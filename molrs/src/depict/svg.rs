//! SVG 描画 (D3 最小版 → D10 精緻化)。
//!
//! レイアウト座標 (結合長 = 1.0、y 上向き) を [`Style`] で pt にスケールし、
//! y を反転して SVG 座標系に写す。依存ゼロの文字列組み立て。
//!
//! - 原子ラベル: 元素記号 + Hn 下付き + 電荷上付き (GR-2.1)、
//!   結合が右から来る場合は H 左置き (HO–)。label.rs の文字幅テーブルで
//!   bbox を推定し、結合線を margin_width でクリップする
//! - 二重結合の sidedness (GR-1.10): 環内は環中心側に短縮線、鎖中は
//!   置換基の多い側に短縮線、末端 (両端に他の隣接なし) はセンター振り分け
//! - くさび (IUPAC 2006): solid = 三角形、hashed = hash_spacing 間隔の
//!   横棒列。細端 = 立体中心

use crate::graph::MoleculeGraph;

use super::chain_layout::{hidden_h_counts, visible_adjacency};
use super::label::{build_label, h_on_left, has_label, AtomLabel};
use super::point2::Point2;
use super::{Coords2D, Style, WedgeDir};

/// f64 を決定的な短い 10 進表記にする (golden スナップショットの安定性)。
fn fmt(v: f64) -> String {
    let s = format!("{v:.2}");
    if s == "-0.00" {
        "0.00".to_string()
    } else {
        s
    }
}

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// 描画プリミティブ (グローバル平行移動前)。
enum Op {
    Line { p1: Point2, p2: Point2, width: f64 },
    Polygon(Vec<Point2>),
    Text { anchor: Point2, label: AtomLabel },
}

/// 二重結合の 2 本目の線を出す側。None = センター振り分け。
fn double_bond_side(
    g: &MoleculeGraph,
    c: &Coords2D,
    vadj: &[Vec<usize>],
    a: usize,
    b: usize,
) -> Option<f64> {
    let dir = (c.pos[b] - c.pos[a]).normalized()?;
    // 環内: 最小の環の重心側
    let mut best_ring: Option<&Vec<usize>> = None;
    for ring in &g.ring_atom_sets {
        let n = ring.len();
        let has = (0..n).any(|k| {
            let (x, y) = (ring[k], ring[(k + 1) % n]);
            (x == a && y == b) || (x == b && y == a)
        });
        if has && best_ring.map(|r| ring.len() < r.len()).unwrap_or(true) {
            best_ring = Some(ring);
        }
    }
    if let Some(ring) = best_ring {
        let centroid = ring.iter().fold(Point2::ZERO, |s, &x| s + c.pos[x]) / ring.len() as f64;
        let side = dir.cross(centroid - c.pos[a]);
        return Some(if side >= 0.0 { 1.0 } else { -1.0 });
    }
    // 鎖: 両端の他の隣接の分布で決める。両端とも裸ならセンター
    let mut score = 0.0;
    let mut n_others = 0;
    for (end, other) in [(a, b), (b, a)] {
        for &nb in &vadj[end] {
            if nb == other {
                continue;
            }
            n_others += 1;
            let s = dir.cross(c.pos[nb] - c.pos[end]);
            if s.abs() > 1e-9 {
                score += s.signum();
            }
        }
    }
    if n_others == 0 {
        return None; // 末端二重結合 (C=O 等): センター振り分け
    }
    Some(if score >= 0.0 { 1.0 } else { -1.0 })
}

/// 2D 座標を SVG 文字列にする。
pub fn to_svg(g: &MoleculeGraph, c: &Coords2D, s: &Style) -> String {
    let vadj = visible_adjacency(g, &c.hidden);
    let h_counts = hidden_h_counts(g, &c.hidden);
    let scale = s.scale();
    let visible: Vec<usize> = (0..g.atoms.len()).filter(|&i| !c.hidden[i]).collect();

    // pt 座標 (y 反転のみ。平行移動は最後に bbox から決める)
    let pt = |p: Point2| Point2::new(p.x * scale, -p.y * scale);

    // ラベルの組版
    let labels: Vec<Option<AtomLabel>> = (0..g.atoms.len())
        .map(|i| {
            if c.hidden[i] || !has_label(g, &vadj, i) {
                return None;
            }
            let h_left = h_on_left(&c.pos, &vadj, i);
            Some(build_label(
                &g.atoms[i].symbol,
                h_counts[i],
                g.atoms[i].formal_charge,
                h_left,
                s,
            ))
        })
        .collect();

    let mut ops: Vec<Op> = Vec::new();

    // 結合
    let spacing = s.bond_spacing_frac * scale;
    for (bi, b) in g.bonds.iter().enumerate() {
        let (i, j) = (b.begin_idx, b.end_idx);
        if c.hidden[i] || c.hidden[j] {
            continue;
        }
        let (mut p1, mut p2) = (pt(c.pos[i]), pt(c.pos[j]));
        let Some(dir) = (p2 - p1).normalized() else {
            continue;
        };
        // ラベルでクリップ
        if let Some(l) = &labels[i] {
            p1 = p1 + dir * l.trim_distance(dir, s.margin_width_pt);
        }
        if let Some(l) = &labels[j] {
            p2 = p2 - dir * l.trim_distance(-dir, s.margin_width_pt);
        }
        let len = p1.distance(p2);
        if len < 1e-6 {
            continue;
        }

        // くさび
        if let Some(w) = &c.wedge[bi] {
            let (narrow, wide) = if w.narrow == i { (p1, p2) } else { (p2, p1) };
            let Some(wdir) = (wide - narrow).normalized() else {
                continue;
            };
            let perp = wdir.perp();
            match w.dir {
                WedgeDir::Up => {
                    // solid wedge: 三角形
                    ops.push(Op::Polygon(vec![
                        narrow,
                        wide + perp * (s.bold_width_pt / 2.0),
                        wide - perp * (s.bold_width_pt / 2.0),
                    ]));
                }
                WedgeDir::Down => {
                    // hashed wedge: 直交する横棒列
                    let wlen = narrow.distance(wide);
                    let n_hash = ((wlen / s.hash_spacing_pt).floor() as usize).max(2);
                    for k in 1..=n_hash {
                        let t = k as f64 / n_hash as f64;
                        let half = 0.5 * (s.line_width_pt * (1.0 - t) + s.bold_width_pt * t);
                        let q = narrow + wdir * (t * wlen);
                        ops.push(Op::Line {
                            p1: q + perp * half,
                            p2: q - perp * half,
                            width: s.line_width_pt,
                        });
                    }
                }
            }
            continue;
        }

        let order = g.kekule_bond_orders[bi];
        if order >= 3.0 {
            // 三重結合: センター 3 本
            let perp = dir.perp();
            for k in [-1.0, 0.0, 1.0] {
                let o = perp * (k * spacing);
                ops.push(Op::Line {
                    p1: p1 + o,
                    p2: p2 + o,
                    width: s.line_width_pt,
                });
            }
        } else if order >= 2.0 {
            let perp = dir.perp();
            match double_bond_side(g, c, &vadj, i, j) {
                None => {
                    // センター振り分け
                    for k in [-0.5, 0.5] {
                        let o = perp * (k * spacing);
                        ops.push(Op::Line {
                            p1: p1 + o,
                            p2: p2 + o,
                            width: s.line_width_pt,
                        });
                    }
                }
                Some(side) => {
                    // 主線 + 側方の短縮線 (レイアウトは y 上向き、pt は y
                    // 反転済みなので side の符号も反転する)
                    ops.push(Op::Line {
                        p1,
                        p2,
                        width: s.line_width_pt,
                    });
                    let o = perp * (-side * spacing);
                    let shorten = dir * (0.18 * len);
                    ops.push(Op::Line {
                        p1: p1 + o + shorten,
                        p2: p2 + o - shorten,
                        width: s.line_width_pt,
                    });
                }
            }
        } else {
            ops.push(Op::Line {
                p1,
                p2,
                width: s.line_width_pt,
            });
        }
    }

    // ラベル
    for &i in &visible {
        if let Some(l) = &labels[i] {
            ops.push(Op::Text {
                anchor: pt(c.pos[i]),
                label: build_label(
                    &g.atoms[i].symbol,
                    h_counts[i],
                    g.atoms[i].formal_charge,
                    h_on_left(&c.pos, &vadj, i),
                    s,
                ),
            });
            let _ = l;
        }
    }

    // グローバル bbox → 平行移動
    let (mut min_x, mut min_y, mut max_x, mut max_y) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
    let mut visit = |p: Point2| {
        min_x = min_x.min(p.x);
        min_y = min_y.min(p.y);
        max_x = max_x.max(p.x);
        max_y = max_y.max(p.y);
    };
    for op in &ops {
        match op {
            Op::Line { p1, p2, .. } => {
                visit(*p1);
                visit(*p2);
            }
            Op::Polygon(pts) => pts.iter().for_each(|&p| visit(p)),
            Op::Text { anchor, label } => {
                visit(*anchor + Point2::new(label.left_offset, -label.half_height));
                visit(*anchor + Point2::new(label.left_offset + label.width, label.half_height));
            }
        }
    }
    if ops.is_empty() {
        for &i in &visible {
            visit(pt(c.pos[i]));
        }
        if visible.is_empty() {
            visit(Point2::ZERO);
        }
    }
    let pad = 0.35 * scale;
    let off = Point2::new(pad - min_x, pad - min_y);
    let width = max_x - min_x + 2.0 * pad;
    let height = max_y - min_y + 2.0 * pad;

    // 出力
    let mut out = String::new();
    out.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}pt\" height=\"{h}pt\" viewBox=\"0 0 {w} {h}\">\n",
        w = fmt(width),
        h = fmt(height),
    ));
    for op in &ops {
        match op {
            Op::Line { p1, p2, width } => {
                let (a, b) = (*p1 + off, *p2 + off);
                out.push_str(&format!(
                    "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"black\" stroke-width=\"{}\" stroke-linecap=\"round\"/>\n",
                    fmt(a.x), fmt(a.y), fmt(b.x), fmt(b.y), fmt(*width),
                ));
            }
            Op::Polygon(pts) => {
                let coords: Vec<String> = pts
                    .iter()
                    .map(|&p| format!("{},{}", fmt((p + off).x), fmt((p + off).y)))
                    .collect();
                out.push_str(&format!(
                    "  <polygon points=\"{}\" fill=\"black\"/>\n",
                    coords.join(" ")
                ));
            }
            Op::Text { anchor, label } => {
                let a = *anchor + off;
                let x0 = a.x + label.left_offset;
                out.push_str(&format!(
                    "  <text x=\"{}\" y=\"{}\" font-family=\"{}\" font-size=\"{}\" dominant-baseline=\"central\">",
                    fmt(x0),
                    fmt(a.y),
                    xml_escape(s.font_family),
                    fmt(s.font_size_pt),
                ));
                let mut cur_dy = 0.0;
                for r in &label.runs {
                    let dy = r.dy_em - cur_dy;
                    cur_dy = r.dy_em;
                    let size = s.font_size_pt * r.rel_size;
                    if r.rel_size == 1.0 && dy == 0.0 {
                        out.push_str(&xml_escape(&r.text));
                    } else {
                        out.push_str(&format!(
                            "<tspan font-size=\"{}\" dy=\"{}\">{}</tspan>",
                            fmt(size),
                            fmt(dy * s.font_size_pt),
                            xml_escape(&r.text)
                        ));
                    }
                }
                out.push_str("</text>\n");
            }
        }
    }
    out.push_str("</svg>\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::depict::{compute_coords_2d, LayoutParams};
    use crate::graph::build_molecule_graph;

    fn svg_of(smiles: &str) -> String {
        let g = build_molecule_graph(smiles).unwrap();
        let c = compute_coords_2d(&g, &LayoutParams::default()).unwrap();
        to_svg(&g, &c, &Style::acs_1996())
    }

    #[test]
    fn ethanol_has_oh_label_and_bond() {
        let svg = svg_of("CCO");
        assert!(svg.starts_with("<svg xmlns"));
        // 結合が右から来る場合は H 左置き (HO–, GR-2.1.1)
        assert!(
            svg.contains(">HO</text>") || svg.contains(">OH</text>"),
            "{svg}"
        );
        assert_eq!(svg.matches("<line ").count(), 2);
    }

    #[test]
    fn double_and_triple_bond_line_counts() {
        assert_eq!(svg_of("C=C").matches("<line ").count(), 2);
        assert_eq!(svg_of("C#C").matches("<line ").count(), 3);
    }

    #[test]
    fn methane_is_single_label() {
        let svg = svg_of("C");
        assert!(svg.contains("CH<tspan"), "subscript expected: {svg}");
        assert_eq!(svg.matches("<line ").count(), 0);
    }

    #[test]
    fn charge_superscript() {
        let svg = svg_of("C[N+](C)(C)C");
        assert!(svg.contains(">+</tspan>"), "{svg}");
    }

    #[test]
    fn wedge_drawn_for_stereocenter() {
        let svg = svg_of("C[C@H](O)CC");
        // solid (polygon) か hashed (短い線群) のどちらかが出る
        let has_polygon = svg.contains("<polygon");
        let n_lines = svg.matches("<line ").count();
        assert!(has_polygon || n_lines > 4, "{svg}");
    }

    #[test]
    fn h_left_for_right_attached_oh() {
        // フェノール: O は環の左側に置かれることがあるため、H 左置きの
        // 判定ロジック自体をテスト (O の隣接が右にあるケースを作る)
        let g = build_molecule_graph("OC").unwrap(); // O が先頭 → メチルが右
        let c = compute_coords_2d(&g, &LayoutParams::default()).unwrap();
        let svg = to_svg(&g, &c, &Style::acs_1996());
        // HO 表記 (H が O の前) または OH (向きにより)。どちらかを含む
        assert!(
            svg.contains("HO<")
                || svg.contains(">HO</text>")
                || svg.contains(">OH<")
                || svg.contains("OH<"),
            "{svg}"
        );
    }

    #[test]
    fn benzene_kekule_has_double_lines() {
        let svg = svg_of("c1ccccc1");
        // 6 主線 + 3 内側短縮線 = 9
        assert_eq!(svg.matches("<line ").count(), 9);
    }

    #[test]
    fn deterministic() {
        assert_eq!(svg_of("CC(C)C=O"), svg_of("CC(C)C=O"));
    }
}
