//! SVG 描画 (D3: 最小版)。
//!
//! レイアウト座標 (結合長 = 1.0、y 上向き) を [`Style`] で pt にスケールし、
//! y を反転して SVG 座標系に写す。依存ゼロの文字列組み立て
//! (conformer/molblock.rs の流儀)。
//!
//! D3 の範囲: 結合線 (単/二重/三重)、ヘテロ原子・電荷・孤立炭素のラベル
//! (プレーンテキスト)、ラベル付き原子での結合線の短縮 (素朴版)。
//! ラベルの bbox クリップ・下付き数字・H 左置き・二重結合の sidedness・
//! くさび描画は D10。

use crate::graph::MoleculeGraph;

use super::chain_layout::{hidden_h_counts, visible_adjacency};
use super::point2::Point2;
use super::{Coords2D, Style};

/// f64 を決定的な短い 10 進表記にする (golden スナップショットの安定性)。
fn fmt(v: f64) -> String {
    let s = format!("{v:.2}");
    // "-0.00" → "0.00"
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

/// 原子にラベルを描くか (D3 版: 非炭素、電荷つき、または可視結合を持たない C)。
fn has_label(g: &MoleculeGraph, vadj: &[Vec<usize>], i: usize) -> bool {
    let a = &g.atoms[i];
    a.symbol != "C" || a.formal_charge != 0 || vadj[i].is_empty()
}

/// D3 版ラベル文字列: 元素記号 + Hn + 電荷 (プレーンテキスト)。
fn label_text(g: &MoleculeGraph, h_counts: &[u8], i: usize) -> String {
    let a = &g.atoms[i];
    let mut s = a.symbol.clone();
    match h_counts[i] {
        0 => {}
        1 => s.push('H'),
        n => {
            s.push('H');
            s.push_str(&n.to_string());
        }
    }
    match a.formal_charge {
        0 => {}
        1 => s.push('+'),
        -1 => s.push('-'),
        c if c > 0 => s.push_str(&format!("{c}+")),
        c => s.push_str(&format!("{}-", -c)),
    }
    s
}

/// 結合 1 本を線分列として書き出す。二重/三重は平行線 (センター振り分け)。
fn push_bond_lines(
    out: &mut String,
    p1: Point2,
    p2: Point2,
    order: f64,
    spacing: f64,
    line_width: f64,
) {
    let Some(dir) = (p2 - p1).normalized() else {
        return;
    };
    let perp = dir.perp();
    let offsets: &[f64] = if order >= 3.0 {
        &[-1.0, 0.0, 1.0]
    } else if order >= 2.0 {
        &[-0.5, 0.5]
    } else {
        &[0.0]
    };
    for &k in offsets {
        let o = perp * (k * spacing);
        let (a, b) = (p1 + o, p2 + o);
        out.push_str(&format!(
            "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"black\" stroke-width=\"{}\" stroke-linecap=\"round\"/>\n",
            fmt(a.x), fmt(a.y), fmt(b.x), fmt(b.y), fmt(line_width),
        ));
    }
}

/// 2D 座標を SVG 文字列にする。
pub fn to_svg(g: &MoleculeGraph, c: &Coords2D, s: &Style) -> String {
    let vadj = visible_adjacency(g, &c.hidden);
    let h_counts = hidden_h_counts(g, &c.hidden);
    let scale = s.scale();

    let visible: Vec<usize> = (0..g.atoms.len()).filter(|&i| !c.hidden[i]).collect();

    // レイアウト座標の範囲 → pt 座標系 (y 反転 + パディング)
    let (mut min_x, mut max_x, mut min_y, mut max_y) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for &i in &visible {
        min_x = min_x.min(c.pos[i].x);
        max_x = max_x.max(c.pos[i].x);
        min_y = min_y.min(c.pos[i].y);
        max_y = max_y.max(c.pos[i].y);
    }
    if visible.is_empty() {
        min_x = 0.0;
        max_x = 0.0;
        min_y = 0.0;
        max_y = 0.0;
    }
    let pad = scale * 0.75;
    let width = (max_x - min_x) * scale + 2.0 * pad;
    let height = (max_y - min_y) * scale + 2.0 * pad;
    let to_pt = |p: Point2| -> Point2 {
        Point2::new((p.x - min_x) * scale + pad, (max_y - p.y) * scale + pad)
    };

    let mut out = String::new();
    out.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}pt\" height=\"{h}pt\" viewBox=\"0 0 {w} {h}\">\n",
        w = fmt(width),
        h = fmt(height),
    ));

    // 結合線 (可視原子間のみ)。ラベル付き端点では線を短縮する (素朴版:
    // フォント高の半分ぶん)
    let trim = 0.55 * s.font_size_pt;
    let spacing = s.bond_spacing_frac * scale;
    for (bi, b) in g.bonds.iter().enumerate() {
        if c.hidden[b.begin_idx] || c.hidden[b.end_idx] {
            continue;
        }
        // 芳香族 (1.5) はケクレ形で描く (GR-6: 内円は非推奨)
        let order = g.kekule_bond_orders[bi];
        let mut p1 = to_pt(c.pos[b.begin_idx]);
        let mut p2 = to_pt(c.pos[b.end_idx]);
        let Some(dir) = (p2 - p1).normalized() else {
            continue;
        };
        if has_label(g, &vadj, b.begin_idx) {
            p1 = p1 + dir * trim;
        }
        if has_label(g, &vadj, b.end_idx) {
            p2 = p2 - dir * trim;
        }
        push_bond_lines(&mut out, p1, p2, order, spacing, s.line_width_pt);
    }

    // 原子ラベル
    for &i in &visible {
        if !has_label(g, &vadj, i) {
            continue;
        }
        let p = to_pt(c.pos[i]);
        out.push_str(&format!(
            "  <text x=\"{}\" y=\"{}\" font-family=\"{}\" font-size=\"{}\" text-anchor=\"middle\" dominant-baseline=\"central\">{}</text>\n",
            fmt(p.x),
            fmt(p.y),
            xml_escape(s.font_family),
            fmt(s.font_size_pt),
            xml_escape(&label_text(g, &h_counts, i)),
        ));
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
        assert!(svg.contains(">OH</text>"));
        assert_eq!(svg.matches("<line ").count(), 2); // C-C, C-O
    }

    #[test]
    fn double_and_triple_bond_line_counts() {
        assert_eq!(svg_of("C=C").matches("<line ").count(), 2);
        assert_eq!(svg_of("C#C").matches("<line ").count(), 3);
    }

    #[test]
    fn methane_is_single_label() {
        let svg = svg_of("C");
        assert!(svg.contains(">CH4</text>"));
        assert_eq!(svg.matches("<line ").count(), 0);
    }

    #[test]
    fn charge_label() {
        let svg = svg_of("C[N+](C)(C)C");
        assert!(svg.contains(">N+</text>"));
    }

    #[test]
    fn deterministic() {
        assert_eq!(svg_of("CC(C)C=O"), svg_of("CC(C)C=O"));
    }
}
