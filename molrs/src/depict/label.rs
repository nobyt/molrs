//! 原子ラベルの組版と bbox 推定 (D10)。
//!
//! フォントエンジンなしで Helvetica の AFM 近似文字幅テーブルを使い、
//! ラベル (元素記号 + Hn 下付き + 電荷上付き、GR-2.1) の実寸を見積もる。
//! H は結合が右から来る場合に左置き (HO–, H₂N–)。bbox は結合線の
//! クリップ (margin_width) に使う。

use super::point2::Point2;
use super::Style;
use crate::graph::MoleculeGraph;

/// Helvetica 近似文字幅 (フォントサイズ比)。
fn char_width(c: char) -> f64 {
    match c {
        'i' | 'j' | 'l' => 0.222,
        'I' => 0.278,
        'f' | 't' => 0.278,
        'r' => 0.333,
        'c' | 's' => 0.5,
        'a' | 'b' | 'd' | 'e' | 'g' | 'h' | 'n' | 'o' | 'p' | 'q' | 'u' | 'v' | 'x' | 'y' | 'z'
        | 'k' | 'w' | 'm' => 0.556,
        '0'..='9' => 0.556,
        '+' | '-' => 0.584,
        'A' | 'B' | 'E' | 'K' | 'P' | 'S' | 'V' | 'X' | 'Y' => 0.667,
        'C' | 'D' | 'H' | 'N' | 'R' | 'U' => 0.722,
        'F' | 'T' | 'Z' => 0.611,
        'G' | 'O' | 'Q' => 0.778,
        'L' => 0.556,
        'M' => 0.833,
        'W' => 0.944,
        _ => 0.6,
    }
}

fn text_width(s: &str) -> f64 {
    s.chars().map(char_width).sum()
}

/// ラベルの 1 パーツ (通常 / 下付き / 上付き)。
pub(crate) struct Run {
    pub text: String,
    /// フォントサイズ比 (下付き・上付きは 0.7)
    pub rel_size: f64,
    /// ベースラインシフト (em、正 = 下)
    pub dy_em: f64,
}

/// 組版済みラベル。
pub(crate) struct AtomLabel {
    pub runs: Vec<Run>,
    /// 全体幅 (pt)
    pub width: f64,
    /// 左端の x オフセット (原子中心からの相対 pt; 主元素記号の中心 =
    /// 原子中心になるように置く)
    pub left_offset: f64,
    /// bbox 半高 (pt)
    pub half_height: f64,
}

/// 原子にラベルを描くか: 非炭素、電荷つき、可視結合を持たない C。
pub(crate) fn has_label(g: &MoleculeGraph, vadj: &[Vec<usize>], i: usize) -> bool {
    let a = &g.atoms[i];
    a.symbol != "C" || a.formal_charge != 0 || vadj[i].is_empty()
}

/// H を左に置くか: 可視隣接の平均方向が右向き (結合が右から来る)。
pub(crate) fn h_on_left(pos: &[Point2], vadj: &[Vec<usize>], i: usize) -> bool {
    if vadj[i].is_empty() {
        return false;
    }
    let mut sum = Point2::ZERO;
    for &nb in &vadj[i] {
        if let Some(u) = (pos[nb] - pos[i]).normalized() {
            sum = sum + u;
        }
    }
    sum.x > 1e-6
}

const SUB_SIZE: f64 = 0.7;
const SUB_DY: f64 = 0.25;
const SUP_DY: f64 = -0.4;

/// ラベルを組版する。
pub(crate) fn build_label(
    symbol: &str,
    n_h: u8,
    charge: i8,
    h_left: bool,
    style: &Style,
) -> AtomLabel {
    let mut runs: Vec<Run> = Vec::new();
    let push = |runs: &mut Vec<Run>, text: String, rel: f64, dy: f64| {
        runs.push(Run {
            text,
            rel_size: rel,
            dy_em: dy,
        });
    };

    let h_part = |runs: &mut Vec<Run>| {
        if n_h >= 1 {
            push(runs, "H".into(), 1.0, 0.0);
            if n_h >= 2 {
                push(runs, n_h.to_string(), SUB_SIZE, SUB_DY);
            }
        }
    };

    if h_left {
        h_part(&mut runs);
    }
    push(&mut runs, symbol.to_string(), 1.0, 0.0);
    if !h_left {
        h_part(&mut runs);
    }
    match charge {
        0 => {}
        1 => push(&mut runs, "+".into(), SUB_SIZE, SUP_DY),
        -1 => push(&mut runs, "-".into(), SUB_SIZE, SUP_DY),
        c if c > 0 => push(&mut runs, format!("{c}+"), SUB_SIZE, SUP_DY),
        c => push(&mut runs, format!("{}-", -c), SUB_SIZE, SUP_DY),
    }

    // 幅と主記号中心のオフセットを計算
    let fs = style.font_size_pt;
    let mut x = 0.0;
    let mut sym_center = 0.0;
    // 主記号 run のインデックス: h_left なら H(+数字) の後
    let sym_run_idx = if h_left {
        if n_h >= 2 {
            2
        } else if n_h >= 1 {
            1
        } else {
            0
        }
    } else {
        0
    };
    for (k, r) in runs.iter().enumerate() {
        let w = text_width(&r.text) * fs * r.rel_size;
        if k == sym_run_idx {
            sym_center = x + w / 2.0;
        }
        x += w;
    }
    AtomLabel {
        runs,
        width: x,
        left_offset: -sym_center,
        half_height: 0.6 * fs,
    }
}

impl AtomLabel {
    /// 原子中心 (pt) から方向 dir へ出る結合線のクリップ距離。
    /// bbox を margin ぶん膨らませた楕円で近似する。
    pub(crate) fn trim_distance(&self, dir: Point2, margin: f64) -> f64 {
        // bbox: x ∈ [left_offset, left_offset+width], y ∈ ±half_height
        // 中心非対称なので方向側の半幅を使う
        let hw = if dir.x >= 0.0 {
            (self.left_offset + self.width).abs()
        } else {
            self.left_offset.abs()
        }
        .max(0.3 * self.half_height)
            + margin;
        let hh = self.half_height + margin;
        let (dx, dy) = (dir.x, dir.y);
        let denom = (hh * dx).powi(2) + (hw * dy).powi(2);
        if denom < 1e-12 {
            return hw.min(hh);
        }
        hw * hh / denom.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn widths_monotone() {
        let s = Style::acs_1996();
        let l1 = build_label("O", 1, 0, false, &s); // OH
        let l2 = build_label("N", 2, 0, false, &s); // NH2
        assert!(l2.width > l1.width);
        assert!(l1.width > 0.0);
    }

    #[test]
    fn h_left_layout() {
        let s = Style::acs_1996();
        let l = build_label("O", 1, 0, true, &s); // HO
        assert_eq!(l.runs[0].text, "H");
        assert_eq!(l.runs[1].text, "O");
        // 主記号 O の中心が原点 → 左オフセットは H の幅ぶん大きい
        assert!(l.left_offset < -0.5 * char_width('O') * s.font_size_pt * 0.9);
    }

    #[test]
    fn charge_superscript() {
        let s = Style::acs_1996();
        let l = build_label("N", 0, 1, false, &s); // N+
        assert_eq!(l.runs.last().unwrap().text, "+");
        assert!(l.runs.last().unwrap().dy_em < 0.0);
    }

    #[test]
    fn trim_distance_reasonable() {
        let s = Style::acs_1996();
        let l = build_label("O", 1, 0, false, &s);
        let t_right = l.trim_distance(Point2::new(1.0, 0.0), s.margin_width_pt);
        let t_up = l.trim_distance(Point2::new(0.0, 1.0), s.margin_width_pt);
        assert!(t_right > t_up * 0.8); // OH は横に長い
        assert!(t_right < s.font_size_pt * 2.5);
    }
}
