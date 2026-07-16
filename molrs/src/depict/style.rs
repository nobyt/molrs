//! 描画スタイル (D1)。
//!
//! ChemDraw のドキュメント設定と同じ語彙を持つ。レイアウト
//! ([`compute_coords_2d`](super::compute_coords_2d)) は結合長 = 1.0 の
//! 無次元座標を生成し、描画時に本構造体がスケール・線幅・フォントを決める
//! (レイアウトとスタイルの直交性)。
//!
//! プリセット値の出典 (RUST_2D_PLAN.md に詳細):
//! - ACS: pubs.acs.org "Preparing Graphics" (ACS Document 1996 設定)
//! - Nature: nature.com/documents/nr-chemical-structures-guide.pdf
//! - RSC: rsc.org Chemical Science 投稿ガイド
//! - Wiley: 公式の数値一覧が確認できず、ChemDraw「Wiley Document」スタイル
//!   シート = ACS 同値との二次情報に基づく (**要確認**)

/// 1 インチ = 72 pt。cm 由来の値は 1 cm = 72/2.54 pt で換算。
pub const PT_PER_INCH: f64 = 72.0;

#[derive(Debug, Clone, PartialEq)]
pub struct Style {
    /// 結合長 (pt)。ChemDraw "fixed length"
    pub bond_length_pt: f64,
    /// 二重結合の 2 本目の線との間隔 (結合長に対する比)。ChemDraw "bond spacing"
    pub bond_spacing_frac: f64,
    /// 結合線の太さ (pt)。ChemDraw "line width"
    pub line_width_pt: f64,
    /// 太線・くさびの最大幅 (pt)。ChemDraw "bold width"
    pub bold_width_pt: f64,
    /// 破線くさびの線間隔 (pt)。ChemDraw "hash spacing"
    pub hash_spacing_pt: f64,
    /// 原子ラベル周囲の空白 (pt)。結合線はこの距離でクリップされる。
    /// ChemDraw "margin width"
    pub margin_width_pt: f64,
    /// 原子ラベルのフォントファミリ (SVG font-family 値)
    pub font_family: &'static str,
    /// 原子ラベルのフォントサイズ (pt)
    pub font_size_pt: f64,
    /// 単段組の図幅上限 (in)。超過してもエラーにはせず利用側の判断材料
    pub max_width_single_col_in: Option<f64>,
    /// 二段組の図幅上限 (in)
    pub max_width_double_col_in: Option<f64>,
}

impl Style {
    /// IUPAC 2008 勧告準拠の既定スタイル。勧告は絶対寸法を定めない
    /// (GR-0.6: 媒体に応じて可読なら可) ため、寸法は ACS 1996 と同値を使う。
    pub fn iupac_default() -> Style {
        Style {
            bond_length_pt: 14.4,
            bond_spacing_frac: 0.18,
            line_width_pt: 0.6,
            bold_width_pt: 2.0,
            hash_spacing_pt: 2.5,
            margin_width_pt: 1.6,
            font_family: "Helvetica, Arial, sans-serif",
            font_size_pt: 10.0,
            max_width_single_col_in: None,
            max_width_double_col_in: None,
        }
    }

    /// ACS Document 1996 (J. Am. Chem. Soc. 等 ACS 全誌)。
    pub fn acs_1996() -> Style {
        Style {
            max_width_single_col_in: Some(3.25),
            max_width_double_col_in: Some(7.0),
            ..Style::iupac_default()
        }
    }

    /// Nature Research 系 (Style guide for chemical structures)。
    /// cm 規定値の pt 換算: fixed 0.381 cm、line 0.021 cm、bold 0.055 cm、
    /// hash 0.06 cm、margin 0.042 cm。ラベルは 6 pt。
    pub fn nature() -> Style {
        const PT_PER_CM: f64 = PT_PER_INCH / 2.54;
        Style {
            bond_length_pt: 0.381 * PT_PER_CM, // 10.80 pt
            bond_spacing_frac: 0.18,
            line_width_pt: 0.021 * PT_PER_CM,   // 0.595 pt
            bold_width_pt: 0.055 * PT_PER_CM,   // 1.559 pt
            hash_spacing_pt: 0.06 * PT_PER_CM,  // 1.701 pt
            margin_width_pt: 0.042 * PT_PER_CM, // 1.191 pt
            font_family: "Arial, Helvetica, sans-serif",
            font_size_pt: 6.0,
            max_width_single_col_in: None,
            max_width_double_col_in: None,
        }
    }

    /// RSC (Chemical Science 等)。bond 12.2 pt、二重結合間隔 20%、
    /// 線 0.5 pt、太線/くさび 1.6 pt、hash 1.8 pt、ラベル 7 pt。
    /// margin は規定がないため ACS 値を流用。
    pub fn rsc() -> Style {
        Style {
            bond_length_pt: 12.2,
            bond_spacing_frac: 0.20,
            line_width_pt: 0.5,
            bold_width_pt: 1.6,
            hash_spacing_pt: 1.8,
            margin_width_pt: 1.6,
            font_family: "Arial, Helvetica, sans-serif",
            font_size_pt: 7.0,
            max_width_single_col_in: None,
            max_width_double_col_in: None,
        }
    }

    /// Wiley (Angewandte Chemie 等)。公式一覧が確認できないため
    /// ChemDraw「Wiley Document」= ACS 同値との情報に基づく (要確認)。
    pub fn wiley() -> Style {
        Style::iupac_default()
    }

    /// レイアウト単位 (結合長 = 1.0) → pt のスケール係数。
    pub fn scale(&self) -> f64 {
        self.bond_length_pt
    }
}

impl Default for Style {
    fn default() -> Style {
        Style::iupac_default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn acs_values() {
        let s = Style::acs_1996();
        assert_eq!(s.bond_length_pt, 14.4);
        assert_eq!(s.bond_spacing_frac, 0.18);
        assert_eq!(s.line_width_pt, 0.6);
        assert_eq!(s.bold_width_pt, 2.0);
        assert_eq!(s.hash_spacing_pt, 2.5);
        assert_eq!(s.margin_width_pt, 1.6);
        assert_eq!(s.font_size_pt, 10.0);
        assert_eq!(s.max_width_single_col_in, Some(3.25));
        assert_eq!(s.max_width_double_col_in, Some(7.0));
    }

    #[test]
    fn nature_values() {
        let s = Style::nature();
        assert!((s.bond_length_pt - 10.7999).abs() < 1e-3);
        assert!((s.line_width_pt - 0.5953).abs() < 1e-3);
        assert!((s.bold_width_pt - 1.5591).abs() < 1e-3);
        assert!((s.hash_spacing_pt - 1.7008).abs() < 1e-3);
        assert!((s.margin_width_pt - 1.1906).abs() < 1e-3);
        assert_eq!(s.font_size_pt, 6.0);
    }

    #[test]
    fn rsc_values() {
        let s = Style::rsc();
        assert_eq!(s.bond_length_pt, 12.2);
        assert_eq!(s.bond_spacing_frac, 0.20);
        assert_eq!(s.line_width_pt, 0.5);
        assert_eq!(s.bold_width_pt, 1.6);
        assert_eq!(s.hash_spacing_pt, 1.8);
        assert_eq!(s.font_size_pt, 7.0);
    }

    #[test]
    fn wiley_matches_acs_drawing_values() {
        let w = Style::wiley();
        let a = Style::acs_1996();
        assert_eq!(w.bond_length_pt, a.bond_length_pt);
        assert_eq!(w.line_width_pt, a.line_width_pt);
    }

    #[test]
    fn default_is_iupac() {
        assert_eq!(Style::default(), Style::iupac_default());
    }
}
