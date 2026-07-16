//! SMILES パーサ (S1.1)。
//!
//! SMILES 文字列を [`ParsedMolecule`] (原子リスト + 結合リスト) に変換する。
//! 原子価チェック・暗黙 H 計算・芳香族認識は S1.2/S1.3 の担当で、ここでは行わない。

mod parser;

pub use parser::parse_smiles;

/// SMILES 中の四面体キラリティ表記。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Chirality {
    /// `@` (反時計回り)
    Anticlockwise,
    /// `@@` (時計回り)
    Clockwise,
}

/// パース直後の原子。記号は正規の大文字小文字 ("C", "Cl", "Se", ...) に正規化済み。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomSpec {
    /// 元素記号 (正規化済み)。ワイルドカードは "*"。
    pub symbol: String,
    /// 入力で芳香族小文字表記だったか
    pub aromatic: bool,
    /// 同位体 (`[13C]` → Some(13))
    pub isotope: Option<u16>,
    /// 形式電荷
    pub charge: i8,
    /// 角括弧内の明示 H 数 (`[CH3]` → Some(3))。角括弧原子で H 指定なしは Some(0)、
    /// 有機サブセット原子 (暗黙 H を原子価から計算する) は None。
    pub explicit_h: Option<u8>,
    /// `@` / `@@`
    pub chirality: Option<Chirality>,
    /// アトムクラス (`[CH4:2]` → Some(2))
    pub atom_class: Option<u32>,
    /// 角括弧原子だったか
    pub bracket: bool,
}

/// SMILES 上の結合種別。Elided (省略) の解決は S1.2/S1.3 で行う。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondKind {
    /// 記号なし (単結合または芳香族結合として後段で解決)
    Elided,
    /// `-`
    Single,
    /// `=`
    Double,
    /// `#`
    Triple,
    /// `$`
    Quadruple,
    /// `:`
    Aromatic,
    /// `/` (cis/trans 用方向付き単結合)
    Up,
    /// `\`
    Down,
}

/// 環閉じ結合の付帯情報。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RingClosure {
    /// 環番号 (RDKit 互換の結合順: 環閉じ結合は末尾に環番号順で並ぶ)
    pub num: u16,
    /// 開き側の数字に結合次数記号 (-, =, #, $, :) が付いていたか。
    /// RDKit はこの場合のみ結合を (開き側, 閉じ側) の向きで格納する
    /// (それ以外は (閉じ側, 開き側))。方向記号 / \ は次数扱いしない。
    pub opened_with_order: bool,
}

/// パース直後の結合。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BondSpec {
    /// 原子インデックス (SMILES 出現順)。環閉じ結合では a = 開き側, b = 閉じ側で、
    /// `kind` の Up/Down は a→b 向きの意味を持つ。
    pub a: usize,
    pub b: usize,
    pub kind: BondKind,
    /// 環閉じ数字由来なら Some。
    pub ring_closure: Option<RingClosure>,
}

/// SMILES パース結果。
///
/// `neighbor_order[i]` は原子 i に接続する結合のインデックスを
/// **SMILES 出現順** に並べたもの。環閉じ結合は数字が現れた位置に置かれる。
/// `@`/`@@` の解釈 (S1.7) はこの順序に依存する。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsedMolecule {
    pub atoms: Vec<AtomSpec>,
    pub bonds: Vec<BondSpec>,
    pub neighbor_order: Vec<Vec<usize>>,
}
