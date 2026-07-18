//! SMILES/分子グラフ → 標準 InChI (`InChI=1S/…`) と InChIKey。
//!
//! IUPAC 公式 InChI (C 実装) とビット完全一致を目標とする。計画:
//! RUST_INCHI_PLAN.md。
//!
//! v1 の対象: 骨格層 (式・接続 `c`・水素 `h`・電荷 `q`/`p`) と InChIKey。
//! 立体 (`b`/`t`/`m`/`s`)・同位体 (`i`)・一般の互変異性正規化は未対応。
//!
//! InChIKey は標準 InChI 文字列の SHA-256 ハッシュ (base-26 符号化) なので、
//! 依存クレートゼロを保つため SHA-256 を自前実装している ([`sha256`])。

pub(crate) mod formula;
pub mod number;
pub mod sha256;

use crate::graph::MoleculeGraph;

/// InChI 生成のエラー。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InchiError {
    /// 不正な SMILES ([`inchi_of`] 経由のみ)
    InvalidSmiles(String),
    /// v1 で未対応の構造クラス (実装の進行に応じて縮小)
    Unsupported(String),
}

impl std::fmt::Display for InchiError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            InchiError::InvalidSmiles(s) => write!(f, "invalid SMILES: {s}"),
            InchiError::Unsupported(s) => write!(f, "InChI unsupported: {s}"),
        }
    }
}

impl std::error::Error for InchiError {}

/// 分子グラフの Hill 式層を返す (I2)。式のみが必要な場合の公開 API。
pub fn formula(g: &MoleculeGraph) -> String {
    formula::formula_layer(g)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    #[test]
    fn formula_public_api() {
        let g = build_molecule_graph("CC(=O)O").unwrap();
        assert_eq!(formula(&g), "C2H4O2");
    }
}
