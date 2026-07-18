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

pub mod base26;
pub(crate) mod formula;
pub(crate) mod layers;
pub mod number;
pub mod sha256;

pub use base26::inchi_key_from_string;

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

/// 標準 InChI (`InChI=1S/…`) を生成する (I4、v1 範囲)。
///
/// v1 は単一成分・中性 (電荷 q/p・立体・同位体・多成分・有機金属は未対応)。
/// 適用範囲外は [`InchiError::Unsupported`] を返す。
pub fn to_inchi(g: &MoleculeGraph) -> Result<String, InchiError> {
    let comps = layers::build_components(g);
    if comps.len() != 1 {
        return Err(InchiError::Unsupported("multi-component (v2)".into()));
    }
    // 電荷を持つ分子は q/p 層が要る (v1 未対応)。可動群で中和される
    // 負電荷は許容 (p 層は今後)。ここでは全原子中性のみ通す。
    if g.atoms.iter().any(|a| a.formal_charge != 0) {
        return Err(InchiError::Unsupported(
            "charged (q/p layer pending)".into(),
        ));
    }

    let formula = formula::formula_layer(g);
    let c = layers::connection_layer(&comps[0]);
    let h = layers::hydrogen_layer(&comps[0]);

    let mut s = format!("InChI=1S/{formula}");
    if !c.is_empty() {
        s.push_str("/c");
        s.push_str(&c);
    }
    if !h.is_empty() {
        s.push_str("/h");
        s.push_str(&h);
    }
    Ok(s)
}

/// SMILES から標準 InChI を生成する便利関数。
pub fn inchi_of(smiles: &str) -> Result<String, InchiError> {
    let g = crate::graph::build_molecule_graph(smiles)
        .map_err(|e| InchiError::InvalidSmiles(e.to_string()))?;
    to_inchi(&g)
}

/// 分子グラフの InChIKey を生成する (I5、v1 範囲)。
pub fn to_inchi_key(g: &MoleculeGraph) -> Result<String, InchiError> {
    let inchi = to_inchi(g)?;
    Ok(inchi_key_from_string(&inchi))
}

/// SMILES から InChIKey を生成する便利関数。
pub fn inchi_key_of(smiles: &str) -> Result<String, InchiError> {
    let inchi = inchi_of(smiles)?;
    Ok(inchi_key_from_string(&inchi))
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
