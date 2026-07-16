//! smiles2iupac: SMILES から IUPAC 2013 優先名を生成するライブラリ (Rust 移植版)。
//!
//! Python 実装 (`src/smiles2iupac/`) の移植。移植計画は RUST_PORT_PLAN.md を参照。
//!
//! 現在の対応範囲 (S2.2 パイプライン骨格):
//! - 保留名テーブル (立体署名付き複合キーで照合)
//! - 直鎖アルカン (C1〜C30)
//!
//! 未対応の構造は [`NameError::Unsupported`] を返す (移植の進行に応じて拡大)。

pub mod constants;
pub mod functional_group;
pub mod retained;

use molrs::graph::{build_molecule_graph, MoleculeGraph};
pub use molrs::ChemError;

/// 命名エラー。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NameError {
    /// 不正な SMILES
    InvalidSmiles(String),
    /// この構造はまだ移植されていない (移植完了までの暫定エラー)
    Unsupported(String),
}

impl std::fmt::Display for NameError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NameError::InvalidSmiles(s) => write!(f, "Invalid SMILES: {s}"),
            NameError::Unsupported(s) => write!(f, "Unsupported structure: {s}"),
        }
    }
}

impl std::error::Error for NameError {}

/// SMILES 文字列を IUPAC 2013 優先名に変換する。
///
/// Python 版 `smiles_to_iupac()` に対応。
pub fn smiles_to_iupac(smiles: &str) -> Result<String, NameError> {
    let smiles = smiles.trim();
    if smiles.contains('.') {
        // 多成分 (塩・混合物) は S5.8 で対応
        return Err(NameError::Unsupported(format!("multicomponent: {smiles}")));
    }
    let graph =
        build_molecule_graph(smiles).map_err(|e| NameError::InvalidSmiles(e.to_string()))?;

    // 保留名テーブル (Python _try_retained_name 相当)
    if let Some(name) = retained::try_retained_name(&graph) {
        return Ok(name.to_string());
    }

    // 直鎖アルカン (パイプライン骨格の最小ケース)
    if let Some(name) = try_linear_alkane(&graph) {
        return Ok(name);
    }

    Err(NameError::Unsupported(smiles.to_string()))
}

/// 直鎖アルカン (メタン〜トリアコンタン) の命名。
/// 全原子が中性の炭素、環なし、分岐なし、単結合のみの場合のみ Some。
fn try_linear_alkane(g: &MoleculeGraph) -> Option<String> {
    let n_kept = g.parser_to_graph.iter().flatten().count();
    if n_kept == 0 || n_kept > 30 {
        return None;
    }
    let heavy_degree = |i: usize| g.adjacency[i].iter().filter(|&&x| x < n_kept).count();
    for i in 0..n_kept {
        let a = &g.atoms[i];
        if a.symbol != "C" || a.formal_charge != 0 || a.in_ring || a.is_aromatic {
            return None;
        }
        if heavy_degree(i) > 2 {
            return None;
        }
    }
    // 全結合が単結合
    for b in &g.bonds {
        if b.begin_idx < n_kept && b.end_idx < n_kept && b.bond_order != 1.0 {
            return None;
        }
    }
    Some(format!("{}ane", constants::CHAIN_PREFIX[n_kept]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_alkanes() {
        assert_eq!(smiles_to_iupac("C").unwrap(), "methane");
        assert_eq!(smiles_to_iupac("CC").unwrap(), "ethane");
        assert_eq!(smiles_to_iupac("CCC").unwrap(), "propane");
        assert_eq!(smiles_to_iupac("C(C)C").unwrap(), "propane"); // 表記違い
        assert_eq!(
            smiles_to_iupac("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC").unwrap(),
            "triacontane"
        );
    }

    #[test]
    fn not_linear_alkanes() {
        assert!(matches!(
            smiles_to_iupac("CC(C)C"),
            Err(NameError::Unsupported(_)) // 分岐 (S3 で対応)
        ));
        assert!(matches!(
            smiles_to_iupac("C1CC1"),
            Err(NameError::Unsupported(_)) // 環 (S4 で対応)
        ));
        assert!(matches!(
            smiles_to_iupac("C=C"),
            Err(NameError::Unsupported(_)) // アルケンは保留名/S3 経由
        ));
    }

    #[test]
    fn invalid_smiles() {
        assert!(matches!(
            smiles_to_iupac("hello"),
            Err(NameError::InvalidSmiles(_))
        ));
    }

    #[test]
    fn constants_table_counts() {
        // 生成テーブルの件数が Python 側と一致すること
        assert_eq!(constants::FUNCTIONAL_GROUPS.len(), 165);
        assert_eq!(constants::CHAIN_PREFIX.len(), 31);
        assert_eq!(constants::RETAINED_NAMES_RAW.len(), 205);
        assert_eq!(constants::CHAIN_PREFIX[1], "meth");
        assert_eq!(constants::CHAIN_PREFIX[30], "triacont");
        assert_eq!(constants::MULTIPLIER[2], "di");
        assert_eq!(constants::halogen_prefix("Cl"), Some("chloro"));
        assert_eq!(constants::suffix_of("carboxylic_acid"), Some("oic acid"));
    }
}
