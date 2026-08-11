//! molrs: smiles2iupac 用の最小ケモインフォマティクスコア (RDKit 代替層)。
//!
//! 実装予定 (RUST_PORT_PLAN.md Phase 1):
//! - S1.1 SMILES パーサ
//! - S1.2 分子グラフ構築 (原子価モデル・暗黙 H)
//! - S1.3 ケクレ化・芳香族認識
//! - S1.4 SSSR 環認識
//! - S1.5 正規 SMILES 生成・フラグメント分解
//! - S1.6 部分構造マッチ (VF2)
//! - S1.7 CIP 立体化学

mod aromaticity;
pub mod canon;
pub mod conformer;
pub mod depict;
#[cfg(feature = "editing")]
pub mod edit;
pub mod elements;
pub mod geometry;
pub mod graph;
pub mod inchi;
pub mod rings;
pub mod smarts;
pub mod smiles;
mod stereo;
pub mod substructure;

/// molrs 全体のエラー型。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ChemError {
    /// 不正な SMILES (パースエラー・原子価超過など)
    InvalidSmiles(String),
    /// 構造自体は正しいが実装の制限で扱えない (環認識の 128 原子上限など)。
    /// ライブラリがパニックする代わりにこれを返す (I31)。
    Unsupported(String),
}

impl std::fmt::Display for ChemError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ChemError::InvalidSmiles(s) => write!(f, "Invalid SMILES: {s}"),
            ChemError::Unsupported(s) => write!(f, "unsupported structure: {s}"),
        }
    }
}

impl std::error::Error for ChemError {}
