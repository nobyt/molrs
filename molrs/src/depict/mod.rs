//! 2D 構造式描画 (depict)。
//!
//! IUPAC Recommendations 2008 (構造式の図示、GR 規則) と 2006 (立体配置の
//! 図示) に従って 2D 座標を生成し、SVG / 2D MOL として出力する。
//! 計画: RUST_2D_PLAN.md
//!
//! パイプラインは 2 段で、既存の 3D (embed → optimize) と同型:
//! 1. [`compute_coords_2d`]: 結合長 = 1.0 の無次元レイアウト座標を生成
//! 2. [`to_svg`] / [`to_mol_block_2d`]: [`Style`] でスケールして描画
//!
//! スタイル ([`Style`]) は ChemDraw 語彙のパラメータ集合で、
//! IUPAC 既定 / ACS 1996 / Nature / RSC / Wiley のプリセットを持つ。

pub(crate) mod chain_layout;
pub mod point2;
pub mod style;

pub use point2::Point2;
pub use style::Style;

use crate::graph::MoleculeGraph;

/// レイアウトパラメータ。
#[derive(Debug, Clone)]
pub struct LayoutParams {
    /// 乱数シード (衝突解消のタイブレークのみに使用。同一入力 → 同一座標)
    pub seed: u64,
    /// 衝突解消の反復上限
    pub max_collision_iters: usize,
    /// 最終手段としての結合伸長を許可するか
    pub allow_bond_stretch: bool,
}

impl Default for LayoutParams {
    fn default() -> LayoutParams {
        LayoutParams {
            seed: 0xD2D2_2D2D,
            max_collision_iters: 40,
            allow_bond_stretch: true,
        }
    }
}

/// くさび結合の向き。narrow 端は常に `bond.begin_idx` 側。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WedgeDir {
    /// solid wedge (手前)
    Up,
    /// hashed wedge (奥)
    Down,
}

/// 2D レイアウト結果。
///
/// `pos` は全原子分 (隠し H 含む)。隠し H は親重原子の位置を持つ
/// (NaN を避け、MOL 出力等での事故を防ぐ)。
#[derive(Debug, Clone)]
pub struct Coords2D {
    /// 原子ごとのレイアウト座標 (結合長 = 1.0 単位)
    pub pos: Vec<Point2>,
    /// 描画時に隠す原子 (ラベルに畳まれた H)
    pub hidden: Vec<bool>,
    /// 結合ごとのくさび指定 (graph.bonds と同順)
    pub wedge: Vec<Option<WedgeDir>>,
}

/// 2D レイアウトのエラー。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DepictError {
    /// レイアウト不能 (詳細メッセージ付き)
    LayoutFailed(String),
    /// まだ移植されていない構造クラス (実装の進行に応じて縮小)
    Unsupported(String),
}

impl std::fmt::Display for DepictError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DepictError::LayoutFailed(s) => write!(f, "2D layout failed: {s}"),
            DepictError::Unsupported(s) => write!(f, "2D layout unsupported: {s}"),
        }
    }
}

impl std::error::Error for DepictError {}

/// 2D レイアウト座標を生成する (結合長 = 1.0 単位)。
///
/// 現在の対応範囲 (実装の進行に応じて拡大):
/// - 無環・単一フラグメント (D2)
/// - 環系は D4-D6、フラグメント並置は D7、くさびは D9 で対応
pub fn compute_coords_2d(
    g: &MoleculeGraph,
    _params: &LayoutParams,
) -> Result<Coords2D, DepictError> {
    let hidden = chain_layout::hidden_h_flags(g);
    if !g.ring_atom_sets.is_empty() {
        return Err(DepictError::Unsupported("ring systems (D4-D6)".into()));
    }
    let vadj = chain_layout::visible_adjacency(g, &hidden);

    // 単一フラグメント確認 (複数は D7)
    let visible: Vec<usize> = (0..g.atoms.len()).filter(|&i| !hidden[i]).collect();
    if let Some(&first) = visible.first() {
        let mut seen = vec![false; g.atoms.len()];
        let mut stack = vec![first];
        seen[first] = true;
        while let Some(v) = stack.pop() {
            for &nb in &vadj[v] {
                if !seen[nb] {
                    seen[nb] = true;
                    stack.push(nb);
                }
            }
        }
        if visible.iter().any(|&i| !seen[i]) {
            return Err(DepictError::Unsupported("multiple fragments (D7)".into()));
        }
    }

    let mut pos = chain_layout::layout_acyclic(g, &hidden, &vadj)?;
    chain_layout::enforce_ez(g, &mut pos, &hidden, &vadj);

    // 隠し H は親重原子の位置に置く (NaN 回避)
    for i in 0..g.atoms.len() {
        if hidden[i] {
            if let Some(&parent) = g.adjacency[i].iter().find(|&&nb| !hidden[nb]) {
                pos[i] = pos[parent];
            }
        }
    }

    let wedge = vec![None; g.bonds.len()];
    Ok(Coords2D { pos, hidden, wedge })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn layout_params_default_is_deterministic() {
        assert_eq!(LayoutParams::default().seed, LayoutParams::default().seed);
        assert!(LayoutParams::default().allow_bond_stretch);
    }

    #[test]
    fn error_display() {
        let e = DepictError::LayoutFailed("x".into());
        assert_eq!(e.to_string(), "2D layout failed: x");
    }
}
