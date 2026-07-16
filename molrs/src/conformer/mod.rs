//! 3D 配座生成 (RUST_3D_PLAN.md)。距離幾何法 + 知識項。
//!
//! - C2: [`params`] 理想幾何パラメータ
//! - C3: [`bounds`] 境界行列
//! - C4/C5: 埋め込みと誤差最小化 (リトライつき)
//! - C6: [`stereo3d`] キラル体積・平面知識項と 3D 立体検証
//! - C7: [`uff`] UFF 力場による仕上げ最小化
//! - C8: [`embed_molecule`] 公開 API、[`molblock`] V2000 出力

pub mod bounds;
pub(crate) mod embed;
pub(crate) mod exp_torsions;
pub(crate) mod minimize;
pub mod molblock;
pub mod params;
pub(crate) mod stereo3d;
pub(crate) mod torsion_lib;
pub(crate) mod uff;

use crate::geometry::{SeededRng, Vec3};

/// 距離違反の最大値がこの値未満なら埋め込み成功とみなす (Å)。
pub(crate) const ACCEPT_MAX_VIOLATION: f64 = 0.08;

/// 埋め込み + 最小化をリトライつきで行う内部ヘルパ。
/// 成功時は (座標, 最大距離違反)。試行ごとにシードを変えて再サンプリングする
/// (距離幾何は初期距離行列しだいで局所解に捕まるため、リトライが本質的に必要)。
pub(crate) fn embed_and_refine(
    bm: &bounds::BoundsMatrix,
    volumes: &[minimize::VolumeConstraint],
    seed: u64,
    max_attempts: u32,
    max_iter: usize,
) -> Option<(Vec<Vec3>, f64)> {
    let field = minimize::ErrorField {
        bounds: bm,
        volumes: volumes.to_vec(),
    };
    for attempt in 0..max_attempts {
        let mut rng = SeededRng::new(
            seed.wrapping_add(attempt as u64)
                .wrapping_mul(0x9E3779B97F4A7C15)
                .wrapping_add(attempt as u64),
        );
        let Some(mut coords) = embed::embed_from_bounds(bm, &mut rng) else {
            continue;
        };
        // 鏡映レスキュー: キラル体積 (符号固定) の違反が鏡映で減るなら反転してから最小化
        let chiral_violation = |cs: &[Vec3]| -> f64 {
            volumes
                .iter()
                .filter(|vc| vc.lower > 0.0 || vc.upper < 0.0)
                .map(|vc| {
                    let (v, _) = minimize::signed_volume_of(cs, &vc.atoms);
                    (vc.lower - v).max(v - vc.upper).max(0.0)
                })
                .sum()
        };
        if chiral_violation(&coords) > 0.0 {
            let mirrored: Vec<Vec3> = coords.iter().map(|c| Vec3::new(-c.x, c.y, c.z)).collect();
            if chiral_violation(&mirrored) < chiral_violation(&coords) {
                coords = mirrored;
            }
        }
        minimize::minimize(&field, &mut coords, max_iter);
        // 最大距離違反
        let mut maxv: f64 = 0.0;
        for i in 0..coords.len() {
            for j in (i + 1)..coords.len() {
                let d = coords[i].distance(coords[j]);
                let v = (d - bm.upper(i, j)).max(bm.lower(i, j) - d).max(0.0);
                maxv = maxv.max(v);
            }
        }
        // 体積拘束の違反
        let ok_volumes = volumes.iter().all(|vc| {
            let (v, _) = minimize::signed_volume_of(&coords, &vc.atoms);
            v >= vc.lower - 0.15 && v <= vc.upper + 0.15
        });
        if maxv < ACCEPT_MAX_VIOLATION && ok_volumes {
            return Some((coords, maxv));
        }
    }
    None
}

/// 配座生成のパラメータ。
#[derive(Debug, Clone)]
pub struct EmbedParams {
    /// 乱数シード (同一シード → 同一座標が API 契約)
    pub seed: u64,
    /// 距離行列の再サンプリング試行回数の上限
    pub max_attempts: u32,
    /// UFF 力場による仕上げ最小化 (C7) を行うか
    pub optimize: bool,
    /// ETKDG 実験トーション選好 (C10) を UFF 段に適用するか
    pub use_exp_torsions: bool,
}

impl Default for EmbedParams {
    fn default() -> Self {
        EmbedParams {
            seed: 0xC0FFEE,
            max_attempts: 30,
            optimize: true,
            use_exp_torsions: true,
        }
    }
}

/// 生成された配座 (全原子分の座標、付加 H を含む)。
#[derive(Debug, Clone)]
pub struct Conformer {
    pub coords: Vec<Vec3>,
}

/// 配座生成の失敗。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EmbedError {
    /// 試行回数内に受理可能な座標が得られなかった。
    /// 再現用にシードを含む。
    Failed { seed: u64, attempts: u32 },
}

impl std::fmt::Display for EmbedError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EmbedError::Failed { seed, attempts } => {
                write!(
                    f,
                    "3D embedding failed after {attempts} attempts (seed {seed})"
                )
            }
        }
    }
}

impl std::error::Error for EmbedError {}

/// SMILES 由来の分子グラフから 3D 座標を生成する (距離幾何 + 知識項)。
///
/// RDKit の `EmbedMolecule` (KDG 相当) に対応する。UFF 仕上げ (C7) は未実装。
pub fn embed_molecule(
    g: &crate::graph::MoleculeGraph,
    params: &EmbedParams,
) -> Result<Conformer, EmbedError> {
    let bm = bounds::build_bounds(g);
    let volumes = stereo3d::build_volume_constraints(g);
    let iters = minimize::default_iterations(g.atoms.len());
    match embed_and_refine(&bm, &volumes, params.seed, params.max_attempts, iters) {
        Some((mut coords, _)) => {
            // UFF 仕上げ (C7): 最適化後に有限性と立体保存を検証し、
            // 問題があれば距離幾何の座標のまま返す
            if params.optimize {
                if let Some(field) = uff::build_uff(g, params.use_exp_torsions) {
                    let mut opt = coords.clone();
                    uff::optimize(&field, &mut opt, iters);
                    let finite = opt
                        .iter()
                        .all(|c| c.x.is_finite() && c.y.is_finite() && c.z.is_finite());
                    let stereo_ok = stereo3d::verify_atom_stereo(g, &opt)
                        .iter()
                        .all(|&(idx, code)| g.atoms[idx].chiral_tag == Some(code))
                        && stereo3d::verify_bond_stereo(g, &opt)
                            .iter()
                            .all(|&(ei, code)| g.bonds[ei].stereo == Some(code));
                    if finite && stereo_ok {
                        coords = opt;
                    }
                }
            }
            Ok(Conformer { coords })
        }
        None => Err(EmbedError::Failed {
            seed: params.seed,
            attempts: params.max_attempts,
        }),
    }
}

/// 実験トーションの収集結果 (検証ゲート用)。
pub fn exp_torsions_of(g: &crate::graph::MoleculeGraph) -> Vec<exp_torsions::ExpTorsion> {
    exp_torsions::collect_exp_torsions(g)
}

/// 各結合の UFF 平衡長 (検証ゲート用)。型付け不能な分子は None。
pub fn uff_bond_rest_lengths(g: &crate::graph::MoleculeGraph) -> Option<Vec<f64>> {
    uff::bond_rest_lengths(g)
}

/// 3D 座標が入力の立体指定 (R/S, E/Z) を全て再現しているか検証する。
pub fn verify_stereo_3d(g: &crate::graph::MoleculeGraph, conf: &Conformer) -> bool {
    stereo3d::verify_atom_stereo(g, &conf.coords)
        .iter()
        .all(|&(idx, code)| g.atoms[idx].chiral_tag == Some(code))
        && stereo3d::verify_bond_stereo(g, &conf.coords)
            .iter()
            .all(|&(ei, code)| g.bonds[ei].stereo == Some(code))
}
