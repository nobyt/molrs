//! 距離誤差の最小化 (RUST_3D_PLAN.md C5)。
//!
//! 誤差関数:
//! - 距離項: 全原子対で境界違反の 2 乗 (d > upper → (d−upper)²、d < lower → (lower−d)²)
//! - キラル体積項: 4 近傍の符号付き体積を目標範囲に (C6 が拘束を供給)
//! - 平面項: sp2 まわり 4 点の符号付き体積 → 0 (C6 が供給)
//!
//! 最小化は Polak–Ribière 共役勾配 + Armijo バックトラッキング直線探索。
//! 解析勾配を使う。

use crate::conformer::bounds::BoundsMatrix;
use crate::geometry::Vec3;

/// キラル体積拘束: 4 点の符号付き体積 (p1−p4)·((p2−p4)×(p3−p4)) を
/// [lower, upper] に収める。
#[derive(Debug, Clone)]
pub(crate) struct VolumeConstraint {
    pub atoms: [usize; 4],
    pub lower: f64,
    pub upper: f64,
    pub weight: f64,
}

/// 誤差関数の定義一式。
pub(crate) struct ErrorField<'a> {
    pub bounds: &'a BoundsMatrix,
    /// キラル体積 (符号固定) と平面 (lower=upper=0) の両方をこの形で持つ
    pub volumes: Vec<VolumeConstraint>,
}

/// 4 点の符号付き体積とその勾配。
pub(crate) fn signed_volume_of(p: &[Vec3], atoms: &[usize; 4]) -> (f64, [Vec3; 4]) {
    let (a, b, c, d) = (p[atoms[0]], p[atoms[1]], p[atoms[2]], p[atoms[3]]);
    let u1 = a - d;
    let u2 = b - d;
    let u3 = c - d;
    let v = u1.dot(u2.cross(u3));
    // ∂v/∂a = u2×u3, ∂v/∂b = u3×u1, ∂v/∂c = u1×u2, ∂v/∂d = −(和)
    let ga = u2.cross(u3);
    let gb = u3.cross(u1);
    let gc = u1.cross(u2);
    let gd = -(ga + gb + gc);
    (v, [ga, gb, gc, gd])
}

impl ErrorField<'_> {
    /// エネルギーと勾配 (座標の負方向が降下方向)。
    pub(crate) fn energy_and_grad(&self, coords: &[Vec3]) -> (f64, Vec<Vec3>) {
        let n = coords.len();
        let mut e = 0.0;
        let mut grad = vec![Vec3::ZERO; n];

        // 距離項
        for i in 0..n {
            for j in (i + 1)..n {
                let diff = coords[i] - coords[j];
                let d = diff.norm().max(1e-8);
                let lo = self.bounds.lower(i, j);
                let up = self.bounds.upper(i, j);
                let viol = if d > up {
                    d - up
                } else if d < lo {
                    d - lo // 負
                } else {
                    continue;
                };
                e += viol * viol;
                // ∂e/∂coords_i = 2 viol * (diff / d)
                let gi = diff * (2.0 * viol / d);
                grad[i] = grad[i] + gi;
                grad[j] = grad[j] - gi;
            }
        }

        // 体積項
        for vc in &self.volumes {
            let (v, gs) = signed_volume_of(coords, &vc.atoms);
            let viol = if v > vc.upper {
                v - vc.upper
            } else if v < vc.lower {
                v - vc.lower
            } else {
                continue;
            };
            e += vc.weight * viol * viol;
            let scale = 2.0 * vc.weight * viol;
            for (t, &ai) in vc.atoms.iter().enumerate() {
                grad[ai] = grad[ai] + gs[t] * scale;
            }
        }
        (e, grad)
    }
}

/// 距離幾何の誤差最小化に適した反復数の既定値。
pub(crate) fn default_iterations(n_atoms: usize) -> usize {
    (200 * n_atoms).max(1000)
}

/// Polak–Ribière 共役勾配 + Armijo バックトラッキング。
/// 収束後 (または反復上限で) の最終エネルギーを返す。
pub(crate) fn minimize(field: &ErrorField<'_>, coords: &mut [Vec3], max_iter: usize) -> f64 {
    minimize_with(&|c: &[Vec3]| field.energy_and_grad(c), coords, max_iter)
}

/// 任意のエネルギー関数に対する共役勾配最小化 (UFF (C7) と共用)。
pub(crate) fn minimize_with<F>(f: &F, coords: &mut [Vec3], max_iter: usize) -> f64
where
    F: Fn(&[Vec3]) -> (f64, Vec<Vec3>),
{
    let n = coords.len();
    if n == 0 {
        return 0.0;
    }
    let (mut e, mut g) = f(coords);
    let mut dir: Vec<Vec3> = g.iter().map(|&v| -v).collect();
    let mut g_norm_sq: f64 = g.iter().map(|v| v.norm_sq()).sum();

    for iter in 0..max_iter {
        if g_norm_sq < 1e-8 {
            break;
        }
        // 直線探索 (Armijo)
        let slope: f64 = g.iter().zip(&dir).map(|(gv, dv)| gv.dot(*dv)).sum();
        if slope >= 0.0 {
            // 降下方向でなくなったら最急降下にリセット
            for (d, gv) in dir.iter_mut().zip(&g) {
                *d = -*gv;
            }
            continue;
        }
        let mut alpha = 1.0;
        let mut trial: Vec<Vec3> = coords.to_vec();
        let mut e_new = e;
        let mut accepted = false;
        for _ in 0..40 {
            for i in 0..n {
                trial[i] = coords[i] + dir[i] * alpha;
            }
            let (et, _) = f(&trial);
            if et <= e + 1e-4 * alpha * slope {
                e_new = et;
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if !accepted {
            break; // これ以上下がらない
        }
        coords.copy_from_slice(&trial);

        let (e2, g_new) = f(coords);
        debug_assert!(e2 <= e + 1e-9, "energy must not increase");
        e = e_new.min(e2);

        // Polak–Ribière β
        let g_new_norm_sq: f64 = g_new.iter().map(|v| v.norm_sq()).sum();
        let mut num = 0.0;
        for (gn, go) in g_new.iter().zip(&g) {
            num += gn.dot(*gn - *go);
        }
        let beta = (num / g_norm_sq.max(1e-12)).max(0.0);
        for i in 0..n {
            dir[i] = -g_new[i] + dir[i] * beta;
        }
        // 周期的リスタート
        if (iter + 1) % (3 * n) == 0 {
            for (d, gv) in dir.iter_mut().zip(&g_new) {
                *d = -*gv;
            }
        }
        g = g_new;
        g_norm_sq = g_new_norm_sq;
    }
    e
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformer::bounds::build_bounds;
    use crate::graph::build_molecule_graph;

    /// 埋め込み + 距離最小化のスモーク。結合長が理想に近づくこと。
    #[test]
    fn minimization_restores_bond_lengths() {
        for smi in ["CC", "CCO", "c1ccccc1", "C1CCCCC1", "CC(C)C", "CCCCCC"] {
            let g = build_molecule_graph(smi).unwrap();
            let bm = build_bounds(&g);
            let iters = default_iterations(g.atoms.len());
            let (coords, maxv) =
                crate::conformer::embed_and_refine(&bm, &[], 11, 10, iters).expect("refine");
            assert!(maxv < crate::conformer::ACCEPT_MAX_VIOLATION);
            // 全結合が境界の ±0.1 Å 以内
            for b in &g.bonds {
                let d = coords[b.begin_idx].distance(coords[b.end_idx]);
                let lo = bm.lower(b.begin_idx, b.end_idx);
                let up = bm.upper(b.begin_idx, b.end_idx);
                assert!(
                    d > lo - 0.1 && d < up + 0.1,
                    "{smi}: bond ({},{}) = {d:.3}, bounds [{lo:.3},{up:.3}]",
                    b.begin_idx,
                    b.end_idx
                );
            }
        }
    }

    /// ベンゼンの粗い平面性。面外変位は距離に二次でしか効かないため、
    /// 境界のみでは強く拘束できない (厳密な平面性は C6 の平面項が担保する)。
    #[test]
    fn benzene_is_roughly_planar_from_bounds_alone() {
        let g = build_molecule_graph("c1ccccc1").unwrap();
        let bm = build_bounds(&g);
        let iters = default_iterations(g.atoms.len());
        let (coords, _) =
            crate::conformer::embed_and_refine(&bm, &[], 3, 10, iters).expect("refine");
        // 環炭素 6 つの最小二乗平面からの RMS
        let ring: Vec<Vec3> = coords[..6].to_vec();
        let centroid = ring.iter().fold(Vec3::ZERO, |a, &b| a + b) / 6.0;
        // 慣性テンソルの最小固有ベクトル = 法線
        let mut m = [0.0f64; 9];
        for p in &ring {
            let r = *p - centroid;
            let v = [r.x, r.y, r.z];
            for a in 0..3 {
                for b in 0..3 {
                    m[a * 3 + b] += v[a] * v[b];
                }
            }
        }
        let (_, vecs) = crate::geometry::jacobi_eigen(&m, 3);
        let normal = Vec3::new(vecs[2][0], vecs[2][1], vecs[2][2]);
        let rms = (ring
            .iter()
            .map(|p| (*p - centroid).dot(normal).powi(2))
            .sum::<f64>()
            / 6.0)
            .sqrt();
        assert!(rms < 0.3, "benzene planarity rms = {rms}");
    }

    /// E/Z 境界ピンが座標に反映される。
    #[test]
    fn ez_geometry_realized() {
        // E-2-ブテン: 末端炭素は遠い
        let g = build_molecule_graph("C/C=C/C").unwrap();
        let bm = build_bounds(&g);
        let iters = default_iterations(g.atoms.len());
        let (coords, _) =
            crate::conformer::embed_and_refine(&bm, &[], 5, 10, iters).expect("refine");
        let d_e = coords[0].distance(coords[3]);
        assert!(d_e > 3.5, "E-butene C1-C4 = {d_e}");

        // Z-2-ブテン: 近い
        let g = build_molecule_graph("C/C=C\\C").unwrap();
        let bm = build_bounds(&g);
        let (coords, _) =
            crate::conformer::embed_and_refine(&bm, &[], 5, 10, iters).expect("refine");
        let d_z = coords[0].distance(coords[3]);
        assert!(d_z < 3.3, "Z-butene C1-C4 = {d_z}");
    }
}
