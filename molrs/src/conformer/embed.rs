//! 距離幾何の埋め込み (RUST_3D_PLAN.md C4)。
//!
//! 境界行列の範囲内で距離行列をランダムにサンプリングし、
//! 古典的 MDS (二重中心化した計量行列の固有分解) で 3D 初期座標を得る。
//! サンプリングされた距離集合は一般に 3 次元に埋め込み可能とは限らないので、
//! ここで得るのはあくまで初期値であり、C5 の誤差最小化で整える。

use crate::conformer::bounds::BoundsMatrix;
use crate::geometry::{jacobi_eigen, SeededRng, Vec3};

/// 境界内サンプリング + 計量行列埋め込みで初期座標を作る。
/// 最大固有値が正でない (縮退した) 場合は None。
pub(crate) fn embed_from_bounds(bm: &BoundsMatrix, rng: &mut SeededRng) -> Option<Vec<Vec3>> {
    let n = bm.n;
    if n == 0 {
        return Some(Vec::new());
    }
    if n == 1 {
        return Some(vec![Vec3::ZERO]);
    }

    // 距離の 2 乗行列をサンプリング
    let mut d2 = vec![0.0f64; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = rng.uniform(bm.lower(i, j), bm.upper(i, j));
            d2[i * n + j] = d * d;
            d2[j * n + i] = d * d;
        }
    }

    // 古典的 MDS: G = -1/2 J D² J (二重中心化)
    let mut row_mean = vec![0.0f64; n];
    let mut total = 0.0f64;
    for i in 0..n {
        let s: f64 = (0..n).map(|j| d2[i * n + j]).sum();
        row_mean[i] = s / n as f64;
        total += s;
    }
    let total_mean = total / (n * n) as f64;
    let mut gmat = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            gmat[i * n + j] = -0.5 * (d2[i * n + j] - row_mean[i] - row_mean[j] + total_mean);
        }
    }

    let (vals, vecs) = jacobi_eigen(&gmat, n);
    if vals[0] <= 1e-9 {
        return None; // 完全に縮退 (通常は起きない)
    }
    // 上位 3 固有対から座標。負の固有値は 0 扱い (3D に潰す)
    let s0 = vals[0].max(0.0).sqrt();
    let s1 = if n > 1 { vals[1].max(0.0).sqrt() } else { 0.0 };
    let s2 = if n > 2 { vals[2].max(0.0).sqrt() } else { 0.0 };
    let coords = (0..n)
        .map(|i| {
            Vec3::new(
                s0 * vecs[0][i],
                if n > 1 { s1 * vecs[1][i] } else { 0.0 },
                if n > 2 { s2 * vecs[2][i] } else { 0.0 },
            )
        })
        .collect();
    Some(coords)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformer::bounds::build_bounds;
    use crate::graph::build_molecule_graph;

    #[test]
    fn embedding_produces_finite_coords() {
        for smi in ["CC", "c1ccccc1", "CCCC", "C1CCCCC1"] {
            let g = build_molecule_graph(smi).unwrap();
            let bm = build_bounds(&g);
            let mut rng = SeededRng::new(1);
            let coords = embed_from_bounds(&bm, &mut rng).expect("embeds");
            assert_eq!(coords.len(), g.atoms.len());
            for c in &coords {
                assert!(c.x.is_finite() && c.y.is_finite() && c.z.is_finite());
            }
            // 原子が全部同一点に潰れていない
            let spread: f64 = coords
                .iter()
                .map(|c| c.distance(coords[0]))
                .fold(0.0, f64::max);
            assert!(spread > 0.5, "{smi}: spread = {spread}");
        }
    }

    #[test]
    fn embedding_is_deterministic() {
        let g = build_molecule_graph("CCO").unwrap();
        let bm = build_bounds(&g);
        let a = embed_from_bounds(&bm, &mut SeededRng::new(7)).unwrap();
        let b = embed_from_bounds(&bm, &mut SeededRng::new(7)).unwrap();
        assert_eq!(a.len(), b.len());
        for (p, q) in a.iter().zip(&b) {
            assert_eq!(p, q);
        }
    }
}
