//! 幾何・数値計算の基盤 (RUST_3D_PLAN.md C1)。
//!
//! 3D 配座生成 (距離幾何法 + UFF) が必要とする最小限の数値計算を
//! 依存クレートなしで提供する:
//!
//! - [`Vec3`]: 3 次元ベクトル
//! - [`jacobi_eigen`]: 対称行列の固有分解 (巡回 Jacobi 法)。
//!   距離幾何の計量行列埋め込み (C4) と Kabsch 重ね合わせが使う
//! - [`kabsch`] / [`kabsch_rmsd`]: 点群の最適重ね合わせ (Horn の四元数法)。
//!   **真回転のみ** (鏡映を含まない) なので、エナンチオマーは重ならない —
//!   立体保存の検証 (C6/C9) にはこの性質が必要
//! - [`SeededRng`]: 再現可能な乱数 (xorshift64*)。配座生成の決定性は
//!   API 契約 (同一シード → 同一座標) なので、std の HashMap 等に依存しない

use std::ops::{Add, Div, Mul, Neg, Sub};

// ---------------------------------------------------------------------------
// Vec3
// ---------------------------------------------------------------------------

/// 3 次元ベクトル。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const ZERO: Vec3 = Vec3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    /// 内積。
    pub fn dot(self, o: Vec3) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// 外積。
    pub fn cross(self, o: Vec3) -> Vec3 {
        Vec3 {
            x: self.y * o.z - self.z * o.y,
            y: self.z * o.x - self.x * o.z,
            z: self.x * o.y - self.y * o.x,
        }
    }

    pub fn norm_sq(self) -> f64 {
        self.dot(self)
    }

    pub fn norm(self) -> f64 {
        self.norm_sq().sqrt()
    }

    /// 単位ベクトル。ゼロベクトルは None。
    pub fn normalized(self) -> Option<Vec3> {
        let n = self.norm();
        if n < 1e-300 {
            None
        } else {
            Some(self / n)
        }
    }

    pub fn distance(self, o: Vec3) -> f64 {
        (self - o).norm()
    }
}

impl Add for Vec3 {
    type Output = Vec3;
    fn add(self, o: Vec3) -> Vec3 {
        Vec3::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}

impl Sub for Vec3 {
    type Output = Vec3;
    fn sub(self, o: Vec3) -> Vec3 {
        Vec3::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

impl Neg for Vec3 {
    type Output = Vec3;
    fn neg(self) -> Vec3 {
        Vec3::new(-self.x, -self.y, -self.z)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Vec3;
    fn mul(self, s: f64) -> Vec3 {
        Vec3::new(self.x * s, self.y * s, self.z * s)
    }
}

impl Div<f64> for Vec3 {
    type Output = Vec3;
    fn div(self, s: f64) -> Vec3 {
        Vec3::new(self.x / s, self.y / s, self.z / s)
    }
}

/// 3×3 行列 (行優先)。回転行列の表現に使う。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mat3(pub [[f64; 3]; 3]);

impl Mat3 {
    pub const IDENTITY: Mat3 = Mat3([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);

    /// ベクトルへの適用 (行列 × 列ベクトル)。
    pub fn apply(&self, v: Vec3) -> Vec3 {
        let m = &self.0;
        Vec3::new(
            m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z,
        )
    }

    pub fn det(&self) -> f64 {
        let m = &self.0;
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }
}

// ---------------------------------------------------------------------------
// Jacobi 固有分解
// ---------------------------------------------------------------------------

/// 対称行列 (row-major, n×n) の固有分解 (巡回 Jacobi 法)。
///
/// 返り値: `(固有値, 固有ベクトル)` — 固有値は**降順**、
/// `eigenvectors[k]` (長さ n) が `eigenvalues[k]` に対応する正規直交ベクトル。
///
/// 距離幾何 (C4) は上位 3 固有対だけを使うが、簡単のため全対を返す。
/// 分子サイズ (n ≤ 128) では計算量 O(n³/sweep) は問題にならない。
///
/// # Panics
/// `mat.len() != n * n` のとき。非対称行列を渡した場合の結果は未規定
/// (デバッグビルドでは対称性を検査する)。
pub fn jacobi_eigen(mat: &[f64], n: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    assert_eq!(mat.len(), n * n, "matrix size mismatch");
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    #[cfg(debug_assertions)]
    for p in 0..n {
        for q in (p + 1)..n {
            debug_assert!(
                (mat[p * n + q] - mat[q * n + p]).abs() < 1e-9,
                "jacobi_eigen expects a symmetric matrix"
            );
        }
    }

    let idx = |r: usize, c: usize| r * n + c;
    let mut a = mat.to_vec();
    // v は固有ベクトルを「列」に持つ (v[k*n + j] = j 番目の固有ベクトルの成分 k)
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[idx(i, i)] = 1.0;
    }

    const MAX_SWEEPS: usize = 100;
    for _ in 0..MAX_SWEEPS {
        // 非対角成分の二乗和で収束判定
        let mut off = 0.0;
        for p in 0..n {
            for q in (p + 1)..n {
                off += a[idx(p, q)] * a[idx(p, q)];
            }
        }
        if off < 1e-24 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[idx(p, q)];
                if apq.abs() < 1e-300 {
                    continue;
                }
                // 回転角の決定 (Numerical Recipes 流の安定な式)
                let theta = (a[idx(q, q)] - a[idx(p, p)]) / (2.0 * apq);
                // f64::signum(+0.0) == 1.0 なので theta == 0 でも t = 1 (45°) になる
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;

                // A ← Jᵀ A J (行と列の両方に回転を適用)
                for k in 0..n {
                    let akp = a[idx(k, p)];
                    let akq = a[idx(k, q)];
                    a[idx(k, p)] = c * akp - s * akq;
                    a[idx(k, q)] = s * akp + c * akq;
                }
                for k in 0..n {
                    let apk = a[idx(p, k)];
                    let aqk = a[idx(q, k)];
                    a[idx(p, k)] = c * apk - s * aqk;
                    a[idx(q, k)] = s * apk + c * aqk;
                }
                // V ← V J
                for k in 0..n {
                    let vkp = v[idx(k, p)];
                    let vkq = v[idx(k, q)];
                    v[idx(k, p)] = c * vkp - s * vkq;
                    v[idx(k, q)] = s * vkp + c * vkq;
                }
            }
        }
    }

    // 固有値 = 対角成分。降順に並べ替え
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| {
        a[idx(j, j)]
            .partial_cmp(&a[idx(i, i)])
            .expect("eigenvalues are finite")
    });
    let eigenvalues: Vec<f64> = order.iter().map(|&i| a[idx(i, i)]).collect();
    let eigenvectors: Vec<Vec<f64>> = order
        .iter()
        .map(|&j| (0..n).map(|k| v[idx(k, j)]).collect())
        .collect();
    (eigenvalues, eigenvectors)
}

// ---------------------------------------------------------------------------
// Kabsch 重ね合わせ (Horn の四元数法)
// ---------------------------------------------------------------------------

/// 重ね合わせの結果。`q ≈ rotation · (p − p_centroid) + q_centroid` の意味で
/// 点群 P を点群 Q に重ねる。
#[derive(Debug, Clone)]
pub struct Superposition {
    /// 真回転行列 (det = +1; 鏡映は含まない)
    pub rotation: Mat3,
    pub p_centroid: Vec3,
    pub q_centroid: Vec3,
    pub rmsd: f64,
}

impl Superposition {
    /// P 側の点を Q 系に写す。
    pub fn transform(&self, p: Vec3) -> Vec3 {
        self.rotation.apply(p - self.p_centroid) + self.q_centroid
    }
}

/// 点群 P を点群 Q に最適重ね合わせする (Horn, 1987 の四元数法)。
///
/// 対応点は同じインデックス同士。**真回転のみ**を許すため、
/// キラルな点群とその鏡像は RMSD が 0 にならない。
///
/// # Panics
/// 点数が一致しない、または点数が 0 のとき。
pub fn kabsch(p: &[Vec3], q: &[Vec3]) -> Superposition {
    assert_eq!(p.len(), q.len(), "point count mismatch");
    assert!(!p.is_empty(), "empty point sets");
    let n = p.len() as f64;

    let p_c = p.iter().fold(Vec3::ZERO, |a, &b| a + b) / n;
    let q_c = q.iter().fold(Vec3::ZERO, |a, &b| a + b) / n;

    // 共分散成分 S_ab = Σ (p−p̄)_a (q−q̄)_b
    let mut s = [[0.0f64; 3]; 3];
    for (&pi, &qi) in p.iter().zip(q) {
        let a = pi - p_c;
        let b = qi - q_c;
        let av = [a.x, a.y, a.z];
        let bv = [b.x, b.y, b.z];
        for (r, &ar) in av.iter().enumerate() {
            for (c, &bc) in bv.iter().enumerate() {
                s[r][c] += ar * bc;
            }
        }
    }
    let (sxx, sxy, sxz) = (s[0][0], s[0][1], s[0][2]);
    let (syx, syy, syz) = (s[1][0], s[1][1], s[1][2]);
    let (szx, szy, szz) = (s[2][0], s[2][1], s[2][2]);

    // Horn の 4×4 対称行列。最大固有値の固有ベクトルが最適四元数
    #[rustfmt::skip]
    let k = [
        sxx + syy + szz, syz - szy,        szx - sxz,        sxy - syx,
        syz - szy,       sxx - syy - szz,  sxy + syx,        szx + sxz,
        szx - sxz,       sxy + syx,        -sxx + syy - szz, syz + szy,
        sxy - syx,       szx + sxz,        syz + szy,        -sxx - syy + szz,
    ];
    let (_, vecs) = jacobi_eigen(&k, 4);
    let qv = &vecs[0]; // 最大固有値に対応 (降順ソート済み)
    let (w, x, y, z) = (qv[0], qv[1], qv[2], qv[3]);

    // 四元数 → 回転行列
    let rotation = Mat3([
        [
            w * w + x * x - y * y - z * z,
            2.0 * (x * y - w * z),
            2.0 * (x * z + w * y),
        ],
        [
            2.0 * (x * y + w * z),
            w * w - x * x + y * y - z * z,
            2.0 * (y * z - w * x),
        ],
        [
            2.0 * (x * z - w * y),
            2.0 * (y * z + w * x),
            w * w - x * x - y * y + z * z,
        ],
    ]);

    let mut sq_sum = 0.0;
    for (&pi, &qi) in p.iter().zip(q) {
        let moved = rotation.apply(pi - p_c) + q_c;
        sq_sum += (moved - qi).norm_sq();
    }
    Superposition {
        rotation,
        p_centroid: p_c,
        q_centroid: q_c,
        rmsd: (sq_sum / n).sqrt(),
    }
}

/// 最適重ね合わせ後の RMSD のみを返す。
pub fn kabsch_rmsd(p: &[Vec3], q: &[Vec3]) -> f64 {
    kabsch(p, q).rmsd
}

// ---------------------------------------------------------------------------
// シード付き乱数 (xorshift64*)
// ---------------------------------------------------------------------------

/// 再現可能な擬似乱数生成器 (xorshift64*)。
///
/// 配座生成の決定性 (同一シード → 同一座標) を API 契約にするための
/// 自前実装。暗号用途には使わないこと。
#[derive(Debug, Clone)]
pub struct SeededRng {
    state: u64,
}

impl SeededRng {
    /// 任意のシードから初期化する (0 も可 — splitmix64 で状態を作るため)。
    pub fn new(seed: u64) -> Self {
        // splitmix64 で 0 やビット偏りのあるシードをほぐす
        let mut z = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        SeededRng {
            state: if z == 0 { 0x9E37_79B9_7F4A_7C15 } else { z },
        }
    }

    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// [0, 1) の一様乱数。
    pub fn next_f64(&mut self) -> f64 {
        // 上位 53 ビットを仮数に
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// [lo, hi) の一様乱数。
    pub fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- Vec3 ----

    #[test]
    fn vec3_basics() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, -5.0, 6.0);
        assert_eq!(a + b, Vec3::new(5.0, -3.0, 9.0));
        assert_eq!(a - b, Vec3::new(-3.0, 7.0, -3.0));
        assert_eq!(a * 2.0, Vec3::new(2.0, 4.0, 6.0));
        assert!((a.dot(b) - 12.0).abs() < 1e-12); // 4 - 10 + 18
                                                  // 外積は両方に直交
        let c = a.cross(b);
        assert!(c.dot(a).abs() < 1e-12);
        assert!(c.dot(b).abs() < 1e-12);
        // 既知値: x × y = z
        assert_eq!(
            Vec3::new(1.0, 0.0, 0.0).cross(Vec3::new(0.0, 1.0, 0.0)),
            Vec3::new(0.0, 0.0, 1.0)
        );
        assert!((Vec3::new(3.0, 4.0, 0.0).norm() - 5.0).abs() < 1e-12);
        assert!(Vec3::ZERO.normalized().is_none());
        let u = a.normalized().unwrap();
        assert!((u.norm() - 1.0).abs() < 1e-12);
    }

    // ---- Jacobi ----

    #[test]
    fn jacobi_known_2x2() {
        // [[2,1],[1,2]] → 固有値 3, 1
        let (vals, vecs) = jacobi_eigen(&[2.0, 1.0, 1.0, 2.0], 2);
        assert!((vals[0] - 3.0).abs() < 1e-10);
        assert!((vals[1] - 1.0).abs() < 1e-10);
        // 固有値 3 のベクトルは (1,1)/√2 (符号は任意)
        let v = &vecs[0];
        assert!((v[0].abs() - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-10);
        assert!((v[0] - v[1]).abs() < 1e-10);
    }

    #[test]
    fn jacobi_diagonal_passthrough() {
        let (vals, _) = jacobi_eigen(&[5.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 7.0], 3);
        assert!((vals[0] - 7.0).abs() < 1e-12);
        assert!((vals[1] - 5.0).abs() < 1e-12);
        assert!((vals[2] + 2.0).abs() < 1e-12);
    }

    /// A v = λ v・正規直交性・再構成 A = Σ λ v vᵀ を、決定的に生成した
    /// ランダム対称行列で検査する。
    #[test]
    fn jacobi_random_symmetric_properties() {
        let mut rng = SeededRng::new(42);
        for &n in &[1usize, 2, 3, 5, 8, 16, 40] {
            // 対称行列を作る
            let mut m = vec![0.0; n * n];
            for i in 0..n {
                for j in i..n {
                    let x = rng.uniform(-10.0, 10.0);
                    m[i * n + j] = x;
                    m[j * n + i] = x;
                }
            }
            let (vals, vecs) = jacobi_eigen(&m, n);
            // 降順
            for k in 1..n {
                assert!(vals[k - 1] >= vals[k] - 1e-9);
            }
            for k in 0..n {
                // A v ≈ λ v
                for i in 0..n {
                    let av: f64 = (0..n).map(|j| m[i * n + j] * vecs[k][j]).sum();
                    assert!(
                        (av - vals[k] * vecs[k][i]).abs() < 1e-7,
                        "n={n} k={k}: residual too large"
                    );
                }
                // 正規直交
                for l in k..n {
                    let dot: f64 = (0..n).map(|j| vecs[k][j] * vecs[l][j]).sum();
                    let expect = if k == l { 1.0 } else { 0.0 };
                    assert!((dot - expect).abs() < 1e-9, "n={n}: orthonormality");
                }
            }
            // 再構成
            for i in 0..n {
                for j in 0..n {
                    let re: f64 = (0..n).map(|k| vals[k] * vecs[k][i] * vecs[k][j]).sum();
                    assert!((re - m[i * n + j]).abs() < 1e-7, "n={n}: reconstruction");
                }
            }
        }
    }

    // ---- Kabsch ----

    fn test_points() -> Vec<Vec3> {
        // 非平面・非対称な点群 (キラル)
        vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.5, 0.0, 0.0),
            Vec3::new(0.0, 2.0, 0.0),
            Vec3::new(0.0, 0.0, 3.0),
            Vec3::new(1.0, 1.0, 0.5),
        ]
    }

    #[test]
    fn kabsch_identity() {
        let p = test_points();
        assert!(kabsch_rmsd(&p, &p) < 1e-10);
    }

    #[test]
    fn kabsch_recovers_rotation_and_translation() {
        let p = test_points();
        // z 軸まわり 90° 回転 + 並進 (5, -3, 2)
        let q: Vec<Vec3> = p
            .iter()
            .map(|v| Vec3::new(-v.y + 5.0, v.x - 3.0, v.z + 2.0))
            .collect();
        let sup = kabsch(&p, &q);
        assert!(sup.rmsd < 1e-10, "rmsd = {}", sup.rmsd);
        // 真回転であること
        assert!((sup.rotation.det() - 1.0).abs() < 1e-9);
        // transform が実際に写像すること
        for (&pi, &qi) in p.iter().zip(&q) {
            assert!(sup.transform(pi).distance(qi) < 1e-9);
        }
    }

    #[test]
    fn kabsch_does_not_superpose_mirror_image() {
        // 鏡像 (x → −x) は真回転では重ならない (キラル検証の要)
        let p = test_points();
        let mirrored: Vec<Vec3> = p.iter().map(|v| Vec3::new(-v.x, v.y, v.z)).collect();
        let sup = kabsch(&p, &mirrored);
        assert!(
            sup.rmsd > 0.3,
            "mirror image should not superpose: rmsd = {}",
            sup.rmsd
        );
        assert!((sup.rotation.det() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn kabsch_noisy_pair_rmsd_value() {
        // 1 点だけ既知量ずらす: RMSD = δ/√N とはならず最適回転後の値になるが、
        // 少なくとも 0 < rmsd < δ の範囲に入る
        let p = test_points();
        let mut q = p.clone();
        q[4] = q[4] + Vec3::new(0.3, 0.0, 0.0);
        let r = kabsch_rmsd(&p, &q);
        assert!(r > 1e-3 && r < 0.3, "rmsd = {r}");
    }

    // ---- SeededRng ----

    #[test]
    fn rng_reproducible() {
        let mut a = SeededRng::new(12345);
        let mut b = SeededRng::new(12345);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
        // 異なるシードは異なる列 (先頭 8 個がすべて一致することはまずない)
        let mut c = SeededRng::new(12346);
        let mut a2 = SeededRng::new(12345);
        let same = (0..8).filter(|_| a2.next_u64() == c.next_u64()).count();
        assert!(same < 8);
        // シード 0 も動く
        let mut z = SeededRng::new(0);
        assert_ne!(z.next_u64(), 0);
    }

    #[test]
    fn rng_uniform_range_and_mean() {
        let mut rng = SeededRng::new(7);
        let n = 10_000;
        let mut sum = 0.0;
        for _ in 0..n {
            let x = rng.uniform(2.0, 4.0);
            assert!((2.0..4.0).contains(&x));
            sum += x;
        }
        let mean = sum / n as f64;
        assert!((mean - 3.0).abs() < 0.05, "mean = {mean}");
    }
}
