//! 距離境界行列 (RUST_3D_PLAN.md C3)。
//!
//! 全原子対 (付加 H 含む) の距離の下界・上界を組み立てる:
//!
//! - 1-2 (結合): 理想結合長 ± 0.01 Å
//! - 1-3 (結合角): 余弦定理。角は混成から、小員環 (4/5 員) と
//!   同一環内は環内角で補正
//! - 1-4 (トーション): cis 距離 〜 trans 距離。E/Z 指定された二重結合の
//!   まわりは該当側に固定
//! - 芳香環内: 平面正多角形の頂点間距離に固定 (平面性の主要因)
//! - デフォルト: 下界 = vdW 半径和 × 0.7、上界 = 大
//! - 三角不等式スムージング (Floyd–Warshall 型) 後、逆転した対は
//!   下界を上界にクランプ (縮合多環で vdW 下界が幾何的に成立しない場合)

use crate::conformer::params::{
    ideal_angle_for, ideal_bond_length, perceive_hybridization, vdw_radius, Hybridization,
};
use crate::graph::MoleculeGraph;
use crate::stereo::cip_ranks;

/// 距離境界行列 (対称)。
pub struct BoundsMatrix {
    pub n: usize,
    lower: Vec<f64>,
    upper: Vec<f64>,
}

impl BoundsMatrix {
    fn new(n: usize, default_upper: f64) -> Self {
        BoundsMatrix {
            n,
            lower: vec![0.0; n * n],
            upper: vec![default_upper; n * n],
        }
    }

    pub fn lower(&self, i: usize, j: usize) -> f64 {
        self.lower[i * self.n + j]
    }

    pub fn upper(&self, i: usize, j: usize) -> f64 {
        self.upper[i * self.n + j]
    }

    fn set(&mut self, i: usize, j: usize, lo: f64, up: f64) {
        let lo = lo.max(0.0);
        self.lower[i * self.n + j] = lo;
        self.lower[j * self.n + i] = lo;
        self.upper[i * self.n + j] = up;
        self.upper[j * self.n + i] = up;
    }

    /// 既存より狭い場合のみ更新する。
    fn tighten(&mut self, i: usize, j: usize, lo: f64, up: f64) {
        let cur_lo = self.lower(i, j);
        let cur_up = self.upper(i, j);
        self.set(i, j, cur_lo.max(lo), cur_up.min(up));
    }
}

/// スピロ等の既定角範囲 (±2°)。
fn return_default_angle(elem_angle: f64) -> (f64, f64) {
    (
        elem_angle - 2.0_f64.to_radians(),
        elem_angle + 2.0_f64.to_radians(),
    )
}

/// 余弦定理: 2 辺と挟角から対辺。
fn third_side(r1: f64, r2: f64, angle: f64) -> f64 {
    (r1 * r1 + r2 * r2 - 2.0 * r1 * r2 * angle.cos()).sqrt()
}

/// 1-4 距離の cis / trans 値。
/// r_ij, r_jk, r_kl: 結合長、theta_j = 角 i-j-k、theta_k = 角 j-k-l。
fn dist14(r_ij: f64, r_jk: f64, r_kl: f64, theta_j: f64, theta_k: f64) -> (f64, f64) {
    // j を原点、k を +x に置く。i は上半平面
    let i = (r_ij * theta_j.cos(), r_ij * theta_j.sin());
    let lx = r_jk - r_kl * theta_k.cos();
    let ly = r_kl * theta_k.sin();
    let cis = ((i.0 - lx).powi(2) + (i.1 - ly).powi(2)).sqrt();
    let trans = ((i.0 - lx).powi(2) + (i.1 + ly).powi(2)).sqrt();
    (cis, trans)
}

/// 与えられた辺長と目標内角に近い「閉じた平面多角形」の頂点座標を解く。
/// 内角を変数として、閉包誤差 |Σ s_i e^{iφ_i}|² + λ Σ(θ−target)² を
/// 勾配降下で最小化する (決定的・数十反復で十分)。
fn planar_ring_coords(sides: &[f64], targets: &[f64]) -> Vec<(f64, f64)> {
    let m = sides.len();
    // 内角の初期値: 目標角を (m−2)π に正規化
    let want: f64 = (m as f64 - 2.0) * std::f64::consts::PI;
    let sum_t: f64 = targets.iter().sum();
    let mut theta: Vec<f64> = targets
        .iter()
        .map(|t| t + (want - sum_t) / m as f64)
        .collect();

    let closure = |theta: &[f64]| -> (f64, f64, Vec<(f64, f64)>) {
        // 進行方向: 外角 = π − θ で曲がりながら辺を進む
        let mut phi = 0.0f64;
        let (mut x, mut y) = (0.0f64, 0.0f64);
        let mut pts = Vec::with_capacity(m);
        for i in 0..m {
            pts.push((x, y));
            x += sides[i] * phi.cos();
            y += sides[i] * phi.sin();
            phi += std::f64::consts::PI - theta[(i + 1) % m];
        }
        (x, y, pts)
    };

    // 数値勾配で閉包誤差を潰す (変数 m 個、コスト無視できる規模)
    let lambda = 0.02;
    for _ in 0..200 {
        let (cx, cy, _) = closure(&theta);
        let err = cx * cx + cy * cy;
        if err < 1e-10 {
            break;
        }
        let h = 1e-6;
        let mut grad = vec![0.0f64; m];
        for i in 0..m {
            let mut tp = theta.clone();
            tp[i] += h;
            let (px, py, _) = closure(&tp);
            let ep = px * px + py * py + lambda * (tp[i] - targets[i]).powi(2)
                - lambda * (theta[i] - targets[i]).powi(2);
            grad[i] = (ep - err) / h;
        }
        let gnorm: f64 = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
        if gnorm < 1e-12 {
            break;
        }
        let step = (0.1 / gnorm).min(0.05);
        for i in 0..m {
            theta[i] -= step * grad[i];
        }
    }
    closure(&theta).2
}

/// 境界行列を構築する。
pub fn build_bounds(g: &MoleculeGraph) -> BoundsMatrix {
    let n = g.atoms.len();
    let hyb = perceive_hybridization(g);

    // 上界のデフォルト: 全結合長の和 (分子の最大伸長)
    let total_len: f64 = g
        .bonds
        .iter()
        .map(|b| {
            ideal_bond_length(
                &g.atoms[b.begin_idx].symbol,
                &g.atoms[b.end_idx].symbol,
                b.bond_order,
            )
        })
        .sum::<f64>()
        .max(10.0);
    let mut bm = BoundsMatrix::new(n, total_len);

    // デフォルト下界: vdW 和 × 0.7 (i < j の非対角のみ)
    for i in 0..n {
        for j in (i + 1)..n {
            let lo = 0.7 * (vdw_radius(&g.atoms[i].symbol) + vdw_radius(&g.atoms[j].symbol));
            bm.set(i, j, lo, bm.upper(i, j));
        }
        // 対角は 0
        bm.lower[i * n + i] = 0.0;
        bm.upper[i * n + i] = 0.0;
    }

    // 結合長テーブル (再利用するのでキャッシュ)
    let bond_len = |i: usize, j: usize, order: f64| {
        ideal_bond_length(&g.atoms[i].symbol, &g.atoms[j].symbol, order)
    };
    let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n]; // (相手, order)
    for b in &g.bonds {
        adj[b.begin_idx].push((b.end_idx, b.bond_order));
        adj[b.end_idx].push((b.begin_idx, b.bond_order));
    }

    // ---- 1-2 ----
    for b in &g.bonds {
        let len = bond_len(b.begin_idx, b.end_idx, b.bond_order);
        bm.set(b.begin_idx, b.end_idx, len - 0.01, len + 0.01);
    }

    // ---- 1-3 ----
    // 方針: 実現可能性優先。角度が不確かな状況 (環内・小員環の環外置換基) は
    // 広い角度範囲で距離範囲に変換し、最終的な角の品質は UFF (C7) に任せる。
    // 中心 j の基準角は元素補正つき (S/P 系は四面体角より狭い)。
    let ring_interior = |i: usize, j: usize, k: usize| -> Option<f64> {
        let mut best: Option<usize> = None;
        for ring in &g.ring_atom_sets {
            if ring.contains(&i) && ring.contains(&j) && ring.contains(&k) {
                best = Some(best.map_or(ring.len(), |b: usize| b.min(ring.len())));
            }
        }
        best.map(|m| (m as f64 - 2.0) * std::f64::consts::PI / m as f64)
    };
    // j を含む最小環 (サイズと、原子が環に属するか判定するための環そのもの)
    let smallest_ring_of = |j: usize| -> Option<&Vec<usize>> {
        g.ring_atom_sets
            .iter()
            .filter(|r| r.contains(&j))
            .min_by_key(|r| r.len())
    };
    // 橋かけ環系の原子 (3 原子以上を共有する環ペア = ノルボルナン等)。
    // ここだけは環内角が大きく歪む (橋頭 93° 等) ため広い範囲を使う
    let mut bridged = vec![false; n];
    for (ri, ra) in g.ring_atom_sets.iter().enumerate() {
        for rb in g.ring_atom_sets.iter().skip(ri + 1) {
            let shared = ra.iter().filter(|a| rb.contains(a)).count();
            if shared >= 3 {
                for &a in ra.iter().chain(rb.iter()) {
                    bridged[a] = true;
                }
            }
        }
    }
    for j in 0..n {
        let nbrs = &adj[j];
        let elem_angle = ideal_angle_for(
            &g.atoms[j].symbol,
            g.atoms[j].is_aromatic,
            hyb[j],
            nbrs.len(),
        );
        for a in 0..nbrs.len() {
            for b in (a + 1)..nbrs.len() {
                let (i, o_ij) = nbrs[a];
                let (k, o_jk) = nbrs[b];
                let r1 = bond_len(i, j, o_ij);
                let r2 = bond_len(j, k, o_jk);
                let (lo_ang, hi_ang) = match ring_interior(i, j, k) {
                    Some(interior) => {
                        if bridged[j] {
                            // 橋かけ環: 橋頭 93° 等を実現可能に保つ広い範囲
                            (
                                interior.min(elem_angle) - 16.0_f64.to_radians(),
                                interior.max(elem_angle) + 5.0_f64.to_radians(),
                            )
                        } else {
                            // 通常の環: 後段の「環の平面多角形ピン」が正しい距離を
                            // 与えるので、ここは広めの整合範囲だけ置く
                            (
                                interior.min(elem_angle) - 8.0_f64.to_radians(),
                                interior.max(elem_angle) + 4.0_f64.to_radians(),
                            )
                        }
                    }
                    None => {
                        // 環外ペア: 小員環 (3/4 員) 中心では角が大きく歪む
                        match smallest_ring_of(j) {
                            Some(ring) if ring.len() <= 4 => {
                                let m = ring.len();
                                let in_i = ring.contains(&i);
                                let in_k = ring.contains(&k);
                                let deg: Option<f64> = if in_i != in_k {
                                    // 片方が環内: sp2 は平面の残り (360−interior)/2
                                    Some(if hyb[j] == Hybridization::Sp2 {
                                        if m == 3 {
                                            150.0
                                        } else {
                                            135.0
                                        }
                                    } else if m == 3 {
                                        118.0
                                    } else {
                                        113.0
                                    })
                                } else if !in_i && !in_k {
                                    // 両方環外 (CH2 の H-H など)
                                    Some(if m == 3 { 115.0 } else { 111.0 })
                                } else {
                                    // 両方環内だが同一環に j と乗らない (スピロ等) → 既定
                                    None
                                };
                                match deg {
                                    Some(d) => {
                                        let c = d.to_radians();
                                        (c - 5.0_f64.to_radians(), c + 5.0_f64.to_radians())
                                    }
                                    None => return_default_angle(elem_angle),
                                }
                            }
                            _ => (
                                elem_angle - 2.0_f64.to_radians(),
                                elem_angle + 2.0_f64.to_radians(),
                            ),
                        }
                    }
                };
                let d_lo = third_side(r1, r2, lo_ang) - 0.02;
                let d_hi = third_side(r1, r2, hi_ang) + 0.02;
                bm.tighten(i, k, d_lo, d_hi);
            }
        }
    }

    // ---- 芳香環内: 閉じた平面多角形の頂点間距離に固定 ----
    // 辺長は実結合長、内角は元素別目標角 (チオフェン S は 92° 等) を
    // 初期値として多角形が閉じるように解く。これで縮合ヘテロ芳香環でも
    // サンプリングが整合的になる (ピンを外すと乱数サンプルの質が落ちて
    // かえって埋め込み失敗が増える)
    for ring in &g.ring_atom_sets {
        let m = ring.len();
        if m < 5 || !ring.iter().all(|&a| g.atoms[a].is_aromatic) {
            continue;
        }
        let sides: Vec<f64> = (0..m)
            .map(|t| bond_len(ring[t], ring[(t + 1) % m], 1.5))
            .collect();
        // targets[i] = 頂点 ring[i] の内角 (辺 i-1 と辺 i の間 → 上の closure の
        // 定義では pts[i] の角は theta[i])
        let targets: Vec<f64> = ring
            .iter()
            .map(|&a| {
                ideal_angle_for(
                    &g.atoms[a].symbol,
                    g.atoms[a].is_aromatic,
                    hyb[a],
                    adj[a].len(),
                )
            })
            .collect();
        let pts = planar_ring_coords(&sides, &targets);
        for a in 0..m {
            for b in (a + 1)..m {
                let k = (b - a).min(m - (b - a));
                if k >= 2 {
                    let (dx, dy) = (pts[a].0 - pts[b].0, pts[a].1 - pts[b].1);
                    let d = (dx * dx + dy * dy).sqrt();
                    bm.tighten(ring[a], ring[b], d - 0.06, d + 0.06);
                }
            }
        }
    }

    // ---- 非芳香族環の 1-3 ピン (平面多角形ソルバ) ----
    // 環を元素別目標角で閉じた平面多角形として解き、隣々接 (k=2) の距離を
    // 固定する。パッカリング (椅子形等) は 1-3 距離をほぼ変えないため
    // タイトで安全。1-4 以遠は自由のまま (柔軟性を保つ)。
    // 橋かけ環は除外 (上の広い範囲に任せる)。
    for ring in &g.ring_atom_sets {
        let m = ring.len();
        if m < 4 || ring.iter().any(|&a| bridged[a]) || ring.iter().all(|&a| g.atoms[a].is_aromatic)
        {
            continue;
        }
        let sides: Vec<f64> = (0..m)
            .map(|t| {
                let (u, v) = (ring[t], ring[(t + 1) % m]);
                let order = adj[u]
                    .iter()
                    .find(|&&(x, _)| x == v)
                    .map(|&(_, o)| o)
                    .unwrap_or(1.0);
                bond_len(u, v, order)
            })
            .collect();
        let targets: Vec<f64> = ring
            .iter()
            .map(|&a| {
                ideal_angle_for(
                    &g.atoms[a].symbol,
                    g.atoms[a].is_aromatic,
                    hyb[a],
                    adj[a].len(),
                )
            })
            .collect();
        let pts = planar_ring_coords(&sides, &targets);
        for a in 0..m {
            let b = (a + 2) % m;
            if m == 4 && a >= 2 {
                break; // 4 員環の対角は 2 本だけ
            }
            let (dx, dy) = (pts[a].0 - pts[b].0, pts[a].1 - pts[b].1);
            let d = (dx * dx + dy * dy).sqrt();
            bm.tighten(ring[a], ring[b], d - 0.08, d + 0.08);
        }
    }

    // ---- 1-4 ----
    // E/Z 指定二重結合の高位置換基 (CIP ランク最高の隣接) を特定
    let n_kept = g.parser_to_graph.iter().flatten().count();
    let has_stereo_bonds = g.bonds.iter().any(|b| b.stereo.is_some());
    let ranks = if has_stereo_bonds {
        cip_ranks(g)
    } else {
        Vec::new()
    };

    for b in &g.bonds {
        let (j, k) = (b.begin_idx, b.end_idx);
        for &(i, o_ij) in &adj[j] {
            if i == k {
                continue;
            }
            for &(l, o_kl) in &adj[k] {
                if l == j || l == i {
                    continue;
                }
                // i-j-k-l の 1-4 対
                if adj[i].iter().any(|&(x, _)| x == l) {
                    continue; // i-l が結合している (小員環) → 1-2 が優先
                }
                // i と l が共通の隣接を持つ = 実質 1-3 (5 員環の橋頭対など)。
                // トーション由来の cis 下限を課すと環の実距離と矛盾する
                if adj[i]
                    .iter()
                    .any(|&(x, _)| adj[l].iter().any(|&(y, _)| x == y))
                {
                    continue;
                }
                let theta_j = ideal_angle_for(
                    &g.atoms[j].symbol,
                    g.atoms[j].is_aromatic,
                    hyb[j],
                    adj[j].len(),
                );
                let theta_k = ideal_angle_for(
                    &g.atoms[k].symbol,
                    g.atoms[k].is_aromatic,
                    hyb[k],
                    adj[k].len(),
                );
                let (cis, trans) = dist14(
                    bond_len(i, j, o_ij),
                    bond_len(j, k, b.bond_order),
                    bond_len(k, l, o_kl),
                    theta_j,
                    theta_k,
                );
                match b.stereo {
                    Some(ez) => {
                        // E/Z の基準は両側の CIP 最高位置換基
                        let hi_j = adj[j]
                            .iter()
                            .filter(|&&(x, _)| x != k && x < n_kept)
                            .max_by_key(|&&(x, _)| ranks[x])
                            .map(|&(x, _)| x);
                        let hi_k = adj[k]
                            .iter()
                            .filter(|&&(x, _)| x != j && x < n_kept)
                            .max_by_key(|&&(x, _)| ranks[x])
                            .map(|&(x, _)| x);
                        let (Some(hi_j), Some(hi_k)) = (hi_j, hi_k) else {
                            bm.tighten(i, l, cis - 0.1, trans + 0.1);
                            continue;
                        };
                        // (i,l) が高位対と同側か
                        let same_as_reference = (i == hi_j) == (l == hi_k);
                        let want_cis = match ez {
                            'Z' => same_as_reference,
                            _ => !same_as_reference, // 'E'
                        };
                        if want_cis {
                            bm.tighten(i, l, cis - 0.05, cis + 0.05);
                        } else {
                            bm.tighten(i, l, trans - 0.05, trans + 0.05);
                        }
                    }
                    None => {
                        bm.tighten(i, l, cis - 0.1, trans + 0.1);
                    }
                }
            }
        }
    }

    // ---- 三角不等式スムージング ----
    for k in 0..n {
        for i in 0..n {
            if i == k {
                continue;
            }
            for j in 0..n {
                if j == k || j == i {
                    continue;
                }
                let u = bm.upper(i, k) + bm.upper(k, j);
                if bm.upper(i, j) > u {
                    bm.upper[i * n + j] = u;
                }
                let l = (bm.lower(i, k) - bm.upper(k, j)).max(bm.lower(j, k) - bm.upper(k, i));
                if bm.lower(i, j) < l {
                    bm.lower[i * n + j] = l;
                }
            }
        }
    }
    // 対称性の復元 (スムージングは片側ずつ更新するため)
    for i in 0..n {
        for j in (i + 1)..n {
            let lo = bm.lower(i, j).max(bm.lower(j, i));
            let up = bm.upper(i, j).min(bm.upper(j, i));
            // 逆転した対 (縮合多環の vdW デフォルトなど) は下界を譲る
            let lo = lo.min(up);
            bm.set(i, j, lo, up);
        }
    }
    bm
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    fn bounds_for(smiles: &str) -> (MoleculeGraph, BoundsMatrix) {
        let g = build_molecule_graph(smiles).expect("valid");
        let bm = build_bounds(&g);
        (g, bm)
    }

    fn assert_valid(bm: &BoundsMatrix) {
        for i in 0..bm.n {
            for j in 0..bm.n {
                assert!(
                    bm.lower(i, j) <= bm.upper(i, j) + 1e-9,
                    "bounds inverted at ({i},{j}): {} > {}",
                    bm.lower(i, j),
                    bm.upper(i, j)
                );
            }
        }
    }

    #[test]
    fn ethane_bounds() {
        let (_, bm) = bounds_for("CC");
        assert_valid(&bm);
        // C-C 結合
        assert!((bm.lower(0, 1) - 1.525).abs() < 0.02);
        assert!((bm.upper(0, 1) - 1.545).abs() < 0.02);
        // ジェミナル H-H (原子 2,3 は C0 の H): 約 1.78 Å
        let d = third_side(1.089, 1.089, 109.471_f64.to_radians());
        assert!(bm.lower(2, 3) <= d && d <= bm.upper(2, 3));
    }

    #[test]
    fn benzene_bounds_are_planar_polygon() {
        let (_, bm) = bounds_for("c1ccccc1");
        assert_valid(&bm);
        // オルト = 結合長 1.394
        assert!(bm.lower(0, 1) > 1.37 && bm.upper(0, 1) < 1.42);
        // メタ (1-3) = √3 × 1.394 ≈ 2.414
        let meta = 3.0f64.sqrt() * 1.394;
        assert!(
            bm.lower(0, 2) <= meta && meta <= bm.upper(0, 2),
            "meta: [{}, {}] vs {meta}",
            bm.lower(0, 2),
            bm.upper(0, 2)
        );
        assert!(
            bm.upper(0, 2) - bm.lower(0, 2) < 0.2,
            "meta should be tight"
        );
        // パラ (1-4) = 2 × 1.394 ≈ 2.788 (芳香環固定)
        let para = 2.0 * 1.394;
        assert!(bm.lower(0, 3) <= para && para <= bm.upper(0, 3));
        assert!(
            bm.upper(0, 3) - bm.lower(0, 3) < 0.2,
            "para should be tight"
        );
    }

    #[test]
    fn butane_torsion_range() {
        let (_, bm) = bounds_for("CCCC");
        assert_valid(&bm);
        // C1-C4: cis ≈ 2.5 Å 〜 trans ≈ 3.9 Å の範囲
        assert!(bm.lower(0, 3) > 2.2 && bm.lower(0, 3) < 2.8);
        assert!(bm.upper(0, 3) > 3.7 && bm.upper(0, 3) < 4.2);
    }

    #[test]
    fn ez_pinning() {
        // trans-2-ブテン (E): 末端 C1-C4 は trans 距離に固定
        let (_, bm) = bounds_for("C/C=C/C");
        assert_valid(&bm);
        assert!(
            bm.lower(0, 3) > 3.5,
            "E-butene C1-C4 should be pinned trans: [{}, {}]",
            bm.lower(0, 3),
            bm.upper(0, 3)
        );
        // cis-2-ブテン (Z): C1-C4 は cis 距離に固定
        let (_, bm) = bounds_for("C/C=C\\C");
        assert_valid(&bm);
        assert!(
            bm.upper(0, 3) < 3.4,
            "Z-butene C1-C4 should be pinned cis: [{}, {}]",
            bm.lower(0, 3),
            bm.upper(0, 3)
        );
    }

    #[test]
    fn cyclohexane_and_fused_rings_stay_consistent() {
        for smi in [
            "C1CCCCC1",
            "C1CCC2CCCCC2C1", // デカリン
            "C1CC2CCC1CC2",   // ビシクロ[2.2.2]
            "c1ccc2ccccc2c1", // ナフタレン
            "C1CC1",          // シクロプロパン
            "C1CCC1",         // シクロブタン
        ] {
            let (_, bm) = bounds_for(smi);
            assert_valid(&bm);
        }
    }
}
