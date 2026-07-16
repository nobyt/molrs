//! 前処理 (隠し H) と鎖 (無環分子) のレイアウト (D2)。
//!
//! - 最長鎖 (木の直径パス) を水平ジグザグ (±30° 交互) に置く (GR-3.2)
//! - 隣接結合は 120° 分離、sp 中心 (三重結合・累積二重結合) は 180° (GR-4.1)
//! - 4 結合原子は 90°×4 を使用 (GR-4.1.3 の許容形)
//! - E/Z 二重結合は cip_ranks で幾何側を強制 (GR: 例外なく正しい幾何)
//!
//! すべて結合長 = 1.0、角度は 30° の倍数。決定的 (乱数不使用)。

use std::collections::VecDeque;
use std::f64::consts::PI;

use crate::conformer::params::{perceive_hybridization, Hybridization};
use crate::graph::MoleculeGraph;
use crate::stereo::cip_ranks;

use super::point2::{snap_angle_30, Point2};
use super::DepictError;

const DEG30: f64 = PI / 6.0;
const DEG60: f64 = PI / 3.0;
const DEG90: f64 = PI / 2.0;
const DEG120: f64 = 2.0 * PI / 3.0;

/// 描画時に隠す H (重原子に結合した明示 H ノード)。
/// 重原子隣接を持たない H ([H+]、[H][H] など) は隠さない。
pub(crate) fn hidden_h_flags(g: &MoleculeGraph) -> Vec<bool> {
    g.atoms
        .iter()
        .map(|a| {
            a.symbol == "H"
                && g.adjacency[a.idx]
                    .iter()
                    .any(|&nb| g.atoms[nb].symbol != "H")
        })
        .collect()
}

/// 可視原子のみの隣接リスト (順序は元の隣接順を保存)。
pub(crate) fn visible_adjacency(g: &MoleculeGraph, hidden: &[bool]) -> Vec<Vec<usize>> {
    (0..g.atoms.len())
        .map(|i| {
            if hidden[i] {
                Vec::new()
            } else {
                g.adjacency[i]
                    .iter()
                    .copied()
                    .filter(|&nb| !hidden[nb])
                    .collect()
            }
        })
        .collect()
}

/// 可視原子の重原子ごとの隠し H 数 (ラベル表示用)。
pub(crate) fn hidden_h_counts(g: &MoleculeGraph, hidden: &[bool]) -> Vec<u8> {
    (0..g.atoms.len())
        .map(|i| g.adjacency[i].iter().filter(|&&nb| hidden[nb]).count() as u8)
        .collect()
}

/// BFS で start から最遠の可視原子を返す (タイは最小インデックス)。
/// parent マップも返す。
fn bfs_farthest(start: usize, vadj: &[Vec<usize>], n: usize) -> (usize, Vec<Option<usize>>) {
    let mut dist = vec![usize::MAX; n];
    let mut parent = vec![None; n];
    let mut queue = VecDeque::new();
    dist[start] = 0;
    queue.push_back(start);
    let mut far = start;
    while let Some(v) = queue.pop_front() {
        if dist[v] > dist[far] {
            far = v;
        }
        for &nb in &vadj[v] {
            if dist[nb] == usize::MAX {
                dist[nb] = dist[v] + 1;
                parent[nb] = Some(v);
                queue.push_back(nb);
            }
        }
    }
    (far, parent)
}

/// sp 中心 (直線で描く原子) か。
fn is_linear(hyb: &[Hybridization], i: usize) -> bool {
    hyb[i] == Hybridization::Sp
}

/// v の未配置の子に方向を割り当てる。角度は v から子へ向かう向き。
pub(crate) fn assign_child_angles(
    placed_dirs: &[f64],
    n_children: usize,
    linear: bool,
) -> Vec<f64> {
    match placed_dirs.len() {
        0 => {
            // 孤立原子から開始 (直径パス端が孤立分岐の場合)。最初の子は +30°、
            // 以降は再帰で埋まる
            (0..n_children).map(|k| DEG30 + k as f64 * DEG120).collect()
        }
        1 => {
            let d = placed_dirs[0];
            let cont = d + PI; // 鎖の続行方向
            if linear && n_children == 1 {
                return vec![cont];
            }
            match n_children {
                1 => {
                    // ジグザグ: d±120° のうち水平に近い側 (GR-3.2: ±30° 最大化)
                    let c1 = d + DEG120;
                    let c2 = d - DEG120;
                    vec![pick_horizontal(c1, c2)]
                }
                2 => vec![d + DEG120, d - DEG120],
                3 => vec![cont + DEG90, cont, cont - DEG90], // 90° 十字 (GR-4.1.3)
                _ => fallback_slots(placed_dirs, n_children),
            }
        }
        2 => {
            // 既存 2 結合の逆二等分方向を基準に配置
            let u = Point2::from_angle(placed_dirs[0]) + Point2::from_angle(placed_dirs[1]);
            match u.normalized() {
                Some(v) => {
                    let b = snap_angle_30((-v).angle());
                    match n_children {
                        1 => vec![b],
                        2 => vec![b + DEG60, b - DEG60],
                        _ => fallback_slots(placed_dirs, n_children),
                    }
                }
                // 既存 2 結合が正反対 (直線貫通): 垂直に出す (90° 十字, GR-4.1.3)
                None => {
                    let b = snap_angle_30(placed_dirs[0] + DEG90);
                    match n_children {
                        1 => vec![b],
                        2 => vec![b, b + PI],
                        _ => fallback_slots(placed_dirs, n_children),
                    }
                }
            }
        }
        _ => fallback_slots(placed_dirs, n_children),
    }
}

/// 水平に近い方 (|sin| が小さい方) を選ぶ。タイは x 正方向 → 小さい角度。
fn pick_horizontal(c1: f64, c2: f64) -> f64 {
    let s1 = c1.sin().abs();
    let s2 = c2.sin().abs();
    if (s1 - s2).abs() > 1e-9 {
        return if s1 < s2 { c1 } else { c2 };
    }
    let x1 = c1.cos();
    let x2 = c2.cos();
    if (x1 - x2).abs() > 1e-9 {
        return if x1 > x2 { c1 } else { c2 };
    }
    if normalize_angle(c1) <= normalize_angle(c2) {
        c1
    } else {
        c2
    }
}

fn normalize_angle(a: f64) -> f64 {
    let mut a = a % (2.0 * PI);
    if a < 0.0 {
        a += 2.0 * PI;
    }
    a
}

/// 30° スロットから「既存方向との最小角距離が最大」のものを順に選ぶ。
/// タイは |sin| が小さい方 → 角度が小さい方。
pub(crate) fn fallback_slots(placed_dirs: &[f64], n_children: usize) -> Vec<f64> {
    let mut dirs: Vec<f64> = placed_dirs.to_vec();
    let mut out = Vec::with_capacity(n_children);
    for _ in 0..n_children {
        let mut best: Option<(f64, f64, f64, f64)> = None; // (-mindist, |sin|, angle, slot)
        for k in 0..12 {
            let slot = k as f64 * DEG30;
            let mind = dirs
                .iter()
                .map(|&d| {
                    let mut diff = (slot - d).abs() % (2.0 * PI);
                    if diff > PI {
                        diff = 2.0 * PI - diff;
                    }
                    diff
                })
                .fold(f64::INFINITY, f64::min);
            let key = (-mind, slot.sin().abs(), normalize_angle(slot), slot);
            if best.is_none()
                || (key.0, key.1, key.2) < (best.unwrap().0, best.unwrap().1, best.unwrap().2)
            {
                best = Some(key);
            }
        }
        let slot = best.unwrap().3;
        dirs.push(slot);
        out.push(slot);
    }
    out
}

/// 無環分子 (可視原子が木) のレイアウト。単一フラグメント前提。
pub(crate) fn layout_acyclic(
    g: &MoleculeGraph,
    hidden: &[bool],
    vadj: &[Vec<usize>],
) -> Result<Vec<Point2>, DepictError> {
    let n = g.atoms.len();
    let visible: Vec<usize> = (0..n).filter(|&i| !hidden[i]).collect();
    let mut pos = vec![Point2::ZERO; n];
    let mut placed = vec![false; n];
    if visible.is_empty() {
        return Err(DepictError::LayoutFailed("no visible atoms".into()));
    }

    let hyb = perceive_hybridization(g);
    // 直進させる原子: sp 中心、または 4 結合で全隣接が末端
    // (GR-4.1.3: 90° 十字が推奨されるケース)
    let straight_through = |i: usize| -> bool {
        is_linear(&hyb, i) || (vadj[i].len() == 4 && vadj[i].iter().all(|&c| vadj[c].len() == 1))
    };

    // 最長鎖 = 木の直径パス (可視部分グラフ上)
    let (u, _) = bfs_farthest(visible[0], vadj, n);
    let (w, parent) = bfs_farthest(u, vadj, n);
    let mut path = vec![w];
    while let Some(p) = parent[*path.last().unwrap()] {
        path.push(p);
    }
    path.reverse(); // u → w

    // 主鎖を水平ジグザグに配置。最初の結合は +30°、以降は 120° 交互。
    // sp 中心は直進。
    placed[path[0]] = true;
    let mut angle = DEG30;
    for i in 1..path.len() {
        let prev = path[i - 1];
        if i >= 2 && !straight_through(prev) {
            // 前の結合から ±120° のうち水平側 = ジグザグ交互
            angle = pick_horizontal(angle + DEG60, angle - DEG60);
        }
        pos[path[i]] = pos[prev] + Point2::from_angle(angle);
        placed[path[i]] = true;
    }

    // 分岐を BFS で配置 (主鎖順 → 子の順)
    let mut queue: VecDeque<usize> = path.iter().copied().collect();
    while let Some(v) = queue.pop_front() {
        let unplaced: Vec<usize> = vadj[v].iter().copied().filter(|&c| !placed[c]).collect();
        if unplaced.is_empty() {
            continue;
        }
        let placed_dirs: Vec<f64> = vadj[v]
            .iter()
            .copied()
            .filter(|&c| placed[c])
            .map(|c| (pos[c] - pos[v]).angle())
            .collect();
        let angles = assign_child_angles(&placed_dirs, unplaced.len(), is_linear(&hyb, v));
        for (&c, &a) in unplaced.iter().zip(angles.iter()) {
            pos[c] = pos[v] + Point2::from_angle(a);
            placed[c] = true;
            queue.push_back(c);
        }
    }

    if visible.iter().any(|&i| !placed[i]) {
        return Err(DepictError::LayoutFailed(
            "acyclic layout: disconnected visible atoms".into(),
        ));
    }
    Ok(pos)
}

/// E/Z の参照隣接 (CIP 最高位の可視隣接、結合相手を除く)。
fn ez_reference(
    g: &MoleculeGraph,
    hidden: &[bool],
    ranks: &[usize],
    atom: usize,
    other: usize,
) -> Option<usize> {
    g.adjacency[atom]
        .iter()
        .copied()
        .filter(|&nb| nb != other && !hidden[nb])
        .max_by_key(|&nb| (ranks[nb], usize::MAX - nb))
}

/// 2D 座標から二重結合の E/Z を再導出する。参照隣接が両端にあれば Some。
pub(crate) fn derive_ez(
    g: &MoleculeGraph,
    pos: &[Point2],
    hidden: &[bool],
    ranks: &[usize],
    bond_idx: usize,
) -> Option<char> {
    let b = &g.bonds[bond_idx];
    let (a1, a2) = (b.begin_idx, b.end_idx);
    let r1 = ez_reference(g, hidden, ranks, a1, a2)?;
    let r2 = ez_reference(g, hidden, ranks, a2, a1)?;
    let axis = pos[a2] - pos[a1];
    let s1 = axis.cross(pos[r1] - pos[a1]);
    let s2 = axis.cross(pos[r2] - pos[a2]);
    if s1.abs() < 1e-9 || s2.abs() < 1e-9 {
        return None; // 参照が軸上 (幾何が退化)
    }
    Some(if s1 * s2 > 0.0 { 'Z' } else { 'E' })
}

/// 木レイアウトの E/Z を強制する。違反時は end 側の部分木を結合軸で鏡映。
pub(crate) fn enforce_ez(
    g: &MoleculeGraph,
    pos: &mut [Point2],
    hidden: &[bool],
    vadj: &[Vec<usize>],
) {
    let ranks = cip_ranks(g);
    for (bi, b) in g.bonds.iter().enumerate() {
        let Some(want) = b.stereo else { continue };
        let Some(got) = derive_ez(g, pos, hidden, &ranks, bi) else {
            continue;
        };
        if got == want {
            continue;
        }
        // end 側部分木 (この結合自体を通らない到達集合) を軸で鏡映。
        // begin まで到達できる場合はこの結合が環内にある → 鏡映不能なので
        // スキップ (環内 E/Z は環レイアウトに従う)
        let (a1, a2) = (b.begin_idx, b.end_idx);
        let mut in_subtree = vec![false; g.atoms.len()];
        let mut stack = vec![a2];
        in_subtree[a2] = true;
        while let Some(v) = stack.pop() {
            for &nb in &vadj[v] {
                if (v == a2 && nb == a1) || (v == a1 && nb == a2) {
                    continue; // この結合の辺は通らない
                }
                if !in_subtree[nb] {
                    in_subtree[nb] = true;
                    stack.push(nb);
                }
            }
        }
        if in_subtree[a1] {
            continue; // 環内二重結合
        }
        let q = pos[a1];
        let Some(u) = (pos[a2] - pos[a1]).normalized() else {
            continue;
        };
        for (i, &inside) in in_subtree.iter().enumerate() {
            if inside && !hidden[i] {
                let v = pos[i] - q;
                pos[i] = q + u * (2.0 * v.dot(u)) - v;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    fn layout(smiles: &str) -> (MoleculeGraph, Vec<Point2>, Vec<bool>) {
        let g = build_molecule_graph(smiles).unwrap();
        let hidden = hidden_h_flags(&g);
        let vadj = visible_adjacency(&g, &hidden);
        let mut pos = layout_acyclic(&g, &hidden, &vadj).unwrap();
        enforce_ez(&g, &mut pos, &hidden, &vadj);
        (g, pos, hidden)
    }

    /// 可視結合が全て長さ 1、可視原子の結合間角度が全て 30° の倍数であること。
    fn assert_geometry(g: &MoleculeGraph, pos: &[Point2], hidden: &[bool]) {
        for b in &g.bonds {
            if hidden[b.begin_idx] || hidden[b.end_idx] {
                continue;
            }
            let d = pos[b.begin_idx].distance(pos[b.end_idx]);
            assert!((d - 1.0).abs() < 1e-6, "bond length {d}");
        }
        for i in 0..g.atoms.len() {
            if hidden[i] {
                continue;
            }
            let dirs: Vec<f64> = g.adjacency[i]
                .iter()
                .filter(|&&nb| !hidden[nb])
                .map(|&nb| (pos[nb] - pos[i]).angle())
                .collect();
            for a in &dirs {
                let steps = a / (PI / 6.0);
                assert!(
                    (steps - steps.round()).abs() < 1e-6,
                    "angle {} not multiple of 30°",
                    a.to_degrees()
                );
            }
            // 重なり (同方向) 禁止
            for (k, a) in dirs.iter().enumerate() {
                for b in &dirs[k + 1..] {
                    let mut diff = (a - b).abs() % (2.0 * PI);
                    if diff > PI {
                        diff = 2.0 * PI - diff;
                    }
                    assert!(diff > PI / 6.0 - 1e-6, "bonds overlap at atom {i}");
                }
            }
        }
    }

    #[test]
    fn butane_zigzag() {
        let (g, pos, hidden) = layout("CCCC");
        assert_geometry(&g, &pos, &hidden);
        // 主鎖はジグザグ: C0→C1 が +30°、C1→C2 が -30°、C2→C3 が +30°
        // (y 座標が 0, 0.5, 0, 0.5 ではなく 0, s, 0, s... の交互)
        let ys: Vec<f64> = (0..4).map(|i| pos[i].y).collect();
        assert!((ys[0] - ys[2]).abs() < 1e-9);
        assert!((ys[1] - ys[3]).abs() < 1e-9);
        assert!((ys[0] - ys[1]).abs() > 0.4);
        // 直線ではなく 120° 折れ
        let v1 = pos[1] - pos[0];
        let v2 = pos[2] - pos[1];
        assert!((v1.dot(v2) - 0.5).abs() < 1e-9); // cos60° = 0.5
    }

    #[test]
    fn isobutane_branch() {
        let (g, pos, hidden) = layout("CC(C)C");
        assert_geometry(&g, &pos, &hidden);
    }

    #[test]
    fn neopentane_cross() {
        let (g, pos, hidden) = layout("CC(C)(C)C");
        assert_geometry(&g, &pos, &hidden);
        // 中心炭素 (idx 1) の 4 結合が互いに 90° 以上
        let c = 1;
        let dirs: Vec<Point2> = g.adjacency[c]
            .iter()
            .filter(|&&nb| !hidden[nb])
            .map(|&nb| pos[nb] - pos[c])
            .collect();
        assert_eq!(dirs.len(), 4);
        for i in 0..4 {
            for j in i + 1..4 {
                assert!(dirs[i].dot(dirs[j]) < 1e-6, "angle < 90°");
            }
        }
    }

    #[test]
    fn allene_linear() {
        let (g, pos, hidden) = layout("C=C=C");
        assert_geometry(&g, &pos, &hidden);
        let v1 = pos[1] - pos[0];
        let v2 = pos[2] - pos[1];
        assert!((v1.dot(v2) - 1.0).abs() < 1e-9, "allene must be collinear");
    }

    #[test]
    fn butyne_linear() {
        let (g, pos, hidden) = layout("CC#CC");
        assert_geometry(&g, &pos, &hidden);
        // C1, C2 は sp: 3 結合 C0-C1-C2-C3 が全て一直線
        let v1 = (pos[1] - pos[0]).normalized().unwrap();
        let v2 = (pos[2] - pos[1]).normalized().unwrap();
        let v3 = (pos[3] - pos[2]).normalized().unwrap();
        assert!((v1.dot(v2) - 1.0).abs() < 1e-9);
        assert!((v2.dot(v3) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn ez_butene() {
        for (smi, want) in [("C/C=C/C", 'E'), ("C/C=C\\C", 'Z')] {
            let (g, pos, hidden) = layout(smi);
            assert_geometry(&g, &pos, &hidden);
            let ranks = cip_ranks(&g);
            let bi = g.bonds.iter().position(|b| b.bond_order == 2.0).unwrap();
            assert_eq!(g.bonds[bi].stereo, Some(want), "graph tag for {smi}");
            assert_eq!(
                derive_ez(&g, &pos, &hidden, &ranks, bi),
                Some(want),
                "layout geometry for {smi}"
            );
        }
    }

    #[test]
    fn ez_diene_both_preserved() {
        // (2E,4Z)-ヘキサジエン様: 複数 E/Z の同時保存
        let (g, pos, hidden) = layout("C/C=C/C=C\\C");
        assert_geometry(&g, &pos, &hidden);
        let ranks = cip_ranks(&g);
        for (bi, b) in g.bonds.iter().enumerate() {
            if let Some(want) = b.stereo {
                assert_eq!(derive_ez(&g, &pos, &hidden, &ranks, bi), Some(want));
            }
        }
    }

    #[test]
    fn water_single_visible_atom() {
        let (g, pos, hidden) = layout("O");
        let vis: Vec<usize> = (0..g.atoms.len()).filter(|&i| !hidden[i]).collect();
        assert_eq!(vis.len(), 1);
        assert_eq!(pos[vis[0]], Point2::ZERO);
    }

    #[test]
    fn hidden_h_semantics() {
        let g = build_molecule_graph("CO").unwrap();
        let hidden = hidden_h_flags(&g);
        let counts = hidden_h_counts(&g, &hidden);
        // C に 3H、O に 1H が畳まれる
        assert_eq!(counts[0], 3);
        assert_eq!(counts[1], 1);
        // 全 H が隠される
        let n_hidden = hidden.iter().filter(|&&h| h).count();
        assert_eq!(n_hidden, 4);
    }
}
