//! 環系のレイアウト (D4-D6)。
//!
//! - 環系 = 辺 (2 原子以上) を共有する SSSR 環の連結成分。
//!   スピロ (1 原子共有) は別系として扱い place.rs が接続する
//! - 単環は正多角形 (GR-3.3: 五員環 108° 等は正多角形で自動充足)。
//!   大員環 (≥9) も正多角形 = 内角均等
//! - 縮合環は共有辺への正多角形貼り付け (ナフタレン等)
//! - 貼り付け不能な系 (橋かけ等) は「外周サイクル = 正多角形 + 内部原子 =
//!   反復重心配置 (Tutte 流)」で解く。橋結合の長さ 1.0 逸脱は許容
//!   (GR も橋かけの長い結合を許容する)
//!
//! すべて決定的。

use std::collections::VecDeque;
use std::f64::consts::PI;

use crate::graph::MoleculeGraph;

use super::point2::Point2;
use super::DepictError;

/// 辺共有で連結した環の集まり (ローカル座標を持つ)。
pub(crate) struct RingSystem {
    /// g.ring_atom_sets へのインデックス
    pub rings: Vec<usize>,
    /// 系に属する原子 (出現順、重複なし)
    pub atoms: Vec<usize>,
}

/// 2 つの環の共有原子数。
fn shared_atoms(r1: &[usize], r2: &[usize]) -> usize {
    r1.iter().filter(|a| r2.contains(a)).count()
}

/// SSSR 環を辺共有 (共有原子 ≥2) で連結成分にまとめる。
pub(crate) fn ring_systems(g: &MoleculeGraph) -> Vec<RingSystem> {
    let rings = &g.ring_atom_sets;
    let nr = rings.len();
    let mut comp = vec![usize::MAX; nr];
    let mut n_comp = 0;
    for start in 0..nr {
        if comp[start] != usize::MAX {
            continue;
        }
        let mut queue = VecDeque::from([start]);
        comp[start] = n_comp;
        while let Some(r) = queue.pop_front() {
            for r2 in 0..nr {
                if comp[r2] == usize::MAX && shared_atoms(&rings[r], &rings[r2]) >= 2 {
                    comp[r2] = n_comp;
                    queue.push_back(r2);
                }
            }
        }
        n_comp += 1;
    }
    (0..n_comp)
        .map(|c| {
            let members: Vec<usize> = (0..nr).filter(|&r| comp[r] == c).collect();
            let mut atoms = Vec::new();
            for &r in &members {
                for &a in &rings[r] {
                    if !atoms.contains(&a) {
                        atoms.push(a);
                    }
                }
            }
            RingSystem {
                rings: members,
                atoms,
            }
        })
        .collect()
}

/// 結合 (i, j) が系内のいずれかの環の辺か。
fn ring_edge_count(g: &MoleculeGraph, sys: &RingSystem, i: usize, j: usize) -> usize {
    sys.rings
        .iter()
        .filter(|&&r| {
            let ring = &g.ring_atom_sets[r];
            let n = ring.len();
            (0..n).any(|k| {
                let (a, b) = (ring[k], ring[(k + 1) % n]);
                (a == i && b == j) || (a == j && b == i)
            })
        })
        .count()
}

/// 系のローカル座標 (結合長 = 1.0)。返り値は (原子, 座標) のリスト。
pub(crate) fn layout_ring_system(
    g: &MoleculeGraph,
    sys: &RingSystem,
) -> Result<Vec<(usize, Point2)>, DepictError> {
    if let Some(res) = try_polygon_attach(g, sys) {
        return Ok(res);
    }
    perimeter_layout(g, sys)
}

/// 正多角形の頂点列を作る。最初の辺は p0→p1、以降 turn (±2π/n) ずつ曲がる。
fn walk_polygon(p0: Point2, p1: Point2, n: usize, turn: f64) -> Vec<Point2> {
    let mut pts = vec![p0, p1];
    let mut dir = p1 - p0;
    for _ in 2..n {
        dir = dir.rotated(turn);
        let next = *pts.last().unwrap() + dir;
        pts.push(next);
    }
    pts
}

/// 共有辺貼り付けによる縮合環レイアウト。橋かけ等で破綻したら None。
fn try_polygon_attach(g: &MoleculeGraph, sys: &RingSystem) -> Option<Vec<(usize, Point2)>> {
    let n_atoms_total = g.atoms.len();
    let mut pos = vec![Point2::ZERO; n_atoms_total];
    let mut placed = vec![false; n_atoms_total];

    // BFS 順 (辺共有隣接)。開始環は最小の環インデックス。
    let mut order: Vec<usize> = Vec::new();
    let mut visited = vec![false; sys.rings.len()];
    let mut queue = VecDeque::from([0usize]);
    visited[0] = true;
    while let Some(k) = queue.pop_front() {
        order.push(k);
        for (k2, vis) in visited.iter_mut().enumerate() {
            if !*vis
                && shared_atoms(
                    &g.ring_atom_sets[sys.rings[k]],
                    &g.ring_atom_sets[sys.rings[k2]],
                ) >= 2
            {
                *vis = true;
                queue.push_back(k2);
            }
        }
    }

    for (oi, &k) in order.iter().enumerate() {
        let ring = &g.ring_atom_sets[sys.rings[k]];
        let n = ring.len();
        let turn = 2.0 * PI / n as f64;
        if oi == 0 {
            // 最初の環: 正多角形。重心を原点、最初の頂点を上 (90°) に。
            let r = 0.5 / (PI / n as f64).sin();
            for (i, &a) in ring.iter().enumerate() {
                let theta = PI / 2.0 + i as f64 * 2.0 * PI / n as f64;
                pos[a] = Point2::from_angle(theta) * r;
                placed[a] = true;
            }
            continue;
        }
        // 配置済み原子との共有をチェック
        let placed_in_ring: Vec<usize> = ring.iter().copied().filter(|&a| placed[a]).collect();
        if placed_in_ring.len() != 2 {
            return None; // 橋かけ・多重共有 → 外周フォールバック
        }
        // 共有 2 原子が環順で隣接し、距離 ≈ 1 であること
        let idx = (0..n).find(|&i| placed[ring[i]] && placed[ring[(i + 1) % n]])?;
        let (a, b) = (ring[idx], ring[(idx + 1) % n]);
        if (pos[a].distance(pos[b]) - 1.0).abs() > 1e-6 {
            return None;
        }
        // 既配置系の重心から遠い側に貼る
        let centroid_placed = {
            let placed_atoms: Vec<Point2> = sys
                .atoms
                .iter()
                .filter(|&&x| placed[x])
                .map(|&x| pos[x])
                .collect();
            placed_atoms.iter().fold(Point2::ZERO, |s, &p| s + p) / placed_atoms.len() as f64
        };
        let cand1 = walk_polygon(pos[a], pos[b], n, turn);
        let cand2 = walk_polygon(pos[a], pos[b], n, -turn);
        let centroid_of = |pts: &[Point2]| -> Point2 {
            pts.iter().fold(Point2::ZERO, |s, &p| s + p) / pts.len() as f64
        };
        let pts = if centroid_of(&cand1).distance(centroid_placed)
            >= centroid_of(&cand2).distance(centroid_placed)
        {
            cand1
        } else {
            cand2
        };
        // 環順 idx, idx+1, ... に沿って割り当て。既配置原子と衝突したら破綻
        for (step, &p) in pts.iter().enumerate() {
            let atom = ring[(idx + step) % n];
            if placed[atom] {
                if pos[atom].distance(p) > 1e-6 {
                    return None;
                }
                continue;
            }
            // 既存原子との重なりチェック
            for &x in &sys.atoms {
                if placed[x] && pos[x].distance(p) < 0.3 {
                    return None;
                }
            }
            pos[atom] = p;
            placed[atom] = true;
        }
    }

    if sys.atoms.iter().any(|&a| !placed[a]) {
        return None;
    }
    Some(sys.atoms.iter().map(|&a| (a, pos[a])).collect())
}

/// 外周サイクル = 正多角形 + 内部原子 = 反復重心配置。
fn perimeter_layout(
    g: &MoleculeGraph,
    sys: &RingSystem,
) -> Result<Vec<(usize, Point2)>, DepictError> {
    // 系内の環結合と、その環所属数
    let mut edges: Vec<(usize, usize, usize)> = Vec::new(); // (i, j, count)
    for (ai, &i) in sys.atoms.iter().enumerate() {
        for &j in &sys.atoms[ai + 1..] {
            let cnt = ring_edge_count(g, sys, i, j);
            if cnt > 0 {
                edges.push((i, j, cnt));
            }
        }
    }
    // 外周候補 = 所属数 1 の結合
    let boundary: Vec<(usize, usize)> = edges
        .iter()
        .filter(|&&(_, _, c)| c == 1)
        .map(|&(i, j, _)| (i, j))
        .collect();

    let cycle = extract_single_cycle(&boundary);

    let (fixed_cycle, interior): (Vec<usize>, Vec<usize>) = match cycle {
        Some(cy) => {
            let interior: Vec<usize> = sys
                .atoms
                .iter()
                .copied()
                .filter(|a| !cy.contains(a))
                .collect();
            (cy, interior)
        }
        None => {
            // 最終フォールバック: 最大の環を固定し、残りを重心配置
            let &rmax = sys
                .rings
                .iter()
                .max_by_key(|&&r| (g.ring_atom_sets[r].len(), usize::MAX - r))
                .ok_or_else(|| DepictError::LayoutFailed("empty ring system".into()))?;
            let cy = g.ring_atom_sets[rmax].clone();
            let interior: Vec<usize> = sys
                .atoms
                .iter()
                .copied()
                .filter(|a| !cy.contains(a))
                .collect();
            (cy, interior)
        }
    };

    // 外周を正多角形に
    let m = fixed_cycle.len();
    let r = 0.5 / (PI / m as f64).sin();
    let n_total = g.atoms.len();
    let mut pos = vec![Point2::ZERO; n_total];
    for (i, &a) in fixed_cycle.iter().enumerate() {
        let theta = PI / 2.0 + i as f64 * 2.0 * PI / m as f64;
        pos[a] = Point2::from_angle(theta) * r;
    }

    // 内部原子: 系内隣接の平均を反復 (Tutte 流バリセントリック)
    // 初期値: わずかに散らした原点近傍 (対称性による退化の回避、決定的)
    for (k, &a) in interior.iter().enumerate() {
        let t = k as f64 * 2.399963; // 黄金角
        pos[a] = Point2::from_angle(t) * 0.01;
    }
    let in_sys = |x: usize| sys.atoms.contains(&x);
    for _ in 0..300 {
        for &a in &interior {
            let nbrs: Vec<usize> = g.adjacency[a]
                .iter()
                .copied()
                .filter(|&x| in_sys(x))
                .collect();
            if nbrs.is_empty() {
                continue;
            }
            let mean = nbrs.iter().fold(Point2::ZERO, |s, &x| s + pos[x]) / nbrs.len() as f64;
            pos[a] = mean;
        }
    }

    // 退化チェック (同一位置に潰れた原子)
    for (ai, &a) in sys.atoms.iter().enumerate() {
        for &b in &sys.atoms[ai + 1..] {
            if pos[a].distance(pos[b]) < 0.05 {
                return Err(DepictError::LayoutFailed(format!(
                    "ring system layout degenerate at atoms {a},{b}"
                )));
            }
        }
    }
    Ok(sys.atoms.iter().map(|&a| (a, pos[a])).collect())
}

/// 辺集合が単一サイクルを成すならその原子列を返す。
fn extract_single_cycle(edges: &[(usize, usize)]) -> Option<Vec<usize>> {
    if edges.is_empty() {
        return None;
    }
    let mut adj: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for &(i, j) in edges {
        adj.entry(i).or_default().push(j);
        adj.entry(j).or_default().push(i);
    }
    if adj.values().any(|v| v.len() != 2) {
        return None;
    }
    let start = *adj.keys().min()?;
    let mut cycle = vec![start];
    let mut prev = start;
    let mut cur = adj[&start][0];
    while cur != start {
        cycle.push(cur);
        let nbrs = &adj[&cur];
        let next = if nbrs[0] == prev { nbrs[1] } else { nbrs[0] };
        prev = cur;
        cur = next;
    }
    if cycle.len() == adj.len() {
        Some(cycle)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    fn layout(smiles: &str) -> (MoleculeGraph, Vec<Vec<(usize, Point2)>>) {
        let g = build_molecule_graph(smiles).unwrap();
        let systems = ring_systems(&g);
        let layouts = systems
            .iter()
            .map(|s| layout_ring_system(&g, s).unwrap())
            .collect();
        (g, layouts)
    }

    fn assert_ring_bonds_unit(g: &MoleculeGraph, coords: &[(usize, Point2)], tol: f64) {
        let pos_of = |a: usize| coords.iter().find(|(x, _)| *x == a).map(|(_, p)| *p);
        for ring in &g.ring_atom_sets {
            let n = ring.len();
            for k in 0..n {
                let (a, b) = (ring[k], ring[(k + 1) % n]);
                if let (Some(pa), Some(pb)) = (pos_of(a), pos_of(b)) {
                    let d = pa.distance(pb);
                    assert!((d - 1.0).abs() < tol, "ring bond {a}-{b} length {d}");
                }
            }
        }
    }

    #[test]
    fn benzene_regular_hexagon() {
        let (g, layouts) = layout("c1ccccc1");
        assert_eq!(layouts.len(), 1);
        assert_ring_bonds_unit(&g, &layouts[0], 1e-9);
        // 重心からの距離 = 外接円半径 1.0
        let c = layouts[0].iter().fold(Point2::ZERO, |s, &(_, p)| s + p) / 6.0;
        for &(_, p) in &layouts[0] {
            assert!((p.distance(c) - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn cyclopentane_regular() {
        let (g, layouts) = layout("C1CCCC1");
        assert_ring_bonds_unit(&g, &layouts[0], 1e-9);
    }

    #[test]
    fn macrocycle_regular() {
        let (g, layouts) = layout("C1CCCCCCCCCCC1"); // 12 員環
        assert_ring_bonds_unit(&g, &layouts[0], 1e-9);
    }

    #[test]
    fn naphthalene_fused() {
        let (g, layouts) = layout("c1ccc2ccccc2c1");
        assert_eq!(layouts.len(), 1);
        assert_eq!(layouts[0].len(), 10);
        assert_ring_bonds_unit(&g, &layouts[0], 1e-9);
        // 全原子ペアが十分離れている
        for (i, &(_, p)) in layouts[0].iter().enumerate() {
            for &(_, q) in &layouts[0][i + 1..] {
                assert!(p.distance(q) > 0.9);
            }
        }
    }

    #[test]
    fn anthracene_linear_fusion() {
        let (g, layouts) = layout("c1ccc2cc3ccccc3cc2c1");
        assert_eq!(layouts[0].len(), 14);
        assert_ring_bonds_unit(&g, &layouts[0], 1e-9);
    }

    #[test]
    fn indole_five_six() {
        let (g, layouts) = layout("c1ccc2[nH]ccc2c1");
        assert_eq!(layouts[0].len(), 9);
        assert_ring_bonds_unit(&g, &layouts[0], 1e-9);
    }

    #[test]
    fn norbornane_bridged() {
        // ビシクロ[2.2.1]ヘプタン: 貼り付け不能 → 外周 6 員環 + 橋原子
        let (_g, layouts) = layout("C1CC2CCC1C2");
        assert_eq!(layouts[0].len(), 7);
        // 外周結合は 1.0、橋結合は逸脱可。重なりなしのみ確認
        for (i, &(_, p)) in layouts[0].iter().enumerate() {
            for &(_, q) in &layouts[0][i + 1..] {
                assert!(p.distance(q) > 0.3, "atoms too close: {p:?} {q:?}");
            }
        }
    }

    #[test]
    fn adamantane_cage_fallback() {
        let (_, layouts) = layout("C1C2CC3CC1CC(C2)C3");
        assert_eq!(layouts[0].len(), 10);
        for (i, &(_, p)) in layouts[0].iter().enumerate() {
            for &(_, q) in &layouts[0][i + 1..] {
                assert!(p.distance(q) > 0.05);
            }
        }
    }

    #[test]
    fn spiro_gives_two_systems() {
        // スピロ[4.5]デカン: 1 原子共有 → 2 系
        let (_, layouts) = layout("C1CCC2(CC1)CCCC2");
        assert_eq!(layouts.len(), 2);
    }

    #[test]
    fn steroid_skeleton() {
        // ゴナン (ステロイド骨格): 6-6-6-5 縮合
        let (g, layouts) = layout("C1CCC2C(C1)CCC3C2CCC4C3CCC4");
        assert_eq!(layouts.len(), 1);
        assert_eq!(layouts[0].len(), 17);
        assert_ring_bonds_unit(&g, &layouts[0], 1e-6);
    }
}
