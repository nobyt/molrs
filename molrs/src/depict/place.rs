//! 分子全体の組み立て (D4-D6: 環系 + 鎖の統合)。
//!
//! 最大の環系をアンカーとして置き、BFS で鎖原子・次の環系 (結合接続 /
//! スピロ接続) を配置していく。環置換基の方向は既存結合の外角二等分
//! (GR-4.2)。E/Z は非環二重結合について部分木鏡映で強制する。

use std::collections::VecDeque;

use crate::conformer::params::{perceive_hybridization, Hybridization};
use crate::graph::MoleculeGraph;

use super::chain_layout::{assign_child_angles, enforce_ez, layout_acyclic};
use super::collide::resolve_collisions;
use super::point2::{snap_angle_30, Point2};
use super::ring_layout::{layout_ring_system, ring_systems, RingSystem};
use super::DepictError;

/// 分子全体のレイアウト (D7: 複数フラグメントは横並び配置)。
pub(crate) fn layout_molecule(
    g: &MoleculeGraph,
    hidden: &[bool],
    vadj: &[Vec<usize>],
    params: &super::LayoutParams,
) -> Result<Vec<Point2>, DepictError> {
    let n = g.atoms.len();
    let mut pos = vec![Point2::ZERO; n];

    // 可視原子の連結成分 (フラグメント) を列挙 (最小原子番号順)
    let mut frag_id = vec![usize::MAX; n];
    let mut n_frags = 0;
    for start in 0..n {
        if hidden[start] || frag_id[start] != usize::MAX {
            continue;
        }
        let mut stack = vec![start];
        frag_id[start] = n_frags;
        while let Some(v) = stack.pop() {
            for &nb in &vadj[v] {
                if frag_id[nb] == usize::MAX {
                    frag_id[nb] = n_frags;
                    stack.push(nb);
                }
            }
        }
        n_frags += 1;
    }

    let mut x_cursor = 0.0;
    for f in 0..n_frags {
        let frag: Vec<usize> = (0..n).filter(|&i| frag_id[i] == f).collect();
        // フラグメント外を隠した仮想グラフでレイアウト
        let hidden_f: Vec<bool> = (0..n).map(|i| hidden[i] || frag_id[i] != f).collect();
        let vadj_f = crate::depict::chain_layout::visible_adjacency(g, &hidden_f);
        let mut fpos = layout_fragment(g, &hidden_f, &vadj_f)?;
        enforce_ez(g, &mut fpos, &hidden_f, &vadj_f);
        resolve_collisions(g, &mut fpos, &hidden_f, &vadj_f, &frag, params);
        orient_fragment(g, &mut fpos, &frag, &hidden_f, &vadj_f);

        // 横並び: 左端を x_cursor に、y は重心 0 に揃える
        let (mut min_x, mut max_x, mut sum_y) = (f64::MAX, f64::MIN, 0.0);
        for &a in &frag {
            min_x = min_x.min(fpos[a].x);
            max_x = max_x.max(fpos[a].x);
            sum_y += fpos[a].y;
        }
        let cy = sum_y / frag.len() as f64;
        for &a in &frag {
            pos[a] = Point2::new(fpos[a].x - min_x + x_cursor, fpos[a].y - cy);
        }
        x_cursor += (max_x - min_x) + 1.5;
    }
    Ok(pos)
}

/// 単一フラグメントのレイアウト (環の有無で分岐)。
fn layout_fragment(
    g: &MoleculeGraph,
    hidden: &[bool],
    vadj: &[Vec<usize>],
) -> Result<Vec<Point2>, DepictError> {
    let has_ring_atom = (0..g.atoms.len()).any(|i| !hidden[i] && g.atoms[i].in_ring);
    if !has_ring_atom {
        layout_acyclic(g, hidden, vadj)
    } else {
        layout_with_rings(g, hidden, vadj)
    }
}

/// フラグメントの全体配向 (GR-3.2): 30° 刻み 12 回転 × 鏡映 2 の 24 候補から
/// 「水平 ±30° 帯の結合数」最大の配向を選ぶ。タイは横長 bbox → 小さい候補
/// インデックス。回転・鏡映とも E/Z の相対幾何を保存する。
fn orient_fragment(
    g: &MoleculeGraph,
    pos: &mut [Point2],
    frag: &[usize],
    hidden: &[bool],
    _vadj: &[Vec<usize>],
) {
    use std::f64::consts::PI;
    let centroid = frag.iter().fold(Point2::ZERO, |s, &a| s + pos[a]) / frag.len() as f64;
    let bonds: Vec<(usize, usize)> = g
        .bonds
        .iter()
        .filter(|b| !hidden[b.begin_idx] && !hidden[b.end_idx])
        .map(|b| (b.begin_idx, b.end_idx))
        .filter(|&(i, j)| frag.contains(&i) && frag.contains(&j))
        .collect();
    if bonds.is_empty() {
        return;
    }
    let mut best: Option<(i64, i64, usize)> = None; // (-score, -width_milli, cand_idx)
    let mut best_transform = (0.0, false);
    for cand in 0..24 {
        let k = cand % 12;
        let mirror = cand >= 12;
        let theta = k as f64 * PI / 6.0;
        let tf = |p: Point2| -> Point2 {
            let mut q = (p - centroid).rotated(theta);
            if mirror {
                q.y = -q.y;
            }
            q
        };
        let score = bonds
            .iter()
            .filter(|&&(i, j)| {
                let d = tf(pos[j]) - tf(pos[i]);
                let len = d.norm();
                len > 1e-9 && (d.y / len).abs() <= 0.5 + 1e-9 // |sin| ≤ sin30°
            })
            .count() as i64;
        let (mut min_x, mut max_x) = (f64::MAX, f64::MIN);
        for &a in frag {
            let q = tf(pos[a]);
            min_x = min_x.min(q.x);
            max_x = max_x.max(q.x);
        }
        let width_milli = ((max_x - min_x) * 1000.0) as i64;
        let key = (-score, -width_milli, cand);
        if best.is_none() || key < best.unwrap() {
            best = Some(key);
            best_transform = (theta, mirror);
        }
    }
    let (theta, mirror) = best_transform;
    for &a in frag {
        let mut q = (pos[a] - centroid).rotated(theta);
        if mirror {
            q.y = -q.y;
        }
        pos[a] = q + centroid;
    }
}

fn layout_with_rings(
    g: &MoleculeGraph,
    hidden: &[bool],
    vadj: &[Vec<usize>],
) -> Result<Vec<Point2>, DepictError> {
    let n = g.atoms.len();
    // アクティブ (可視) な原子を含む環系のみ対象 (フラグメントマスク対応)
    let mut systems = ring_systems(g);
    systems.retain(|s| s.atoms.iter().all(|&a| !hidden[a]));
    let layouts: Vec<Vec<(usize, Point2)>> = systems
        .iter()
        .map(|s| layout_ring_system(g, s))
        .collect::<Result<_, _>>()?;

    let hyb = perceive_hybridization(g);
    let mut pos = vec![Point2::ZERO; n];
    let mut placed = vec![false; n];
    let mut sys_placed = vec![false; systems.len()];

    // アンカー: 原子数最大の系 (タイは最小の環インデックス)
    let anchor = (0..systems.len())
        .max_by_key(|&i| (systems[i].atoms.len(), usize::MAX - i))
        .ok_or_else(|| DepictError::LayoutFailed("no ring systems".into()))?;
    for &(a, p) in &layouts[anchor] {
        pos[a] = p;
        placed[a] = true;
    }
    sys_placed[anchor] = true;

    let mut queue: VecDeque<usize> = systems[anchor].atoms.iter().copied().collect();
    while let Some(v) = queue.pop_front() {
        // 1. v をスピロ原子として共有する未配置の系
        for (si, sys) in systems.iter().enumerate() {
            if !sys_placed[si] && sys.atoms.contains(&v) {
                attach_spiro(g, &layouts[si], sys, v, &mut pos, &mut placed, vadj);
                sys_placed[si] = true;
                for &a in &sys.atoms {
                    queue.push_back(a);
                }
            }
        }
        // 2. 未配置の可視隣接
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
        let angles = assign_child_angles(&placed_dirs, unplaced.len(), hyb[v] == Hybridization::Sp);
        for (&c, &ang) in unplaced.iter().zip(angles.iter()) {
            if placed[c] {
                continue; // 環系接続で同時に置かれた場合
            }
            // c が未配置の環系に属するなら系ごと接続
            if let Some(si) =
                (0..systems.len()).find(|&si| !sys_placed[si] && systems[si].atoms.contains(&c))
            {
                attach_by_bond(&layouts[si], c, pos[v], ang, &mut pos, &mut placed);
                sys_placed[si] = true;
                for &a in &systems[si].atoms {
                    queue.push_back(a);
                }
            } else {
                pos[c] = pos[v] + Point2::from_angle(ang);
                placed[c] = true;
                queue.push_back(c);
            }
        }
        // v 自身を再度キューに (環系接続で新たな隣接が生じた場合に備える)
    }

    let visible: Vec<usize> = (0..n).filter(|&i| !hidden[i]).collect();
    if visible.iter().any(|&i| !placed[i]) {
        return Err(DepictError::LayoutFailed(
            "ring layout: unreached visible atoms".into(),
        ));
    }
    Ok(pos)
}

/// 環系 S を結合 (parent → entry) で接続する。entry は S の原子。
/// S の重心が結合方向の先に来るように回転する。
fn attach_by_bond(
    layout: &[(usize, Point2)],
    entry: usize,
    parent_pos: Point2,
    angle: f64,
    pos: &mut [Point2],
    placed: &mut [bool],
) {
    let p_entry = layout
        .iter()
        .find(|(a, _)| *a == entry)
        .map(|(_, p)| *p)
        .expect("entry atom in ring system layout");
    let centroid = layout.iter().fold(Point2::ZERO, |s, &(_, p)| s + p) / layout.len() as f64;
    let target = parent_pos + Point2::from_angle(angle);
    let rot = angle - (centroid - p_entry).angle();
    for &(a, p) in layout {
        if !placed[a] {
            pos[a] = (p - p_entry).rotated(rot) + target;
            placed[a] = true;
        }
    }
}

/// 環系 S をスピロ原子 v で接続する。系の重心が v の既存結合の
/// 反対側に来るように回転する。
fn attach_spiro(
    _g: &MoleculeGraph,
    layout: &[(usize, Point2)],
    _sys: &RingSystem,
    v: usize,
    pos: &mut [Point2],
    placed: &mut [bool],
    vadj: &[Vec<usize>],
) {
    let p_v = layout
        .iter()
        .find(|(a, _)| *a == v)
        .map(|(_, p)| *p)
        .expect("spiro atom in ring system layout");
    let centroid = layout.iter().fold(Point2::ZERO, |s, &(_, p)| s + p) / layout.len() as f64;
    // 既存結合方向の平均の反対 = 空いている方向
    let placed_dirs: Vec<Point2> = vadj[v]
        .iter()
        .copied()
        .filter(|&c| placed[c] && c != v)
        .map(|c| {
            (pos[c] - pos[v])
                .normalized()
                .unwrap_or(Point2::new(1.0, 0.0))
        })
        .collect();
    let free_dir = if placed_dirs.is_empty() {
        0.0
    } else {
        let sum = placed_dirs.iter().fold(Point2::ZERO, |s, &d| s + d);
        match sum.normalized() {
            Some(u) => snap_angle_30((-u).angle()),
            None => snap_angle_30((placed_dirs[0].perp()).angle()),
        }
    };
    let rot = free_dir - (centroid - p_v).angle();
    for &(a, p) in layout {
        if !placed[a] {
            pos[a] = (p - p_v).rotated(rot) + pos[v];
            placed[a] = true;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::depict::chain_layout::{hidden_h_flags, visible_adjacency};
    use crate::graph::build_molecule_graph;

    fn layout(smiles: &str) -> (MoleculeGraph, Vec<Point2>, Vec<bool>) {
        let g = build_molecule_graph(smiles).unwrap();
        let hidden = hidden_h_flags(&g);
        let vadj = visible_adjacency(&g, &hidden);
        let pos =
            layout_molecule(&g, &hidden, &vadj, &crate::depict::LayoutParams::default()).unwrap();
        (g, pos, hidden)
    }

    fn assert_no_overlap(g: &MoleculeGraph, pos: &[Point2], hidden: &[bool], min_d: f64) {
        let vis: Vec<usize> = (0..g.atoms.len()).filter(|&i| !hidden[i]).collect();
        for (k, &i) in vis.iter().enumerate() {
            for &j in &vis[k + 1..] {
                let d = pos[i].distance(pos[j]);
                assert!(d > min_d, "atoms {i},{j} too close: {d}");
            }
        }
    }

    /// 非環結合は長さ 1.0 (環内結合は ring_layout 側で担保)。
    fn assert_chain_bonds_unit(g: &MoleculeGraph, pos: &[Point2], hidden: &[bool]) {
        for (bi, b) in g.bonds.iter().enumerate() {
            if hidden[b.begin_idx] || hidden[b.end_idx] {
                continue;
            }
            let in_ring = g.ring_atom_sets.iter().any(|ring| {
                let n = ring.len();
                (0..n).any(|k| {
                    let (x, y) = (ring[k], ring[(k + 1) % n]);
                    (x == b.begin_idx && y == b.end_idx) || (x == b.end_idx && y == b.begin_idx)
                })
            });
            if in_ring {
                continue;
            }
            let d = pos[b.begin_idx].distance(pos[b.end_idx]);
            assert!((d - 1.0).abs() < 1e-6, "chain bond {bi} length {d}");
        }
    }

    #[test]
    fn toluene() {
        let (g, pos, hidden) = layout("Cc1ccccc1");
        assert_chain_bonds_unit(&g, &pos, &hidden);
        assert_no_overlap(&g, &pos, &hidden, 0.9);
        // メチル基は環外角二等分方向: 環重心から見て置換原子の延長上
        let ring_centroid = (1..7).fold(Point2::ZERO, |s, i| s + pos[i]) / 6.0;
        let ipso = pos[1];
        let methyl = pos[0];
        let out = (ipso - ring_centroid).normalized().unwrap();
        let sub = (methyl - ipso).normalized().unwrap();
        assert!(out.dot(sub) > 0.99, "substituent not on exterior bisector");
    }

    #[test]
    fn cyclohexanone_exocyclic_double() {
        let (g, pos, hidden) = layout("O=C1CCCCC1");
        assert_chain_bonds_unit(&g, &pos, &hidden);
        assert_no_overlap(&g, &pos, &hidden, 0.9);
    }

    #[test]
    fn biphenyl_two_systems() {
        let (g, pos, hidden) = layout("c1ccc(-c2ccccc2)cc1");
        assert_chain_bonds_unit(&g, &pos, &hidden);
        assert_no_overlap(&g, &pos, &hidden, 0.8);
    }

    #[test]
    fn spiro_decane() {
        let (g, pos, hidden) = layout("C1CCC2(CC1)CCCC2");
        assert_no_overlap(&g, &pos, &hidden, 0.5);
    }

    #[test]
    fn norbornane_no_overlap() {
        let (g, pos, hidden) = layout("C1CC2CCC1C2");
        assert_no_overlap(&g, &pos, &hidden, 0.3);
    }

    #[test]
    fn phenol_and_styrene() {
        let (g, pos, hidden) = layout("Oc1ccccc1");
        assert_chain_bonds_unit(&g, &pos, &hidden);
        assert_no_overlap(&g, &pos, &hidden, 0.9);
        let (g2, pos2, hidden2) = layout("C=Cc1ccccc1");
        assert_chain_bonds_unit(&g2, &pos2, &hidden2);
        assert_no_overlap(&g2, &pos2, &hidden2, 0.8);
    }

    #[test]
    fn salt_fragments_side_by_side() {
        let (g, pos, hidden) = layout("[Na+].[Cl-]");
        let vis: Vec<usize> = (0..g.atoms.len()).filter(|&i| !hidden[i]).collect();
        assert_eq!(vis.len(), 2);
        // 重ならず、横に並ぶ
        assert!(pos[vis[0]].distance(pos[vis[1]]) > 1.0);
        assert!((pos[vis[0]].y - pos[vis[1]].y).abs() < 1e-9);
    }

    #[test]
    fn orientation_maximizes_horizontal_band() {
        // 直鎖ヘキサン: 全結合が水平 ±30° 帯に入る
        let (g, pos, hidden) = layout("CCCCCC");
        for b in &g.bonds {
            if hidden[b.begin_idx] || hidden[b.end_idx] {
                continue;
            }
            let d = pos[b.end_idx] - pos[b.begin_idx];
            assert!(
                (d.y / d.norm()).abs() <= 0.5 + 1e-9,
                "bond outside horizontal band"
            );
        }
    }

    #[test]
    fn deterministic_two_runs() {
        let (_, p1, _) = layout("CC(=O)Oc1ccccc1C(=O)O.[Na+]");
        let (_, p2, _) = layout("CC(=O)Oc1ccccc1C(=O)O.[Na+]");
        for (a, b) in p1.iter().zip(p2.iter()) {
            assert_eq!(a.x, b.x);
            assert_eq!(a.y, b.y);
        }
    }

    #[test]
    fn ring_ez_not_broken() {
        // 環と E/Z 鎖の混在: E 配置が保存される
        let (g, pos, hidden) = layout("c1ccccc1/C=C/C");
        let ranks = crate::stereo::cip_ranks(&g);
        for (bi, b) in g.bonds.iter().enumerate() {
            if let Some(want) = b.stereo {
                let got = crate::depict::chain_layout::derive_ez(&g, &pos, &hidden, &ranks, bi);
                assert_eq!(got, Some(want));
            }
        }
    }
}
