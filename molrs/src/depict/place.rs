//! 分子全体の組み立て (D4-D6: 環系 + 鎖の統合)。
//!
//! 最大の環系をアンカーとして置き、BFS で鎖原子・次の環系 (結合接続 /
//! スピロ接続) を配置していく。環置換基の方向は既存結合の外角二等分
//! (GR-4.2)。E/Z は非環二重結合について部分木鏡映で強制する。

use std::collections::VecDeque;

use crate::conformer::params::{perceive_hybridization, Hybridization};
use crate::graph::MoleculeGraph;

use super::chain_layout::{assign_child_angles, enforce_ez, layout_acyclic};
use super::point2::{snap_angle_30, Point2};
use super::ring_layout::{layout_ring_system, ring_systems, RingSystem};
use super::DepictError;

/// 単一フラグメントの分子全体レイアウト。
pub(crate) fn layout_molecule(
    g: &MoleculeGraph,
    hidden: &[bool],
    vadj: &[Vec<usize>],
) -> Result<Vec<Point2>, DepictError> {
    let mut pos = if g.ring_atom_sets.is_empty() {
        layout_acyclic(g, hidden, vadj)?
    } else {
        layout_with_rings(g, hidden, vadj)?
    };
    enforce_ez(g, &mut pos, hidden, vadj);
    Ok(pos)
}

fn layout_with_rings(
    g: &MoleculeGraph,
    hidden: &[bool],
    vadj: &[Vec<usize>],
) -> Result<Vec<Point2>, DepictError> {
    let n = g.atoms.len();
    let systems = ring_systems(g);
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
        let pos = layout_molecule(&g, &hidden, &vadj).unwrap();
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
