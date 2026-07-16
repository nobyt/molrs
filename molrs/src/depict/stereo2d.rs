//! 立体くさびの割当てと 2D 立体検証 (D9)。
//!
//! IUPAC 2006 (立体配置の図示): くさびの細端は立体中心に置く。
//! グラフは CIP ラベル (R/S) を持ちパリティを持たないため、くさびの
//! solid/hashed は「選んだ隣接を ±z に持ち上げ → R/S を再導出 → 入力
//! ラベルと一致する側を採用」で決める (再導出一致方式)。
//!
//! くさび先の選択優先順 (立体中心間・環内結合を可能な限り回避):
//! 1. 末端の非環単結合 (相手が非立体中心)
//! 2. 非末端の非環単結合 (相手が非立体中心)
//! 3. 隠し H の再表示 (空き方向に置いてくさび先にする)
//! 4. 非環単結合 (相手が立体中心)
//! 5. 環内結合 (最終手段)

use crate::geometry::Vec3;
use crate::graph::MoleculeGraph;
use crate::stereo::cip_ranks;

use super::chain_layout::{derive_ez, fallback_slots};
use super::point2::Point2;
use super::{Coords2D, Wedge, WedgeDir};

/// 3D 隣接ベクトル (CIP 順位降順) から R/S を求める。
/// 規約: 最下位置換基を奥に見て 1→2→3 が時計回りなら R。
/// 3 隣接 (孤立電子対) はファントム (単位ベクトル和の逆) を最下位に置く。
fn rs_from_3d(center: Vec3, mut nbrs: Vec<(usize, Vec3)>) -> Option<char> {
    if nbrs.len() < 3 || nbrs.len() > 4 {
        return None;
    }
    // CIP 順位降順 (rank 大 = 高順位)
    nbrs.sort_by_key(|x| std::cmp::Reverse(x.0));
    let v = |p: Vec3| p - center;
    let phantom;
    let v4 = if nbrs.len() == 4 {
        v(nbrs[3].1)
    } else {
        let mut sum = Vec3::ZERO;
        for &(_, p) in &nbrs {
            if let Some(u) = v(p).normalized() {
                sum = sum + u;
            }
        }
        phantom = -sum;
        phantom
    };
    let (v1, v2, v3) = (v(nbrs[0].1), v(nbrs[1].1), v(nbrs[2].1));
    let n = v1.cross(v2) + v2.cross(v3) + v3.cross(v1);
    let d = n.dot(v4);
    if d.abs() < 1e-9 {
        return None; // 退化 (平面的)
    }
    Some(if d < 0.0 { 'S' } else { 'R' })
}

fn lift(p: Point2, z: f64) -> Vec3 {
    Vec3::new(p.x, p.y, z)
}

/// 立体中心 c の R/S を 2D 座標 + くさびから再導出する。
/// c を細端とするくさびが付いた隣接を ±z に持ち上げる。
pub(crate) fn derive_rs(
    g: &MoleculeGraph,
    coords: &Coords2D,
    ranks: &[usize],
    c: usize,
) -> Option<char> {
    let mut nbrs: Vec<(usize, Vec3)> = Vec::new();
    for &nb in &g.adjacency[c] {
        if coords.hidden[nb] {
            continue;
        }
        let bi = bond_index(g, c, nb)?;
        let z = match &coords.wedge[bi] {
            Some(w) if w.narrow == c => match w.dir {
                WedgeDir::Up => 0.5,
                WedgeDir::Down => -0.5,
            },
            _ => 0.0,
        };
        nbrs.push((ranks[nb], lift(coords.pos[nb], z)));
    }
    rs_from_3d(lift(coords.pos[c], 0.0), nbrs)
}

fn bond_index(g: &MoleculeGraph, i: usize, j: usize) -> Option<usize> {
    g.bonds
        .iter()
        .position(|b| (b.begin_idx == i && b.end_idx == j) || (b.begin_idx == j && b.end_idx == i))
}

fn is_ring_bond(g: &MoleculeGraph, i: usize, j: usize) -> bool {
    g.ring_atom_sets.iter().any(|ring| {
        let n = ring.len();
        (0..n).any(|k| {
            let (a, b) = (ring[k], ring[(k + 1) % n]);
            (a == i && b == j) || (a == j && b == i)
        })
    })
}

/// 全立体中心にくさびを割り当てる。必要なら隠し H を再表示する。
pub(crate) fn assign_wedges(g: &MoleculeGraph, coords: &mut Coords2D) {
    let ranks = cip_ranks(g);
    let stereocenters: Vec<usize> = (0..g.atoms.len())
        .filter(|&i| g.atoms[i].chiral_tag.is_some())
        .collect();

    for &c in &stereocenters {
        let want = g.atoms[c].chiral_tag.unwrap();
        // 候補を優先順に列挙: (優先度, 隣接, 結合 idx)
        let mut candidates: Vec<(u8, usize, usize)> = Vec::new();
        for &nb in &g.adjacency[c] {
            let Some(bi) = bond_index(g, c, nb) else {
                continue;
            };
            if coords.wedge[bi].is_some() {
                continue; // 既に他の中心のくさび
            }
            if g.kekule_bond_orders[bi] != 1.0 {
                continue; // くさびは単結合のみ
            }
            let nb_stereo = g.atoms[nb].chiral_tag.is_some();
            let ring = is_ring_bond(g, c, nb);
            let prio = if coords.hidden[nb] {
                2 // 隠し H の再表示
            } else if ring {
                4
            } else if nb_stereo {
                3
            } else {
                let terminal = g.adjacency[nb]
                    .iter()
                    .filter(|&&x| !coords.hidden[x])
                    .count()
                    <= 1;
                if terminal {
                    0
                } else {
                    1
                }
            };
            candidates.push((prio, nb, bi));
        }
        candidates.sort();

        'cand: for &(_prio, nb, bi) in &candidates {
            // 隠し H は空き方向に再表示してから使う
            let mut revealed = false;
            if coords.hidden[nb] {
                let placed_dirs: Vec<f64> = g.adjacency[c]
                    .iter()
                    .copied()
                    .filter(|&x| !coords.hidden[x])
                    .map(|x| (coords.pos[x] - coords.pos[c]).angle())
                    .collect();
                let dir = fallback_slots(&placed_dirs, 1)[0];
                coords.pos[nb] = coords.pos[c] + Point2::from_angle(dir);
                coords.hidden[nb] = false;
                revealed = true;
            }
            for dir in [WedgeDir::Up, WedgeDir::Down] {
                coords.wedge[bi] = Some(Wedge { dir, narrow: c });
                if derive_rs(g, coords, &ranks, c) == Some(want) {
                    break 'cand; // 確定
                }
            }
            // この候補では一致しない → 巻き戻して次の候補へ
            coords.wedge[bi] = None;
            if revealed {
                coords.hidden[nb] = true;
                coords.pos[nb] = coords.pos[c];
            }
        }
        // 全候補で一致しなければくさびなしのまま (verify_stereo_2d /
        // D12 ゲートで検出される)
    }
}

/// 2D 座標 + くさびから全立体 (R/S と E/Z) を再導出し、入力と一致するか。
/// 不一致の原子/結合リストを返す (空 = 合格)。
pub fn verify_stereo_2d(g: &MoleculeGraph, coords: &Coords2D) -> Vec<String> {
    let ranks = cip_ranks(g);
    let mut failures = Vec::new();
    for (i, a) in g.atoms.iter().enumerate() {
        if let Some(want) = a.chiral_tag {
            match derive_rs(g, coords, &ranks, i) {
                Some(got) if got == want => {}
                got => failures.push(format!("atom {i}: want {want}, got {got:?}")),
            }
        }
    }
    for (bi, b) in g.bonds.iter().enumerate() {
        if let Some(want) = b.stereo {
            match derive_ez(g, &coords.pos, &coords.hidden, &ranks, bi) {
                Some(got) if got == want => {}
                got => failures.push(format!("bond {bi}: want {want}, got {got:?}")),
            }
        }
    }
    failures
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::depict::{compute_coords_2d, LayoutParams};
    use crate::graph::build_molecule_graph;

    fn coords_of(smiles: &str) -> (MoleculeGraph, Coords2D) {
        let g = build_molecule_graph(smiles).unwrap();
        let c = compute_coords_2d(&g, &LayoutParams::default()).unwrap();
        (g, c)
    }

    #[test]
    fn rs_from_3d_known_case() {
        // u4 (最下位) を -z に、u1..u3 を xy 面で反時計回り (0°,120°,240°)
        // → 奥の u4 から見ず、手前から見て 1→2→3 反時計回り = S
        let c = Vec3::ZERO;
        let nbrs = vec![
            (4, Vec3::new(1.0, 0.0, 0.0)),
            (3, Vec3::new(-0.5, 0.866, 0.0)),
            (2, Vec3::new(-0.5, -0.866, 0.0)),
            (1, Vec3::new(0.0, 0.0, -1.0)),
        ];
        assert_eq!(rs_from_3d(c, nbrs), Some('S'));
        // u4 を +z (手前) にすると R
        let nbrs_r = vec![
            (4, Vec3::new(1.0, 0.0, 0.0)),
            (3, Vec3::new(-0.5, 0.866, 0.0)),
            (2, Vec3::new(-0.5, -0.866, 0.0)),
            (1, Vec3::new(0.0, 0.0, 1.0)),
        ];
        assert_eq!(rs_from_3d(c, nbrs_r), Some('R'));
    }

    #[test]
    fn simple_stereocenters_roundtrip() {
        for smi in [
            "N[C@@H](C)C(=O)O", // L-alanine (S)
            "N[C@H](C)C(=O)O",  // D-alanine (R)
            "C[C@H](O)CC",      // 2-butanol
            "C[C@@H](O)CC",
            "Br[C@H](Cl)F",         // 全ハロゲン
            "C[C@H]1CCCC[C@@H]1C",  // 環上立体中心 2 つ
            "O[C@H]1CC[C@H](O)CC1", // 1,4-シクロヘキサンジオール
        ] {
            let (g, c) = coords_of(smi);
            let failures = verify_stereo_2d(&g, &c);
            assert!(failures.is_empty(), "{smi}: {failures:?}");
            // 立体中心には必ずくさびが付く
            for (i, a) in g.atoms.iter().enumerate() {
                if a.chiral_tag.is_some() {
                    let has_wedge = c
                        .wedge
                        .iter()
                        .any(|w| w.as_ref().is_some_and(|w| w.narrow == i));
                    assert!(has_wedge, "{smi}: no wedge at stereocenter {i}");
                }
            }
        }
    }

    #[test]
    fn ez_and_rs_combined() {
        let (g, c) = coords_of("C/C=C/[C@H](C)O");
        assert!(verify_stereo_2d(&g, &c).is_empty());
    }

    #[test]
    fn wedge_narrow_end_is_stereocenter() {
        let (g, c) = coords_of("C[C@H](O)CC");
        for w in c.wedge.iter().flatten() {
            assert!(g.atoms[w.narrow].chiral_tag.is_some());
        }
    }
}
