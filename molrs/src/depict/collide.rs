//! 衝突検出と解消 (D8)。
//!
//! エネルギー = 非結合原子対の接近ペナルティ + 結合線交差ペナルティ。
//! 解消は回転可能結合 (非環単結合) まわりの部分木変換の貪欲探索:
//! - 結合軸での部分木鏡映 (GR-4.3 の対称再配分に相当)
//! - ピボット原子まわりの ±30°/±60°/±90° 回転
//!
//! 全変換は 30° 量子化を保存し、適用前に全 E/Z の再導出一致を検査する
//! (違反する変換は棄却)。決定的な貪欲降下で、改善がなくなるか反復上限で
//! 停止する。解消しきれない衝突はゲート (D12) の例外リストで追跡する。

use std::f64::consts::PI;

use crate::graph::MoleculeGraph;
use crate::stereo::cip_ranks;

use super::chain_layout::derive_ez;
use super::point2::Point2;
use super::LayoutParams;

/// 非結合原子対がこの距離未満に近づいたらペナルティ。
const CLASH_DIST: f64 = 0.6;
/// 結合線交差 1 件あたりのペナルティ。
const CROSS_PENALTY: f64 = 0.25;

/// 線分の厳密交差 (端点共有・接触は交差とみなさない)。
fn segments_cross(p1: Point2, p2: Point2, p3: Point2, p4: Point2) -> bool {
    let d1 = (p2 - p1).cross(p3 - p1);
    let d2 = (p2 - p1).cross(p4 - p1);
    let d3 = (p4 - p3).cross(p1 - p3);
    let d4 = (p4 - p3).cross(p2 - p3);
    d1 * d2 < -1e-12 && d3 * d4 < -1e-12
}

/// 現在座標の衝突エネルギー。
pub(crate) fn collision_energy(
    g: &MoleculeGraph,
    pos: &[Point2],
    hidden: &[bool],
    frag: &[usize],
) -> f64 {
    let mut e = 0.0;
    // 非結合原子対の接近
    for (k, &i) in frag.iter().enumerate() {
        for &j in &frag[k + 1..] {
            if g.adjacency[i].contains(&j) {
                continue;
            }
            let d = pos[i].distance(pos[j]);
            if d < CLASH_DIST {
                let x = CLASH_DIST - d;
                e += x * x;
            }
        }
    }
    // 結合線交差 (原子を共有しない可視結合対)
    let bonds: Vec<(usize, usize)> = g
        .bonds
        .iter()
        .filter(|b| !hidden[b.begin_idx] && !hidden[b.end_idx])
        .map(|b| (b.begin_idx, b.end_idx))
        .filter(|&(i, j)| frag.contains(&i) && frag.contains(&j))
        .collect();
    for (k, &(a, b)) in bonds.iter().enumerate() {
        for &(c, d) in &bonds[k + 1..] {
            if a == c || a == d || b == c || b == d {
                continue;
            }
            if segments_cross(pos[a], pos[b], pos[c], pos[d]) {
                e += CROSS_PENALTY;
            }
        }
    }
    e
}

/// 部分木変換の候補。
#[derive(Clone, Copy)]
enum Transform {
    /// 結合軸 (pivot→child) での鏡映
    Reflect,
    /// pivot まわりの回転 (rad)
    Rotate(f64),
}

const TRANSFORMS: [Transform; 7] = [
    Transform::Reflect,
    Transform::Rotate(PI / 6.0),
    Transform::Rotate(-PI / 6.0),
    Transform::Rotate(PI / 3.0),
    Transform::Rotate(-PI / 3.0),
    Transform::Rotate(PI / 2.0),
    Transform::Rotate(-PI / 2.0),
];

fn apply_transform(
    tf: Transform,
    pivot: Point2,
    child: Point2,
    subtree: &[usize],
    pos: &mut [Point2],
) {
    match tf {
        Transform::Reflect => {
            let Some(u) = (child - pivot).normalized() else {
                return;
            };
            for &a in subtree {
                let v = pos[a] - pivot;
                pos[a] = pivot + u * (2.0 * v.dot(u)) - v;
            }
        }
        Transform::Rotate(theta) => {
            for &a in subtree {
                pos[a] = (pos[a] - pivot).rotated(theta) + pivot;
            }
        }
    }
}

/// フラグメント内の衝突を貪欲に解消する。
pub(crate) fn resolve_collisions(
    g: &MoleculeGraph,
    pos: &mut [Point2],
    hidden: &[bool],
    vadj: &[Vec<usize>],
    frag: &[usize],
    params: &LayoutParams,
) {
    let n = g.atoms.len();
    let mut energy = collision_energy(g, pos, hidden, frag);
    if energy < 1e-9 {
        return;
    }
    let ranks = cip_ranks(g);
    let stereo_bonds: Vec<usize> = g
        .bonds
        .iter()
        .enumerate()
        .filter(|(_, b)| b.stereo.is_some())
        .map(|(bi, _)| bi)
        .collect();

    // 回転可能結合: 可視・単結合 (ケクレ次数 1)・非環 (部分木が閉じない)
    struct Rotor {
        pivot: usize,
        child: usize,
        subtree: Vec<usize>,
    }
    let mut rotors: Vec<Rotor> = Vec::new();
    for (bi, b) in g.bonds.iter().enumerate() {
        let (u, v) = (b.begin_idx, b.end_idx);
        if hidden[u] || hidden[v] || g.kekule_bond_orders[bi] != 1.0 {
            continue;
        }
        if !frag.contains(&u) || !frag.contains(&v) {
            continue;
        }
        // v 側の部分木 (この結合の辺を通らない到達集合)
        let mut reach = vec![false; n];
        let mut stack = vec![v];
        reach[v] = true;
        while let Some(x) = stack.pop() {
            for &nb in &vadj[x] {
                if (x == v && nb == u) || (x == u && nb == v) {
                    continue;
                }
                if !reach[nb] {
                    reach[nb] = true;
                    stack.push(nb);
                }
            }
        }
        if reach[u] {
            continue; // 環内結合
        }
        let side_v: Vec<usize> = frag.iter().copied().filter(|&a| reach[a]).collect();
        // 小さい側を動かす
        if side_v.len() * 2 <= frag.len() {
            rotors.push(Rotor {
                pivot: u,
                child: v,
                subtree: side_v,
            });
        } else {
            let side_u: Vec<usize> = frag.iter().copied().filter(|&a| !reach[a]).collect();
            rotors.push(Rotor {
                pivot: v,
                child: u,
                subtree: side_u,
            });
        }
    }

    let mut scratch: Vec<Point2> = pos.to_vec();
    for _ in 0..params.max_collision_iters {
        let mut best: Option<(f64, usize, usize)> = None; // (energy, rotor, tf)
        for (ri, rot) in rotors.iter().enumerate() {
            for (ti, &tf) in TRANSFORMS.iter().enumerate() {
                scratch.copy_from_slice(pos);
                apply_transform(
                    tf,
                    pos[rot.pivot],
                    pos[rot.child],
                    &rot.subtree,
                    &mut scratch,
                );
                // E/Z を壊す変換は棄却
                let ez_ok = stereo_bonds.iter().all(|&bi| {
                    derive_ez(g, &scratch, hidden, &ranks, bi)
                        .map(|got| Some(got) == g.bonds[bi].stereo)
                        .unwrap_or(true)
                });
                if !ez_ok {
                    continue;
                }
                let e = collision_energy(g, &scratch, hidden, frag);
                if e < energy - 1e-9 && (best.is_none() || e < best.unwrap().0 - 1e-12) {
                    best = Some((e, ri, ti));
                }
            }
        }
        let Some((e, ri, ti)) = best else { break };
        let rot = &rotors[ri];
        apply_transform(
            TRANSFORMS[ti],
            pos[rot.pivot],
            pos[rot.child],
            &rot.subtree,
            pos,
        );
        energy = e;
        if energy < 1e-9 {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::depict::chain_layout::{hidden_h_flags, visible_adjacency};
    use crate::graph::build_molecule_graph;

    #[test]
    fn segments_cross_basic() {
        let p = |x: f64, y: f64| Point2::new(x, y);
        assert!(segments_cross(
            p(0.0, 0.0),
            p(1.0, 1.0),
            p(0.0, 1.0),
            p(1.0, 0.0)
        ));
        assert!(!segments_cross(
            p(0.0, 0.0),
            p(1.0, 0.0),
            p(0.0, 1.0),
            p(1.0, 1.0)
        ));
        // 端点共有は交差ではない
        assert!(!segments_cross(
            p(0.0, 0.0),
            p(1.0, 0.0),
            p(0.0, 0.0),
            p(0.0, 1.0)
        ));
    }

    #[test]
    fn crowded_molecule_energy_reduced() {
        // 隣接環位置の長鎖 2 本 (二等分方向が 60° 差 → 先で衝突しがち)
        let g = build_molecule_graph("CCCCCc1ccccc1CCCCC").unwrap();
        let hidden = hidden_h_flags(&g);
        let vadj = visible_adjacency(&g, &hidden);
        let mut pos = crate::depict::place::layout_molecule(
            &g,
            &hidden,
            &vadj,
            &crate::depict::LayoutParams::default(),
        )
        .unwrap();
        let frag: Vec<usize> = (0..g.atoms.len()).filter(|&i| !hidden[i]).collect();
        let e = collision_energy(&g, &pos, &hidden, &frag);
        resolve_collisions(
            &g,
            &mut pos,
            &hidden,
            &vadj,
            &frag,
            &crate::depict::LayoutParams::default(),
        );
        let e2 = collision_energy(&g, &pos, &hidden, &frag);
        assert!(e2 <= e + 1e-12, "energy must not increase: {e} -> {e2}");
    }
}
