//! 立体化学の 3D 拘束と検証 (RUST_3D_PLAN.md C6)。
//!
//! - キラル体積拘束: CIP タグ付き原子について、近傍を CIP ランク降順
//!   (明示 H は最下位) に並べた 4 点の符号付き体積
//!   V = (p1−p4)·((p2−p4)×(p3−p4)) の符号を固定する。R → V < 0, S → V > 0
//!   (正四面体配置での解析値と一致することをテストで固定)。
//! - 平面項: sp2 中心 (3 配位) の [n1, n2, n3, 中心] 体積 → 0。
//!   面外変位は距離に二次でしか効かないため、芳香環の平面性はこの項が担保する。
//!   非芳香族二重結合には置換基クアッド [i, j, k, l] → 0 も追加する。
//! - 検証: 3D 座標から R/S と E/Z を再計算する。既存の CIP ランク
//!   (stereo.rs) を再利用するので、入力 SMILES の立体指定との
//!   全数一致ゲート (C9) がそのまま立体保存の証明になる。
//!
//! 注: アミド N の平面化は未実装 (混成推定が単結合 N を sp3 と見なすため)。
//! 必要になれば params の混成推定を拡張する。

use crate::conformer::minimize::VolumeConstraint;
use crate::conformer::params::{perceive_hybridization, Hybridization};
use crate::geometry::Vec3;
use crate::graph::MoleculeGraph;
use crate::stereo::cip_ranks;

/// CIP ランク (H ノードは最下位 0、kept 原子は rank+1)。
fn extended_ranks(g: &MoleculeGraph) -> Vec<usize> {
    let n_kept = g.parser_to_graph.iter().flatten().count();
    let base = cip_ranks(g);
    (0..g.atoms.len())
        .map(|i| if i < n_kept { base[i] + 1 } else { 0 })
        .collect()
}

/// タグ付きキラル原子の (中心, ランク降順の近傍列)。
/// 3 配位 (リン等) は中心自身を最下位として補う。
fn chiral_quads(g: &MoleculeGraph, ranks: &[usize]) -> Vec<(usize, [usize; 4], char)> {
    let mut out = Vec::new();
    for a in &g.atoms {
        let Some(tag) = a.chiral_tag else { continue };
        let mut nbrs: Vec<usize> = g.adjacency[a.idx].clone();
        if nbrs.len() < 3 || nbrs.len() > 4 {
            continue;
        }
        nbrs.sort_by(|&x, &y| ranks[y].cmp(&ranks[x]));
        let quad = if nbrs.len() == 4 {
            [nbrs[0], nbrs[1], nbrs[2], nbrs[3]]
        } else {
            // 孤立電子対を最下位扱い: 中心を第 4 点に
            [nbrs[0], nbrs[1], nbrs[2], a.idx]
        };
        out.push((a.idx, quad, tag));
    }
    out
}

/// 体積拘束一式を作る (キラル + 平面)。
pub(crate) fn build_volume_constraints(g: &MoleculeGraph) -> Vec<VolumeConstraint> {
    let ranks = extended_ranks(g);
    let hyb = perceive_hybridization(g);
    let mut out = Vec::new();

    // キラル体積 (典型値 |V| ≈ 3〜10 Å³)
    for (_, quad, tag) in chiral_quads(g, &ranks) {
        let (lower, upper) = match tag {
            'R' => (-30.0, -0.5),
            _ => (0.5, 30.0), // 'S'
        };
        out.push(VolumeConstraint {
            atoms: quad,
            lower,
            upper,
            weight: 3.0,
        });
    }

    // sp2 中心の平面化 (3 配位のみ)
    for a in &g.atoms {
        if hyb[a.idx] != Hybridization::Sp2 {
            continue;
        }
        let nbrs = &g.adjacency[a.idx];
        if nbrs.len() != 3 {
            continue;
        }
        out.push(VolumeConstraint {
            atoms: [nbrs[0], nbrs[1], nbrs[2], a.idx],
            lower: -0.1,
            upper: 0.1,
            weight: 10.0,
        });
    }

    // 芳香環の全体平面化: 連続する環内 4 原子のクアッド。
    // 隣接クアッドが 3 原子を共有するため平面が環全体に連鎖する
    // (sp2 中心クアッドは 2 原子しか共有せず蝶番になり得る)
    for ring in &g.ring_atom_sets {
        let m = ring.len();
        if m < 4 || !ring.iter().all(|&a| g.atoms[a].is_aromatic) {
            continue;
        }
        for t in 0..m {
            out.push(VolumeConstraint {
                atoms: [
                    ring[t],
                    ring[(t + 1) % m],
                    ring[(t + 2) % m],
                    ring[(t + 3) % m],
                ],
                lower: -0.05,
                upper: 0.05,
                weight: 10.0,
            });
        }
    }

    // 非芳香族二重結合の置換基クアッド (ねじれ防止)
    for b in &g.bonds {
        if b.bond_order != 2.0 {
            continue;
        }
        let (j, k) = (b.begin_idx, b.end_idx);
        for &i in &g.adjacency[j] {
            if i == k {
                continue;
            }
            for &l in &g.adjacency[k] {
                if l == j || l == i {
                    continue;
                }
                out.push(VolumeConstraint {
                    atoms: [i, j, k, l],
                    lower: -0.15,
                    upper: 0.15,
                    weight: 5.0,
                });
            }
        }
    }
    out
}

fn signed_volume(p: &[Vec3], q: &[usize; 4]) -> f64 {
    let (a, b, c, d) = (p[q[0]], p[q[1]], p[q[2]], p[q[3]]);
    (a - d).dot((b - d).cross(c - d))
}

/// 3D 座標からタグ付き原子の R/S を再計算する。
/// 返り値: (原子 idx, 計算された 'R'/'S')。縮退時は 'R'/'S' の代わりに '?'。
pub(crate) fn verify_atom_stereo(g: &MoleculeGraph, coords: &[Vec3]) -> Vec<(usize, char)> {
    let ranks = extended_ranks(g);
    chiral_quads(g, &ranks)
        .into_iter()
        .map(|(center, quad, _)| {
            let v = signed_volume(coords, &quad);
            let code = if v < -1e-3 {
                'R'
            } else if v > 1e-3 {
                'S'
            } else {
                '?'
            };
            (center, code)
        })
        .collect()
}

/// 3D 座標からタグ付き二重結合の E/Z を再計算する。
/// 返り値: (結合 idx, 計算された 'E'/'Z'/'?')。
pub(crate) fn verify_bond_stereo(g: &MoleculeGraph, coords: &[Vec3]) -> Vec<(usize, char)> {
    let n_kept = g.parser_to_graph.iter().flatten().count();
    let has = g.bonds.iter().any(|b| b.stereo.is_some());
    if !has {
        return Vec::new();
    }
    let ranks = cip_ranks(g);
    let mut out = Vec::new();
    for (ei, b) in g.bonds.iter().enumerate() {
        if b.stereo.is_none() {
            continue;
        }
        let (j, k) = (b.begin_idx, b.end_idx);
        let hi = |center: usize, other: usize| {
            g.adjacency[center]
                .iter()
                .filter(|&&x| x != other && x < n_kept)
                .max_by_key(|&&x| ranks[x])
                .copied()
        };
        let (Some(i), Some(l)) = (hi(j, k), hi(k, j)) else {
            continue;
        };
        let Some(axis) = (coords[k] - coords[j]).normalized() else {
            out.push((ei, '?'));
            continue;
        };
        // 軸に垂直な成分どうしの向き
        let u = {
            let w = coords[i] - coords[j];
            w - axis * w.dot(axis)
        };
        let v = {
            let w = coords[l] - coords[k];
            w - axis * w.dot(axis)
        };
        let d = u.dot(v);
        let code = if d > 1e-6 {
            'Z' // 高位置換基が同じ側
        } else if d < -1e-6 {
            'E'
        } else {
            '?'
        };
        out.push((ei, code));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformer::bounds::build_bounds;
    use crate::conformer::minimize::default_iterations;
    use crate::graph::build_molecule_graph;

    /// 手組みの正四面体で符号規約を固定する:
    /// p4 (最下位) を視線の奥に置き p1→p2→p3 が時計回り = R → V < 0。
    #[test]
    fn volume_sign_convention() {
        let coords = vec![
            Vec3::new(1.0, 0.0, 0.333),     // p1
            Vec3::new(-0.5, -0.866, 0.333), // p2 (上から見て時計回り)
            Vec3::new(-0.5, 0.866, 0.333),  // p3
            Vec3::new(0.0, 0.0, -1.0),      // p4 (奥)
        ];
        let v = signed_volume(&coords, &[0, 1, 2, 3]);
        assert!(v < -0.5, "clockwise (R) arrangement must be negative: {v}");
    }

    fn embed_with_stereo(smiles: &str, seed: u64) -> (MoleculeGraph, Vec<Vec3>) {
        let g = build_molecule_graph(smiles).expect("valid");
        let bm = build_bounds(&g);
        let volumes = build_volume_constraints(&g);
        let iters = default_iterations(g.atoms.len());
        let (coords, _) = crate::conformer::embed_and_refine(&bm, &volumes, seed, 20, iters)
            .unwrap_or_else(|| panic!("{smiles}: embedding failed"));
        (g, coords)
    }

    #[test]
    fn alanine_enantiomers_preserved() {
        for smi in ["N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"] {
            let (g, coords) = embed_with_stereo(smi, 1);
            for (idx, computed) in verify_atom_stereo(&g, &coords) {
                assert_eq!(
                    Some(computed),
                    g.atoms[idx].chiral_tag,
                    "{smi}: atom {idx} stereo not preserved"
                );
            }
        }
    }

    #[test]
    fn ez_preserved() {
        for smi in ["C/C=C/C", "C/C=C\\C", r"CC/C=C\CC", "C/C=C/C=C/C"] {
            let (g, coords) = embed_with_stereo(smi, 2);
            for (ei, computed) in verify_bond_stereo(&g, &coords) {
                assert_eq!(
                    Some(computed),
                    g.bonds[ei].stereo,
                    "{smi}: bond {ei} stereo not preserved"
                );
            }
        }
    }

    #[test]
    fn benzene_strictly_planar_with_constraints() {
        let (_, coords) = embed_with_stereo("c1ccccc1", 3);
        let ring: Vec<Vec3> = coords[..6].to_vec();
        let centroid = ring.iter().fold(Vec3::ZERO, |a, &b| a + b) / 6.0;
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
        assert!(rms < 0.03, "benzene planarity rms = {rms}");
    }

    #[test]
    fn multi_center_diastereomer() {
        // L-トレオニン (2S,3R): 2 中心が同時に保存されること
        let (g, coords) = embed_with_stereo("C[C@@H](O)[C@@H](N)C(=O)O", 4);
        let computed = verify_atom_stereo(&g, &coords);
        assert_eq!(computed.len(), 2);
        for (idx, code) in computed {
            assert_eq!(Some(code), g.atoms[idx].chiral_tag);
        }
    }
}
