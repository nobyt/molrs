//! ETKDG 実験トーションの収集 (C10)。
//! RDKit `getExperimentalTorsions` (TorsionPreferences.cpp) の移植。
//!
//! ライブラリの 365 パターンをファイル順に照合し、**中央結合ごとに
//! 最初にマッチしたパターン**のポテンシャルを採用する (先勝ち)。
//! 2 本以上の結合を共有する環ペア (橋かけ環系) の結合は除外する。

use std::sync::LazyLock;

use crate::conformer::torsion_lib::TORSION_LIB;
use crate::graph::MoleculeGraph;
use crate::smarts::{parse_smarts, smarts_matches, MolView, SmartsPattern};

/// パース済みライブラリエントリ。
struct LibEntry {
    pattern: SmartsPattern,
    /// マップ番号 1..4 に対応するパターン原子インデックス
    map_idx: [usize; 4],
    signs: [i8; 6],
    v: [f64; 6],
}

static PARSED_LIB: LazyLock<Vec<LibEntry>> = LazyLock::new(|| {
    TORSION_LIB
        .iter()
        .map(|&(smarts, signs, v)| {
            let pattern = parse_smarts(smarts)
                .unwrap_or_else(|e| panic!("torsion lib SMARTS fails to parse: {smarts}: {e}"));
            let mut map_idx = [usize::MAX; 4];
            for (i, m) in pattern.atom_maps.iter().enumerate() {
                if let Some(n @ 1..=4) = m {
                    map_idx[(*n - 1) as usize] = i;
                }
            }
            assert!(
                map_idx.iter().all(|&x| x != usize::MAX),
                "torsion lib SMARTS lacks maps 1-4: {smarts}"
            );
            LibEntry {
                pattern,
                map_idx,
                signs,
                v,
            }
        })
        .collect()
});

/// 採用された実験トーション。
#[derive(Debug, Clone)]
pub struct ExpTorsion {
    pub atoms: [usize; 4],
    pub signs: [i8; 6],
    pub v: [f64; 6],
    /// ライブラリ内のパターン番号 (フィクスチャ照合・デバッグ用)
    pub torsion_idx: usize,
}

/// 分子の実験トーション一式を収集する。
pub(crate) fn collect_exp_torsions(g: &MoleculeGraph) -> Vec<ExpTorsion> {
    let n = g.atoms.len();
    if n < 4 {
        return Vec::new();
    }
    let view = MolView::build(g);

    // 結合 idx の逆引き
    let bond_index = |u: usize, v: usize| -> Option<usize> {
        g.bonds.iter().position(|b| {
            (b.begin_idx == u && b.end_idx == v) || (b.begin_idx == v && b.end_idx == u)
        })
    };

    // 環ごとの結合集合
    let ring_bonds: Vec<Vec<usize>> = g
        .ring_atom_sets
        .iter()
        .map(|ring| {
            (0..ring.len())
                .filter_map(|t| bond_index(ring[t], ring[(t + 1) % ring.len()]))
                .collect()
        })
        .collect();

    // 2 本以上の結合を共有する環ペア → 両環の結合を除外
    // (マクロ環 (9 員以上) 同士は除外しない — RDKit MIN_MACROCYCLE_SIZE)
    let mut excluded = vec![false; g.bonds.len()];
    for (ri, ra) in ring_bonds.iter().enumerate() {
        for (rj, rb) in ring_bonds.iter().enumerate().skip(ri + 1) {
            let size_a = g.ring_atom_sets[ri].len();
            let size_b = g.ring_atom_sets[rj].len();
            if size_a >= 9 && size_b >= 9 {
                continue;
            }
            let shared = ra.iter().filter(|b| rb.contains(b)).count();
            if shared > 1 {
                if size_a < 9 {
                    for &b in ra {
                        excluded[b] = true;
                    }
                }
                if size_b < 9 {
                    for &b in rb {
                        excluded[b] = true;
                    }
                }
            }
        }
    }
    // 結合ごとの所属環数
    let mut n_bond_rings = vec![0usize; g.bonds.len()];
    for rb in &ring_bonds {
        for &b in rb {
            n_bond_rings[b] += 1;
        }
    }

    let mut done = vec![false; g.bonds.len()];
    let mut out = Vec::new();
    for (ti, entry) in PARSED_LIB.iter().enumerate() {
        let matches = smarts_matches(&view, &entry.pattern);
        for m in &matches {
            let a1 = m[entry.map_idx[0]];
            let a2 = m[entry.map_idx[1]];
            let a3 = m[entry.map_idx[2]];
            let a4 = m[entry.map_idx[3]];
            let Some(bid) = bond_index(a2, a3) else {
                continue;
            };
            if excluded[bid] || n_bond_rings[bid] > 3 {
                done[bid] = true;
            }
            if !done[bid] {
                done[bid] = true;
                out.push(ExpTorsion {
                    atoms: [a1, a2, a3, a4],
                    signs: entry.signs,
                    v: entry.v,
                    torsion_idx: ti,
                });
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    #[test]
    fn whole_library_parses() {
        assert_eq!(PARSED_LIB.len(), 365);
    }

    #[test]
    fn amide_gets_a_torsion() {
        // N-メチルアセトアミド: アミド C-N 結合に実験トーションが付く
        let g = build_molecule_graph("CC(=O)NC").unwrap();
        let ts = collect_exp_torsions(&g);
        assert!(!ts.is_empty());
        // 中央結合 (atoms[1], atoms[2]) = C(=O)-N
        let has_amide = ts.iter().any(|t| {
            let (a, b) = (t.atoms[1].min(t.atoms[2]), t.atoms[1].max(t.atoms[2]));
            (a, b) == (1, 3)
        });
        assert!(has_amide, "torsions: {ts:?}");
    }

    #[test]
    fn one_torsion_per_bond() {
        let g = build_molecule_graph("CCCCO").unwrap();
        let ts = collect_exp_torsions(&g);
        let mut bonds: Vec<(usize, usize)> = ts
            .iter()
            .map(|t| (t.atoms[1].min(t.atoms[2]), t.atoms[1].max(t.atoms[2])))
            .collect();
        let total = bonds.len();
        bonds.sort_unstable();
        bonds.dedup();
        assert_eq!(bonds.len(), total, "one torsion per central bond");
    }
}
