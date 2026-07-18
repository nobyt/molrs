//! InChI 正準番号付け (I3)。
//!
//! InChI は分子を成分に分け、各成分内の重原子に 1 始まりの正準番号を付ける。
//! 番号は初期色 (元素・次数・H 数・電荷など) を Morgan 系の反復精緻化で
//! 細分し、同値類が残る場合は正準最小化 (各候補を試して最小の接続表を選ぶ)
//! で確定する。
//!
//! 検証: RDKit AuxInfo の `/N:` フィールド (公式の正準番号) と一致させる。
//!
//! 本ファイルは I3 で番号付け本体を実装する。まず成分分解
//! ([`connected_components`]) を提供する (formula.rs が使用)。

use crate::graph::MoleculeGraph;

/// 重原子 (非 H) のインデックス一覧。番号付け本体 (I3) で使用する。
#[allow(dead_code)]
pub(crate) fn heavy_atoms(g: &MoleculeGraph) -> Vec<usize> {
    (0..g.atoms.len())
        .filter(|&i| g.atoms[i].symbol != "H")
        .collect()
}

/// 重原子の連結成分。各成分は原子インデックスの Vec。
/// 成分の順序は最小原子インデックスの昇順 (安定・決定的)。
pub(crate) fn connected_components(g: &MoleculeGraph) -> Vec<Vec<usize>> {
    let n = g.atoms.len();
    let is_heavy = |i: usize| g.atoms[i].symbol != "H";
    let mut comp = vec![usize::MAX; n];
    let mut n_comp = 0;
    for start in 0..n {
        if !is_heavy(start) || comp[start] != usize::MAX {
            continue;
        }
        let mut stack = vec![start];
        comp[start] = n_comp;
        while let Some(v) = stack.pop() {
            for &nb in &g.adjacency[v] {
                if is_heavy(nb) && comp[nb] == usize::MAX {
                    comp[nb] = n_comp;
                    stack.push(nb);
                }
            }
        }
        n_comp += 1;
    }
    let mut comps: Vec<Vec<usize>> = vec![Vec::new(); n_comp];
    for i in 0..n {
        if is_heavy(i) && comp[i] != usize::MAX {
            comps[comp[i]].push(i);
        }
    }
    // 各成分内は昇順、成分自体は最小 idx 昇順 (push 順で既に満たす)
    comps
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    #[test]
    fn components_split() {
        let g = build_molecule_graph("CCO.O").unwrap();
        let comps = connected_components(&g);
        assert_eq!(comps.len(), 2);
        // 重原子のみ (H は含まない)
        for c in &comps {
            for &a in c {
                assert_ne!(g.atoms[a].symbol, "H");
            }
        }
        assert_eq!(comps[0].len(), 3); // C C O
        assert_eq!(comps[1].len(), 1); // O
    }

    #[test]
    fn heavy_atoms_excludes_h() {
        let g = build_molecule_graph("C").unwrap();
        assert_eq!(heavy_atoms(&g).len(), 1);
    }
}
