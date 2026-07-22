//! 主鎖の探索とロカント番号付け (S3.2, Python chain_finder.py の移植)。
//!
//! IUPAC 2013 Blue Book P-44 の選択規則を順に適用して主鎖を決める:
//! 1. 主特性基の炭素を含む最長鎖 2. 最長炭素鎖 3. 多重結合数最多
//! 4. 置換基数最多 5. ロカント集合最小。

use molrs::graph::{get_bond_order, MoleculeGraph};

use crate::constants;
use crate::functional_group::FunctionalGroup;

/// 主鎖: ロカント 1 から始まる順序付き炭素インデックス。
pub struct PrincipalChain {
    pub atom_indices: Vec<usize>,
}

impl PrincipalChain {
    pub fn length(&self) -> usize {
        self.atom_indices.len()
    }
    /// 原子インデックス → ロカント (1 始まり)。鎖外は None。
    pub fn locant_of(&self, atom: usize) -> Option<usize> {
        self.atom_indices.iter().position(|&a| a == atom).map(|i| i + 1)
    }
}

fn is_c(g: &MoleculeGraph, i: usize) -> bool {
    g.atoms[i].symbol == "C"
}

/// 非水素・非炭素の重原子か (置換基の存在判定用)。
fn carbon_indices(g: &MoleculeGraph) -> Vec<usize> {
    (0..g.atoms.len()).filter(|&i| is_c(g, i)).collect()
}

fn non_ring_carbon_indices(g: &MoleculeGraph) -> Vec<usize> {
    (0..g.atoms.len())
        .filter(|&i| is_c(g, i) && !g.atoms[i].in_ring)
        .collect()
}

/// 主鎖を選択してロカント方向を決める。
pub fn find_principal_chain(
    g: &MoleculeGraph,
    principal: Option<&FunctionalGroup>,
) -> PrincipalChain {
    let has_ring_c = g.atoms.iter().any(|a| a.symbol == "C" && a.in_ring);
    let mut c_idxs = if has_ring_c {
        non_ring_carbon_indices(g)
    } else {
        carbon_indices(g)
    };
    if c_idxs.is_empty() {
        // 環のみ (非環炭素なし) — 呼び出し側で環処理に回すため空鎖
        return PrincipalChain { atom_indices: Vec::new() };
    }

    // 主基がニトリル系でなければ末端 C≡N の炭素は cyano 置換基 → 鎖から除外
    let nitrile_types = [
        "nitrile",
        "dinitrile",
        "isocyanate",
        "cyanate",
        "isothiocyanate",
        "thiocyanate",
    ];
    let is_nitrile_principal = principal.is_some_and(|p| nitrile_types.contains(&p.group_type));
    if !is_nitrile_principal {
        let mut cyano_c = Vec::new();
        for &ci in &c_idxs {
            for &nb in &g.adjacency[ci] {
                if g.atoms[nb].symbol == "N" && get_bond_order(g, ci, nb) == 3.0 {
                    let n_heavy = g.adjacency[nb]
                        .iter()
                        .any(|&x| x != ci && g.atoms[x].symbol != "H");
                    if !n_heavy {
                        cyano_c.push(ci);
                        break;
                    }
                }
            }
        }
        if !cyano_c.is_empty() {
            let filtered: Vec<usize> =
                c_idxs.iter().copied().filter(|c| !cyano_c.contains(c)).collect();
            c_idxs = if filtered.is_empty() { cyano_c } else { filtered };
        }
    }

    // 主基に属する炭素
    let required: Vec<usize> = principal
        .map(|p| {
            p.atom_indices
                .iter()
                .copied()
                .filter(|&ai| is_c(g, ai))
                .collect()
        })
        .unwrap_or_default();

    let mut all_paths = enumerate_carbon_paths(g, &c_idxs);
    if !required.is_empty() {
        let filtered: Vec<Vec<usize>> = all_paths
            .iter()
            .filter(|p| required.iter().all(|r| p.contains(r)))
            .cloned()
            .collect();
        if !filtered.is_empty() {
            all_paths = filtered;
        }
    }
    if all_paths.is_empty() {
        return PrincipalChain { atom_indices: c_idxs[..1].to_vec() };
    }

    let best = select_best_path(g, &all_paths);
    let oriented = orient_chain(g, &best, principal);
    PrincipalChain { atom_indices: oriented }
}

/// 炭素のみを通る全単純経路 (端点ペアごとに 1 本) を DFS 列挙する。
fn enumerate_carbon_paths(g: &MoleculeGraph, c_idxs: &[usize]) -> Vec<Vec<usize>> {
    let c_set: std::collections::HashSet<usize> = c_idxs.iter().copied().collect();
    let c_adj: std::collections::HashMap<usize, Vec<usize>> = c_idxs
        .iter()
        .map(|&c| {
            (
                c,
                g.adjacency[c]
                    .iter()
                    .copied()
                    .filter(|nb| c_set.contains(nb))
                    .collect(),
            )
        })
        .collect();

    let mut all_paths: Vec<Vec<usize>> = Vec::new();
    let mut seen_pairs: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();

    fn dfs(
        cur: usize,
        vis: &mut Vec<bool>,
        path: &mut Vec<usize>,
        c_adj: &std::collections::HashMap<usize, Vec<usize>>,
        all_paths: &mut Vec<Vec<usize>>,
        seen: &mut std::collections::HashSet<(usize, usize)>,
    ) {
        path.push(cur);
        vis[cur] = true;
        let mut extended = false;
        for &nb in &c_adj[&cur] {
            if !vis[nb] {
                extended = true;
                dfs(nb, vis, path, c_adj, all_paths, seen);
            }
        }
        if !extended {
            let (a, b) = (path[0], *path.last().unwrap());
            let pair = (a.min(b), a.max(b));
            if seen.insert(pair) {
                all_paths.push(path.clone());
            }
        }
        path.pop();
        vis[cur] = false;
    }

    let n = g.atoms.len();
    for &start in c_idxs {
        let mut vis = vec![false; n];
        let mut path = Vec::new();
        dfs(start, &mut vis, &mut path, &c_adj, &mut all_paths, &mut seen_pairs);
    }
    if all_paths.is_empty() {
        all_paths = c_idxs.iter().map(|&c| vec![c]).collect();
    }
    all_paths
}

fn count_multiple_bonds(g: &MoleculeGraph, path: &[usize]) -> usize {
    (0..path.len().saturating_sub(1))
        .filter(|&i| get_bond_order(g, path[i], path[i + 1]) >= 2.0)
        .count()
}

fn count_substituents(g: &MoleculeGraph, path: &[usize]) -> usize {
    let set: std::collections::HashSet<usize> = path.iter().copied().collect();
    let mut count = 0;
    for &c in path {
        for &nb in &g.adjacency[c] {
            if !set.contains(&nb) && g.atoms[nb].symbol != "H" {
                count += 1;
            }
        }
    }
    count
}

fn compute_locant_set(g: &MoleculeGraph, path: &[usize]) -> Vec<i64> {
    let set: std::collections::HashSet<usize> = path.iter().copied().collect();
    let mut locants = Vec::new();
    for (i, &c) in path.iter().enumerate() {
        let locant = (i + 1) as i64;
        for &nb in &g.adjacency[c] {
            if !set.contains(&nb) && g.atoms[nb].symbol != "H" {
                locants.push(locant);
            }
        }
        if i + 1 < path.len() && get_bond_order(g, c, path[i + 1]) >= 2.0 {
            locants.push(locant);
        }
    }
    locants.sort_unstable();
    locants
}

fn locant_set_for_sorting(g: &MoleculeGraph, path: &[usize]) -> Vec<i64> {
    let fwd = compute_locant_set(g, path);
    let rev: Vec<usize> = path.iter().rev().copied().collect();
    let rev_set = compute_locant_set(g, &rev);
    fwd.min(rev_set)
}

fn select_best_path(g: &MoleculeGraph, paths: &[Vec<usize>]) -> Vec<usize> {
    // sort_key = (length, mb, subst, [-l for l in locant_set]) を最大化
    paths
        .iter()
        .max_by(|a, b| {
            let ka = (
                a.len(),
                count_multiple_bonds(g, a),
                count_substituents(g, a),
                locant_set_for_sorting(g, a).iter().map(|&l| -l).collect::<Vec<_>>(),
            );
            let kb = (
                b.len(),
                count_multiple_bonds(g, b),
                count_substituents(g, b),
                locant_set_for_sorting(g, b).iter().map(|&l| -l).collect::<Vec<_>>(),
            );
            ka.cmp(&kb)
        })
        .cloned()
        .unwrap()
}

fn sub_locants(g: &MoleculeGraph, path: &[usize]) -> Vec<usize> {
    let set: std::collections::HashSet<usize> = path.iter().copied().collect();
    let mut locs = Vec::new();
    for (i, &c) in path.iter().enumerate() {
        for &nb in &g.adjacency[c] {
            if !set.contains(&nb) && g.atoms[nb].symbol != "H" {
                locs.push(i + 1);
            }
        }
    }
    locs.sort_unstable();
    locs
}

fn choose_by_substituent_locants(g: &MoleculeGraph, fwd: &[usize], rev: &[usize]) -> Vec<usize> {
    if sub_locants(g, fwd) <= sub_locants(g, rev) {
        fwd.to_vec()
    } else {
        rev.to_vec()
    }
}

fn mb_locs(g: &MoleculeGraph, path: &[usize]) -> Vec<usize> {
    let mut locs = Vec::new();
    for i in 0..path.len().saturating_sub(1) {
        if get_bond_order(g, path[i], path[i + 1]) >= 2.0 {
            locs.push(i + 1);
        }
    }
    locs.sort_unstable();
    locs
}

fn db_locs(g: &MoleculeGraph, path: &[usize]) -> Vec<usize> {
    let mut locs = Vec::new();
    for i in 0..path.len().saturating_sub(1) {
        if get_bond_order(g, path[i], path[i + 1]) == 2.0 {
            locs.push(i + 1);
        }
    }
    locs.sort_unstable();
    locs
}

fn choose_by_multiple_bond_locants(g: &MoleculeGraph, fwd: &[usize], rev: &[usize]) -> Vec<usize> {
    let (lf, lr) = (mb_locs(g, fwd), mb_locs(g, rev));
    if lf != lr {
        return if lf < lr { fwd.to_vec() } else { rev.to_vec() };
    }
    let (df, dr) = (db_locs(g, fwd), db_locs(g, rev));
    if df != dr {
        return if df < dr { fwd.to_vec() } else { rev.to_vec() };
    }
    choose_by_substituent_locants(g, fwd, rev)
}

fn orient_chain(
    g: &MoleculeGraph,
    path: &[usize],
    principal: Option<&FunctionalGroup>,
) -> Vec<usize> {
    if path.len() == 1 {
        return path.to_vec();
    }
    let fwd = path.to_vec();
    let rev: Vec<usize> = path.iter().rev().copied().collect();

    let gtype = principal.map(|p| p.group_type).unwrap_or("alkane");
    if gtype == "alkane" {
        return choose_by_substituent_locants(g, &fwd, &rev);
    }
    if gtype == "alkene" || gtype == "alkyne" {
        return choose_by_multiple_bond_locants(g, &fwd, &rev);
    }
    let spec = constants::FUNCTIONAL_GROUP_MAP.get(gtype);
    let Some(spec) = spec else {
        return choose_by_substituent_locants(g, &fwd, &rev);
    };
    let grp_carbons: std::collections::HashSet<usize> = principal
        .map(|p| {
            p.atom_indices
                .iter()
                .copied()
                .filter(|&ai| is_c(g, ai))
                .collect()
        })
        .unwrap_or_default();

    if spec.anchor_c1 {
        let both_ends =
            grp_carbons.contains(&path[0]) && grp_carbons.contains(path.last().unwrap());
        if both_ends {
            let (mf, mr) = (mb_locs(g, &fwd), mb_locs(g, &rev));
            if mf < mr {
                return fwd;
            }
            if mr < mf {
                return rev;
            }
            return choose_by_substituent_locants(g, &fwd, &rev);
        }
        if grp_carbons.contains(&path[0]) {
            return fwd;
        }
        if grp_carbons.contains(path.last().unwrap()) {
            return rev;
        }
        return fwd;
    }

    if spec.needs_locant {
        let loc_fwd = fwd
            .iter()
            .position(|c| grp_carbons.contains(c))
            .map(|i| i + 1)
            .unwrap_or(fwd.len());
        let loc_rev = rev
            .iter()
            .position(|c| grp_carbons.contains(c))
            .map(|i| i + 1)
            .unwrap_or(rev.len());
        if loc_fwd != loc_rev {
            return if loc_fwd < loc_rev { fwd } else { rev };
        }
        let (mf, mr) = (mb_locs(g, &fwd), mb_locs(g, &rev));
        if mf != mr {
            return if mf < mr { fwd } else { rev };
        }
        return choose_by_substituent_locants(g, &fwd, &rev);
    }

    choose_by_substituent_locants(g, &fwd, &rev)
}

/// 主鎖上の多重結合ロカント。返り値 (ene, yne)。
pub fn multiple_bond_locants(g: &MoleculeGraph, chain: &PrincipalChain) -> (Vec<usize>, Vec<usize>) {
    let path = &chain.atom_indices;
    let mut ene = Vec::new();
    let mut yne = Vec::new();
    for i in 0..path.len().saturating_sub(1) {
        let bo = get_bond_order(g, path[i], path[i + 1]);
        if bo == 2.0 {
            ene.push(i + 1);
        } else if bo == 3.0 {
            yne.push(i + 1);
        }
    }
    ene.sort_unstable();
    yne.sort_unstable();
    (ene, yne)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::functional_group::detect_groups;
    use molrs::graph::build_molecule_graph;

    fn chain_locants(smiles: &str) -> Vec<usize> {
        let g = build_molecule_graph(smiles).unwrap();
        let groups = detect_groups(&g);
        let principal = crate::functional_group::principal_group(&groups);
        find_principal_chain(&g, principal).atom_indices
    }

    #[test]
    fn straight_chain() {
        assert_eq!(chain_locants("CCCC").len(), 4);
        assert_eq!(chain_locants("CC(C)CC").len(), 4); // longest = 4
    }

    #[test]
    fn acid_anchors_c1() {
        // ブタン酸: COOH 炭素が C1
        let g = build_molecule_graph("CCCC(=O)O").unwrap();
        let groups = detect_groups(&g);
        let principal = crate::functional_group::principal_group(&groups);
        let chain = find_principal_chain(&g, principal);
        assert_eq!(chain.length(), 4);
        // C1 は COOH 炭素 (=O を持つ)
        let c1 = chain.atom_indices[0];
        assert!(g.adjacency[c1].iter().any(|&nb| g.atoms[nb].symbol == "O"
            && get_bond_order(&g, c1, nb) == 2.0));
    }
}
