//! 置換基の命名 (S3.3, Python substituent.py の簡約移植)。
//!
//! v1 対応: ハロゲン、O 系 (hydroxy/oxo/alkoxy)、N 系 (amino/nitro)、
//! S 系 (sulfanyl/alkylsulfanyl)、炭素鎖 (直鎖・分岐、付け根に最小ロカント)。
//! アリール・複雑複合置換基は未対応 (None を返す → 呼び出し側で失敗)。

use std::collections::HashSet;

use molrs::graph::{get_bond_order, MoleculeGraph};

use crate::constants::{halogen_prefix, CHAIN_PREFIX, MULTIPLIER};

/// 主鎖上の炭素の置換基と、アミン主基の N 上の置換基を列挙する。
/// 返り値 (chain_subs (ロカント,名前), n_subs (N番号,名前))。
/// いずれかが命名不能、または全重原子を覆えなければ None。
#[allow(clippy::type_complexity)]
pub fn collect_substituents(
    g: &MoleculeGraph,
    chain: &[usize],
    principal_atoms: &[usize],
) -> Option<(Vec<(usize, String)>, Vec<(usize, String)>)> {
    let chain_set: HashSet<usize> = chain.iter().copied().collect();
    let principal_set: HashSet<usize> = principal_atoms.iter().copied().collect();
    let excluded: HashSet<usize> = chain_set.union(&principal_set).copied().collect();

    let mut consumed: HashSet<usize> = excluded.clone();
    let mut chain_subs: Vec<(usize, String)> = Vec::new();
    for (i, &c) in chain.iter().enumerate() {
        let locant = i + 1;
        for &nb in &g.adjacency[c] {
            if g.atoms[nb].symbol == "H" || chain_set.contains(&nb) || principal_set.contains(&nb) {
                continue;
            }
            if g.atoms[nb].symbol == "C" && get_bond_order(g, c, nb) == 2.0 {
                return None; // exo C=C (alkylidene) 未対応
            }
            let name = name_substituent(g, nb, &excluded, &mut consumed)?;
            chain_subs.push((locant, name));
        }
    }

    // アミン主基の N 上の置換基 (N-methyl 等)
    let n_atoms: Vec<usize> = principal_atoms
        .iter()
        .copied()
        .filter(|&a| g.atoms[a].symbol == "N")
        .collect();
    let mut n_subs: Vec<(usize, String)> = Vec::new();
    for (ni, &n) in n_atoms.iter().enumerate() {
        for &nb in &g.adjacency[n] {
            if g.atoms[nb].symbol == "H" || principal_set.contains(&nb) || chain_set.contains(&nb) {
                continue;
            }
            let name = name_substituent(g, nb, &excluded, &mut consumed)?;
            n_subs.push((ni, name));
        }
    }

    // 完全性チェック
    for i in 0..g.atoms.len() {
        if g.atoms[i].symbol != "H" && !consumed.contains(&i) {
            return None;
        }
    }
    chain_subs.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    Some((chain_subs, n_subs))
}

/// root から始まる置換基名。消費原子を `consumed` に加える。命名不能なら None。
pub fn name_substituent(
    g: &MoleculeGraph,
    root: usize,
    excluded: &HashSet<usize>,
    consumed: &mut HashSet<usize>,
) -> Option<String> {
    let sym = g.atoms[root].symbol.as_str();
    consumed.insert(root);

    if let Some(h) = halogen_prefix(sym) {
        return Some(h.to_string());
    }

    match sym {
        "O" => name_oxygen(g, root, excluded, consumed),
        "S" => name_sulfur(g, root, excluded, consumed),
        "N" => name_nitrogen(g, root, excluded, consumed),
        "C" => name_carbon(g, root, excluded, consumed),
        _ => None,
    }
}

fn has_h(g: &MoleculeGraph, i: usize) -> bool {
    g.adjacency[i].iter().any(|&nb| g.atoms[nb].symbol == "H")
}

fn name_oxygen(
    g: &MoleculeGraph,
    root: usize,
    excluded: &HashSet<usize>,
    consumed: &mut HashSet<usize>,
) -> Option<String> {
    // ペルオキシド (O-O) は未対応
    if g.adjacency[root].iter().any(|&nb| g.atoms[nb].symbol == "O") {
        return None;
    }
    // =O → oxo
    if g.adjacency[root]
        .iter()
        .any(|&nb| get_bond_order(g, root, nb) == 2.0)
    {
        return Some("oxo".to_string());
    }
    if has_h(g, root) {
        return Some("hydroxy".to_string());
    }
    // エーテル -O-C → alkoxy
    let c_other: Vec<usize> = g.adjacency[root]
        .iter()
        .copied()
        .filter(|&nb| !excluded.contains(&nb) && g.atoms[nb].symbol == "C")
        .collect();
    if let Some(&c) = c_other.first() {
        let mut ex = excluded.clone();
        ex.insert(root);
        // アシルオキシ (相手が C=O) は未対応
        let is_acyl = g.adjacency[c]
            .iter()
            .any(|&o| o != root && g.atoms[o].symbol == "O" && get_bond_order(g, c, o) == 2.0);
        if is_acyl {
            return None;
        }
        let alkyl = name_carbon(g, c, &ex, consumed)?;
        return Some(make_oxy_name(&alkyl));
    }
    Some("oxy".to_string())
}

/// alkyl + "oxy"、母音重複や特例 (methyl→methoxy 等) を処理。
fn make_oxy_name(alkyl: &str) -> String {
    match alkyl {
        "methyl" => "methoxy".to_string(),
        "ethyl" => "ethoxy".to_string(),
        "propyl" => "propoxy".to_string(),
        "butyl" => "butoxy".to_string(),
        _ => {
            // 複合名は括弧
            if alkyl.contains('-') || alkyl.chars().any(|c| c.is_ascii_digit()) {
                format!("({alkyl}oxy)")
            } else {
                format!("{alkyl}oxy")
            }
        }
    }
}

fn name_sulfur(
    g: &MoleculeGraph,
    root: usize,
    excluded: &HashSet<usize>,
    consumed: &mut HashSet<usize>,
) -> Option<String> {
    if has_h(g, root) {
        return Some("sulfanyl".to_string());
    }
    let n_oxo = g.adjacency[root]
        .iter()
        .filter(|&&nb| g.atoms[nb].symbol == "O" && get_bond_order(g, root, nb) == 2.0)
        .count();
    let c_other: Vec<usize> = g.adjacency[root]
        .iter()
        .copied()
        .filter(|&nb| !excluded.contains(&nb) && g.atoms[nb].symbol == "C")
        .collect();
    if let Some(&c) = c_other.first() {
        if n_oxo != 0 {
            return None; // sulfinyl/sulfonyl は未対応
        }
        let mut ex = excluded.clone();
        ex.insert(root);
        let alkyl = name_carbon(g, c, &ex, consumed)?;
        let name = format!("{alkyl}sulfanyl");
        if alkyl.contains('-') || alkyl.chars().any(|c| c.is_ascii_digit()) {
            return Some(format!("({alkyl})sulfanyl"));
        }
        return Some(name);
    }
    match n_oxo {
        2 => Some("sulfonyl".to_string()),
        1 => Some("sulfinyl".to_string()),
        _ => Some("sulfanyl".to_string()),
    }
}

fn name_nitrogen(
    g: &MoleculeGraph,
    root: usize,
    excluded: &HashSet<usize>,
    consumed: &mut HashSet<usize>,
) -> Option<String> {
    // ニトロ: N(=O)=O or [N+](=O)[O-]
    let o_double = g.adjacency[root]
        .iter()
        .filter(|&&nb| g.atoms[nb].symbol == "O" && get_bond_order(g, root, nb) == 2.0)
        .count();
    let o_neg = g.adjacency[root]
        .iter()
        .filter(|&&nb| {
            g.atoms[nb].symbol == "O"
                && get_bond_order(g, root, nb) == 1.0
                && g.atoms[nb].formal_charge < 0
        })
        .count();
    if o_double + o_neg >= 2 {
        for &nb in &g.adjacency[root] {
            if g.atoms[nb].symbol == "O" {
                consumed.insert(nb);
            }
        }
        return Some("nitro".to_string());
    }
    // 二重/三重結合を持つ N (イミン/ジアゾ/アジド等) は未対応
    if g.adjacency[root]
        .iter()
        .any(|&nb| get_bond_order(g, root, nb) >= 2.0)
    {
        return None;
    }
    // 非炭素・非水素の重原子隣接 (N-N, N-O 等) は未対応
    if g.adjacency[root].iter().any(|&nb| {
        !excluded.contains(&nb) && !matches!(g.atoms[nb].symbol.as_str(), "H" | "C")
    }) {
        return None;
    }
    // アミノ: -NH2 (無置換のみ対応。N-置換アミノは未対応)
    let c_other: Vec<usize> = g.adjacency[root]
        .iter()
        .copied()
        .filter(|&nb| !excluded.contains(&nb) && g.atoms[nb].symbol == "C")
        .collect();
    if c_other.is_empty() {
        return Some("amino".to_string());
    }
    None
}

/// 炭素置換基 (yl 基)。付け根に最小ロカント。直鎖・分岐対応。
fn name_carbon(
    g: &MoleculeGraph,
    root: usize,
    excluded: &HashSet<usize>,
    consumed: &mut HashSet<usize>,
) -> Option<String> {
    // 芳香環・カルボニル・カルボキシ・シアノ等は未対応
    let a = &g.atoms[root];
    if a.in_ring {
        return None;
    }
    // C=O (アシル)・C≡N (シアノ)・カルボキシは未対応
    for &nb in &g.adjacency[root] {
        let s = g.atoms[nb].symbol.as_str();
        if s == "O" && get_bond_order(g, root, nb) == 2.0 {
            return None;
        }
        if s == "N" && get_bond_order(g, root, nb) == 3.0 {
            return None;
        }
    }

    // 置換基の炭素集合 (root から到達、excluded を通らない、炭素のみ)
    let carbons = collect_sub_carbons(g, root, excluded);
    for &c in &carbons {
        consumed.insert(c);
    }
    if carbons.is_empty() {
        return Some("methyl".to_string());
    }
    // 全て単結合炭素・非環・ヘテロ置換なしの純アルキルのみ対応
    let cset: HashSet<usize> = carbons.iter().copied().collect();
    for &c in &carbons {
        for &nb in &g.adjacency[c] {
            if g.atoms[nb].symbol == "H" {
                continue;
            }
            if cset.contains(&nb) {
                if get_bond_order(g, c, nb) != 1.0 {
                    return None; // 不飽和 yl は未対応
                }
                continue;
            }
            if excluded.contains(&nb) {
                continue; // 付け根の主鎖側
            }
            // 置換基上の更なるヘテロ/炭素分岐 (メチル以外) — 単純分岐のみ許可
            if g.atoms[nb].symbol != "C" {
                return None;
            }
        }
    }

    // 付け根を含む最長鎖を探し、付け根に最小ロカントを与える
    let (chain, attach_locant) = longest_chain_through_root(g, root, &cset)?;
    let n = chain.len();
    let stem = CHAIN_PREFIX.get(n)?;
    // 鎖外の炭素 = メチル分岐
    let chain_set: HashSet<usize> = chain.iter().copied().collect();
    let mut branches: Vec<usize> = Vec::new();
    for (i, &c) in chain.iter().enumerate() {
        for &nb in &g.adjacency[c] {
            if cset.contains(&nb) && !chain_set.contains(&nb) && g.atoms[nb].symbol == "C" {
                branches.push(i + 1);
            }
        }
    }
    // 分岐は全てメチルであることを既に保証 (それ以外は上で None)
    if branches.is_empty() {
        if attach_locant == 1 {
            // 単純直鎖: methyl, ethyl, propyl, butyl...
            return Some(format!("{stem}yl"));
        }
        // 付け根が内部: propan-2-yl 等
        return Some(format!("{stem}an-{attach_locant}-yl"));
    }
    // メチル分岐つき: {locs}-methyl... {stem}-{attach}-yl
    branches.sort_unstable();
    let mult = MULTIPLIER.get(branches.len()).copied().unwrap_or("");
    let loc_str = branches
        .iter()
        .map(|l| l.to_string())
        .collect::<Vec<_>>()
        .join(",");
    if attach_locant == 1 {
        // 付け根が末端: 2-methylpropyl 等 ("an-1" を省略)
        Some(format!("{loc_str}-{mult}methyl{stem}yl"))
    } else {
        Some(format!("{loc_str}-{mult}methyl{stem}an-{attach_locant}-yl"))
    }
}

/// root から到達可能な (excluded を通らない) 炭素集合。
fn collect_sub_carbons(g: &MoleculeGraph, root: usize, excluded: &HashSet<usize>) -> Vec<usize> {
    let mut seen = HashSet::new();
    let mut stack = vec![root];
    let mut out = Vec::new();
    while let Some(c) = stack.pop() {
        if seen.contains(&c) || g.atoms[c].symbol != "C" {
            continue;
        }
        seen.insert(c);
        out.push(c);
        for &nb in &g.adjacency[c] {
            if !excluded.contains(&nb) && !seen.contains(&nb) && g.atoms[nb].symbol == "C" {
                stack.push(nb);
            }
        }
    }
    out
}

/// 付け根 root を含む最長単純鎖と、その中での root の最小ロカントを返す。
fn longest_chain_through_root(
    g: &MoleculeGraph,
    root: usize,
    cset: &HashSet<usize>,
) -> Option<(Vec<usize>, usize)> {
    // cset 内の炭素隣接
    let adj = |c: usize| -> Vec<usize> {
        g.adjacency[c]
            .iter()
            .copied()
            .filter(|nb| cset.contains(nb))
            .collect()
    };
    // root から各末端への最長パスを全探索 (置換基は小さい前提)
    let mut best: Option<Vec<usize>> = None;
    let n = g.atoms.len();
    fn dfs(
        cur: usize,
        vis: &mut Vec<bool>,
        path: &mut Vec<usize>,
        adj: &dyn Fn(usize) -> Vec<usize>,
        best: &mut Option<Vec<usize>>,
    ) {
        path.push(cur);
        vis[cur] = true;
        let mut extended = false;
        for nb in adj(cur) {
            if !vis[nb] {
                extended = true;
                dfs(nb, vis, path, adj, best);
            }
        }
        if !extended && best.as_ref().map(|b| path.len() > b.len()).unwrap_or(true) {
            *best = Some(path.clone());
        }
        path.pop();
        vis[cur] = false;
    }
    // root を必ず含むよう、root を一端とするパスと、root を経由する 2 端結合を考える。
    // 簡略: root から両方向の最長を連結する。
    let neighbors = adj(root);
    if neighbors.len() <= 1 {
        // root 末端: root からの最長パス
        let mut vis = vec![false; n];
        let mut path = Vec::new();
        dfs(root, &mut vis, &mut path, &adj, &mut best);
        let chain = best?;
        return Some((chain, 1));
    }
    // root 内部: 2 本の最長腕を連結、root は接続点
    // 各腕 = root の隣接から先の最長パス (root を除く)
    let mut arms: Vec<Vec<usize>> = Vec::new();
    for &start in &neighbors {
        let mut vis = vec![false; n];
        vis[root] = true;
        let mut path = Vec::new();
        let mut b: Option<Vec<usize>> = None;
        dfs(start, &mut vis, &mut path, &adj, &mut b);
        if let Some(arm) = b {
            arms.push(arm);
        }
    }
    arms.sort_by_key(|a| std::cmp::Reverse(a.len()));
    if arms.len() < 2 {
        let mut chain = vec![root];
        if let Some(a) = arms.first() {
            chain.extend(a);
        }
        return Some((chain, 1));
    }
    // 最長 2 腕: arm0 (逆順) + root + arm1
    let mut chain: Vec<usize> = arms[0].iter().rev().copied().collect();
    chain.push(root);
    chain.extend(&arms[1]);
    // root のロカント = arm0 長 + 1。逆向きなら arm1 長 + 1。小さい方。
    let loc_a = arms[0].len() + 1;
    let loc_b = arms[1].len() + 1;
    if loc_b < loc_a {
        chain.reverse();
        Some((chain, loc_b))
    } else {
        Some((chain, loc_a))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use molrs::graph::build_molecule_graph;

    fn sub(smiles: &str, root: usize, ex: &[usize]) -> Option<String> {
        let g = build_molecule_graph(smiles).unwrap();
        let exs: HashSet<usize> = ex.iter().copied().collect();
        let mut consumed = HashSet::new();
        name_substituent(&g, root, &exs, &mut consumed)
    }

    #[test]
    fn simple_prefixes() {
        // CCO: O(idx2) は末端 hydroxy
        assert_eq!(sub("CO", 1, &[0]), Some("hydroxy".to_string()));
        // ClC: Cl(idx0)
        assert_eq!(sub("ClC", 0, &[1]), Some("chloro".to_string()));
    }

    #[test]
    fn alkyl_groups() {
        // 2-methylpropane の中心から見たメチル: CC(C)C, root=メチル
        assert_eq!(sub("CC(C)C", 0, &[1]), Some("methyl".to_string()));
    }
}
