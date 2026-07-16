//! 官能基の検出と優先順位付け (S3.1)。
//!
//! Python `functional_group.py` (コミット a01eccd) の移植。
//! IUPAC 2013 Blue Book P-65 の seniority order に従い、
//! 分子グラフから官能基を検出して優先順に返す。
//!
//! 検出順・隣接走査順・タイブレークまで Python と一致させる
//! (ゲート: `tests/functional_group_compat.rs` でコーパス全数一致)。

use std::collections::HashSet;

use molrs::graph::{get_bond_order, MoleculeGraph};

use crate::constants::functional_group_priority;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FunctionalGroup {
    /// 'carboxylic_acid', 'alcohol', ...
    pub group_type: &'static str,
    /// 官能基を構成する原子インデックス
    pub atom_indices: Vec<usize>,
    /// FUNCTIONAL_GROUP_PRIORITY の値
    pub priority: i32,
}

fn prio(name: &str) -> i32 {
    functional_group_priority(name)
        .unwrap_or_else(|| panic!("unknown functional group type: {name}"))
}

fn push(groups: &mut Vec<FunctionalGroup>, group_type: &'static str, atom_indices: Vec<usize>) {
    groups.push(FunctionalGroup {
        group_type,
        atom_indices,
        priority: prio(group_type),
    });
}

/// 指定シンボルの隣接原子一覧 (隣接順)。
fn nbrs(g: &MoleculeGraph, idx: usize, symbol: &str) -> Vec<usize> {
    g.adjacency[idx]
        .iter()
        .copied()
        .filter(|&nb| g.atoms[nb].symbol == symbol)
        .collect()
}

/// 指定シンボルに二重結合している最初の隣接原子。
fn double_bonded_to(g: &MoleculeGraph, idx: usize, symbol: &str) -> Option<usize> {
    g.adjacency[idx]
        .iter()
        .copied()
        .find(|&nb| g.atoms[nb].symbol == symbol && get_bond_order(g, idx, nb) == 2.0)
}

fn has_h_neighbor(g: &MoleculeGraph, idx: usize) -> bool {
    g.adjacency[idx].iter().any(|&nb| g.atoms[nb].symbol == "H")
}

fn is_halogen(symbol: &str) -> bool {
    matches!(symbol, "F" | "Cl" | "Br" | "I")
}

// ─── 内部ヘルパー述語 (Python の _is_* / _get_* に対応) ────────────────

fn is_carboxylic_acid(g: &MoleculeGraph, c_idx: usize) -> bool {
    let mut has_double_o = false;
    let mut has_single_oh = false;
    for &nb_idx in &g.adjacency[c_idx] {
        if g.atoms[nb_idx].symbol != "O" {
            continue;
        }
        let bo = get_bond_order(g, c_idx, nb_idx);
        if bo == 2.0 {
            has_double_o = true;
        } else if bo == 1.0 && has_h_neighbor(g, nb_idx) {
            has_single_oh = true;
        }
    }
    has_double_o && has_single_oh
}

fn is_aldehyde(g: &MoleculeGraph, c_idx: usize) -> bool {
    if double_bonded_to(g, c_idx, "O").is_none() {
        return false;
    }
    let has_h = has_h_neighbor(g, c_idx);
    // カルボン酸系の単結合 O が存在してはいけない
    let has_single_o = g.adjacency[c_idx].iter().any(|&nb| {
        g.atoms[nb].symbol == "O" && get_bond_order(g, c_idx, nb) == 1.0
    });
    has_h && !has_single_o
}

fn is_ketone(g: &MoleculeGraph, c_idx: usize) -> bool {
    if double_bonded_to(g, c_idx, "O").is_none() {
        return false;
    }
    if is_carboxylic_acid(g, c_idx) || is_aldehyde(g, c_idx) {
        return false;
    }
    let c_neighbors = nbrs(g, c_idx, "C");
    if c_neighbors.len() >= 2 {
        return true;
    }
    // ketene 型 C=C=O (累積系)
    if c_neighbors.len() == 1 && get_bond_order(g, c_idx, c_neighbors[0]) == 2.0 {
        return true;
    }
    // チオラクトン: S を挟んだ環内 C=O
    if c_neighbors.len() == 1 && g.atoms[c_idx].in_ring {
        let s_nbs = nbrs(g, c_idx, "S");
        if !s_nbs.is_empty()
            && g.ring_atom_sets
                .iter()
                .any(|rt| rt.contains(&c_idx) && rt.contains(&s_nbs[0]))
        {
            return true;
        }
    }
    false
}

fn get_double_bonded_oxygen(g: &MoleculeGraph, c_idx: usize) -> Option<usize> {
    double_bonded_to(g, c_idx, "O")
}

fn get_carbonyl_oxygens(g: &MoleculeGraph, c_idx: usize) -> Vec<usize> {
    nbrs(g, c_idx, "O")
}

fn is_carboxylate(g: &MoleculeGraph, c_idx: usize) -> bool {
    let mut has_double_o = false;
    let mut has_neg_o = false;
    let mut has_single_oh = false;
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        if nb.symbol != "O" {
            continue;
        }
        let bo = get_bond_order(g, c_idx, nb_idx);
        if bo == 2.0 {
            has_double_o = true;
        } else if bo == 1.0 && nb.formal_charge == -1 {
            has_neg_o = true;
        } else if bo == 1.0 && nb.num_hs >= 1 {
            has_single_oh = true;
        }
    }
    has_double_o && has_neg_o && !has_single_oh
}

fn is_nitrile(g: &MoleculeGraph, c_idx: usize) -> bool {
    if g.atoms[c_idx].in_ring {
        return false;
    }
    for &nb_idx in &g.adjacency[c_idx] {
        if g.atoms[nb_idx].symbol == "N" && get_bond_order(g, c_idx, nb_idx) == 3.0 {
            let n_heavy = g.adjacency[nb_idx]
                .iter()
                .any(|&n| n != c_idx && g.atoms[n].symbol != "H");
            if !n_heavy {
                return true;
            }
        }
    }
    false
}

fn get_triple_bonded_nitrogen(g: &MoleculeGraph, c_idx: usize) -> Option<usize> {
    g.adjacency[c_idx]
        .iter()
        .copied()
        .find(|&nb| g.atoms[nb].symbol == "N" && get_bond_order(g, c_idx, nb) == 3.0)
}

fn has_double_bonded_oxygen(g: &MoleculeGraph, c_idx: usize) -> bool {
    get_double_bonded_oxygen(g, c_idx).is_some()
}

fn is_imidic_acid(g: &MoleculeGraph, c_idx: usize) -> bool {
    let mut has_imine_n = false;
    let mut has_oh = false;
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        if nb.symbol == "N" && get_bond_order(g, c_idx, nb_idx) == 2.0 {
            has_imine_n = true;
        } else if nb.symbol == "O"
            && get_bond_order(g, c_idx, nb_idx) == 1.0
            && has_h_neighbor(g, nb_idx)
        {
            has_oh = true;
        }
    }
    has_imine_n && has_oh
}

fn is_imidate_ester(g: &MoleculeGraph, c_idx: usize) -> bool {
    let mut has_imine_n = false;
    let mut has_ester_o = false;
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        if nb.symbol == "N" && get_bond_order(g, c_idx, nb_idx) == 2.0 {
            has_imine_n = true;
        } else if nb.symbol == "O" && get_bond_order(g, c_idx, nb_idx) == 1.0 {
            let o_c_nb = g.adjacency[nb_idx]
                .iter()
                .any(|&n| n != c_idx && g.atoms[n].symbol == "C");
            if o_c_nb && !has_h_neighbor(g, nb_idx) {
                has_ester_o = true;
            }
        }
    }
    has_imine_n && has_ester_o
}

fn get_imidate_atoms(g: &MoleculeGraph, c_idx: usize) -> (Option<usize>, Option<usize>) {
    let mut n_idx = None;
    let mut o_idx = None;
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        if nb.symbol == "N" && get_bond_order(g, c_idx, nb_idx) == 2.0 {
            n_idx = Some(nb_idx);
        } else if nb.symbol == "O" && get_bond_order(g, c_idx, nb_idx) == 1.0 {
            let o_c_nb = g.adjacency[nb_idx]
                .iter()
                .any(|&n| n != c_idx && g.atoms[n].symbol == "C");
            if o_c_nb && !has_h_neighbor(g, nb_idx) {
                o_idx = Some(nb_idx);
            }
        }
    }
    (n_idx, o_idx)
}

/// C=N の C 側に N 以外への二重結合があるか (イソシアネート/カルボジイミド除外)。
fn c_has_other_double(g: &MoleculeGraph, c_idx: usize, n_idx: usize) -> bool {
    g.adjacency[c_idx]
        .iter()
        .any(|&nb2| nb2 != n_idx && get_bond_order(g, c_idx, nb2) == 2.0)
}

fn is_imine(g: &MoleculeGraph, c_idx: usize) -> bool {
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        if nb.symbol == "N" && get_bond_order(g, c_idx, nb_idx) == 2.0 {
            if has_h_neighbor(g, nb_idx) {
                return true;
            }
            // N-置換イミン: N に C のみが付く
            let n_heavy: Vec<usize> = g.adjacency[nb_idx]
                .iter()
                .copied()
                .filter(|&n| n != c_idx && g.atoms[n].symbol != "H")
                .collect();
            let n_c_only = n_heavy.iter().all(|&n| g.atoms[n].symbol == "C");
            if !n_heavy.is_empty() && n_c_only && !c_has_other_double(g, c_idx, nb_idx) {
                return true;
            }
        }
    }
    false
}

fn get_imine_nitrogen(g: &MoleculeGraph, c_idx: usize) -> Option<usize> {
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        if nb.symbol == "N" && get_bond_order(g, c_idx, nb_idx) == 2.0 {
            if has_h_neighbor(g, nb_idx) {
                return Some(nb_idx);
            }
            let n_heavy: Vec<usize> = g.adjacency[nb_idx]
                .iter()
                .copied()
                .filter(|&n| n != c_idx && g.atoms[n].symbol != "H")
                .collect();
            if !n_heavy.is_empty()
                && n_heavy.iter().all(|&n| g.atoms[n].symbol == "C")
                && !c_has_other_double(g, c_idx, nb_idx)
            {
                return Some(nb_idx);
            }
        }
    }
    None
}

/// C=N-OH の N を見つけて C 隣接数を返す (オキシム判定の共通部)。
/// `require_exo_n`: N が環外であることを要求 (ketoxime)。
fn oxime_c_count(g: &MoleculeGraph, c_idx: usize, require_exo_n: bool) -> Option<usize> {
    for &nb_idx in &g.adjacency[c_idx] {
        if g.atoms[nb_idx].symbol != "N" || get_bond_order(g, c_idx, nb_idx) != 2.0 {
            continue;
        }
        if require_exo_n && g.atoms[nb_idx].in_ring {
            continue;
        }
        for &n_nb_idx in &g.adjacency[nb_idx] {
            if g.atoms[n_nb_idx].symbol == "O"
                && get_bond_order(g, nb_idx, n_nb_idx) == 1.0
                && has_h_neighbor(g, n_nb_idx)
            {
                let c_count = g.adjacency[c_idx]
                    .iter()
                    .filter(|&&nb2| g.atoms[nb2].symbol == "C")
                    .count();
                return Some(c_count);
            }
        }
    }
    None
}

fn is_ketoxime(g: &MoleculeGraph, c_idx: usize) -> bool {
    oxime_c_count(g, c_idx, true).is_some_and(|c| c >= 2)
}

fn is_aldoxime(g: &MoleculeGraph, c_idx: usize) -> bool {
    if g.atoms[c_idx].in_ring {
        return false;
    }
    oxime_c_count(g, c_idx, false).is_some_and(|c| c < 2)
}

fn get_oxime_nitrogen(g: &MoleculeGraph, c_idx: usize) -> Option<usize> {
    for &nb_idx in &g.adjacency[c_idx] {
        if g.atoms[nb_idx].symbol != "N" || get_bond_order(g, c_idx, nb_idx) != 2.0 {
            continue;
        }
        for &n_nb_idx in &g.adjacency[nb_idx] {
            if g.atoms[n_nb_idx].symbol == "O"
                && get_bond_order(g, nb_idx, n_nb_idx) == 1.0
                && has_h_neighbor(g, n_nb_idx)
            {
                return Some(nb_idx);
            }
        }
    }
    None
}

fn is_carbonate(g: &MoleculeGraph, c_idx: usize) -> bool {
    if g.atoms[c_idx].in_ring {
        return false;
    }
    if get_double_bonded_oxygen(g, c_idx).is_none() {
        return false;
    }
    if !nbrs(g, c_idx, "C").is_empty() {
        return false;
    }
    let mut single_o_with_c = 0usize;
    let mut single_o_oh = 0usize;
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        if nb.symbol != "O" || get_bond_order(g, c_idx, nb_idx) != 1.0 {
            continue;
        }
        if g.adjacency[nb_idx]
            .iter()
            .any(|&on| on != c_idx && g.atoms[on].symbol == "C")
        {
            single_o_with_c += 1;
        } else if has_h_neighbor(g, nb_idx) {
            single_o_oh += 1;
        }
    }
    single_o_with_c == 2 || (single_o_with_c == 1 && single_o_oh == 1)
}

fn is_anhydride(g: &MoleculeGraph, c_idx: usize) -> bool {
    if g.atoms[c_idx].in_ring {
        return false;
    }
    if get_double_bonded_oxygen(g, c_idx).is_none() {
        return false;
    }
    for &nb_idx in &g.adjacency[c_idx] {
        if g.atoms[nb_idx].symbol != "O" || get_bond_order(g, c_idx, nb_idx) != 1.0 {
            continue;
        }
        for &o_nb_idx in &g.adjacency[nb_idx] {
            if o_nb_idx == c_idx {
                continue;
            }
            if g.atoms[o_nb_idx].symbol == "C"
                && get_double_bonded_oxygen(g, o_nb_idx).is_some()
            {
                return true;
            }
        }
    }
    false
}

fn is_ester(g: &MoleculeGraph, c_idx: usize) -> bool {
    let mut has_double_o = false;
    let mut has_single_o_c = false;
    for &nb_idx in &g.adjacency[c_idx] {
        if g.atoms[nb_idx].symbol != "O" {
            continue;
        }
        let bo = get_bond_order(g, c_idx, nb_idx);
        if bo == 2.0 {
            has_double_o = true;
        } else if bo == 1.0
            && g.adjacency[nb_idx]
                .iter()
                .any(|&o_nb| g.atoms[o_nb].symbol == "C" && o_nb != c_idx)
        {
            has_single_o_c = true;
        }
    }
    has_double_o && has_single_o_c
}

fn get_double_bonded_sulfur(g: &MoleculeGraph, c_idx: usize) -> Option<usize> {
    double_bonded_to(g, c_idx, "S")
}

fn has_double_bonded_sulfur(g: &MoleculeGraph, c_idx: usize) -> bool {
    get_double_bonded_sulfur(g, c_idx).is_some()
}

fn is_o_thiocarbamate(g: &MoleculeGraph, c_idx: usize) -> bool {
    if get_double_bonded_sulfur(g, c_idx).is_none() {
        return false;
    }
    let mut has_n = false;
    let mut has_o_c = false;
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        if nb.symbol == "N" && !nb.in_ring && get_bond_order(g, c_idx, nb_idx) == 1.0 {
            has_n = true;
        }
        if nb.symbol == "O"
            && get_bond_order(g, c_idx, nb_idx) == 1.0
            && g.adjacency[nb_idx]
                .iter()
                .any(|&x| x != c_idx && g.atoms[x].symbol == "C")
        {
            has_o_c = true;
        }
    }
    has_n && has_o_c
}

fn is_s_carbamothioate(g: &MoleculeGraph, c_idx: usize) -> bool {
    if get_double_bonded_oxygen(g, c_idx).is_none() {
        return false;
    }
    let mut has_n = false;
    let mut has_s_c = false;
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        if nb.symbol == "N" && !nb.in_ring && get_bond_order(g, c_idx, nb_idx) == 1.0 {
            has_n = true;
        }
        if nb.symbol == "S"
            && get_bond_order(g, c_idx, nb_idx) == 1.0
            && g.adjacency[nb_idx]
                .iter()
                .any(|&x| x != c_idx && g.atoms[x].symbol == "C")
        {
            has_s_c = true;
        }
    }
    has_n && has_s_c
}

fn is_s_carbamodithioate(g: &MoleculeGraph, c_idx: usize) -> bool {
    let Some(s_double) = get_double_bonded_sulfur(g, c_idx) else {
        return false;
    };
    let mut has_n = false;
    let mut has_s_c = false;
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        if nb.symbol == "N" && !nb.in_ring && get_bond_order(g, c_idx, nb_idx) == 1.0 {
            has_n = true;
        }
        if nb.symbol == "S"
            && nb_idx != s_double
            && get_bond_order(g, c_idx, nb_idx) == 1.0
            && g.adjacency[nb_idx]
                .iter()
                .any(|&x| x != c_idx && g.atoms[x].symbol == "C")
        {
            has_s_c = true;
        }
    }
    has_n && has_s_c
}

fn is_thioamide(g: &MoleculeGraph, c_idx: usize) -> bool {
    if get_double_bonded_sulfur(g, c_idx).is_none() {
        return false;
    }
    let c_in_ring = g.atoms[c_idx].in_ring;
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        if nb.symbol == "N" {
            if c_in_ring && nb.in_ring {
                continue; // thiolactam は環命名系が扱う
            }
            if get_bond_order(g, c_idx, nb_idx) == 1.0 {
                return true;
            }
        }
    }
    false
}

fn get_thioamide_nitrogen(g: &MoleculeGraph, c_idx: usize) -> Option<usize> {
    let c_in_ring = g.atoms[c_idx].in_ring;
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        if nb.symbol == "N" {
            if c_in_ring && nb.in_ring {
                continue;
            }
            if get_bond_order(g, c_idx, nb_idx) == 1.0 {
                return Some(nb_idx);
            }
        }
    }
    None
}

fn get_double_bonded_selenium(g: &MoleculeGraph, c_idx: usize) -> Option<usize> {
    double_bonded_to(g, c_idx, "Se")
}

fn get_double_bonded_tellurium(g: &MoleculeGraph, c_idx: usize) -> Option<usize> {
    double_bonded_to(g, c_idx, "Te")
}

/// C(=X)-N (環外 N, 単結合) パターン (seleno/telluramide 共通)。
fn has_exo_single_n(g: &MoleculeGraph, c_idx: usize) -> bool {
    g.adjacency[c_idx].iter().any(|&nb_idx| {
        let nb = &g.atoms[nb_idx];
        nb.symbol == "N" && !nb.in_ring && get_bond_order(g, c_idx, nb_idx) == 1.0
    })
}

fn is_selenoamide(g: &MoleculeGraph, c_idx: usize) -> bool {
    get_double_bonded_selenium(g, c_idx).is_some() && has_exo_single_n(g, c_idx)
}

fn is_telluramide(g: &MoleculeGraph, c_idx: usize) -> bool {
    get_double_bonded_tellurium(g, c_idx).is_some() && has_exo_single_n(g, c_idx)
}

fn is_carbamate(g: &MoleculeGraph, c_idx: usize) -> bool {
    let mut has_double_o = false;
    let mut has_single_o_c = false;
    let mut has_single_n = false;
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        if nb.symbol == "O" {
            let bo = get_bond_order(g, c_idx, nb_idx);
            if bo == 2.0 {
                has_double_o = true;
            } else if bo == 1.0
                && g.adjacency[nb_idx]
                    .iter()
                    .any(|&o_nb| g.atoms[o_nb].symbol == "C" && o_nb != c_idx)
            {
                has_single_o_c = true;
            }
        } else if nb.symbol == "N" && get_bond_order(g, c_idx, nb_idx) == 1.0 {
            has_single_n = true;
        }
    }
    has_double_o && has_single_o_c && has_single_n
}

fn get_imine_double_bonded_nitrogen(g: &MoleculeGraph, c_idx: usize) -> Option<usize> {
    double_bonded_to(g, c_idx, "N")
}

fn is_kethydrazone(g: &MoleculeGraph, c_idx: usize) -> bool {
    let Some(n1_idx) = get_imine_double_bonded_nitrogen(g, c_idx) else {
        return false;
    };
    if g.atoms[c_idx].in_ring && g.atoms[n1_idx].in_ring {
        return false;
    }
    for &nb in &g.adjacency[n1_idx] {
        if nb == c_idx {
            continue;
        }
        if g.atoms[nb].symbol == "N" {
            if get_bond_order(g, n1_idx, nb) != 1.0 {
                continue;
            }
            let c_nbrs = nbrs(g, c_idx, "C").len();
            let h_nbrs = nbrs(g, c_idx, "H").len();
            return c_nbrs >= 2 || (c_nbrs == 1 && h_nbrs == 0);
        }
    }
    false
}

fn is_aldhydrazone(g: &MoleculeGraph, c_idx: usize) -> bool {
    let Some(n1_idx) = get_imine_double_bonded_nitrogen(g, c_idx) else {
        return false;
    };
    if g.atoms[c_idx].in_ring && g.atoms[n1_idx].in_ring {
        return false;
    }
    for &nb in &g.adjacency[n1_idx] {
        if nb == c_idx {
            continue;
        }
        if g.atoms[nb].symbol == "N" {
            if get_bond_order(g, n1_idx, nb) != 1.0 {
                continue;
            }
            let h_nbrs = nbrs(g, c_idx, "H").len();
            let c_nbrs = nbrs(g, c_idx, "C").len();
            return h_nbrs >= 1 && c_nbrs <= 1;
        }
    }
    false
}

fn get_hydrazone_nitrogen(g: &MoleculeGraph, c_idx: usize) -> Option<usize> {
    get_imine_double_bonded_nitrogen(g, c_idx)
}

fn is_chloroformate(g: &MoleculeGraph, c_idx: usize) -> bool {
    if get_double_bonded_oxygen(g, c_idx).is_none() {
        return false;
    }
    let mut has_halide = false;
    let mut has_ester_o = false;
    let mut has_c_neighbor = false;
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        if is_halogen(&nb.symbol) {
            has_halide = true;
        } else if nb.symbol == "O" && get_bond_order(g, c_idx, nb_idx) == 1.0 {
            if g.adjacency[nb_idx]
                .iter()
                .any(|&n| n != c_idx && g.atoms[n].symbol == "C")
            {
                has_ester_o = true;
            }
        } else if nb.symbol == "C" {
            has_c_neighbor = true;
        }
    }
    has_halide && has_ester_o && !has_c_neighbor
}

fn is_acid_halide(g: &MoleculeGraph, c_idx: usize) -> bool {
    if get_double_bonded_oxygen(g, c_idx).is_none() {
        return false;
    }
    g.adjacency[c_idx]
        .iter()
        .any(|&nb| is_halogen(&g.atoms[nb].symbol))
}

fn is_acyl_azide(g: &MoleculeGraph, c_idx: usize) -> bool {
    if get_double_bonded_oxygen(g, c_idx).is_none() {
        return false;
    }
    for &nb_idx in &g.adjacency[c_idx] {
        if g.atoms[nb_idx].symbol == "N" && get_bond_order(g, c_idx, nb_idx) == 1.0 {
            let n_has_n_double = g.adjacency[nb_idx].iter().any(|&n2| {
                n2 != c_idx
                    && g.atoms[n2].symbol == "N"
                    && get_bond_order(g, nb_idx, n2) == 2.0
            });
            if n_has_n_double {
                return true;
            }
        }
    }
    false
}

fn is_amide(g: &MoleculeGraph, c_idx: usize) -> bool {
    get_double_bonded_oxygen(g, c_idx).is_some() && has_exo_single_n(g, c_idx)
}

fn get_amide_nitrogen(g: &MoleculeGraph, c_idx: usize) -> Option<usize> {
    g.adjacency[c_idx].iter().copied().find(|&nb_idx| {
        let nb = &g.atoms[nb_idx];
        nb.symbol == "N" && !nb.in_ring && get_bond_order(g, c_idx, nb_idx) == 1.0
    })
}

fn is_carbamic_acid(g: &MoleculeGraph, c_idx: usize) -> bool {
    let mut has_imine_o = false;
    let mut has_n = false;
    let mut has_oh = false;
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        let bo = get_bond_order(g, c_idx, nb_idx);
        if nb.symbol == "O" && bo == 2.0 {
            has_imine_o = true;
        } else if nb.symbol == "O" && bo == 1.0 {
            if has_h_neighbor(g, nb_idx) {
                has_oh = true;
            }
        } else if nb.symbol == "N" && bo == 1.0 {
            has_n = true;
        }
    }
    has_imine_o && has_n && has_oh
}

fn get_carbamic_oh(g: &MoleculeGraph, c_idx: usize) -> Option<usize> {
    g.adjacency[c_idx].iter().copied().find(|&nb_idx| {
        g.atoms[nb_idx].symbol == "O"
            && get_bond_order(g, c_idx, nb_idx) == 1.0
            && has_h_neighbor(g, nb_idx)
    })
}

fn is_amidine(g: &MoleculeGraph, c_idx: usize) -> bool {
    let mut has_imine_n = false;
    let mut has_amine_n = false;
    for &nb_idx in &g.adjacency[c_idx] {
        if g.atoms[nb_idx].symbol != "N" {
            continue;
        }
        let bo = get_bond_order(g, c_idx, nb_idx);
        if bo == 2.0 {
            // =N が他に二重結合を持たないこと (除: C=N=O, C=N=C)
            let has_other_db = g.adjacency[nb_idx]
                .iter()
                .any(|&n2| n2 != c_idx && get_bond_order(g, nb_idx, n2) == 2.0);
            if !has_other_db {
                has_imine_n = true;
            }
        } else if bo == 1.0 {
            // -N が二重結合を持たないこと
            let has_other_db = g.adjacency[nb_idx]
                .iter()
                .any(|&n2| n2 != c_idx && get_bond_order(g, nb_idx, n2) >= 2.0);
            if !has_other_db {
                has_amine_n = true;
            }
        }
    }
    has_imine_n && has_amine_n
}

fn get_amidine_nitrogens(g: &MoleculeGraph, c_idx: usize) -> (Option<usize>, Option<usize>) {
    let mut n_imine = None;
    let mut n_amine = None;
    for &nb_idx in &g.adjacency[c_idx] {
        if g.atoms[nb_idx].symbol != "N" {
            continue;
        }
        let bo = get_bond_order(g, c_idx, nb_idx);
        if bo == 2.0 {
            n_imine = Some(nb_idx);
        } else if bo == 1.0 {
            n_amine = Some(nb_idx);
        }
    }
    (n_imine, n_amine)
}

/// 末端 C≡N + 単結合 X (cyanate: X=O, thiocyanate: X=S) の共通判定。
fn is_pseudohalide_ester(g: &MoleculeGraph, c_idx: usize, x_symbol: &str) -> bool {
    let mut has_triple_n = false;
    let mut has_single_x = false;
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        let bo = get_bond_order(g, c_idx, nb_idx);
        if nb.symbol == "N" && bo == 3.0 {
            let n_heavy = g.adjacency[nb_idx]
                .iter()
                .any(|&n| n != c_idx && g.atoms[n].symbol != "H");
            if !n_heavy {
                has_triple_n = true;
            }
        } else if nb.symbol == x_symbol && bo == 1.0 {
            has_single_x = true;
        }
    }
    has_triple_n && has_single_x
}

fn is_cyanate(g: &MoleculeGraph, c_idx: usize) -> bool {
    is_pseudohalide_ester(g, c_idx, "O")
}

fn get_cyanate_oxygen(g: &MoleculeGraph, c_idx: usize) -> Option<usize> {
    g.adjacency[c_idx]
        .iter()
        .copied()
        .find(|&nb| g.atoms[nb].symbol == "O" && get_bond_order(g, c_idx, nb) == 1.0)
}

fn is_thiocyanate(g: &MoleculeGraph, c_idx: usize) -> bool {
    is_pseudohalide_ester(g, c_idx, "S")
}

fn get_thiocyanate_sulfur(g: &MoleculeGraph, c_idx: usize) -> Option<usize> {
    g.adjacency[c_idx]
        .iter()
        .copied()
        .find(|&nb| g.atoms[nb].symbol == "S" && get_bond_order(g, c_idx, nb) == 1.0)
}

fn is_carbodiimide(g: &MoleculeGraph, c_idx: usize) -> bool {
    if g.atoms[c_idx].in_ring {
        return false;
    }
    let n_double: Vec<usize> = g.adjacency[c_idx]
        .iter()
        .copied()
        .filter(|&nb| g.atoms[nb].symbol == "N" && get_bond_order(g, c_idx, nb) == 2.0)
        .collect();
    if n_double.len() != 2 {
        return false;
    }
    for &n_idx in &n_double {
        let has_c_sgl = g.adjacency[n_idx].iter().any(|&nb| {
            nb != c_idx && g.atoms[nb].symbol == "C" && get_bond_order(g, n_idx, nb) == 1.0
        });
        if !has_c_sgl {
            return false;
        }
    }
    true
}

fn is_semicarbazone_or_thio(g: &MoleculeGraph, c_idx: usize) -> Option<&'static str> {
    let n1_idx = get_imine_double_bonded_nitrogen(g, c_idx)?;

    // N1 の隣の N2 (単結合)
    let n2_idx = g.adjacency[n1_idx].iter().copied().find(|&nb| {
        nb != c_idx && g.atoms[nb].symbol == "N" && get_bond_order(g, n1_idx, nb) == 1.0
    })?;

    // N2 の隣の C2 (単結合)
    let c2_idx = g.adjacency[n2_idx].iter().copied().find(|&nb| {
        nb != n1_idx && g.atoms[nb].symbol == "C" && get_bond_order(g, n2_idx, nb) == 1.0
    })?;

    let has_carbonyl = get_double_bonded_oxygen(g, c2_idx).is_some();
    let has_thio = get_double_bonded_sulfur(g, c2_idx).is_some();
    if !(has_carbonyl || has_thio) {
        return None;
    }
    let has_n_h = g.adjacency[c2_idx].iter().any(|&nb| {
        nb != n2_idx
            && g.atoms[nb].symbol == "N"
            && get_bond_order(g, c2_idx, nb) == 1.0
            && has_h_neighbor(g, nb)
    });
    if !has_n_h {
        return None;
    }

    let h_count = nbrs(g, c_idx, "H").len();
    let c_count = g.adjacency[c_idx]
        .iter()
        .filter(|&&nb| nb != n1_idx && g.atoms[nb].symbol == "C")
        .count();

    let aldehyde_like = h_count >= 1 && c_count <= 1;
    if has_carbonyl {
        Some(if aldehyde_like { "aldsemicarbazone" } else { "semicarbazone" })
    } else {
        Some(if aldehyde_like { "aldthiosemicarbazone" } else { "thiosemicarbazone" })
    }
}

fn is_peroxyacid(g: &MoleculeGraph, c_idx: usize) -> bool {
    if get_double_bonded_oxygen(g, c_idx).is_none() {
        return false;
    }
    for &nb_idx in &g.adjacency[c_idx] {
        if g.atoms[nb_idx].symbol != "O" || get_bond_order(g, c_idx, nb_idx) != 1.0 {
            continue;
        }
        let o2 = g.adjacency[nb_idx]
            .iter()
            .copied()
            .find(|&n| n != c_idx && g.atoms[n].symbol == "O");
        let Some(o2_idx) = o2 else { continue };
        if has_h_neighbor(g, o2_idx) {
            return true;
        }
    }
    false
}

fn is_peroxyester(g: &MoleculeGraph, c_idx: usize) -> bool {
    if get_double_bonded_oxygen(g, c_idx).is_none() {
        return false;
    }
    for &nb_idx in &g.adjacency[c_idx] {
        if g.atoms[nb_idx].symbol != "O" || get_bond_order(g, c_idx, nb_idx) != 1.0 {
            continue;
        }
        let o2 = g.adjacency[nb_idx]
            .iter()
            .copied()
            .find(|&n| n != c_idx && g.atoms[n].symbol == "O");
        let Some(o2_idx) = o2 else { continue };
        // O2 に C が隣接すること (アシルペルオキシドは除外)
        for &alkyl_c in &g.adjacency[o2_idx] {
            if alkyl_c == nb_idx || g.atoms[alkyl_c].symbol != "C" {
                continue;
            }
            if get_double_bonded_oxygen(g, alkyl_c).is_some() {
                continue;
            }
            return true;
        }
    }
    false
}

/// C(=X)-N1-N2 (単結合鎖) パターンのヒドラジド骨格判定。
/// `require_exo_n1`: N1 が環外であることを要求 (thio/seleno 版)。
fn hydrazide_like(g: &MoleculeGraph, c_idx: usize, require_exo_n1: bool) -> bool {
    for &nb_idx in &g.adjacency[c_idx] {
        let nb = &g.atoms[nb_idx];
        if nb.symbol != "N" || get_bond_order(g, c_idx, nb_idx) != 1.0 {
            continue;
        }
        if require_exo_n1 && nb.in_ring {
            continue;
        }
        for &n2_idx in &g.adjacency[nb_idx] {
            if n2_idx == c_idx || g.atoms[n2_idx].symbol != "N" {
                continue;
            }
            if get_bond_order(g, nb_idx, n2_idx) != 1.0 {
                continue;
            }
            // N2 が C への二重結合を持たないこと (semicarbazone 除外)
            let n2_dbl_c = g.adjacency[n2_idx].iter().any(|&nb3| {
                nb3 != nb_idx
                    && g.atoms[nb3].symbol == "C"
                    && get_bond_order(g, n2_idx, nb3) == 2.0
            });
            if n2_dbl_c {
                continue;
            }
            return true;
        }
    }
    false
}

fn is_hydrazide(g: &MoleculeGraph, c_idx: usize) -> bool {
    get_double_bonded_oxygen(g, c_idx).is_some() && hydrazide_like(g, c_idx, false)
}

fn is_thiohydrazide(g: &MoleculeGraph, c_idx: usize) -> bool {
    get_double_bonded_sulfur(g, c_idx).is_some() && hydrazide_like(g, c_idx, true)
}

fn is_selenohydrazide(g: &MoleculeGraph, c_idx: usize) -> bool {
    get_double_bonded_selenium(g, c_idx).is_some() && hydrazide_like(g, c_idx, true)
}

fn is_sulfonyl_azide(g: &MoleculeGraph, s_idx: usize, n1_idx: usize) -> bool {
    g.adjacency[n1_idx].iter().any(|&n2| {
        n2 != s_idx && g.atoms[n2].symbol == "N" && get_bond_order(g, n1_idx, n2) == 2.0
    })
}

fn is_sulfonohydrazide(g: &MoleculeGraph, s_idx: usize, n1_idx: usize) -> bool {
    g.adjacency[n1_idx].iter().any(|&n2| {
        n2 != s_idx && g.atoms[n2].symbol == "N" && get_bond_order(g, n1_idx, n2) == 1.0
    })
}

// ─── メイン検出ループ ────────────────────────────────────────────

/// MoleculeGraph を走査し、全官能基を検出して優先順位の高い順に返す。
/// 最初の要素が principal characteristic group になる。
pub fn detect_groups(g: &MoleculeGraph) -> Vec<FunctionalGroup> {
    let mut groups: Vec<FunctionalGroup> = Vec::new();

    // --- 炭素中心の官能基パターン検出 ---
    for (idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "C" {
            continue;
        }

        // 酸無水物: C(=O)-O-C(=O) — ester より先に判定
        if is_anhydride(g, idx) {
            let mut o_link: Option<usize> = None;
            let mut c2: Option<usize> = None;
            'outer: for &nb_idx in &g.adjacency[idx] {
                if g.atoms[nb_idx].symbol == "O" && get_bond_order(g, idx, nb_idx) == 1.0 {
                    for &o_nb_idx in &g.adjacency[nb_idx] {
                        if o_nb_idx == idx {
                            continue;
                        }
                        if g.atoms[o_nb_idx].symbol == "C"
                            && get_double_bonded_oxygen(g, o_nb_idx).is_some()
                        {
                            o_link = Some(nb_idx);
                            c2 = Some(o_nb_idx);
                            break;
                        }
                    }
                }
                if o_link.is_some() {
                    break 'outer;
                }
            }
            let mut indices = vec![idx];
            indices.extend(o_link);
            indices.extend(c2);
            push(&mut groups, "anhydride", indices);
            continue;
        }

        // 炭酸エステル: RO-C(=O)-OR — エステルより先に判定
        if is_carbonate(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_oxygen(g, idx));
            for &nb_idx in &g.adjacency[idx] {
                if g.atoms[nb_idx].symbol == "O" && get_bond_order(g, idx, nb_idx) == 1.0 {
                    indices.push(nb_idx);
                }
            }
            push(&mut groups, "carbonate", indices);
            continue;
        }

        // カルバメート: N-C(=O)-O-R — ester/amide より先に判定
        if is_carbamate(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_oxygen(g, idx));
            for &nb_idx in &g.adjacency[idx] {
                let sym = &g.atoms[nb_idx].symbol;
                if (sym == "O" || sym == "N") && get_bond_order(g, idx, nb_idx) == 1.0 {
                    indices.push(nb_idx);
                }
            }
            push(&mut groups, "carbamate", indices);
            continue;
        }

        // クロロホルメート: Cl-C(=O)-O-R — ester より先に判定
        if is_chloroformate(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_oxygen(g, idx));
            for &nb_idx in &g.adjacency[idx] {
                if g.atoms[nb_idx].symbol == "O" && get_bond_order(g, idx, nb_idx) == 1.0 {
                    indices.push(nb_idx);
                    break;
                }
            }
            push(&mut groups, "chloroformate", indices);
            continue;
        }

        // ペルオキシ酸: C(=O)-O-O-H — ester より先に判定
        if is_peroxyacid(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_oxygen(g, idx));
            for &nb_idx in &g.adjacency[idx] {
                if g.atoms[nb_idx].symbol == "O" && get_bond_order(g, idx, nb_idx) == 1.0 {
                    indices.push(nb_idx);
                    for &o2_nb in &g.adjacency[nb_idx] {
                        if o2_nb != idx && g.atoms[o2_nb].symbol == "O" {
                            indices.push(o2_nb);
                        }
                    }
                    break;
                }
            }
            push(&mut groups, "peroxyacid", indices);
            continue;
        }

        // ペルオキシエステル: C(=O)-O-O-C — ester より先に判定
        if is_peroxyester(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_oxygen(g, idx));
            for &nb_idx in &g.adjacency[idx] {
                if g.atoms[nb_idx].symbol == "O" && get_bond_order(g, idx, nb_idx) == 1.0 {
                    indices.push(nb_idx);
                    for &o2_nb in &g.adjacency[nb_idx] {
                        if o2_nb != idx && g.atoms[o2_nb].symbol == "O" {
                            indices.push(o2_nb);
                        }
                    }
                    break;
                }
            }
            push(&mut groups, "peroxy_ester", indices);
            continue;
        }

        // エステル: C(=O)-O-R
        if is_ester(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_oxygen(g, idx));
            for &nb_idx in &g.adjacency[idx] {
                if g.atoms[nb_idx].symbol == "O" && get_bond_order(g, idx, nb_idx) == 1.0 {
                    indices.push(nb_idx);
                    break;
                }
            }
            push(&mut groups, "ester", indices);
            continue;
        }

        // 酸ハライド: C(=O)-X
        if is_acid_halide(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_oxygen(g, idx));
            push(&mut groups, "acid_halide", indices);
            continue;
        }

        // カルバミン酸: RnN-C(=O)-OH
        if is_carbamic_acid(g, idx) {
            let n_idx_ca = g.adjacency[idx].iter().copied().find(|&nb| {
                g.atoms[nb].symbol == "N" && get_bond_order(g, idx, nb) == 1.0
            });
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_oxygen(g, idx));
            indices.extend(n_idx_ca);
            indices.extend(get_carbamic_oh(g, idx));
            push(&mut groups, "carbamic_acid", indices);
            continue;
        }

        // S-アルキルカルバモジチオアート: N-C(=S)-S-R — チオアミドより先
        if is_s_carbamodithioate(g, idx) {
            let s_double = get_double_bonded_sulfur(g, idx);
            let n_idx = get_thioamide_nitrogen(g, idx);
            let s_ester = g.adjacency[idx].iter().copied().find(|&nb| {
                g.atoms[nb].symbol == "S"
                    && Some(nb) != s_double
                    && get_bond_order(g, idx, nb) == 1.0
                    && g.adjacency[nb]
                        .iter()
                        .any(|&x| x != idx && g.atoms[x].symbol == "C")
            });
            let mut indices = vec![idx];
            indices.extend(s_double);
            indices.extend(n_idx);
            indices.extend(s_ester);
            push(&mut groups, "s_carbamodithioate", indices);
            continue;
        }

        // O-アルキルチオカルバメート: N-C(=S)-O-R — チオアミドより先
        if is_o_thiocarbamate(g, idx) {
            let s_idx = get_double_bonded_sulfur(g, idx);
            let n_idx = get_thioamide_nitrogen(g, idx);
            let o_idx = g.adjacency[idx].iter().copied().find(|&nb| {
                g.atoms[nb].symbol == "O"
                    && get_bond_order(g, idx, nb) == 1.0
                    && g.adjacency[nb]
                        .iter()
                        .any(|&x| x != idx && g.atoms[x].symbol == "C")
            });
            let mut indices = vec![idx];
            indices.extend(s_idx);
            indices.extend(n_idx);
            indices.extend(o_idx);
            push(&mut groups, "o_thiocarbamate", indices);
            continue;
        }

        // チオヒドラジド: C(=S)-NH-NH₂ — チオアミドより先に検出
        if is_thiohydrazide(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_sulfur(g, idx));
            indices.extend(get_thioamide_nitrogen(g, idx));
            push(&mut groups, "thiohydrazide", indices);
            continue;
        }

        // セレノヒドラジド: C(=[Se])-NH-NH₂ — セレノアミドより先に検出
        if is_selenohydrazide(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_selenium(g, idx));
            indices.extend(get_thioamide_nitrogen(g, idx));
            push(&mut groups, "selenohydrazide", indices);
            continue;
        }

        // チオアミド: C(=S)-NR₂ — アミン検出より先に処理
        if is_thioamide(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_sulfur(g, idx));
            indices.extend(get_thioamide_nitrogen(g, idx));
            push(&mut groups, "thioamide", indices);
            continue;
        }

        // セレノアミド: C(=[Se])-NR₂
        if is_selenoamide(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_selenium(g, idx));
            indices.extend(get_thioamide_nitrogen(g, idx));
            push(&mut groups, "selenoamide", indices);
            continue;
        }

        // テルラミド: C(=[Te])-NR₂
        if is_telluramide(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_tellurium(g, idx));
            indices.extend(get_thioamide_nitrogen(g, idx));
            push(&mut groups, "telluramide", indices);
            continue;
        }

        // チオカルボン酸 — チオアルデヒドより先に判定
        if !atom.in_ring {
            let o2_tc: Vec<usize> = g.adjacency[idx]
                .iter()
                .copied()
                .filter(|&nb| g.atoms[nb].symbol == "O" && get_bond_order(g, idx, nb) == 2.0)
                .collect();
            let s2_tc: Vec<usize> = g.adjacency[idx]
                .iter()
                .copied()
                .filter(|&nb| g.atoms[nb].symbol == "S" && get_bond_order(g, idx, nb) == 2.0)
                .collect();
            let s1_tc: Vec<usize> = g.adjacency[idx]
                .iter()
                .copied()
                .filter(|&nb| {
                    g.atoms[nb].symbol == "S"
                        && get_bond_order(g, idx, nb) == 1.0
                        && (g.atoms[nb].num_hs >= 1 || has_h_neighbor(g, nb))
                })
                .collect();
            let o1_tc: Vec<usize> = g.adjacency[idx]
                .iter()
                .copied()
                .filter(|&nb| {
                    g.atoms[nb].symbol == "O"
                        && get_bond_order(g, idx, nb) == 1.0
                        && (g.atoms[nb].num_hs >= 1 || has_h_neighbor(g, nb))
                })
                .collect();
            if !o2_tc.is_empty() && !s1_tc.is_empty() {
                // C(=O)-SH: thioic S-acid
                let mut indices = vec![idx];
                indices.extend(&o2_tc);
                indices.extend(&s1_tc);
                push(&mut groups, "thioic_s_acid", indices);
                continue;
            }
            if !s2_tc.is_empty() && !o1_tc.is_empty() {
                // C(=S)-OH: thioic O-acid
                let mut indices = vec![idx];
                indices.extend(&s2_tc);
                indices.extend(&o1_tc);
                push(&mut groups, "thioic_o_acid", indices);
                continue;
            }
            if !s2_tc.is_empty() && !s1_tc.is_empty() {
                // C(=S)-SH: dithioic acid
                let mut indices = vec![idx];
                indices.extend(&s2_tc);
                indices.extend(&s1_tc);
                push(&mut groups, "dithioic_acid", indices);
                continue;
            }
        }

        // チオアルデヒド / チオケトン: C=S (exocyclic S のみ)
        let s_idx_tk = get_double_bonded_sulfur(g, idx);
        if let Some(s_tk) = s_idx_tk {
            if !g.atoms[s_tk].in_ring {
                // N または O への二重結合がないこと (isothiocyanate 等を除外)
                let has_other_double = g.adjacency[idx].iter().any(|&nb| {
                    matches!(g.atoms[nb].symbol.as_str(), "N" | "O")
                        && nb != s_tk
                        && get_bond_order(g, idx, nb) == 2.0
                });
                if !has_other_double {
                    // O-チオエステル / S-ジチオエステル
                    let ether_os: Vec<usize> = g.adjacency[idx]
                        .iter()
                        .copied()
                        .filter(|&nb| {
                            g.atoms[nb].symbol == "O"
                                && nb != s_tk
                                && !has_h_neighbor(g, nb)
                                && g.adjacency[nb]
                                    .iter()
                                    .any(|&onc| onc != idx && g.atoms[onc].symbol == "C")
                        })
                        .collect();
                    let thioether_ss: Vec<usize> = g.adjacency[idx]
                        .iter()
                        .copied()
                        .filter(|&nb| {
                            g.atoms[nb].symbol == "S"
                                && nb != s_tk
                                && !has_h_neighbor(g, nb)
                                && g.adjacency[nb]
                                    .iter()
                                    .any(|&snc| snc != idx && g.atoms[snc].symbol == "C")
                        })
                        .collect();
                    // carbonothioate / carbonodithioate — 中心 C に直接 C なし
                    let has_c_nbr = g.adjacency[idx]
                        .iter()
                        .any(|&nb| g.atoms[nb].symbol == "C" && nb != s_tk);
                    if !has_c_nbr {
                        if ether_os.len() == 2 {
                            let mut indices = vec![idx, s_tk];
                            indices.extend(&ether_os);
                            push(&mut groups, "carbonothioate", indices);
                            continue;
                        }
                        if thioether_ss.len() == 2 {
                            let mut indices = vec![idx, s_tk];
                            indices.extend(&thioether_ss);
                            push(&mut groups, "carbonodithioate", indices);
                            continue;
                        }
                    }
                    if !ether_os.is_empty() {
                        let o_ester_idx = ether_os[0];
                        let alkyl_cs: Vec<usize> = g.adjacency[o_ester_idx]
                            .iter()
                            .copied()
                            .filter(|&nb| nb != idx && g.atoms[nb].symbol == "C")
                            .collect();
                        let mut indices = vec![o_ester_idx, idx];
                        indices.extend(&alkyl_cs);
                        push(&mut groups, "o_thioester", indices);
                        continue;
                    }
                    if !thioether_ss.is_empty() {
                        let s_ester_idx = thioether_ss[0];
                        let alkyl_cs: Vec<usize> = g.adjacency[s_ester_idx]
                            .iter()
                            .copied()
                            .filter(|&nb| nb != idx && g.atoms[nb].symbol == "C")
                            .collect();
                        let mut indices = vec![s_ester_idx, idx];
                        indices.extend(&alkyl_cs);
                        push(&mut groups, "s_dithioate_ester", indices);
                        continue;
                    }
                    let c_single = g.adjacency[idx]
                        .iter()
                        .filter(|&&nb| {
                            nb != s_tk
                                && g.atoms[nb].symbol == "C"
                                && matches!(get_bond_order(g, idx, nb), 1.0 | 1.5)
                        })
                        .count();
                    let h_on_c = has_h_neighbor(g, idx);
                    // thioaldehyde: H あり (R-CH=S); H なし + c_single==0 は thioketene 型
                    if c_single <= 1 && !atom.in_ring && h_on_c {
                        push(&mut groups, "thioaldehyde", vec![idx, s_tk]);
                        continue;
                    } else if c_single >= 2 || (!atom.in_ring && !h_on_c) {
                        push(&mut groups, "thioketone", vec![idx, s_tk]);
                        continue;
                    }
                }
            }
        }

        // ヒドラジド: C(=O)-NH-NH₂ — アミドより先に判定
        if is_hydrazide(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_oxygen(g, idx));
            indices.extend(get_amide_nitrogen(g, idx));
            push(&mut groups, "hydrazide", indices);
            continue;
        }

        // アシルアジド: C(=O)-N=[N+]=[N-] — アミドより先に検出
        if is_acyl_azide(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_oxygen(g, idx));
            push(&mut groups, "acyl_azide", indices);
            continue;
        }

        // S-アルキルカルバモチオアート: N-C(=O)-S-R — アミドより先
        if is_s_carbamothioate(g, idx) {
            let s_idx = g.adjacency[idx].iter().copied().find(|&nb| {
                g.atoms[nb].symbol == "S"
                    && get_bond_order(g, idx, nb) == 1.0
                    && g.adjacency[nb]
                        .iter()
                        .any(|&x| x != idx && g.atoms[x].symbol == "C")
            });
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_oxygen(g, idx));
            indices.extend(get_amide_nitrogen(g, idx));
            indices.extend(s_idx);
            push(&mut groups, "s_carbamothioate", indices);
            continue;
        }

        // アミド: C(=O)-NR₂
        if is_amide(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_oxygen(g, idx));
            indices.extend(get_amide_nitrogen(g, idx));
            push(&mut groups, "amide", indices);
            continue;
        }

        // カルボキシレートアニオン: C(=O)[O-]
        if is_carboxylate(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_carbonyl_oxygens(g, idx));
            push(&mut groups, "carboxylate", indices);
            continue;
        }

        // カルボン酸: C(=O)O-H
        if is_carboxylic_acid(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_carbonyl_oxygens(g, idx));
            push(&mut groups, "carboxylic_acid", indices);
            continue;
        }

        // シアン酸エステル: O-C≡N — ニトリルより先に判定
        if is_cyanate(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_triple_bonded_nitrogen(g, idx));
            indices.extend(get_cyanate_oxygen(g, idx));
            push(&mut groups, "cyanate", indices);
            continue;
        }

        // チオシアン酸エステル: S-C≡N — ニトリルより先に判定
        if is_thiocyanate(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_triple_bonded_nitrogen(g, idx));
            indices.extend(get_thiocyanate_sulfur(g, idx));
            push(&mut groups, "thiocyanate", indices);
            continue;
        }

        // ニトリル: C≡N (末端)
        if is_nitrile(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_triple_bonded_nitrogen(g, idx));
            push(&mut groups, "nitrile", indices);
            continue;
        }

        // アルデヒド: C(=O)H
        if is_aldehyde(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_oxygen(g, idx));
            push(&mut groups, "aldehyde", indices);
            continue;
        }

        // ケトン: C(=O)C
        if is_ketone(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_double_bonded_oxygen(g, idx));
            push(&mut groups, "ketone", indices);
            continue;
        }

        // カルボジイミド: R-N=C=N-R
        if is_carbodiimide(g, idx) {
            let mut indices = vec![idx];
            indices.extend(g.adjacency[idx].iter().copied().filter(|&nb| {
                g.atoms[nb].symbol == "N" && get_bond_order(g, idx, nb) == 2.0
            }));
            push(&mut groups, "carbodiimide", indices);
            continue;
        }

        // アミジン: C(=N-H)(N-H) — imine より先に検出
        if is_amidine(g, idx) {
            let (n_imine_idx, n_amine_idx) = get_amidine_nitrogens(g, idx);
            // 全原子が環内にある場合は環命名系に委譲
            let all_in_ring = g.atoms[idx].in_ring
                && n_imine_idx.is_none_or(|n| g.atoms[n].in_ring)
                && n_amine_idx.is_none_or(|n| g.atoms[n].in_ring);
            if !all_in_ring {
                let mut indices = vec![idx];
                indices.extend(n_imine_idx);
                indices.extend(n_amine_idx);
                push(&mut groups, "amidine", indices);
                continue;
            }
        }

        // イミド酸: C(=N)(O-H) — imidate_ester より先に検出
        if is_imidic_acid(g, idx) {
            let n_idx = g.adjacency[idx].iter().copied().find(|&nb| {
                g.atoms[nb].symbol == "N" && get_bond_order(g, idx, nb) == 2.0
            });
            let mut indices = vec![idx];
            indices.extend(n_idx);
            push(&mut groups, "imidic_acid", indices);
            continue;
        }

        // イミデートエステル: C(=N)(O-R) — imine より先に検出
        if is_imidate_ester(g, idx) {
            let (n_idx, o_idx) = get_imidate_atoms(g, idx);
            let mut indices = vec![idx];
            indices.extend(n_idx);
            indices.extend(o_idx);
            push(&mut groups, "imidate_ester", indices);
            continue;
        }

        // イミン: C=N-H / C=N-R
        if is_imine(g, idx) {
            let n_idx = get_imine_nitrogen(g, idx);
            // C=N 両原子が環内にある場合は環命名系に委譲
            let both_in_ring =
                g.atoms[idx].in_ring && n_idx.is_some_and(|n| g.atoms[n].in_ring);
            if !both_in_ring {
                let mut indices = vec![idx];
                indices.extend(n_idx);
                push(&mut groups, "imine", indices);
                continue;
            }
        }

        // オキシム: C=N-OH
        if is_ketoxime(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_oxime_nitrogen(g, idx));
            push(&mut groups, "ketoxime", indices);
            continue;
        }

        if is_aldoxime(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_oxime_nitrogen(g, idx));
            push(&mut groups, "aldoxime", indices);
            continue;
        }

        // セミカルバゾン — hydrazone より先に判定
        if let Some(sc) = is_semicarbazone_or_thio(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_hydrazone_nitrogen(g, idx));
            push(&mut groups, sc, indices);
            continue;
        }

        // ヒドラゾン: C=N-NH₂ / C=N-NHR
        if is_kethydrazone(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_hydrazone_nitrogen(g, idx));
            push(&mut groups, "kethydrazone", indices);
            continue;
        }

        if is_aldhydrazone(g, idx) {
            let mut indices = vec![idx];
            indices.extend(get_hydrazone_nitrogen(g, idx));
            push(&mut groups, "aldhydrazone", indices);
            continue;
        }
    }

    // アルコール: C-O-H (カルボニル系 C を除く)
    let carbonyl_carbons: HashSet<usize> =
        groups.iter().map(|fg| fg.atom_indices[0]).collect();
    for (o_idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "O" {
            continue;
        }
        let c_neighbors = nbrs(g, o_idx, "C");
        let h_neighbors = nbrs(g, o_idx, "H");
        if !c_neighbors.is_empty() && !h_neighbors.is_empty() {
            let c_idx = c_neighbors[0];
            if !carbonyl_carbons.contains(&c_idx) {
                push(&mut groups, "alcohol", vec![c_idx, o_idx]);
            }
        }
    }

    // ヒドロペルオキシド: C-O1-O2-H
    let mut seen_peroxy: HashSet<(usize, usize)> = HashSet::new();
    for (o1_idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "O" {
            continue;
        }
        let o1_c_nbrs = nbrs(g, o1_idx, "C");
        let o1_o_nbrs = nbrs(g, o1_idx, "O");
        if o1_c_nbrs.is_empty() || o1_o_nbrs.is_empty() {
            continue;
        }
        let o2_idx = o1_o_nbrs[0];
        if !has_h_neighbor(g, o2_idx) {
            continue;
        }
        // O1 に隣接する C が carbonyl C ならペルオキシ酸 → 除外
        let c_idx = o1_c_nbrs[0];
        if get_double_bonded_oxygen(g, c_idx).is_some() {
            continue;
        }
        let pair = (o1_idx.min(o2_idx), o1_idx.max(o2_idx));
        if !seen_peroxy.insert(pair) {
            continue;
        }
        push(&mut groups, "hydroperoxide", vec![c_idx, o1_idx, o2_idx]);
    }

    // 有機ペルオキシド: C-O-O-C (H なし)
    let mut seen_peroxide: HashSet<(usize, usize)> = HashSet::new();
    for (o1_idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "O" {
            continue;
        }
        let o1_c_nbrs = nbrs(g, o1_idx, "C");
        let o1_o_nbrs = nbrs(g, o1_idx, "O");
        if o1_c_nbrs.is_empty() || o1_o_nbrs.is_empty() {
            continue;
        }
        let o2_idx = o1_o_nbrs[0];
        let o2_c_nbrs: Vec<usize> = g.adjacency[o2_idx]
            .iter()
            .copied()
            .filter(|&nb| nb != o1_idx && g.atoms[nb].symbol == "C")
            .collect();
        if o2_c_nbrs.is_empty() {
            continue;
        }
        if has_h_neighbor(g, o1_idx) || has_h_neighbor(g, o2_idx) {
            continue;
        }
        let c1_idx = o1_c_nbrs[0];
        let c2_idx = o2_c_nbrs[0];
        if get_double_bonded_oxygen(g, c1_idx).is_some()
            || get_double_bonded_oxygen(g, c2_idx).is_some()
        {
            continue;
        }
        let pair = (o1_idx.min(o2_idx), o1_idx.max(o2_idx));
        if !seen_peroxide.insert(pair) {
            continue;
        }
        push(&mut groups, "peroxide", vec![c1_idx, o1_idx, o2_idx, c2_idx]);
    }

    // チオール / スルホキシド / スルホン / スルホンアミド: S 原子を走査
    for (s_idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "S" || atom.in_ring {
            continue;
        }
        let h_neighbors = nbrs(g, s_idx, "H");
        let c_neighbors = nbrs(g, s_idx, "C");
        let n_neighbors = nbrs(g, s_idx, "N");
        let o_double: Vec<usize> = g.adjacency[s_idx]
            .iter()
            .copied()
            .filter(|&nb| {
                g.atoms[nb].symbol == "O"
                    && (get_bond_order(g, s_idx, nb) == 2.0
                        || (atom.formal_charge == 1
                            && g.atoms[nb].formal_charge == -1
                            && get_bond_order(g, s_idx, nb) == 1.0))
            })
            .collect();
        let o_single: Vec<usize> = g.adjacency[s_idx]
            .iter()
            .copied()
            .filter(|&nb| {
                g.atoms[nb].symbol == "O"
                    && get_bond_order(g, s_idx, nb) == 1.0
                    && !(atom.formal_charge == 1 && g.atoms[nb].formal_charge == -1)
            })
            .collect();
        let o_single_oh: Vec<usize> = o_single
            .iter()
            .copied()
            .filter(|&nb| has_h_neighbor(g, nb))
            .collect();
        let halogen_neighbors: Vec<usize> = g.adjacency[s_idx]
            .iter()
            .copied()
            .filter(|&nb| is_halogen(&g.atoms[nb].symbol))
            .collect();

        // チェーン 1: スルホニルハライド / スルホン酸 / スルホン酸アニオン
        if o_double.len() == 2
            && c_neighbors.len() == 1
            && halogen_neighbors.len() == 1
            && n_neighbors.is_empty()
            && o_single_oh.is_empty()
        {
            let mut indices = vec![s_idx];
            indices.extend(&o_double);
            indices.extend(&c_neighbors);
            indices.extend(&halogen_neighbors);
            push(&mut groups, "sulfonyl_chloride", indices);
        } else if o_double.len() == 2
            && c_neighbors.len() == 1
            && !o_single_oh.is_empty()
            && n_neighbors.is_empty()
        {
            let mut indices = vec![s_idx];
            indices.extend(&o_double);
            indices.extend(&c_neighbors);
            indices.extend(&o_single_oh);
            push(&mut groups, "sulfonic_acid", indices);
        } else if o_double.len() == 2
            && c_neighbors.len() == 1
            && n_neighbors.is_empty()
            && o_single_oh.is_empty()
        {
            // スルホン酸アニオン: C-S(=O)₂-[O⁻]
            let o_single_neg: Vec<usize> = o_single
                .iter()
                .copied()
                .filter(|&nb| g.atoms[nb].formal_charge == -1)
                .collect();
            if !o_single_neg.is_empty() {
                let mut indices = vec![s_idx];
                indices.extend(&o_double);
                indices.extend(&c_neighbors);
                indices.extend(&o_single_neg);
                push(&mut groups, "sulfonate_anion", indices);
            }
        }

        // チェーン 2: スルフィニル系〜ポリスルフィド
        if o_double.len() == 1
            && c_neighbors.len() == 1
            && halogen_neighbors.len() == 1
            && n_neighbors.is_empty()
            && o_single_oh.is_empty()
        {
            // スルフィニルハライド: C-S(=O)-X
            let mut indices = vec![s_idx];
            indices.extend(&o_double);
            indices.extend(&c_neighbors);
            indices.extend(&halogen_neighbors);
            push(&mut groups, "sulfinyl_chloride", indices);
        } else if o_double.len() == 1
            && c_neighbors.len() == 1
            && !o_single_oh.is_empty()
            && n_neighbors.is_empty()
        {
            // スルフィン酸: C-S(=O)-OH
            let mut indices = vec![s_idx];
            indices.extend(&o_double);
            indices.extend(&c_neighbors);
            indices.extend(&o_single_oh);
            push(&mut groups, "sulfinic_acid", indices);
        } else if halogen_neighbors.len() == 1
            && c_neighbors.len() == 1
            && o_double.is_empty()
            && h_neighbors.is_empty()
            && o_single_oh.is_empty()
            && n_neighbors.is_empty()
        {
            // スルフェニルハライド: C-S-X
            let mut indices = vec![s_idx];
            indices.extend(&c_neighbors);
            indices.extend(&halogen_neighbors);
            push(&mut groups, "sulfenyl_halide", indices);
        } else if o_single_oh.len() == 1
            && c_neighbors.len() == 1
            && o_double.is_empty()
            && h_neighbors.is_empty()
        {
            // スルフェン酸: C-S-OH
            let mut indices = vec![s_idx];
            indices.extend(&c_neighbors);
            indices.extend(&o_single_oh);
            push(&mut groups, "sulfenic_acid", indices);
        } else if o_single.len() == 1
            && o_single_oh.is_empty()
            && c_neighbors.len() == 1
            && o_double.is_empty()
            && h_neighbors.is_empty()
            && halogen_neighbors.is_empty()
            && n_neighbors.is_empty()
        {
            // スルフェン酸エステル: C-S-O-C
            let o_ester_idx = o_single[0];
            let has_o_ester_c = g.adjacency[o_ester_idx]
                .iter()
                .any(|&nb| nb != s_idx && g.atoms[nb].symbol == "C");
            if has_o_ester_c {
                let mut indices = vec![s_idx];
                indices.extend(&c_neighbors);
                indices.push(o_ester_idx);
                push(&mut groups, "sulfenate_ester", indices);
            }
        } else if !h_neighbors.is_empty() && c_neighbors.len() == 1 && o_double.is_empty() {
            // チオール: C-SH
            push(&mut groups, "thiol", vec![c_neighbors[0], s_idx]);
        } else if o_double.len() == 2 && n_neighbors.len() == 2 && c_neighbors.is_empty() {
            // スルファミド: H2N-S(=O)₂-NH2
            let mut indices = vec![s_idx];
            indices.extend(&o_double);
            indices.extend(&n_neighbors);
            push(&mut groups, "sulfamide", indices);
        } else if o_double.len() == 2
            && n_neighbors.len() == 1
            && c_neighbors.is_empty()
            && !o_single_oh.is_empty()
        {
            // スルファミン酸: H2N-S(=O)₂-OH
            let mut indices = vec![s_idx];
            indices.extend(&o_double);
            indices.extend(&n_neighbors);
            indices.extend(&o_single_oh);
            push(&mut groups, "sulfamic_acid", indices);
        } else if o_double.len() == 2
            && n_neighbors.len() == 1
            && c_neighbors.len() == 1
            && is_sulfonyl_azide(g, s_idx, n_neighbors[0])
        {
            // スルホニルアジド: C-S(=O)₂-N=[N+]=[N-]
            let n_idx_sa = n_neighbors[0];
            let n2_nbrs_sa: Vec<usize> = g.adjacency[n_idx_sa]
                .iter()
                .copied()
                .filter(|&nb| {
                    nb != s_idx
                        && g.atoms[nb].symbol == "N"
                        && get_bond_order(g, n_idx_sa, nb) == 2.0
                })
                .collect();
            let mut n3_nbrs_sa: Vec<usize> = Vec::new();
            for &n2_sa in &n2_nbrs_sa {
                n3_nbrs_sa.extend(g.adjacency[n2_sa].iter().copied().filter(|&nb| {
                    nb != n_idx_sa && g.atoms[nb].symbol == "N"
                }));
            }
            let mut indices = vec![s_idx];
            indices.extend(&o_double);
            indices.extend(&c_neighbors);
            indices.push(n_idx_sa);
            indices.extend(&n2_nbrs_sa);
            indices.extend(&n3_nbrs_sa);
            push(&mut groups, "sulfonyl_azide", indices);
        } else if o_double.len() == 2
            && n_neighbors.len() == 1
            && c_neighbors.len() == 1
            && is_sulfonohydrazide(g, s_idx, n_neighbors[0])
        {
            // スルホノヒドラジド: C-S(=O)₂-NH-NH₂
            let n1_idx = n_neighbors[0];
            let n2_nbrs: Vec<usize> = g.adjacency[n1_idx]
                .iter()
                .copied()
                .filter(|&nb| {
                    nb != s_idx
                        && g.atoms[nb].symbol == "N"
                        && get_bond_order(g, n1_idx, nb) == 1.0
                })
                .collect();
            let mut indices = vec![s_idx];
            indices.extend(&o_double);
            indices.extend(&c_neighbors);
            indices.push(n1_idx);
            indices.extend(&n2_nbrs);
            push(&mut groups, "sulfonohydrazide", indices);
        } else if o_double.len() == 2 && !n_neighbors.is_empty() && c_neighbors.len() == 1 {
            // スルホンアミド: C-S(=O)₂-N
            let mut indices = vec![s_idx];
            indices.extend(&o_double);
            indices.extend(&c_neighbors);
            indices.extend(&n_neighbors);
            push(&mut groups, "sulfonamide", indices);
        } else if o_double.len() == 1
            && c_neighbors.len() == 1
            && n_neighbors.len() == 1
            && is_sulfonohydrazide(g, s_idx, n_neighbors[0])
        {
            // スルフィニルヒドラジド: C-S(=O)-NH-NH₂
            let n1_idx = n_neighbors[0];
            let n2_nbrs: Vec<usize> = g.adjacency[n1_idx]
                .iter()
                .copied()
                .filter(|&nb| {
                    nb != s_idx
                        && g.atoms[nb].symbol == "N"
                        && get_bond_order(g, n1_idx, nb) == 1.0
                })
                .collect();
            let mut indices = vec![s_idx];
            indices.extend(&o_double);
            indices.extend(&c_neighbors);
            indices.push(n1_idx);
            indices.extend(&n2_nbrs);
            push(&mut groups, "sulfinylhydrazide", indices);
        } else if o_double.len() == 1 && c_neighbors.len() == 1 && !n_neighbors.is_empty() {
            // スルフィナミド: C-S(=O)-N
            let mut indices = vec![s_idx];
            indices.extend(&o_double);
            indices.extend(&c_neighbors);
            indices.extend(&n_neighbors);
            push(&mut groups, "sulfinamide", indices);
        } else if o_double.len() == 2
            && c_neighbors.len() == 1
            && o_single_oh.is_empty()
            && n_neighbors.is_empty()
        {
            // スルホン酸エステル: C-S(=O)₂-O-C — sulfone より先に
            let o_ester: Vec<usize> = o_single
                .iter()
                .copied()
                .filter(|&nb| {
                    !o_single_oh.contains(&nb)
                        && g.adjacency[nb]
                            .iter()
                            .any(|&occ| occ != s_idx && g.atoms[occ].symbol == "C")
                })
                .collect();
            if !o_ester.is_empty() {
                let mut indices = vec![s_idx];
                indices.extend(&o_double);
                indices.extend(&c_neighbors);
                indices.push(o_ester[0]);
                push(&mut groups, "sulfonate_ester", indices);
            }
        } else if o_double.len() == 1
            && c_neighbors.len() == 1
            && o_single_oh.is_empty()
            && n_neighbors.is_empty()
        {
            // スルフィン酸エステル: C-S(=O)-O-C — sulfoxide より先に
            let o_ester: Vec<usize> = o_single
                .iter()
                .copied()
                .filter(|&nb| {
                    !o_single_oh.contains(&nb)
                        && g.adjacency[nb]
                            .iter()
                            .any(|&occ| occ != s_idx && g.atoms[occ].symbol == "C")
                })
                .collect();
            if !o_ester.is_empty() {
                let mut indices = vec![s_idx];
                indices.extend(&o_double);
                indices.extend(&c_neighbors);
                indices.push(o_ester[0]);
                push(&mut groups, "sulfinate_ester", indices);
            }
        } else if o_double.len() == 2 && c_neighbors.len() == 2 {
            // スルホン: C-S(=O)₂-C
            let mut indices = vec![s_idx];
            indices.extend(&o_double);
            indices.extend(&c_neighbors);
            push(&mut groups, "sulfone", indices);
        } else if o_double.len() == 1 && c_neighbors.len() == 2 {
            // スルホキシド: C-S(=O)-C
            let mut indices = vec![s_idx];
            indices.extend(&o_double);
            indices.extend(&c_neighbors);
            push(&mut groups, "sulfoxide", indices);
        } else if h_neighbors.is_empty()
            && c_neighbors.len() == 2
            && o_double.is_empty()
            && n_neighbors.is_empty()
        {
            // チオエステル vs チオエーテル
            let carbonyl_cs: Vec<usize> = c_neighbors
                .iter()
                .copied()
                .filter(|&c| has_double_bonded_oxygen(g, c))
                .collect();
            if !carbonyl_cs.is_empty() {
                // チオラクトン (環状チオエステル) は ring ketone として扱う
                let same_ring = g
                    .ring_atom_sets
                    .iter()
                    .any(|rt| rt.contains(&s_idx) && rt.contains(&carbonyl_cs[0]));
                if !same_ring {
                    let alkyl_cs: Vec<usize> = c_neighbors
                        .iter()
                        .copied()
                        .filter(|c| !carbonyl_cs.contains(c))
                        .collect();
                    let mut indices = vec![s_idx, carbonyl_cs[0]];
                    indices.extend(&alkyl_cs);
                    push(&mut groups, "thioester", indices);
                }
            } else {
                let mut indices = vec![s_idx];
                indices.extend(&c_neighbors);
                push(&mut groups, "sulfide", indices);
            }
        } else if h_neighbors.is_empty()
            && c_neighbors.len() == 1
            && o_double.is_empty()
            && n_neighbors.is_empty()
        {
            // ジスルフィド / トリスルフィド / テトラスルフィド: C-Sn-C
            let s_neighbors = nbrs(g, s_idx, "S");
            if s_neighbors.len() == 1 {
                let mut s_chain = vec![s_idx, s_neighbors[0]];
                loop {
                    let tail = *s_chain.last().unwrap();
                    let next_s: Vec<usize> = g.adjacency[tail]
                        .iter()
                        .copied()
                        .filter(|&nb| g.atoms[nb].symbol == "S" && !s_chain.contains(&nb))
                        .collect();
                    if next_s.len() == 1 {
                        s_chain.push(next_s[0]);
                    } else {
                        break;
                    }
                }
                let end_s = *s_chain.last().unwrap();
                let end_c = nbrs(g, end_s, "C");
                if !end_c.is_empty() && s_idx < end_s {
                    let gtype = match s_chain.len() {
                        3 => "trisulfide",
                        4 => "tetrasulfide",
                        _ => "disulfide",
                    };
                    let mut indices = s_chain.clone();
                    indices.extend(&c_neighbors);
                    indices.extend(&end_c);
                    push(&mut groups, gtype, indices);
                }
            }
        }
    }

    // セレノール / セレニド / テルロール / テルリド / セレン酸: Se/Te 原子を走査
    for (se_idx, atom) in g.atoms.iter().enumerate() {
        if !matches!(atom.symbol.as_str(), "Se" | "Te") || atom.in_ring {
            continue;
        }
        let h_neighbors = nbrs(g, se_idx, "H");
        let c_neighbors = nbrs(g, se_idx, "C");
        let se_neighbors = nbrs(g, se_idx, &atom.symbol);
        let o_neighbors = nbrs(g, se_idx, "O");
        let o_double: Vec<usize> = o_neighbors
            .iter()
            .copied()
            .filter(|&nb| get_bond_order(g, se_idx, nb) == 2.0)
            .collect();
        let o_single: Vec<usize> = o_neighbors
            .iter()
            .copied()
            .filter(|&nb| get_bond_order(g, se_idx, nb) == 1.0)
            .collect();
        let o_single_oh: Vec<usize> = o_single
            .iter()
            .copied()
            .filter(|&nb| has_h_neighbor(g, nb))
            .collect();
        let is_se = atom.symbol == "Se";

        if c_neighbors.len() == 1 && o_double.len() == 2 && !o_single_oh.is_empty() {
            // セレノン酸 / テルロン酸: R-Se(=O)2-OH
            let gtype = if is_se { "selenonic_acid" } else { "telluronic_acid" };
            let mut indices = vec![se_idx];
            indices.extend(&c_neighbors);
            indices.extend(&o_double);
            indices.push(o_single_oh[0]);
            push(&mut groups, gtype, indices);
        } else if c_neighbors.len() == 1 && o_double.len() == 1 && !o_single_oh.is_empty() {
            // セレニン酸 / テルリン酸: R-Se(=O)-OH
            let gtype = if is_se { "seleninic_acid" } else { "tellurinic_acid" };
            let mut indices = vec![se_idx];
            indices.extend(&c_neighbors);
            indices.extend(&o_double);
            indices.push(o_single_oh[0]);
            push(&mut groups, gtype, indices);
        } else if c_neighbors.len() == 1
            && !o_single_oh.is_empty()
            && o_double.is_empty()
            && h_neighbors.is_empty()
        {
            // セレネン酸 / テルレン酸: R-Se-OH
            let gtype = if is_se { "selenenic_acid" } else { "tellurenic_acid" };
            let mut indices = vec![se_idx];
            indices.extend(&c_neighbors);
            indices.push(o_single_oh[0]);
            push(&mut groups, gtype, indices);
        } else if !h_neighbors.is_empty() && c_neighbors.len() == 1 {
            // セレノール / テルロール: C-SeH / C-TeH
            let gtype = if is_se { "selenol" } else { "tellurol" };
            push(&mut groups, gtype, vec![c_neighbors[0], se_idx]);
        } else if o_double.len() == 2
            && c_neighbors.len() == 2
            && se_neighbors.is_empty()
            && h_neighbors.is_empty()
        {
            // セレノン / テルロン: C-Se(=O)₂-C
            let gtype = if is_se { "selenone" } else { "telluride" };
            let mut indices = vec![se_idx];
            indices.extend(&o_double);
            indices.extend(&c_neighbors);
            push(&mut groups, gtype, indices);
        } else if o_double.len() == 1
            && c_neighbors.len() == 2
            && se_neighbors.is_empty()
            && h_neighbors.is_empty()
        {
            // セレノキシド / テルロキシド: C-Se(=O)-C
            let gtype = if is_se { "selenoxide" } else { "telluride" };
            let mut indices = vec![se_idx];
            indices.extend(&o_double);
            indices.extend(&c_neighbors);
            push(&mut groups, gtype, indices);
        } else if c_neighbors.len() == 2
            && se_neighbors.is_empty()
            && h_neighbors.is_empty()
            && o_double.is_empty()
        {
            // セレニド / テルリド: C-Se-C / C-Te-C
            let gtype = if is_se { "selenide" } else { "telluride" };
            let mut indices = vec![se_idx];
            indices.extend(&c_neighbors);
            push(&mut groups, gtype, indices);
        } else if c_neighbors.len() == 1 && se_neighbors.len() == 1 && h_neighbors.is_empty() {
            // ジセレニド / ジテルリド: C-Se-Se-C / C-Te-Te-C
            let se2_idx = se_neighbors[0];
            let se2_c: Vec<usize> = g.adjacency[se2_idx]
                .iter()
                .copied()
                .filter(|&nb| nb != se_idx && g.atoms[nb].symbol == "C")
                .collect();
            if !se2_c.is_empty() && se_idx < se2_idx {
                let gtype = if is_se { "diselenide" } else { "ditelluride" };
                let mut indices = vec![se_idx, se2_idx];
                indices.extend(&c_neighbors);
                indices.extend(&se2_c);
                push(&mut groups, gtype, indices);
            }
        }
    }

    // アミン検出: 環外 N を走査 (アミドの N は除外)
    for (n_idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "N" || atom.in_ring {
            continue;
        }
        if atom.formal_charge == 1 {
            // N+ で O/N 隣接なし (純アンモニウム) は ammonium 検出へ
            let has_o = g.adjacency[n_idx].iter().any(|&nb| g.atoms[nb].symbol == "O");
            let has_n = g.adjacency[n_idx].iter().any(|&nb| g.atoms[nb].symbol == "N");
            if !has_o && !has_n {
                continue;
            }
        }
        let h_neighbors = nbrs(g, n_idx, "H");
        let c_neighbors = nbrs(g, n_idx, "C");

        // イソシアネート / イソチオシアネート: R-N=C=O/S
        let c_dbl: Vec<usize> = c_neighbors
            .iter()
            .copied()
            .filter(|&nb| get_bond_order(g, n_idx, nb) == 2.0)
            .collect();
        let c_sgl: Vec<usize> = c_neighbors
            .iter()
            .copied()
            .filter(|&nb| get_bond_order(g, n_idx, nb) == 1.0)
            .collect();
        if c_dbl.len() == 1 && c_sgl.len() == 1 && h_neighbors.is_empty() {
            let iso_c = c_dbl[0];
            for &iso_c_nb in &g.adjacency[iso_c] {
                if get_bond_order(g, iso_c, iso_c_nb) == 2.0 {
                    let sym = &g.atoms[iso_c_nb].symbol;
                    if sym == "O" {
                        push(&mut groups, "isocyanate", vec![c_sgl[0], n_idx, iso_c]);
                        break;
                    } else if sym == "S" {
                        push(
                            &mut groups,
                            "isothiocyanate",
                            vec![c_sgl[0], n_idx, iso_c, iso_c_nb],
                        );
                        break;
                    }
                }
            }
            continue;
        }

        // アジド (zwitterion形): R-[N-]-[N+]≡N
        if atom.formal_charge == -1 && c_neighbors.len() == 1 && h_neighbors.is_empty() {
            let n_sgl_pos: Vec<usize> = g.adjacency[n_idx]
                .iter()
                .copied()
                .filter(|&nb| {
                    g.atoms[nb].symbol == "N"
                        && g.atoms[nb].formal_charge == 1
                        && get_bond_order(g, n_idx, nb) == 1.0
                })
                .collect();
            if n_sgl_pos.len() == 1 {
                let n2_idx = n_sgl_pos[0];
                let n3 = g.adjacency[n2_idx].iter().copied().find(|&nb| {
                    nb != n_idx
                        && g.atoms[nb].symbol == "N"
                        && get_bond_order(g, n2_idx, nb) == 3.0
                });
                if let Some(n3_idx) = n3 {
                    push(
                        &mut groups,
                        "azide",
                        vec![c_neighbors[0], n_idx, n2_idx, n3_idx],
                    );
                    continue;
                }
            }
        }

        // アジド: R-N=[N+]=[N-]
        let n_dbl: Vec<usize> = g.adjacency[n_idx]
            .iter()
            .copied()
            .filter(|&nb| g.atoms[nb].symbol == "N" && get_bond_order(g, n_idx, nb) == 2.0)
            .collect();
        if n_dbl.len() == 1 && c_neighbors.len() == 1 && h_neighbors.is_empty() {
            let n2_idx = n_dbl[0];
            let n3 = g.adjacency[n2_idx]
                .iter()
                .copied()
                .find(|&nb| nb != n_idx && g.atoms[nb].symbol == "N");
            if let Some(n3_idx) = n3 {
                push(
                    &mut groups,
                    "azide",
                    vec![c_neighbors[0], n_idx, n2_idx, n3_idx],
                );
                continue;
            }
        }

        // ニトロソ: R-N=O
        let o_all = nbrs(g, n_idx, "O");
        let o_dbl: Vec<usize> = o_all
            .iter()
            .copied()
            .filter(|&nb| get_bond_order(g, n_idx, nb) == 2.0)
            .collect();
        if o_all.len() == 1 && o_dbl.len() == 1 && c_neighbors.len() == 1 && h_neighbors.is_empty()
        {
            push(&mut groups, "nitroso", vec![c_neighbors[0], n_idx, o_dbl[0]]);
            continue;
        }

        if h_neighbors.len() >= 2 && c_neighbors.len() == 1 {
            // 第一級アミン
            let c_idx = c_neighbors[0];
            if has_double_bonded_oxygen(g, c_idx) || has_double_bonded_sulfur(g, c_idx) {
                continue;
            }
            push(&mut groups, "amine", vec![c_idx, n_idx]);
        } else if h_neighbors.len() == 1 && c_neighbors.len() == 2 {
            // 第二級アミン
            if c_neighbors.iter().any(|&c| has_double_bonded_oxygen(g, c))
                || c_neighbors.iter().any(|&c| has_double_bonded_sulfur(g, c))
            {
                continue;
            }
            let mut indices = vec![n_idx];
            indices.extend(&c_neighbors);
            push(&mut groups, "amine", indices);
        } else if h_neighbors.is_empty() && c_neighbors.len() == 3 {
            // 第三級アミン
            if c_neighbors.iter().any(|&c| has_double_bonded_oxygen(g, c))
                || c_neighbors.iter().any(|&c| has_double_bonded_sulfur(g, c))
            {
                continue;
            }
            let mut indices = vec![n_idx];
            indices.extend(&c_neighbors);
            push(&mut groups, "amine", indices);
        } else {
            // N-ハロアミン (N に C + ハロゲンが付く場合)
            let has_n_halogen = g.adjacency[n_idx]
                .iter()
                .any(|&nb| is_halogen(&g.atoms[nb].symbol));
            if has_n_halogen
                && !c_neighbors.is_empty()
                && !c_neighbors.iter().any(|&c| has_double_bonded_oxygen(g, c))
            {
                let mut indices = vec![n_idx];
                indices.extend(&c_neighbors);
                push(&mut groups, "amine", indices);
            }
        }
    }

    // アルケン: C=C (非芳香族、非環状のみ)
    let mut seen_double: HashSet<(usize, usize)> = HashSet::new();
    for (a_idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "C" || atom.in_ring {
            continue;
        }
        for &neighbor_idx in &g.adjacency[a_idx] {
            let neighbor = &g.atoms[neighbor_idx];
            if neighbor.symbol != "C" || neighbor.in_ring {
                continue;
            }
            if get_bond_order(g, a_idx, neighbor_idx) == 2.0 {
                let pair = (a_idx.min(neighbor_idx), a_idx.max(neighbor_idx));
                if seen_double.insert(pair) {
                    push(&mut groups, "alkene", vec![pair.0, pair.1]);
                }
            }
        }
    }

    // アルキン: C≡C (非環状のみ)
    let mut seen_triple: HashSet<(usize, usize)> = HashSet::new();
    for (a_idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "C" || atom.in_ring {
            continue;
        }
        for &neighbor_idx in &g.adjacency[a_idx] {
            let neighbor = &g.atoms[neighbor_idx];
            if neighbor.symbol != "C" || neighbor.in_ring {
                continue;
            }
            if get_bond_order(g, a_idx, neighbor_idx) == 3.0 {
                let pair = (a_idx.min(neighbor_idx), a_idx.max(neighbor_idx));
                if seen_triple.insert(pair) {
                    push(&mut groups, "alkyne", vec![pair.0, pair.1]);
                }
            }
        }
    }

    // リン化合物検出
    for (p_idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "P" || atom.in_ring {
            continue;
        }
        let c_neighbors = nbrs(g, p_idx, "C");
        let o_neighbors = nbrs(g, p_idx, "O");
        let o_double: Vec<usize> = o_neighbors
            .iter()
            .copied()
            .filter(|&nb| get_bond_order(g, p_idx, nb) == 2.0)
            .collect();
        let o_single: Vec<usize> = o_neighbors
            .iter()
            .copied()
            .filter(|&nb| get_bond_order(g, p_idx, nb) == 1.0)
            .collect();
        let o_single_oh: Vec<usize> = o_single
            .iter()
            .copied()
            .filter(|&nb| has_h_neighbor(g, nb))
            .collect();
        let o_ester_p: Vec<usize> = o_single
            .iter()
            .copied()
            .filter(|&nb| {
                !o_single_oh.contains(&nb)
                    && g.adjacency[nb]
                        .iter()
                        .any(|&occ| occ != p_idx && g.atoms[occ].symbol == "C")
            })
            .collect();

        if o_double.len() == 1 && o_single_oh.len() >= 2 && c_neighbors.len() == 1 {
            // ホスホン酸: R-P(=O)(OH)2
            let mut indices = vec![p_idx];
            indices.extend(&c_neighbors);
            indices.extend(&o_double);
            indices.extend(&o_single_oh);
            push(&mut groups, "phosphonic_acid", indices);
        } else if o_double.len() == 1
            && c_neighbors.len() == 1
            && !o_ester_p.is_empty()
            && !o_single_oh.is_empty()
        {
            // ホスホン酸部分エステル: R-P(=O)(OR')(OH)
            let mut indices = vec![p_idx];
            indices.extend(&c_neighbors);
            indices.extend(&o_double);
            indices.extend(&o_ester_p);
            indices.extend(&o_single_oh);
            push(&mut groups, "phosphonate_halfester", indices);
        } else if o_double.len() == 1 && !o_single_oh.is_empty() && c_neighbors.len() == 2 {
            // ホスフィン酸: R2P(=O)(OH)
            let mut indices = vec![p_idx];
            indices.extend(&c_neighbors);
            indices.extend(&o_double);
            indices.extend(&o_single_oh);
            push(&mut groups, "phosphinic_acid", indices);
        } else if o_double.len() == 1
            && o_single_oh.len() == 1
            && c_neighbors.len() == 1
            && o_ester_p.is_empty()
        {
            // ホスフィン酸 (モノアルキル): R-PH(=O)(OH)
            let mut indices = vec![p_idx];
            indices.extend(&c_neighbors);
            indices.extend(&o_double);
            indices.extend(&o_single_oh);
            push(&mut groups, "phosphinic_acid", indices);
        } else if o_double.is_empty() && o_single_oh.len() >= 2 && c_neighbors.len() == 1 {
            // ホスホナス酸: R-P(OH)2
            let mut indices = vec![p_idx];
            indices.extend(&c_neighbors);
            indices.extend(&o_single_oh[..2]);
            push(&mut groups, "phosphonous_acid", indices);
        } else if o_double.is_empty() && o_single_oh.len() == 1 && !c_neighbors.is_empty() {
            // ホスフィナス酸: R_n-PH_{2-n}(OH)
            let mut indices = vec![p_idx];
            indices.extend(&c_neighbors);
            indices.extend(&o_single_oh);
            push(&mut groups, "phosphinous_acid", indices);
        } else if !c_neighbors.is_empty() && o_neighbors.is_empty() {
            // ホスファン: R_n-PH_{3-n}
            let mut indices = vec![p_idx];
            indices.extend(&c_neighbors);
            push(&mut groups, "phosphane", indices);
        } else if o_double.len() == 1 && c_neighbors.is_empty() {
            // ホスフェートエステル: (RO)_n P(=O)(OH)_{3-n}
            let o_ester: Vec<usize> = o_single
                .iter()
                .copied()
                .filter(|&nb| {
                    !o_single_oh.contains(&nb)
                        && g.adjacency[nb]
                            .iter()
                            .any(|&occ| occ != p_idx && g.atoms[occ].symbol == "C")
                })
                .collect();
            if !o_ester.is_empty() {
                let mut indices = vec![p_idx];
                indices.extend(&o_double);
                indices.extend(&o_ester);
                indices.extend(&o_single_oh);
                push(&mut groups, "phosphate_ester", indices);
            }
        } else if o_double.len() == 1 && c_neighbors.len() == 1 && o_single_oh.is_empty() {
            // ホスホネートエステル: R-P(=O)(OR)2
            let o_ester: Vec<usize> = o_single
                .iter()
                .copied()
                .filter(|&nb| {
                    g.adjacency[nb]
                        .iter()
                        .any(|&occ| occ != p_idx && g.atoms[occ].symbol == "C")
                })
                .collect();
            if !o_ester.is_empty() {
                let mut indices = vec![p_idx];
                indices.extend(&c_neighbors);
                indices.extend(&o_double);
                indices.extend(&o_ester);
                push(&mut groups, "phosphonate_ester", indices);
            }
        } else if o_double.is_empty() && c_neighbors.len() == 1 && o_single_oh.is_empty() {
            // ホスホノチオアートエステル: R-P(=S)(OR)2
            let s_double_p: Vec<usize> = g.adjacency[p_idx]
                .iter()
                .copied()
                .filter(|&nb| g.atoms[nb].symbol == "S" && get_bond_order(g, p_idx, nb) == 2.0)
                .collect();
            if !s_double_p.is_empty() {
                let o_ester: Vec<usize> = o_single
                    .iter()
                    .copied()
                    .filter(|&nb| {
                        g.adjacency[nb]
                            .iter()
                            .any(|&occ| occ != p_idx && g.atoms[occ].symbol == "C")
                    })
                    .collect();
                if !o_ester.is_empty() {
                    let mut indices = vec![p_idx];
                    indices.extend(&c_neighbors);
                    indices.extend(&s_double_p);
                    indices.extend(&o_ester);
                    // Python 同様 phosphonate_ester の優先度を使う
                    groups.push(FunctionalGroup {
                        group_type: "phosphonothioate_ester",
                        atom_indices: indices,
                        priority: prio("phosphonate_ester"),
                    });
                }
            }
        } else if o_double.len() == 1 && c_neighbors.len() == 2 && o_single_oh.is_empty() {
            // ホスフィネートエステル: R2-P(=O)(OR)
            let o_ester: Vec<usize> = o_single
                .iter()
                .copied()
                .filter(|&nb| {
                    g.adjacency[nb]
                        .iter()
                        .any(|&occ| occ != p_idx && g.atoms[occ].symbol == "C")
                })
                .collect();
            if !o_ester.is_empty() {
                let mut indices = vec![p_idx];
                indices.extend(&c_neighbors);
                indices.extend(&o_double);
                indices.extend(&o_ester);
                push(&mut groups, "phosphinate_ester", indices);
            }
        } else if o_double.len() == 1 && c_neighbors.len() >= 3 && o_single.is_empty() {
            // ホスフィンオキシド: R3P=O
            let mut indices = vec![p_idx];
            indices.extend(&c_neighbors[..3]);
            indices.extend(&o_double);
            push(&mut groups, "phosphine_oxide", indices);
        } else if o_double.is_empty()
            && c_neighbors.is_empty()
            && o_single.len() >= 3
            && o_single_oh.is_empty()
        {
            // 亜リン酸トリエステル: (RO)3P
            let o_ester: Vec<usize> = o_single
                .iter()
                .copied()
                .filter(|&nb| {
                    g.adjacency[nb]
                        .iter()
                        .any(|&occ| occ != p_idx && g.atoms[occ].symbol == "C")
                })
                .collect();
            if o_ester.len() >= 3 {
                let mut indices = vec![p_idx];
                indices.extend(&o_ester[..3]);
                push(&mut groups, "phosphite_ester", indices);
            }
        }
    }

    // ヒ素化合物検出
    for (as_idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "As" || atom.in_ring {
            continue;
        }
        let c_neighbors = nbrs(g, as_idx, "C");
        let o_neighbors = nbrs(g, as_idx, "O");
        let o_double: Vec<usize> = o_neighbors
            .iter()
            .copied()
            .filter(|&nb| get_bond_order(g, as_idx, nb) == 2.0)
            .collect();
        let o_single_oh: Vec<usize> = o_neighbors
            .iter()
            .copied()
            .filter(|&nb| get_bond_order(g, as_idx, nb) == 1.0 && has_h_neighbor(g, nb))
            .collect();

        if o_double.len() == 1 && o_single_oh.len() >= 2 && c_neighbors.len() == 1 {
            // ヒ素酸 (arsonic): R-As(=O)(OH)2
            let mut indices = vec![as_idx];
            indices.extend(&c_neighbors);
            indices.extend(&o_double);
            indices.extend(&o_single_oh[..2]);
            push(&mut groups, "arsonic_acid", indices);
        } else if o_double.len() == 1 && !o_single_oh.is_empty() && c_neighbors.len() >= 2 {
            // アルシン酸 (arsinic): R2As(=O)(OH)
            let mut indices = vec![as_idx];
            indices.extend(&c_neighbors);
            indices.extend(&o_double);
            indices.push(o_single_oh[0]);
            push(&mut groups, "arsinic_acid", indices);
        } else if o_double.is_empty() && o_single_oh.len() >= 2 && c_neighbors.len() == 1 {
            // 亜ヒ酸 (arsonous): R-As(OH)2
            let mut indices = vec![as_idx];
            indices.extend(&c_neighbors);
            indices.extend(&o_single_oh[..2]);
            push(&mut groups, "arsonous_acid", indices);
        } else if o_double.is_empty() && o_single_oh.len() == 1 && !c_neighbors.is_empty() {
            // 亜アルシン酸 (arsinous): R_n-As-OH
            let mut indices = vec![as_idx];
            indices.extend(&c_neighbors);
            indices.extend(&o_single_oh);
            push(&mut groups, "arsinous_acid", indices);
        } else if !c_neighbors.is_empty() && o_neighbors.is_empty() {
            // ヒ化水素 (arsane): R_n-AsH_{3-n}
            let mut indices = vec![as_idx];
            indices.extend(&c_neighbors);
            push(&mut groups, "arsane_org", indices);
        }
    }

    // 有機水銀化合物検出
    for (hg_idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "Hg" || atom.in_ring {
            continue;
        }
        let c_neighbors = nbrs(g, hg_idx, "C");
        if !c_neighbors.is_empty() {
            let mut indices = vec![hg_idx];
            indices.extend(&c_neighbors);
            push(&mut groups, "organomercury", indices);
        }
    }

    // ホウ素化合物検出
    for (b_idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "B" || atom.in_ring {
            continue;
        }
        let c_neighbors = nbrs(g, b_idx, "C");
        let o_neighbors = nbrs(g, b_idx, "O");
        let o_single_oh: Vec<usize> = o_neighbors
            .iter()
            .copied()
            .filter(|&nb| get_bond_order(g, b_idx, nb) == 1.0 && has_h_neighbor(g, nb))
            .collect();

        if o_single_oh.len() >= 2 && c_neighbors.len() == 1 {
            // ボロン酸: R-B(OH)2
            let mut indices = vec![b_idx];
            indices.extend(&c_neighbors);
            indices.extend(&o_single_oh);
            push(&mut groups, "boronic_acid", indices);
        } else if o_single_oh.len() == 1 && c_neighbors.len() == 2 {
            // ボリン酸: R2B(OH)
            let mut indices = vec![b_idx];
            indices.extend(&c_neighbors);
            indices.extend(&o_single_oh);
            push(&mut groups, "borinic_acid", indices);
        } else if c_neighbors.len() == 1 && !o_neighbors.is_empty() && o_single_oh.is_empty() {
            // ボロン酸エステル: R-B(OR')2
            let o_ester: Vec<usize> = o_neighbors
                .iter()
                .copied()
                .filter(|&nb| {
                    g.adjacency[nb]
                        .iter()
                        .any(|&occ| occ != b_idx && g.atoms[occ].symbol == "C")
                })
                .collect();
            if o_ester.len() >= 2 {
                let mut indices = vec![b_idx];
                indices.extend(&c_neighbors);
                indices.extend(&o_ester[..2]);
                push(&mut groups, "boronate_ester", indices);
            }
        } else if !c_neighbors.is_empty() && o_neighbors.is_empty() {
            // ボラン: R_n-BH_{3-n}
            let mut indices = vec![b_idx];
            indices.extend(&c_neighbors);
            push(&mut groups, "borane_org", indices);
        } else if c_neighbors.is_empty() && o_neighbors.len() >= 3 {
            // トリアルコキシボラン: B(OR)3
            let o_ester: Vec<usize> = o_neighbors
                .iter()
                .copied()
                .filter(|&nb| {
                    g.adjacency[nb]
                        .iter()
                        .any(|&occ| occ != b_idx && g.atoms[occ].symbol == "C")
                })
                .collect();
            if o_ester.len() >= 3 {
                let mut indices = vec![b_idx];
                indices.extend(&o_ester[..3]);
                push(&mut groups, "borate_ester", indices);
            }
        }
    }

    // ケイ素化合物検出
    let mut seen_si: HashSet<usize> = HashSet::new();
    for (si_idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "Si" || atom.in_ring {
            continue;
        }
        if seen_si.contains(&si_idx) {
            continue;
        }
        let c_neighbors = nbrs(g, si_idx, "C");
        let o_neighbors = nbrs(g, si_idx, "O");
        let n_neighbors = nbrs(g, si_idx, "N");
        let si_neighbors = nbrs(g, si_idx, "Si");
        // Si-O-Si ジシロキサン検出
        let mut si_o_si: Option<(usize, usize)> = None;
        for &o_idx in &o_neighbors {
            if has_h_neighbor(g, o_idx) {
                continue;
            }
            let si2 = g.adjacency[o_idx]
                .iter()
                .copied()
                .find(|&nb| nb != si_idx && g.atoms[nb].symbol == "Si");
            if let Some(si2) = si2 {
                si_o_si = Some((o_idx, si2));
                break;
            }
        }
        // Si-N-Si ジシラザン検出
        let mut si_n_si: Option<(usize, usize)> = None;
        if si_o_si.is_none() {
            for &n_idx in &n_neighbors {
                let si2 = g.adjacency[n_idx]
                    .iter()
                    .copied()
                    .find(|&nb| nb != si_idx && g.atoms[nb].symbol == "Si");
                if let Some(si2) = si2 {
                    si_n_si = Some((n_idx, si2));
                    break;
                }
            }
        }
        // Si-Si ジシラン検出
        let mut si_si: Option<usize> = None;
        if si_o_si.is_none() && si_n_si.is_none() && !si_neighbors.is_empty() {
            si_si = si_neighbors
                .iter()
                .copied()
                .find(|nb| !seen_si.contains(nb));
        }
        if let Some((o_idx, si2_idx)) = si_o_si {
            seen_si.insert(si_idx);
            seen_si.insert(si2_idx);
            let c2 = nbrs(g, si2_idx, "C");
            let mut indices = vec![si_idx, si2_idx, o_idx];
            indices.extend(&c_neighbors);
            indices.extend(&c2);
            push(&mut groups, "disiloxane_org", indices);
        } else if let Some((n_idx, si2_idx)) = si_n_si {
            seen_si.insert(si_idx);
            seen_si.insert(si2_idx);
            let c2 = nbrs(g, si2_idx, "C");
            let mut indices = vec![si_idx, si2_idx, n_idx];
            indices.extend(&c_neighbors);
            indices.extend(&c2);
            push(&mut groups, "disilazane_org", indices);
        } else if let Some(si2_idx) = si_si {
            seen_si.insert(si_idx);
            seen_si.insert(si2_idx);
            let c2 = nbrs(g, si2_idx, "C");
            let mut indices = vec![si_idx, si2_idx];
            indices.extend(&c_neighbors);
            indices.extend(&c2);
            push(&mut groups, "disilane_org", indices);
        } else {
            let o_oh: Vec<usize> = o_neighbors
                .iter()
                .copied()
                .filter(|&nb| has_h_neighbor(g, nb))
                .collect();
            let o_or: Vec<usize> = o_neighbors
                .iter()
                .copied()
                .filter(|&nb| {
                    !o_oh.contains(&nb)
                        && g.adjacency[nb].iter().any(|&occ| g.atoms[occ].symbol == "C")
                })
                .collect();
            if !c_neighbors.is_empty() && !o_oh.is_empty() {
                // シラノール/ジオール/トリオール: R_n Si(OH)_{4-n}
                let mut indices = vec![si_idx];
                indices.extend(&c_neighbors);
                indices.extend(&o_oh);
                push(&mut groups, "silanol_org", indices);
            } else if !o_or.is_empty() {
                // シリルエーテル / テトラアルコキシシラン
                let mut indices = vec![si_idx];
                indices.extend(&c_neighbors);
                indices.extend(&o_or);
                push(&mut groups, "silyl_ether_org", indices);
            } else if !c_neighbors.is_empty() {
                let mut indices = vec![si_idx];
                indices.extend(&c_neighbors);
                push(&mut groups, "silane_org", indices);
            }
        }
    }

    // ゲルマン・スタンナン検出
    for (central_idx, atom) in g.atoms.iter().enumerate() {
        if !matches!(atom.symbol.as_str(), "Ge" | "Sn") || atom.in_ring {
            continue;
        }
        let c_neighbors = nbrs(g, central_idx, "C");
        if !c_neighbors.is_empty() {
            let gtype = if atom.symbol == "Ge" { "germane_org" } else { "stannane_org" };
            let mut indices = vec![central_idx];
            indices.extend(&c_neighbors);
            push(&mut groups, gtype, indices);
        }
    }

    // ビスマス・アンチモン・鉛 有機水素化物
    for (central_idx, atom) in g.atoms.iter().enumerate() {
        if !matches!(atom.symbol.as_str(), "Bi" | "Sb" | "Pb") || atom.in_ring {
            continue;
        }
        let c_neighbors = nbrs(g, central_idx, "C");
        if !c_neighbors.is_empty() {
            let gtype = match atom.symbol.as_str() {
                "Bi" => "bismuthane_org",
                "Sb" => "stibane_org",
                _ => "plumbane_org",
            };
            let mut indices = vec![central_idx];
            indices.extend(&c_neighbors);
            push(&mut groups, gtype, indices);
        }
    }

    // イソシアニド検出 ([C-]#[N+]-R)
    for (n_idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "N" || atom.in_ring || atom.formal_charge != 1 {
            continue;
        }
        let cn_neighbors: Vec<usize> = g.adjacency[n_idx]
            .iter()
            .copied()
            .filter(|&nb| {
                g.atoms[nb].symbol == "C"
                    && g.atoms[nb].formal_charge == -1
                    && get_bond_order(g, n_idx, nb) == 3.0
            })
            .collect();
        if !cn_neighbors.is_empty() {
            let c_alkyl: Vec<usize> = g.adjacency[n_idx]
                .iter()
                .copied()
                .filter(|&nb| !cn_neighbors.contains(&nb) && g.atoms[nb].symbol == "C")
                .collect();
            let h_neighbors = nbrs(g, n_idx, "H");
            if !c_alkyl.is_empty() || !h_neighbors.is_empty() {
                let mut indices = vec![n_idx];
                indices.extend(&cn_neighbors);
                indices.extend(if c_alkyl.is_empty() { &h_neighbors } else { &c_alkyl });
                push(&mut groups, "isocyanide", indices);
            }
        }
    }

    // アンモニウムイオン検出
    for (n_idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "N" || atom.in_ring || atom.formal_charge != 1 {
            continue;
        }
        // O 隣接あり → N-oxide or nitro / N 隣接あり → diazo/azo
        if g.adjacency[n_idx].iter().any(|&nb| g.atoms[nb].symbol == "O")
            || g.adjacency[n_idx].iter().any(|&nb| g.atoms[nb].symbol == "N")
        {
            continue;
        }
        let c_neighbors = nbrs(g, n_idx, "C");
        // イソシアニド R-N≡C は ammonium ではない
        if c_neighbors.iter().any(|&cn| {
            g.atoms[cn].formal_charge == -1 && get_bond_order(g, n_idx, cn) == 3.0
        }) {
            continue;
        }
        if !c_neighbors.is_empty() {
            let mut indices = vec![n_idx];
            indices.extend(&c_neighbors);
            push(&mut groups, "ammonium", indices);
        }
    }

    // ホスホニウム / スルホニウム / アルソニウム検出
    for (het_idx, atom) in g.atoms.iter().enumerate() {
        if !matches!(atom.symbol.as_str(), "P" | "S" | "As") || atom.formal_charge != 1 {
            continue;
        }
        let c_nbs = nbrs(g, het_idx, "C");
        let gtype = match (atom.symbol.as_str(), c_nbs.len()) {
            ("P", 4) => "phosphanium",
            ("As", 4) => "arsonium",
            ("S", 3) => "sulfonium",
            _ => continue,
        };
        let mut indices = vec![het_idx];
        indices.extend(&c_nbs);
        push(&mut groups, gtype, indices);
    }

    // ニトレート/ニトライトエステル検出
    for (n_idx, atom) in g.atoms.iter().enumerate() {
        if atom.symbol != "N" || atom.in_ring {
            continue;
        }
        let o_double: Vec<usize> = g.adjacency[n_idx]
            .iter()
            .copied()
            .filter(|&nb| g.atoms[nb].symbol == "O" && get_bond_order(g, n_idx, nb) == 2.0)
            .collect();
        let o_single: Vec<usize> = g.adjacency[n_idx]
            .iter()
            .copied()
            .filter(|&nb| g.atoms[nb].symbol == "O" && get_bond_order(g, n_idx, nb) == 1.0)
            .collect();
        if atom.formal_charge == 1 && o_double.len() == 1 && o_single.len() == 2 {
            // Nitrate ester: R-O-[N+](=O)[O-]
            let o_alkyl: Vec<usize> = o_single
                .iter()
                .copied()
                .filter(|&ob| g.adjacency[ob].iter().any(|&c| g.atoms[c].symbol == "C"))
                .collect();
            let o_neg: Vec<usize> = o_single
                .iter()
                .copied()
                .filter(|&ob| g.atoms[ob].formal_charge == -1)
                .collect();
            if o_alkyl.len() == 1 && !o_neg.is_empty() {
                let mut indices = vec![n_idx, o_alkyl[0]];
                indices.extend(&o_double);
                indices.push(o_neg[0]);
                push(&mut groups, "nitrate_ester", indices);
            }
        } else if atom.formal_charge == 0 && o_double.len() == 1 && o_single.len() == 1 {
            // Nitrite ester: R-O-N=O
            let o_alkyl_n: Vec<usize> = o_single
                .iter()
                .copied()
                .filter(|&ob| g.adjacency[ob].iter().any(|&c| g.atoms[c].symbol == "C"))
                .collect();
            if o_alkyl_n.len() == 1 {
                let mut indices = vec![n_idx, o_alkyl_n[0]];
                indices.extend(&o_double);
                push(&mut groups, "nitrite_ester", indices);
            }
        }
    }

    // アルカン: 上記なし
    if groups.is_empty() {
        push(&mut groups, "alkane", vec![]);
    }

    // 優先順位の高い順に安定ソート
    groups.sort_by_key(|fg| std::cmp::Reverse(fg.priority));

    // 同一最高優先度グループを diol/dione/dioic_acid 等に集約
    aggregate_groups(groups, g)
}

// ─── 集約と主基選択 ─────────────────────────────────────────────

/// 最高優先度の同一 group_type が複数ある場合、diol/dione/dioic_acid 等に集約する。
pub fn aggregate_groups(
    groups: Vec<FunctionalGroup>,
    g: &MoleculeGraph,
) -> Vec<FunctionalGroup> {
    if groups.is_empty() {
        return groups;
    }

    let top_type = groups[0].group_type;
    let top_priority = groups[0].priority;
    let top_groups: Vec<&FunctionalGroup> = groups
        .iter()
        .filter(|fg| fg.group_type == top_type)
        .collect();

    if top_groups.len() <= 1 {
        return groups;
    }

    // carboxylic_acid のみ: 芳香族環直結と側鎖が混在する場合はマージしない
    if top_type == "carboxylic_acid" {
        let is_ring_cooh = |fg: &FunctionalGroup| -> bool {
            let anchor_c = fg.atom_indices[0];
            g.adjacency[anchor_c].iter().any(|&nb| {
                g.atoms[nb].symbol == "C" && g.atoms[nb].in_ring && g.atoms[nb].is_aromatic
            })
        };
        let n_ring = top_groups.iter().filter(|fg| is_ring_cooh(fg)).count();
        let n_chain = top_groups.len() - n_ring;
        if n_ring > 0 && n_chain > 0 && n_ring == 1 {
            return groups; // 集約しない (ring COOH が 1 つのみ)
        }
    }

    // ketone 特例: 両方の C=O が同一芳香族環に直結している場合はマージしない
    if top_type == "ketone" && top_groups.len() == 2 {
        let adj_arom_ring_idxs = |c: usize| -> Vec<usize> {
            g.ring_atom_sets
                .iter()
                .enumerate()
                .filter(|(_, rs)| {
                    g.adjacency[c].iter().any(|&nb| {
                        g.atoms[nb].symbol == "C"
                            && g.atoms[nb].is_aromatic
                            && rs.contains(&nb)
                    })
                })
                .map(|(i, _)| i)
                .collect()
        };
        let rs0 = adj_arom_ring_idxs(top_groups[0].atom_indices[0]);
        let rs1 = adj_arom_ring_idxs(top_groups[1].atom_indices[0]);
        if !rs0.is_empty() && !rs1.is_empty() && rs0.iter().any(|s| rs1.contains(s)) {
            return groups; // 同一芳香族環直結: マージしない
        }
    }

    let multi_type = match (top_type, top_groups.len()) {
        ("carboxylic_acid", 2) => "dioic_acid",
        ("alcohol", 2) => "diol",
        ("alcohol", 3) => "triol",
        ("alcohol", 4) => "tetraol",
        ("alcohol", 5) => "pentaol",
        ("alcohol", 6) => "hexaol",
        ("ketone", 2) => "dione",
        ("ketone", 3) => "trione",
        ("ketone", 4) => "tetraone",
        ("ketone", 5) => "pentaone",
        ("aldehyde", 2) => "dial",
        ("ester", 2) => "diester",
        ("amine", 2) => "diamine",
        ("amine", 3) => "triamine",
        ("acid_halide", 2) => "diacid_halide",
        ("nitrile", 2) => "dinitrile",
        ("carboxylate", 2) => "dicarboxylate",
        ("thiol", 2) => "dithiol",
        ("amide", 2) => "diamide",
        ("sulfonic_acid", 2) => "disulfonic_acid",
        ("sulfonamide", 2) => "disulfonamide",
        ("isocyanate", 2) => "diisocyanate",
        ("isothiocyanate", 2) => "diisothiocyanate",
        ("imine", 2) => "diimine",
        ("thioamide", 2) => "dithioamide",
        ("selenoamide", 2) => "diselenoamide",
        ("amidine", 2) => "diamidine",
        _ => return groups, // 3+ ketones 等の未対応ケースはそのまま返す
    };

    // アミドのマージ特例: 同一 N 原子を共有する 2 つのアミド (N-アシルアミド) は
    // diamide に集約しない
    if top_type == "amide" && top_groups.len() == 2 {
        let n_of = |fg: &FunctionalGroup| -> Option<usize> {
            fg.atom_indices
                .iter()
                .copied()
                .find(|&ai| g.atoms[ai].symbol == "N")
        };
        let n0 = n_of(top_groups[0]);
        if n0.is_some() && n0 == n_of(top_groups[1]) {
            return groups;
        }
    }

    let mut seen: HashSet<usize> = HashSet::new();
    let mut merged_atoms: Vec<usize> = Vec::new();
    for fg in &top_groups {
        for &ai in &fg.atom_indices {
            if seen.insert(ai) {
                merged_atoms.push(ai);
            }
        }
    }

    let merged = FunctionalGroup {
        group_type: multi_type,
        atom_indices: merged_atoms,
        priority: top_priority,
    };
    let mut result = vec![merged];
    result.extend(groups.iter().filter(|fg| fg.group_type != top_type).cloned());
    result
}

/// 官能基リストの最優先グループを返す (alkane の場合は None)。
pub fn principal_group(groups: &[FunctionalGroup]) -> Option<&FunctionalGroup> {
    let top = groups.first()?;
    if top.group_type == "alkane" {
        None
    } else {
        Some(top)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use molrs::graph::build_molecule_graph;

    fn types_of(smiles: &str) -> Vec<&'static str> {
        let g = build_molecule_graph(smiles).unwrap();
        detect_groups(&g).iter().map(|fg| fg.group_type).collect()
    }

    #[test]
    fn basic_groups() {
        assert_eq!(types_of("CC(=O)O")[0], "carboxylic_acid");
        assert_eq!(types_of("CCO")[0], "alcohol");
        assert_eq!(types_of("CC(=O)C")[0], "ketone");
        assert_eq!(types_of("CC=O")[0], "aldehyde");
        assert_eq!(types_of("CC#N")[0], "nitrile");
        assert_eq!(types_of("CCN")[0], "amine");
        assert_eq!(types_of("CC(=O)OC")[0], "ester");
        assert_eq!(types_of("CC(=O)N")[0], "amide");
        assert_eq!(types_of("CCS")[0], "thiol");
        assert_eq!(types_of("C=C")[0], "alkene");
        assert_eq!(types_of("C#C")[0], "alkyne");
        assert_eq!(types_of("CC")[0], "alkane");
    }

    #[test]
    fn aggregation() {
        assert_eq!(types_of("OCCO")[0], "diol");
        assert_eq!(types_of("OC(=O)CCC(=O)O")[0], "dioic_acid");
        assert_eq!(types_of("CC(=O)CC(=O)C")[0], "dione");
    }

    #[test]
    fn principal_group_alkane_is_none() {
        let g = build_molecule_graph("CC").unwrap();
        let groups = detect_groups(&g);
        assert!(principal_group(&groups).is_none());
    }

    #[test]
    fn sulfur_groups() {
        assert_eq!(types_of("CS(=O)(=O)O")[0], "sulfonic_acid");
        assert_eq!(types_of("CS(=O)(=O)N")[0], "sulfonamide");
        assert_eq!(types_of("CS(=O)C")[0], "sulfoxide");
        assert_eq!(types_of("CS(=O)(=O)C")[0], "sulfone");
        assert_eq!(types_of("CSC")[0], "sulfide");
        assert_eq!(types_of("CSSC")[0], "disulfide");
    }
}
