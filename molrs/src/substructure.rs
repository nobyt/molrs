//! 部分構造マッチ (S1.6)。RDKit `GetSubstructMatches(uniquify=False)` 相当。
//!
//! heterocycle_handler.py の 2 つの使用形態を再現する:
//! - `MolQuery`: `MolFromSmiles` で作った分子をクエリにする通常形。
//!   原子は原子番号のみ照合 (電荷はクエリ側が非ゼロのときだけ一致要求、
//!   H 数・芳香族フラグは照合しない、`*` はダミー同士のみ)。
//!   結合は次数クラスの完全一致 (クエリも同じサニタイズを通るため
//!   ケクレ表記クエリは芳香族化されて 1.5 同士で一致する)。
//! - SMARTS 形 (`substruct_matches_smarts`): SMILES 文字列を SMARTS として
//!   解釈する形 (`MolFromSmarts(core_smi)`)。SMARTS はサニタイズを通らない
//!   (ケクレ化不能なパターンも有効) ため、パース結果から直接クエリを作る。
//!   小文字/大文字が芳香族/脂肪族の制約、角括弧の H 指定 (`[nH]`) が
//!   総 H 数の制約になる。省略結合は「単結合または芳香族」。`*` は任意原子。
//!
//! 列挙は VF2 系バックトラッキング。返り値の各タプルは
//! クエリ原子インデックス順のターゲット原子インデックス。
//! 列挙順は RDKit と一致させない (呼び出し側は集合として扱う)。

use crate::elements::atomic_number;
use crate::graph::MoleculeGraph;
use crate::smiles::{parse_smiles, BondKind};
use crate::ChemError;

struct QAtom {
    atomic_num: u8,
    /// Some(flag) なら芳香族性を照合 (SMARTS 形のみ)
    aromatic: Option<bool>,
    /// Some(c) なら形式電荷を照合
    charge: Option<i8>,
    /// Some(n) なら総 H 数を照合 (SmartsLike の角括弧原子のみ)
    n_h: Option<u8>,
    /// 任意原子 (SmartsLike の `*`)
    wildcard: bool,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum QBondKind {
    /// 次数クラス完全一致 (0=単, 1=芳香族, 2=二重, 3=三重, 4=四重)
    Exact(u8),
    /// 単結合または芳香族 (SMARTS の省略結合)
    SingleOrAromatic,
}

/// クエリ表現: (原子制約, 隣接リスト)
type QueryRep = (Vec<QAtom>, Vec<Vec<(usize, QBondKind)>>);

/// ターゲット側の前処理済みビュー。
struct TargetView {
    n: usize,
    atomic_num: Vec<u8>,
    aromatic: Vec<bool>,
    charge: Vec<i8>,
    n_h: Vec<u8>,
    /// atom → [(相手, 次数クラス)]
    adj: Vec<Vec<(usize, u8)>>,
}

fn order_class(order: f64) -> u8 {
    if order == 1.5 {
        1
    } else if order == 2.0 {
        2
    } else if order == 3.0 {
        3
    } else if order == 4.0 {
        4
    } else {
        0
    }
}

fn num_kept_atoms(g: &MoleculeGraph) -> usize {
    g.parser_to_graph.iter().flatten().count()
}

fn build_target_view(g: &MoleculeGraph) -> TargetView {
    let n = num_kept_atoms(g);
    let mut adj = vec![Vec::new(); n];
    for b in &g.bonds {
        if b.begin_idx < n && b.end_idx < n {
            let oc = order_class(b.bond_order);
            adj[b.begin_idx].push((b.end_idx, oc));
            adj[b.end_idx].push((b.begin_idx, oc));
        }
    }
    let n_h: Vec<u8> = (0..n)
        .map(|i| g.adjacency[i].iter().filter(|&&x| x >= n).count() as u8)
        .collect();
    TargetView {
        n,
        atomic_num: (0..n).map(|i| g.atoms[i].atomic_num).collect(),
        aromatic: (0..n).map(|i| g.atoms[i].is_aromatic).collect(),
        charge: (0..n).map(|i| g.atoms[i].formal_charge).collect(),
        n_h,
        adj,
    }
}

/// Mol クエリ (サニタイズ済み MoleculeGraph) を制約表現にする。
fn build_mol_query(g: &MoleculeGraph) -> QueryRep {
    let n = num_kept_atoms(g);
    let mut atoms = Vec::with_capacity(n);
    for gi in 0..n {
        let a = &g.atoms[gi];
        atoms.push(QAtom {
            atomic_num: a.atomic_num,
            aromatic: None,
            charge: (a.formal_charge != 0).then_some(a.formal_charge),
            n_h: None,
            wildcard: false, // mol クエリの `*` はダミー同士のみ (atomic_num 0 で照合)
        });
    }
    let mut adj = vec![Vec::new(); n];
    for b in &g.bonds {
        if b.begin_idx >= n || b.end_idx >= n {
            continue;
        }
        let k = QBondKind::Exact(order_class(b.bond_order));
        adj[b.begin_idx].push((b.end_idx, k));
        adj[b.end_idx].push((b.begin_idx, k));
    }
    (atoms, adj)
}

/// SMILES 文字列を SMARTS として解釈したクエリ。サニタイズなし
/// (RDKit MolFromSmarts と同じく、ケクレ化不能なパターンも許す)。
fn build_smarts_query(pattern: &str) -> Result<QueryRep, ChemError> {
    let parsed = parse_smiles(pattern)?;
    let n = parsed.atoms.len();
    let mut atoms = Vec::with_capacity(n);
    for pa in &parsed.atoms {
        atoms.push(QAtom {
            atomic_num: atomic_number(&pa.symbol).unwrap_or(0),
            aromatic: Some(pa.aromatic),
            charge: (pa.charge != 0).then_some(pa.charge),
            // 角括弧の H 指定のみ制約になる ([nH] → 1)。
            // 指定なし角括弧は Some(0) になるが、対象は [se] 等の
            // H を持ち得ない原子なので実害はない
            n_h: pa.explicit_h,
            wildcard: pa.symbol == "*",
        });
    }
    let mut adj = vec![Vec::new(); n];
    for b in &parsed.bonds {
        let k = match b.kind {
            BondKind::Elided => QBondKind::SingleOrAromatic,
            BondKind::Single | BondKind::Up | BondKind::Down => QBondKind::Exact(0),
            BondKind::Aromatic => QBondKind::Exact(1),
            BondKind::Double => QBondKind::Exact(2),
            BondKind::Triple => QBondKind::Exact(3),
            BondKind::Quadruple => QBondKind::Exact(4),
        };
        adj[b.a].push((b.b, k));
        adj[b.b].push((b.a, k));
    }
    Ok((atoms, adj))
}

fn atom_compat(qa: &QAtom, t: &TargetView, ti: usize) -> bool {
    if !qa.wildcard && qa.atomic_num != t.atomic_num[ti] {
        return false;
    }
    if let Some(ar) = qa.aromatic {
        if !qa.wildcard && ar != t.aromatic[ti] {
            return false;
        }
    }
    if let Some(c) = qa.charge {
        if c != t.charge[ti] {
            return false;
        }
    }
    if let Some(h) = qa.n_h {
        if h != t.n_h[ti] {
            return false;
        }
    }
    true
}

fn bond_compat(qk: QBondKind, target_oc: u8) -> bool {
    match qk {
        QBondKind::Exact(k) => k == target_oc,
        QBondKind::SingleOrAromatic => target_oc == 0 || target_oc == 1,
    }
}

/// Mol クエリでの全マッチ列挙 (uniquify=False 相当)。タプルはクエリ原子順。
pub fn substruct_matches(target: &MoleculeGraph, query: &MoleculeGraph) -> Vec<Vec<usize>> {
    let (q_atoms, q_adj) = build_mol_query(query);
    run_vf2(target, &q_atoms, &q_adj)
}

/// SMILES 文字列を SMARTS として解釈した全マッチ列挙。
pub fn substruct_matches_smarts(
    target: &MoleculeGraph,
    pattern: &str,
) -> Result<Vec<Vec<usize>>, ChemError> {
    let (q_atoms, q_adj) = build_smarts_query(pattern)?;
    Ok(run_vf2(target, &q_atoms, &q_adj))
}

fn run_vf2(
    target: &MoleculeGraph,
    q_atoms: &[QAtom],
    q_adj: &[Vec<(usize, QBondKind)>],
) -> Vec<Vec<usize>> {
    let t = build_target_view(target);
    let nq = q_atoms.len();
    if nq == 0 || nq > t.n {
        return Vec::new();
    }

    // クエリ原子の探索順: 既出原子に隣接するものを優先 (連結クエリなら常に隣接)
    let mut order = Vec::with_capacity(nq);
    let mut placed = vec![false; nq];
    while order.len() < nq {
        let next = (0..nq)
            .filter(|&i| !placed[i])
            .find(|&i| q_adj[i].iter().any(|&(j, _)| placed[j]))
            .or_else(|| (0..nq).find(|&i| !placed[i]))
            .expect("unplaced atom exists");
        placed[next] = true;
        order.push(next);
    }

    const MAX_MATCHES: usize = 100_000;
    let mut mapping = vec![usize::MAX; nq]; // query atom → target atom
    let mut used = vec![false; t.n];
    let mut results = Vec::new();

    #[allow(clippy::too_many_arguments)]
    fn backtrack(
        depth: usize,
        order: &[usize],
        q_atoms: &[QAtom],
        q_adj: &[Vec<(usize, QBondKind)>],
        t: &TargetView,
        mapping: &mut [usize],
        used: &mut [bool],
        results: &mut Vec<Vec<usize>>,
    ) {
        if results.len() >= MAX_MATCHES {
            return;
        }
        if depth == order.len() {
            results.push(mapping.to_vec());
            return;
        }
        let qi = order[depth];
        // マップ済みのクエリ隣接があればその像の隣接だけを候補にする
        let anchor = q_adj[qi]
            .iter()
            .find(|&&(j, _)| mapping[j] != usize::MAX)
            .map(|&(j, _)| mapping[j]);
        let candidates: Vec<usize> = match anchor {
            Some(ta) => t.adj[ta].iter().map(|&(v, _)| v).collect(),
            None => (0..t.n).collect(),
        };
        'cand: for ti in candidates {
            if used[ti] || !atom_compat(&q_atoms[qi], t, ti) {
                continue;
            }
            // qi とマップ済み隣接の全結合がターゲットに存在し互換であること
            for &(qj, qk) in &q_adj[qi] {
                let tj = mapping[qj];
                if tj == usize::MAX {
                    continue;
                }
                let Some(&(_, toc)) = t.adj[ti].iter().find(|&&(v, _)| v == tj) else {
                    continue 'cand;
                };
                if !bond_compat(qk, toc) {
                    continue 'cand;
                }
            }
            mapping[qi] = ti;
            used[ti] = true;
            backtrack(depth + 1, order, q_atoms, q_adj, t, mapping, used, results);
            mapping[qi] = usize::MAX;
            used[ti] = false;
        }
    }
    backtrack(
        0,
        &order,
        q_atoms,
        q_adj,
        &t,
        &mut mapping,
        &mut used,
        &mut results,
    );
    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    fn mol_matches(target: &str, query: &str) -> Vec<Vec<usize>> {
        let t = build_molecule_graph(target).expect("target");
        let q = build_molecule_graph(query).expect("query");
        let mut m = substruct_matches(&t, &q);
        m.sort();
        m
    }

    fn smarts_matches(target: &str, pattern: &str) -> Vec<Vec<usize>> {
        let t = build_molecule_graph(target).expect("target");
        let mut m = substruct_matches_smarts(&t, pattern).expect("pattern");
        m.sort();
        m
    }

    #[test]
    fn benzene_in_naphthalene() {
        // 12 自己同型 × 2 環 = 24 マッチ
        assert_eq!(mol_matches("c1ccc2ccccc2c1", "c1ccccc1").len(), 24);
    }

    #[test]
    fn charge_semantics() {
        // 中性クエリは荷電ターゲットにマッチ、逆は不可 (RDKit 実測)
        assert_eq!(mol_matches("CC(=O)[O-]", "CC(=O)O").len(), 1);
        assert_eq!(mol_matches("CC(=O)O", "CC(=O)[O-]").len(), 0);
    }

    #[test]
    fn wildcard_semantics() {
        // mol モードのダミーは何にもマッチしない
        assert_eq!(mol_matches("CC", "*C").len(), 0);
        // SMARTS モードでは任意原子
        assert_eq!(smarts_matches("CC", "*C").len(), 2);
    }

    #[test]
    fn smarts_nh_constraint() {
        // [nH] は総 H 数を制約する: N-メチルインドールには不一致
        let indole_q = "c1ccc2[nH]ccc2c1";
        assert_eq!(smarts_matches("Cn1ccc2ccccc21", indole_q).len(), 0);
        assert!(!smarts_matches("c1ccc2[nH]ccc2c1", indole_q).is_empty());
        // mol モードは H を照合しないので N-メチルにもマッチ
        assert!(!mol_matches("Cn1ccc2ccccc21", indole_q).is_empty());
        // ケクレ化不能なパターンも SMARTS としては有効
        assert!(smarts_matches("Cn1ccc2ccccc21", "c1ccc2nccc2c1").len() == 1);
    }

    #[test]
    fn aliphatic_vs_aromatic() {
        // ピペリジンクエリはピリジンにマッチしない (結合次数が異なる)
        assert_eq!(mol_matches("c1ccncc1", "C1CCNCC1").len(), 0);
        // SMARTS モードでは大文字 = 脂肪族制約
        assert_eq!(smarts_matches("c1ccncc1", "C1CCNCC1").len(), 0);
    }
}
