//! 電荷正規化 (`/q`・`/p` 層、I10)。
//!
//! InChI は分子を可能な限り中性化してから骨格層を計算する:
//! - 負電荷のヘテロ酸点 (カルボキシラート等) はプロトン付加で中性化 → `/p` 減
//! - プロトン付き塩基点 (アンモニウム等) はプロトン除去で中性化 → `/p` 増
//! - 中性化できない電荷 (四級 N・金属イオン) は残余電荷として `/q`
//!
//! [`neutralize`] は中性化したグラフのクローンと (q, p) を返す。重原子の
//! インデックスは保存する (0..n_heavy) ため `parsed`/`parser_to_graph` と
//! CIP ランク・立体タグは有効なまま使える。

use crate::graph::{AtomInfo, BondInfo, MoleculeGraph};
use std::collections::HashMap;

/// 負電荷を中性化 (プロトン付加) できる元素。ハロゲン化物イオンも含む
/// (InChI は [Cl-] を HCl/p-1 とする)。
fn is_protonatable(sym: &str) -> bool {
    matches!(sym, "N" | "O" | "S" | "Se" | "Te" | "F" | "Cl" | "Br" | "I")
}

/// 陽イオンから脱プロトンして中性化できる元素 (塩基点)。ハロゲンは対象外。
fn is_deprotonatable(sym: &str) -> bool {
    matches!(sym, "N" | "O" | "S" | "Se" | "Te")
}

/// 中性化したグラフと (q, p)。q = 残余電荷合計、p = 除去 - 付加 プロトン数。
pub(crate) fn neutralize(g: &MoleculeGraph) -> (MoleculeGraph, i32, i32) {
    let n_heavy = g.atoms.iter().filter(|a| a.symbol != "H").count();
    // 重原子が先頭に連続していることを前提 (build_molecule_graph の不変条件)
    let heavy_contiguous = (0..n_heavy).all(|i| g.atoms[i].symbol != "H");
    if !heavy_contiguous {
        return (g.clone(), 0, 0);
    }

    // 各重原子の現在の H 数
    let cur_h = |i: usize| {
        g.adjacency[i]
            .iter()
            .filter(|&&x| g.atoms[x].symbol == "H")
            .count() as i32
    };

    let mut new_charge = vec![0i8; n_heavy];
    let mut final_h = vec![0i32; n_heavy];
    let mut n_add = 0i32;
    let mut n_remove = 0i32;
    let mut q = 0i32;

    // 隣接に逆符号の電荷を持つ原子 (イリド/N-オキシド/ニトロ/アジド等の
    // 電荷分離) はプロトン化しない — InChI は共有結合の中性形で扱う。
    let has_opposite_charged_neighbor = |i: usize| {
        let ci = g.atoms[i].formal_charge;
        g.adjacency[i].iter().any(|&nb| {
            (g.atoms[nb].formal_charge < 0) != (ci < 0) && g.atoms[nb].formal_charge != 0
        })
    };

    for i in 0..n_heavy {
        let a = &g.atoms[i];
        let h = cur_h(i);
        let ch = a.formal_charge as i32;
        if ch != 0 && has_opposite_charged_neighbor(i) {
            // 電荷分離 (ネット中性の zwitterion) → 触らない
            final_h[i] = h;
            new_charge[i] = a.formal_charge;
            q += ch;
        } else if ch < 0 && is_protonatable(&a.symbol) {
            // 負電荷 → プロトン付加で中性化
            let add = -ch;
            n_add += add;
            final_h[i] = h + add;
            new_charge[i] = 0;
        } else if ch > 0 && is_deprotonatable(&a.symbol) && h > 0 {
            // プロトン付き陽イオン → 除去で中性化 (除去可能な H まで)
            let rem = ch.min(h);
            n_remove += rem;
            final_h[i] = h - rem;
            let residual = ch - rem;
            new_charge[i] = residual as i8;
            q += residual;
        } else {
            // 中性化不能 (四級 N・金属など)
            final_h[i] = h;
            new_charge[i] = a.formal_charge;
            q += ch;
        }
    }

    let p = n_remove - n_add;
    if n_add == 0 && n_remove == 0 && q == g.atoms.iter().map(|a| a.formal_charge as i32).sum() {
        // 変化なし (かつ元から電荷なし) ならクローン省略のため元を返す
        if q == 0 {
            return (g.clone(), 0, 0);
        }
    }

    // 中性化グラフを再構築: 重原子 (電荷調整) + final_h に基づく H ノード
    let mut atoms: Vec<AtomInfo> = Vec::with_capacity(n_heavy);
    for (i, &nc) in new_charge.iter().enumerate() {
        let mut a = g.atoms[i].clone();
        a.formal_charge = nc;
        atoms.push(a);
    }
    // 重原子-重原子結合を保持
    let mut bonds: Vec<BondInfo> = Vec::new();
    let mut kekule: Vec<f64> = Vec::new();
    for (bi, b) in g.bonds.iter().enumerate() {
        if b.begin_idx < n_heavy && b.end_idx < n_heavy {
            bonds.push(b.clone());
            kekule.push(g.kekule_bond_orders[bi]);
        }
    }
    // H ノードを再付加
    for (i, &fh) in final_h.iter().enumerate() {
        for _ in 0..fh {
            let h_idx = atoms.len();
            atoms.push(AtomInfo {
                idx: h_idx,
                symbol: "H".into(),
                atomic_num: 1,
                is_aromatic: false,
                in_ring: false,
                num_hs: 0,
                chiral_tag: None,
                formal_charge: 0,
            });
            bonds.push(BondInfo {
                begin_idx: i,
                end_idx: h_idx,
                bond_order: 1.0,
                stereo: None,
            });
            kekule.push(1.0);
        }
    }
    // 金属切断で生じた孤立 H (どの重原子にも結合しない) は独立成分として
    // 保持する必要がある (`[BiH3]` → `Bi.3H`)。上のループは重原子に結合した
    // H しか再生成しないため、ここで補う (I20)。
    let mut h_remap: HashMap<usize, usize> = HashMap::new();
    for old in 0..g.atoms.len() {
        let is_lone_h = g.atoms[old].symbol == "H"
            && !g.adjacency[old].iter().any(|&nb| g.atoms[nb].symbol != "H");
        if !is_lone_h {
            continue;
        }
        let h_idx = atoms.len();
        h_remap.insert(old, h_idx);
        atoms.push(AtomInfo {
            idx: h_idx,
            symbol: "H".into(),
            atomic_num: 1,
            is_aromatic: false,
            in_ring: false,
            num_hs: 0,
            chiral_tag: None,
            formal_charge: 0,
        });
    }
    // 孤立 H 同士の結合 (水素分子 `[H][H]`) は成分としてまとめる必要があるので
    // 保持する。金属水素化物由来の H は互いに結合していないので何も増えない。
    for (bi, b) in g.bonds.iter().enumerate() {
        if let (Some(&i), Some(&j)) = (h_remap.get(&b.begin_idx), h_remap.get(&b.end_idx)) {
            bonds.push(BondInfo {
                begin_idx: i,
                end_idx: j,
                bond_order: b.bond_order,
                stereo: None,
            });
            kekule.push(g.kekule_bond_orders[bi]);
        }
    }
    // idx を振り直し (重原子は不変、H は末尾)
    for (i, a) in atoms.iter_mut().enumerate() {
        a.idx = i;
    }
    let mut adjacency = vec![Vec::new(); atoms.len()];
    let mut bond_orders = HashMap::new();
    for b in &bonds {
        adjacency[b.begin_idx].push(b.end_idx);
        adjacency[b.end_idx].push(b.begin_idx);
        bond_orders.insert(
            (b.begin_idx.min(b.end_idx), b.begin_idx.max(b.end_idx)),
            b.bond_order,
        );
    }

    let ng = MoleculeGraph {
        atoms,
        bonds,
        adjacency,
        bond_orders,
        ring_atom_sets: g.ring_atom_sets.clone(),
        kekule_bond_orders: kekule,
        parsed: g.parsed.clone(),
        parser_to_graph: g.parser_to_graph.clone(),
    };
    (ng, q, p)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    fn qp(smiles: &str) -> (i32, i32, String) {
        let g = build_molecule_graph(smiles).unwrap();
        let (ng, q, p) = neutralize(&g);
        (q, p, super::super::formula::formula_layer(&ng))
    }

    #[test]
    fn carboxylate_adds_h() {
        let (q, p, f) = qp("CC(=O)[O-]");
        assert_eq!((q, p), (0, -1));
        assert_eq!(f, "C2H4O2"); // 中性酸の式
    }

    #[test]
    fn dicarboxylate() {
        let (q, p, _) = qp("O=C([O-])[O-]");
        assert_eq!((q, p), (0, -2));
    }

    #[test]
    fn ammonium_removes_h() {
        let (q, p, f) = qp("C[NH3+]");
        assert_eq!((q, p), (0, 1));
        assert_eq!(f, "CH5N"); // 中性アミン
    }

    #[test]
    fn quaternary_keeps_charge() {
        let (q, p, f) = qp("C[N+](C)(C)C");
        assert_eq!((q, p), (1, 0));
        assert_eq!(f, "C4H12N");
    }

    #[test]
    fn neutral_unchanged() {
        let (q, p, f) = qp("CCO");
        assert_eq!((q, p), (0, 0));
        assert_eq!(f, "C2H6O");
    }
}
