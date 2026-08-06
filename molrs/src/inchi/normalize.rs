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

/// 重原子の連結成分 id (重原子 idx → 成分 id)。電荷の中性化は成分ごとに
/// 行うため ([`neutralize`])、`inchi::number::connected_components` を待たずに
/// ここで求める必要がある。
fn heavy_components(g: &MoleculeGraph, n_heavy: usize) -> Vec<usize> {
    let mut comp = vec![usize::MAX; n_heavy];
    let mut next = 0usize;
    for start in 0..n_heavy {
        if comp[start] != usize::MAX {
            continue;
        }
        let mut stack = vec![start];
        comp[start] = next;
        while let Some(a) = stack.pop() {
            for &nb in &g.adjacency[a] {
                if nb < n_heavy && comp[nb] == usize::MAX {
                    comp[nb] = next;
                    stack.push(nb);
                }
            }
        }
        next += 1;
    }
    comp
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

    // 隣接に逆符号の電荷を持つ原子 (イリド/N-オキシド/ニトロ/アジド等の
    // 電荷分離) はプロトン化しない — InChI は共有結合の中性形で扱う。
    let has_opposite_charged_neighbor = |i: usize| {
        let ci = g.atoms[i].formal_charge;
        g.adjacency[i].iter().any(|&nb| {
            (g.atoms[nb].formal_charge < 0) != (ci < 0) && g.atoms[nb].formal_charge != 0
        })
    };

    // 成分ごとのプロトン移動予算 (I31)。
    //
    // InChI が中性化するのは「各原子」ではなく「**各成分の正味電荷**」。
    // 分子内塩 (アセチルカルニチン `CC(=O)OC(CC(=O)[O-])C[N+](C)(C)C` の
    // ように四級 N+ とカルボキシラートが同じ成分にある) は正味 0 なので
    // **プロトンを一切動かさない** — 実 InChI は `C9H17NO4` (H 17 個) で
    // `/q` も `/p` も付かない。原子ごとに中性化すると O- をプロトン化して
    // H 18 個になり、四級 N+ が中和できず `/q+1/p-1` が付いてしまう。
    //
    // 正味 +1 の側 (`CC(=O)OC(CC(=O)O)C[N+](C)(C)C`) は COOH から 1 個
    // 外して正味 0 にし `/p+1`。基準状態はどちらも同じ双性イオンになる。
    //
    // 酢酸ナトリウム `[Na+].[O-]C(=O)C` のような塩は成分が分かれているので
    // 従来どおり: 酢酸イオン成分は正味 -1 → プロトン付加で `/p-1`、Na 成分は
    // 正味 +1 で外せる H がなく `/q+1`。
    let comp_of = heavy_components(g, n_heavy);
    let n_comps = comp_of.iter().copied().max().map_or(0, |m| m + 1);
    let mut budget = vec![0i32; n_comps];
    for i in 0..n_heavy {
        budget[comp_of[i]] += g.atoms[i].formal_charge as i32;
    }

    for i in 0..n_heavy {
        let a = &g.atoms[i];
        let h = cur_h(i);
        let ch = a.formal_charge as i32;
        let b = &mut budget[comp_of[i]];
        if ch != 0 && has_opposite_charged_neighbor(i) {
            // 電荷分離 (ニトロ・N-オキシド等) → 触らない
            final_h[i] = h;
            new_charge[i] = a.formal_charge;
        } else if ch < 0 && is_protonatable(&a.symbol) && *b < 0 {
            // 負電荷 → プロトン付加で中性化 (成分の正味電荷を超えない範囲で)
            let add = (-ch).min(-*b);
            n_add += add;
            *b += add;
            final_h[i] = h + add;
            new_charge[i] = (ch + add) as i8;
        } else if ch > 0 && is_deprotonatable(&a.symbol) && h > 0 && *b > 0 {
            // プロトン付き陽イオン → 除去で中性化
            let rem = ch.min(h).min(*b);
            n_remove += rem;
            *b -= rem;
            final_h[i] = h - rem;
            new_charge[i] = (ch - rem) as i8;
        } else {
            // 中性化不能 (四級 N・金属など)、または成分が既に正味中性
            final_h[i] = h;
            new_charge[i] = a.formal_charge;
        }
    }

    // 第 2 パス (I31): 成分にまだ正味の陽電荷が残っていて、外せる**陽イオンの**
    // H が無い場合は、**中性の酸点** (カルボン酸等の O-H) からプロトンを外して
    // 釣り合わせる。
    //
    // カルニチン `C[N+](C)(C)CC(CC(=O)O)O` がこれで、四級 N+ は脱プロトン
    // できないが実 InChI は COOH から 1 個外した双性イオンを基準にして
    // `C7H15NO3/…/p+1` を出す。外さないと `/q+1` になってしまう。
    //
    // **酸性の O-H/S-H だけ**が対象。隣の重原子 (C/N/S/P) が O/S へ二重結合を
    // 持つもの — カルボン酸・スルホン酸・リン酸・ヒドロキサム酸など。
    //
    // 単なるアルコールを外してはいけない (I35 で判明)。コリン
    // `C[N+](C)(C)CCO` の実 InChI は `C5H14NO/…/q+1` で、OH の H を保ったまま
    // `/q` に残す。ここを「任意の O-H」にしていると `/p+1` になってしまう。
    let is_acidic_oh = |i: usize| -> bool {
        g.adjacency[i].iter().any(|&c| {
            matches!(g.atoms[c].symbol.as_str(), "C" | "N" | "S" | "P")
                && g.bonds.iter().enumerate().any(|(bi, b)| {
                    (b.begin_idx == c || b.end_idx == c)
                        && g.kekule_bond_orders[bi] == 2.0
                        && matches!(
                            g.atoms[if b.begin_idx == c {
                                b.end_idx
                            } else {
                                b.begin_idx
                            }]
                            .symbol
                            .as_str(),
                            "O" | "S"
                        )
                })
        })
    };
    for (c, left) in budget.iter_mut().enumerate().take(n_comps) {
        while *left > 0 {
            let pick = (0..n_heavy)
                .filter(|&i| {
                    comp_of[i] == c
                        && new_charge[i] == 0
                        && final_h[i] > 0
                        && matches!(g.atoms[i].symbol.as_str(), "O" | "S")
                        && is_acidic_oh(i)
                })
                .min_by_key(|&i| (g.atoms[i].symbol != "O", i));
            let Some(i) = pick else { break };
            final_h[i] -= 1;
            new_charge[i] = -1;
            n_remove += 1;
            *left -= 1;
        }
    }
    // q は最終的な電荷の総和 (2 パス目で変わるのでここで数え直す)
    let q: i32 = new_charge.iter().map(|&c| c as i32).sum();

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
