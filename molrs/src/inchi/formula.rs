//! InChI 式層 (Hill 式) の生成 (I2)。
//!
//! Hill システム: 炭素があれば C, H, 残りをアルファベット順。炭素がなければ
//! H を含め全元素をアルファベット順。多成分は各成分の式をソートし、同一式は
//! 数係数でまとめる (`3H2O`, `2C2H6O`, `ClH.Na`)。
//!
//! 注意: InChI は式を計算する前に電荷正規化 (プロトン移動) を行う
//! (`[Na+].[Cl-]` → `ClH.Na`)。本モジュールは与えられたグラフの組成を
//! そのまま数えるため、正規化で組成が変わる分子は normalize.rs (I4) を
//! 通した後に呼ぶこと。

use crate::graph::MoleculeGraph;

use super::number::connected_components;

/// 1 成分の元素数え上げ → Hill 式文字列。
fn component_formula(g: &MoleculeGraph, atoms: &[usize]) -> String {
    use std::collections::BTreeMap;
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for &a in atoms {
        *counts.entry(g.atoms[a].symbol.as_str()).or_insert(0) += 1;
        // 重原子に結合した明示 H ノードを数える
        for &nb in &g.adjacency[a] {
            if g.atoms[nb].symbol == "H" {
                *counts.entry("H").or_insert(0) += 1;
            }
        }
    }
    let mut out = String::new();
    let mut emit = |sym: &str, n: usize| {
        out.push_str(sym);
        if n > 1 {
            out.push_str(&n.to_string());
        }
    };
    let has_c = counts.contains_key("C");
    if has_c {
        if let Some(&n) = counts.get("C") {
            emit("C", n);
        }
        if let Some(&n) = counts.get("H") {
            emit("H", n);
        }
        for (&sym, &n) in &counts {
            if sym != "C" && sym != "H" {
                emit(sym, n);
            }
        }
    } else {
        // 炭素なし: H を含め全てアルファベット順 (BTreeMap の順)
        for (&sym, &n) in &counts {
            emit(sym, n);
        }
    }
    out
}

/// 分子全体の式層 (先頭の `InChI=1S/` と最初の `/` の間の部分)。
pub(crate) fn formula_layer(g: &MoleculeGraph) -> String {
    let comps = connected_components(g);
    // 各成分の式を計算し、辞書順にソート
    let mut formulas: Vec<String> = comps
        .iter()
        .map(|atoms| component_formula(g, atoms))
        .collect();
    formulas.sort();

    // 連続する同一式を数係数でまとめる
    let mut out = String::new();
    let mut i = 0;
    while i < formulas.len() {
        let mut j = i + 1;
        while j < formulas.len() && formulas[j] == formulas[i] {
            j += 1;
        }
        let count = j - i;
        if !out.is_empty() {
            out.push('.');
        }
        if count > 1 {
            out.push_str(&count.to_string());
        }
        out.push_str(&formulas[i]);
        i = j;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    fn formula(smiles: &str) -> String {
        formula_layer(&build_molecule_graph(smiles).unwrap())
    }

    #[test]
    fn single_component_hill() {
        assert_eq!(formula("CC(=O)O"), "C2H4O2");
        assert_eq!(formula("c1ccccc1"), "C6H6");
        assert_eq!(formula("FC(F)F"), "CHF3");
        assert_eq!(formula("c1cc[nH]c1"), "C4H5N");
        assert_eq!(formula("C"), "CH4");
    }

    #[test]
    fn no_carbon_alphabetical() {
        assert_eq!(formula("O"), "H2O");
        assert_eq!(formula("[NH4+]"), "H4N"); // 電荷は q 層、式は組成のみ
    }

    #[test]
    fn multi_component_collapse_and_sort() {
        assert_eq!(formula("O.O.O"), "3H2O");
        assert_eq!(formula("CCO.CCO"), "2C2H6O");
    }
}
