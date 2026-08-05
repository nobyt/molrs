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

/// 成分の並び順キー (I20)。
///
/// 実 InChI の多成分順序をコーパス 32 例から導出した規則:
/// **炭素を含む成分が先 → 重原子数の昇順 → H 数の降順 → 式の辞書順昇順**。
///
/// - 炭素優先は Hill 式の思想と同じ: 硫酸ナトリウムが `2Na.H2O4S` (Na が先)
///   なのに安息香酸カリウムが `C7H6O2.K` (有機が先) になるのは、後者だけが
///   炭素を含むため。重原子数だけでは説明できない。
/// - 重原子数**昇順**: `2Na.H2O4S` は Na (重原子 1) が H2O4S (重原子 5) より
///   先に来る。
/// - H 数降順は重原子数が並んだときの決定打: `CH3.ClH.Hg` の ClH (H 1 個) が
///   Hg (H 0 個) より先。
///
/// 金属から切り離された孤立 H 成分は重原子 0 個でこの比較に載らないため、
/// 呼び出し側 ([`formula_layer`]) が常に末尾へ追加する (`Bi.3H`)。
pub(crate) fn component_sort_key(
    g: &MoleculeGraph,
    atoms: &[usize],
) -> (bool, usize, std::cmp::Reverse<usize>, String) {
    let has_carbon = atoms.iter().any(|&a| g.atoms[a].symbol == "C");
    let h = atoms
        .iter()
        .map(|&a| {
            g.adjacency[a]
                .iter()
                .filter(|&&nb| g.atoms[nb].symbol == "H")
                .count()
        })
        .sum::<usize>();
    (
        !has_carbon,
        atoms.len(),
        std::cmp::Reverse(h),
        component_formula(g, atoms),
    )
}

/// 1 成分の元素数え上げ → Hill 式文字列。
pub(crate) fn component_formula(g: &MoleculeGraph, atoms: &[usize]) -> String {
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
    let mut comps = connected_components(g);
    // 成分順序は c/h/q 層と共通の規則で決める (I20)
    comps.sort_by_key(|atoms| component_sort_key(g, atoms));
    let mut formulas: Vec<String> = comps
        .iter()
        .map(|atoms| component_formula(g, atoms))
        .collect();
    // 重原子を含まない H だけの成分 (金属から切り離された H、水素分子) は
    // 重原子比較の対象外なので常に末尾に置く
    for size in super::disconnect::hydrogen_component_sizes(g) {
        formulas.push(super::disconnect::hydrogen_component_formula(size));
    }

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
