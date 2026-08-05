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
/// 実 InChI の多成分順序 (I29):
/// **炭素数の降順 → 重原子数の降順 → H 数の降順 → 式の辞書順昇順**。
///
/// - 炭素数**降順**が主キー。PubChem 実データ 863 件の多成分分子で
///   例外なく成立する。安息香酸カリウム `C7H6O2.K` も硫酸ナトリウム
///   `2Na.H2O4S` (どちらも炭素数 0) もこれで説明でき、`C6H5.C2H4O2.Hg` の
///   ように炭素数の大きい成分が先に来るのが本質。
/// - 重原子数**降順**が第 2 キー: `C10H11O.C5H5.Fe`、`C9H12O.3CO.Fe`。
/// - H 数降順は重原子数が並んだときの決定打: `CH3.ClH.Hg` の ClH (H 1 個) が
///   Hg (H 0 個) より先。
///
/// 金属から切り離された孤立 H 成分は重原子 0 個でこの比較に載らないため、
/// 呼び出し側 ([`formula_layer`]) が常に末尾へ追加する (`Bi.3H`)。
///
/// # 既知の残差 (無機塩)
///
/// アルカリ金属塩など「単原子カチオン + 多原子アニオン」で、実 InChI は
/// カチオンを先に置くことがある (`2Na.H2O4S`、`5Na.H3O4P.H2O3S`、
/// `Cu.N2O4.2NO3`)。一方 `FH.O3Si.2Zn` は単原子の Zn が最後に来るので
/// 「単原子金属を先頭」という規則では説明できず、電荷層との関係も含めて
/// 未解明。PubChem 863 件中 33 件 (3.8%) がこの系統で残る。
///
/// 旧実装は重原子数**昇順**だったが、これはリポジトリ内コーパスの
/// 多成分 32 例 (`2Na.H2O4S` 系が多数) に過適合したもので、PubChem 実データ
/// では 33.7% しか再現できなかった。本規則は 96.2%。
pub(crate) fn component_sort_key(
    g: &MoleculeGraph,
    atoms: &[usize],
) -> (
    std::cmp::Reverse<usize>,
    std::cmp::Reverse<usize>,
    std::cmp::Reverse<usize>,
    String,
) {
    use std::cmp::Reverse;
    let n_c = atoms.iter().filter(|&&a| g.atoms[a].symbol == "C").count();
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
        Reverse(n_c),
        Reverse(atoms.len()),
        Reverse(h),
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
