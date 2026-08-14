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

/// [`component_sort_key`] の主キー 1 要素分: Hill 式 (H を除く) を要素記号
/// ごとに分解したときの 1 トークン。
///
/// 実 InChI `ichimake.c::GetElementAndCount`/`CompareHillFormulasNoH` の
/// 移植: 炭素は常に他のどの元素よりも小さい特別扱い (`Carbon`)、式が尽きた
/// 側は他のどの元素よりも大きい番兵 (`End`) — 「式が短い方が (辞書順で)
/// 後ろ」という C 側のコメント通りの挙動になる。
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum ElemKey {
    Carbon,
    Elem(String),
    End,
}

/// Hill 式文字列 (H を含む) を `CompareHillFormulasNoH` と同じ規則で
/// トークン列に分解する: 元素記号+個数を先頭から読み、H だけは完全に
/// 読み飛ばす (比較にも位置にも影響しない)。個数は降順 (`Reverse`) — 同じ
/// 元素なら原子数が多い方が先。
fn hill_no_h_key(formula: &str) -> Vec<(ElemKey, std::cmp::Reverse<u64>)> {
    let bytes = formula.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        let start = i;
        i += 1;
        if i < bytes.len() && bytes[i].is_ascii_lowercase() {
            i += 1;
        }
        let sym = &formula[start..i];
        let dstart = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        let count: u64 = if dstart == i {
            1
        } else {
            formula[dstart..i].parse().unwrap_or(1)
        };
        if sym == "H" {
            continue;
        }
        let key = if sym == "C" {
            ElemKey::Carbon
        } else {
            ElemKey::Elem(sym.to_string())
        };
        out.push((key, std::cmp::Reverse(count)));
    }
    out.push((ElemKey::End, std::cmp::Reverse(0)));
    out
}

/// 成分の並び順キー (I20、I53 で全面書き直し)。
///
/// 実 InChI `ichimake.c::CompINChI2` の主キーそのままの移植:
/// **Hill 式 (H を除く) のトークン単位比較 → 重原子数の降順 → H 数の降順 →
/// 式の辞書順昇順**。
///
/// トークン単位比較 (`hill_no_h_key`) が本質: 炭素含有成分が常に先頭
/// (炭素は他のどの元素よりも「小さい」特別値)、以降は残り元素を出現順
/// (Hill 式なのでアルファベット順) に記号→個数の順で比較し、記号が違えば
/// 記号のアルファベット順(小さい方が先)、同じ記号なら個数の多い方が先。
/// H は完全に読み飛ばす (位置にも比較にも数えない) — `2H3N` (アンモニウム)
/// が `H2O3S2` より先に来るのはこれで説明できる (H を除くと N と O3S2、
/// N < O)。式が短く尽きた側は「無限に大きい元素」が続くとみなし後ろに回る
/// (`CH2O` が `CH2` より先)。
///
/// 旧実装 (炭素数の単純な降順が主キー) は C ソースを読まずに PubChem
/// 863 件から逆算したもので、`C10H6ClNO` vs `C10H7NO2` のような
/// 同炭素数・同重原子数の成分順序を取り違えていた。
pub(crate) type ComponentSortKey = (
    Vec<(ElemKey, std::cmp::Reverse<u64>)>,
    std::cmp::Reverse<usize>,
    std::cmp::Reverse<usize>,
    String,
);

pub(crate) fn component_sort_key(g: &MoleculeGraph, atoms: &[usize]) -> ComponentSortKey {
    use std::cmp::Reverse;
    let formula = component_formula(g, atoms);
    let h = atoms
        .iter()
        .map(|&a| {
            g.adjacency[a]
                .iter()
                .filter(|&&nb| g.atoms[nb].symbol == "H")
                .count()
        })
        .sum::<usize>();
    (hill_no_h_key(&formula), Reverse(atoms.len()), Reverse(h), formula)
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
