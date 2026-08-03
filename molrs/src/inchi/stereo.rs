//! InChI 立体層 `/b`・`/t`・`/m`・`/s` の生成 (I6/I7)。
//!
//! molrs のグラフは E/Z (CIP 基準) と R/S (CIP) を持つが、InChI の立体記述は
//! **正準番号** 基準のパリティで表す。ここでは正準番号 (number.rs) と CIP
//! ランク (stereo.rs) を使って InChI パリティに変換する。
//!
//! v1 の対象: 非環・単純な四面体中心と二重結合。累積二重結合・sp3 以外の
//! 立体・立体中心の相互依存 (擬似不斉) は未対応。

use crate::graph::MoleculeGraph;
use crate::stereo::cip_ranks;

use super::layers::Component;

/// 原子 → 正準番号 (1..=n) を単一成分の Component から引く。
fn canon_of(comp: &Component, orig: usize) -> usize {
    comp.num[orig]
}

/// 二重結合の一端 `end` (相手 `other`) の非水素・非相手隣接一覧。
fn ez_neighbors(g: &MoleculeGraph, end: usize, other: usize) -> Vec<usize> {
    g.adjacency[end]
        .iter()
        .copied()
        .filter(|&nb| nb != other && g.atoms[nb].symbol != "H")
        .collect()
}

/// `/b` 層 (先頭の `b` は含めない)。空なら空文字列。
pub(crate) fn double_bond_layer(g: &MoleculeGraph, comp: &Component) -> String {
    let ranks = cip_ranks(g);
    let mut entries: Vec<(usize, usize, char)> = Vec::new();

    for b in &g.bonds {
        let Some(ez) = b.stereo else { continue };
        let (a, c) = (b.begin_idx, b.end_idx);
        // 成分内かつ番号付け済み
        if comp.num.get(a).copied().unwrap_or(0) == 0 || comp.num.get(c).copied().unwrap_or(0) == 0
        {
            continue;
        }
        let na = ez_neighbors(g, a, c);
        let nc = ez_neighbors(g, c, a);
        if na.is_empty() || nc.is_empty() {
            continue;
        }
        // CIP 最高隣接 (molrs の E/Z の基準)
        let cip_hi = |ns: &[usize]| *ns.iter().max_by_key(|&&x| ranks[x]).unwrap();
        // 正準番号最高隣接 (InChI の基準)
        let canon_hi = |ns: &[usize]| *ns.iter().max_by_key(|&&x| canon_of(comp, x)).unwrap();
        let flip_a = cip_hi(&na) != canon_hi(&na);
        let flip_c = cip_hi(&nc) != canon_hi(&nc);
        let opposite_cip = ez == 'E'; // E = CIP 最高が反対側
        let opposite_canon = opposite_cip ^ flip_a ^ flip_c;
        let parity = if opposite_canon { '+' } else { '-' };

        let (ca, cc) = (canon_of(comp, a), canon_of(comp, c));
        let (hi, lo) = (ca.max(cc), ca.min(cc));
        entries.push((hi, lo, parity));
    }

    // 未定義の立体源性二重結合 (`?`): 定義済みの立体二重結合が 1 つでも
    // ある場合のみ、立体源性だが SMILES で構成が指定されていない二重結合
    // (代表例: イミン/ヒドラゾンの C=N) を `hi-lo?` として併記する
    // (実 InChI/RDKit の挙動。未定義のみの分子では /b 層自体が省略される)。
    if !entries.is_empty() {
        let listed: std::collections::HashSet<(usize, usize)> =
            entries.iter().map(|&(hi, lo, _)| (hi, lo)).collect();
        // 可動 H 群のメンバーを端点に持つ二重結合 (アミジン/グアニジンの
        // C=N 等) は、二重結合の位置自体が互変異性で動くため対象外。
        let tgroup = super::number::tautomer_group_members(g);
        for (bi, b) in g.bonds.iter().enumerate() {
            if b.stereo.is_some() || g.kekule_bond_orders[bi] != 2.0 {
                continue;
            }
            let (a, c) = (b.begin_idx, b.end_idx);
            if g.atoms[a].is_aromatic || g.atoms[c].is_aromatic {
                continue;
            }
            if tgroup[a] || tgroup[c] {
                continue;
            }
            // 環内二重結合は対象外 (InChI は小環の E/Z を持たない)
            if g.ring_atom_sets
                .iter()
                .any(|ring| ring.contains(&a) && ring.contains(&c))
            {
                continue;
            }
            if comp.num.get(a).copied().unwrap_or(0) == 0
                || comp.num.get(c).copied().unwrap_or(0) == 0
            {
                continue;
            }
            // 累積二重結合 (アレン等) の中心は対象外
            let has_other_multiple = |end: usize| {
                g.bonds.iter().enumerate().any(|(bj, bb)| {
                    bj != bi
                        && (bb.begin_idx == end || bb.end_idx == end)
                        && g.kekule_bond_orders[bj] >= 2.0
                })
            };
            if has_other_multiple(a) || has_other_multiple(c) {
                continue;
            }
            // 両端がそれぞれ立体源性 (2 つの置換基 — H・孤立電子対を含む —
            // が互いに異なる) であること。
            let end_stereogenic = |end: usize, other: usize| -> bool {
                let sym = g.atoms[end].symbol.as_str();
                let heavy = ez_neighbors(g, end, other);
                let n_h = g.adjacency[end]
                    .iter()
                    .filter(|&&x| g.atoms[x].symbol == "H")
                    .count();
                match (sym, heavy.len()) {
                    // C: 置換基 2 個が異なるとき (CIP ランクで判定)。
                    // =CH2 (H 2 個) は非立体源性。
                    ("C", 2) => ranks[heavy[0]] != ranks[heavy[1]],
                    ("C", 1) => n_h == 1,
                    // N: 孤立電子対が 2 つ目の置換基となるため、置換基
                    // (重原子または H) が 1 個なら常に立体源性。
                    ("N", 1) => n_h == 0,
                    ("N", 0) => n_h == 1,
                    _ => false,
                }
            };
            if !end_stereogenic(a, c) || !end_stereogenic(c, a) {
                continue;
            }
            let (ca, cc) = (canon_of(comp, a), canon_of(comp, c));
            let (hi, lo) = (ca.max(cc), ca.min(cc));
            if !listed.contains(&(hi, lo)) {
                entries.push((hi, lo, '?'));
            }
        }
    }

    entries.sort_unstable();
    entries
        .iter()
        .map(|(hi, lo, p)| format!("{hi}-{lo}{p}"))
        .collect::<Vec<_>>()
        .join(",")
}

/// 置換のパリティ (偶=0/奇=1)。src, dst は同一要素の並び。
fn perm_parity(src: &[usize], dst: &[usize]) -> usize {
    let mut work = src.to_vec();
    let mut sw = 0usize;
    for i in 0..dst.len() {
        if work[i] != dst[i] {
            let j = (i + 1..work.len()).find(|&j| work[j] == dst[i]).unwrap();
            work.swap(i, j);
            sw += 1;
        }
    }
    sw % 2
}

/// 四面体中心 1 個の生パリティ ('+'/'-')。四面体中心でなければ None。
/// raw = '-' iff (rs_bit XOR perm(CIP昇順→正準昇順)) == 1。
/// H は正準番号最大・CIP 最下位の仮想隣接として 4 番目に加える。
fn tetra_raw_parity(
    g: &MoleculeGraph,
    comp: &Component,
    ranks: &[usize],
    center: usize,
) -> Option<char> {
    let rs = g.atoms[center].chiral_tag?;
    let heavy: Vec<usize> = g.adjacency[center]
        .iter()
        .copied()
        .filter(|&nb| g.atoms[nb].symbol != "H")
        .collect();
    let n_h = g.adjacency[center]
        .iter()
        .filter(|&&nb| g.atoms[nb].symbol == "H")
        .count();
    // 4 heavy か 3 heavy + 1 H のみ対応 (孤立電子対中心は v2)
    let mut items: Vec<(i64, i64, usize)> = heavy
        .iter()
        .map(|&nb| (canon_of(comp, nb) as i64, ranks[nb] as i64, nb))
        .collect();
    if heavy.len() == 3 && n_h == 1 {
        items.push((i64::MAX, -1, usize::MAX)); // 仮想 H
    } else if heavy.len() != 4 {
        return None;
    }
    // CIP 昇順・正準昇順の原子 id 列
    let mut cip_asc = items.clone();
    cip_asc.sort_by_key(|&(_, r, id)| (r, id as i64));
    let mut canon_asc = items.clone();
    canon_asc.sort_by_key(|&(c, _, id)| (c, id as i64));
    let src: Vec<usize> = cip_asc.iter().map(|&(_, _, id)| id).collect();
    let dst: Vec<usize> = canon_asc.iter().map(|&(_, _, id)| id).collect();
    let pp = perm_parity(&src, &dst);
    let rs_bit = if rs == 'R' { 1 } else { 0 };
    Some(if (rs_bit ^ pp) == 1 { '-' } else { '+' })
}

/// `/t`・`/m`・`/s` 層。返り値は (t本体, mありなら Some(m文字), sありなら Some(s文字))。
/// 立体中心がなければ (空, None, None)。
pub(crate) fn tetrahedral_layers(
    g: &MoleculeGraph,
    comp: &Component,
) -> (String, Option<char>, Option<char>) {
    let ranks = cip_ranks(g);
    let mut centers: Vec<(usize, char)> = Vec::new();
    for &orig in &comp.inv {
        if let Some(raw) = tetra_raw_parity(g, comp, &ranks, orig) {
            centers.push((canon_of(comp, orig), raw));
        }
    }
    if centers.is_empty() {
        return (String::new(), None, None);
    }
    centers.sort_unstable();
    // /m 正規化: 最初の中心の生パリティが '+' なら全反転して '-' に、m=1。
    let invert = centers[0].1 == '+';
    let m = if invert { '1' } else { '0' };
    let t = centers
        .iter()
        .map(|&(c, p)| {
            let shown = if invert {
                if p == '+' {
                    '-'
                } else {
                    '+'
                }
            } else {
                p
            };
            format!("{c}{shown}")
        })
        .collect::<Vec<_>>()
        .join(",");
    (t, Some(m), Some('1'))
}

#[cfg(test)]
mod tests {
    use super::super::layers::build_components;
    use super::*;
    use crate::graph::build_molecule_graph;

    fn b_layer(smiles: &str) -> String {
        let g = build_molecule_graph(smiles).unwrap();
        let comps = build_components(&g);
        double_bond_layer(&g, &comps[0])
    }

    fn tms(smiles: &str) -> String {
        let g = build_molecule_graph(smiles).unwrap();
        let comps = build_components(&g);
        let (t, m, s) = tetrahedral_layers(&g, &comps[0]);
        let mut out = String::new();
        if !t.is_empty() {
            out.push_str(&format!("t{t}"));
        }
        if let Some(m) = m {
            out.push_str(&format!("/m{m}"));
        }
        if let Some(s) = s {
            out.push_str(&format!("/s{s}"));
        }
        out
    }

    #[test]
    fn tetrahedral() {
        assert_eq!(tms("C[C@H](N)C(=O)O"), "t2-/m0/s1");
        assert_eq!(tms("C[C@@H](N)C(=O)O"), "t2-/m1/s1");
        assert_eq!(tms("C[C@H](O)CC"), "t4-/m0/s1");
        assert_eq!(tms("C[C@@H](O)CC"), "t4-/m1/s1");
        assert_eq!(tms("F[C@H](Cl)Br"), "t1-/m0/s1");
        assert_eq!(tms("F[C@@H](Cl)Br"), "t1-/m1/s1");
        assert_eq!(tms("C[C@H](O)[C@@H](O)C"), "t3-,4-/m0/s1");
    }

    #[test]
    fn simple_ez() {
        assert_eq!(b_layer("C/C=C/C"), "4-3+"); // E
        assert_eq!(b_layer("C/C=C\\C"), "4-3-"); // Z
        assert_eq!(b_layer("F/C=C/Cl"), "2-1+");
        assert_eq!(b_layer("CC/C=C/CC"), "6-5+");
    }

    #[test]
    fn conjugated_diene() {
        assert_eq!(b_layer("C/C=C/C=C/C"), "5-3+,6-4+");
    }

    #[test]
    fn no_stereo() {
        assert_eq!(b_layer("CCCC"), "");
        assert_eq!(b_layer("C=C"), "");
    }

    #[test]
    fn undefined_stereogenic_cn_gets_question_mark() {
        // I18: 定義済みの立体二重結合がある分子では、立体源性だが未定義の
        // C=N (イミン/ヒドラゾン) を `?` として併記する。
        assert_eq!(b_layer("C/C=C/C(=N)OC"), "4-3+,6-5?");
        assert_eq!(b_layer("C/C=C/CC(=NN)C"), "4-3+,8-6?");
        // 未定義のみの分子では /b 層自体を省略する。
        assert_eq!(b_layer("CC(=NN)C"), "");
        // 可動 H 群のメンバー (アミジンの C=N) は対象外。
        assert_eq!(b_layer("C/C=C/C(=N)N"), "3-2+");
    }
}
