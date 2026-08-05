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
    // 重原子 4 個の中心 (第四級) は 1 個反転する (I30)。
    //
    // R/S 側 (`crate::stereo::assign_atom_chiral_codes`) は隣接 3 個 + 暗黙 H の
    // ときだけ `n_swaps += 1` の補正を入れる。こちらの仮想 H (CIP 最下位 =
    // 並びの先頭、正準番号最大 = 並びの末尾) を 4 要素の置換に含めると、その
    // 補正とちょうど打ち消し合って 3 重原子 + H の場合は正しくなるが、
    // 重原子 4 個の場合は打ち消す相手がないぶん符号が 1 つずれる。
    //
    // 最小再現: `C[C@](N)(O)CC` は実 InChI が `/t4-/m0/s1`、補正なしだと
    // `/m1` になる。`C[C@H](N)C(=O)O` (3 重原子 + H) は補正なしで正しい。
    let quaternary = usize::from(heavy.len() == 4);
    Some(if (rs_bit ^ pp ^ quaternary) == 1 {
        '-'
    } else {
        '+'
    })
}

/// **立体源性だが SMILES で配置が指定されていない**四面体中心か (I30)。
///
/// 実 InChI はこの中心を `/t` に `12?` として列挙する。molrs は
/// [`tetra_raw_parity`] が `chiral_tag` 無しで `None` を返すため、従来は
/// 完全に脱落していた (PubChem 実データで 239 件の不一致の原因)。
///
/// 置換基が互いに異なるかは **1-WL 精緻化色** ([`refined_colors`]) で見る。
/// CIP ランクは同順位をほとんど作らない全順序に近く、これで判定すると sp3
/// 炭素が軒並み立体源性になってしまう (実測 94.66% → 61.76% に悪化した)。
/// 1-WL 色は対称な枝 (イソプロピルの 2 つのメチル等) に同じ色を与えるので、
/// 「4 つの置換基が構成的に区別できるか」を正しく表す。
fn is_undefined_tetra(
    g: &MoleculeGraph,
    comp: &Component,
    colors: &[usize],
    center: usize,
) -> bool {
    let a = &g.atoms[center];
    if a.chiral_tag.is_some() || a.is_aromatic {
        return false;
    }
    // 四面体になりうる中心元素。炭素以外は実 InChI が扱う範囲に絞る
    // (第 14 族と、孤立電子対を置換基に数えない荷電 N/P/S)。
    let ok_elem = match a.symbol.as_str() {
        "C" | "Si" | "Ge" | "Sn" => true,
        "N" | "P" | "S" => a.formal_charge > 0,
        _ => false,
    };
    if !ok_elem {
        return false;
    }
    let heavy: Vec<usize> = g.adjacency[center]
        .iter()
        .copied()
        .filter(|&nb| g.atoms[nb].symbol != "H")
        .collect();
    let n_h = g.adjacency[center].len() - heavy.len();
    // 4 heavy か 3 heavy + 1 H のみ (二重結合を持つ中心は接続数が 4 未満)
    if !(heavy.len() == 4 && n_h == 0 || heavy.len() == 3 && n_h == 1) {
        return false;
    }
    // 重原子置換基の 1-WL 色が互いに相異なること (H は 1 個までなので
    // 重原子とは必ず区別できる)
    let mut cs: Vec<usize> = heavy.iter().map(|&nb| colors[canon_of(comp, nb)]).collect();
    cs.sort_unstable();
    let before = cs.len();
    cs.dedup();
    cs.len() == before
}

/// 正準番号空間の色 (自己同型が保たなければならない原子の不変量) を
/// 1-WL 精緻化で求める。返り値は 1-indexed (添字 0 は未使用)。
fn refined_colors(g: &MoleculeGraph, comp: &Component) -> Vec<usize> {
    let n = comp.inv.len();
    // 初期色: 元素・電荷・固定 H 数・次数 (正準番号付けが使う不変量と同じ粒度)
    let mut keys: Vec<String> = (0..=n)
        .map(|c| {
            if c == 0 {
                return String::new();
            }
            let a = &g.atoms[comp.inv[c - 1]];
            format!(
                "{}|{}|{}|{}",
                a.symbol,
                a.formal_charge,
                comp.fixed_h[c],
                comp.adj[c].len()
            )
        })
        .collect();
    let mut color = assign_ids(&keys);
    for _ in 0..n {
        let next_keys: Vec<String> = (0..=n)
            .map(|c| {
                if c == 0 {
                    return String::new();
                }
                let mut nb: Vec<usize> = comp.adj[c].iter().map(|&x| color[x]).collect();
                nb.sort_unstable();
                format!("{}|{nb:?}", color[c])
            })
            .collect();
        let next = assign_ids(&next_keys);
        if next == color {
            break;
        }
        color = next;
        keys = next_keys;
    }
    let _ = keys;
    color
}

/// 文字列キー列 → 辞書順に基づく連番 ID 列。
fn assign_ids(keys: &[String]) -> Vec<usize> {
    let mut uniq: Vec<&String> = keys.iter().collect();
    uniq.sort_unstable();
    uniq.dedup();
    keys.iter()
        .map(|k| uniq.binary_search(&k).unwrap())
        .collect()
}

/// 分子がアキラル (メソ体など、鏡像が自分自身と一致する) かどうか。
///
/// 実 InChI は構造とその鏡像を両方正準化し、結果が同一なら `/m`・`/s` を
/// 出力しない (メソ酒石酸 `.../t1-,2+` に `/m` が付かないのはこのため)。
/// 鏡像は全ての生パリティを反転したものなので、条件は「グラフの自己同型 π が
/// 存在して、π で番号を付け替えた割り当てが全反転した割り当てと一致する」。
///
/// π で付け替えると、中心 c のパリティは隣接の正準昇順が π によって並べ替わる
/// 分だけ反転する (`new(π(c)) = old(c) XOR σ_c`、σ_c はその置換のパリティ)。
/// CIP ランクと R/S は自己同型で不変なのでこの項だけが効く。したがって
/// 求める条件は全ての中心で `old(c) XOR σ_c == !old(π(c))`。
///
/// 探索は色クラス内のバックトラッキング。ノード数に上限を設け、超えた場合は
/// 「アキラルと証明できなかった」= false を返す (安全側 — `/m` は出力される)。
fn is_achiral(g: &MoleculeGraph, comp: &Component, centers: &[(usize, char)]) -> bool {
    let n = comp.inv.len();
    let color = refined_colors(g, comp);
    // 中心の正準番号 → パリティ (0/1)
    let mut parity: Vec<Option<usize>> = vec![None; n + 1];
    for &(c, p) in centers {
        parity[c] = Some(usize::from(p == '-'));
    }
    // 各中心の置換基正準番号 (H は仮想的に最大値)
    let subs = |c: usize| -> Vec<usize> {
        let mut v = comp.adj[c].clone();
        let orig = comp.inv[c - 1];
        let n_h = g.adjacency[orig]
            .iter()
            .filter(|&&x| g.atoms[x].symbol == "H")
            .count();
        if v.len() == 3 && n_h == 1 {
            v.push(usize::MAX);
        }
        v
    };

    let mut pi = vec![0usize; n + 1];
    let mut used = vec![false; n + 1];
    let mut budget = 200_000usize;
    // 局所的な再帰ヘルパ。状態をまとめた構造体を切るより引数のままの方が
    // 読みやすいので、引数個数の lint はここでは無効化する。
    #[allow(clippy::too_many_arguments)]
    fn search(
        c: usize,
        n: usize,
        color: &[usize],
        adj: &[Vec<usize>],
        pi: &mut Vec<usize>,
        used: &mut Vec<bool>,
        budget: &mut usize,
        check: &dyn Fn(&[usize]) -> bool,
    ) -> bool {
        if c > n {
            return check(pi);
        }
        for cand in 1..=n {
            if *budget == 0 {
                return false;
            }
            *budget -= 1;
            if used[cand] || color[cand] != color[c] {
                continue;
            }
            // 既に割り当て済みの原子との隣接関係が保たれるか
            if !(1..c).all(|k| adj[c].contains(&k) == adj[cand].contains(&pi[k])) {
                continue;
            }
            pi[c] = cand;
            used[cand] = true;
            if search(c + 1, n, color, adj, pi, used, budget, check) {
                return true;
            }
            used[cand] = false;
        }
        false
    }

    let check = |pi: &[usize]| -> bool {
        for &(c, _) in centers {
            let img = pi[c];
            let (Some(pc), Some(pimg)) = (parity[c], parity[img]) else {
                return false;
            };
            let s = subs(c);
            let mut a = s.clone();
            a.sort_unstable();
            let mut b = s.clone();
            b.sort_by_key(|&x| if x == usize::MAX { usize::MAX } else { pi[x] });
            let sigma = perm_parity(&a, &b);
            if (pc ^ sigma) != 1 - pimg {
                return false;
            }
        }
        true
    };

    search(
        1,
        n,
        &color,
        &comp.adj,
        &mut pi,
        &mut used,
        &mut budget,
        &check,
    )
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
    // 未定義の立体源性中心 (`?`): **定義済みの中心が 1 つでもある場合のみ**
    // 併記する (I30)。全て未定義の分子では `/t` 層自体が省略される —
    // `CC(CN)O` の実 InChI に `/t3?` が付かないのがこれで、`/b` 側の未定義
    // 二重結合とまったく同じ規則。
    let mut defined = centers.clone();
    defined.sort_unstable();
    let colors = refined_colors(g, comp);
    for &orig in &comp.inv {
        if is_undefined_tetra(g, comp, &colors, orig) {
            centers.push((canon_of(comp, orig), '?'));
        }
    }
    centers.sort_unstable();
    // /m 正規化: 最初の**定義済み**中心の生パリティが '+' なら全反転して
    // '-' に、m=1。未定義中心 (`?`) は反転の対象にも基準にもならない。
    let invert = defined[0].1 == '+';
    // メソ体 (鏡像が自分自身と一致) では実 InChI は /m・/s を出さない。
    let m = if is_achiral(g, comp, &defined) {
        None
    } else {
        Some(if invert { '1' } else { '0' })
    };
    let t = centers
        .iter()
        .map(|&(c, p)| {
            let shown = if invert && p != '?' {
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
    let s = m.map(|_| '1');
    (t, m, s)
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

    /// メソ体 (鏡像が自分自身と一致) は `/t` だけで `/m`・`/s` を出さない。
    /// 対称に等価な 2 中心が逆パリティを持ち、両者を入れ替える自己同型が
    /// 全パリティ反転を実現する場合。
    #[test]
    fn meso_compounds_omit_m_and_s() {
        // メソ酒石酸
        assert_eq!(tms("OC(=O)[C@@H](O)[C@@H](O)C(=O)O"), "t1-,2+");
        // cis-シクロヘキサン-1,2-ジオール
        assert_eq!(tms("O[C@@H]1CCCC[C@@H]1O"), "t5-,6+");
        assert_eq!(tms("[C@@H](F)([C@H](F)Cl)Cl"), "t1-,2+");
    }

    /// 逆に、光学活性な (対称でない、あるいは同符号の) 分子では /m・/s が
    /// 残ること — メソ判定が広すぎないことの回帰テスト。
    #[test]
    fn chiral_compounds_keep_m_and_s() {
        // (R,R)/(S,S)-酒石酸は光学活性
        assert_eq!(tms("OC(=O)[C@@H](O)[C@H](O)C(=O)O"), "t1-,2-/m0/s1");
        assert_eq!(tms("O[C@@H]1CCCC[C@H]1O"), "t5-,6-/m1/s1");
        assert_eq!(tms("C[C@H](N)C(=O)O"), "t2-/m0/s1");
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
