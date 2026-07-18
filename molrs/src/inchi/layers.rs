//! InChI の接続層 `c`・水素層 `h` の直列化 (I4)。
//!
//! 正準番号 (number.rs) を前提に、公式 InChI と同形式の文字列を組み立てる。
//! v1 は単一成分・中性〜単純電荷を対象 (q/p は charge.rs)。

use crate::graph::MoleculeGraph;

use super::number::{canonical_numbering, tautomer_group_members};

/// 1 成分の正準情報 (番号付け結果から構築)。
pub(crate) struct Component {
    /// canonical 番号 (1..=n) → 元の原子 idx
    pub inv: Vec<usize>,
    /// 元の原子 idx → canonical 番号 (1..=n)。q/p 層 (I4 続き) で使用予定。
    #[allow(dead_code)]
    pub num: Vec<usize>,
    /// canonical 隣接 (num a → [num b, ...]、a のみ)。1-indexed。
    pub adj: Vec<Vec<usize>>,
    /// canonical 番号 → 固定 H 数
    pub fixed_h: Vec<u8>,
    /// 可動 H 群: (mobile H 数, [canonical 番号...])
    pub mobile: Vec<(u8, Vec<usize>)>,
}

/// 分子を成分ごとに正準情報へ分解する。
pub(crate) fn build_components(g: &MoleculeGraph) -> Vec<Component> {
    let tgroup = tautomer_group_members(g);
    let numbering = canonical_numbering(g);
    let heavy_deg = |i: usize| {
        g.adjacency[i]
            .iter()
            .filter(|&&x| g.atoms[x].symbol != "H")
            .count()
    };
    let n_h_of = |i: usize| {
        g.adjacency[i]
            .iter()
            .filter(|&&x| g.atoms[x].symbol == "H")
            .count() as u8
    };

    numbering
        .iter()
        .map(|inv| {
            let n = inv.len();
            // num: 元 idx → canonical (1-indexed)。inv[c-1] = 元 idx
            let max_idx = g.atoms.len();
            let mut num = vec![0usize; max_idx];
            for (ci, &orig) in inv.iter().enumerate() {
                num[orig] = ci + 1;
            }
            // canonical 隣接 (成分内のみ)
            let mut adj = vec![Vec::new(); n + 1];
            for (ci, &orig) in inv.iter().enumerate() {
                let cnum = ci + 1;
                for &nb in &g.adjacency[orig] {
                    if g.atoms[nb].symbol != "H" && num[nb] != 0 && num[nb] != cnum {
                        adj[cnum].push(num[nb]);
                    }
                }
                adj[cnum].sort_unstable();
                adj[cnum].dedup();
            }
            // 固定 H (t-group メンバーは 0) と可動群
            let mut fixed_h = vec![0u8; n + 1];
            for (ci, &orig) in inv.iter().enumerate() {
                if !tgroup[orig] {
                    fixed_h[ci + 1] = n_h_of(orig);
                }
            }
            // 可動群: 同一中心に結合した t-group 端点をまとめる
            let mobile = collect_mobile_groups(g, inv, &num, &tgroup, &heavy_deg, &n_h_of);

            Component {
                inv: inv.clone(),
                num,
                adj,
                fixed_h,
                mobile,
            }
        })
        .collect()
}

/// 可動 H 群を成分内で収集する。各群 = ある中心の t-group 端点集合。
/// mobile H 数 = 端点上の H 総数 + 中和される負電荷数。
fn collect_mobile_groups(
    g: &MoleculeGraph,
    inv: &[usize],
    num: &[usize],
    tgroup: &[bool],
    heavy_deg: &impl Fn(usize) -> usize,
    n_h_of: &impl Fn(usize) -> u8,
) -> Vec<(u8, Vec<usize>)> {
    let comp_set: std::collections::HashSet<usize> = inv.iter().copied().collect();
    let mut groups: Vec<(u8, Vec<usize>)> = Vec::new();
    let mut used = std::collections::HashSet::new();
    // 中心をなめ、その t-group 端点をまとめる
    for &center in inv {
        let endpoints: Vec<usize> = g.adjacency[center]
            .iter()
            .copied()
            .filter(|&nb| {
                tgroup[nb]
                    && comp_set.contains(&nb)
                    && heavy_deg(nb) == 1
                    && g.atoms[nb].symbol != "H"
            })
            .collect();
        if endpoints.len() < 2 || endpoints.iter().any(|e| used.contains(e)) {
            continue;
        }
        let n_hydro: u8 = endpoints.iter().map(|&e| n_h_of(e)).sum();
        let n_neg = endpoints
            .iter()
            .filter(|&&e| g.atoms[e].formal_charge < 0)
            .count() as u8;
        let mobile_h = n_hydro + n_neg;
        let mut nums: Vec<usize> = endpoints.iter().map(|&e| num[e]).collect();
        nums.sort_unstable();
        for &e in &endpoints {
            used.insert(e);
        }
        groups.push((mobile_h, nums));
    }
    // 出力順: 群の最小 canonical 番号昇順
    groups.sort_by_key(|(_, nums)| nums[0]);
    groups
}

/// c 層の本体 (先頭の `c` は含めない)。単一成分。
///
/// 2 段構成: (1) atom 1 から「最小未訪問隣接優先」の DFS で全域木と後退辺を
/// 決める。(2) 各原子で「後退辺 (昇順) → 木の子 (昇順)」を列挙し、最後の
/// 項目はインライン、他は括弧。連結子 `-` は直前が `)` の項目には付けない。
pub(crate) fn connection_layer(comp: &Component) -> String {
    let n = comp.inv.len();
    if n <= 1 {
        return String::new();
    }
    // 開始原子 = 最小番号の末端 (次数 1) 原子。末端がなければ (純環系) atom 1。
    let start = (1..=n).find(|&c| comp.adj[c].len() == 1).unwrap_or(1);

    // (1) DFS で全域木・discovery 時刻を決める
    let mut discovery = vec![usize::MAX; n + 1];
    let mut parent = vec![0usize; n + 1];
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); n + 1];
    let mut clock = 0usize;
    // 明示スタックで再帰 DFS を模倣 (深い分子でのスタック溢れ回避)
    let mut stack: Vec<usize> = vec![start];
    discovery[start] = clock;
    clock += 1;
    parent[start] = 0;
    while let Some(&a) = stack.last() {
        // a の最小未訪問隣接を 1 つ取る
        let next = comp.adj[a]
            .iter()
            .copied()
            .filter(|&x| discovery[x] == usize::MAX)
            .min();
        match next {
            Some(c) => {
                discovery[c] = clock;
                clock += 1;
                parent[c] = a;
                children[a].push(c);
                stack.push(c);
            }
            None => {
                stack.pop();
            }
        }
    }
    // 部分木サイズ (子の出力順に使う): 葉が先、幹 (最大部分木) が最後
    let mut subtree = vec![1usize; n + 1];
    // discovery 降順 = 葉から根へ (後順)
    let mut order: Vec<usize> = (1..=n).filter(|&c| discovery[c] != usize::MAX).collect();
    order.sort_by_key(|&c| std::cmp::Reverse(discovery[c]));
    for &c in &order {
        for &ch in &children[c] {
            subtree[c] += subtree[ch];
        }
    }

    // 各ノードが後退辺を持つか (出力上の「純末端」判定に使う)。
    // 純末端 = 木の子なし かつ 後退辺なし → カンマ結合してよい。
    let mut has_back = vec![false; n + 1];
    for c in 1..=n {
        if discovery[c] == usize::MAX {
            continue;
        }
        has_back[c] = comp.adj[c]
            .iter()
            .any(|&x| x != parent[c] && !children[c].contains(&x) && discovery[x] < discovery[c]);
    }
    let terminal: Vec<bool> = (0..=n)
        .map(|c| c >= 1 && c <= n && children[c].is_empty() && !has_back[c])
        .collect();

    // (2) 直列化
    let mut out = start.to_string();
    serialize_node(
        comp, start, &discovery, &parent, &children, &subtree, &terminal, &mut out,
    );
    out
}

#[allow(clippy::too_many_arguments)]
fn serialize_node(
    comp: &Component,
    a: usize,
    discovery: &[usize],
    parent: &[usize],
    children: &[Vec<usize>],
    subtree: &[usize],
    terminal: &[bool],
    out: &mut String,
) {
    // 後退辺: 親でも木の子でもない隣接で、より早く発見された相手 (a 側で 1 度)
    let mut back: Vec<usize> = comp.adj[a]
        .iter()
        .copied()
        .filter(|&x| x != parent[a] && !children[a].contains(&x) && discovery[x] < discovery[a])
        .collect();
    back.sort_unstable();
    // 子は (部分木サイズ, 番号) 昇順 → 幹 (最大部分木) が最後 = インライン
    let mut kids = children[a].clone();
    kids.sort_by_key(|&c| (subtree[c], c));

    // インラインは最後の 1 個 (子があれば最大部分木の子、なければ最後の後退辺)
    let inline_child = kids.pop(); // None なら子なし
    let mut prev_was_paren = false;

    // 後退辺: 葉 (子なし) なのでまとめて 1 つのカンマ括弧にする
    // (子がインラインになる場合)。後退辺のみで子がないときは最後をインライン。
    let inline_back = if inline_child.is_none() {
        back.pop()
    } else {
        None
    };
    if !back.is_empty() {
        out.push('(');
        out.push_str(
            &back
                .iter()
                .map(|b| b.to_string())
                .collect::<Vec<_>>()
                .join(","),
        );
        out.push(')');
        prev_was_paren = true;
    }

    // 残りの子 (インライン以外): 純末端 (木の子も後退辺もなし) は 1 つの
    // カンマ括弧に統合、それ以外は各自の括弧で再帰。純末端を先に。
    let leaves: Vec<usize> = kids.iter().copied().filter(|&c| terminal[c]).collect();
    let branches: Vec<usize> = kids.iter().copied().filter(|&c| !terminal[c]).collect();
    if !leaves.is_empty() {
        out.push('(');
        out.push_str(
            &leaves
                .iter()
                .map(|c| c.to_string())
                .collect::<Vec<_>>()
                .join(","),
        );
        out.push(')');
        prev_was_paren = true;
    }
    for &c in &branches {
        out.push('(');
        out.push_str(&c.to_string());
        serialize_node(comp, c, discovery, parent, children, subtree, terminal, out);
        out.push(')');
        prev_was_paren = true;
    }

    // インライン (最後)
    if let Some(c) = inline_child {
        if !prev_was_paren {
            out.push('-');
        }
        out.push_str(&c.to_string());
        serialize_node(comp, c, discovery, parent, children, subtree, terminal, out);
    } else if let Some(b) = inline_back {
        if !prev_was_paren {
            out.push('-');
        }
        out.push_str(&b.to_string());
    }
}

/// h 層の本体 (先頭の `h` は含めない)。単一成分。
pub(crate) fn hydrogen_layer(comp: &Component) -> String {
    let n = comp.inv.len();
    // 固定 H: 数ごとにグループ化 (数昇順)、原子番号は範囲圧縮
    let mut by_count: std::collections::BTreeMap<u8, Vec<usize>> =
        std::collections::BTreeMap::new();
    for c in 1..=n {
        let h = comp.fixed_h[c];
        if h > 0 {
            by_count.entry(h).or_default().push(c);
        }
    }
    let mut parts: Vec<String> = Vec::new();
    for (count, atoms) in &by_count {
        let list = compress_ranges(atoms);
        if *count == 1 {
            parts.push(format!("{list}H"));
        } else {
            parts.push(format!("{list}H{count}"));
        }
    }
    // 可動群
    for (mh, nums) in &comp.mobile {
        let list = nums
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",");
        if *mh == 1 {
            parts.push(format!("(H,{list})"));
        } else {
            parts.push(format!("(H{mh},{list})"));
        }
    }
    parts.join(",")
}

/// 昇順の番号列を InChI の範囲表記に圧縮する (`1,2,3,5` → `1-3,5`)。
fn compress_ranges(nums: &[usize]) -> String {
    if nums.is_empty() {
        return String::new();
    }
    let mut out = Vec::new();
    let mut start = nums[0];
    let mut prev = nums[0];
    for &x in &nums[1..] {
        if x == prev + 1 {
            prev = x;
        } else {
            out.push(if start == prev {
                start.to_string()
            } else {
                format!("{start}-{prev}")
            });
            start = x;
            prev = x;
        }
    }
    out.push(if start == prev {
        start.to_string()
    } else {
        format!("{start}-{prev}")
    });
    out.join(",")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    fn layers(smiles: &str) -> (String, String) {
        let g = build_molecule_graph(smiles).unwrap();
        let comps = build_components(&g);
        assert_eq!(comps.len(), 1, "multi-component in test");
        (connection_layer(&comps[0]), hydrogen_layer(&comps[0]))
    }

    #[test]
    fn connection_simple() {
        assert_eq!(layers("CCO").0, "1-2-3");
        assert_eq!(layers("CC(=O)O").0, "1-2(3)4");
        assert_eq!(layers("CCOCC").0, "1-3-5-4-2");
        assert_eq!(layers("CC(C)CC").0, "1-4-5(2)3");
    }

    #[test]
    fn connection_rings() {
        assert_eq!(layers("c1ccccc1").0, "1-2-4-6-5-3-1");
        assert_eq!(layers("C1CCCCC1").0, "1-2-4-6-5-3-1");
        assert_eq!(layers("c1ccc2ccccc2c1").0, "1-2-6-10-8-4-3-7-9(10)5-1");
    }

    #[test]
    fn hydrogen_simple() {
        assert_eq!(layers("CCO").1, "3H,2H2,1H3");
        assert_eq!(layers("c1ccccc1").1, "1-6H");
        assert_eq!(layers("C1CCCCC1").1, "1-6H2");
        assert_eq!(layers("CC(C)CC").1, "5H,4H2,1-3H3");
    }

    #[test]
    fn hydrogen_mobile() {
        assert_eq!(layers("CC(=O)O").1, "1H3,(H,3,4)");
        assert_eq!(layers("CC(=O)N").1, "1H3,(H2,3,4)");
        assert_eq!(layers("OCC(=O)O").1, "3H,1H2,(H,4,5)");
    }
}
