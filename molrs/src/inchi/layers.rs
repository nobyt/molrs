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
    /// 可動 H 群: (mobile H 数, 負電荷数, [canonical 番号...])
    pub mobile: Vec<(u8, u8, Vec<usize>)>,
}

/// 分子を成分ごとに正準情報へ分解する。
pub(crate) fn build_components(g: &MoleculeGraph) -> Vec<Component> {
    let tgroup = tautomer_group_members(g);
    let all_groups = super::number::mobile_groups(g); // (端点原子, 可動 H 数)
    let numbering = canonical_numbering(g);
    let n_h_of = |i: usize| {
        g.adjacency[i]
            .iter()
            .filter(|&&x| g.atoms[x].symbol == "H")
            .count() as u8
    };

    let mut comps: Vec<Component> = numbering
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
            // 可動群 (成分内のもの) を canonical 番号列に変換。可動 H 数 0 の
            // 群 (電荷分離ニトロ等、番号付けの等価化には使うが h 層には
            // 出さない) は除外する。
            let mut mobile: Vec<(u8, u8, Vec<usize>)> = all_groups
                .iter()
                .filter(|(eps, mh, _)| *mh > 0 && eps.iter().all(|&e| num[e] != 0))
                .map(|(eps, mh, neg)| {
                    let mut nums: Vec<usize> = eps.iter().map(|&e| num[e]).collect();
                    nums.sort_unstable();
                    (*mh, *neg, nums)
                })
                .collect();
            // 出力順: 群のメンバー数 (端点原子数) 昇順、同数なら最小 canonical
            // 番号昇順 (I19 §3.5、実 InChI で確認: 小さい群が先に来る。
            // 例: アルギニンは (H,11,12) [2 端点] が (H4,8,9,10) [3 端点] より
            // 先)。
            mobile.sort_by_key(|(_, _, nums)| (nums.len(), nums[0]));

            Component {
                inv: inv.clone(),
                num,
                adj,
                fixed_h,
                mobile,
            }
        })
        .collect();
    // 成分順序は式層と共通の規則 (I20)
    comps.sort_by_key(|c| super::formula::component_sort_key(g, &c.inv));
    comps
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
    // 部分木の「描画重み」= 部分木の原子数 + 部分木内で言及される環閉合
    // (後退辺) の数。子の出力順に使う: 軽い枝 (末端 O 等) が先に括弧、
    // 重い幹が最後にインライン。環閉合も 1 項目として数えるのが実 InChI の
    // 挙動 (例: 末端 O の枝 `(7)` は環閉合を含む枝 `6-3` より先)。
    let mut subtree = vec![0usize; n + 1];
    for c in 1..=n {
        if discovery[c] == usize::MAX {
            continue;
        }
        // 自身 1 + この原子で言及する後退辺の数
        let backs = comp.adj[c]
            .iter()
            .filter(|&&x| {
                x != parent[c] && !children[c].contains(&x) && discovery[x] < discovery[c]
            })
            .count();
        subtree[c] = 1 + backs;
    }
    // discovery 降順 = 葉から根へ (後順)
    let mut order: Vec<usize> = (1..=n).filter(|&c| discovery[c] != usize::MAX).collect();
    order.sort_by_key(|&c| std::cmp::Reverse(discovery[c]));
    for &c in &order {
        for &ch in &children[c] {
            subtree[c] += subtree[ch];
        }
    }

    // (2) 直列化
    let mut out = start.to_string();
    serialize_node(
        comp, start, &discovery, &parent, &children, &subtree, &mut out,
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

    // 全項目 = 後退辺 (葉、昇順) ++ 木の子 (部分木/番号昇順)。
    // 最後の 1 個をインライン、それ以外は 1 つのカンマ括弧に入れる
    // (各項目は自身の部分木も直列化する。例: 四級 N の (6-2,7-3)8-4)。
    #[derive(Clone, Copy)]
    enum Item {
        Back(usize),
        Child(usize),
    }
    let mut items: Vec<Item> = Vec::new();
    for &b in &back {
        items.push(Item::Back(b));
    }
    for &c in &kids {
        items.push(Item::Child(c));
    }
    let Some(inline) = items.pop() else {
        return;
    };

    let render = |out: &mut String, it: Item| match it {
        Item::Back(b) => out.push_str(&b.to_string()),
        Item::Child(c) => {
            out.push_str(&c.to_string());
            serialize_node(comp, c, discovery, parent, children, subtree, out);
        }
    };

    let mut prev_was_paren = false;
    if !items.is_empty() {
        out.push('(');
        for (k, &it) in items.iter().enumerate() {
            if k > 0 {
                out.push(',');
            }
            render(out, it);
        }
        out.push(')');
        prev_was_paren = true;
    }
    if !prev_was_paren {
        out.push('-');
    }
    render(out, inline);
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
    // 固定 H 部はカンマ結合
    let fixed = parts.join(",");
    // 可動群はカンマなしで連結 ((H,5,6)(H,7,8))。固定部との間にはカンマ。
    let mut mobile = String::new();
    for (mh, neg, nums) in &comp.mobile {
        let list = nums
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",");
        // `H` + (2 以上なら個数) + (負電荷があれば `-`)。実 InChI の
        // `(H-,10,11,12)` は「可動 H 1 個 + 負電荷 1」、`(H3-,…)` は
        // 「可動 H 3 個 + 負電荷 1」(I32)。
        let count = if *mh == 1 {
            String::new()
        } else {
            mh.to_string()
        };
        let charge = "-".repeat(*neg as usize);
        mobile.push_str(&format!("(H{count}{charge},{list})"));
    }
    match (fixed.is_empty(), mobile.is_empty()) {
        (false, false) => format!("{fixed},{mobile}"),
        (false, true) => fixed,
        (true, false) => mobile,
        (true, true) => String::new(),
    }
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

/// 同位体層 `/i` の本体 (先頭の `i` は含めない)。空なら空文字列 (I37)。
///
/// 2 種類のエントリを正準番号の昇順で `,` 区切りに並べる:
///
/// - **重原子そのものが同位体**: `{番号}{+|-}{標準質量数との差}` (`4+1` = 13C)
/// - **同位体水素が付いている**: `{番号}T{個数}D{個数}` (`5T2`、`12D2`、`1TD`)。
///   個数 1 は省略。T (三重水素) を D より先に書く。
///
/// 同位体情報は `AtomInfo` に無く、SMILES パース結果 (`g.parsed`) にしか
/// 無いので `parser_to_graph` を逆引きして辿る。`[2H]`/`[3H]` は同位体付き
/// なので `build_molecule_graph` の H マージ対象外で、原子ノードとして残る
/// (このため電荷正規化も走らず、添字はそのまま使える)。
pub(crate) fn isotope_layer(g: &MoleculeGraph, comp: &Component) -> String {
    // グラフ原子 idx → パーサ原子 idx
    let mut to_parsed = vec![usize::MAX; g.atoms.len()];
    for (pi, gi) in g.parser_to_graph.iter().enumerate() {
        if let Some(gi) = gi {
            if *gi < to_parsed.len() {
                to_parsed[*gi] = pi;
            }
        }
    }
    let isotope_of = |gi: usize| -> Option<u16> {
        let pi = *to_parsed.get(gi)?;
        g.parsed.atoms.get(pi)?.isotope
    };

    let mut entries: Vec<(usize, String)> = Vec::new();
    for (ci, &orig) in comp.inv.iter().enumerate() {
        let canon = ci + 1;
        // 1) 重原子自身の同位体 (標準質量数との差)
        if let (Some(iso), Some(std)) = (isotope_of(orig), standard_mass(&g.atoms[orig].symbol)) {
            let delta = iso as i32 - std as i32;
            if delta != 0 {
                entries.push((canon, format!("{canon}{delta:+}")));
            }
        }
        // 2) 付いている同位体水素 (D = 2H, T = 3H)
        let (mut d, mut t) = (0usize, 0usize);
        for &nb in &g.adjacency[orig] {
            if g.atoms[nb].symbol != "H" {
                continue;
            }
            match isotope_of(nb) {
                Some(2) => d += 1,
                Some(3) => t += 1,
                _ => {}
            }
        }
        if d > 0 || t > 0 {
            let mut spec = format!("{canon}");
            let mut push = |sym: char, n: usize| {
                if n > 0 {
                    spec.push(sym);
                    if n > 1 {
                        spec.push_str(&n.to_string());
                    }
                }
            };
            push('T', t);
            push('D', d);
            entries.push((canon, spec));
        }
    }
    entries.sort_by_key(|&(c, _)| c);
    entries
        .into_iter()
        .map(|(_, s)| s)
        .collect::<Vec<_>>()
        .join(",")
}

/// InChI が同位体差の基準に使う質量数 (最も存在比の高い同位体)。
fn standard_mass(sym: &str) -> Option<u16> {
    Some(match sym {
        "H" => 1,
        "B" => 11,
        "C" => 12,
        "N" => 14,
        "O" => 16,
        "F" => 19,
        "Si" => 28,
        "P" => 31,
        "S" => 32,
        "Cl" => 35,
        "As" => 75,
        "Se" => 80,
        "Br" => 79,
        "I" => 127,
        _ => return None,
    })
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

    #[test]
    fn hydrogen_mobile_groups_sorted_by_size_ascending() {
        // I19 §3.5: 複数の可動 H 群があるとき、実 InChI はメンバー数
        // (端点原子数) の少ない群を先に出す (最小 canonical 番号の昇順では
        // ない)。アルギニン: グアニジノ基の 3 端点・可動 H4 の群より、
        // カルボキシルの 2 端点・可動 H1 の群が先。
        let g = build_molecule_graph("NC(CCCNC(N)=N)C(=O)O").unwrap();
        let h = crate::inchi::to_inchi(&g).unwrap();
        assert!(h.ends_with(",(H,11,12)(H4,8,9,10)"), "got: {h}");
    }
}
