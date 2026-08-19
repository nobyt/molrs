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
            // 出さない) は除外する。負電荷だけが残った群 (脱プロトン後の
            // カルボキシラート) も同じ — `(H0-,…)` という表記は存在せず、
            // カルニチンの実 InChI は h 層に何も出さない (I38)。
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
    // 成分順序は式層と共通の規則 (I20)。式・原子数・H 数まで一致する構成
    // 異性体どうし (同じ分子式の異なる化合物の混合物) は、それだけでは
    // タイブレークできず安定ソートで入力順のまま残ってしまう。
    //
    // I44 では `/c` の**レンダリング後の文字列**を辞書式に比較していたが、
    // 実 InChI (`ichimake.c` の `CompINChI2`) は接続表を**レンダリング前の
    // 整数配列**として比較する (`nConnTable[i]`、要素ごとに「大きい方が先」、
    // 配列長も「長い方が先」) — 文字列比較は 2 桁の正準番号が混ざる配列で
    // 数値比較と食い違う (例: 文字列では `"12" < "9"`、数値では `12 > 9`)
    // ため、3 成分以上が同時にタイする場合に実 InChI と順序が合わないことが
    // あった (I46 の節を参照)。`connection_layer` が使う `Component::adj`
    // (正準番号ごとの整数隣接リスト) をそのまま数値比較することで、
    // レンダリングを経由しない忠実な再現にした (I47)。
    let total_h = |inv: &[usize]| -> usize { inv.iter().map(|&i| n_h_of(i) as usize).sum() };
    comps.sort_by(|a, b| {
        super::formula::component_sort_key(g, &a.inv)
            .cmp(&super::formula::component_sort_key(g, &b.inv))
            .then_with(|| compare_conn_table_numeric(b, a))
            // `CompINChI2` は conn table が一致した後、合計 H 数を「多い方が
            // 先」で比較する (1792〜1798 行、`num_H2 - num_H1`)。骨格
            // (conn table) が同じ飽和/不飽和対をここで正しく隣接させる。
            .then_with(|| total_h(&b.inv).cmp(&total_h(&a.inv)))
            // 合計 H 数まで一致した場合、`CompINChI2` はさらに原子ごとの H 数
            // (`nNum_H[i]`、可動群メンバーは 0 扱い = `fixed_h` と同じ) を
            // 正準番号順に比較する (1800〜1812 行): 「H 0 の原子が先、非 0
            // 同士なら H が多い方が先」(コメント曰く `N < NH3 < NH2 < NH` の
            // 順、すなわち出力順は N→NH3→NH2→NH)。全く同じ接続表・合計 H を
            // 持つ構造異性体成分どうし (例: インドール互変異性体の混合物)
            // が、この段まで来て初めて順序が決まるケースがある (I80)。
            .then_with(|| per_atom_h_key(a).cmp(&per_atom_h_key(b)))
    });
    comps
}

/// [`build_components`] の並べ替え追加タイブレークが使う、`nNum_H[i]` 相当の
/// 原子ごと優先度キー。「H 0 が最優先、非 0 同士なら H が多い方が優先」を
/// 通常の昇順 `Ord` にそのまま乗るよう符号を反転して表現する
/// (H=0 → `i64::MIN`、H=h(>0) → `-h`)。
fn per_atom_h_key(c: &Component) -> Vec<i64> {
    c.fixed_h[1..]
        .iter()
        .map(|&h| if h == 0 { i64::MIN } else { -(h as i64) })
        .collect()
}

/// [`build_components`] の並べ替えタイブレークが使う、`nConnTable` 相当の
/// 整数列の辞書式比較。「配列が長い方が先」「同じ位置なら値が大きい方が先」
/// という実 InChI (`ichimake.c::CompINChI2`) の規則を数値配列で再現する。
///
/// `nConnTable` 自体は `ichicano.c::UpdateFullLinearCT` (`CT_ATOMID_IS_CURRANK`
/// モード) が正準番号順に構築する: 各原子について「自分自身の正準番号」→
/// 「自分より小さい正準番号を持つ隣接原子 (後退辺) を昇順」の順で積む。
/// I48 時点では自分自身の番号を含めず全隣接 (前方・後退両方) を積んでいたが、
/// これは複数成分が同じ Hill 式でタイし、かつ環閉じ位置の違いで後退辺の
/// 蓄積ペースが成分ごとに異なる場合に実 InChI と順序がずれる (I53 直後の
/// PubChem `c` バケツで確認: 自分自身の番号を挟むことで配列中の対応位置が
/// ずれ、比較結果が反転するケースがあった)。
fn compare_conn_table_numeric(a: &Component, b: &Component) -> std::cmp::Ordering {
    let flat = |c: &Component| -> Vec<usize> {
        let n = c.adj.len() - 1;
        let mut out = Vec::with_capacity(c.adj.iter().map(|v| v.len()).sum::<usize>() + n);
        for i in 1..=n {
            out.push(i);
            out.extend(c.adj[i].iter().copied().filter(|&nb| nb < i));
        }
        out
    };
    let (ta, tb) = (flat(a), flat(b));
    ta.len().cmp(&tb.len()).then_with(|| ta.cmp(&tb))
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

    /// I38: 中性化でプロトンが外れた成分は、可動 H 群が 1 つに併合され
    /// 負電荷 (`-`) を共有する。共役で繋がっていない群どうしも併合される。
    #[test]
    fn mobile_charge_merges_groups() {
        // チアミン三リン酸: アミノピリミジンの N 群と三リン酸の O 群が 1 群に
        let smi = "CC1=C(SC=[N+]1CC2=CN=C(N=C2N)C)CCOP(=O)(O)OP(=O)(O)OP(=O)(O)O";
        let h = crate::inchi::inchi_of(smi).unwrap();
        assert!(
            h.contains("(H5-,13,14,15,17,18,19,20,21,22,23)/p+1"),
            "got: {h}"
        );
        // 単独の群しかなければ併合は no-op。可動 H 0 + 負電荷だけの群は
        // h 層に出さない (カルニチンの実 InChI に `(H0-,…)` は付かない)。
        let h = crate::inchi::inchi_of("C[N+](C)(C)CC(CC(=O)O)O").unwrap();
        assert_eq!(
            h,
            "InChI=1S/C7H15NO3/c1-8(2,3)5-6(9)4-7(10)11/h6,9H,4-5H2,1-3H3/p+1"
        );
    }

    /// I38: 脱プロトンの対象は酸性 O-H (フェノール等) とカルボニル型の
    /// 可動 H 群 (アミド)。塩基性アミンと単なるアルコールは外さない。
    #[test]
    fn deprotonation_site_acidity() {
        // 一級アミドから外して `(H-,8,10)/p+1`
        let h = crate::inchi::inchi_of("C[N+]1=CC=CC(=C1)C(=O)N").unwrap();
        assert!(h.contains("(H-,8,10)/p+1"), "got: {h}");
        // コリンのアルコールは外さない (I36)
        let h = crate::inchi::inchi_of("C[N+](C)(C)CCO").unwrap();
        assert!(h.ends_with("/q+1"), "got: {h}");
        // アミノチアゾリウムの塩基性 NH2 も外さない
        let h = crate::inchi::inchi_of("CC1=C(SC=[N+]1CC2=CN=C(N=C2N)C)CCO").unwrap();
        assert!(h.contains("(H2,13,14,15)/q+1"), "got: {h}");
    }

    /// I39: 酸対 (中心に二重結合 O/S ≥1 かつ 供与体 O/S ≥1) に **S を含む**
    /// 場合、N 端点は末端 NH2 (heavy_deg==1) に限らず二級 NHR も可動群に
    /// 含める。ジチオカルバミン酸型 `R-NH-C(=S)-S-`、チオカルバミン酸型
    /// `R-NH-C(=O)-S-H` はどちらも N が二級でも実 InChI は N を可動群に
    /// 含める — 通常のカルバミン酸 (酸対が O,O のみ) の N はアミド性の
    /// まま固定されるのと対照的 (PubChem 実データで確認)。
    #[test]
    fn secondary_n_joins_thio_acid_pair() {
        let h = crate::inchi::inchi_of("CNC(=S)[S-]").unwrap();
        assert!(h.contains("(H2,3,4,5)/p-1"), "got: {h}");
        let h = crate::inchi::inchi_of("CNC(=O)S").unwrap();
        assert!(h.contains("(H2,3,4,5)"), "got: {h}");
        // 対照: 酸対が O,O のみ (通常のカルバミン酸) なら N は固定のまま
        let h = crate::inchi::inchi_of("CNC(=O)O").unwrap();
        assert!(h.contains("(H,4,5)"), "got: {h}");
        assert!(!h.contains(",3,4,5)"), "N should stay fixed: {h}");
    }
}
