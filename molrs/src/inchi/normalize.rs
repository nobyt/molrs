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

/// 負電荷を中性化 (プロトン付加) できるか (I41)。
///
/// `metal_locked` に含まれる原子 (金属との異方開裂で電荷を得た O/N/S 等) は
/// 隣接原子の種類によらずプロトン化しない (`COCCO[Hg]` の実 InChI は
/// `/q-1;+1` で O をプロトン化しない)。ハロゲン化物イオンはこの対象外
/// (`C[Hg]Cl` → HCl/p-1 は検証済み、[`super::disconnect::disconnect_metals`]
/// が `metal_locked` に含めていない)。
///
/// ハロゲン化物イオンは**孤立** (重原子への結合なし) のときだけ対象
/// (InChI は孤立 `[Cl-]` を HCl/p-1 とする)。共有結合したハロゲン
/// (超原子価錯体のブロモニウム/クロロニウム型アニオン) はプロトン化しない。
///
/// O/S/Se/Te は、隣接する重原子が全て**炭素**か、二重結合 O/S を持つ
/// **酸性の S/P 中心** (スルホン酸・リン酸型) のときだけプロトン化する。
/// 隣が N (ヒドロキシルアミン型 `N-O-` 結合)、あるいは酸性でない中心
/// (亜ヒ酸 `[O-][As]([O-])[O-]` の As 等) はプロトン化しない
/// (PubChem 実データで確認)。重原子への結合が無い孤立イオンは常にプロトン化。
fn is_protonatable(
    g: &MoleculeGraph,
    atom: usize,
    metal_locked: &std::collections::HashSet<usize>,
) -> bool {
    if metal_locked.contains(&atom) {
        return false;
    }
    let sym = g.atoms[atom].symbol.as_str();
    if matches!(sym, "F" | "Cl" | "Br" | "I") {
        return !g.adjacency[atom]
            .iter()
            .any(|&nb| g.atoms[nb].symbol != "H");
    }
    if sym == "N" {
        // アミド型 N⁻ (隣接炭素が酸性中心、すなわち O/S へ実二重結合を
        // 持つ) だけプロトン化する。カルボキシラート O⁻ と違い、N⁻ は
        // 隣が「ただの炭素」というだけでは対象にしない — 単純なアルキル
        // アミンの脱プロトン化アニオン (`CC[N-]C`) は炭素だけに結合して
        // いても実 InChI は `/q-1` のまま (I41 で `is_protonatable` から
        // N を撤廃した理由)。一方でアミド型 (`C(=O)[NH-]`) は炭素の先に
        // 実二重結合 O/S があるため、下の `acidic_center` 判定 (炭素の
        // ショートカットを使わない版) でだけ拾える。
        return g.adjacency[atom]
            .iter()
            .filter(|&&nb| g.atoms[nb].symbol != "H")
            .all(|&nb| acidic_center(g, nb));
    }
    if !matches!(sym, "O" | "S" | "Se" | "Te") {
        return false;
    }
    g.adjacency[atom]
        .iter()
        .filter(|&&nb| g.atoms[nb].symbol != "H")
        .all(|&nb| protonation_neighbor_ok(g, nb, sym == "O"))
}

/// 隣が炭素で、かつ酸性 O⁻ の判定に使われるものと同じ意味で「酸性化する」
/// 炭素かどうか (芳香族 = フェノール型、二重結合を持つ = エノール/カルボ
/// ニル型)。単なる sp3 アルキル炭素は対象外 (I58、下記参照)。
fn is_acidic_carbon(g: &MoleculeGraph, c: usize) -> bool {
    g.atoms[c].is_aromatic
        || g.bonds.iter().enumerate().any(|(bi, b)| {
            (b.begin_idx == c || b.end_idx == c) && g.kekule_bond_orders[bi] > 1.0
        })
}

/// [`is_protonatable`] の隣接原子側の判定: 炭素、または二重結合 O/S を
/// 持つ**酸性中心**だけを許す。元素を限定しない — スルホン酸・リン酸型
/// (S/P) だけでなく、亜硝酸 `N(=O)[O-]`・過塩素酸 `[O-]Cl(=O)(=O)=O`・
/// 亜セレン酸 `[O-][Se](=O)[O-]`・ヒ酸 `O[As](=O)([O-])[O-]` のように
/// N/Cl/Se/As が中心でも同じ「二重結合 O/S を持つ」条件で酸性となる。
/// 二重結合を持たない中心 (亜ヒ酸の As、ヒドロキシルアミン型の N-O 結合の
/// N 等) は対象外。
///
/// 隣が O 単結合 (ペルオキシ酸型: `C(=O)O[O-]`、過ギ酸アニオンの末端 O)
/// のときは、その先 (隣の隣) に酸性中心があれば許可する (I54)。ヒドロキシ
/// ルアミン型の N-O は対象外のまま (この分岐は隣が O のときだけ発火し、
/// N のケースには触れない)。
///
/// `owner_is_o` (I58): O⁻ 自身が持つ「隣が炭素なら常に許可」のショート
/// カットは、**単なる sp3 アルキル炭素** (アルコキシド `CC[O-]`・
/// `CC(C)(C)[O-]`) には適用しない — 実 InChI はこれらを孤立イオンでも
/// プロトン化せず `/q-1` のまま残す (`inchi-1` で確認: `CC[O-].[Na+]` は
/// `/q-1;+1` で `/p` なし)。カルボキシラート・フェノラート・エノラート
/// (炭素が芳香族または二重結合を持つ) は従来どおりプロトン化する。
/// **S/Se/Te には適用しない** — チオラート `CC[S-]` は sp3 アルキル炭素
/// 隣接でも実 InChI がプロトン化する (`CC[S-].[Na+]` → `/p-1`)、O だけが
/// 例外的に非酸性アルコキシドを区別する。
fn protonation_neighbor_ok(g: &MoleculeGraph, nb: usize, owner_is_o: bool) -> bool {
    if g.atoms[nb].symbol == "C" {
        if owner_is_o {
            return is_acidic_carbon(g, nb);
        }
        return true;
    }
    if acidic_center(g, nb) {
        return true;
    }
    if g.atoms[nb].symbol == "O" {
        return g.adjacency[nb]
            .iter()
            .filter(|&&x| g.atoms[x].symbol != "H")
            .any(|&x| acidic_center(g, x));
    }
    false
}

/// [`protonation_neighbor_ok`] から「隣が炭素なら常に許す」ショートカットを
/// 除いたもの — 実二重結合 O/S を持つ酸性中心かどうかだけを見る。N⁻ の
/// プロトン化可否判定 ([`is_protonatable`]) はこちらを直接使う。
fn acidic_center(g: &MoleculeGraph, nb: usize) -> bool {
    g.bonds.iter().enumerate().any(|(bi, b)| {
        (b.begin_idx == nb || b.end_idx == nb)
            && g.kekule_bond_orders[bi] == 2.0
            && matches!(
                g.atoms[if b.begin_idx == nb {
                    b.end_idx
                } else {
                    b.begin_idx
                }]
                .symbol
                .as_str(),
                "O" | "S"
            )
    })
}

/// 陽イオンから脱プロトンして中性化できるか (塩基点)。ハロゲンは対象外。
///
/// N は例外がある: 唯一の重原子隣接が O のとき (`[NH3+]O`・`[NH3+]OC` の
/// ような N-O-R 型ヒドロキシルアミニウム) は脱プロトンしない。これは
/// [`protonation_neighbor_ok`] が O⁻ 側で「隣が N (ヒドロキシルアミン型
/// N-O 結合) ならプロトン化しない」としているのと同じ N-O 境界を逆方向
/// (N 側が陽電荷) にも適用したもの — `inchi-1` で確認: `[NH3+]O` は
/// `H4NO/q+1` のまま (`/p` を使わない) だが、`[NH3+]N` (ヒドラジニウム)・
/// `C[NH2+]O` (N が O 以外にも重原子隣接を持つ) はどちらも通常どおり
/// `/p+1` で脱プロトンする。
fn is_deprotonatable(g: &MoleculeGraph, atom: usize) -> bool {
    let sym = g.atoms[atom].symbol.as_str();
    if !matches!(sym, "N" | "O" | "S" | "Se" | "Te") {
        return false;
    }
    if sym == "N" {
        let mut heavy_neighbors = g.adjacency[atom]
            .iter()
            .copied()
            .filter(|&nb| g.atoms[nb].symbol != "H");
        if let Some(only) = heavy_neighbors.next() {
            if heavy_neighbors.next().is_none() && g.atoms[only].symbol == "O" {
                return false;
            }
        }
    }
    true
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
/// `metal_locked` は [`super::disconnect::disconnect_metals`] が返す、
/// 金属由来の電荷を恒久として扱うべき原子集合 ([`is_protonatable`] 参照)。
pub(crate) fn neutralize(
    g: &MoleculeGraph,
    metal_locked: &std::collections::HashSet<usize>,
) -> (MoleculeGraph, i32, i32) {
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

    // アミジニウム/グアニジニウム型の陽電荷リレー (I41)。
    //
    // `[N+]1=C(...)N` (環内 N+ が二重結合、環外 NH2 が単結合) のような
    // 中心は、電荷を持つ N+ 自身は H を持たないが、可動 H 群のもう一方の
    // 端点 (環外 NH2) が H を持つ。実 InChI はこの群全体を 1 つのアミジン
    // として中性化する — N+ の電荷と NH2 の 1 個の H が**同時に**消え、
    // どちらの原子にも電荷は残らない (中性のアミジン `N=C-NH` になる)。
    //
    // これは I38 の「酸点をまたぐ可動負電荷」(無関係な外部の陽イオンが
    // 別の酸性 O-H から借りる) とは逆方向で、**群自身が電荷を持つ**場合
    // だけに限る。無関係な陽イオンが塩基性のアミジン/グアニジン基から
    // 借りることはない (I38 で確認済み、チアミンのアミノピリミジンが
    // `/q+1` のまま残るのがその証拠) ので、「群のメンバー自身が正電荷」
    // という条件だけで安全に区別できる。
    //
    // ```
    // smi  CCC(=O)C[N+]1=C(OC2=C1CCCC2)N
    // want …h12H,…/p+1                     ← NH2 が NH になり N+ は消える
    // ```
    let mut relay_donor: HashMap<usize, usize> = HashMap::new();
    for (eps, mh, _) in super::number::mobile_groups(g) {
        if mh == 0 {
            continue;
        }
        let Some(&donor) = eps.iter().find(|&&e| cur_h(e) > 0) else {
            continue;
        };
        for &e in &eps {
            if e != donor && g.atoms[e].formal_charge > 0 {
                relay_donor.insert(e, donor);
            }
        }
    }
    let mut relay_debt = vec![0i32; n_heavy];

    for i in 0..n_heavy {
        let a = &g.atoms[i];
        let h = cur_h(i);
        let ch = a.formal_charge as i32;
        let b = &mut budget[comp_of[i]];
        if ch != 0 && has_opposite_charged_neighbor(i) {
            // 電荷分離 (ニトロ・N-オキシド等) → 触らない
            final_h[i] = h;
            new_charge[i] = a.formal_charge;
        } else if ch < 0 && is_protonatable(g, i, metal_locked) && *b < 0 {
            // 負電荷 → プロトン付加で中性化 (成分の正味電荷を超えない範囲で)
            let add = (-ch).min(-*b);
            n_add += add;
            *b += add;
            final_h[i] = h + add;
            new_charge[i] = (ch + add) as i8;
        } else if ch > 0 && is_deprotonatable(g, i) && h > 0 && *b > 0 {
            // プロトン付き陽イオン → 除去で中性化
            let rem = ch.min(h).min(*b);
            n_remove += rem;
            *b -= rem;
            final_h[i] = h - rem;
            new_charge[i] = (ch - rem) as i8;
        } else if ch > 0 && h == 0 && *b > 0 && relay_donor.contains_key(&i) {
            // アミジニウム/グアニジニウム型リレー: 自分自身は H を持たないが
            // 可動 H 群のドナー端点が H を持つ。電荷はここで消し、H の除去は
            // ドナー側に借りとして記録し、ループ後にまとめて適用する。
            let rem = ch.min(*b);
            n_remove += rem;
            *b -= rem;
            final_h[i] = h;
            new_charge[i] = (ch - rem) as i8;
            relay_debt[relay_donor[&i]] += rem;
        } else {
            // 中性化不能 (四級 N・金属など)、または成分が既に正味中性
            final_h[i] = h;
            new_charge[i] = a.formal_charge;
        }
    }
    // リレーの借りをドナー側に適用 (訪問順に依存しないよう、ここでまとめて)
    for (i, &debt) in relay_debt.iter().enumerate() {
        final_h[i] -= debt;
    }

    // 第 2 パス (I31): 成分にまだ正味の陽電荷が残っていて、外せる**陽イオンの**
    // H が無い場合は、**中性の酸点** (カルボン酸等の O-H) からプロトンを外して
    // 釣り合わせる。
    //
    // カルニチン `C[N+](C)(C)CC(CC(=O)O)O` がこれで、四級 N+ は脱プロトン
    // できないが実 InChI は COOH から 1 個外した双性イオンを基準にして
    // `C7H15NO3/…/p+1` を出す。外さないと `/q+1` になってしまう。
    //
    // **酸性の O-H/S-H** が第 1 候補。単なるアルコールを外してはいけない
    // (I35 で判明)。コリン `C[N+](C)(C)CCO` の実 InChI は `C5H14NO/…/q+1` で、
    // OH の H を保ったまま `/q` に残す。ここを「任意の O-H」にしていると
    // `/p+1` になってしまう。
    //
    // I36 では「隣の重原子が O/S へ二重結合を持つ」(カルボン酸・スルホン酸・
    // リン酸型) に限っていたが、実 InChI はフェノール・エノール・ヒドロキサム酸の
    // N-O-H からも外す。糖・アルコールとの境界は「**飽和炭素に付いた O-H か
    // どうか**」— 芳香族炭素 (フェノール)、二重結合を持つ炭素 (エノール)、
    // 炭素ですらない相手 (N-OH) はいずれも酸性側に入る (I38)。
    //
    // ```text
    // smi  CN1CCC2=CC(=…)…[N+](CCC…)(C)C…   ← 四級 N+ + カテコール
    // want …/h6-11,18-21,28-29H,…,(H-,40,41)/p+1   ← フェノールから外す
    // ```
    let is_acidic_oh = |i: usize| -> bool {
        g.adjacency[i].iter().any(|&c| {
            let a = &g.atoms[c];
            if a.symbol == "H" {
                return false;
            }
            // ホウ酸型 (B-OH) は「非炭素」の例外に含めない — ボロン酸
            // `B(O)O` の O-H は無関係な陽電荷 (第四級アンモニウム等) の
            // 相手にはならず、実 InChI は中性のまま `/q+1` に残す
            // (I57、`c1ccc(C[n+]...)cc1B(O)O` で確認)。
            //
            // 隣が**それ自身が陽電荷を持つ N** (ヒドロキシルアミニウム型
            // `[NH3+]O`) のときも除外する — ヒドロキサム酸 `C(=O)NO` の
            // ように「中性の N を介して離れた四級 N+ から借りる」場合とは
            // 違い、この O-H はその N 自身の電荷を中和する相手そのものに
            // なってしまう。[`is_deprotonatable`] が N 側で同じ N-O 境界を
            // 理由に脱プロトンを拒否しているのと対称に、O 側でも同じ
            // 「その N 自身の電荷」を借りない (I61、`inchi-1` で確認:
            // `[NH3+]O` は `/p` を使わず `H4NO/q+1` のまま)。
            (a.symbol != "C" && a.symbol != "B" && !(a.symbol == "N" && a.formal_charge > 0))
                || a.is_aromatic
                || g.bonds.iter().enumerate().any(|(bi, b)| {
                    (b.begin_idx == c || b.end_idx == c) && g.kekule_bond_orders[bi] > 1.0
                })
        })
    };
    // 第 2 候補: **カルボニル型の可動 H 群** (アミド・イミド・尿素) の H。
    // ニコチンアミド `C[N+]1=CC=CC(=C1)C(=O)N` の実 InChI は `(H-,8,10)/p+1` で、
    // 一級アミドから外している。酸点らしい O-H が無くても四級 N+ は中和される。
    //
    // 「可動 H 群ならどれでも」にすると外し過ぎる。**塩基性の**アミン群 —
    // アミノピリミジン (チアミン `CC1=C(SC=[N+]1CC2=CN=C(N=C2N)C)CCO` の実 InChI は
    // `/q+1`)、アミジン、アミノチアゾール、グアニジン — は外さない。一級
    // スルホンアミド `S(=O)(=O)N` も外さない。境界は「群に**炭素に付いた**
    // O/S 端点があるか」で、これがカルボニル (C=O) を持つ群だけを選び出す
    // (スルホンアミドの O は S に付いているので落ちる)。
    //
    // 併合後は群全体で 1 個の負電荷を共有するので、群内のどのメンバーを
    // 選んでも結果は同じ。
    let mut acidic_taut: Option<Vec<bool>> = None;
    for (c, left) in budget.iter_mut().enumerate().take(n_comps) {
        while *left > 0 {
            let has_h = |i: usize| comp_of[i] == c && new_charge[i] == 0 && final_h[i] > 0;
            let pick = (0..n_heavy)
                .filter(|&i| {
                    has_h(i) && matches!(g.atoms[i].symbol.as_str(), "O" | "S") && is_acidic_oh(i)
                })
                .min_by_key(|&i| (g.atoms[i].symbol != "O", i))
                .or_else(|| {
                    let tm = acidic_taut.get_or_insert_with(|| {
                        let mut m = vec![false; n_heavy];
                        for (eps, _, _) in super::number::mobile_groups(g) {
                            let carbonyl = eps.iter().any(|&e| {
                                matches!(g.atoms[e].symbol.as_str(), "O" | "S")
                                    && g.adjacency[e].iter().any(|&c| g.atoms[c].symbol == "C")
                            });
                            if carbonyl {
                                for e in eps {
                                    m[e] = true;
                                }
                            }
                        }
                        m
                    });
                    (0..n_heavy).find(|&i| has_h(i) && tm[i])
                });
            let Some(i) = pick else { break };
            final_h[i] -= 1;
            new_charge[i] = -1;
            n_remove += 1;
            *left -= 1;
        }
    }
    // 孤立した `[H+]` (どの重原子にも結合しない、正電荷を持つ H 単体) は
    // 実 InChI では「どこかの酸性中心が実際にプロトン化されるための実体」
    // として直接消費される — 上のパスで数えた `n_add` (中性化のために
    // "仮想的に" 足したプロトン数) のうち、実際にこの孤立 H+ の数だけは
    // 実在の原子で賄えるので `/p` に計上しない。孤立 H+ を消費せず単なる
    // 別成分として残すと `.H` という余分な成分と `+1;/p-1` が付いてしまう
    // (`[H+].c1ccccc1N.[Cl-]` の実 InChI は `C6H8ClN` 相当を 2 成分
    // `C6H7N.ClH` にまとめるだけで `/q`・`/p` を一切出さない)。
    //
    // 消費しきれない余り (孤立 H+ の数 > n_add) はそのまま孤立成分として
    // 残す (下の孤立 H 保持ループ) — 対応する酸性中心が無い場合に相当し、
    // 実 InChI での挙動は未検証だが少なくとも退行はしない安全側の扱い。
    let mut consumed_lone_h_plus: std::collections::HashSet<usize> = std::collections::HashSet::new();
    if n_add > 0 {
        for old in 0..g.atoms.len() {
            if consumed_lone_h_plus.len() as i32 >= n_add {
                break;
            }
            let is_lone_h = g.atoms[old].symbol == "H"
                && !g.adjacency[old].iter().any(|&nb| g.atoms[nb].symbol != "H");
            if is_lone_h && g.atoms[old].formal_charge == 1 {
                consumed_lone_h_plus.insert(old);
            }
        }
        n_add -= consumed_lone_h_plus.len() as i32;
    }

    // q は最終的な電荷の総和 (2 パス目で変わるのでここで数え直す)
    let q: i32 = new_charge.iter().map(|&c| c as i32).sum();

    let p = n_remove - n_add;
    if n_add == 0
        && n_remove == 0
        && consumed_lone_h_plus.is_empty()
        && q == g.atoms.iter().map(|a| a.formal_charge as i32).sum()
    {
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
    // H ノードを再付加。`crate::stereo::build_ctx` は「明示的に解析された原子
    // (`parser_to_graph` に対応あり) はグラフ先頭の連続範囲を占める」という
    // 不変条件に依存している (同位体を持ちうる原子だけを CIP ランクのノードに
    // 含め、残りの暗黙 H は `n_h` カウントに畳み込むため)。そのため単純に
    // 「重原子の出現順」で再構築すると、明示同位体 H (`[2H]` 等) が暗黙 H の
    // 間に紛れ込み、その不変条件と `parser_to_graph` の対応 (`h_remap`) の
    // 両方が壊れる。ここでは 2 パスに分け、明示原子を先に (元の解析順を保った
    // まま) 重原子の直後へ詰め、暗黙 H は必ず末尾に回す。
    let mut h_remap: HashMap<usize, usize> = HashMap::new();
    let mut remaining: Vec<i32> = final_h.clone();
    // パス 1: 明示的に解析された H (同位体標識の有無を問わない) を、元の解析順
    // のまま先頭寄りに詰め直す。`fh` を超えた分 (脱プロトン等で消えた原子) は
    // 対応なしのまま置き去りにする — `h_remap` に入らないので、後段の
    // `parser_to_graph` 再構築で自動的に `None` になる。
    for slot in &g.parser_to_graph {
        let Some(old_gi) = *slot else { continue };
        if old_gi >= g.atoms.len() || g.atoms[old_gi].symbol != "H" {
            continue;
        }
        let Some(&owner) = g.adjacency[old_gi].iter().find(|&&nb| nb < n_heavy) else {
            continue; // 孤立 H (重原子に結合しない) は別扱い
        };
        if remaining[owner] <= 0 {
            continue;
        }
        remaining[owner] -= 1;
        let h_idx = atoms.len();
        h_remap.insert(old_gi, h_idx);
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
            begin_idx: owner,
            end_idx: h_idx,
            bond_order: 1.0,
            stereo: None,
        });
        kekule.push(1.0);
    }
    // 金属切断で生じた孤立 H (どの重原子にも結合しない) は独立成分として
    // 保持する必要がある (`[BiH3]` → `Bi.3H`)。パス 1 は重原子に結合した H
    // しか再生成しないため、ここで補う (I20)。これも「明示的に解析された
    // 原子」の一種なので、暗黙 H を埋めるパス 2 より**前**に置いて、明示原子の
    // 連続範囲 (`build_ctx` の不変条件) を維持する。
    for old in 0..g.atoms.len() {
        let is_lone_h = g.atoms[old].symbol == "H"
            && !g.adjacency[old].iter().any(|&nb| g.atoms[nb].symbol != "H");
        if !is_lone_h || consumed_lone_h_plus.contains(&old) {
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
            // 金属水素化物切断由来の孤立 H は元々中性だが、`[H-]` のように
            // 入力 SMILES 自体が孤立 H に電荷を持たせている場合もある
            // (`[H-].C(=O)O[O-].[K+]` 等) — 電荷を保持しないと /q 層から
            // 丸ごと消えてしまう。
            formal_charge: g.atoms[old].formal_charge,
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
    // パス 2: 残り (真に新規のプロトン、または暗黙 H) を重原子順に追加。明示
    // 原子 (パス 1 + 孤立 H) より必ず後ろに来るようにする。
    for (i, &rem) in remaining.iter().enumerate() {
        for _ in 0..rem.max(0) {
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

    // parser_to_graph も H の再番号に合わせて付け替える。重原子はインデックス
    // 不変なのでそのまま、H は `h_remap` にあれば新番号、無ければ (プロトン化/
    // 脱プロトンで消えた原子) 対応なしとして None にする — さもないと
    // 同位体層 (`/i`) が古い番号のまま別の H を指してしまう (I63)。
    let parser_to_graph: Vec<Option<usize>> = g
        .parser_to_graph
        .iter()
        .map(|old| {
            old.and_then(|old_gi| {
                if old_gi < n_heavy {
                    Some(old_gi)
                } else {
                    h_remap.get(&old_gi).copied()
                }
            })
        })
        .collect();

    let ng = MoleculeGraph {
        atoms,
        bonds,
        adjacency,
        bond_orders,
        ring_atom_sets: g.ring_atom_sets.clone(),
        kekule_bond_orders: kekule,
        parsed: g.parsed.clone(),
        parser_to_graph,
    };
    (ng, q, p)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    fn qp(smiles: &str) -> (i32, i32, String) {
        let g = build_molecule_graph(smiles).unwrap();
        let (ng, q, p) = neutralize(&g, &std::collections::HashSet::new());
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

    /// I41: 負電荷の隣が N (ヒドロキシルアミン型 N-O 結合) や非酸性中心
    /// (亜ヒ酸の As) だとプロトン化しない。共有結合したハロゲン化物
    /// (ブロモニウム型) も同様。一方、酸性中心が S/P 以外 (亜硝酸の N、
    /// 過塩素酸の Cl、亜セレン酸の Se、ヒ酸の As) でも二重結合 O/S さえ
    /// あれば通常どおりプロトン化する。
    #[test]
    fn protonation_neighbor_matters() {
        // ヒドロキシルアミン型: N-O 結合の O はプロトン化しない
        let (q, _, f) = qp("N[O-]");
        assert_eq!(q, -1, "N[O-]: {f}");
        // 亜ヒ酸: As が二重結合を持たない → プロトン化しない
        let (q, _, _) = qp("[O-][As]([O-])[O-]");
        assert_eq!(q, -3);
        // ヒ酸: As が二重結合 O を持つ → 通常どおりプロトン化
        let (q, _, f) = qp("O[As](=O)([O-])[O-]");
        assert_eq!(q, 0, "arsenic acid: {f}");
        assert_eq!(f, "AsH3O4");
        // 亜硝酸・過塩素酸も同様に酸性中心として扱う (中心 N/Cl 自身に
        // 二重結合 O/S があれば、その中心が S/P でなくてもプロトン化する)
        let (q, _, f) = qp("N(=O)[O-]");
        assert_eq!(q, 0, "nitrite: {f}");
        let (q, _, f) = qp("[O-]Cl(=O)(=O)=O");
        assert_eq!(q, 0, "perchlorate: {f}");
    }

    /// I41: N⁻ (脱プロトン化アミン) は炭素に結合していても通常プロトン化しない
    /// — カルボキシラート/スルホン酸型の O⁻/S⁻ とは異なり、実 InChI は
    /// アミド/アミン型アニオンを恒久電荷のまま残す。
    #[test]
    fn amine_anion_stays_charged() {
        let (q, _, _) = qp("CC[NH-]");
        assert_eq!(q, -1);
    }

    /// I41: アミジニウム/グアニジニウム型リレー。環内 N+ 自身は H を持たない
    /// が、可動 H 群のもう一方の端点 (環外 NH2) から 1 個を代わりに除去し、
    /// N+ の電荷も同時に消える (`/p+1`、恒久電荷は残らない)。
    #[test]
    fn amidinium_relay_removes_group_h() {
        let g = build_molecule_graph("CCC(=O)C[N+]1=C(OC2=C1CCCC2)N").unwrap();
        let h = crate::inchi::to_inchi(&g).unwrap();
        assert!(h.ends_with("/p+1"), "got: {h}");
        assert!(!h.contains("/q"), "got: {h}");
    }

    /// I41: 金属との異方開裂で電荷を得た O/N/S 等 (ハロゲン以外) は
    /// プロトン化せず恒久電荷として残す (`COCCO[Hg]` の実 InChI は
    /// `/q-1;+1`)。ハロゲン化物 (`C[Hg]Cl`) は対象外で通常どおり
    /// プロトン化する。
    #[test]
    fn metal_locked_heteroatom_stays_charged() {
        let h = crate::inchi::inchi_of("COCCO[Hg]").unwrap();
        assert!(h.contains("/q-1;+1"), "got: {h}");
        let h = crate::inchi::inchi_of("C[Hg]Cl").unwrap();
        assert!(h.contains("/p-1"), "got: {h}");
    }

    /// I61: ヒドロキシルアミニウム型 (N が唯一の重原子隣接として O だけを
    /// 持つ陽イオン `[NH3+]O`・`[NH3+]OC`) は脱プロトンせず `/q` のまま
    /// 残す。N が O 以外の重原子 (C・N) にも結合していれば通常どおり
    /// `/p` で脱プロトンする。
    #[test]
    fn hydroxylammonium_stays_charged() {
        let h = crate::inchi::inchi_of("[NH3+]O").unwrap();
        assert_eq!(h, "InChI=1S/H4NO/c1-2/h2H,1H3/q+1");
        let h = crate::inchi::inchi_of("[NH3+]OC").unwrap();
        assert_eq!(h, "InChI=1S/CH6NO/c1-3-2/h1-2H3/q+1");
        // N が O 以外にも重原子隣接を持てば脱プロトンする
        let h = crate::inchi::inchi_of("C[NH2+]O").unwrap();
        assert!(h.ends_with("/p+1"), "got: {h}");
        let h = crate::inchi::inchi_of("[NH3+]N").unwrap();
        assert!(h.ends_with("/p+1"), "got: {h}");
    }
}

