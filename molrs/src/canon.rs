//! 正規 SMILES 生成とフラグメント分解 (S1.5)。
//!
//! RDKit の正規形との文字列一致は**要求しない** (RUST_PORT_PLAN.md 基本方針 5)。
//! 要件は 2 つ:
//! - 決定性: 同じ入力グラフから常に同じ文字列
//! - 不変性: 同一分子のどの SMILES 表記から構築しても同じ正規形
//!
//! アルゴリズム: Morgan/CANGEN 系の不変量反復精緻化で原子を順位付けし、
//! 同値類が残る場合は各メンバーを強制的に先頭に置いて分岐し、
//! 得られる SMILES 文字列の辞書順最小を採用する (健全な正準化)。
//!
//! 立体化学 (@/@@, /\) は未出力 (S1.7 で対応予定)。
//!
//! 2 つの書き出しモード:
//! - Strict: 完全な分子の正規形。再パースで同一グラフに戻る
//!   (電荷・同位体・推論不能な H 数は角括弧で明示)
//! - Lenient: 部分構造キー用 (RDKit `MolFragmentToSmiles` 相当)。
//!   有機サブセット原子は H 数が合わなくても裸で書く (置換位置の
//!   コアパターンが独立分子のパターンと同じ文字列になる)

use crate::graph::{aromatic_takes_double_bond, default_valences, MoleculeGraph};

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum WriteMode {
    Strict,
    Lenient,
}

/// 書き出し用の原子情報。
struct CAtom {
    symbol: String,
    aromatic: bool,
    charge: i8,
    isotope: Option<u16>,
    n_h: u8,
    atomic_num: u8,
    /// CIP コード (0=なし, 1=R, 2=S)。タイブレーク不変量にのみ使い、
    /// 出力文字列そのものは変えない (立体は複合キー側で表現する)。
    /// メソ体などで対称な 2 中心の出力位置が入力表記に依らず安定する。
    stereo_ord: u8,
}

/// 書き出し用のローカルグラフ (subset 内のインデックスに再番号済み)。
struct CMol {
    atoms: Vec<CAtom>,
    /// (a, b, order_class) order_class: 0=単, 1=芳香族, 2=二重, 3=三重, 4=四重
    bonds: Vec<(usize, usize, u8)>,
    adj: Vec<Vec<(usize, usize)>>, // atom → [(相手, bond_idx)]
}

fn order_class(order: f64) -> u8 {
    // f64 は match パターンにできないため段階比較
    if order == 1.5 {
        1
    } else if order == 2.0 {
        2
    } else if order == 3.0 {
        3
    } else if order == 4.0 {
        4
    } else {
        0
    }
}

/// グラフの「保持原子」(付加 H を除く) の数。
fn num_kept_atoms(g: &MoleculeGraph) -> usize {
    g.parser_to_graph.iter().flatten().count()
}

/// MoleculeGraph の部分集合から CMol を作る。
fn build_cmol(g: &MoleculeGraph, subset: &[usize]) -> CMol {
    let n_kept = num_kept_atoms(g);
    // graph idx → parser idx (同位体の取得用)
    let mut graph_to_parser = vec![usize::MAX; n_kept];
    for (pi, slot) in g.parser_to_graph.iter().enumerate() {
        if let Some(gi) = slot {
            graph_to_parser[*gi] = pi;
        }
    }

    let mut local = vec![usize::MAX; n_kept];
    let mut atoms = Vec::with_capacity(subset.len());
    for (li, &gi) in subset.iter().enumerate() {
        local[gi] = li;
        let a = &g.atoms[gi];
        // 付加 H (idx >= n_kept) の数 = この原子の H 数
        let n_h = g.adjacency[gi].iter().filter(|&&x| x >= n_kept).count() as u8;
        atoms.push(CAtom {
            symbol: a.symbol.clone(),
            aromatic: a.is_aromatic,
            charge: a.formal_charge,
            isotope: g.parsed.atoms[graph_to_parser[gi]].isotope,
            n_h,
            atomic_num: a.atomic_num,
            stereo_ord: match a.chiral_tag {
                Some('R') => 1,
                Some('S') => 2,
                _ => 0,
            },
        });
    }

    let mut bonds = Vec::new();
    for b in &g.bonds {
        if b.begin_idx >= n_kept || b.end_idx >= n_kept {
            continue; // 付加 H への結合
        }
        let (la, lb) = (local[b.begin_idx], local[b.end_idx]);
        if la == usize::MAX || lb == usize::MAX {
            continue; // subset 外
        }
        bonds.push((la, lb, order_class(b.bond_order)));
    }
    let mut adj = vec![Vec::new(); atoms.len()];
    for (ei, &(a, b, _)) in bonds.iter().enumerate() {
        adj[a].push((b, ei));
        adj[b].push((a, ei));
    }
    CMol { atoms, bonds, adj }
}

/// 不変量の反復精緻化。ranks はクラス id (0 始まり)。クラス数を返す。
fn refine(mol: &CMol, ranks: &mut [usize]) -> usize {
    let n = mol.atoms.len();
    let mut n_classes = ranks.iter().max().map_or(0, |m| m + 1);
    loop {
        // key = (自クラス, ソート済み [(結合クラス, 相手クラス)])
        let mut keys: Vec<(usize, Vec<(u8, usize)>)> = Vec::with_capacity(n);
        for i in 0..n {
            let mut nbrs: Vec<(u8, usize)> = mol.adj[i]
                .iter()
                .map(|&(j, ei)| (mol.bonds[ei].2, ranks[j]))
                .collect();
            nbrs.sort_unstable();
            keys.push((ranks[i], nbrs));
        }
        let mut sorted: Vec<&(usize, Vec<(u8, usize)>)> = keys.iter().collect();
        sorted.sort();
        sorted.dedup();
        let new_n = sorted.len();
        for i in 0..n {
            ranks[i] = sorted.binary_search(&&keys[i]).expect("key exists");
        }
        if new_n == n_classes {
            return n_classes;
        }
        n_classes = new_n;
    }
}

/// 初期不変量から正準順位を作り、同値類が残れば分岐して
/// 辞書順最小の SMILES を返す。第 2 要素は原子の出力順 (ローカル idx)。
fn canonical_component_smiles(mol: &CMol, mode: WriteMode) -> (String, Vec<usize>) {
    let n = mol.atoms.len();
    let mut init: Vec<(u8, bool, i8, u8, usize, u16, u8)> = Vec::with_capacity(n);
    for (i, a) in mol.atoms.iter().enumerate() {
        init.push((
            a.atomic_num,
            a.aromatic,
            a.charge,
            a.n_h,
            mol.adj[i].len(),
            a.isotope.unwrap_or(0),
            a.stereo_ord,
        ));
    }
    let mut sorted: Vec<_> = init.clone();
    sorted.sort_unstable();
    sorted.dedup();
    let mut ranks: Vec<usize> = init
        .iter()
        .map(|k| sorted.binary_search(k).expect("key"))
        .collect();

    let mut budget = 2000usize; // 分岐の安全上限
    resolve(mol, &mut ranks, mode, &mut budget)
}

fn resolve(
    mol: &CMol,
    ranks: &mut [usize],
    mode: WriteMode,
    budget: &mut usize,
) -> (String, Vec<usize>) {
    let n = mol.atoms.len();
    let n_classes = refine(mol, ranks);
    if n_classes == n || *budget == 0 {
        return write_smiles(mol, ranks, mode);
    }
    // 最小ランクの非単一クラスを見つけ、各メンバーで分岐
    let mut class_size = vec![0usize; n_classes];
    for &r in ranks.iter() {
        class_size[r] += 1;
    }
    let target = (0..n_classes)
        .find(|&c| class_size[c] > 1)
        .expect("tied class exists");
    let members: Vec<usize> = (0..n).filter(|&i| ranks[i] == target).collect();

    let mut best: Option<(String, Vec<usize>)> = None;
    for &m in &members {
        if *budget == 0 {
            break;
        }
        *budget -= 1;
        // m を先頭に強制: 他の全クラスを 1 つずらして m を単独クラスに
        let mut forced: Vec<usize> = ranks.iter().map(|&r| r * 2 + 1).collect();
        forced[m] = target * 2;
        let s = resolve(mol, &mut forced, mode, budget);
        if best.as_ref().is_none_or(|b| s.0 < b.0) {
            best = Some(s);
        }
    }
    best.expect("at least one branch")
}

/// 原子 1 つ分のトークンを書く。
fn write_atom(a: &CAtom, mol: &CMol, i: usize, mode: WriteMode) -> String {
    let organic = matches!(
        a.symbol.as_str(),
        "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
    );
    let aromatic_organic = matches!(a.symbol.as_str(), "B" | "C" | "N" | "O" | "P" | "S");
    let sym = if a.aromatic {
        a.symbol.to_lowercase()
    } else {
        a.symbol.clone()
    };

    let mut needs_bracket = a.charge != 0
        || a.isotope.is_some()
        || a.symbol == "H"
        || !organic
        || (a.aromatic && !aromatic_organic);

    if !needs_bracket && a.symbol != "*" {
        // H 数が SMILES リーダの推論と一致しないなら角括弧が必要。
        // Lenient (部分構造キー) では有機サブセットの H 不一致を無視するが、
        // 芳香族 N/P (ピロール型 [nH]) のみ H の有無がパターンの意味を
        // 変えるため常に判定する (芳香族 C は裸のまま)。
        let check = match mode {
            WriteMode::Strict => true,
            WriteMode::Lenient => a.aromatic && matches!(a.symbol.as_str(), "N" | "P"),
        };
        if check {
            let mut base = 0usize;
            for &(_, ei) in &mol.adj[i] {
                base += match mol.bonds[ei].2 {
                    2 => 2,
                    3 => 3,
                    4 => 4,
                    _ => 1, // 単結合・芳香族
                };
            }
            let v = if a.aromatic && aromatic_takes_double_bond(&a.symbol) {
                base + 1
            } else {
                base
            };
            let inferred = default_valences(&a.symbol)
                .iter()
                .find(|&&t| t as usize >= v)
                .map(|&t| t as usize - v);
            let inferred = match inferred {
                Some(h) => h,
                None if a.aromatic => 0,
                None => usize::MAX, // 原子価超過: 角括弧で明示するしかない
            };
            if inferred != a.n_h as usize {
                needs_bracket = true;
            }
        }
    }

    if !needs_bracket {
        return sym;
    }
    let mut s = String::from("[");
    if let Some(iso) = a.isotope {
        s.push_str(&iso.to_string());
    }
    s.push_str(&sym);
    if a.n_h == 1 {
        s.push('H');
    } else if a.n_h > 1 {
        s.push('H');
        s.push_str(&a.n_h.to_string());
    }
    match a.charge {
        0 => {}
        1 => s.push('+'),
        -1 => s.push('-'),
        c if c > 0 => s.push_str(&format!("+{c}")),
        c => s.push_str(&format!("-{}", -c)),
    }
    s.push(']');
    s
}

/// 結合トークン (省略可能なら空文字)。
fn bond_token(mol: &CMol, ei: usize) -> &'static str {
    let (a, b, oc) = mol.bonds[ei];
    match oc {
        0 => {
            // 芳香族原子同士の単結合は '-' を明示 (省略すると芳香族結合扱いになる)
            if mol.atoms[a].aromatic && mol.atoms[b].aromatic {
                "-"
            } else {
                ""
            }
        }
        1 => "", // 芳香族結合は省略
        2 => "=",
        3 => "#",
        4 => "$",
        _ => unreachable!(),
    }
}

/// 与えられた完全順位で SMILES を書き出す。第 2 要素は原子の出力順。
fn write_smiles(mol: &CMol, ranks: &[usize], mode: WriteMode) -> (String, Vec<usize>) {
    let n = mol.atoms.len();
    let root = (0..n).min_by_key(|&i| ranks[i]).expect("nonempty");

    // DFS 木と後退辺 (環閉じ) を決める。子は順位昇順。(分子は小さいので再帰で十分)
    let mut visited = vec![false; n];
    let mut ring_bonds_at: Vec<Vec<usize>> = vec![Vec::new(); n]; // atom → 環閉じ bond idx
    let mut tree_children: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n]; // (child, bond)
    {
        fn dfs(
            u: usize,
            parent_bond: Option<usize>,
            mol: &CMol,
            ranks: &[usize],
            visited: &mut [bool],
            tree_children: &mut [Vec<(usize, usize)>],
            ring_bonds_at: &mut [Vec<usize>],
        ) {
            let mut nbrs: Vec<(usize, usize)> = mol.adj[u].clone();
            nbrs.sort_by_key(|&(v, _)| ranks[v]);
            for (v, ei) in nbrs {
                if Some(ei) == parent_bond {
                    continue;
                }
                if visited[v] {
                    // 後退辺: 両端に記録 (u 側が閉じ、v 側が開き)
                    if !ring_bonds_at[u].contains(&ei) {
                        ring_bonds_at[u].push(ei);
                        ring_bonds_at[v].push(ei);
                    }
                } else {
                    visited[v] = true;
                    tree_children[u].push((v, ei));
                    dfs(
                        v,
                        Some(ei),
                        mol,
                        ranks,
                        visited,
                        tree_children,
                        ring_bonds_at,
                    );
                }
            }
        }
        visited[root] = true;
        dfs(
            root,
            None,
            mol,
            ranks,
            &mut visited,
            &mut tree_children,
            &mut ring_bonds_at,
        );
    }

    // 環閉じ番号の割当て: 書き出し順に開き、閉じたら番号を再利用
    let mut digit_of_bond: Vec<Option<u16>> = vec![None; mol.bonds.len()];
    let mut used_digits = vec![false; 100];
    let mut out = String::new();
    let mut order: Vec<usize> = Vec::with_capacity(n);

    #[allow(clippy::too_many_arguments)]
    fn write_rec(
        u: usize,
        mol: &CMol,
        mode: WriteMode,
        tree_children: &[Vec<(usize, usize)>],
        ring_bonds_at: &[Vec<usize>],
        digit_of_bond: &mut [Option<u16>],
        used_digits: &mut [bool],
        out: &mut String,
        order: &mut Vec<usize>,
    ) {
        order.push(u);
        out.push_str(&write_atom(&mol.atoms[u], mol, u, mode));
        // 環閉じ数字
        for &ei in &ring_bonds_at[u] {
            match digit_of_bond[ei] {
                None => {
                    // 開き: 最小の未使用番号
                    let d = (1..100).find(|&d| !used_digits[d]).expect("digit") as u16;
                    used_digits[d as usize] = true;
                    digit_of_bond[ei] = Some(d);
                    out.push_str(bond_token(mol, ei));
                    push_digit(out, d);
                }
                Some(d) => {
                    used_digits[d as usize] = false;
                    out.push_str(bond_token(mol, ei));
                    push_digit(out, d);
                }
            }
        }
        // 子: 最後以外は括弧
        let children = &tree_children[u];
        for (k, &(v, ei)) in children.iter().enumerate() {
            let last = k == children.len() - 1;
            if !last {
                out.push('(');
            }
            out.push_str(bond_token(mol, ei));
            write_rec(
                v,
                mol,
                mode,
                tree_children,
                ring_bonds_at,
                digit_of_bond,
                used_digits,
                out,
                order,
            );
            if !last {
                out.push(')');
            }
        }
    }
    write_rec(
        root,
        mol,
        mode,
        &tree_children,
        &ring_bonds_at,
        &mut digit_of_bond,
        &mut used_digits,
        &mut out,
        &mut order,
    );
    (out, order)
}

fn push_digit(out: &mut String, d: u16) {
    if d < 10 {
        out.push_str(&d.to_string());
    } else {
        out.push('%');
        out.push_str(&format!("{d:02}"));
    }
}

/// 保持原子 (付加 H を除く) の連結成分。原子は昇順、成分は最小原子順。
/// (RDKit `GetMolFrags` 相当)
pub fn get_fragments(g: &MoleculeGraph) -> Vec<Vec<usize>> {
    let n_kept = num_kept_atoms(g);
    let mut comp = vec![usize::MAX; n_kept];
    let mut n_comp = 0;
    for start in 0..n_kept {
        if comp[start] != usize::MAX {
            continue;
        }
        let mut stack = vec![start];
        comp[start] = n_comp;
        while let Some(u) = stack.pop() {
            for &v in &g.adjacency[u] {
                if v < n_kept && comp[v] == usize::MAX {
                    comp[v] = n_comp;
                    stack.push(v);
                }
            }
        }
        n_comp += 1;
    }
    let mut frags = vec![Vec::new(); n_comp];
    for a in 0..n_kept {
        frags[comp[a]].push(a);
    }
    frags
}

/// 分子全体の正規 SMILES。多成分は各成分の正規形を辞書順に '.' 結合。
pub fn to_canonical_smiles(g: &MoleculeGraph) -> String {
    to_canonical_smiles_with_order(g).0
}

/// 正規 SMILES と、その出力順に並んだグラフ原子インデックス。
/// 立体署名付き複合キー (canonical_key) の位置決めに使う。
pub fn to_canonical_smiles_with_order(g: &MoleculeGraph) -> (String, Vec<usize>) {
    let mut parts: Vec<(String, Vec<usize>)> = get_fragments(g)
        .iter()
        .map(|frag| {
            let (s, local_order) =
                canonical_component_smiles(&build_cmol(g, frag), WriteMode::Strict);
            (s, local_order.into_iter().map(|li| frag[li]).collect())
        })
        .collect();
    parts.sort();
    let mut smi = String::new();
    let mut order = Vec::new();
    for (i, (s, o)) in parts.into_iter().enumerate() {
        if i > 0 {
            smi.push('.');
        }
        smi.push_str(&s);
        order.extend(o);
    }
    (smi, order)
}

/// 原子部分集合の正規 SMILES (部分構造キー用、RDKit `MolFragmentToSmiles` 相当)。
/// 部分集合が非連結なら成分ごとの正規形を辞書順に '.' 結合。
pub fn fragment_smiles(g: &MoleculeGraph, atoms: &[usize]) -> String {
    // subset 内での連結成分に分割
    let in_subset: std::collections::HashSet<usize> = atoms.iter().copied().collect();
    let n_kept = num_kept_atoms(g);
    let mut seen: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut parts = Vec::new();
    for &start in atoms {
        if seen.contains(&start) {
            continue;
        }
        let mut compo = vec![start];
        seen.insert(start);
        let mut stack = vec![start];
        while let Some(u) = stack.pop() {
            for &v in &g.adjacency[u] {
                if v < n_kept && in_subset.contains(&v) && !seen.contains(&v) {
                    seen.insert(v);
                    compo.push(v);
                    stack.push(v);
                }
            }
        }
        compo.sort_unstable();
        parts.push(canonical_component_smiles(&build_cmol(g, &compo), WriteMode::Lenient).0);
    }
    parts.sort();
    parts.join(".")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    fn canon(s: &str) -> String {
        to_canonical_smiles(&build_molecule_graph(s).expect("valid"))
    }

    #[test]
    fn invariance_basic() {
        assert_eq!(canon("CCO"), canon("OCC"));
        assert_eq!(canon("C(C)C"), canon("CCC"));
        assert_eq!(canon("c1ccccc1"), canon("C1=CC=CC=C1"));
        assert_eq!(canon("c1ccccc1C"), canon("Cc1ccccc1"));
        assert_eq!(canon("OC(=O)c1ccccc1"), canon("c1ccccc1C(O)=O"));
        assert_eq!(canon("[Na+].[Cl-]"), canon("[Cl-].[Na+]"));
        assert_eq!(canon("N1CC1"), canon("C1NC1"));
    }

    #[test]
    fn idempotency_basic() {
        for s in [
            "CCO",
            "c1ccccc1",
            "O=c1cccc[nH]1",
            "O=C1C=CC(=O)C=C1",
            "Cn1c(=S)oc2ccccc21",
            "[NH4+].[Cl-]",
            "[2H]OC",
            "F/C=C/F", // 立体は落ちるが冪等ではあること
            "CC(C)(C)c1ccccc1",
        ] {
            let c1 = canon(s);
            let c2 = canon(&c1);
            assert_eq!(c1, c2, "not idempotent for {s}: {c1} -> {c2}");
        }
    }

    #[test]
    fn distinct_molecules_differ() {
        assert_ne!(canon("CCO"), canon("CCC"));
        assert_ne!(canon("c1ccccc1"), canon("C1CCCCC1"));
        // N-メチルインドール vs インドール (nH の有無)
        assert_ne!(canon("Cn1ccc2ccccc21"), canon("c1ccc2[nH]ccc2c1"));
    }

    #[test]
    fn fragments() {
        let g = build_molecule_graph("CC(=O)[O-].[Na+]").expect("valid");
        let frags = get_fragments(&g);
        assert_eq!(frags.len(), 2);
        assert_eq!(frags[0], vec![0, 1, 2, 3]);
        assert_eq!(frags[1], vec![4]);
    }

    #[test]
    fn fragment_smiles_matches_standalone() {
        // ベンゼン環の部分集合キー = 独立ベンゼンの正規形 (Lenient の要件)
        let g = build_molecule_graph("c1ccccc1CCN").expect("valid");
        let ring: Vec<usize> = (0..6).collect();
        let standalone = build_molecule_graph("c1ccccc1").expect("valid");
        let key_frag = fragment_smiles(&g, &ring);
        let key_full = fragment_smiles(&standalone, &[0, 1, 2, 3, 4, 5]);
        assert_eq!(key_frag, key_full);

        // ピリジン核でも同様
        let g = build_molecule_graph("Clc1ccncc1").expect("valid");
        let ring: Vec<usize> = (1..7).collect();
        let standalone = build_molecule_graph("c1ccncc1").expect("valid");
        assert_eq!(
            fragment_smiles(&g, &ring),
            fragment_smiles(&standalone, &[0, 1, 2, 3, 4, 5])
        );
    }
}
