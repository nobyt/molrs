//! InChI 正準番号付け (I3)。
//!
//! InChI は分子を成分に分け、各成分内の重原子に 1 始まりの正準番号を付ける。
//! 番号は初期色 (元素・次数・H 数・電荷など) を Morgan 系の反復精緻化で
//! 細分し、同値類が残る場合は正準最小化 (各候補を試して最小の接続表を選ぶ)
//! で確定する。
//!
//! 検証: RDKit AuxInfo の `/N:` フィールド (公式の正準番号) と一致させる。
//!
//! 本ファイルは I3 で番号付け本体を実装する。まず成分分解
//! ([`connected_components`]) を提供する (formula.rs が使用)。

use crate::graph::MoleculeGraph;

/// 重原子 (非 H) のインデックス一覧。番号付け本体 (I3) で使用する。
#[allow(dead_code)]
pub(crate) fn heavy_atoms(g: &MoleculeGraph) -> Vec<usize> {
    (0..g.atoms.len())
        .filter(|&i| g.atoms[i].symbol != "H")
        .collect()
}

/// 重原子の連結成分。各成分は原子インデックスの Vec。
/// 成分の順序は最小原子インデックスの昇順 (安定・決定的)。
pub(crate) fn connected_components(g: &MoleculeGraph) -> Vec<Vec<usize>> {
    let n = g.atoms.len();
    let is_heavy = |i: usize| g.atoms[i].symbol != "H";
    let mut comp = vec![usize::MAX; n];
    let mut n_comp = 0;
    for start in 0..n {
        if !is_heavy(start) || comp[start] != usize::MAX {
            continue;
        }
        let mut stack = vec![start];
        comp[start] = n_comp;
        while let Some(v) = stack.pop() {
            for &nb in &g.adjacency[v] {
                if is_heavy(nb) && comp[nb] == usize::MAX {
                    comp[nb] = n_comp;
                    stack.push(nb);
                }
            }
        }
        n_comp += 1;
    }
    let mut comps: Vec<Vec<usize>> = vec![Vec::new(); n_comp];
    for i in 0..n {
        if is_heavy(i) && comp[i] != usize::MAX {
            comps[comp[i]].push(i);
        }
    }
    // 各成分内は昇順、成分自体は最小 idx 昇順 (push 順で既に満たす)
    comps
}

/// 成分内の 1 原子の番号付け用データ。
///
/// InChI の正準番号付けは**単純グラフ** (結合次数なし) で行う — c 層に
/// 結合次数が現れないことと対応する (C#CCC#N が恒等番号になるのは
/// このため)。H 数はトポロジー精緻化の後の第 2 段でのみ効く (プロペンの
/// CH2 が CH3 より先)。可動 H 群 (カルボキシル等) のメンバーは固定 H を
/// 持たない扱いで等価化する (マロン酸で実証)。
struct NAtom {
    /// 元素順位: 炭素 = (0, "")、他は (1, symbol) でアルファベット順
    elem_key: (u8, String),
    degree: usize,
    /// 固定 H 数 (可動 H 群メンバーは 0)
    n_h: u8,
    /// 電荷 (可動群で相殺されるものは 0)
    charge: i8,
    /// 可動 H 群 (t-group) のメンバーか
    in_tgroup: bool,
    /// 成分内ローカル隣接 (ローカル idx)
    nbrs: Vec<usize>,
}

/// 可動 H 群 (標準互変異性の 1,3-シフト) のメンバー原子を検出する。
///
/// v1 の規則: ある中心原子に結合した末端 (重原子次数 1) の O/S/Se/Te/N を
/// 端点集合とし、端点が 2 個以上・二重結合端点が 1 個以上・(H を持つ端点
/// または負電荷端点が 1 個以上) のとき、その端点集合を可動群とする。
/// カルボキシル/スルホン酸/アミド/カルボキシラート/ニトロを覆う。
/// 環内 (イミダゾール等) の長距離互変異性は未対応 (v2)。
pub(crate) fn tautomer_group_members(g: &MoleculeGraph) -> Vec<bool> {
    let n = g.atoms.len();
    let mut member = vec![false; n];
    let heavy_deg = |i: usize| {
        g.adjacency[i]
            .iter()
            .filter(|&&x| g.atoms[x].symbol != "H")
            .count()
    };
    let n_h = |i: usize| {
        g.adjacency[i]
            .iter()
            .filter(|&&x| g.atoms[x].symbol == "H")
            .count()
    };
    for center in 0..n {
        if g.atoms[center].symbol == "H" {
            continue;
        }
        let mut endpoints = Vec::new();
        for &nb in &g.adjacency[center] {
            let sym = g.atoms[nb].symbol.as_str();
            if matches!(sym, "O" | "S" | "Se" | "Te" | "N") && heavy_deg(nb) == 1 {
                endpoints.push(nb);
            }
        }
        if endpoints.len() < 2 {
            continue;
        }
        let n_double = endpoints
            .iter()
            .filter(|&&e| {
                let key = (center.min(e), center.max(e));
                g.bond_orders.get(&key).copied().unwrap_or(1.0) >= 2.0
            })
            .count();
        let n_hydrogens: usize = endpoints.iter().map(|&e| n_h(e)).sum();
        let n_neg = endpoints
            .iter()
            .filter(|&&e| g.atoms[e].formal_charge < 0)
            .count();
        if n_double >= 1 && (n_hydrogens >= 1 || n_neg >= 1) {
            for &e in &endpoints {
                member[e] = true;
            }
        }
    }
    member
}

/// 成分 (重原子ローカル集合) の番号付けデータを作る。
fn build_natoms(g: &MoleculeGraph, atoms: &[usize], tgroup: &[bool]) -> Vec<NAtom> {
    let mut local = std::collections::HashMap::new();
    for (li, &gi) in atoms.iter().enumerate() {
        local.insert(gi, li);
    }
    atoms
        .iter()
        .map(|&gi| {
            let a = &g.atoms[gi];
            let in_tgroup = tgroup[gi];
            let n_h = if in_tgroup {
                0 // 可動 H はメンバー間で等価化 (群レベルで h 層に出る)
            } else {
                g.adjacency[gi]
                    .iter()
                    .filter(|&&x| g.atoms[x].symbol == "H")
                    .count() as u8
            };
            let charge = if in_tgroup && a.formal_charge < 0 {
                0 // 群内の負電荷はプロトン除去 (p 層) に正規化される
            } else {
                a.formal_charge
            };
            let mut nbrs = Vec::new();
            for &nb in &g.adjacency[gi] {
                if let Some(&lj) = local.get(&nb) {
                    nbrs.push(lj);
                }
            }
            nbrs.sort_unstable();
            let elem_key = if a.symbol == "C" {
                (0u8, String::new())
            } else {
                (1u8, a.symbol.clone())
            };
            NAtom {
                elem_key,
                degree: nbrs.len(),
                n_h,
                charge,
                in_tgroup,
                nbrs,
            }
        })
        .collect()
}

/// 隣接ランクによる反復精緻化。ranks はクラス id (0 始まり)。クラス数を返す。
fn refine(atoms: &[NAtom], ranks: &mut [usize]) -> usize {
    let n = atoms.len();
    let mut n_classes = ranks.iter().max().map_or(0, |m| m + 1);
    loop {
        let mut keys: Vec<(usize, Vec<usize>)> = Vec::with_capacity(n);
        for (i, a) in atoms.iter().enumerate() {
            let mut nb: Vec<usize> = a.nbrs.iter().map(|&j| ranks[j]).collect();
            nb.sort_unstable();
            keys.push((ranks[i], nb));
        }
        let mut sorted: Vec<&(usize, Vec<usize>)> = keys.iter().collect();
        sorted.sort();
        sorted.dedup();
        let new_n = sorted.len();
        for (i, r) in ranks.iter_mut().enumerate() {
            *r = sorted.binary_search(&&keys[i]).expect("key exists");
        }
        if new_n == n_classes {
            return n_classes;
        }
        n_classes = new_n;
    }
}

/// (辺シグネチャ, 番号付け) — 最小化の候補。
type Candidate = (Vec<(usize, usize)>, Vec<usize>);

/// 番号付け → InChI 接続表の比較キー。
/// InChI の c 層は「原子 k (2..n) ごとに、より小さい番号の隣接」を並べる。
/// 辺 (j,k) j<k は k のグループに現れるため、比較順は (大きい端点, 小さい端点)。
fn edge_signature(atoms: &[NAtom], numbering: &[usize]) -> Vec<(usize, usize)> {
    let mut edges = Vec::new();
    for (i, a) in atoms.iter().enumerate() {
        for &j in &a.nbrs {
            if i < j {
                let (u, v) = (numbering[i], numbering[j]);
                edges.push((u.max(v), u.min(v)));
            }
        }
    }
    edges.sort_unstable();
    edges
}

/// 精緻化されたランクから正準番号を確定する。同値類が残る場合は
/// 各メンバーで分岐し、辺シグネチャが辞書順最小となる番号付けを採る。
fn resolve(atoms: &[NAtom], ranks: &[usize], budget: &mut usize, best: &mut Option<Candidate>) {
    let n = atoms.len();
    let mut ranks = ranks.to_vec();
    let n_classes = refine(atoms, &mut ranks);

    if n_classes == n {
        // 全原子が一意ランク → ランク昇順に 0..n を割り当て
        let numbering = ranks.clone(); // rank i (0..n) = canonical番号
        let sig = edge_signature(atoms, &numbering);
        if best.as_ref().map(|(s, _)| &sig < s).unwrap_or(true) {
            *best = Some((sig, numbering));
        }
        return;
    }
    if *budget == 0 {
        return;
    }

    let mut class_size = vec![0usize; n_classes];
    for &r in &ranks {
        class_size[r] += 1;
    }
    let target = (0..n_classes)
        .find(|&c| class_size[c] > 1)
        .expect("tied class");
    let members: Vec<usize> = (0..n).filter(|&i| ranks[i] == target).collect();
    for &m in &members {
        if *budget == 0 {
            break;
        }
        *budget -= 1;
        // m を target クラスの先頭に固定 (rank を 1 つ下げる)
        let mut branched = ranks.clone();
        branched[m] = target; // 保持
        for (i, r) in branched.iter_mut().enumerate() {
            if i != m && ranks[i] >= target {
                *r += 1;
            }
        }
        resolve(atoms, &branched, budget, best);
    }
}

/// 色キー列からランク (0 始まりクラス id) を作る。
fn ranks_from_keys<K: Ord>(keys: &[K]) -> Vec<usize> {
    let mut sorted: Vec<&K> = keys.iter().collect();
    sorted.sort();
    sorted.dedup();
    keys.iter()
        .map(|k| sorted.binary_search(&k).expect("key"))
        .collect()
}

/// 成分の正準番号 (ローカル idx → 0 始まり正準番号)。
///
/// InChI の多段正準化 (実測から再構成):
/// 1. トポロジーのみ (元素 + 次数、結合次数なし) で精緻化
/// 2. 固定 H 数 (+t-group フラグ) を加えて精緻化 (H 昇順)
/// 3. 電荷を加えて精緻化
/// 4. 残る同値類は分岐し、InChI 接続表 (edge_signature) 最小の番号を採用
fn number_component(atoms: &[NAtom]) -> Vec<usize> {
    // 段 1: (元素, 次数)
    let keys1: Vec<(&(u8, String), usize)> =
        atoms.iter().map(|a| (&a.elem_key, a.degree)).collect();
    let mut ranks = ranks_from_keys(&keys1);
    refine(atoms, &mut ranks);
    // 段 2: + (t-group, 固定 H)
    let keys2: Vec<(usize, bool, u8)> = atoms
        .iter()
        .enumerate()
        .map(|(i, a)| (ranks[i], a.in_tgroup, a.n_h))
        .collect();
    ranks = ranks_from_keys(&keys2);
    refine(atoms, &mut ranks);
    // 段 3: + 電荷
    let keys3: Vec<(usize, i8)> = atoms
        .iter()
        .enumerate()
        .map(|(i, a)| (ranks[i], a.charge))
        .collect();
    ranks = ranks_from_keys(&keys3);
    refine(atoms, &mut ranks);

    let mut budget = 5000usize;
    let mut best: Option<Candidate> = None;
    resolve(atoms, &ranks, &mut budget, &mut best);
    best.map(|(_, num)| num).unwrap_or_else(|| {
        let mut r = ranks.clone();
        refine(atoms, &mut r);
        r
    })
}

/// 分子全体の正準番号付け。返り値は成分ごとに
/// `canonical番号 (1 始まり) → 元の原子インデックス (0 始まり)` のベクタ。
/// RDKit AuxInfo `/N:` と同じ形式 (成分順は connected_components 準拠)。
pub fn canonical_numbering(g: &MoleculeGraph) -> Vec<Vec<usize>> {
    let tgroup = tautomer_group_members(g);
    connected_components(g)
        .iter()
        .map(|atoms| {
            let natoms = build_natoms(g, atoms, &tgroup);
            let numbering = number_component(&natoms); // local idx → 0-based canon
                                                       // canon番号 → 元の原子 idx
            let mut inv = vec![0usize; atoms.len()];
            for (li, &cn) in numbering.iter().enumerate() {
                inv[cn] = atoms[li];
            }
            inv
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    /// canonical_numbering を AuxInfo /N: 形式 (1 始まり元 idx) に変換。
    fn numbering_1based(smiles: &str) -> Vec<Vec<usize>> {
        let g = build_molecule_graph(smiles).unwrap();
        canonical_numbering(&g)
            .iter()
            .map(|comp| comp.iter().map(|&i| i + 1).collect())
            .collect()
    }

    #[test]
    fn linear_and_simple() {
        // RDKit AuxInfo /N: と一致すべき既知ケース
        assert_eq!(numbering_1based("CCO"), vec![vec![1, 2, 3]]);
        assert_eq!(numbering_1based("CCN"), vec![vec![1, 2, 3]]);
        assert_eq!(numbering_1based("CCCC"), vec![vec![1, 4, 2, 3]]);
        // 単純グラフ (結合次数無視): C#CCC#N は恒等番号
        assert_eq!(numbering_1based("C#CCC#N"), vec![vec![1, 2, 3, 4, 5]]);
        // H はトポロジー後の第 2 段 (昇順): プロペンは CH2 が先
        assert_eq!(numbering_1based("CC=C"), vec![vec![3, 1, 2]]);
        assert_eq!(numbering_1based("CCC=C"), vec![vec![4, 1, 3, 2]]);
        assert_eq!(
            numbering_1based("C(/C=C/C)CC"),
            vec![vec![4, 6, 3, 5, 2, 1]]
        );
    }

    #[test]
    fn mobile_h_symmetrization() {
        // カルボキシル O は可動 H 群で等価化 (マロン酸で OH/=O が交互になる)
        assert_eq!(numbering_1based("CC(=O)O"), vec![vec![1, 2, 3, 4]]);
        assert_eq!(numbering_1based("OCC(=O)O"), vec![vec![2, 3, 1, 4, 5]]);
        assert_eq!(
            numbering_1based("OC(=O)CC(=O)O"),
            vec![vec![4, 2, 5, 1, 3, 6, 7]]
        );
        // スルホン酸: 3 つの O が同一群
        assert_eq!(numbering_1based("CS(=O)(=O)O"), vec![vec![1, 3, 4, 5, 2]]);
    }

    #[test]
    fn components_split() {
        let g = build_molecule_graph("CCO.O").unwrap();
        let comps = connected_components(&g);
        assert_eq!(comps.len(), 2);
        // 重原子のみ (H は含まない)
        for c in &comps {
            for &a in c {
                assert_ne!(g.atoms[a].symbol, "H");
            }
        }
        assert_eq!(comps[0].len(), 3); // C C O
        assert_eq!(comps[1].len(), 1); // O
    }

    #[test]
    fn heavy_atoms_excludes_h() {
        let g = build_molecule_graph("C").unwrap();
        assert_eq!(heavy_atoms(&g).len(), 1);
    }
}
