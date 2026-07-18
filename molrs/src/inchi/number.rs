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
struct NAtom {
    /// 元素順位: 炭素 = (0, "")、他は (1, symbol) でアルファベット順
    elem_key: (u8, String),
    degree: usize,
    n_h: u8,
    charge: i8,
    /// 成分内ローカル隣接 (ローカル idx, 結合次数クラス)
    nbrs: Vec<(usize, u8)>,
}

fn order_class(order: f64) -> u8 {
    if order == 3.0 {
        3
    } else if order == 2.0 {
        2
    } else {
        1
    }
}

/// 成分 (重原子ローカル集合) の番号付けデータを作る。
fn build_natoms(g: &MoleculeGraph, atoms: &[usize]) -> Vec<NAtom> {
    let mut local = std::collections::HashMap::new();
    for (li, &gi) in atoms.iter().enumerate() {
        local.insert(gi, li);
    }
    atoms
        .iter()
        .map(|&gi| {
            let a = &g.atoms[gi];
            let n_h = g.adjacency[gi]
                .iter()
                .filter(|&&x| g.atoms[x].symbol == "H")
                .count() as u8;
            let mut nbrs = Vec::new();
            for &nb in &g.adjacency[gi] {
                if let Some(&lj) = local.get(&nb) {
                    let ord = *g.bond_orders.get(&(gi.min(nb), gi.max(nb))).unwrap_or(&1.0);
                    nbrs.push((lj, order_class(ord)));
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
                degree: g.adjacency[gi]
                    .iter()
                    .filter(|&&x| local.contains_key(&x))
                    .count(),
                n_h,
                charge: a.formal_charge,
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
        let mut keys: Vec<(usize, Vec<(u8, usize)>)> = Vec::with_capacity(n);
        for (i, a) in atoms.iter().enumerate() {
            let mut nb: Vec<(u8, usize)> = a.nbrs.iter().map(|&(j, o)| (o, ranks[j])).collect();
            nb.sort_unstable();
            keys.push((ranks[i], nb));
        }
        let mut sorted: Vec<&(usize, Vec<(u8, usize)>)> = keys.iter().collect();
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

/// 番号付け → ソート済み辺リスト (min,max) の正準番号ペア。最小化の対象。
fn edge_signature(atoms: &[NAtom], numbering: &[usize]) -> Vec<(usize, usize)> {
    let mut edges = Vec::new();
    for (i, a) in atoms.iter().enumerate() {
        for &(j, _) in &a.nbrs {
            if i < j {
                let (u, v) = (numbering[i], numbering[j]);
                edges.push((u.min(v), u.max(v)));
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

/// 成分の正準番号 (ローカル idx → 0 始まり正準番号)。
fn number_component(atoms: &[NAtom]) -> Vec<usize> {
    // 初期色: (元素キー, 次数, H 数, 電荷)
    let mut init: Vec<(&(u8, String), usize, u8, i8)> = atoms
        .iter()
        .map(|a| (&a.elem_key, a.degree, a.n_h, a.charge))
        .collect();
    let mut sorted = init.clone();
    sorted.sort();
    sorted.dedup();
    let mut ranks: Vec<usize> = init
        .iter()
        .map(|k| sorted.binary_search(k).expect("init key"))
        .collect();
    init.clear();

    let mut budget = 5000usize;
    let mut best: Option<Candidate> = None;
    resolve(atoms, &ranks, &mut budget, &mut best);
    ranks = best.map(|(_, num)| num).unwrap_or_else(|| {
        // フォールバック: 単純精緻化順
        let mut r = ranks.clone();
        refine(atoms, &mut r);
        r
    });
    ranks
}

/// 分子全体の正準番号付け。返り値は成分ごとに
/// `canonical番号 (1 始まり) → 元の原子インデックス (0 始まり)` のベクタ。
/// RDKit AuxInfo `/N:` と同じ形式 (成分順は connected_components 準拠)。
pub fn canonical_numbering(g: &MoleculeGraph) -> Vec<Vec<usize>> {
    connected_components(g)
        .iter()
        .map(|atoms| {
            let natoms = build_natoms(g, atoms);
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
        // 対称な可動 H 群 (カルボキシルの 2 酸素等) を含む分子は
        // normalize (I4) で酸素を等価化してからでないと一致しない。
        // OCC(=O)O 等はそこで対応する。
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
