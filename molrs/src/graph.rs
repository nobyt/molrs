//! 分子グラフ構築 (S1.2)。
//!
//! Python 版 `molecule_analyzer.py` の `MoleculeGraph` を再現する。
//! RDKit の `MolFromSmiles` → `AddHs` 相当の処理:
//!
//! 1. 素の `[H]` 原子 (同位体・電荷・クラスなし、次数 1、単結合) は
//!    隣接重原子の H 数にマージして原子ノードから除去 (RDKit `RemoveHs` 相当)
//! 2. 有機サブセット原子の暗黙 H を原子価テーブルから計算
//! 3. 全ての H を明示的な原子ノードとして末尾に付加 (`AddHs` 相当、重原子順)
//! 4. 環閉じ結合は RDKit と同じく結合リスト末尾に環番号順で置き、向きは
//!    開き側に次数記号があれば (開き側, 閉じ側)、なければ (閉じ側, 開き側)
//!
//! `num_hs` は Python 版 (AddHs 後の `GetTotalNumHs()`) と同じく常に 0。
//! H の情報は明示 H 原子ノードと adjacency が持つ。
//!
//! 芳香族認識 (ケクレ入力の 1.5 化) は S1.3、`ring_atom_sets` (SSSR) は S1.4 で埋める。

use std::collections::HashMap;

use crate::aromaticity::{kekulize, perceive_aromaticity, AromAtom, AromBond};
use crate::elements::atomic_number;
use crate::rings::symmetrized_sssr;
use crate::smiles::{parse_smiles, BondKind, ParsedMolecule};
use crate::stereo::assign_stereochemistry;
use crate::ChemError;

/// Python 版 AtomInfo 相当。
#[derive(Debug, Clone, PartialEq)]
pub struct AtomInfo {
    pub idx: usize,
    pub symbol: String,
    pub atomic_num: u8,
    pub is_aromatic: bool,
    pub in_ring: bool,
    /// Python 互換: AddHs 後は常に 0 (H は原子ノードとして持つ)
    pub num_hs: u8,
    /// CIP コード 'R'/'S' (S1.7 で割り当て)
    pub chiral_tag: Option<char>,
    pub formal_charge: i8,
}

/// Python 版 BondInfo 相当。
#[derive(Debug, Clone, PartialEq)]
pub struct BondInfo {
    pub begin_idx: usize,
    pub end_idx: usize,
    /// 1.0 / 1.5 / 2.0 / 3.0
    pub bond_order: f64,
    /// 'E'/'Z' (S1.7 で割り当て)
    pub stereo: Option<char>,
}

/// Python 版 MoleculeGraph 相当。
#[derive(Debug, Clone)]
pub struct MoleculeGraph {
    pub atoms: Vec<AtomInfo>,
    pub bonds: Vec<BondInfo>,
    /// atom_idx → 隣接 atom_idx (結合追加順)
    pub adjacency: Vec<Vec<usize>>,
    /// (i, j) i<j → bond_order
    pub bond_orders: HashMap<(usize, usize), f64>,
    /// SSSR (S1.4 で構築)
    pub ring_atom_sets: Vec<Vec<usize>>,
    /// ケクレ化後の結合次数 (bonds と同順)。芳香族結合は 1.5 ではなく
    /// ケクレ形の 1.0/2.0 を持つ (2D 描画のケクレ表示用, D4)
    pub kekule_bond_orders: Vec<f64>,
    /// 元の SMILES パース結果 (S1.7 立体化学で使用)
    pub parsed: ParsedMolecule,
    /// パーサ原子 idx → グラフ原子 idx (マージされた H は None)
    pub parser_to_graph: Vec<Option<usize>>,
}

/// 2 原子間の結合次数。結合がなければ 0.0。(Python `get_bond_order` 相当)
pub fn get_bond_order(graph: &MoleculeGraph, i: usize, j: usize) -> f64 {
    let key = (i.min(j), i.max(j));
    graph.bond_orders.get(&key).copied().unwrap_or(0.0)
}

/// RDKit デフォルトの原子価リスト (有機サブセット元素のみ; 暗黙 H 計算対象)。
pub(crate) fn default_valences(symbol: &str) -> &'static [u8] {
    match symbol {
        "B" => &[3],
        "C" => &[4],
        "N" => &[3],
        "O" => &[2],
        "P" => &[3, 5],
        "S" => &[2, 4, 6],
        "F" | "Cl" | "Br" | "I" => &[1],
        _ => unreachable!("implicit H is only computed for organic-subset atoms"),
    }
}

/// 芳香族環内で二重結合を 1 本受け持つ元素 (暗黙 H 計算で原子価 +1)。
/// O/S/Se/Te は芳香環内で常に単結合 2 本なので +1 しない。
pub(crate) fn aromatic_takes_double_bond(symbol: &str) -> bool {
    matches!(symbol, "C" | "B" | "N" | "P" | "As")
}

/// SMILES から MoleculeGraph を構築する。(Python `build_molecule_graph` 相当)
pub fn build_molecule_graph(smiles: &str) -> Result<MoleculeGraph, ChemError> {
    let parsed = parse_smiles(smiles)?;
    let n = parsed.atoms.len();

    // ---- 1. 素の [H] のマージ判定 (RDKit RemoveHs 相当) ----
    let mut degree = vec![0usize; n];
    for b in &parsed.bonds {
        degree[b.a] += 1;
        degree[b.b] += 1;
    }
    let mut merged = vec![false; n];
    for (i, a) in parsed.atoms.iter().enumerate() {
        if a.symbol == "H"
            && a.isotope.is_none()
            && a.charge == 0
            && a.explicit_h.unwrap_or(0) == 0
            && a.atom_class.is_none()
            && degree[i] == 1
        {
            let bond = parsed
                .bonds
                .iter()
                .find(|b| b.a == i || b.b == i)
                .expect("degree 1");
            let other = if bond.a == i { bond.b } else { bond.a };
            if matches!(bond.kind, BondKind::Elided | BondKind::Single)
                && parsed.atoms[other].symbol != "H"
            {
                merged[i] = true;
            }
        }
    }

    // ---- 2. 重原子の再インデックス ----
    let mut parser_to_graph = vec![None; n];
    let mut graph_to_parser = Vec::new();
    for i in 0..n {
        if !merged[i] {
            parser_to_graph[i] = Some(graph_to_parser.len());
            graph_to_parser.push(i);
        }
    }
    let n_heavy = graph_to_parser.len();

    let mut merged_h = vec![0u8; n_heavy];
    for b in &parsed.bonds {
        for (h, other) in [(b.a, b.b), (b.b, b.a)] {
            if merged[h] {
                merged_h[parser_to_graph[other].expect("heavy")] += 1;
            }
        }
    }

    // ---- 3. 重原子間の結合 ----
    // RDKit 互換の結合順: 鎖結合を出現順に並べ、環閉じ結合は全体の末尾に
    // 環番号順 (同番号再利用時は閉じた順) で置く。向きは (閉じ側, 開き側)。
    struct KeptBond {
        a: usize,
        b: usize,
        kind: BondKind,
    }
    let mut kept_bonds: Vec<KeptBond> = Vec::new();
    let mut closure_bonds: Vec<(u16, KeptBond)> = Vec::new();
    for b in &parsed.bonds {
        if merged[b.a] || merged[b.b] {
            continue;
        }
        let x = parser_to_graph[b.a].expect("heavy");
        let y = parser_to_graph[b.b].expect("heavy");
        match b.ring_closure {
            Some(rc) => {
                // 開き側に次数記号があれば (開き側, 閉じ側)、なければ (閉じ側, 開き側)
                let (ba, bb) = if rc.opened_with_order { (x, y) } else { (y, x) };
                closure_bonds.push((
                    rc.num,
                    KeptBond {
                        a: ba,
                        b: bb,
                        kind: b.kind,
                    },
                ));
            }
            None => kept_bonds.push(KeptBond {
                a: x,
                b: y,
                kind: b.kind,
            }),
        }
    }
    closure_bonds.sort_by_key(|&(num, _)| num); // 安定ソート: 同番号は閉じた順のまま
    kept_bonds.extend(closure_bonds.into_iter().map(|(_, kb)| kb));

    // ---- 4. 橋 (bridge) 検出 → 環メンバーシップ ----
    let edges: Vec<(usize, usize)> = kept_bonds.iter().map(|b| (b.a, b.b)).collect();
    let is_bridge = find_bridges(n_heavy, &edges);
    let mut atom_in_ring = vec![false; n_heavy];
    for (ei, &(a, b)) in edges.iter().enumerate() {
        if !is_bridge[ei] {
            atom_in_ring[a] = true;
            atom_in_ring[b] = true;
        }
    }

    // ---- 5. H 数の確定 ----
    let mut num_hs_total = vec![0u8; n_heavy];
    for (gi, nh) in num_hs_total.iter_mut().enumerate() {
        let a = &parsed.atoms[graph_to_parser[gi]];
        *nh = if let Some(e) = a.explicit_h {
            // 角括弧原子: 明示 H のみ (暗黙 H なし)
            e + merged_h[gi]
        } else if a.symbol == "*" {
            merged_h[gi]
        } else {
            // 有機サブセット: 原子価テーブルから暗黙 H を計算。
            // 芳香族結合は 1 と数え、二重結合を受け持つ元素は +1 する
            // (RDKit はケクレ化後に計算するが、有効な SMILES では同じ結果になる)。
            let mut base = merged_h[gi] as usize;
            for kb in &kept_bonds {
                if kb.a == gi || kb.b == gi {
                    base += match kb.kind {
                        BondKind::Double => 2,
                        BondKind::Triple => 3,
                        BondKind::Quadruple => 4,
                        _ => 1,
                    };
                }
            }
            let v = if a.aromatic && aromatic_takes_double_bond(&a.symbol) {
                base + 1
            } else {
                base
            };
            let valences = default_valences(&a.symbol);
            match valences.iter().find(|&&t| t as usize >= v) {
                Some(&t) => merged_h[gi] + (t as usize - v) as u8,
                None if a.aromatic => merged_h[gi],
                None => {
                    return Err(ChemError::InvalidSmiles(format!(
                        "valence {v} for atom {} ({}) exceeds permitted in {smiles:?}",
                        gi, a.symbol
                    )));
                }
            }
        };
    }

    // ---- 6. ケクレ化と芳香族認識 (S1.3) ----
    let arom_atoms: Vec<AromAtom<'_>> = (0..n_heavy)
        .map(|gi| {
            let a = &parsed.atoms[graph_to_parser[gi]];
            AromAtom {
                symbol: &a.symbol,
                charge: a.charge,
                input_aromatic: a.aromatic,
                num_hs: num_hs_total[gi],
            }
        })
        .collect();
    let mut arom_bonds: Vec<AromBond> = kept_bonds
        .iter()
        .enumerate()
        .map(|(ei, kb)| {
            let both_aromatic = parsed.atoms[graph_to_parser[kb.a]].aromatic
                && parsed.atoms[graph_to_parser[kb.b]].aromatic;
            let (order, aromatic_candidate) = match kb.kind {
                BondKind::Single | BondKind::Up | BondKind::Down => (1.0, false),
                BondKind::Double => (2.0, false),
                BondKind::Triple => (3.0, false),
                BondKind::Quadruple => (4.0, false),
                BondKind::Aromatic => (1.0, true),
                // 省略結合: 両端が芳香族表記なら芳香族候補
                BondKind::Elided => (1.0, both_aromatic),
            };
            AromBond {
                a: kb.a,
                b: kb.b,
                order,
                aromatic_candidate,
                in_ring: !is_bridge[ei],
            }
        })
        .collect();
    kekulize(&arom_atoms, &mut arom_bonds)?;
    // RDKit と同じく対称化 SSSR を環情報として使う (芳香族認識にも同じものを渡す)
    let sssr_rings = symmetrized_sssr(n_heavy, &edges);
    let (atom_arom, bond_arom) = perceive_aromaticity(&arom_atoms, &arom_bonds, &sssr_rings);

    // ---- 7. AtomInfo / BondInfo の組み立て ----
    let mut atoms: Vec<AtomInfo> = Vec::with_capacity(n_heavy);
    for gi in 0..n_heavy {
        let a = &parsed.atoms[graph_to_parser[gi]];
        atoms.push(AtomInfo {
            idx: gi,
            symbol: a.symbol.clone(),
            atomic_num: atomic_number(&a.symbol).unwrap_or(0),
            is_aromatic: atom_arom[gi],
            in_ring: atom_in_ring[gi],
            num_hs: 0,
            chiral_tag: None,
            formal_charge: a.charge,
        });
    }

    let mut bonds: Vec<BondInfo> = Vec::new();
    let mut kekule_bond_orders: Vec<f64> = Vec::new();
    for (ei, ab) in arom_bonds.iter().enumerate() {
        bonds.push(BondInfo {
            begin_idx: ab.a,
            end_idx: ab.b,
            bond_order: if bond_arom[ei] { 1.5 } else { ab.order },
            stereo: None,
        });
        kekule_bond_orders.push(ab.order);
    }

    // ---- 7. 明示 H 原子の付加 (AddHs 相当: 重原子順に末尾へ) ----
    for (gi, &nh) in num_hs_total.iter().enumerate() {
        for _ in 0..nh {
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
                begin_idx: gi,
                end_idx: h_idx,
                bond_order: 1.0,
                stereo: None,
            });
            kekule_bond_orders.push(1.0);
        }
    }

    // ---- 8. adjacency / bond_orders ----
    let mut adjacency = vec![Vec::new(); atoms.len()];
    let mut bond_orders = HashMap::new();
    for b in &bonds {
        adjacency[b.begin_idx].push(b.end_idx);
        adjacency[b.end_idx].push(b.begin_idx);
        let key = (b.begin_idx.min(b.end_idx), b.begin_idx.max(b.end_idx));
        bond_orders.insert(key, b.bond_order);
    }

    let mut graph = MoleculeGraph {
        atoms,
        bonds,
        adjacency,
        bond_orders,
        ring_atom_sets: sssr_rings,
        kekule_bond_orders,
        parsed,
        parser_to_graph,
    };
    // ---- 9. 立体化学 (S1.7): CIP R/S と E/Z の割当て ----
    assign_stereochemistry(&mut graph);
    Ok(graph)
}

/// Tarjan の橋検出。edges[i] が橋なら結果の [i] が true。
fn find_bridges(n: usize, edges: &[(usize, usize)]) -> Vec<bool> {
    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n]; // (相手, edge_idx)
    for (ei, &(a, b)) in edges.iter().enumerate() {
        adj[a].push((b, ei));
        adj[b].push((a, ei));
    }
    let mut disc = vec![usize::MAX; n];
    let mut low = vec![usize::MAX; n];
    let mut is_bridge = vec![false; edges.len()];
    let mut timer = 0usize;

    // 反復 DFS (長鎖でのスタックオーバーフロー回避)
    for start in 0..n {
        if disc[start] != usize::MAX {
            continue;
        }
        // (node, parent_edge, 隣接リストの走査位置)
        let mut stack: Vec<(usize, usize, usize)> = vec![(start, usize::MAX, 0)];
        disc[start] = timer;
        low[start] = timer;
        timer += 1;
        while let Some(&mut (u, pe, ref mut it)) = stack.last_mut() {
            if *it < adj[u].len() {
                let (v, ei) = adj[u][*it];
                *it += 1;
                if ei == pe {
                    continue;
                }
                if disc[v] == usize::MAX {
                    disc[v] = timer;
                    low[v] = timer;
                    timer += 1;
                    stack.push((v, ei, 0));
                } else {
                    low[u] = low[u].min(disc[v]);
                }
            } else {
                stack.pop();
                if let Some(&mut (p, _, _)) = stack.last_mut() {
                    low[p] = low[p].min(low[u]);
                    if low[u] > disc[p] {
                        is_bridge[pe] = true;
                    }
                }
            }
        }
    }
    is_bridge
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g(s: &str) -> MoleculeGraph {
        build_molecule_graph(s).unwrap_or_else(|e| panic!("{s}: {e}"))
    }

    fn h_neighbors(graph: &MoleculeGraph, i: usize) -> usize {
        graph.adjacency[i]
            .iter()
            .filter(|&&j| graph.atoms[j].symbol == "H")
            .count()
    }

    #[test]
    fn ethanol() {
        let m = g("CCO");
        assert_eq!(m.atoms.len(), 9); // C,C,O + 6H
        assert_eq!(h_neighbors(&m, 0), 3);
        assert_eq!(h_neighbors(&m, 1), 2);
        assert_eq!(h_neighbors(&m, 2), 1);
        assert!(m.atoms.iter().all(|a| a.num_hs == 0));
        // H は末尾に重原子順で付加される
        assert_eq!(m.bonds[2].begin_idx, 0); // C0 の最初の H
        assert_eq!(m.bonds[2].end_idx, 3);
    }

    #[test]
    fn plain_h_merged() {
        // [H]C(=O)C: 素の [H] は C にマージされ、その後 AddHs で復元される
        let m = g("[H]C(=O)C");
        assert_eq!(m.atoms.len(), 7); // C,O,C + 4H
        assert_eq!(m.atoms[0].symbol, "C");
        assert_eq!(h_neighbors(&m, 0), 1); // アルデヒド C
        assert_eq!(h_neighbors(&m, 2), 3);
        assert_eq!(m.bonds[0].bond_order, 2.0); // C=O
        assert_eq!(m.parser_to_graph[0], None); // マージされた H
    }

    #[test]
    fn deuterium_kept() {
        // [2H] は原子として残る
        let m = g("[2H]OC");
        assert_eq!(m.atoms[0].symbol, "H");
        assert_eq!(h_neighbors(&m, 1), 1); // O の H は重水素のみ
        assert_eq!(h_neighbors(&m, 2), 3);
    }

    #[test]
    fn benzene_aromatic_input() {
        let m = g("c1ccccc1");
        assert_eq!(m.atoms.len(), 12);
        assert!(m.atoms[..6].iter().all(|a| a.is_aromatic && a.in_ring));
        assert!(m.bonds[..6].iter().all(|b| b.bond_order == 1.5));
        for i in 0..6 {
            assert_eq!(h_neighbors(&m, i), 1);
        }
        // 環閉じ結合は (閉じ側, 開き側)
        assert_eq!((m.bonds[5].begin_idx, m.bonds[5].end_idx), (5, 0));
    }

    #[test]
    fn aromatic_heteroatom_h_counts() {
        // チオフェン: 芳香族 S は H を持たない
        let m = g("c1ccsc1");
        assert_eq!(h_neighbors(&m, 3), 0);
        // フラン O
        let m = g("c1ccoc1");
        assert_eq!(h_neighbors(&m, 3), 0);
        // ピリジン n
        let m = g("c1ccncc1");
        assert_eq!(h_neighbors(&m, 3), 0);
        // ピロール [nH]
        let m = g("c1cc[nH]c1");
        assert_eq!(h_neighbors(&m, 3), 1);
        // ナフタレン縮合炭素 (次数 3)
        let m = g("c1ccc2ccccc2c1");
        let fused: Vec<usize> = (0..10)
            .filter(|&i| {
                m.adjacency[i]
                    .iter()
                    .filter(|&&j| m.atoms[j].symbol != "H")
                    .count()
                    == 3
            })
            .collect();
        assert_eq!(fused.len(), 2);
        for i in fused {
            assert_eq!(h_neighbors(&m, i), 0);
        }
    }

    #[test]
    fn biphenyl_bridge_is_single() {
        let m = g("c1ccc(cc1)c1ccccc1");
        // 架橋結合 (3, 6) は環外なので 1.0
        assert_eq!(get_bond_order(&m, 3, 6), 1.0);
        assert!(!m.bonds.iter().any(|b| b.bond_order == 1.5
            && ((b.begin_idx == 3 && b.end_idx == 6) || (b.begin_idx == 6 && b.end_idx == 3))));
    }

    #[test]
    fn charges_and_ions() {
        let m = g("[NH4+].[Cl-]");
        assert_eq!(m.atoms[0].formal_charge, 1);
        assert_eq!(m.atoms[1].formal_charge, -1);
        assert_eq!(h_neighbors(&m, 0), 4);
        assert_eq!(h_neighbors(&m, 1), 0);
    }

    #[test]
    fn hypervalent_sulfur() {
        let m = g("CS(=O)(=O)O"); // メタンスルホン酸: S は 6 価
        assert_eq!(h_neighbors(&m, 1), 0);
        assert_eq!(h_neighbors(&m, 4), 1); // OH
                                           // 中性 N の 4 価はエラー
        assert!(build_molecule_graph("CN(C)(C)C").is_err());
    }

    #[test]
    fn kekule_input_h_counts() {
        // ケクレ表記ベンゼン: 芳香族認識前でも H 数は一致する
        let m = g("C1=CC=CC=C1");
        for i in 0..6 {
            assert_eq!(h_neighbors(&m, i), 1);
        }
        assert!(m.atoms[..6].iter().all(|a| a.in_ring));
    }

    #[test]
    fn aromaticity_perception() {
        // ケクレ表記ベンゼン → 芳香族認識
        let m = g("C1=CC=CC=C1");
        assert!(m.atoms[..6].iter().all(|a| a.is_aromatic));
        assert!(m.bonds[..6].iter().all(|b| b.bond_order == 1.5));

        // キノン: 芳香族表記でも非芳香族にケクレ化される
        let m = g("O=c1ccc(=O)cc1");
        assert!(m.atoms.iter().all(|a| !a.is_aromatic));
        assert!(!m.bonds.iter().any(|b| b.bond_order == 1.5));

        // 2-ピリドン: 芳香族のまま (nH 2e + c(=O) 0e + 4c 4e = 6)
        let m = g("O=c1cccc[nH]1");
        assert!(m.atoms[1..7].iter().all(|a| a.is_aromatic));

        // N-メチルベンゾオキサゾール-2-チオン (Phase 853 系)
        let m = g("Cn1c(=S)oc2ccccc21");
        assert!(m.atoms[1].is_aromatic); // n
        assert!(m.atoms[2].is_aromatic); // c(=S)
        assert!(!m.atoms[0].is_aromatic && !m.atoms[3].is_aromatic); // CH3, S

        // シクロブタジエン: 4n 電子 → 非芳香族
        let m = g("C1=CC=C1");
        assert!(m.atoms.iter().all(|a| !a.is_aromatic));

        // フルベン: 環外 C=C の炭素は候補外 → 非芳香族
        let m = g("C=C1C=CC=C1");
        assert!(m.atoms.iter().all(|a| !a.is_aromatic));

        // アズレン: 単環では 4n だがユニオン (10 原子) で芳香族
        let m = g("c1ccc2cccc2cc1");
        assert!(m.atoms[..10].iter().all(|a| a.is_aromatic));

        // 環に入っていない芳香族原子はエラー
        assert!(build_molecule_graph("cc").is_err());
        // 5 員環全炭素芳香族はケクレ化不能
        assert!(build_molecule_graph("c1cccc1").is_err());
    }

    #[test]
    fn ring_membership() {
        let m = g("C1CC1CC"); // シクロプロパン + エチル
        assert!(m.atoms[0].in_ring && m.atoms[1].in_ring && m.atoms[2].in_ring);
        assert!(!m.atoms[3].in_ring && !m.atoms[4].in_ring);
        // スピロ環: 全原子が環内
        let m = g("C1CC12CC2");
        assert!(m.atoms[..5].iter().all(|a| a.in_ring));
    }
}
