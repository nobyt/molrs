//! 立体化学 (S1.7): RDKit レガシー CIP 実装の忠実な移植。
//!
//! Python 版は `Chem.AssignStereochemistry(cleanIt=True, force=True)` の
//! `_CIPCode` (R/S) と `GetStereo()` (E/Z) を読む。これは新 CIP ラベラではなく
//! レガシー実装 (Chirality.cpp, Release_2023_09) なので、それをそのまま移植する:
//!
//! 1. `buildCIPInvariants` + `iterateCIPRanks`: 原子番号・同位体質量差から
//!    初期不変量を作り、隣接ランクを結合次数の 2 倍回複製した列で反復精緻化
//! 2. `AdjustAtomChiralityFlags` (SmilesParseOps.cpp): SMILES 出現順の
//!    近傍リストと分子の結合リスト順のパリティ差でキラルタグを補正
//! 3. `assignAtomChiralCodes`: 正当な四面体中心か検査し、CIP ランク昇順への
//!    置換パリティ (+ 3 隣接 1H の補正) で R/S を決定
//! 4. `assignBondStereoCodes`: 方向結合を二重結合起点向きに正規化し、
//!    両側の最高ランク置換基の向きが同じなら Z、違えば E

use std::collections::HashMap;

use crate::graph::MoleculeGraph;
use crate::smiles::{BondKind, Chirality};

/// 主同位体の質量数 (RDKit getMostCommonIsotope 相当; 使用元素のみ)。
fn most_common_isotope(atomic_num: u8) -> i32 {
    match atomic_num {
        1 => 1,
        3 => 7,
        5 => 11,
        6 => 12,
        7 => 14,
        8 => 16,
        9 => 19,
        11 => 23,
        12 => 24,
        13 => 27,
        14 => 28,
        15 => 31,
        16 => 32,
        17 => 35,
        19 => 39,
        20 => 40,
        26 => 56,
        30 => 64,
        32 => 74,
        33 => 75,
        34 => 80,
        35 => 79,
        50 => 120,
        51 => 121,
        52 => 130,
        53 => 127,
        80 => 202,
        82 => 208,
        83 => 209,
        _ => 0,
    }
}

/// getTwiceBondType 相当 (結合次数の 2 倍; 芳香族 = 3)。
fn twice_bond_type(order: f64) -> u32 {
    if order == 1.5 {
        3
    } else if order == 2.0 {
        4
    } else if order == 3.0 {
        6
    } else if order == 4.0 {
        8
    } else {
        2
    }
}

/// Rankers::rankVect 相当 (稠密ランク、同値は同ランク)。
fn rank_vect<T: Ord>(vals: &[T], ranks: &mut [usize]) {
    let n = vals.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| vals[a].cmp(&vals[b]));
    let mut curr = 0usize;
    let mut last = idx[0];
    for &i in &idx {
        if vals[i] == vals[last] {
            ranks[i] = curr;
        } else {
            curr += 1;
            ranks[i] = curr;
            last = i;
        }
    }
}

/// 2 つの同要素列の置換パリティ (countSwapsToInterconvert 相当、偶奇のみ)。
fn swap_parity(reference: &[usize], probe: &[usize]) -> usize {
    debug_assert_eq!(reference.len(), probe.len());
    let mut work: Vec<usize> = probe.to_vec();
    let mut swaps = 0usize;
    for i in 0..reference.len() {
        if work[i] != reference[i] {
            let j = (i + 1..work.len())
                .find(|&j| work[j] == reference[i])
                .expect("same elements");
            work.swap(i, j);
            swaps += 1;
        }
    }
    swaps
}

/// 立体化学コンテキスト (kept 原子のみ)。
struct Ctx {
    n: usize,
    atomic_num: Vec<u8>,
    isotope: Vec<Option<u16>>,
    n_h: Vec<u8>,
    /// atom → グラフ結合 id のリスト (graph.bonds 順 = RDKit の結合リスト順)
    atom_bonds: Vec<Vec<usize>>,
    /// グラフ結合 id → (begin, end, order)
    bonds: Vec<(usize, usize, f64)>,
    /// 結合が属する最小環サイズ (環外は None)
    min_bond_ring: Vec<Option<usize>>,
    /// 原子の環結合数 (橋頭判定用)
    ring_bond_count: Vec<usize>,
    /// 原子が 3 員環に属するか
    in_ring3: Vec<bool>,
}

fn build_ctx(g: &MoleculeGraph) -> Ctx {
    let n = g.parser_to_graph.iter().flatten().count();
    let bonds: Vec<(usize, usize, f64)> = g
        .bonds
        .iter()
        .filter(|b| b.begin_idx < n && b.end_idx < n)
        .map(|b| (b.begin_idx, b.end_idx, b.bond_order))
        .collect();
    let mut atom_bonds = vec![Vec::new(); n];
    for (ei, &(a, b, _)) in bonds.iter().enumerate() {
        atom_bonds[a].push(ei);
        atom_bonds[b].push(ei);
    }
    // 環メンバーシップ (対称化 SSSR から)
    let mut min_bond_ring = vec![None; bonds.len()];
    let mut ring_bond_count = vec![0usize; n];
    let mut in_ring3 = vec![false; n];
    let bond_between = |u: usize, v: usize| -> Option<usize> {
        bonds
            .iter()
            .position(|&(a, b, _)| (a == u && b == v) || (a == v && b == u))
    };
    let mut counted = vec![false; bonds.len()];
    for ring in &g.ring_atom_sets {
        for k in 0..ring.len() {
            let (u, v) = (ring[k], ring[(k + 1) % ring.len()]);
            if ring.len() == 3 {
                in_ring3[u] = true;
            }
            if let Some(bi) = bond_between(u, v) {
                let cur = min_bond_ring[bi];
                if cur.is_none() || ring.len() < cur.unwrap() {
                    min_bond_ring[bi] = Some(ring.len());
                }
                if !counted[bi] {
                    counted[bi] = true;
                    ring_bond_count[u] += 1;
                    ring_bond_count[v] += 1;
                }
            }
        }
    }
    // graph idx → parser idx (同位体)
    let mut graph_to_parser = vec![usize::MAX; n];
    for (pi, slot) in g.parser_to_graph.iter().enumerate() {
        if let Some(gi) = slot {
            graph_to_parser[*gi] = pi;
        }
    }
    Ctx {
        n,
        atomic_num: (0..n).map(|i| g.atoms[i].atomic_num).collect(),
        isotope: (0..n)
            .map(|i| g.parsed.atoms[graph_to_parser[i]].isotope)
            .collect(),
        n_h: (0..n)
            .map(|i| g.adjacency[i].iter().filter(|&&x| x >= n).count() as u8)
            .collect(),
        atom_bonds,
        bonds,
        min_bond_ring,
        ring_bond_count,
        in_ring3,
    }
}

/// buildCIPInvariants + iterateCIPRanks 相当。
fn assign_cip_ranks(ctx: &Ctx) -> Vec<usize> {
    let n = ctx.n;
    if n == 0 {
        return Vec::new();
    }
    // 初期不変量: (原子番号 % 128) << 10 | 質量フィールド、さらに << 10 (map 番号 0)
    let invars: Vec<u64> = (0..n)
        .map(|i| {
            let num = (ctx.atomic_num[i] as u64) % 128;
            let mut mass: i64 = 0;
            if let Some(iso) = ctx.isotope[i] {
                mass = iso as i64 - most_common_isotope(ctx.atomic_num[i]) as i64;
                if mass >= 0 {
                    mass += 1;
                }
            }
            mass += 512;
            let mass = if mass < 0 { 0 } else { (mass % 1024) as u64 };
            // 最後の << 10 はアトムマップ番号フィールド (常に 0)
            ((num << 10) | mass) << 10
        })
        .collect();

    let mut ranks = vec![0usize; n];
    rank_vect(&invars, &mut ranks);

    // cipEntries[i] = [原子番号, rank] で開始
    let mut entries: Vec<Vec<i64>> = (0..n)
        .map(|i| vec![ctx.atomic_num[i] as i64, ranks[i] as i64])
        .collect();

    let max_its = n / 2 + 1;
    let mut num_its = 0usize;
    let mut last_num_ranks: i64 = -1;
    let mut num_ranks = ranks.iter().max().unwrap() + 1;
    let mut counts = vec![0u32; n];

    while num_ranks < n
        && num_its < max_its
        && (last_num_ranks < 0 || (last_num_ranks as usize) < num_ranks)
    {
        let mut longest = 0usize;
        for (i, entry) in entries.iter_mut().enumerate() {
            let mut nbr_idxs: Vec<usize> = Vec::with_capacity(ctx.atom_bonds[i].len());
            for &ei in &ctx.atom_bonds[i] {
                let (a, b, order) = ctx.bonds[ei];
                let nbr = if a == i { b } else { a };
                nbr_idxs.push(nbr);
                // キラルリン化合物の特例 (二重結合先が P で次数 3/4 → 重み 1)
                let weight = if order == 2.0
                    && ctx.atomic_num[nbr] == 15
                    && (ctx.atom_bonds[nbr].len() == 3 || ctx.atom_bonds[nbr].len() == 4)
                {
                    1
                } else {
                    twice_bond_type(order)
                };
                counts[nbr] += weight;
            }
            // ランク降順に、重み回数だけ (rank+1) を複製して追加
            nbr_idxs.sort_by(|&x, &y| ranks[y].cmp(&ranks[x]));
            for nbr in nbr_idxs {
                let c = counts[nbr] as usize;
                for _ in 0..c {
                    entry.push(ranks[nbr] as i64 + 1);
                }
                counts[nbr] = 0;
            }
            // H 1 つにつき 0 を追加
            for _ in 0..ctx.n_h[i] {
                entry.push(0);
            }
            longest = longest.max(entry.len());
        }
        // -1 でパディングして同じ長さに
        for e in entries.iter_mut() {
            e.resize(longest, -1);
        }

        last_num_ranks = num_ranks as i64;
        rank_vect(&entries, &mut ranks);
        num_ranks = ranks.iter().max().unwrap() + 1;

        if last_num_ranks as usize != num_ranks {
            for (i, e) in entries.iter_mut().enumerate() {
                e[num_its + 1] = ranks[i] as i64;
                e.truncate(num_its + 2);
            }
        }
        num_its += 1;
    }
    ranks
}

/// SMILES 出現順のキラルタグを分子の結合リスト順基準に補正する
/// (AdjustAtomChiralityFlags 相当)。返り値: 補正済みタグ (graph atom idx →)。
fn adjusted_chiral_tags(g: &MoleculeGraph, ctx: &Ctx) -> HashMap<usize, Chirality> {
    let n = ctx.n;
    // parser bond idx → graph bond idx
    let mut bond_map: HashMap<(usize, usize), usize> = HashMap::new();
    for (ei, &(a, b, _)) in ctx.bonds.iter().enumerate() {
        bond_map.insert((a.min(b), a.max(b)), ei);
    }

    let mut tags = HashMap::new();
    for (pi, pa) in g.parsed.atoms.iter().enumerate() {
        let Some(chir) = pa.chirality else { continue };
        let Some(gi) = g.parser_to_graph[pi] else {
            continue;
        };
        if gi >= n {
            continue;
        }
        // SMILES 出現順のグラフ結合 id 列
        let mut smiles_order: Vec<usize> = Vec::new();
        let mut ok = true;
        for &pb in &g.parsed.neighbor_order[pi] {
            let b = &g.parsed.bonds[pb];
            let (pa_i, pb_i) = (b.a, b.b);
            let (Some(x), Some(y)) = (g.parser_to_graph[pa_i], g.parser_to_graph[pb_i]) else {
                ok = false; // マージされた H が近傍 (コーパスには出現しない)
                break;
            };
            match bond_map.get(&(x.min(y), x.max(y))) {
                Some(&ei) => smiles_order.push(ei),
                None => {
                    ok = false;
                    break;
                }
            }
        }
        if !ok {
            continue;
        }
        let mut n_swaps = swap_parity(&ctx.atom_bonds[gi], &smiles_order);

        // chiralAtomNeedsTagInversion 相当
        let degree = smiles_order.len();
        let explicit_h = pa.explicit_h.unwrap_or(0);
        let is_first = !g
            .parsed
            .bonds
            .iter()
            .any(|b| b.b == pi && b.ring_closure.is_none());
        // この原子に環閉じ数字が付いた数 (開き側・閉じ側どちらも)
        let n_closures = g.parsed.neighbor_order[pi]
            .iter()
            .filter(|&&pb| g.parsed.bonds[pb].ring_closure.is_some())
            .count();
        let has_fourth_valence = explicit_h == 1;
        let unsaturated = g.parsed.neighbor_order[pi].iter().any(|&pb| {
            matches!(
                g.parsed.bonds[pb].kind,
                BondKind::Double | BondKind::Triple | BondKind::Quadruple | BondKind::Aromatic
            )
        });
        if degree == 3
            && ((is_first && explicit_h == 1)
                || (!has_fourth_valence && n_closures == 1 && !unsaturated))
        {
            n_swaps += 1;
        }

        let tag = if n_swaps % 2 == 1 {
            match chir {
                Chirality::Anticlockwise => Chirality::Clockwise,
                Chirality::Clockwise => Chirality::Anticlockwise,
            }
        } else {
            chir
        };
        tags.insert(gi, tag);
    }
    tags
}

/// assignAtomChiralCodes 相当。atoms[gi].chiral_tag に 'R'/'S' を設定。
fn assign_atom_chiral_codes(g: &mut MoleculeGraph, ctx: &Ctx, ranks: &[usize]) {
    let tags = adjusted_chiral_tags(g, ctx);
    for (&gi, &tag) in &tags {
        let nz_degree = ctx.atom_bonds[gi].len();
        let tnz_degree = nz_degree + ctx.n_h[gi] as usize;
        let anum = ctx.atomic_num[gi];

        // 正当な四面体中心の判定 (isAtomPotentialChiralCenter 相当)
        let legal =
            if !(3..=4).contains(&tnz_degree) || (nz_degree < 3 && !(anum == 15 || anum == 33)) {
                false
            } else if nz_degree == 3 {
                if ctx.n_h[gi] == 1 {
                    true // protium 隣接は MolFromSmiles 後には存在しない
                } else {
                    match anum {
                        7 => ctx.in_ring3[gi] || ctx.ring_bond_count[gi] >= 3,
                        15 | 33 => true,
                        16 | 34 => {
                            let val: f64 = ctx.atom_bonds[gi]
                                .iter()
                                .map(|&ei| ctx.bonds[ei].2)
                                .sum::<f64>()
                                + ctx.n_h[gi] as f64;
                            val == 4.0 || (val == 3.0 && g.atoms[gi].formal_charge == 1)
                        }
                        _ => false,
                    }
                }
            } else {
                true
            };
        if !legal {
            continue;
        }
        // 隣接ランクの重複チェック
        let nbr_ranks: Vec<usize> = ctx.atom_bonds[gi]
            .iter()
            .map(|&ei| {
                let (a, b, _) = ctx.bonds[ei];
                ranks[if a == gi { b } else { a }]
            })
            .collect();
        let mut seen = nbr_ranks.clone();
        seen.sort_unstable();
        seen.dedup();
        if seen.len() != nbr_ranks.len() {
            continue; // hasDupes
        }

        // CIP ランク昇順に並べ替えるのに必要なスワップ数
        let mut sorted_bonds: Vec<(usize, usize)> = ctx.atom_bonds[gi]
            .iter()
            .map(|&ei| {
                let (a, b, _) = ctx.bonds[ei];
                (ranks[if a == gi { b } else { a }], ei)
            })
            .collect();
        sorted_bonds.sort_unstable();
        let probe: Vec<usize> = sorted_bonds.iter().map(|&(_, ei)| ei).collect();
        let mut n_swaps = swap_parity(&ctx.atom_bonds[gi], &probe);
        if probe.len() == 3 && ctx.n_h[gi] == 1 {
            n_swaps += 1;
        }
        let mut tag = tag;
        if n_swaps % 2 == 1 {
            tag = match tag {
                Chirality::Anticlockwise => Chirality::Clockwise,
                Chirality::Clockwise => Chirality::Anticlockwise,
            };
        }
        g.atoms[gi].chiral_tag = Some(match tag {
            Chirality::Anticlockwise => 'S',
            Chirality::Clockwise => 'R',
        });
    }
}

/// assignBondStereoCodes 相当。bonds[ei].stereo に 'E'/'Z' を設定。
fn assign_bond_stereo_codes(g: &mut MoleculeGraph, ctx: &Ctx, ranks: &[usize]) {
    // 方向結合 (グラフ結合 id → begin→end 向きの Up/Down)
    let mut dirs: Vec<Option<BondKind>> = vec![None; ctx.bonds.len()];
    {
        let mut bond_map: HashMap<(usize, usize), usize> = HashMap::new();
        for (ei, &(a, b, _)) in ctx.bonds.iter().enumerate() {
            bond_map.insert((a.min(b), a.max(b)), ei);
        }
        for b in &g.parsed.bonds {
            if !matches!(b.kind, BondKind::Up | BondKind::Down) {
                continue;
            }
            let (Some(x), Some(y)) = (g.parser_to_graph[b.a], g.parser_to_graph[b.b]) else {
                continue;
            };
            let Some(&ei) = bond_map.get(&(x.min(y), x.max(y))) else {
                continue;
            };
            // parsed は a→b 向き。グラフ結合の begin→end と逆なら反転
            let (gb, _, _) = ctx.bonds[ei];
            let kind = if gb == x {
                b.kind
            } else {
                match b.kind {
                    BondKind::Up => BondKind::Down,
                    BondKind::Down => BondKind::Up,
                    k => k,
                }
            };
            dirs[ei] = Some(kind);
        }
    }

    // findAtomNeighborDirHelper 相当
    let neighbor_dirs = |atom: usize, ref_ei: usize| -> Vec<(usize, Option<BondKind>)> {
        let mut nbrs: Vec<(usize, Option<BondKind>)> = Vec::new();
        let mut seen_dir = false;
        for &ei in &ctx.atom_bonds[atom] {
            if ei == ref_ei {
                continue;
            }
            let (a, b, _) = ctx.bonds[ei];
            let mut dir = dirs[ei];
            if let Some(d) = dir {
                seen_dir = true;
                // 原子から見て逆向きなら反転
                if atom != a {
                    dir = Some(match d {
                        BondKind::Up => BondKind::Down,
                        BondKind::Down => BondKind::Up,
                        k => k,
                    });
                }
            }
            nbrs.push((if a == atom { b } else { a }, dir));
        }
        if !seen_dir {
            return Vec::new();
        }
        if nbrs.len() == 2 && ranks[nbrs[0].0] == ranks[nbrs[1].0] {
            return Vec::new(); // 置換基が同一 → 立体なし
        }
        // 片側だけ方向がある場合は逆向きを補完
        if nbrs[0].1.is_none() && nbrs.len() > 1 {
            nbrs[0].1 = Some(match nbrs[1].1.expect("one dir set") {
                BondKind::Up => BondKind::Down,
                _ => BondKind::Up,
            });
        } else if nbrs.len() > 1 && nbrs[1].1.is_none() {
            nbrs[1].1 = Some(match nbrs[0].1.expect("one dir set") {
                BondKind::Up => BondKind::Down,
                _ => BondKind::Up,
            });
        }
        nbrs
    };

    for ei in 0..ctx.bonds.len() {
        let (beg, end, order) = ctx.bonds[ei];
        if order != 2.0 {
            continue;
        }
        // 8 員未満の環内二重結合は対象外
        if let Some(sz) = ctx.min_bond_ring[ei] {
            if sz < 8 {
                continue;
            }
        }
        let (db, de) = (ctx.atom_bonds[beg].len(), ctx.atom_bonds[end].len());
        if !(2..=3).contains(&db) || !(2..=3).contains(&de) {
            continue;
        }
        let beg_nbrs = neighbor_dirs(beg, ei);
        let end_nbrs = neighbor_dirs(end, ei);
        if beg_nbrs.is_empty() || end_nbrs.is_empty() {
            continue;
        }
        // 各側の最高ランク近傍を選ぶ
        let pick = |nbrs: &[(usize, Option<BondKind>)]| -> (usize, BondKind) {
            if nbrs.len() == 1 || ranks[nbrs[0].0] > ranks[nbrs[1].0] {
                (nbrs[0].0, nbrs[0].1.expect("dir"))
            } else {
                (nbrs[1].0, nbrs[1].1.expect("dir"))
            }
        };
        // 両隣接が同方向なら矛盾 → 立体なし
        let conflicting =
            |nbrs: &[(usize, Option<BondKind>)]| nbrs.len() == 2 && nbrs[0].1 == nbrs[1].1;
        if conflicting(&beg_nbrs) || conflicting(&end_nbrs) {
            continue;
        }
        let (_, beg_dir) = pick(&beg_nbrs);
        let (_, end_dir) = pick(&end_nbrs);
        g.bonds[ei].stereo = Some(if beg_dir == end_dir { 'Z' } else { 'E' });
    }
}

/// CIP ランクを返す (kept 原子のみ)。
/// 配座生成 (conformer) が E/Z 幾何拘束の置換基選択に使う。
pub(crate) fn cip_ranks(g: &MoleculeGraph) -> Vec<usize> {
    let ctx = build_ctx(g);
    assign_cip_ranks(&ctx)
}

/// 立体化学の割当て (build_molecule_graph の最終段で呼ぶ)。
pub(crate) fn assign_stereochemistry(g: &mut MoleculeGraph) {
    let has_chirality = g.parsed.atoms.iter().any(|a| a.chirality.is_some());
    let has_dirs = g
        .parsed
        .bonds
        .iter()
        .any(|b| matches!(b.kind, BondKind::Up | BondKind::Down));
    if !has_chirality && !has_dirs {
        return;
    }
    let ctx = build_ctx(g);
    let ranks = assign_cip_ranks(&ctx);
    if has_chirality {
        assign_atom_chiral_codes(g, &ctx, &ranks);
    }
    if has_dirs {
        assign_bond_stereo_codes(g, &ctx, &ranks);
    }
}

#[cfg(test)]
mod tests {
    use crate::graph::build_molecule_graph;

    fn cip(smiles: &str) -> Vec<(usize, char)> {
        let g = build_molecule_graph(smiles).expect("valid");
        g.atoms
            .iter()
            .filter_map(|a| a.chiral_tag.map(|c| (a.idx, c)))
            .collect()
    }

    fn ez(smiles: &str) -> Vec<char> {
        let g = build_molecule_graph(smiles).expect("valid");
        g.bonds.iter().filter_map(|b| b.stereo).collect()
    }

    #[test]
    fn alanine() {
        // L-アラニン = (S)、D-アラニン = (R)
        assert_eq!(cip("N[C@@H](C)C(=O)O"), vec![(1, 'S')]);
        assert_eq!(cip("N[C@H](C)C(=O)O"), vec![(1, 'R')]);
    }

    #[test]
    fn smiles_start_chirality() {
        // 先頭原子 + 明示 H の特例 (RDKit 実測):
        // [C@@H](Cl)(F)Br は S で、Cl[C@H](F)Br と同じ分子 (Cl[C@@H] は R)
        assert_eq!(cip("[C@@H](Cl)(F)Br"), vec![(0, 'S')]);
        assert_eq!(cip("Cl[C@@H](F)Br"), vec![(1, 'R')]);
        assert_eq!(cip("Cl[C@H](F)Br"), vec![(1, 'S')]);
    }

    #[test]
    fn ring_closure_chirality() {
        // 環閉じ数字を含むキラル中心 (RDKit 実測: S)
        assert_eq!(cip("C[C@H]1CCCO1"), vec![(1, 'S')]);
    }

    #[test]
    fn non_stereocenter_cleaned() {
        // 同一置換基 2 つ → CIP コードなし
        assert!(cip("C[C@H](C)O").is_empty());
    }

    #[test]
    fn double_bond_stereo() {
        assert_eq!(ez("F/C=C/F"), vec!['E']);
        assert_eq!(ez("F/C=C\\F"), vec!['Z']);
        assert_eq!(ez("C/C=C\\C"), vec!['Z']);
        // 方向結合なし → 立体なし
        assert!(ez("CC=CC").is_empty());
        // 小さい環内の二重結合は対象外
        assert!(ez("C1/C=C\\CC1").is_empty());
    }
}
