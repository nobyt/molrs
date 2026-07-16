//! 環認識 (S1.4): RDKit の findSSSR + symmetrizeSSSR の忠実な移植。
//!
//! RDKit の `MolFromSmiles` はサニタイズで `symmetrizeSSSR` を実行するため、
//! `GetRingInfo().AtomRings()` は SSSR に「同サイズで置換可能な余剰環」を
//! 末尾に加えた対称化 SSSR になる (ビシクロ[2.2.2]オクタンは 3 環、
//! キュバンは 6 環)。環リストの順序と環内の原子順は Figueras 系 BFS の
//! 探索順そのものなので、Python 版命名ロジックとの互換のため
//! RDKit (Release_2023_09) の FindRings.cpp をアルゴリズムごと移植する。
//!
//! 前提: `bonds` は RDKit 互換順 (graph.rs が構築する順) であること。
//! 隣接リストの走査順が結果の原子順を決めるため、結合リスト順は本質的。

use std::collections::{BTreeMap, BTreeSet};

const NIL: usize = usize::MAX;

/// 環の原子集合ビットマスク (RDKit RINGINVAR 相当; 数値比較の順序も一致)
type RingInvar = u128;

fn ring_invariant(ring: &[usize]) -> RingInvar {
    let mut m: RingInvar = 0;
    for &a in ring {
        m |= 1 << a;
    }
    m
}

struct Graph {
    n_atoms: usize,
    n_bonds: usize,
    bonds: Vec<(usize, usize)>,
    /// atom → [(相手, bond_idx)] 結合リスト順
    adj: Vec<Vec<(usize, usize)>>,
}

impl Graph {
    fn new(n_atoms: usize, bonds: &[(usize, usize)]) -> Self {
        assert!(n_atoms <= 128, "molecule too large for ring perception");
        let mut adj = vec![Vec::new(); n_atoms];
        for (ei, &(a, b)) in bonds.iter().enumerate() {
            adj[a].push((b, ei));
            adj[b].push((a, ei));
        }
        Graph {
            n_atoms,
            n_bonds: bonds.len(),
            bonds: bonds.to_vec(),
            adj,
        }
    }

    fn bond_between(&self, u: usize, v: usize) -> Option<usize> {
        self.adj[u]
            .iter()
            .find(|&&(w, _)| w == v)
            .map(|&(_, ei)| ei)
    }

    /// 環 (原子列) → 結合インデックス列 (RDKit convertToBonds 相当)
    fn to_bond_ring(&self, ring: &[usize]) -> Vec<usize> {
        let mut br = Vec::with_capacity(ring.len());
        for i in 0..ring.len() - 1 {
            br.push(self.bond_between(ring[i], ring[i + 1]).expect("ring bond"));
        }
        br.push(
            self.bond_between(ring[ring.len() - 1], ring[0])
                .expect("ring bond"),
        );
        br
    }
}

/// RDKit `getMolFrags` 相当: 連結成分ごとの原子リスト (昇順)。
fn mol_frags(g: &Graph) -> Vec<Vec<usize>> {
    let mut comp = vec![NIL; g.n_atoms];
    let mut n_comp = 0;
    for start in 0..g.n_atoms {
        if comp[start] != NIL {
            continue;
        }
        let mut stack = vec![start];
        comp[start] = n_comp;
        while let Some(u) = stack.pop() {
            for &(v, _) in &g.adj[u] {
                if comp[v] == NIL {
                    comp[v] = n_comp;
                    stack.push(v);
                }
            }
        }
        n_comp += 1;
    }
    let mut frags = vec![Vec::new(); n_comp];
    for a in 0..g.n_atoms {
        frags[comp[a]].push(a);
    }
    frags
}

/// RDKit `trimBonds` 相当。
fn trim_bonds(
    cand: usize,
    g: &Graph,
    changed: &mut BTreeSet<usize>,
    atom_degrees: &mut [i32],
    active_bonds: &mut [bool],
) {
    for &(other, ei) in &g.adj[cand] {
        if !active_bonds[ei] {
            continue;
        }
        if atom_degrees[other] <= 2 {
            changed.insert(other);
        }
        active_bonds[ei] = false;
        atom_degrees[other] -= 1;
        atom_degrees[cand] -= 1;
    }
}

/// RDKit `markUselessD2s` 相当 (再帰)。
fn mark_useless_d2s(
    root: usize,
    g: &Graph,
    forb: &mut [bool],
    atom_degrees: &[i32],
    active_bonds: &[bool],
) {
    for &(other, ei) in &g.adj[root] {
        if !active_bonds[ei] {
            continue;
        }
        if !forb[other] && atom_degrees[other] == 2 {
            forb[other] = true;
            mark_useless_d2s(other, g, forb, atom_degrees, active_bonds);
        }
    }
}

/// RDKit `pickD2Nodes` 相当。
fn pick_d2_nodes(
    g: &Graph,
    cur_frag: &[usize],
    atom_degrees: &[i32],
    active_bonds: &[bool],
) -> Vec<usize> {
    let mut d2nodes = Vec::new();
    let mut forb = vec![false; g.n_atoms];
    loop {
        let mut root = NIL;
        for &axci in cur_frag {
            if atom_degrees[axci] == 2 && !forb[axci] {
                root = axci;
                d2nodes.push(axci);
                forb[axci] = true;
                break;
            }
        }
        if root == NIL {
            break;
        }
        mark_useless_d2s(root, g, &mut forb, atom_degrees, active_bonds);
    }
    d2nodes
}

/// RDKit `BFSWorkspace::smallestRingsBfs` 相当。
/// root を含む最小環を全て見つける (経路復元順も一致させる)。
fn smallest_rings_bfs(
    g: &Graph,
    root: usize,
    rings: &mut Vec<Vec<usize>>,
    active_bonds: &[bool],
    forbidden: Option<&[usize]>,
) -> usize {
    const WHITE: u8 = 0;
    const GRAY: u8 = 1;
    const BLACK: u8 = 2;
    let mut done = vec![WHITE; g.n_atoms];
    if let Some(forb) = forbidden {
        for &i in forb {
            done[i] = BLACK;
        }
    }
    let mut parents = vec![NIL; g.n_atoms];
    let mut depths = vec![0usize; g.n_atoms];
    let mut bfsq = std::collections::VecDeque::new();
    bfsq.push_back(root);

    let mut cur_size = usize::MAX;
    while let Some(curr) = bfsq.pop_front() {
        done[curr] = BLACK;
        let depth = depths[curr] + 1;
        if depth > cur_size {
            break;
        }
        for &(nbr, ei) in &g.adj[curr] {
            if !active_bonds[ei] {
                continue;
            }
            if done[nbr] == BLACK || parents[curr] == nbr {
                continue;
            }
            if done[nbr] == WHITE {
                parents[nbr] = curr;
                done[nbr] = GRAY;
                depths[nbr] = depth;
                bfsq.push_back(nbr);
            } else {
                // 別経路で到達済み → 環閉合の可能性。2 経路を縫い合わせる
                let mut ring = vec![nbr];
                // forwards path (nbr の祖先、root は含めない)
                let mut parent = parents[nbr];
                while parent != NIL && parent != root {
                    ring.push(parent);
                    parent = parents[parent];
                }
                // backwards path (curr の祖先を先頭に挿入、root を含む)
                ring.insert(0, curr);
                parent = parents[curr];
                while parent != NIL {
                    // 最小共通祖先が root でなければ環ではない
                    if ring.contains(&parent) {
                        ring.clear();
                        break;
                    }
                    ring.insert(0, parent);
                    parent = parents[parent];
                }
                if ring.len() > 1 {
                    if ring.len() <= cur_size {
                        cur_size = ring.len();
                        rings.push(ring);
                    } else {
                        return rings.len();
                    }
                }
            }
        }
    }
    rings.len()
}

/// findRingsD2nodes 内で環を res に登録する共通処理。
fn add_ring_if_new(
    g: &Graph,
    ring: &[usize],
    res: &mut Vec<Vec<usize>>,
    invars: &mut BTreeSet<RingInvar>,
    ring_bonds: Option<&mut Vec<bool>>,
    ring_atoms: Option<&mut Vec<bool>>,
) -> bool {
    let invr = ring_invariant(ring);
    if invars.contains(&invr) {
        return false;
    }
    res.push(ring.to_vec());
    invars.insert(invr);
    if let (Some(rb), Some(ra)) = (ring_bonds, ring_atoms) {
        for i in 0..ring.len() - 1 {
            let bi = g.bond_between(ring[i], ring[i + 1]).expect("ring bond");
            rb[bi] = true;
            ra[ring[i]] = true;
        }
        let bi = g
            .bond_between(ring[0], ring[ring.len() - 1])
            .expect("ring bond");
        rb[bi] = true;
        ra[ring[ring.len() - 1]] = true;
    }
    true
}

/// RDKit `findSSSRforDupCands` 相当。
#[allow(clippy::too_many_arguments)]
fn find_sssr_for_dup_cands(
    g: &Graph,
    res: &mut Vec<Vec<usize>>,
    invars: &mut BTreeSet<RingInvar>,
    dup_map: &BTreeMap<usize, Vec<usize>>,
    dup_d2_cands: &BTreeMap<RingInvar, Vec<usize>>,
    atom_degrees: &[i32],
    active_bonds: &[bool],
) {
    for dup_cands in dup_d2_cands.values() {
        if dup_cands.len() <= 1 {
            continue;
        }
        let mut nrings: Vec<Vec<usize>> = Vec::new();
        let mut min_siz = usize::MAX;
        for &dup_cand in dup_cands {
            let mut degrees_copy = atom_degrees.to_vec();
            let mut active_copy = active_bonds.to_vec();
            let mut changed = BTreeSet::new();
            for &dni in dup_map.get(&dup_cand).expect("duplicate in dupMap") {
                trim_bonds(dni, g, &mut changed, &mut degrees_copy, &mut active_copy);
            }
            let mut srings = Vec::new();
            smallest_rings_bfs(g, dup_cand, &mut srings, &active_copy, None);
            for sring in srings {
                min_siz = min_siz.min(sring.len());
                nrings.push(sring);
            }
        }
        for nring in &nrings {
            if nring.len() == min_siz {
                add_ring_if_new(g, nring, res, invars, None, None);
            }
        }
    }
}

/// RDKit `findRingsD2nodes` 相当。
#[allow(clippy::too_many_arguments)]
fn find_rings_d2_nodes(
    g: &Graph,
    res: &mut Vec<Vec<usize>>,
    invars: &mut BTreeSet<RingInvar>,
    d2nodes: &[usize],
    atom_degrees: &mut [i32],
    active_bonds: &mut [bool],
    ring_bonds: &mut Vec<bool>,
    ring_atoms: &mut Vec<bool>,
) {
    let mut dup_d2_cands: BTreeMap<RingInvar, Vec<usize>> = BTreeMap::new();
    let mut dup_map: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    let mut node_invars: BTreeMap<usize, Vec<RingInvar>> = BTreeMap::new();

    for &cand in d2nodes {
        let mut srings = Vec::new();
        smallest_rings_bfs(g, cand, &mut srings, active_bonds, None);
        let srings_empty = srings.is_empty();
        for nring in &srings {
            let invr = ring_invariant(nring);
            add_ring_if_new(g, nring, res, invars, Some(ring_bonds), Some(ring_atoms));

            node_invars.entry(cand).or_default().push(invr);
            for (&node, invs) in &node_invars {
                if node != cand && invs.contains(&invr) {
                    dup_map.entry(cand).or_default().push(node);
                    dup_map.entry(node).or_default().push(cand);
                }
            }
            dup_d2_cands.entry(invr).or_default().push(cand);
        }

        // 環が見つからなかった場合のみ、この場で刈り取る (Issue 134 対応)
        if srings_empty {
            let mut changed: BTreeSet<usize> = BTreeSet::from([cand]);
            while let Some(&c) = changed.iter().next() {
                changed.remove(&c);
                trim_bonds(c, g, &mut changed, atom_degrees, active_bonds);
            }
        }
    }

    find_sssr_for_dup_cands(
        g,
        res,
        invars,
        &dup_map,
        &dup_d2_cands,
        atom_degrees,
        active_bonds,
    );
}

/// RDKit `findRingsD3Node` 相当。
fn find_rings_d3_node(
    g: &Graph,
    res: &mut Vec<Vec<usize>>,
    invars: &mut BTreeSet<RingInvar>,
    cand: usize,
    active_bonds: &[bool],
) {
    let mut srings = Vec::new();
    let nsmall = smallest_rings_bfs(g, cand, &mut srings, active_bonds, None);
    for nring in &srings {
        add_ring_if_new(g, nring, res, invars, None, None);
    }
    if nsmall >= 3 {
        return;
    }

    // 次数 3 ノードの活性な隣接 3 つ
    let nbrs: Vec<usize> = g.adj[cand]
        .iter()
        .filter(|&&(_, ei)| active_bonds[ei])
        .map(|&(v, _)| v)
        .collect();
    let (n1, n2, n3) = (nbrs[0], nbrs[1], nbrs[2]);

    if nsmall == 2 {
        // 2 環に共通する隣接原子 f を禁止して 3 つ目の環を探す
        let f = [n1, n2, n3]
            .into_iter()
            .find(|&x| srings[0].contains(&x) && srings[1].contains(&x))
            .expect("third ring not found");
        let mut trings = Vec::new();
        smallest_rings_bfs(g, cand, &mut trings, active_bonds, Some(&[f]));
        for nring in &trings {
            add_ring_if_new(g, nring, res, invars, None, None);
        }
    } else if nsmall == 1 {
        // 見つかった環に入っていない隣接 f1/f2 それぞれを含む環を探す
        let (f1, f2) = if !srings[0].contains(&n1) {
            (n2, n3)
        } else if !srings[0].contains(&n2) {
            (n1, n3)
        } else {
            (n1, n2)
        };
        let mut trings = Vec::new();
        smallest_rings_bfs(g, cand, &mut trings, active_bonds, Some(&[f2]));
        for nring in &trings {
            add_ring_if_new(g, nring, res, invars, None, None);
        }
        trings.clear();
        smallest_rings_bfs(g, cand, &mut trings, active_bonds, Some(&[f1]));
        for nring in &trings {
            add_ring_if_new(g, nring, res, invars, None, None);
        }
    }
}

/// RDKit `removeExtraRings` 相当。余剰環を除去して返す (symmetrize 用)。
fn remove_extra_rings(g: &Graph, res: &mut Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    // サイズで安定ソート (libstdc++ std::sort は小規模列で挿入ソート = 安定)
    res.sort_by_key(|r| r.len());

    let brings: Vec<Vec<usize>> = res.iter().map(|r| g.to_bond_ring(r)).collect();
    let bit_brings: Vec<u128> = brings
        .iter()
        .map(|br| {
            let mut m: u128 = 0;
            for &bi in br {
                m |= 1 << bi;
            }
            m
        })
        .collect();

    let n = res.len();
    let mut avail = vec![true; n];
    let mut keep = vec![false; n];
    let mut munion: u128 = 0;

    for i in 0..n {
        if bit_brings[i] & !munion == 0 {
            avail[i] = false;
        }
        if !avail[i] {
            continue;
        }
        munion |= bit_brings[i];
        keep[i] = true;

        let mut consider = vec![false; n];
        for j in i + 1..n {
            if avail[j] && brings[j].len() == brings[i].len() {
                consider[j] = true;
            }
        }
        while consider.iter().any(|&c| c) {
            let mut best_j = i + 1;
            let mut best_overlap: i64 = -1;
            let mut j = i + 1;
            while j < n && brings[j].len() == brings[i].len() {
                if consider[j] && avail[j] {
                    let overlap = (bit_brings[j] & munion).count_ones() as i64;
                    if overlap > best_overlap {
                        best_overlap = overlap;
                        best_j = j;
                    }
                }
                j += 1;
            }
            consider[best_j] = false;
            if bit_brings[best_j] & !munion == 0 {
                avail[best_j] = false;
            } else {
                keep[best_j] = true;
                avail[best_j] = false;
                munion |= bit_brings[best_j];
            }
        }
    }

    let mut extras = Vec::new();
    let temp = std::mem::take(res);
    for (i, ring) in temp.into_iter().enumerate() {
        if keep[i] {
            res.push(ring);
        } else {
            extras.push(ring);
        }
    }
    extras
}

/// RDKit `_atomSearchBFS` 相当 (縮合環系フォールバック用)。
fn atom_search_bfs(
    g: &Graph,
    start: usize,
    end: usize,
    ring_atoms: &[bool],
    invars: &BTreeSet<RingInvar>,
) -> Option<Vec<usize>> {
    let mut bfsq = std::collections::VecDeque::new();
    bfsq.push_back(vec![start]);
    while let Some(tv) = bfsq.pop_front() {
        let curr = *tv.last().expect("path");
        for &(nbr, _) in &g.adj[curr] {
            if nbr == end {
                if curr != start {
                    let mut nv = tv.clone();
                    nv.push(nbr);
                    if !invars.contains(&ring_invariant(&nv)) {
                        return Some(nv);
                    }
                }
            } else if ring_atoms[nbr] && !tv.contains(&nbr) {
                let mut nv = tv.clone();
                nv.push(nbr);
                bfsq.push_back(nv);
            }
        }
    }
    None
}

/// RDKit `findRingConnectingAtoms` 相当。
fn find_ring_connecting_atoms(
    g: &Graph,
    bond: (usize, usize),
    res: &mut Vec<Vec<usize>>,
    invars: &mut BTreeSet<RingInvar>,
    ring_bonds: &mut Vec<bool>,
    ring_atoms: &mut Vec<bool>,
) -> bool {
    if let Some(nring) = atom_search_bfs(g, bond.0, bond.1, ring_atoms, invars) {
        add_ring_if_new(g, &nring, res, invars, Some(ring_bonds), Some(ring_atoms))
    } else {
        false
    }
}

/// RDKit `fastFindRings` の DFS 相当 (最後のフォールバック)。
fn fast_find_rings(g: &Graph) -> Vec<Vec<usize>> {
    fn dfs(
        g: &Graph,
        atom: usize,
        colors: &mut [u8],
        order: &mut Vec<usize>,
        res: &mut Vec<Vec<usize>>,
        from: Option<usize>,
    ) {
        colors[atom] = 1;
        order.push(atom);
        for &(nbr, _) in &g.adj[atom] {
            if colors[nbr] == 0 {
                if g.adj[nbr].len() < 2 {
                    colors[nbr] = 2;
                } else {
                    dfs(g, nbr, colors, order, res, Some(atom));
                }
            } else if colors[nbr] == 1 {
                if let Some(f) = from {
                    if nbr != f {
                        let mut cycle = Vec::new();
                        let last = order.iter().rposition(|&x| x == atom).expect("in order");
                        for k in (0..=last).rev() {
                            if order[k] == nbr {
                                break;
                            }
                            cycle.push(order[k]);
                        }
                        cycle.push(nbr);
                        res.push(cycle);
                    }
                }
            }
        }
        colors[atom] = 2;
        order.pop();
    }

    let mut res = Vec::new();
    let mut colors = vec![0u8; g.n_atoms];
    for i in 0..g.n_atoms {
        if colors[i] != 0 {
            continue;
        }
        if g.adj[i].len() < 2 {
            colors[i] = 2;
            continue;
        }
        let mut order = Vec::new();
        dfs(g, i, &mut colors, &mut order, &mut res, None);
    }
    res
}

/// RDKit `findSSSR` 相当。返り値: (SSSR, 余剰環)。
fn find_sssr(g: &Graph) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
    let mut res: Vec<Vec<usize>> = Vec::new();
    let mut invars: BTreeSet<RingInvar> = BTreeSet::new();
    let mut active_bonds = vec![true; g.n_bonds];
    let mut ring_bonds = vec![false; g.n_bonds];
    let mut ring_atoms = vec![false; g.n_atoms];
    let mut atom_degrees: Vec<i32> = (0..g.n_atoms).map(|a| g.adj[a].len() as i32).collect();
    // extraRings プロパティ相当 (フラグメントごとに上書きされる RDKit の挙動を再現)
    let mut extras_prop: Vec<Vec<usize>> = Vec::new();

    for cur_frag in mol_frags(g) {
        if cur_frag.len() < 3 {
            continue;
        }
        let mut frag_res: Vec<Vec<usize>> = Vec::new();
        let mut changed: BTreeSet<usize> = BTreeSet::new();
        let mut nbnds = 0usize;
        for &a in &cur_frag {
            let deg = atom_degrees[a];
            nbnds += deg as usize;
            if deg < 2 {
                changed.insert(a);
            }
        }
        nbnds /= 2;
        if (nbnds as i64 - cur_frag.len() as i64 + 1) < 1 {
            continue;
        }

        let mut done_ats = vec![false; g.n_atoms];
        let mut n_atoms_done = 0usize;
        while n_atoms_done < cur_frag.len() {
            while let Some(&cand) = changed.iter().next() {
                changed.remove(&cand);
                if !done_ats[cand] {
                    done_ats[cand] = true;
                    n_atoms_done += 1;
                    trim_bonds(cand, g, &mut changed, &mut atom_degrees, &mut active_bonds);
                }
            }
            let d2nodes = pick_d2_nodes(g, &cur_frag, &atom_degrees, &active_bonds);
            if !d2nodes.is_empty() {
                find_rings_d2_nodes(
                    g,
                    &mut frag_res,
                    &mut invars,
                    &d2nodes,
                    &mut atom_degrees,
                    &mut active_bonds,
                    &mut ring_bonds,
                    &mut ring_atoms,
                );
                for &d2i in &d2nodes {
                    if !done_ats[d2i] {
                        done_ats[d2i] = true;
                        n_atoms_done += 1;
                    }
                    trim_bonds(d2i, g, &mut changed, &mut atom_degrees, &mut active_bonds);
                }
            } else if n_atoms_done < cur_frag.len() {
                let Some(&cand) = cur_frag.iter().find(|&&a| atom_degrees[a] == 3) else {
                    break;
                };
                find_rings_d3_node(g, &mut frag_res, &mut invars, cand, &active_bonds);
                done_ats[cand] = true;
                n_atoms_done += 1;
                trim_bonds(cand, g, &mut changed, &mut atom_degrees, &mut active_bonds);
            }
        }

        let nexpt = nbnds as i64 - cur_frag.len() as i64 + 1;
        let mut ssiz = frag_res.len() as i64;

        if ssiz < nexpt {
            // 高度に縮合した環系のフォールバック (RDKit Issue 3514824)。
            // RDKit と同じく、走査する結合インデックスの範囲は
            // このフラグメントの結合数 nbnds まで、という挙動を再現する。
            let scan = nbnds.min(g.n_bonds);
            let mut dead_bonds = vec![false; g.n_bonds];
            loop {
                let possible = (0..scan).find(|&i| {
                    !ring_bonds[i] && !dead_bonds[i] && {
                        let (a, b) = g.bonds[i];
                        ring_atoms[a] && ring_atoms[b]
                    }
                });
                let Some(bi) = possible else { break };
                let found = find_ring_connecting_atoms(
                    g,
                    g.bonds[bi],
                    &mut frag_res,
                    &mut invars,
                    &mut ring_bonds,
                    &mut ring_atoms,
                );
                if !found {
                    dead_bonds[bi] = true;
                }
            }
            ssiz = frag_res.len() as i64;
            if ssiz < nexpt {
                // 近似アルゴリズムへ切替 (RDKit fastFindRings フォールバック)
                return (fast_find_rings(g), Vec::new());
            }
        }
        if ssiz > nexpt {
            extras_prop = remove_extra_rings(g, &mut frag_res);
        }
        res.extend(frag_res);
    }
    (res, extras_prop)
}

/// RDKit `symmetrizeSSSR` 相当: `MolFromSmiles` 後の
/// `GetRingInfo().AtomRings()` と同じ環リストを返す。
pub fn symmetrized_sssr(n_atoms: usize, bonds: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let g = Graph::new(n_atoms, bonds);
    let (sssrs, extras) = find_sssr(&g);
    let mut res = sssrs.clone();
    if extras.is_empty() {
        return res;
    }

    let bond_sssrs: Vec<Vec<usize>> = sssrs.iter().map(|r| g.to_bond_ring(r)).collect();
    let mut bond_counts = vec![0i32; g.n_bonds];
    for r in &bond_sssrs {
        for &b in r {
            bond_counts[b] += 1;
        }
    }

    // 各余剰環について、同サイズの SSSR 環と置換可能なら追加
    for extra_atom_ring in &extras {
        let extra_ring = g.to_bond_ring(extra_atom_ring);
        for ring in &bond_sssrs {
            if ring.len() != extra_ring.len() {
                continue;
            }
            let mut share_bond = false;
            let mut replaces_all_unique_bonds = true;
            for &bond_id in ring {
                let bond_count = bond_counts[bond_id];
                if bond_count == 1 || !share_bond {
                    if extra_ring.contains(&bond_id) {
                        share_bond = true;
                    } else if bond_count == 1 {
                        replaces_all_unique_bonds = false;
                    }
                }
            }
            if share_bond && replaces_all_unique_bonds {
                res.push(extra_atom_ring.clone());
                break;
            }
        }
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn benzene_ring_order() {
        // c1ccccc1 の結合順: (0,1)..(4,5), 閉じ (5,0)
        let bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
        let rings = symmetrized_sssr(6, &bonds);
        assert_eq!(rings, vec![vec![0, 5, 4, 3, 2, 1]]); // RDKit 実測値
    }

    #[test]
    fn pyrrole_ring_order() {
        let bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)];
        let rings = symmetrized_sssr(5, &bonds);
        assert_eq!(rings, vec![vec![0, 1, 2, 3, 4]]); // RDKit 実測値
    }

    #[test]
    fn bicyclo222_symmetrized() {
        // BrC1CC2CCC1CC2 相当の環部分 (原子 1..8): RDKit は 3 環を返す
        let bonds = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 8),
            (6, 1),
            (8, 3),
        ];
        let rings = symmetrized_sssr(9, &bonds);
        assert_eq!(rings.len(), 3);
        assert!(rings.iter().all(|r| r.len() == 6));
    }

    #[test]
    fn norbornane_no_extra() {
        // BrC1CC2CCC1C2: 5,5 の 2 環のみ (6 員環は追加されない)
        let bonds = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (6, 1),
            (7, 3),
        ];
        let rings = symmetrized_sssr(8, &bonds);
        assert_eq!(rings.len(), 2);
        assert!(rings.iter().all(|r| r.len() == 5));
    }

    #[test]
    fn chain_no_rings() {
        assert!(symmetrized_sssr(3, &[(0, 1), (1, 2)]).is_empty());
        assert!(symmetrized_sssr(0, &[]).is_empty());
    }
}
