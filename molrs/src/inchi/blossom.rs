//! 一般グラフの交互到達可能性 (Edmonds のブロッサム法、I14)。
//!
//! IUPAC 公式 InChI (`ichi_bns.c`) の可動 H 判定は Kocay–Stone の
//! balanced network flow (Edmonds ブロッサム法を容量付きに一般化したもの)
//! だが、molrs が扱う「1 個の原子は高々 1 本の二重結合を持つ」という
//! 標準的な有機分子の Kekule 構造では容量は常に 0/1 であり、通常の
//! (容量なし) 一般グラフマッチングのブロッサム法と数学的に等価になる。
//! 本モジュールはその標準形 (O(V^3) の古典アルゴリズム) を実装する。
//!
//! 用途: ある露出 (未マッチ) 頂点 `root` を可動 H の供与体とみなし、
//! そこから「交互パス」(自由辺 → マッチ辺 → 自由辺 → …) で到達できる
//! 「マッチ済み」頂点の集合を求める。これらは全て「この供与体の H が
//! 移動した先」になり得る受容体候補である。奇数長サイクル (5 員芳香
//! ヘテロ環など) を跨ぐ探索では素朴な交互 BFS では正しい可否判定が
//! できず、ブロッサム (奇閉路) の収縮が必須。

use std::collections::VecDeque;

/// 交互到達可能性判定用のグラフ (無向、頂点は `0..n` の整数)。
pub(crate) struct MatchGraph {
    n: usize,
    adj: Vec<Vec<usize>>,
}

impl MatchGraph {
    pub(crate) fn new(n: usize) -> Self {
        MatchGraph {
            n,
            adj: vec![Vec::new(); n],
        }
    }

    pub(crate) fn add_edge(&mut self, a: usize, b: usize) {
        if a == b {
            return;
        }
        if !self.adj[a].contains(&b) {
            self.adj[a].push(b);
        }
        if !self.adj[b].contains(&a) {
            self.adj[b].push(a);
        }
    }

    /// `a`-`b` 間に辺があるか。
    #[cfg(test)]
    fn has_edge(&self, a: usize, b: usize) -> bool {
        self.adj[a].contains(&b)
    }

    /// `root` (露出頂点) から交互パスで到達できる「マッチ済み」頂点の集合。
    /// `matched[v] = Some(u)` は v-u が現在マッチ (二重結合) であることを
    /// 示す (対称、`matched[u] == Some(v)`)。root 自身は結果に含まない。
    ///
    /// 標準の一般グラフ最大マッチング探索 (findPath) を、最初の露出頂点で
    /// 打ち切らずに最後まで実行し、「自由辺で到達しマッチ済みだった」頂点
    /// (= その辺で H を受け取れる受容体候補) を全て収集する形に変形した版。
    pub(crate) fn alternating_reachable(
        &self,
        matched: &[Option<usize>],
        root: usize,
    ) -> Vec<usize> {
        let n = self.n;
        let mut used = vec![false; n];
        // p[v] = v に到達したときの交互木上の親頂点 (自由辺の相手側)。
        let mut parent: Vec<Option<usize>> = vec![None; n];
        // base[v] = v が現在どのブロッサムに属しているか (代表頂点)。
        let mut base: Vec<usize> = (0..n).collect();

        let mut odd_labeled = vec![false; n];
        let mut odd_reached: Vec<usize> = Vec::new();

        used[root] = true;
        let mut queue = VecDeque::new();
        queue.push_back(root);

        while let Some(v) = queue.pop_front() {
            for &to in &self.adj[v] {
                // 同一ブロッサム内、または v 自身のマッチ辺は自由辺として
                // 辿らない (マッチ辺はブロッサム収縮や `to` からの折り返し
                // でのみ暗黙に使う)。
                if base[v] == base[to] || matched[v] == Some(to) {
                    continue;
                }
                if to == root || (matched[to].is_some() && parent[matched[to].unwrap()].is_some()) {
                    // 奇閉路 (ブロッサム) を検出 — 収縮する。奇閉路の性質上、
                    // ブロッサムに吸収された各頂点は「回転」によって
                    // どれもが未マッチ (露出) 側になり得るため、吸収された
                    // 頂点は (root 自身を除き) 到達可能な新規露出候補として
                    // 記録する — 通常の (非ブロッサム) 発見経路と同様に扱う。
                    let curbase = lca(matched, &parent, &base, v, to);
                    let mut in_blossom = vec![false; n];
                    mark_path(matched, &mut parent, &base, &mut in_blossom, v, curbase, to);
                    mark_path(matched, &mut parent, &base, &mut in_blossom, to, curbase, v);
                    for i in 0..n {
                        if in_blossom[base[i]] {
                            base[i] = curbase;
                            if i != root && !odd_labeled[i] {
                                odd_labeled[i] = true;
                                odd_reached.push(i);
                            }
                            if !used[i] {
                                used[i] = true;
                                queue.push_back(i);
                            }
                        }
                    }
                } else if parent[to].is_none() {
                    parent[to] = Some(v);
                    // 記録すべきは to 自身ではなく to の「現在のマッチ相手」
                    // (partner)。v-to に新たに二重結合を作ると to は二重結合
                    // が 1 本余るため、to は元々の相手 partner との結合を
                    // 手放す — つまり H を新たに持つ側になるのは partner の
                    // 方 (偶数長の交互パスで partner が露出頂点になる)。
                    // to が既に露出頂点 (matched[to] が None) の場合は
                    // 「2 つの独立した供与体を繋ぐ奇数長の増加路」に相当し、
                    // これは (両者とも H を失って新たな二重結合を得るだけで
                    // 系内の H 総数が減ってしまうため) 互変異性の移動としては
                    // 無効 — 記録も継続探索もしない。
                    if let Some(partner) = matched[to] {
                        if !odd_labeled[partner] {
                            odd_labeled[partner] = true;
                            odd_reached.push(partner);
                        }
                        if !used[partner] {
                            used[partner] = true;
                            queue.push_back(partner);
                        }
                    }
                }
            }
        }
        odd_reached
    }

    /// `matched` を初期マッチングとして増加路を繰り返し探し、最大マッチングの
    /// 辺数を返す (I24)。`skip[v]` が真の頂点はグラフから取り除かれたものとして
    /// 扱う (`matched` 側にも残っていてはならない — 呼び出し側で外すこと)。
    ///
    /// [`Self::alternating_reachable`] が「1 個の供与体から出る 1 本の交互路」
    /// しか見ないのに対し、こちらは完全マッチングの存在判定に使う。t-group の
    /// 仮想ハブ頂点を含むグラフでは全頂点がマッチ済みになるため「露出頂点から
    /// 探索」が使えず、「辺 (a, hub) を含む最大マッチングが存在するか」
    /// = 「`G − a − hub` の最大マッチングが 1 本だけ小さいか」で判定する。
    pub(crate) fn max_matching(&self, matched: &mut [Option<usize>], skip: &[bool]) -> usize {
        for root in 0..self.n {
            if skip[root] || matched[root].is_some() {
                continue;
            }
            if let Some((last, parent)) = self.find_augmenting_path(matched, skip, root) {
                // 増加路に沿ってマッチを反転する。
                let mut v = last;
                loop {
                    let pv = parent[v].expect("augmenting path vertex must have parent");
                    let ppv = matched[pv];
                    matched[v] = Some(pv);
                    matched[pv] = Some(v);
                    match ppv {
                        Some(next) => v = next,
                        None => break,
                    }
                }
            }
        }
        matched.iter().filter(|m| m.is_some()).count() / 2
    }

    /// 露出頂点 `root` から交互木を伸ばし、別の露出頂点で終わる増加路を探す
    /// (Edmonds のブロッサム法、標準形)。見つかれば (終点, 親配列) を返す。
    fn find_augmenting_path(
        &self,
        matched: &[Option<usize>],
        skip: &[bool],
        root: usize,
    ) -> Option<(usize, Vec<Option<usize>>)> {
        let n = self.n;
        let mut used = vec![false; n];
        let mut parent: Vec<Option<usize>> = vec![None; n];
        let mut base: Vec<usize> = (0..n).collect();

        used[root] = true;
        let mut queue = VecDeque::new();
        queue.push_back(root);

        while let Some(v) = queue.pop_front() {
            for &to in &self.adj[v] {
                if skip[to] || base[v] == base[to] || matched[v] == Some(to) {
                    continue;
                }
                if to == root || (matched[to].is_some() && parent[matched[to].unwrap()].is_some()) {
                    let curbase = lca(matched, &parent, &base, v, to);
                    let mut in_blossom = vec![false; n];
                    mark_path(matched, &mut parent, &base, &mut in_blossom, v, curbase, to);
                    mark_path(matched, &mut parent, &base, &mut in_blossom, to, curbase, v);
                    for i in 0..n {
                        if in_blossom[base[i]] {
                            base[i] = curbase;
                            if !used[i] {
                                used[i] = true;
                                queue.push_back(i);
                            }
                        }
                    }
                } else if parent[to].is_none() {
                    parent[to] = Some(v);
                    match matched[to] {
                        None => return Some((to, parent)),
                        Some(partner) => {
                            used[partner] = true;
                            queue.push_back(partner);
                        }
                    }
                }
            }
        }
        None
    }
}

/// v, to から交互木を根 (露出頂点) に向かって遡り、最初に共有する
/// ブロッサム基点 (最近共通祖先) を求める。
fn lca(
    matched: &[Option<usize>],
    parent: &[Option<usize>],
    base: &[usize],
    a: usize,
    b: usize,
) -> usize {
    let n = base.len();
    let mut on_path_a = vec![false; n];
    let mut x = a;
    loop {
        x = base[x];
        on_path_a[x] = true;
        match matched[x] {
            None => break,
            Some(m) => x = parent[m].expect("matched vertex on tree must have parent"),
        }
    }
    let mut y = b;
    loop {
        y = base[y];
        if on_path_a[y] {
            return y;
        }
        y = parent[matched[y].expect("unmatched vertex reached without hitting lca")]
            .expect("matched vertex on tree must have parent");
    }
}

/// v から基点 `b` まで交互パスを遡りながら、通過した頂点を全てブロッサム
/// (`in_blossom`) としてマークし、親ポインタを (子側優先で) 張り替える。
#[allow(clippy::too_many_arguments)]
fn mark_path(
    matched: &[Option<usize>],
    parent: &mut [Option<usize>],
    base: &[usize],
    in_blossom: &mut [bool],
    mut v: usize,
    b: usize,
    mut child: usize,
) {
    while base[v] != b {
        in_blossom[base[v]] = true;
        let m = matched[v].expect("odd vertex in blossom path must be matched");
        in_blossom[base[m]] = true;
        parent[v] = Some(child);
        child = m;
        v = parent[m].expect("matched vertex on tree must have parent");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// matched ペアの Vec から matched[] 配列を組み立てるテストヘルパ。
    fn matching(n: usize, pairs: &[(usize, usize)]) -> Vec<Option<usize>> {
        let mut m = vec![None; n];
        for &(a, b) in pairs {
            m[a] = Some(b);
            m[b] = Some(a);
        }
        m
    }

    #[test]
    fn simple_path() {
        // 0(露出)-1=2-3=4 (辺 0-1,1-2,2-3,3-4、マッチは 1-2,3-4)
        let mut g = MatchGraph::new(5);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        g.add_edge(2, 3);
        g.add_edge(3, 4);
        let m = matching(5, &[(1, 2), (3, 4)]);
        let mut reached = g.alternating_reachable(&m, 0);
        reached.sort_unstable();
        // 0-1 に新たな二重結合ができると 1 は元の相手 2 との結合を手放す
        // (2 が新たに露出) → さらに 2-3 で継続すると 3 は 4 との結合を
        // 手放す (4 が新たに露出)。「新たに露出しうる」頂点 = 2, 4
        // (1,3 自体は経路の中継点であって、そこが露出になるわけではない)。
        assert_eq!(reached, vec![2, 4]);
    }

    #[test]
    fn five_cycle_blossom() {
        // 5 員環 (奇閉路): 0-1-2-3-4-0。マッチは 1-2, 3-4。0 は露出。
        // 奇閉路の性質上、この環は「どの 1 頂点でも回転により未マッチに
        // なりうる」ブロッサムを成す。0 以外の全頂点 {1,2,3,4} が新たな
        // 露出候補として届く。
        let mut g = MatchGraph::new(5);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        g.add_edge(2, 3);
        g.add_edge(3, 4);
        g.add_edge(4, 0);
        let m = matching(5, &[(1, 2), (3, 4)]);
        let mut reached = g.alternating_reachable(&m, 0);
        reached.sort_unstable();
        assert_eq!(reached, vec![1, 2, 3, 4]);
    }

    #[test]
    fn blossom_enables_far_reach() {
        // 5 員環 0-1-2-3-4-0 (マッチ 1-2, 3-4、0 露出) に、2 から外部の
        // 露出していない鎖 2-5=6 (5-6 がマッチ) をぶら下げる。
        // 素朴な前方のみの BFS (ブロッサムなし) だと、0->1(match2)->2 に
        // 到達した時点で 2 は「マッチ済みの相手 1 経由」で処理済みになり、
        // 2 からさらに 5 へ自由辺で進む一般化ができない実装だと 6 に
        // 届かない。ブロッサム法では 0-1-2-3-4-0 の奇閉路を収縮した上で
        // 2 (ブロッサムの一部) からの辺 2-5 も探索対象になり、5 が 6 との
        // 結合を手放して 6 が新たに露出する経路まで届く。
        let mut g = MatchGraph::new(7);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        g.add_edge(2, 3);
        g.add_edge(3, 4);
        g.add_edge(4, 0);
        g.add_edge(2, 5);
        g.add_edge(5, 6);
        let m = matching(7, &[(1, 2), (3, 4), (5, 6)]);
        let mut reached = g.alternating_reachable(&m, 0);
        reached.sort_unstable();
        assert_eq!(reached, vec![1, 2, 3, 4, 6]);
    }

    #[test]
    fn no_edges_no_reach() {
        let g = MatchGraph::new(3);
        let m = matching(3, &[]);
        assert!(g.alternating_reachable(&m, 0).is_empty());
    }

    #[test]
    fn max_matching_augments_odd_path() {
        // 0-1-2-3-4 のパス。初期マッチが 1-2 だけでも最大 (2 本) まで増える。
        let mut g = MatchGraph::new(5);
        for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 4)] {
            g.add_edge(a, b);
        }
        let mut m = matching(5, &[(1, 2)]);
        assert_eq!(g.max_matching(&mut m, &[false; 5]), 2);
    }

    #[test]
    fn max_matching_handles_blossom() {
        // 5 員環 0-1-2-3-4-0 に尻尾 0-5。完全マッチング (3 本) が存在するが、
        // 素朴な交互 BFS では奇閉路の収縮なしに見つけられない。
        let mut g = MatchGraph::new(6);
        for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (0, 5)] {
            g.add_edge(a, b);
        }
        let mut m = matching(6, &[(1, 2)]);
        assert_eq!(g.max_matching(&mut m, &[false; 6]), 3);
    }

    #[test]
    fn max_matching_skips_removed_vertices() {
        // 5 頂点の完全マッチングから 1 頂点を除くと奇数頂点になり 1 本減る。
        let mut g = MatchGraph::new(6);
        for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)] {
            g.add_edge(a, b);
        }
        let mut skip = [false; 6];
        skip[0] = true;
        let mut m = vec![None; 6];
        assert_eq!(g.max_matching(&mut m, &skip), 2);
        assert!(m[0].is_none());
    }

    #[test]
    fn add_edge_is_symmetric_and_dedups() {
        let mut g = MatchGraph::new(2);
        g.add_edge(0, 1);
        g.add_edge(1, 0);
        assert!(g.has_edge(0, 1));
        assert!(g.has_edge(1, 0));
        assert_eq!(g.adj[0].len(), 1);
        assert_eq!(g.adj[1].len(), 1);
    }
}
