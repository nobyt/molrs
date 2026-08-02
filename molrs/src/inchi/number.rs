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

use super::blossom;
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

fn is_hetero(sym: &str) -> bool {
    matches!(sym, "O" | "S" | "Se" | "Te" | "N")
}

fn n_h_of(g: &MoleculeGraph, i: usize) -> usize {
    g.adjacency[i]
        .iter()
        .filter(|&&x| g.atoms[x].symbol == "H")
        .count()
}

/// Kekule 化済みの結合次数 (芳香族結合も 1/2 の実値。`g.bond_orders` は
/// 芳香族を 1.5 で保持するため、可動 H 判定 (二重結合受容体の検出) には
/// こちらを使う)。
fn kekule_order_map(g: &MoleculeGraph) -> std::collections::HashMap<(usize, usize), f64> {
    let mut m = std::collections::HashMap::with_capacity(g.bonds.len());
    for (bi, b) in g.bonds.iter().enumerate() {
        let key = (b.begin_idx.min(b.end_idx), b.begin_idx.max(b.end_idx));
        m.insert(key, g.kekule_bond_orders[bi]);
    }
    m
}

/// 原子 `atom` の重原子結合次数の合計 (`chem_bonds_valence` 相当、H は
/// 含まない)。
fn chem_bonds_valence(
    g: &MoleculeGraph,
    kekule: &std::collections::HashMap<(usize, usize), f64>,
    atom: usize,
) -> f64 {
    g.adjacency[atom]
        .iter()
        .filter(|&&x| g.atoms[x].symbol != "H")
        .map(|&nb| {
            kekule
                .get(&(atom.min(nb), atom.max(nb)))
                .copied()
                .unwrap_or(1.0)
        })
        .sum()
}

/// 原子が (集約的に見て)「受容体」(二重結合を 1 本持ち、それを手放せば H を
/// 受け取れる) かどうか。結合次数合計 − 次数 = 1、すなわちどこか 1 本
/// 二重結合を持っている、を判定する (どの隣接原子との結合が二重かは問わない)。
fn is_acceptor_agg(
    g: &MoleculeGraph,
    kekule: &std::collections::HashMap<(usize, usize), f64>,
    atom: usize,
) -> bool {
    let deg = g.adjacency[atom]
        .iter()
        .filter(|&&x| g.atoms[x].symbol != "H")
        .count() as f64;
    (chem_bonds_valence(g, kekule, atom) - deg - 1.0).abs() < 1e-9
}

/// `center` から見て `nb` が受容体端点として使えるか。
///
/// `nb` が芳香族なら、IUPAC 公式 InChI (`ichitaut.c` の
/// `nGetEndpointInfo`) と同様に **どの隣接原子との間の結合が二重か** では
/// なく [`is_acceptor_agg`] (nb 自身の価数余裕) だけで判定する。芳香環内
/// では「どの結合が Kekule 二重か」は環の中で任意に選んだ 1 通りの共鳴
/// 構造に過ぎず、環内のどのヘテロ原子も「(そのヘテロ原子自身の) 二重結合を
/// どこかに持っている」なら常に受容体候補になり得る。center-nb の特定の
/// 結合が二重かどうかで判定すると (旧実装のバグ)、たまたま Kekule 化で
/// 選ばれた向きにだけ依存してしまい、縮環 (インダゾール/アザインドール型)
/// で片方の異性体だけ正しく判定できない非対称な誤りが生じる。
///
/// `nb` が非芳香族なら、Kekule 選択の曖昧さがない (孤立した二重結合や
/// 鎖状の官能基は 1 通りにしか描けない) ので、`is_acceptor_agg` を使うと
/// **無関係な遠くの二重結合** (例: スルホキシド S=O から見て、同じ S に
/// 単結合している全く別の -OH) まで誤って受容体扱いしてしまう。この場合は
/// center-nb の当該結合が実際に二重結合かどうかで判定する (旧実装のまま)。
fn is_acceptor_from(
    g: &MoleculeGraph,
    kekule: &std::collections::HashMap<(usize, usize), f64>,
    center: usize,
    nb: usize,
) -> bool {
    if g.atoms[nb].is_aromatic {
        is_acceptor_agg(g, kekule, nb)
    } else {
        let key = (center.min(nb), center.max(nb));
        kekule.get(&key).copied().unwrap_or(1.0) >= 2.0
    }
}

/// 単純な union-find (可動 H 群のブリッジ統合、I13 で使用)。
struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
        }
    }
    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[ra] = rb;
        }
    }
}

/// 中心原子 1 個の直接ヘテロ端点から「種」となる可動 H 群を検出する
/// (旧 I9 アルゴリズム本体)。返り値は (種群のリスト, 意図的除外された原子)。
/// 可動 H 数は [`mobile_groups`] 側で統合後にまとめて計算する。
///
/// 規則: ある中心原子に、ヘテロ原子 (O/S/Se/Te/N) が結合し、そのうち
/// 少なくとも 1 つが二重結合 (受容体)、単結合のものは H を持つか負電荷
/// (供与体) のとき、それらヘテロ原子端点で 1 群を作る。端点は末端でなくて
/// よい (アミド/ラクタムの N は次数 2)。カルボン酸・スルホン酸・アミド・
/// アミジン・グアニジン・尿素・ラクタムを覆う。
///
/// 2 番目の返り値「意図的除外」は、ある中心で有効な多端点候補
/// (endpoints.len()>=2 && has_double) が成立したにもかかわらず、酸除外規則
/// (カルバミン酸の N 等) で `chosen` から落とされた原子。[`mobile_groups`]
/// はこれらを孤立供与体としてブリッジ探索の起点に**しない** (探索すると
/// 除外規則を回避して間接的に元の酸素対と再結合してしまうため)。単に
/// どの中心の候補にもならなかった原子 (トリアゾール環の中心 N など) は
/// ここに含まれず、孤立供与体として正常にブリッジ探索の起点になる。
fn seed_groups(
    g: &MoleculeGraph,
    kekule: &std::collections::HashMap<(usize, usize), f64>,
) -> (Vec<Vec<usize>>, std::collections::HashSet<usize>) {
    let n = g.atoms.len();
    let mut used = vec![false; n];
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut excluded: std::collections::HashSet<usize> = std::collections::HashSet::new();
    // N=N-N(H)- (トリアゾール/テトラゾール環) のように中心自体が N のことも
    // ある。この場合も N 端点を許す (中心が O/S のような組み合わせは想定外
    // なので許可しない)。
    let center_is_c_or_n = |c: usize| matches!(g.atoms[c].symbol.as_str(), "C" | "N");
    let heavy_deg = |i: usize| {
        g.adjacency[i]
            .iter()
            .filter(|&&x| g.atoms[x].symbol != "H")
            .count()
    };
    for center in 0..n {
        if g.atoms[center].symbol == "H" {
            continue;
        }
        // 中心の二重結合 O 数 (非炭素中心の N 端点可否に使う)
        let n_double_o = g.adjacency[center]
            .iter()
            .filter(|&&nb| {
                g.atoms[nb].symbol == "O"
                    && kekule
                        .get(&(center.min(nb), center.max(nb)))
                        .copied()
                        .unwrap_or(1.0)
                        >= 2.0
            })
            .count();
        // 中心に結合したヘテロ原子端点を分類 (受容体 = 二重結合、供与体 = H/負電荷)。
        // N が端点になるのは中心が炭素・N (トリアゾール/テトラゾールの
        // N=N-N(H)- 等)、または (末端 N かつ 中心が二重結合 O ≥2 = スルホニル
        // 級) のとき。一級スルホンアミド NH2 は可動、二級・スルフィンアミドは
        // 非可動。
        let mut endpoints: Vec<usize> = Vec::new();
        let mut has_double = false;
        for &nb in &g.adjacency[center] {
            let sym = g.atoms[nb].symbol.as_str();
            if !is_hetero(sym) {
                continue;
            }
            if sym == "N" && !center_is_c_or_n(center) && !(heavy_deg(nb) == 1 && n_double_o >= 2) {
                continue;
            }
            if is_acceptor_from(g, kekule, center, nb) {
                endpoints.push(nb);
                has_double = true;
            } else if n_h_of(g, nb) >= 1 || g.atoms[nb].formal_charge < 0 {
                endpoints.push(nb);
            }
        }
        if !has_double || endpoints.len() < 2 {
            continue;
        }
        // O/S 端点だけで酸系 (二重 O/S ≥1 かつ 供与体 O/S ≥1) を成すなら、
        // N を除外して酸のみを群とする (カルバミン酸は O,O のみで N は固定)。
        let os_ep: Vec<usize> = endpoints
            .iter()
            .copied()
            .filter(|&e| g.atoms[e].symbol != "N")
            .collect();
        let os_double = os_ep
            .iter()
            .any(|&e| is_acceptor_from(g, kekule, center, e));
        let os_donor = os_ep
            .iter()
            .any(|&e| n_h_of(g, e) >= 1 || g.atoms[e].formal_charge < 0);
        let chosen: Vec<usize> = if os_double && os_donor && os_ep.len() >= 2 {
            for &e in &endpoints {
                if !os_ep.contains(&e) {
                    excluded.insert(e);
                }
            }
            os_ep
        } else {
            endpoints
        };
        if chosen.len() < 2 {
            continue;
        }
        // この中心だけでは H/負電荷を持つ端点がなくてもよい (例: 縮環の
        // 共有原子が両環のヘテロ原子と直接隣接するが、そのどちらも局所的
        // には H を持たない場合)。真の供与体は別の中心や孤立供与体経由の
        // ブリッジ探索で後から合流し、最終的な群の妥当性判定
        // (`total_h + total_neg > 0`) は [`mobile_groups`] 側で全ブリッジ
        // 後にまとめて行う。
        if chosen.iter().any(|&e| used[e]) {
            continue;
        }
        for &e in &chosen {
            used[e] = true;
        }
        groups.push(chosen);
    }
    (groups, excluded)
}

/// 可動 H 群 (互変異性) を検出する。返り値は (端点原子集合, 可動 H 数)。
///
/// I9 (単一中心の星型検出、[`seed_groups`]) に加え、I13/I14 で環をまたぐ
/// 多中心の互変異性も扱う。IUPAC 公式 InChI (`ichi_bns.c`) の可動 H 判定は
/// Kocay–Stone の容量付きバランスドネットワークフローだが、molrs が扱う
/// 「1 原子は高々二重結合 1 本」という通常の有機分子の Kekule 構造では
/// 容量は常に 0/1 であり、標準的な一般グラフマッチングの**ブロッサム法**
/// ([`blossom`]) と数学的に等価になる。芳香環の Kekule 構造を「マッチング」
/// とみなし、各供与体 (種の未マッチメンバー・孤立供与体) から交互到達可能な
/// 「マッチ済み原子の相手」(= H を新たに持ちうる位置) を求め union-find で
/// 統合する。単純な前方のみの交互探索 (旧 I13) は 5 員芳香ヘテロ環のような
/// 奇閉路をまたぐ判定を原理的に誤るため (ブロッサム収縮が必要)、単一 SSSR
/// 環へのロックという場当たり的な回避策が要っていたが、ブロッサム法は
/// これを正しく扱うため不要になった。
/// カルバミン酸系 (O,O 酸対に対して N を除外する規則) 等、種の段階での
/// 除外はそのままブリッジ探索にも伝播する (除外された N は種のどの
/// マッチ済みメンバーからも自由辺 1 本では到達できないため)。
pub(crate) fn mobile_groups(g: &MoleculeGraph) -> Vec<(Vec<usize>, u8)> {
    let n = g.atoms.len();
    let kekule = kekule_order_map(g);
    let (seeds, excluded) = seed_groups(g, &kekule);

    // 縮環の共有原子 (2 つ以上の SSSR 環に属する原子) は、探索グラフ上では
    // 環ごとに別頂点に分裂させる (頂点分割)。ある環の辺だけを通って共有
    // 原子に到達した探索が、そのまま別の環の辺へ「乗り換え」られないように
    // するため — インダゾール/アザインドール型で、ピロール型 N-H が縮環
    // 越しにピリジン型 N まで誤って到達しないための歯止め。芳香環をまたぐ
    // 正当な互変異性 (縮環の共有原子自体が両環のヘテロ原子と直接隣接する
    // 場合) は `seed_groups` の局所検出で既に拾われるため、ブリッジ探索側で
    // 縮環をまたぐ必要はない。ブロッサム法自体 ([`blossom`]) は変更しない
    // 汎用グラフアルゴリズムのまま — この頂点分割はグラフ構築側だけの工夫。
    let mut ring_membership: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (ri, ring) in g.ring_atom_sets.iter().enumerate() {
        for &a in ring {
            if a < n {
                ring_membership[a].push(ri);
            }
        }
    }
    let mut vertex_of: std::collections::HashMap<(usize, Option<usize>), usize> =
        std::collections::HashMap::new();
    let mut vertex_atom: Vec<usize> = Vec::new();
    #[allow(clippy::needless_range_loop)]
    for atom in 0..n {
        if g.atoms[atom].symbol == "H" {
            continue;
        }
        if ring_membership[atom].len() >= 2 {
            for &r in &ring_membership[atom] {
                let id = vertex_atom.len();
                vertex_of.insert((atom, Some(r)), id);
                vertex_atom.push(atom);
            }
        } else {
            let id = vertex_atom.len();
            vertex_of.insert((atom, None), id);
            vertex_atom.push(atom);
        }
    }
    let clones_of = |atom: usize| -> Vec<usize> {
        if ring_membership[atom].len() >= 2 {
            ring_membership[atom]
                .iter()
                .map(|&r| vertex_of[&(atom, Some(r))])
                .collect()
        } else {
            vec![vertex_of[&(atom, None)]]
        }
    };
    // 縮環の共有原子は指定した環の分身頂点、それ以外はその原子唯一の頂点。
    let clone_in = |atom: usize, ring: usize| -> usize {
        if ring_membership[atom].len() >= 2 {
            vertex_of[&(atom, Some(ring))]
        } else {
            vertex_of[&(atom, None)]
        }
    };
    let shared_rings = |u: usize, v: usize| -> Vec<usize> {
        ring_membership[u]
            .iter()
            .copied()
            .filter(|r| ring_membership[v].contains(r))
            .collect()
    };

    // ブリッジ探索用グラフ: 分子結合のうち、少なくとも一方の端が芳香族の
    // ものだけを辺として採用する (孤立アルケンやスルホニル中心越しの
    // 橋渡しを防ぐ — 二級スルホンアミド N やジヒドロピリジン環 N-H を
    // 誤って可動化しないための歯止め)。種メンバー (アミド N 等、非芳香族)
    // から芳香環への入口はこの条件で自然に含まれる。
    let mut graph = blossom::MatchGraph::new(vertex_atom.len());
    let mut matched: Vec<Option<usize>> = vec![None; vertex_atom.len()];
    for b in &g.bonds {
        let (u, v) = (b.begin_idx, b.end_idx);
        if g.atoms[u].symbol == "H" || g.atoms[v].symbol == "H" {
            continue;
        }
        let shared = shared_rings(u, v);
        let bo = kekule.get(&(u.min(v), u.max(v))).copied().unwrap_or(1.0);
        if bo >= 2.0 {
            // マッチ (二重結合) は共有環ごとの分身どうしを対応付ける
            // (縮環の共有辺なら両環の分身ペアそれぞれに設定)。
            if shared.is_empty() {
                let (cu, cv) = (clones_of(u)[0], clones_of(v)[0]);
                matched[cu] = Some(cv);
                matched[cv] = Some(cu);
            } else {
                for &r in &shared {
                    let (cu, cv) = (clone_in(u, r), clone_in(v, r));
                    matched[cu] = Some(cv);
                    matched[cv] = Some(cu);
                }
            }
        }
        if !(g.atoms[u].is_aromatic || g.atoms[v].is_aromatic) {
            continue;
        }
        if shared.is_empty() {
            for cu in clones_of(u) {
                for cv in clones_of(v) {
                    graph.add_edge(cu, cv);
                }
            }
        } else {
            for &r in &shared {
                graph.add_edge(clone_in(u, r), clone_in(v, r));
            }
        }
    }

    let mut uf = UnionFind::new(n);
    let mut members: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut roots: Vec<usize> = Vec::new();
    for grp in &seeds {
        for &m in &grp[1..] {
            uf.union(grp[0], m);
        }
        for &m in grp {
            members.insert(m);
            roots.push(m);
        }
    }
    // 孤立供与体 (どの中心の候補にもならなかった H/負電荷ヘテロ原子、例:
    // トリアゾール環の中心 N や環外アミノ基の N) もブリッジ探索の起点に
    // 加える。ただし `excluded` (酸除外規則で意図的に落とされた原子) は除く。
    // 非芳香族の環メンバーは除外する (ジヒドロピリジン環の N-H のように
    // 環内で隣の芳香環に単結合している場合でも、非芳香族環自体のメンバー
    // としては可動 H 対象にならない — 環外の -NH- 置換基や芳香環自身の
    // ヘテロ原子とは扱いが異なる)。
    for i in 0..n {
        if excluded.contains(&i) || members.contains(&i) {
            continue;
        }
        let a = &g.atoms[i];
        if a.in_ring && !a.is_aromatic {
            continue;
        }
        if is_hetero(a.symbol.as_str()) && (n_h_of(g, i) >= 1 || a.formal_charge < 0) {
            members.insert(i);
            roots.push(i);
        }
    }

    for root in roots {
        // ブロッサム探索の起点になれるのは未マッチ (露出) の供与体分身の
        // みなので、root の分身のうち未マッチなものそれぞれから探索する。
        for &root_clone in &clones_of(root) {
            if matched[root_clone].is_some() {
                continue;
            }
            for reached_clone in graph.alternating_reachable(&matched, root_clone) {
                let reached = vertex_atom[reached_clone];
                uf.union(root, reached);
                if is_hetero(g.atoms[reached].symbol.as_str()) {
                    members.insert(reached);
                }
            }
        }
    }

    // 統合後の成分ごとにヘテロ原子端点のみを集約 (炭素は経路の中継のみ)。
    let mut by_root: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();
    for &i in &members {
        if is_hetero(g.atoms[i].symbol.as_str()) {
            by_root.entry(uf.find(i)).or_default().push(i);
        }
    }
    let mut groups: Vec<(Vec<usize>, u8)> = Vec::new();
    for (_, mut grp) in by_root {
        if grp.len() < 2 {
            continue;
        }
        grp.sort_unstable();
        let total_h: usize = grp.iter().map(|&e| n_h_of(g, e)).sum();
        let total_neg = grp
            .iter()
            .filter(|&&e| g.atoms[e].formal_charge < 0)
            .count();
        if total_h + total_neg == 0 {
            continue;
        }
        groups.push((grp, (total_h + total_neg) as u8));
    }
    groups.sort_by_key(|(m, _)| m[0]);
    groups
}

/// 可動 H 群のメンバー原子の bool マップ (番号付けの等価化に使う)。
pub(crate) fn tautomer_group_members(g: &MoleculeGraph) -> Vec<bool> {
    let n = g.atoms.len();
    let mut member = vec![false; n];
    for (eps, _) in mobile_groups(g) {
        for e in eps {
            member[e] = true;
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
    fn mobile_h_aromatic_ring() {
        // 芳香族環 (Kekule 結合次数が bond_orders では 1.5 に潰れる) でも
        // 可動 H 中心を検出できること (I12: kekule_bond_orders を使うよう修正)。
        // 2-methyl-4,5,6,7-tetrahydro-1H-benzimidazole: 環内の N=C-N(H) が
        // 芳香族認識されても中心 C の 2 ヘテロ端点として群になるべき。
        let g = build_molecule_graph("CC1=NC2=C(N1)CCCC2").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].1, 1); // 可動 H 数 = 1
        assert_eq!(groups[0].0.len(), 2); // 端点 = 環内 N 2 個
        for &e in &groups[0].0 {
            assert_eq!(g.atoms[e].symbol, "N");
        }
    }

    #[test]
    fn mobile_h_bridges_through_aromatic_ring() {
        // I13: 単一中心の星型検出では届かない、芳香環越しの多中心互変異性。
        // ピリジン置換アミド: アミド N-H が環を挟んで対側のピリジン N まで
        // 到達し、カルボニル O・アミド N・環 N の 3 端点 1 群になるべき。
        let g = build_molecule_graph("CC(=O)Nc1ccncc1").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].1, 1);
        assert_eq!(groups[0].0.len(), 3);

        // 2-アミノピリミジン: 環外アミノ N (孤立供与体) から環の両側の N まで
        // 橋渡しして 3 端点 1 群になるべき。
        let g = build_molecule_graph("CNc1ncccn1").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 3);

        // トリアゾール環 (中心原子自体が N): 3 つの環内 N が 1 群になるべき。
        let g = build_molecule_graph("C1CCc2[nH]nnc2C1").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 3);
    }

    #[test]
    fn mobile_h_does_not_over_bridge() {
        // 二級スルホンアミド N は非可動 (孤立鎖越しに他の中心と誤結合しない)。
        let g = build_molecule_graph("CNS(=O)(=O)C").unwrap();
        assert!(mobile_groups(&g).is_empty());

        // カルバミン酸: O,O 対のみが可動群になり、置換 N は含まれない。
        let g = build_molecule_graph("CNC(=O)O").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 2);
        for &e in &groups[0].0 {
            assert_eq!(g.atoms[e].symbol, "O");
        }

        // 非芳香族環 (ジヒドロピリジン部分) の N-H は、隣接する芳香族ピリジン
        // 環に単結合しているだけでは可動化しない。
        let g = build_molecule_graph("C1=Cc2ccncc2CN1").unwrap();
        assert!(mobile_groups(&g).is_empty());

        // 縮環 (アザインドール型、ピロール N-H) はピリジン型 N まで
        // 到達しない (縮環の共有原子は両環のヘテロ原子と直接隣接しない)。
        let g = build_molecule_graph("Cc1c[nH]c2cccnc12").unwrap();
        assert!(mobile_groups(&g).is_empty());
    }

    #[test]
    fn mobile_h_fused_ring_isomers_discriminate_correctly() {
        // I14: 縮合ヘテロ二環の可否は「縮環の共有原子が両環のヘテロ原子と
        // 直接隣接するか」で決まる (IUPAC 公式 InChI ソース
        // ichitaut.c/ichi_bns.c を参照して確認、実機 inchi-1 でも検証済み)。
        // トポロジー的にほぼ同一なピラゾロピリジンの位置異性体で、正しく
        // 異なる判定になることを確認する。

        // 縮環の共有原子がピラゾールの非 NH 側 N とピリジン N の両方に
        // 直接隣接 → 縮環をまたいで 3 端点 1 群になる。
        let g = build_molecule_graph("Cc1[nH]nc2ncccc12").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 3);

        // ピリジン N が「もう一方」の縮環共有原子に隣接 (ピラゾールの N とは
        // 直接隣接しない) → 縮環をまたがず、ピラゾール環内の 2 端点のみ。
        let g = build_molecule_graph("Cc1[nH]nc2cccnc12").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 2);
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
