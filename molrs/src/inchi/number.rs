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

/// 原子の負電荷が「電荷分離 (zwitterion) で固定されている」か。
///
/// ニトロ基 (`N+(=O)O-`) や芳香族 N-オキシド (`n+ - O-`) の O- は、隣接する
/// 正電荷原子と対をなす電荷分離型で、実 InChI ではこの負電荷を「除去可能な
/// 可動プロトン」として扱わない (I11 の zwitterion 中性化スキップと同じ思想)。
/// 該当する O- は可動 H 群のメンバーにはなり得る (硝酸 `O=[N+]([O-])O` は
/// 3 つの O を 1 群とし可動 H = 1 だが、その 1 は O-H 由来で O- 由来ではない)
/// が、可動 H 数には数えない。
fn is_locked_zwitterion_neg(g: &MoleculeGraph, atom: usize) -> bool {
    g.atoms[atom].formal_charge < 0
        && g.adjacency[atom]
            .iter()
            .any(|&nb| g.atoms[nb].symbol != "H" && g.atoms[nb].formal_charge > 0)
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
///
/// I19 §3.4 では、`nb` が橋頭ヘテロ原子を含む環のメンバーのときに集約判定を
/// 無効化していた (環の Kekule 構造が橋頭原子で固定されており、隣接原子の
/// 実在二重結合が問題の結合位置まで移動できるとは限らないため)。これは
/// 「再 Kekule 化の連鎖が本当に成立するか」を局所情報で近似する回避策で、
/// I22 の厳密な検証 ([`hub_allowed_endpoints`]) がその判定を正確に
/// 行えるようになったため撤去した。撤去により、橋頭 N を持つ縮環
/// (`Oc1ccn2ccnc2n1` 型) で本来より小さい群しか作れていなかった 16 件が
/// 正しくなる (99.37% → 99.58%)。
fn is_acceptor_from(
    g: &MoleculeGraph,
    kekule: &std::collections::HashMap<(usize, usize), f64>,
    center: usize,
    nb: usize,
) -> bool {
    if g.atoms[nb].is_aromatic {
        is_acceptor_agg(g, kekule, nb)
    } else {
        // 三重結合 (ニトリル等) は「二重結合を 1 本手放せる」受容体ではない
        // ので、正確に二重結合 (次数 2) のときだけ受容体とみなす。
        let key = (center.min(nb), center.max(nb));
        kekule.get(&key).copied().unwrap_or(1.0) == 2.0
    }
}

/// 原子の重原子次数 (`chem_bonds_valence` と対にして使う、H は含まない)。
fn heavy_degree(g: &MoleculeGraph, atom: usize) -> usize {
    g.adjacency[atom]
        .iter()
        .filter(|&&x| g.atoms[x].symbol != "H")
        .count()
}

/// 原子が「既存の (実 Kekule) 二重結合を最低 1 本持つ」か。
/// 実 InChI の `MAX_AT_FLOW(atom) = chem_bonds_valence − 次数` が正、に相当。
fn has_own_double(
    g: &MoleculeGraph,
    kekule: &std::collections::HashMap<(usize, usize), f64>,
    atom: usize,
) -> bool {
    chem_bonds_valence(g, kekule, atom) > heavy_degree(g, atom) as f64 + 1e-9
}

/// 原子が可動 H 探索グラフの頂点として「余裕」を持つか (I16、balanced
/// network flow の `MAX_AT_FLOW` に相当)。
///
/// 実 InChI (`ichi_bns.c`) では各原子は `MAX_AT_FLOW(atom) = chem_bonds_valence
/// − 次数` の容量を持ち、これが 0 の原子はどの結合辺にも一切マッチング
/// (二重結合の付け替え) できない — ただし互変異性の端点として認識された
/// 原子 (現在 H/負電荷を持つ供与体) は、別途 t-group hub 経由の専用辺で
/// この容量に +1 されるため実質的に参加できる。
///
/// molrs では t-group hub を作らない代わりに、この「+1」を直接ここで
/// モデル化する: 既存の二重結合を持つ原子 (素の `MAX_AT_FLOW` > 0) は常に
/// 参加可、既存の二重結合を持たなくても H/負電荷を持つヘテロ原子 (真の
/// 供与体候補) は参加可。逆に、既存の二重結合を持たず H/負電荷も持たない
/// 原子 (例: 縮環系の橋頭ヘテロ原子で価数が環内結合だけで使い切られている
/// 場合、インドリジン型の 3 配位 N) は、探索グラフのどの辺にも参加させない
/// — 実 InChI で全ての接続辺の容量が `min(MAX_AT_FLOW(i), MAX_AT_FLOW(j), 2)
/// = 0` になるのと同じ結果になる。単純なピロール型 N-H (H を持つので参加可)
/// との違いはここで区別される。
fn has_search_slack(
    g: &MoleculeGraph,
    kekule: &std::collections::HashMap<(usize, usize), f64>,
    atom: usize,
) -> bool {
    has_own_double(g, kekule, atom)
        || (is_hetero(g.atoms[atom].symbol.as_str())
            && (n_h_of(g, atom) >= 1 || g.atoms[atom].formal_charge < 0))
}

/// 環の全メンバーが探索スラックを持つ (= 環全体で結合を組み替えられる) か。
///
/// 実 InChI の `nGet15TautIn6MembAltRing` が要求する「交互環」に相当する判定。
/// 「環内の二重結合をちょうど 3 本数える」のは Kekule 化の任意性に左右される
/// (縮環では共有辺の二重結合をどちらの環が「借りる」かで変わる) ため、
/// 各メンバーが [`has_search_slack`] を持つかで判定する — ナフタレンのように
/// 環内二重結合が 2 本しかない側の環でも、全原子は隣の環との共有辺で
/// 二重結合を持つので交互環と判定される。可動 H 群の端点そのもの
/// (キノロンの環内 N-H など、二重結合を持たない代わりに H を持つ) も
/// スラックありなので環を塞がない。
///
/// 逆にクマリンのピラノン環はラクトンの O が二重結合も H も持たない
/// (孤立電子対で芳香族性に寄与しているだけ) ためここで弾かれる。
fn is_alternating_ring(
    g: &MoleculeGraph,
    kekule: &std::collections::HashMap<(usize, usize), f64>,
    ring: &[usize],
) -> bool {
    ring.iter().all(|&a| has_search_slack(g, kekule, a))
}

/// 2 つの端点 `a`, `b` が実 InChI の認める互変異性シフト位置関係にあるか。
///
/// 実 InChI (`ichitaut.c`) が標準で検出するのは次の 2 種類だけで、
/// 任意長の共役経路を互変異性とはみなさない:
///
/// - **1,3-シフト**: `a` と `b` が共通の隣接原子 (中心) を持つ。
///   単一中心の星型パターンそのもので、[`seed_groups`] が拾う範囲と重なる。
/// - **1,5-シフト**: `a-x-y-z-b` の 4 結合経路で、中間の 3 原子 `x,y,z` が
///   同一の**交互環** ([`is_alternating_ring`]) に載っているとき。
///   4-ヒドロキシピリジンの環外 O ↔ 環内 N、イサチンのケト O ↔ アミド N が
///   これ (前者は 6 員環、後者は縮環 5 員環を経由する)。
///
/// ブロッサム法の交互到達判定 ([`mobile_groups`] のブリッジ探索) は
/// 「再 Kekule 化が成立するか」という**必要条件**しか見ないため、実 InChI が
/// 互変異性とみなさない長距離・偶数長のシフトまで拾ってしまう。例:
///
/// - `O=Nc1ccc[nH]1` のニトロソ O ↔ ピロール N-H は 1,4 (5 員環の奇閉路を
///   ブロッサム収縮すると到達できてしまう)
/// - 4-ヒドロキシクマリンの環外 OH ↔ ケト O は 1,5 だが、経路の 3 原子が
///   載る環がラクトン O のせいで交互環でない
///
/// 距離だけでは分離できない (同距離で正反対の正解を持つ分子対が実在する)
/// ため、交互環という環の条件と組み合わせるのが要点。
fn taut_shift_ok(
    g: &MoleculeGraph,
    kekule: &std::collections::HashMap<(usize, usize), f64>,
    a: usize,
    b: usize,
) -> bool {
    let heavy = |i: usize| g.atoms[i].symbol != "H";
    // 1,2: 直接隣接 (ピラゾール/トリアゾール環の N-N など)。実 InChI も
    // インダゾールの 1H/2H を 1 群にする。
    if g.adjacency[a].contains(&b) {
        return true;
    }
    // 1,3: 共通の中心原子を持つ
    if g.adjacency[a]
        .iter()
        .any(|&c| heavy(c) && c != b && g.adjacency[b].contains(&c))
    {
        return true;
    }
    // 1,5: 中間 3 原子が同一の交互環に載る
    let alt: Vec<&Vec<usize>> = g
        .ring_atom_sets
        .iter()
        .filter(|r| is_alternating_ring(g, kekule, r))
        .collect();
    if alt.is_empty() {
        return false;
    }
    for &x in &g.adjacency[a] {
        if !heavy(x) || x == b {
            continue;
        }
        for &y in &g.adjacency[x] {
            if !heavy(y) || y == a || y == b {
                continue;
            }
            for &z in &g.adjacency[y] {
                if !heavy(z) || z == a || z == x || z == b {
                    continue;
                }
                if !g.adjacency[z].contains(&b) {
                    continue;
                }
                if alt
                    .iter()
                    .any(|r| r.contains(&x) && r.contains(&y) && r.contains(&z))
                {
                    return true;
                }
            }
        }
    }
    false
}

/// 容量ゼロの原子 (橋頭ヘテロ原子など、[`has_search_slack`] が false) を
/// 1 つでも含む SSSR 環に属する原子を全て「無効化」した集合を返す
/// (I19 §3.4)。
///
/// 実 InChI で確認された挙動: 橋頭ヘテロ原子 (例: インドリジン型の 3 配位
/// N、二重結合を一切持たない) を含む環は、その環自体だけでなく、その環に
/// 融合しているもう一方の環も含めて、可動 H 検出の対象から完全に除外
/// される。これは単に「橋頭原子越しの橋渡しを禁止する」以上に強い制約で、
/// 橋頭原子を全く経由しない**局所的な単一中心パターン**さえも無効化する
/// (例: `Oc1cn2cccc2cn1` — 環外 O は環内の隣接 N と直接の 1,3 パターンを
/// 成すが、その N が属する環がもう一方の環の橋頭 N と縮環しているだけで
/// 可動化されない)。コーパス実測 (RDKit 経由の実 InChI 出力) で確認した
/// 56 件以上の橋頭型縮環分子から導出した規則。
fn poisoned_ring_atoms(
    g: &MoleculeGraph,
    kekule: &std::collections::HashMap<(usize, usize), f64>,
) -> std::collections::HashSet<usize> {
    let n = g.atoms.len();
    let heavy_degree = |i: usize| {
        g.adjacency[i]
            .iter()
            .filter(|&&x| g.atoms[x].symbol != "H")
            .count()
    };
    let mut poisoned = std::collections::HashSet::new();
    for ring in &g.ring_atom_sets {
        // 非芳香族の環メンバー (縮環した飽和側の CH2 等) は、そもそも
        // Kekule/芳香族系に参加していないので二重結合を持たなくて当然
        // (単なる sp3 炭素)。橋頭「容量ゼロ」の概念は芳香族系の中でのみ
        // 意味を持つため、判定は芳香族原子に限る。
        //
        // さらに、次数 3 以上 (2 個以上の環に属する真の縮環共有原子=橋頭)
        // に限定する。次数 2 の通常の環ヘテロ原子 (フラン/チオフェン型の
        // O/S など、H を持たない単純な芳香環ヘテロ原子) も容量スラック 0
        // になるが、これは「橋頭」ではなく単に自分自身の 2 本の環内結合が
        // 常に単結合という局所的事実に過ぎない (それはブリッジ探索の辺条件
        // `has_search_slack` で個別に既に正しく除外されている) — 環全体を
        // 毒すると、その環自体が持つ正当な互変異性 (例: ベンゾオキサゾロン
        // の N-H/C=O) まで壊してしまう。次数 3 以上の橋頭だけが「その環の
        // Kekule 構造を一意に固定してしまい環全体の交互性を奪う」という
        // より強い効果を持つ (I19 §3.4、コーパス実測で確認)。
        let has_bridgehead_zero_slack = ring.iter().any(|&a| {
            a < n
                && g.atoms[a].symbol != "H"
                && g.atoms[a].is_aromatic
                && heavy_degree(a) >= 3
                && !has_search_slack(g, kekule, a)
        });
        if has_bridgehead_zero_slack {
            for &a in ring {
                if a < n && g.atoms[a].is_aromatic {
                    poisoned.insert(a);
                }
            }
        }
    }
    poisoned
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
        // 中心原子自身が価数スラック 0 (橋頭ヘテロ原子など、二重結合を
        // 一切持たない) の場合、中心にはなれない (I19 §3.4 追加修正)。
        // 中心の隣接ヘテロ原子がそれぞれ独立した実在の二重結合を持つ
        // (集約判定で「受容体」に見える) だけで、実際には中心自身がその
        // どちらとも二重結合を作れない (全ての結合が恒久的に単結合) のに
        // 誤って星型と判定してしまうケースがあった (例: `Oc1cc2cccnn2n1`
        // の橋頭 N が、両隣の独立した N=C 二重結合を持つ 2 つの N を
        // 「受容体端点」として誤って星型を形成し、本来無関係な環外 O の
        // 群に合流していた)。通常の正当な中心 (カルボン酸・アミド等) は
        // 必ず自分自身の実在二重結合を持つため、この条件で除外されない。
        if g.atoms[center].symbol == "H" || !has_search_slack(g, kekule, center) {
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
                        == 2.0
            })
            .count();
        // 中心に結合したヘテロ原子端点を分類 (受容体 = 二重結合、供与体 = H/負電荷)。
        // N が端点になるのは中心が炭素・N (トリアゾール/テトラゾールの
        // N=N-N(H)- 等)、または (末端 N かつ 中心が二重結合 O ≥2 = スルホニル
        // 級) のとき。一級スルホンアミド NH2 は可動、二級・スルフィンアミドは
        // 非可動。
        let mut endpoints: Vec<usize> = Vec::new();
        let mut has_double = false;
        // 中心自身が「単結合で繋がる端点」を最低 1 つ持つこと。中心の全端点が
        // (中心から見て) 実際に二重結合の場合 (例: N=C=N のような累積二重
        // 結合の中心) は、各端点が芳香族の集約判定で「受容体」に見えても
        // 単純な単結合↔二重結合の交換を表す星型パターンにならないため除外
        // する (スルホン酸 S(=O)(=O)OH は S からの単結合が -OH に残っている
        // ので該当しない)。
        let mut has_single_from_center = false;
        for &nb in &g.adjacency[center] {
            let sym = g.atoms[nb].symbol.as_str();
            if !is_hetero(sym) {
                continue;
            }
            // 末端 N (NH2) を非 C/N 中心の端点として許すのは、中心が
            // 高酸化状態のヘテロ原子であることを示す二重結合 O の数が
            // 十分なとき: S は 2 個 (スルホンアミド、スルフィンアミドは
            // 1 個で非可動) だが、P はリン酸トリアミド等で 1 個の P=O
            // でも一級アミドが可動になる (I19 §3.1)。
            let center_sym = g.atoms[center].symbol.as_str();
            let n_double_o_ok = n_double_o >= 2 || (center_sym == "P" && n_double_o >= 1);
            if sym == "N" && !center_is_c_or_n(center) && !(heavy_deg(nb) == 1 && n_double_o_ok) {
                continue;
            }
            let bond_from_center_is_double = kekule
                .get(&(center.min(nb), center.max(nb)))
                .copied()
                .unwrap_or(1.0)
                == 2.0;
            if !bond_from_center_is_double {
                has_single_from_center = true;
            }
            if is_acceptor_from(g, kekule, center, nb) {
                endpoints.push(nb);
                has_double = true;
            } else if n_h_of(g, nb) >= 1 || g.atoms[nb].formal_charge < 0 {
                endpoints.push(nb);
            }
        }
        if !has_double || !has_single_from_center || endpoints.len() < 2 {
            continue;
        }
        // O/S 端点だけで酸系 (二重 O/S ≥1 かつ 供与体 O/S ≥1) を成すなら、
        // N を除外して酸のみを群とする。ただし例外として、末端の一級 NH2
        // (heavy_deg == 1) は、中心が炭素で酸対が純粋に酸素だけ (カルバミン
        // 酸型、アミド性の N は非可動) でない限り除外しない — ジチオ
        // カルバミン酸 `NC(=S)S` (酸対に S を含む) やスルファミン酸
        // `NS(=O)(=O)O` (中心が S) の NH2 は酸対と一緒に可動になるが、
        // 通常のカルバミン酸 `NC(=O)O` (中心が C、酸対が O,O のみ) の NH2
        // はアミド性のまま固定される (I19 §3.2/§3.3、実 InChI で確認)。
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
        let carbamic_acid_like =
            g.atoms[center].symbol == "C" && os_ep.iter().all(|&e| g.atoms[e].symbol == "O");
        let is_primary_n =
            |e: usize| g.atoms[e].symbol == "N" && heavy_deg(e) == 1 && !carbamic_acid_like;
        let chosen: Vec<usize> = if os_double && os_donor && os_ep.len() >= 2 {
            for &e in &endpoints {
                if !os_ep.contains(&e) && !is_primary_n(e) {
                    excluded.insert(e);
                }
            }
            endpoints
                .iter()
                .copied()
                .filter(|&e| os_ep.contains(&e) || is_primary_n(e))
                .collect()
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
        //
        // 複数の中心が端点を共有してもよい (例: 縮環の 2 つの共有原子が
        // それぞれ別のヘテロ原子と直接隣接するとき、両方が独立した種を
        // 形成し、共有する原子経由で [`mobile_groups`] の union-find が
        // まとめて 1 つの群に統合する)。ここで早期に「既に使用済み」として
        // 一方を弾くと、その中心のもう一方の端点が種にも孤立供与体にも
        // ならず group から完全に脱落してしまう。
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
/// 検出済みの候補群から**偽陽性の端点を落とすフィルタ** (I22 → I24)。
/// 実 InChI (`ichi_bns.c`) と同じく、各 t-group に**仮想のハブ頂点**を置いた
/// マッチンググラフ上で「その端点が可動 H を保持できるか」を厳密に判定する。
///
/// 何を端点と認めるか (元素・価数・電荷の条件) は実 InChI 固有の規則で、
/// それを緩く取ると群を作りすぎる。そこで端点の認定は従来どおり
/// [`seed_groups`] / ブリッジ探索に任せ、その結果に対して「本当にその位置に
/// H を置いた Kekule 構造が成立するか」だけをこの関数で検証する。
///
/// # グラフの作り
///
/// スルホニル S やホスホリル P のように**二重結合を 2 本以上**持つ原子が
/// あるため、1 原子 1 マッチの単純マッチングでは Kekule 構造を表現できない。
/// 容量の分だけ原子を複製する標準的な次数制約部分グラフ → マッチングの帰着で
/// 対応する。容量は「実在の二重結合の本数 + その原子が現に持つ可動 H 数」で、
/// 後者が実 InChI の「端点はハブ経由で `MAX_AT_FLOW` が +1 される」に当たる。
///
/// 各候補群には可動 H 数 `k` だけのハブ頂点を置き、群の全端点の全複製と結ぶ。
/// ハブは**群ごとに別々**に作る (1 個の大域ハブにすると群をまたいだ H 移動を
/// 許してしまう)。どの候補群にも属さない供与体 (非互変異性の N-H など) には
/// 自分だけに繋がる専用ハブを与え、その H が動かないよう固定する。
/// この構成で初期マッチング (二重結合 + 「現に H を持つ端点 ↔ ハブ」) は
/// **完全マッチング**になり、「端点 a が H を持てる」⟺「辺 (a, ハブ) を含む
/// 最大マッチングが存在する」⟺「`G − a − hub` の最大マッチングが 1 本だけ
/// 小さい」で判定できる。
///
/// # なぜ単一の交互パス判定 (I22) では足りないか
///
/// I22 は「供与体 d を根とする M-交互パスで受容体 a に到達できるか」だけを
/// 見ていた。これは **H が 1 個だけ動く**場合には厳密に正しいが、ポリアザ
/// 縮環 (`[H]Oc1cc2c[nH]nc2nn1` 等) では**2 個の可動 H が別々の群で同時に
/// 動く**ことで初めて互変異性が成立する。このとき対称差は 2 本の独立した
/// 交互パスに分かれ、単一パス判定では永久に繋がらない。ハブを置くと 2 本が
/// ハブを経由して 1 本の交互閉路に繋がり、正しく判定できる。
/// # 判定基準は「H 数が変わりうること」
///
/// 端点が t-group のメンバーである条件は、そこに載る H の数が**現状と違う値も
/// 取れる**こと。すなわち今より 1 個少ない配置か、1 個多い配置のどちらかが
/// 妥当な Kekule 構造として成立すること。「受け取れるか」だけを見ると元々 H を
/// 持つ端点が自明に真になり、固定 O-H が誤って群に取り込まれる
/// (`[H]Oc1nc[nH]c2ccnc1-2` 型の 6-5 縮環)。逆に「手放せるか」だけを見ると、
/// 容量の都合で H を手放せない端点 (アミジン `CC(=N)N` の =NH は、相方の NH2 が
/// 既に容量一杯なので必ず H を 1 個持つ) を落としてしまう — この =NH は
/// 「H を 2 個持つ」側へは動けるので可動である。
fn hub_allowed_endpoints(g: &MoleculeGraph, cands: &[Vec<usize>]) -> Vec<Vec<bool>> {
    let n = g.atoms.len();
    let kekule = kekule_order_map(g);
    let in_pi = |i: usize| g.atoms[i].symbol != "H" && has_search_slack(g, &kekule, i);
    let bo = |u: usize, v: usize| kekule.get(&(u.min(v), u.max(v))).copied().unwrap_or(1.0);
    // 探索グラフに載る隣接原子との二重結合だけを数える (載らない相手との
    // 結合には辺が張られないので、複製を用意しても永久にマッチできない)。
    let n_double = |i: usize| {
        g.adjacency[i]
            .iter()
            .filter(|&&nb| in_pi(nb) && bo(i, nb) == 2.0)
            .count()
    };
    // その原子が現に保持している可動 H 相当の数 (H + 負電荷)。
    let donor_h = |i: usize| -> usize {
        if is_hetero(g.atoms[i].symbol.as_str()) {
            n_h_of(g, i) + usize::from(g.atoms[i].formal_charge < 0)
        } else {
            0
        }
    };

    // --- 頂点の割り当て (原子の複製 → 群ハブ → 専用ハブ) ---
    let mut n_vert = 0usize;
    let mut take = |k: usize| -> Vec<usize> {
        let v: Vec<usize> = (n_vert..n_vert + k).collect();
        n_vert += k;
        v
    };
    let mut clones: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (i, slot) in clones.iter_mut().enumerate() {
        if in_pi(i) {
            // 容量 0 (三重結合しか持たないニトリル炭素など) でも 1 頂点は置き、
            // 経路の中継点として使えるようにする。
            *slot = take((n_double(i) + donor_h(i)).max(1));
        }
    }
    let mut group_of: Vec<Option<usize>> = vec![None; n];
    for (gi, grp) in cands.iter().enumerate() {
        for &a in grp {
            group_of[a] = Some(gi);
        }
    }
    let hubs: Vec<Vec<usize>> = cands
        .iter()
        .map(|grp| take(grp.iter().filter(|&&a| in_pi(a)).map(|&a| donor_h(a)).sum()))
        .collect();
    let mut private_hubs: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (i, slot) in private_hubs.iter_mut().enumerate() {
        if in_pi(i) && group_of[i].is_none() {
            *slot = take(donor_h(i));
        }
    }
    let n_vert = n_vert;

    // --- グラフ構築 ---
    // `limit` を渡すと、その原子だけ辺を張る複製を 1 個減らす:
    // - `Fewer`: ハブに繋ぐ複製を 1 個減らす → 保持できる可動 H が今より 1 個
    //   少ない構造しか作れないグラフ
    // - `More`: 結合に繋ぐ複製を 1 個減らす → 二重結合が今より 1 本少ない
    //   = 可動 H が今より 1 個多い構造しか作れないグラフ
    //
    // 原子の複製どうしは (ハブ辺を除いて) 隣接関係が同一なので、「先頭/末尾の
    // 何個を繋ぐか」で一般性を失わずに上限を課せる。二重結合は先頭の複製から
    // 埋めるので、結合側は先頭・ハブ側は末尾を使う。初期マッチングは
    // 「実在の二重結合」と「現に H を持つ端点 ↔ ハブ」からなり、`limit` が
    // None なら完全マッチングになる。
    let build = |limit: Option<(usize, Shift)>| -> (blossom::MatchGraph, Vec<Option<usize>>) {
        let bond_clones = |i: usize| -> &[usize] {
            match limit {
                Some((a, Shift::More)) if a == i => &clones[i][..n_double(i) - 1],
                _ => &clones[i],
            }
        };
        let hub_clones = |i: usize| -> &[usize] {
            match limit {
                Some((a, Shift::Fewer)) if a == i => {
                    &clones[i][clones[i].len() - (donor_h(i) - 1)..]
                }
                _ => &clones[i],
            }
        };
        let mut graph = blossom::MatchGraph::new(n_vert);
        let mut matched: Vec<Option<usize>> = vec![None; n_vert];
        for b in &g.bonds {
            let (u, v) = (b.begin_idx, b.end_idx);
            if !in_pi(u) || !in_pi(v) {
                continue;
            }
            for &cu in bond_clones(u) {
                for &cv in bond_clones(v) {
                    graph.add_edge(cu, cv);
                }
            }
            if bo(u, v) == 2.0 {
                // 空いている複製どうしを 1 組だけ対応付ける
                let cu = bond_clones(u)
                    .iter()
                    .copied()
                    .find(|&c| matched[c].is_none());
                let cv = bond_clones(v)
                    .iter()
                    .copied()
                    .find(|&c| matched[c].is_none());
                if let (Some(cu), Some(cv)) = (cu, cv) {
                    matched[cu] = Some(cv);
                    matched[cv] = Some(cu);
                }
            }
        }
        let mut connect = |members: &[usize], hub: &[usize]| {
            let mut free = hub.iter().copied();
            for &a in members {
                let ports = hub_clones(a);
                for &ca in ports {
                    for &h in hub {
                        graph.add_edge(ca, h);
                    }
                }
                for _ in 0..donor_h(a).min(ports.len()) {
                    let (Some(h), Some(ca)) = (
                        free.next(),
                        ports.iter().copied().find(|&c| matched[c].is_none()),
                    ) else {
                        continue;
                    };
                    matched[ca] = Some(h);
                    matched[h] = Some(ca);
                }
            }
        };
        for (gi, grp) in cands.iter().enumerate() {
            connect(grp, &hubs[gi]);
        }
        for (i, hub) in private_hubs.iter().enumerate() {
            if !hub.is_empty() {
                connect(&[i], hub);
            }
        }
        (graph, matched)
    };

    let no_skip = vec![false; n_vert];
    let mu = {
        let (graph, matched) = build(None);
        graph.max_matching(&mut matched.clone(), &no_skip)
    };
    // 制限付きグラフでも同じサイズの最大マッチングが取れるか
    // (= その制限を満たす妥当な Kekule 構造 + H 配置が存在するか)。
    let feasible = |a: usize, shift: Shift| -> bool {
        let (graph, matched) = build(Some((a, shift)));
        graph.max_matching(&mut matched.clone(), &no_skip) == mu
    };

    cands
        .iter()
        .map(|grp| {
            grp.iter()
                .map(|&a| {
                    if !in_pi(a) {
                        return false;
                    }
                    (donor_h(a) >= 1 && feasible(a, Shift::Fewer))
                        || (n_double(a) >= 1 && feasible(a, Shift::More))
                })
                .collect()
        })
        .collect()
}

/// [`hub_allowed_endpoints`] の制限方向 (端点の可動 H 数を今より 1 個
/// 減らす / 増やす)。
#[derive(Clone, Copy, PartialEq, Eq)]
enum Shift {
    Fewer,
    More,
}

pub(crate) fn mobile_groups(g: &MoleculeGraph) -> Vec<(Vec<usize>, u8)> {
    let n = g.atoms.len();
    let kekule = kekule_order_map(g);
    // I19 §3.4: 橋頭ヘテロ原子を含む環 (とその融合相手の環) は、複数原子を
    // またぐブリッジ探索 (辺・孤立供与体起点) の対象から除外する
    // (poisoned_ring_atoms 参照)。単一中心の局所パターン
    // ([`seed_groups`]) 自体は橋頭原子の存在に関わらず有効だが、
    // `is_acceptor_from` の集約判定 (芳香族隣接原子が「どこかに」実在
    // 二重結合を持てば受容体とみなす) は、その隣接原子が橋頭原子を含む
    // 環のメンバーの場合は信頼できない (環の Kekule 構造が橋頭原子で
    // 一意に固定されているため、隣接原子の実在二重結合が問題の結合位置に
    // 移動できるとは限らない) ので、`seed_groups` にも渡して集約判定を
    // 局所的に無効化する。
    let poisoned = poisoned_ring_atoms(g, &kekule);
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

    // ブリッジ探索用グラフ: 両端が芳香族の辺 (環内結合)、または片方が芳香族で
    // もう片方がヘテロ原子の辺 (環とその環外ヘテロ端点、例: アミド N やケト
    // 基の O) だけを採用する。孤立アルケンやスルホニル中心越しの橋渡しを
    // 防ぐと同時に (二級スルホンアミド N やジヒドロピリジン環 N-H を誤って
    // 可動化しない)、芳香環が「無関係な環外炭素」(例: カルボン酸の C) へ
    // 逆向きに抜けてしまうのも防ぐ — 環外側がヘテロ原子でなければ、その先の
    // 全く無関係な官能基 (別のヘテロ原子群) まで誤って橋渡ししてしまうため。
    // 種メンバー (アミド N 等、非芳香族ヘテロ原子) から芳香環への入口は
    // この条件で自然に含まれる。
    let mut graph = blossom::MatchGraph::new(vertex_atom.len());
    let mut matched: Vec<Option<usize>> = vec![None; vertex_atom.len()];
    for b in &g.bonds {
        let (u, v) = (b.begin_idx, b.end_idx);
        if g.atoms[u].symbol == "H" || g.atoms[v].symbol == "H" {
            continue;
        }
        let shared = shared_rings(u, v);
        let bo = kekule.get(&(u.min(v), u.max(v))).copied().unwrap_or(1.0);
        if bo == 2.0 {
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
        // 非芳香族側が「カルボニル型受容体炭素」(ヘテロ原子への二重結合を
        // 持つ C) で、かつ結合が環内 (共有環あり) の場合も辺として許可する
        // (I18)。イサチン/フタルイミド型の縮環 5 員環で、芳香環を経由した
        // ビニロガスな経路がケト基の O まで届くために必要。環外置換基
        // (ベンゼン環上のカルボン酸の C など) は共有環を持たないため影響
        // せず、無関係な官能基の誤結合は起きない。
        let is_ring_acceptor_c = |x: usize| -> bool {
            !shared.is_empty()
                && g.atoms[x].symbol == "C"
                && g.adjacency[x].iter().any(|&nb| {
                    is_hetero(g.atoms[nb].symbol.as_str())
                        && kekule.get(&(x.min(nb), x.max(nb))).copied().unwrap_or(1.0) == 2.0
                })
        };
        let edge_ok = match (g.atoms[u].is_aromatic, g.atoms[v].is_aromatic) {
            // I17: 両端が芳香族でも、その結合が「どの環にも属さない」
            // (共有環なし) ビアリール連結結合 (別々の芳香環を単結合でつなぐ)
            // は除外する。この結合はどの Kekule 構造でも単結合のままで、
            // 二重結合の付け替え (互変異性経路) に参加できないため、環をまたぐ
            // 誤った橋渡し (例: `Nc1ccc(-c2ccncc2)cc1` のアミノ基→対側環の
            // ピリジン N) を生む。縮環 (共有原子を持つ) の結合は shared が
            // 空でないので影響しない。
            (true, true) => !shared.is_empty(),
            (true, false) => is_hetero(g.atoms[v].symbol.as_str()) || is_ring_acceptor_c(v),
            (false, true) => is_hetero(g.atoms[u].symbol.as_str()) || is_ring_acceptor_c(u),
            (false, false) => false,
        };
        if !edge_ok {
            continue;
        }
        // I16: 容量ゼロの原子 (橋頭ヘテロ原子など) はどの辺にも参加させない
        // (balanced network flow の MAX_AT_FLOW=0 相当、has_search_slack 参照)。
        if !has_search_slack(g, &kekule, u) || !has_search_slack(g, &kekule, v) {
            continue;
        }
        // I19 §3.4: 橋頭ヘテロ原子を含む環 (とその融合相手の環) に属する
        // 原子はどの辺にも参加させない (poisoned_ring_atoms 参照)。
        if poisoned.contains(&u) || poisoned.contains(&v) {
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

    // 縮環の共有原子は、実際の Kekule 二重結合が「別の」環に割り当たって
    // いることがある (どちらの環が二重結合を「借りる」かは Kekule 化の
    // 任意の選択に過ぎない)。この場合、その原子の「この環では」の分身は
    // 二重結合を持たないように見えてしまい、環内のブロッサム探索がそこで
    // 途切れる — 例えば環内の別のヘテロ原子 (トリアゾール/ピラゾール環の
    // 「中心」寄りの N 等) がブリッジで正しく拾えなくなる。対処として、
    // 同じ環内で隣接し合う「その環では二重結合を持たないが実際にはどこかに
    // 二重結合を持つ」原子どうし (縮環の共有辺の両端であることが多い) を、
    // その環内限定の仮想的なマッチとして対応付ける — 環自体は依然として
    // 閉じた奇閉路として探索でき、ブロッサム収縮で環内の全ヘテロ原子が
    // 正しく到達可能になる。実際の外側の環との接続は、共有原子のもう一方の
    // 分身 (実結合が属する側) に既に正しく設定済みなので、こちらは変更しない。
    let has_real_double = |a: usize| -> bool {
        g.adjacency[a].iter().any(|&nb| {
            g.atoms[nb].symbol != "H"
                && kekule.get(&(a.min(nb), a.max(nb))).copied().unwrap_or(1.0) == 2.0
        })
    };
    // 孤児は 2 個とは限らない (アクリドン型の中央環では、環内の全炭素の
    // 実二重結合が外側の環や環外 =O に割り当たり、4 個以上になる)。環に
    // 沿って隣接する孤児どうしを貪欲にペアリングする: まず「孤児隣接が
    // 1 個だけ」の端点から確定し (パスの端)、残りが全て次数 2 (環全体が
    // 孤児) なら最小番号から開始する。奇数長のパス等でペアにできない
    // 孤児は未マッチのまま残す (従来どおり探索の行き止まりになるだけ)。
    for (ri, ring) in g.ring_atom_sets.iter().enumerate() {
        let mut orphans: Vec<usize> = ring
            .iter()
            .copied()
            .filter(|&a| a < n && g.atoms[a].symbol != "H")
            .filter(|&a| matched[clone_in(a, ri)].is_none() && has_real_double(a))
            .collect();
        if orphans.len() < 2 {
            continue;
        }
        orphans.sort_unstable();
        let orphan_set: std::collections::HashSet<usize> = orphans.iter().copied().collect();
        let mut unpaired: std::collections::HashSet<usize> = orphan_set.clone();
        loop {
            let deg = |a: usize, unpaired: &std::collections::HashSet<usize>| {
                g.adjacency[a]
                    .iter()
                    .filter(|&&nb| orphan_set.contains(&nb) && unpaired.contains(&nb))
                    .count()
            };
            // パスの端 (孤児隣接 1) を優先、なければ環 (全て次数 2) の最小番号
            let pick = orphans
                .iter()
                .copied()
                .filter(|a| unpaired.contains(a))
                .find(|&a| deg(a, &unpaired) == 1)
                .or_else(|| {
                    orphans
                        .iter()
                        .copied()
                        .filter(|a| unpaired.contains(a))
                        .find(|&a| deg(a, &unpaired) >= 1)
                });
            let Some(a) = pick else { break };
            let Some(b) = g.adjacency[a]
                .iter()
                .copied()
                .filter(|nb| orphan_set.contains(nb) && unpaired.contains(nb))
                .min()
            else {
                break;
            };
            unpaired.remove(&a);
            unpaired.remove(&b);
            let (ca, cb) = (clone_in(a, ri), clone_in(b, ri));
            matched[ca] = Some(cb);
            matched[cb] = Some(ca);
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
    //
    // 縮環系の橋頭ヘテロ原子 (例: インドリジン型の 3 配位 N で二重結合を
    // 一切持たない) を含む環のメンバーは孤立供与体の起点にもしない
    // (poisoned_ring_atoms、I19 §3.4)。
    for i in 0..n {
        if excluded.contains(&i) || members.contains(&i) || poisoned.contains(&i) {
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
                // I19 §3.3: 縮環の共有原子 (分身を持つ) は、環ごとに独立した
                // 別々の探索から「たまたま両方とも通過する」ことがある —
                // その原子自身はヘテロ原子でない単なる経由点であっても、
                // 原子 ID ベースの union-find では同一原子として扱われて
                // しまい、本来無関係な 2 つの互変異性系 (例: 縮環の一方の
                // 環にある環外フェノール性 OH と、もう一方の環にある独立
                // した環内 N-H 互変異性) が誤って 1 群に併合されてしまう。
                // 頂点分割は探索グラフレベルでは 2 つの環を正しく分離して
                // いるので、union-find もヘテロ原子 (最終的に群メンバーに
                // なりうる原子) への到達だけを併合対象にすれば十分 — 経由
                // 点の炭素まで union すると分割の効果が最終段で失われる。
                if is_hetero(g.atoms[reached].symbol.as_str())
                    && taut_shift_ok(g, &kekule, root, reached)
                {
                    uf.union(root, reached);
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
    // I22/I24: 検出済みの群を「本当にそこで H が可動か」で検証し、動けない
    // メンバーを落とす (偽陽性フィルタ)。頂点分割のない π 全体のグラフに
    // 群ごとの t-group ハブを足し、ブロッサム法の最大マッチングで厳密に
    // 判定する ([`hub_allowed_endpoints`]) — 縮環をまたぐ長い再 Kekule 化の
    // 連鎖が実際には破綻するケース (例: `[H]Oc1c[nH]c2ccnc-2n1` の環外 OH)
    // を、局所パターンだけを見る `seed_groups` では弾けない。
    let mut cands: Vec<Vec<usize>> = by_root.into_values().collect();
    for c in cands.iter_mut() {
        c.sort_unstable();
    }
    cands.sort();
    let allowed = hub_allowed_endpoints(g, &cands);
    let mut groups: Vec<(Vec<usize>, u8)> = Vec::new();
    {
        for (gi, cand) in cands.iter().enumerate() {
            let grp: Vec<usize> = cand
                .iter()
                .copied()
                .zip(allowed[gi].iter())
                .filter(|&(_, &ok)| ok)
                .map(|(a, _)| a)
                .collect();
            if grp.len() < 2 {
                continue;
            }
            let total_h: usize = grp.iter().map(|&e| n_h_of(g, e)).sum();
            let raw_neg = grp
                .iter()
                .filter(|&&e| g.atoms[e].formal_charge < 0)
                .count();
            // H も負電荷も一切ない群は完全に偽陽性なので除外する。
            if total_h + raw_neg == 0 {
                continue;
            }
            // 可動 H 数: 電荷分離 (zwitterion) で固定された負電荷 (ニトロ・
            // N-オキシドの O-) は可動プロトンとして数えない (実 InChI 準拠、
            // is_locked_zwitterion_neg)。可動数が 0 になっても群自体は残す —
            // メンバーは正準番号付けで等価化される (ニトロの 2 つの O は
            // 対称) が、h 層には出力されない (可動 H がないため)。
            let unlocked_neg = grp
                .iter()
                .filter(|&&e| g.atoms[e].formal_charge < 0 && !is_locked_zwitterion_neg(g, e))
                .count();
            groups.push((grp, (total_h + unlocked_neg) as u8));
        }
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
/// 番号付けの比較キー: (c 層, h 層, q 層) の順。InChI が層を並べる順序と同じで、
/// 前の層で差が付かない番号付けだけが次の層で比較される。
type CandKey = (
    Vec<(usize, usize)>,
    Vec<(u8, Vec<usize>)>,
    Vec<(i8, Vec<usize>)>,
);
type Candidate = (CandKey, Vec<usize>);

/// 番号付け → h 層の比較キー。
///
/// h 層は固定 H 数の昇順にグループ化され (`1-3H,4H2,5H3`)、各グループ内は
/// 番号の昇順に並ぶ。**H を持たない原子は現れない**。番号付けは置換なので
/// 「どの H 数がいくつ現れるか」は候補間で不変であり、(H 数, 番号昇順リスト)
/// を H 数昇順に並べた列の辞書順比較が、印字される文字列の比較と一致する。
fn h_signature(atoms: &[NAtom], numbering: &[usize]) -> Vec<(u8, Vec<usize>)> {
    let mut by_k: std::collections::BTreeMap<u8, Vec<usize>> = std::collections::BTreeMap::new();
    for (i, a) in atoms.iter().enumerate() {
        if a.n_h > 0 {
            by_k.entry(a.n_h).or_default().push(numbering[i]);
        }
    }
    by_k.into_iter()
        .map(|(k, mut v)| {
            v.sort_unstable();
            (k, v)
        })
        .collect()
}

/// 番号付け → q (電荷) 層の比較キー。h 層と同じ作りで、電荷 0 の原子は現れない。
fn q_signature(atoms: &[NAtom], numbering: &[usize]) -> Vec<(i8, Vec<usize>)> {
    let mut by_q: std::collections::BTreeMap<i8, Vec<usize>> = std::collections::BTreeMap::new();
    for (i, a) in atoms.iter().enumerate() {
        if a.charge != 0 {
            by_q.entry(a.charge).or_default().push(numbering[i]);
        }
    }
    by_q.into_iter()
        .map(|(q, mut v)| {
            v.sort_unstable();
            (q, v)
        })
        .collect()
}

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
        let key = (
            edge_signature(atoms, &numbering),
            h_signature(atoms, &numbering),
            q_signature(atoms, &numbering),
        );
        if best.as_ref().map(|(k, _)| &key < k).unwrap_or(true) {
            *best = Some((key, numbering));
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
/// 2. t-group フラグを加えて精緻化
/// 3. 残る同値類は分岐し、(c 層, h 層, q 層) が辞書順最小の番号を採用
///
/// **固定 H 数と電荷は精緻化の順序キーには使わない**。これらは「同値類内で
/// どちらが小さい番号か」を per-atom に決めてしまうが、実 InChI は骨格
/// (H なし) の正準化で残った自己同型の中から **h 層・q 層の文字列を最小化
/// する**番号付けを選ぶ。両者は一致しない:
///
/// - per-atom に H 昇順とすると、シアナミド `NC#N` の 2 つの N (骨格上は
///   対称) で H を持たないニトリル N が先になり `h3H2` になってしまう。
///   実 InChI は `h2H2` (NH2 が 2 番)。
/// - 逆に per-atom に「H を持つ方が先」とすると、1-ブチン `C#CCC` で
///   CH2 が 3 番になり `h1H,3H2,2H3`。実 InChI は `h1H,4H2,2H3` —
///   末端の選択 (1↔2) と中間の選択 (3↔4) は骨格の自己同型として**連動**
///   しており、per-atom な基準では表せない。
///
/// 文字列最小化として扱えば両方とも自然に出る。電荷も同様で、層の順序が
/// c → h → q なので、メチルイソシアニド `[C-]#[N+]C` は電荷ではなく先に
/// h 層で決まり `h1H3` (CH3 が 1 番) になる。
fn number_component(atoms: &[NAtom]) -> Vec<usize> {
    // 段 1: (元素, 次数)
    let keys1: Vec<(&(u8, String), usize)> =
        atoms.iter().map(|a| (&a.elem_key, a.degree)).collect();
    let mut ranks = ranks_from_keys(&keys1);
    refine(atoms, &mut ranks);
    // 段 2: + t-group フラグ (可動 H 群のメンバーは非メンバーと区別される)
    let keys2: Vec<(usize, bool)> = atoms
        .iter()
        .enumerate()
        .map(|(i, a)| (ranks[i], a.in_tgroup))
        .collect();
    ranks = ranks_from_keys(&keys2);
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
    fn mobile_h_ring_center_atom_joins_via_virtual_pairing() {
        // 縮環の共有原子の実際の Kekule 二重結合が「別の」環に割り当たって
        // いる場合、その環では単純な交互パスが途切れてしまうため、環内
        // 限定の仮想マッチ (縮環共有原子どうしのペア付け) で環の中心寄りの
        // ヘテロ原子も正しく群に含める。SMILES の書き方 (H がどちらの N に
        // つくか) を入れ替えても対称的に正しく判定されることを確認する。
        let g = build_molecule_graph("Cc1n[nH]c2ncccc12").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 3);
    }

    #[test]
    fn mobile_h_overlapping_seeds_merge_via_shared_endpoint() {
        // 2 つの異なる中心がヘテロ端点を共有する場合、両方が独立した種を
        // 形成し union-find でまとめて 1 つの群に統合されるべき (どちらか
        // 一方だけを "used" として弾くと、共有原子の反対側の端点が
        // 群から完全に脱落する)。
        let g = build_molecule_graph("Cc1nc2cccnc2[nH]1").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 3);
    }

    #[test]
    fn mobile_h_ignores_triple_bonds_and_cumulated_centers() {
        // 三重結合 (ニトリル) は「二重結合を 1 本手放せる」受容体ではない。
        let g = build_molecule_graph("NC#N").unwrap();
        assert!(mobile_groups(&g).is_empty());

        // 累積二重結合 (N=C=N) の中心は、両端点が中心から見て共に二重結合
        // で単結合の "自由な" 端点を持たないため、単純な星型パターンとして
        // 有効ではない。
        let g = build_molecule_graph("N=C=N").unwrap();
        assert!(mobile_groups(&g).is_empty());
        let g = build_molecule_graph("N=C=O").unwrap();
        assert!(mobile_groups(&g).is_empty());
    }

    #[test]
    fn mobile_h_unrelated_functional_groups_stay_separate() {
        // 芳香環でつながっているだけの無関係な 2 つの官能基 (カルボン酸と
        // アセトアミド) は、環外側 (環に対して非芳香族) がヘテロ原子で
        // なければ環をまたいで橋渡ししない — カルボン酸の C はヘテロでは
        // ないので、環がそこへ抜け出す辺を作らない。2 つの別々の群になる
        // べき。
        let g = build_molecule_graph("OC(=O)c1ccc(NC(C)=O)cc1").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 2);
        for grp in &groups {
            assert_eq!(grp.0.len(), 2);
        }
    }

    #[test]
    fn mobile_h_zwitterion_charge_not_counted_as_mobile_proton() {
        // I17: ニトロ基 (N+(=O)O-) の O- は電荷分離で固定された負電荷であり、
        // 実 InChI では可動プロトンとして扱われない。純粋なニトロ (可動 H
        // なし) は h 層に群を出さない — ただし 2 つの O は正準番号付けの
        // 等価化のため群メンバーとしては残る (可動数 0 で返る)。
        let g = build_molecule_graph("CC[N+](=O)[O-]").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 2); // 2 つの O は等価化のためメンバー
        assert_eq!(groups[0].1, 0); // 可動 H 数 0 → h 層には出さない
        assert_eq!(
            crate::inchi::to_inchi(&g).unwrap(),
            "InChI=1S/C2H5NO2/c1-2-3(4)5/h2H2,1H3"
        );

        // 硝酸の zwitterion 形 (O=N+(O-)OH): 3 つの O は 1 群だが可動 H は
        // O-H 由来の 1 のみ (O- は数えない) → (H,2,3,4)。
        let g = build_molecule_graph("O=[N+]([O-])O").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 3);
        assert_eq!(groups[0].1, 1);
    }

    #[test]
    fn mobile_h_does_not_bridge_across_biaryl_bond() {
        // I17: 別々の芳香環を単結合 (ビアリール連結) でつなぐ分子で、一方の
        // 環の環外アミノ基が、対側の環のピリジン N まで誤って橋渡ししない
        // こと。ビアリール結合はどの Kekule 構造でも単結合のままで二重結合の
        // 付け替えに参加できないため、互変異性経路にならない。アミノ基は
        // 固定 (NH2)、群なしが正しい。
        let g = build_molecule_graph("Nc1ccc(-c2ccncc2)cc1").unwrap();
        assert!(mobile_groups(&g).is_empty());

        // 対照: アミノ基が同一のピリジン環に直接ついている場合は従来どおり
        // 橋渡しする (ビアリール除外は縮環・同一環の橋渡しに影響しない)。
        let g = build_molecule_graph("CNc1ccncc1").unwrap();
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

    #[test]
    fn mobile_h_phosphoric_triamide() {
        // I19 §3.1: リン酸トリアミド (P(=O)(NH2)3) は P=O と 3 つの NH2 が
        // 1 群になり、可動 H = 6 (各 NH2 の 2H)。硫黄と異なり P は P=O が
        // 1 個でも一級アミドが可動になる (スルホンアミドは S=O が 2 個
        // 必要、スルフィンアミドは非可動のまま — 既存の
        // `mobile_h_does_not_over_bridge` の二級スルホンアミドケースは
        // 変わらないはず)。
        let g = build_molecule_graph("NP(=O)(N)N").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 4); // O + 3 N
        assert_eq!(groups[0].1, 6);

        // 炭素置換の膦酸アミド (ホスフィン酸アミド、非可動な炭素置換のみ)
        // は群を作らない。
        let g = build_molecule_graph("CP(=O)(C)C").unwrap();
        assert!(mobile_groups(&g).is_empty());
    }

    #[test]
    fn mobile_h_primary_amine_joins_os_acid_group() {
        // I19 §3.2: O/S だけの酸系対に対する N 除外規則
        // (カルバミン酸 `CNC(=O)O` のような**置換された**二級 N を除外する)
        // は、末端の一級 NH2 (中心以外に重原子隣接を持たない) までは除外
        // しない。ジチオカルバミン酸・スルファミン酸は NH2 も酸対と一緒に
        // 可動になるべき (実 InChI で確認)。
        let g = build_molecule_graph("NC(=S)S").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 3); // N + S + S
        assert_eq!(groups[0].1, 3); // NH2 (2) + SH (1)

        let g = build_molecule_graph("NS(=O)(=O)O").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 4); // N + O + O + O
        assert_eq!(groups[0].1, 3); // NH2 (2) + OH (1)

        // 対照: カルバミン酸の N は甲基置換 (二級) なので引き続き除外。
        let g = build_molecule_graph("CNC(=O)O").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 2);
        for &e in &groups[0].0 {
            assert_eq!(g.atoms[e].symbol, "O");
        }

        // 対照: 無置換のカルバミン酸 (NH2 は一級) も、中心が炭素で酸対が
        // 純粋に酸素だけ (アミド性) なので N は除外されたまま。
        // ジチオカルバミン酸 (酸対に S) やスルファミン酸 (中心が S) との
        // 違いはここ。
        let g = build_molecule_graph("NC(=O)O").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 2);
        for &e in &groups[0].0 {
            assert_eq!(g.atoms[e].symbol, "O");
        }
    }

    #[test]
    fn mobile_h_independent_ring_systems_stay_separate() {
        // I19 §3.3: 縮環の共有原子 (分身を持つ) を、環ごとに独立した別々の
        // 探索が「たまたま両方とも通過する」ことがある。その原子自身は
        // ヘテロ原子でない単なる経由点でも、原子 ID ベースで union すると
        // 無関係な 2 つの互変異性系が誤って 1 群に併合されてしまっていた
        // (縮環ピラゾロピリジノール: 環外フェノール性 OH ↔ 環内ピリジン N
        // の系と、ピラゾール環内 N-N-H の系は本来無関係)。到達原子がヘテロ
        // のときだけ union するよう修正し、正しく 2 群に分離することを
        // 確認する (実 InChI で検証済み)。
        let g = build_molecule_graph("Oc1cc2[nH]ncc2cn1").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 2);
        let mut sizes: Vec<usize> = groups.iter().map(|(m, _)| m.len()).collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![2, 2]);

        // 対照: 正当な縮環越しブリッジ (共有原子が両環のヘテロ端点に直接
        // 隣接する場合) は引き続き 1 群に統合される (I16 で検証済みの規則、
        // 退行していないこと)。
        let g = build_molecule_graph("Cc1[nH]nc2ncccc12").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 3);
    }

    #[test]
    fn mobile_h_acridone_vinylogous_amide() {
        // I18: アクリドン (中央 6 員環に C=O と N-H が 1,4)。中央環の全炭素の
        // 実二重結合は外側のベンゾ環や環外 =O に割り当たるため孤児が 4 個
        // 以上になる。環に沿った孤児の貪欲ペアリングで中央環を交互閉路と
        // して探索でき、N-H から C=O の O まで届いて 1 群 (N, O) になる。
        let g = build_molecule_graph("O=c1c2ccccc2[nH]c2ccccc12").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 2);
    }

    #[test]
    fn mobile_h_isatin_reaches_both_carbonyls() {
        // I18: イサチン (非芳香族 5 員環に 2 つの C=O と N-H)。芳香族を経由
        // する経路で環内カルボニル炭素 (ヘテロへの二重結合を持つ受容体 C)
        // まで辺を延ばし、両方の O とアミド N が 1 群 (3 端点) になる。
        let g = build_molecule_graph("O=C1Nc2ccccc2C1=O").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 3);
    }

    /// I22/I24: 検出済みの群を「本当にそこで H が可動か」で厳密検証し、
    /// 動けないメンバーを落とす ([`hub_allowed_endpoints`])。縮環をまたぐ
    /// 長い再 Kekule 化の連鎖が実際には破綻するケースを、局所パターン
    /// だけの seed_groups では弾けなかった。
    #[test]
    fn mobile_h_exact_filter_drops_unreachable_endpoints() {
        // 6-5 縮環アザインドール型 + 環外 OH。局所的には OH-C=N の 1,3 パターンに
        // 見えるが、必要な再 Kekule 化の連鎖が N-H の価数で破綻するため、実 InChI は
        // 可動 H 群を作らない (OH も N-H も固定)。
        let g = build_molecule_graph("[H]Oc1c[nH]c2ccnc-2n1").unwrap();
        assert!(mobile_groups(&g).is_empty());
        let g = build_molecule_graph("[H]Sc1ncc2ccc[nH]c1-2").unwrap();
        assert!(mobile_groups(&g).is_empty());
    }

    /// I24: 2 個の可動 H が**別々の群で同時に**動くことで初めて成立する
    /// 互変異性 (ポリアザ縮環)。単一の交互パス判定 (I22) では、対称差が
    /// 2 本の独立したパスに分かれるため永久に繋がらない。t-group ハブを
    /// 置くと 2 本がハブ経由で 1 本の交互閉路になり正しく判定できる。
    #[test]
    fn mobile_h_two_groups_shift_simultaneously() {
        // `[H]Oc1cc2c[nH]nc2nn1`: 環外 OH と環内 N が 1 群 (H,7,10)、
        // ピラゾール側の 3 つの N がもう 1 群 (H,6,8,9) になる。I22 までは
        // 前者が作れず OH が固定されていた。
        let g = build_molecule_graph("[H]Oc1cc2c[nH]nc2nn1").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups, vec![(vec![0, 9], 1), (vec![5, 6, 8], 1)]);

        // 対照 (同じ 6-5 縮環でも縮環結合が単結合で書かれる型): こちらは
        // 環外 OH が固定されたままでなければならない。局所パターンは同一で、
        // 大域的な Kekule 充足可能性でも区別できない (I24 の棄却案 1)。
        let g = build_molecule_graph("[H]Oc1nc[nH]c2ccnc1-2").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert!(!groups[0].0.contains(&0), "環外 O(0) は群に入らない");
    }

    /// I24: 端点の可動性は「H 数が今と違う値も取れること」で判定する。
    /// 「H を手放せるか」だけを見るとアミジンの =NH を落としてしまう —
    /// 相方の NH2 が容量一杯なので必ず H を 1 個持つが、「H を 2 個持つ」側
    /// (`CC(N)=N`) へは動けるので可動である。
    #[test]
    fn mobile_h_amidine_endpoint_that_can_only_gain() {
        let g = build_molecule_graph("CC(=N)N").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 2);
        assert_eq!(groups[0].1, 3, "可動 H は 3 個");
    }

    /// 厳密検証は、二重結合を 2 本持つスルホニル S を容量 2 の頂点複製で
    /// 扱う。単純マッチング (1 原子 1 マッチ) だと S=O を 1 本しか表現できず、
    /// スルホン酸の 3 つの O が 1 群にならなくなる回帰があった。
    #[test]
    fn mobile_h_sulfonyl_capacity_two_keeps_all_oxygens() {
        let g = build_molecule_graph("CS(=O)(=O)O").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 3, "3 つの O が 1 群");
        let g = build_molecule_graph("C=CS(=O)(=O)N").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 3, "2 つの O と N が 1 群");
    }

    #[test]
    fn mobile_h_bridgehead_ring_blocks_bridge_search_only() {
        // I19 §3.4: 縮環の橋頭ヘテロ原子 (二重結合を一切持たない、次数 3
        // 以上、例: インドリジン型の 3 配位 N) を含む環は、複数原子をまたぐ
        // ブリッジ探索 (辺・孤立供与体起点) の対象から除外される
        // (poisoned_ring_atoms)。実 InChI で確認: このケースは環外フェノール
        // 性 OH がどこにも橋渡ちせず完全に固定される。
        let g = build_molecule_graph("Oc1ccnc2cccn12").unwrap();
        assert!(mobile_groups(&g).is_empty());

        // 対照 1: 橋頭原子を含む環に属していても、単一中心の局所パターン
        // (実在する二重結合を使う直接の O-C=N 型) はブリッジ探索を必要と
        // しないため、橋頭原子の有無に関わらず常に有効であるべき — ここでは
        // O が結合する炭素が環内窒素と実際に二重結合しており (橋頭原子とは
        // 別の窒素)、局所的に 2 端点の群が正しく成立する。
        let g = build_molecule_graph("Oc1cc2ccccn2n1").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 2);

        // 対照 2: 橋頭原子ではない通常の次数 2 環ヘテロ原子 (フラン型の O、
        // H を持たず二重結合もないので価数スラック 0 になるが橋頭ではない)
        // を含む環は毒されない — ベンゾオキサゾロンの N-H/C=O 互変異性が
        // 正しく機能する。
        let g = build_molecule_graph("O=c1[nH]c2ccccc2o1").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 2);

        // 対照 3: 縮環した飽和 (非芳香族) 側に価数スラック 0 の sp3 炭素
        // (単なる CH2、橋頭とは無関係) があっても毒されない — トリアゾール
        // 環の 3 つの N が正しく 1 群になる (I9 の既存テストと同じ分子)。
        let g = build_molecule_graph("C1CCc2[nH]nnc2C1").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 3);
    }

    #[test]
    fn mobile_h_bridgehead_atom_cannot_be_seed_center() {
        // I19 §3.4 追加修正: 橋頭ヘテロ原子自身 (価数スラック 0、全ての
        // 結合が恒久的に単結合) が `seed_groups` の中心として使われ、その
        // 両隣の独立した (橋頭原子とは無関係な実在二重結合を持つ) ヘテロ
        // 原子 2 つを集約判定で「受容体端点」とみなし、誤った星型を形成
        // する抜け穴があった。橋頭原子はどの結合も二重結合にできない
        // (全て恒久的に単結合) ため、中心にはなり得ない。
        //
        // `Oc1cc2cccnn2n1` (トリアゾロピリジン型): 環外フェノール性 OH は
        // 隣接する N と局所的に正しく 2 端点の群を成すべきだが、橋頭 N が
        // 中心となって両隣の N を誤って拾い、無関係な 3 端点の群に併合
        // されていた。
        let g = build_molecule_graph("Oc1cc2cccnn2n1").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 2);

        // 対照: 橋頭原子を持たない通常の縮環系 (I16 で検証済みの規則、
        // 共有原子が両環のヘテロ端点に直接隣接) は退行しない。
        let g = build_molecule_graph("Cc1[nH]nc2ncccc12").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 3);
    }

    #[test]
    fn mobile_h_aggregate_acceptor_disabled_in_bridgehead_ring() {
        // `Oc1cn2cccc2cn1`: 環外フェノール性 OH が結合する炭素の実際の
        // 二重結合は別の隣接炭素にあり、環内のもう一方の隣接 N は独立した
        // 実在二重結合を持つため `is_acceptor_agg` の集約判定では「受容体」に
        // 見える。しかし橋頭 N を含む環では Kekule 構造が固定されており、
        // その二重結合を当該結合位置まで移動させる再 Kekule 化が成立しない
        // ため、群は一切形成されないのが正しい。
        //
        // I19 §3.4 ではこれを「橋頭を含む環のメンバーなら集約判定を無効化」
        // という局所的な近似で実現していたが、I22/I24 の厳密な検証
        // (hub_allowed_endpoints) が同じ結論をより正確に出せるため近似は
        // 撤去した。このテストは撤去後も結論が変わらないことを固定する。
        let g = build_molecule_graph("Oc1cn2cccc2cn1").unwrap();
        assert!(mobile_groups(&g).is_empty());

        // 対照: 橋頭原子を持たない通常の柔軟な芳香環では、集約判定が
        // 引き続き機能する (アミノピリジンの橋渡し、退行しないこと)。
        let g = build_molecule_graph("CNc1ccncc1").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 2);
    }

    #[test]
    fn mobile_h_fused_ring_bridges_when_fusion_atom_touches_both_heteroatoms() {
        // I16: 実機 inchi-1 で検証済みの規則 — 縮環越しの互変異性ブリッジは、
        // 縮環の共有原子 (どちらか一方) が両方の環でそれぞれヘテロ端点に
        // 直接隣接しているときに限って成立する。この分子は共有原子が
        // ピロール型 N-H (環1) とピリジン型 N (環2) の両方に直接隣接する
        // ため、アザインドール型 (mobile_h_does_not_over_bridge 参照、共有
        // 原子が両側のヘテロ原子に同時には隣接しない) とは異なりブリッジが
        // 成立するべき。この規則は seed_groups の「中心原子 1 個 + 直接
        // ヘテロ端点」検出でそのまま自然に捉えられる (環をまたぐ専用の
        // 探索は不要)。
        let g = build_molecule_graph("Cc1c[nH]c2ncccc12").unwrap();
        let groups = mobile_groups(&g);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.len(), 2);
        for &e in &groups[0].0 {
            assert_eq!(g.atoms[e].symbol, "N");
        }
    }
}
