# InChI v1 残課題 実装計画 (I19 以降)

> 対象読者: この計画だけを見て実装を進められること。前提知識は
> [RUST_INCHI_PLAN.md](RUST_INCHI_PLAN.md) の I0–I18 実装状況。
> 現状 (I18, 2026-08-04): フル InChI 一致 **96.18%** (6955/7439)、
> 正準番号付け 99.38%、InChIKey 機構は非立体で完全一致。

---

## 0. 絶対に守るルール (鉄則)

1. **退行させない**。各変更の前後で必ずコーパスゲートを回し、
   `full InChI` の exact 件数が**減らない**ことを確認する。1 件でも減ったら
   その変更は撤回するか、原因を潰してから進める。
2. **1 タスク = 1 コミット**。下の §3 の各タスクは独立している。まとめて
   やらない。1 つ実装 → 検証 → コミット → 次、を厳守。
3. **推測で直さない**。特に可動 H (互変異性) の可否は微妙で、動作中の
   ケース (アミノピリジン `CNc1ccncc1`、アクリドン、イサチン、カルボン酸)
   を壊しやすい。「こう直せば直りそう」で辺やルールを足すと高確率で退行
   する。必ず §1 の手順でコーパス全体を測ってから採否を決める。
4. **必ず回帰テストを追加**する。直した代表分子を `#[test]` で固定
   (`molrs/src/inchi/number.rs` か `stereo.rs` の `mod tests`)。
5. **コミット前に必ず**: `cargo clippy --release --all-targets -- -D warnings`
   と `cargo fmt -p molrs` (差分が出たら適用) をパスさせる。一時的な
   デバッグ用テストファイル (`tests/analyze_inchi.rs` 等) は**必ず削除**
   してからコミットする。
6. **コミットメッセージ規約** (既存の I16–I18 に倣う):
   ```
   I19: <一行要約> - full InChI 96.18% -> XX.XX%

   <本文: 何を・なぜ・根拠 (どの分子で確認したか)>

   Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
   Claude-Session: https://claude.ai/code/session_...
   ```

---

## 1. 検証ワークフロー (これを最初にセットアップする)

### 1-a. コーパスゲート (毎回これで測る)

```bash
cd molrs
# 退行検知用の閾値付きゲート (全メトリクス)
cargo test --release --test inchi_gate -- --nocapture --test-threads=1 2>&1 \
  | grep -E "full InChI|to_inchi_key:|canonical numbering|formula layer"
```

読むべき行:
```
full InChI: produced 7439, 6955 exact (96.18%); ...
canonical numbering: 7398/7453 match (99.38%)
```
`exact` の**件数** (6955) を基準にする。増えれば採用、減れば退行。

> ⚠️ ゲートの `assert!(acc >= 0.960, ...)` の閾値は、改善したら**必ず**
> 実測値に合わせて引き上げる (`molrs/tests/inchi_gate.rs` の 2 箇所:
> `full_inchi_matches_rdkit_where_produced` と
> `to_inchi_key_matches_rdkit_where_produced`)。

### 1-b. 不一致の分類テスト (どのカテゴリが残っているか見る)

`molrs/tests/analyze_inchi.rs` を**一時的に**作る (コミット前に削除):

```rust
//! 一時診断テスト: フル InChI 不一致の分類 (コミット前に削除)。
use std::io::Read;
use std::path::PathBuf;
use flate2::read::GzDecoder;
use molrs::graph::build_molecule_graph;

fn mob(s: &str) -> Vec<String> {
    let mut v = Vec::new();
    let mut rest = s;
    while let Some(i) = rest.find('(') {
        let Some(j) = rest[i..].find(')').map(|j| i + j) else { break };
        v.push(rest[i..=j].to_string());
        rest = &rest[j + 1..];
    }
    v
}

#[test]
#[ignore]
fn categorize_mismatches() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../corpus/inchi_dump.jsonl.gz");
    let mut text = String::new();
    GzDecoder::new(std::fs::File::open(&path).unwrap()).read_to_string(&mut text).unwrap();
    let mut cats: std::collections::HashMap<&str, Vec<String>> = std::collections::HashMap::new();
    for l in text.lines().filter(|l| !l.trim().is_empty()) {
        let v: serde_json::Value = serde_json::from_str(l).unwrap();
        let (smiles, want) = (v["s"].as_str().unwrap(), v["inchi"].as_str().unwrap_or(""));
        if want.is_empty() || want.contains("/i") { continue; }
        let Ok(g) = build_molecule_graph(smiles) else { continue };
        let Ok(got) = molrs::inchi::to_inchi(&g) else { continue };
        if got == want { continue; }
        let (gp, wp): (Vec<_>, Vec<_>) =
            (got.split('/').collect(), want.split('/').collect());
        let cat = if gp.len() != wp.len() {
            "layer-count-differs"
        } else {
            let mut c = "?";
            for (a, b) in gp.iter().zip(wp.iter()) {
                if a != b {
                    c = match b.chars().next() {
                        Some('c') => "c-layer",
                        Some('h') => {
                            let (gm, wm) = (mob(a), mob(b));
                            if gm.len() > wm.len() { "h-extra-mobile-group" }
                            else if gm.len() < wm.len() { "h-missing-mobile-group" }
                            else if gm != wm { "h-mobile-group-content" }
                            else { "h-fixed-h-only" }
                        }
                        Some('b') => "b-layer", Some('t') => "t-layer",
                        _ => "other",
                    };
                    break;
                }
            }
            c
        };
        cats.entry(cat).or_default().push(format!("{smiles} | got {got} | want {want}"));
    }
    let mut keys: Vec<_> = cats.iter().map(|(k, v)| (v.len(), *k)).collect();
    keys.sort_unstable(); keys.reverse();
    println!("TOTAL: {}", keys.iter().map(|(n, _)| n).sum::<usize>());
    for (n, k) in &keys {
        println!("=== {k}: {n}");
        for m in cats[k].iter().take(10) { println!("  {m}"); }
    }
}
```

実行:
```bash
cargo test --release --test analyze_inchi -- --ignored --nocapture 2>&1 \
  | grep -E "^TOTAL|^===|^  "
```

### 1-c. 実機 inchi-1 での正解確認 (可動 H の可否判断に必須)

コーパスの `want` は RDKit 経由の**公式 InChI ライブラリ**の出力なので、
それが正解。ある分子の可動 H がどうなるべきか単発で確認したいときは
コーパスを引く (`corpus/inchi_dump.jsonl.gz` を SMILES で grep) か、
実機 `inchi-1` を使う。実機のビルドと計装 (I16 で使用):

```bash
# scratchpad は session ごとに消えるので都度用意する
# IUPAC-InChI 公式ソースを gh api で取得 (curl はフックで不可)
gh api repos/IUPAC-InChI/InChI/contents/<path> --jq '.content' | base64 -d > <file>
# INCHI-1-SRC/INCHI_EXE/inchi-1/gcc/makefile で make (CMake は遅いので避ける)
# ichitaut.c / ichi_bns.c には getenv("INCHI_DEBUG_TAUT") 計装が可能:
INCHI_DEBUG_TAUT=1 inchi-1 mol.mol out.txt -AuxNone -NoLabels -Snon 2>trace
```
可動 H の可否は `ichitaut.c` の `MarkTautomerGroups` の中心点走査
(`[outer]`/`[inner]` トレース) と `bExistsAltPath` の呼び出し有無で分かる。
**「なぜ RDKit がこの群を出す/出さないか」を推測できないときは必ずこれで
確認する。**

### 1-d. 各タスクの実装サイクル (毎回これを回す)

```
1. cargo build --release                              # コンパイル
2. cargo test --release --lib inchi::                 # 単体テスト (~45 件)
3. §1-a のゲート                                       # コーパス exact 件数
4. 退行があれば §1-b で原因分子を特定 → デバッグ → 修正
5. cargo clippy --release --all-targets -- -D warnings
6. cargo fmt -p molrs
7. 回帰テスト追加 → cargo test --release -p molrs (全 217+ 件)
8. 一時ファイル削除 → コミット
```

---

## 2. 残差カテゴリ (I18 時点、合計 284 件)

| カテゴリ | 件数 | 内容 | 難易度 | 対応タスク |
|---|---|---|---|---|
| h-extra-mobile-group | 103 | 可動 H 群を余分に出す (過剰検出) | 高 | §3.4 |
| h-mobile-group-content | 76 | 群のメンバー/H 数が違う | 中 | §3.2, §3.3 |
| h-missing-mobile-group | 71 | 群が足りない/誤って 1 つに併合 | 中 | §3.1, §3.3 |
| extra-layer-other | 13 | 有機金属の多成分化 | — | 対象外 (v2) |
| formula-or-other | 7 | 有機金属の多成分化 | — | 対象外 (v2) |
| h-fixed-h-only | 6 | 固定 H の位置 (番号付けタイブレーク) | 中 | §3.5 |
| c-layer | 5 | 縮環の正準番号付け直列化 | 高 | §3.6 |

> 有機金属 (`extra-layer-other`, `formula-or-other`, 計 20 件) は
> 金属-炭素結合の切断による**多成分化**が必要で v1 対象外。触らない。

**推奨着手順**: §3.1 → §3.2 → §3.3 → §3.5 → §3.4 → §3.6
(易しく安全なものから。§3.4 と §3.6 は難しいので最後)。

---

## 3. タスク詳細

各タスクの構成: **症状 / 代表分子 / 原因の当たり / 変更箇所 / 検証 / 退行注意**。

### 3.1 リン中心の可動 H 検出 (h-missing の一部)

- **症状**: リン酸アミド等で群が全く出ない。
- **代表分子**:
  ```
  NP(=O)(N)N   got: h1-3H2 (群なし)   want: h(H6,1,2,3,4)
  ```
  (ホスホリックトリアミド。P に =O と 3 つの NH2 → 6 個の可動 H が
  P 周りの 3 N + O = 4 原子で移動する 1 群。)
- **原因の当たり**: `seed_groups` (`number.rs`) の中心原子は現状
  炭素・N のみを許可 (`center_is_c_or_n`)。P 中心が全く見られていない。
- **変更箇所**: `molrs/src/inchi/number.rs` の `seed_groups`。
  - 中心として P を許可する (S も同様の硫黄酸で既に部分対応、§3.2 と合わせて
    確認)。P 中心で「=O 受容体 + NH2/OH 供与体」の星型を作る。
  - `is_hetero` は現状 O/S/Se/Te/N。P は端点ではなく**中心**なので
    `is_hetero` は変えず、中心許可ロジック側を広げる。
- **検証**: `NP(=O)(N)N` が `(H6,1,2,3,4)` になること。コーパス全体で
  リン化合物が退行しないこと (`P` を含む want を grep して数件確認)。
- **退行注意**: ホスフィン・ホスホン酸エステル等、可動 H を持たない P
  化合物に誤って群を付けないこと。

### 3.2 チオアミド・スルファミン酸の NH2 を群に含める (h-content の一部)

- **症状**: N-H が固定扱いされ、隣の O/S の酸群に合流しない。
- **代表分子**:
  ```
  NC(=S)S       got: 2H2,(H,3,4)      want: (H3,2,3,4)
  NS(=O)(=O)O   got: 1H2,(H,2,3,4)    want: (H3,1,2,3,4)
  ```
  (ジチオカルバミン酸: NH2 の 2 H + SH の 1 H が N,S,S 3 原子で移動。
  スルファミン酸: NH2 + OH が N,O,O,O で移動。)
- **原因の当たり**: `seed_groups` の「O/S だけで酸系を成すなら N を除外」
  ルール (カルバミン酸 `CNC(=O)O` 対策) が、ここでは**過剰に**効いて
  NH2 まで落としている。カルバミン酸 (N が置換され H1 個・二級) と、
  一級 NH2 (H2 個) を区別する必要がある。
- **変更箇所**: `number.rs` の `seed_groups` 内、`os_double && os_donor`
  の酸除外分岐。除外対象の N が**一級アミン (末端・NH2)** の場合は
  除外しない、という条件を追加する。
  - ヒント: カルバミン酸 `CNC(=O)O` は N が炭素置換 (heavy_deg(N)≥2)、
    ジチオカルバミン酸 `NC(=S)S`・スルファミン酸 `NS(=O)(=O)O` は
    N が末端 NH2 (heavy_deg(N)==1)。
- **検証**: 上記 2 分子 + カルバミン酸 `CNC(=O)O` が**引き続き** O,O のみ
  (N 除外) のままであること (退行しないこと)。尿素 `NC(=O)N`・
  チオ尿素 `NC(=S)N` も確認。
- **退行注意**: `CNC(=O)O` (カルバミン酸、N 除外が正) を壊さない。

### 3.3 群の併合 vs 分割 (h-missing / h-content の主要部)

- **症状**: 独立した 2 つの可動 H 系を 1 つの大きな群に併合してしまう。
- **代表分子**:
  ```
  Oc1cc2[nH]ncc2cn1   got: (H2,7,8,9,10)   want: (H,7,10)(H,8,9)
  Oc1cc2cn[nH]c2cn1   got: (H2,7,8,9,10)   want: (H,7,10)(H,8,9)
  ```
  (ピラゾール縮環ピリジノン。ピリドン系の (H,7,10) と、ピラゾール
  N-N の (H,8,9) は**別々の**互変異性系。)
- **原因の当たり**: `mobile_groups` の union-find が、共有する炭素経由
  などで本来独立な 2 系を 1 つに繋いでいる。実 InChI は「H 供与体と
  受容体の対応がとれる極大部分系」ごとに分ける。
- **変更箇所**: `number.rs` の `mobile_groups` の union-find 統合ロジック。
  **難所**。まず §1-c の実機 inchi-1 で各群の正確なメンバーを確認し、
  「どの 2 原子が同一群か」の判定条件を掴んでから着手する。安易に
  「距離で切る」等はやらない (アクリドンの 1,4 が壊れる)。
- **検証**: 上記 2 分子が 2 群になること。**アクリドン・イサチン・
  アミノピリジン・カルボン酸が退行しないこと** (これらは 1 群が正)。
- **退行注意**: 高い。§3.1, §3.2, §3.5 を先に済ませてから、時間をかけて。

### 3.4 環外ドナーの過剰検出 (h-extra、最大カテゴリ 103 件)

- **症状**: 芳香/縮環に付いた環外の NH2・OH・N=O を可動化してしまうが、
  実 InChI は固定扱い。
- **代表分子**:
  ```
  Nc1cc2ccccc2[nH]1     got: (H3,9,10)      want: 10H,9H2    (アミノインドール)
  O=c1cc(N)c2ccccc2o1   got: (H2,10,11)     want: 10H2       (アミノクマリン)
  O=c1cc(O)c2ccccc2o1   got: (H,10,11)      want: 10H        (4-ヒドロキシクマリン)
  O=Nc1ccc[nH]1         got: (H,5,7)        want: 5H         (ニトロソピロール)
  ```
  一方、**壊してはいけない** (可動が正) 動作中ケース:
  ```
  CNc1ccncc1            (H,7,8)     アミノピリジン
  CC(=O)Nc1ccncc1       3 端点群    アミドピリジン
  O=c1c2ccccc2[nH]c2ccccc12  (H,14,15)  アクリドン
  ```
- **原因の当たり (要実機確認)**: 実 InChI が認める互変異性経路には
  「共役の型・長さ・端点元素」の条件がある。ヒドロキシクマリンは
  O-H···C=O が炭素 2 個をまたぐ (1,5 的) 経路で、これは認められない
  一方、アミノピリジンの N···環 N (1,4) やアクリドンの N-H···C=O (1,4)
  は認められる。**この境界を実機 inchi-1 の `[outer]/[inner]` 走査で
  1 パターンずつ確認**してからルール化する。
- **変更箇所**: `number.rs` の `mobile_groups` のブリッジ探索 (辺条件
  および根選択)。おそらく「環外供与体から芳香環へ入るときの受容体条件」を
  厳しくする。**ただし §3.3 と同じく高リスク**。
- **検証**: 上記過剰検出 4 分子が固定になり、かつ動作中 3 ケースが退行
  しないこと。この両立が難しいので、必ず §1-a で全体件数を見て純増を確認。
- **退行注意**: 最も高い。**最後に着手**。1 パターン直すごとにコーパス
  全体で純増か確認し、純減や横ばいなら撤回する。§3.4 は「全 103 件を
  一気に」ではなく「実機で確認できた 1 サブパターンずつ」進める。

### 3.5 末端基の固定 H 位置 (h-fixed-h-only、6 件)

- **症状**: 対称な末端で H がどちらの原子に付くか (番号付けタイブレーク)
  が実 InChI と逆。
- **代表分子**:
  ```
  NC#N          got: h3H2   want: h2H2   (シアナミド)
  C[N+]#[C-]    got: h2H3   want: h1H3   (イソシアニド)
  c1c[nH]c2cncc-2n1  got: h...9H  want: h...8H
  ```
- **原因の当たり**: `canonical_numbering` (`number.rs`) のタイブレークが
  実 InChI と一部異なる。H 数を持つ原子の順位付けの詳細。
- **変更箇所**: `number.rs` の `number_component` / `resolve` のタイブレーク
  段。慎重に (番号付けは 99.38% と高精度なので、いじると広範囲に影響)。
- **検証**: 上記 + `canonical numbering` メトリクスが 99.38% を**下回らない**
  こと。
- **退行注意**: 番号付けは c 層・h 層・立体層すべての土台。1 件でも
  numbering が減ったら即撤回。

### 3.6 縮環の正準番号付け直列化 (c-layer、5 件)

- **症状**: 対称性の高い縮環 (アセナフチレン等) で番号の割り当てが
  実 InChI と異なる。
- **代表分子**:
  ```
  C1=Cc2cccc3cccc(c23)C1
  c1cc2c3c(cccc3c1)CCC2
  ```
- **原因の当たり**: `number_component` の正準最小化 (`resolve`) の
  タイブレーク、または開始原子選択。残り 5 件と少なく、個別性が高い。
- **変更箇所**: `number.rs` の `resolve` / `edge_signature`。
- **退行注意**: §3.5 同様、番号付けの土台なので慎重に。5 件と少ないので
  優先度は低い。

---

## 4. 参考: 主要ファイルと関数の地図

| ファイル | 関数 | 役割 |
|---|---|---|
| `number.rs` | `seed_groups` | 中心原子 1 個の星型で可動 H 群の「種」を検出 |
| `number.rs` | `mobile_groups` | 種 + ブロッサム法ブリッジ探索で最終的な可動 H 群 |
| `number.rs` | `is_acceptor_from` | 端点が二重結合受容体か (芳香族は価数スラック) |
| `number.rs` | `has_search_slack` | 原子が探索辺に参加できるか (MAX_AT_FLOW 相当) |
| `number.rs` | `is_locked_zwitterion_neg` | ニトロ等の固定負電荷判定 |
| `number.rs` | `canonical_numbering` / `number_component` / `resolve` | 正準番号付け |
| `layers.rs` | `connection_layer` / `serialize_node` | c 層直列化 (枝の重み順) |
| `layers.rs` | `hydrogen_layer` / `build_components` | h 層・成分構築 |
| `stereo.rs` | `double_bond_layer` | /b 層 (E/Z + 不定 `?`) |
| `stereo.rs` | `tetrahedral_layers` | /t/m/s 層 |
| `blossom.rs` | `alternating_reachable` | 一般グラフ交互到達 (ブリッジ探索の核) |
| `tests/inchi_gate.rs` | — | コーパス退行ゲート (閾値はここ) |

コーパス: `corpus/inchi_dump.jsonl.gz` (JSONL、各行
`{"s": SMILES, "inchi": 正解 InChI, "key": 正解 Key, ...}`)。

---

## 5. スコープ外 (v1 では触らない)

- 多成分 (塩・有機金属の金属-C 切断)、同位体層 `/i`、荷電分子の
  完全対応 → v2。
- §2 の `extra-layer-other` / `formula-or-other` (有機金属 20 件) は
  多成分化が前提なので v1 では放置。
