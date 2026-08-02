# SMILES → InChI / InChIKey (pure Rust) 実装計画

> **v1 実装状況 (2026-08-02 更新)**: I0–I13 実装済み。実測 (コーパス 7,453 分子、
> RDKit オラクル):
> - 式層 (Hill): 99.13% 一致 (残差は電荷正規化・有機金属切断)
> - 正準番号 (AuxInfo /N:): **99.29%** 一致
> - **InChIKey 機構** (`inchi_key_from_string` vs RDKit `InchiToInchiKey`,
>   非立体 InChI 7,245 件): **7245/7245 = 完全一致**。SHA-256 (FIPS 180-4) +
>   base-26 (`ikey_base26.c` 移植、除外トリプレット対応) は公式とビット一致
> - フル InChI 文字列・InChIKey (v1 適用範囲 = 中性・単一成分・立体/同位体
>   なし): **89.11%** が RDKit と完全一致 (I4 時点は 74.17%、詳細な推移は下記)
> - CLI: `cargo run --bin inchi` (stdin SMILES → JSONL)
> - **立体対応済み (I6-I8, 2026-07-19)**: `/b` (E/Z 二重結合)・`/t/m/s`
>   (四面体) 層を実装。四面体パリティは molrs の R/S を再利用し
>   `raw = '-' iff (rs_bit XOR perm(CIP昇順→正準昇順))`、/m は最初の中心を
>   '-' に正規化。InChIKey の minor 文字列は立体層 (先頭 '/' 込み) の
>   **2 重連結** を SHA-256 (公式 ikey.c 準拠)。立体分子 208 件中 151 件が
>   InChI 文字列・キー完全一致 (残りは骨格側の可動 H/電荷が原因)。
> - **可動 H 認識を一般化 (I9, 2026-07-23)**: 中心原子に結合したヘテロ端点
>   (受容体 = 二重結合、供与体 = H/負電荷) で 1,3-互変異性群を作る一般規則に。
>   端点は末端でなくてよい (アミド/ラクタムの二級 N)。N 端点は中心が炭素
>   または (末端 N かつ 中心が二重結合 O ≥2 = スルホニル級) のとき。O/S だけで
>   酸系を成すなら N を除外 (カルバミン酸)。アミド/アミジン/グアニジン/尿素/
>   ラクタム/一級スルホンアミドを新たに正しく処理。
> - **電荷 q/p + c/h 層修正 (I10, 2026-07-23)**: normalize.rs で電荷正規化
>   (負の酸点・ハロゲン化物イオンをプロトン付加で中性化、陽イオン塩基を
>   脱プロトン、残余 → /q、移動 → /p)。c 層の分岐は全非末端項目を 1 つの
>   カンマ括弧に統合する規則へ修正 (四級 N の (6-2,7-3)8-4)。h 層の可動群は
>   カンマなし連結 ((H,5,6)(H,7,8))。単一成分の荷電分子 53/53 完全一致。
>   フル InChI 一致 **79.04%** (5880/7439)。
> - **zwitterion/イリド中性化のスキップ (I11, 2026-07-23)**: 電荷分離した
>   隣接原子 (N-オキシド・ニトロ・アジド等) はプロトン化・脱プロトン化の
>   対象から除外 (InChI は共有結合の中性形のまま扱う)。フル InChI 一致
>   **79.6%**。
> - **芳香環の可動 H 検出修正 (I12, 2026-07-26)**: `mobile_groups` (I9) が
>   結合次数の判定に `g.bond_orders` (芳香族は 1.5 に潰れる) を使っていたため、
>   芳香族として認識される環 (イミダゾール/トリアゾール/テトラゾール/
>   ピリジノン様の縮合環など) では中心原子の二重結合受容体を一切検出できて
>   いなかった。`g.kekule_bond_orders` (Kekule 化済みの実値 1/2) を使うよう
>   修正。単一中心 (中心 1 原子 + ヘテロ端点 2 個) パターンの環内版が新たに
>   拾えるようになり、フル InChI 一致 **79.6% → 84.50%** (6286/7439)。
> - **芳香環をまたぐ多中心互変異性検出 (I13, 2026-08-02)**: `mobile_groups`
>   を「種 (中心 1 原子の星型検出、旧 I9) → ブリッジ統合」の 2 段に拡張。
>   Kekule 構造を「マッチング」とみなし、種の各メンバー・孤立供与体
>   (どの中心の候補にもならなかった H 保持ヘテロ原子、例: トリアゾール環の
>   中心 N) から「自由結合 → 二重結合 (マッチング辺) で反転」を繰り返す
>   交互パス探索 (union-find で統合) を実装。誤爆抑制のため
>   (a) 中継原子は芳香族限定 (孤立アルケン・スルホニル中心越しの誤統合を防止)、
>   (b) 環内に入ったら単一の SSSR 環にロックし、縮環の共有原子越しに
>   別の環へ「乗り換える」ことを禁止 (インダゾール/アザインドール型で
>   ピロール N-H とピリジン N を誤って同一群にしない)、
>   (c) 酸除外規則 (カルバミン酸で N を落とす等) で意図的に除外された原子は
>   孤立供与体として再導入しない、の 3 点を実装。
>   ピリジン置換アミド・2-アミノピリミジン・キナゾリンジオン・
>   トリアゾール/テトラゾール環等が新たに正しく群化。
>   フル InChI 一致 **84.50% → 89.11%** (6629/7439)。
> - **I14 調査 (2026-08-02、実装は見送り)**: 残る主因である縮合ヘテロ二環
>   (インダゾール/アザインドール系) の可否判定を、IUPAC 公式 InChI ソース
>   (`github.com/IUPAC-InChI/InChI`, v1.07, `INCHI-1-SRC/INCHI_BASE/src/`)
>   を直接参照して解明・移植を試みた。
>   - **判明した実装原理**: 実際の可動 H 判定は単純なグラフ探索ではなく
>     **フローネットワークの実行可能性判定**。`ichitaut.c` の
>     `nGetEndpointInfo()` が原子ごとに供与体/受容体を「価数のスラック」
>     (中性価数と実結合次数和の差) で分類し (記号による特別扱いではなく
>     純粋に算術)、`RegisterEndPoints`/`FindAccessibleEndPoints` が
>     `ichi_bns.c` の **`bExistsAnyAltPath` → `bExistsAltPath` →
>     `RunBalancedNetworkSearch`** を呼ぶ。これは `CreateTGroupInBnStruct`
>     が構築する **単一の架空「t-group ハブ頂点」** (分子内の全適格ヘテロ
>     原子が容量 `neutral_valence - num_bonds`、流量 `min(num_H, cap)` で
>     直結) を含む **balanced network (最大流型) 構造**上の s-t 到達可能性
>     判定であり、2 原子間の直接パス探索ではなく分子全体のフロー整合性を
>     一括で問う問題になっている。`nGet15TautIn6MembAltRing` 等の環サイズ
>     別関数群は候補対を高速に絞る事前フィルタに過ぎず、最終的な可否は
>     常にこのフロー判定に帰着する。
>   - **実装を試みて判明した限界**: このフロー判定の本質 (前方: 自由結合で
>     新たに二重結合を作る/後方: 既存の二重結合を手放す、の 2 状態交互探索)
>     を単純化し、I13 の単一環ロック機構を撤去してこの 2 状態探索に置換
>     したところ、コーパス全体でフル InChI 一致率が **89.11% → 80.76%**
>     まで悪化した (過剰な群化: 本来分離すべき分子まで誤って結合)。単純な
>     2 頂点間到達可能性では、ハブ頂点の容量制約 (分子全体で同時に
>     整合するフロー配分が必要という制約) を反映できておらず、局所的には
>     「パスが存在する」ように見えても実際には無効な組み合わせを大量に
>     受理してしまうことが原因と判断。この変更は**リバート済み** (I13 の
>     単一環ロック・芳香族中継限定ヒューリスティックのままコミット
>     `69f8afb` を維持)。
> - v1 の非対応 (→ `Unsupported`): 多成分 (塩)・同位体 (i)・有機金属切断。
>   残る文字列不一致の主因は**縮合した 2 芳香ヘテロ環にまたがる互変異性**
>   (インダゾール/アザインドール系)。正確な再現には `CreateTGroupInBnStruct`
>   相当の容量付きバランスドネットワーク (ハブ頂点 + 各ヘテロ原子の価数
>   スラック容量 + 分子全体での流量整合性) を `ichi_bns.c` の
>   `RunBalancedNetworkSearch` を参照しつつ正しく移植する必要がある —
>   単純な 2 状態到達可能性では代替不可 (逐次拡張、I14 候補、要再挑戦)。


関連文書: [RUST_PORT_HANDOFF.md](RUST_PORT_HANDOFF.md) (開発サイクル)、
[RUST_3D_PLAN.md](RUST_3D_PLAN.md)・[RUST_2D_PLAN.md](RUST_2D_PLAN.md) (完了済み)。

## 目的

`MoleculeGraph` から **標準 InChI (InChI=1S/…)** と **InChIKey** を生成する。
IUPAC 公式 InChI (C 実装) と**ビット完全一致**させる — InChIKey は InChI 文字列の
SHA-256 ハッシュなので、公式キー (PubChem 等) と一致させるには InChI 文字列を
公式と完全一致させることが前提。molrs は依存クレートゼロのため **SHA-256 を
自前実装** する (FIPS 180-4)。

**v1 の対象**: 骨格層 — 式 (Hill)・接続 `c`・水素 `h` (共通 O/N 系の可動 H 群
`(H,n,m)` を含む)・電荷 `q`/`p`・InChIKey。立体 (`b`/`t`/`m`/`s`)・同位体 (`i`)・
一般の互変異性正規化は v2 以降。

## 差分検証の武器 (RDKit venv)

- `Chem.MolToInchi` — フル InChI 文字列
- `Chem.MolToInchiKey` — InChIKey
- `inchi.MolToInchiAndAuxInfo` の AuxInfo `/N:` — **公式の正準番号付け**。
  最難関の番号を文字列生成と切り離して単独ゲートできる
- `inchi.InchiToInchiKey` — 任意 InChI 文字列 → キー。SHA-256/base-26 の
  キー機構を InChI 生成と独立に検証できる

## モジュール構成 (`molrs::inchi`)

```
molrs/src/inchi/
  mod.rs        # 公開 API・InchiError・統括
  sha256.rs     # FIPS 180-4 SHA-256
  base26.rs     # InChIKey base-26 + キーレイアウト
  formula.rs    # Hill 式・多成分
  number.rs     # InChI 正準番号付け
  normalize.rs  # 電荷 q/p・可動 H 群認識
  layers.rs     # c/h/q/p 層直列化
```

公開 API:
```rust
pub fn to_inchi(g: &MoleculeGraph) -> Result<String, InchiError>;
pub fn to_inchi_key(g: &MoleculeGraph) -> Result<String, InchiError>;
pub fn inchi_key_from_string(inchi: &str) -> String;
pub fn inchi_of(smiles: &str) -> Result<String, InchiError>;
pub fn inchi_key_of(smiles: &str) -> Result<String, InchiError>;
```

## ステップ (I0–I5)

- **I0**: 本計画保存
- **I1**: sha256 (NIST ベクタ) + base26 + InChIKey レイアウト
  (`AAAAAAAAAAAAAA-BBBBBBBBSA-N`: 骨格 14 + 副層 8 + フラグ `SA` + プロトン化 1)。
  ゲート: `inchi_key_from_string` == RDKit `InchiToInchiKey`
- **I2**: Hill 式。ゲート: 式層一致
- **I3**: 正準番号付け (初期色 → Morgan 精緻化 → 正準最小化)。最難関。
  ゲート: AuxInfo `/N:` 一致
- **I4**: c/h/q/p 層 + `to_inchi`。ゲート: フル文字列一致 (v1 範囲)、全体一致率
- **I5**: `to_inchi_key` 結線 + CLI + `tests/inchi_gate.rs`。文書化

## リスク

| リスク | 対策 |
|---|---|
| 正準番号の公式完全一致 (最大) | AuxInfo `/N:` を単独ゲート化、文字列前に確定 |
| 可動 H (酢酸すら `(H,3,4)`) | v1 は共通 O/N 酸系に限定、逸脱は適用外として追跡 |
| コーパス全数一致は非現実的 | v1 範囲で 100%、全体は一致率報告 (単調増加) |
| SHA-256/base-26 実装ミス | NIST ベクタ + `InchiToInchiKey` 二重オラクル |

## 検証

フィクスチャは smiles2iupac `tools/` (RDKit venv) で AuxInfo `/N:`・フル InChI・
InChIKey をコーパスから採取し `MOLRS_ROOT/corpus/` へ。`tests/inchi_gate.rs` が
(a) 番号一致 (b) v1 範囲のフル文字列一致 (c) キー一致 (d) 全体一致率 を検査。
キー機構は `inchi_key_from_string` == `InchiToInchiKey` で恒久検証。
