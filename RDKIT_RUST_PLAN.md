# RDKit 機能の Rust 実装計画 (molrs 汎用化)

`molrs` は smiles2iupac が必要とする RDKit 機能のサブセットを、
**7,453 分子の全数一致ゲート**付きで Rust 実装したものである (RUST_PORT_PLAN.md Phase 1)。
本計画はこれを一般用途のケモインフォマティクスライブラリへ拡張するためのロードマップ。
各ステップは Sonnet クラスが 1 セッションで完了できる粒度に分割してある。

## 0. 現在の資産 (実装済み・検証済み)

| 機能 | RDKit 対応物 | 検証 |
|---|---|---|
| SMILES パーサ (OpenSMILES コア) | `MolFromSmiles` | コーパス 7,453 件 |
| 原子価モデル・暗黙 H・AddHs/RemoveHs 相当 | `updatePropertyCache` 等 | RDKit ダンプ全数一致 |
| ケクレ化 + 芳香族認識 (デフォルトモデル) | `Kekulize`/`setAromaticity` | 同上 |
| 対称化 SSSR (環順・原子順まで一致) | `GetSymmSSSR`/`GetRingInfo` | 同上 |
| 正規 SMILES (RDKit 非互換・健全) | `MolToSmiles` | 冪等 + 異表記 498 群 |
| フラグメント SMILES / 成分分解 | `MolFragmentToSmiles`/`GetMolFrags` | ユニットテスト |
| VF2 部分構造 (mol / SMILES-as-SMARTS) | `GetSubstructMatches` | フィクスチャ 2,995 対 |
| レガシー CIP (R/S, E/Z) | `AssignStereochemistry` | RDKit ダンプ全数一致 |

**検証手法も資産である**: 「RDKit を Python から叩いて正解フィクスチャを生成し、
Rust 統合テストで全数比較する」パイプライン (tools/*.py + corpus/*.gz) は
以降の全ステップでそのまま使える。

## 1. 方針

1. **差分テスト駆動**: 新機能は必ず RDKit フィクスチャとの機械比較ゲートを先に作る。
   互換を要求する機能 (記述子・フィンガープリント等) は値の完全一致、
   互換不要の機能 (正規 SMILES 等) は性質 (冪等・不変性) をゲートにする
2. **クレート分離**: smiles2iupac への影響を避けるため、汎用化は
   `molrs` の API を保ったまま進める。破壊的変更は smiles2iupac の
   全ゲート再実行とセットでのみ許す
3. スコープ外 (当面): 量子化学、描画の美観、Python バインディング。
   **3D 配座生成と UFF 力場は独立計画 [RUST_3D_PLAN.md](RUST_3D_PLAN.md) を参照**
   (距離幾何法 + UFF、C1〜C10 の 10 ステップ)

## 2. フェーズ構成

### Phase R1: ライブラリとしての土台固め

#### R1.1 スケール制約の除去
現在 `u128` ビットマスク前提で **128 原子上限** (rings.rs, canon.rs ほか)。
可変長ビットセット (自作 64bit ブロック or `fixedbitset`) に置換し、
上限を撤廃。巨大分子 (ポリマー等) のパース・環認識のプロパティテストを追加。

#### R1.2 正規 SMILES への立体出力
`@`/`@@` と `/`,`\` の書き出し。書き手側のパリティ計算は
stereo.rs の `AdjustAtomChiralityFlags` 逆変換として実装
(読み手の規約は実装済みなので「自分の読み手で読み戻して一致」が完全なゲートになる)。
完了条件: コーパス立体分子でラウンドトリップ (canon → parse → CIP 再計算) の
CIP コードが元と全数一致。

#### R1.3 エラー・API 整備と公開準備
`ChemError` の細分化 (パースエラー位置情報の構造化)、`Molecule` ファサード型
(graph 直触りを減らす)、rustdoc、`cargo publish --dry-run` 通過。
ベンチマーク (criterion): SMILES 1 万件/秒オーダーの確認。

#### R1.4 ファジング
`cargo-fuzz` で SMILES パーサとサニタイズにファザーを掛け、
panic ゼロを確認 (任意バイト列 → Err で返ること)。
コーパス + PubChem キャッシュを種にする。

### Phase R2: SMARTS フルパーサ

現状は「SMILES 文字列の SMARTS 的再解釈」のみ。一般 SMARTS には:

#### R2.1 SMARTS 字句・論理式
原子プリミティブ (`#n`, `X`, `D`, `H`, `R`, `r`, `v`, `+/-`, `a`, `A`) と
論理演算子 (`!`, `&`, `,`, `;`) の AST。結合式 (`-`, `=`, `~`, `@`, `/`,`\`)。
ゲート: RDKit `MolFromSmarts` が受理/拒否するパターン集で挙動一致。

#### R2.2 マッチング統合と再帰 SMARTS
VF2 の原子/結合述語を SMARTS 式評価に差し替え。`$(...)` 再帰 SMARTS。
ゲート: 代表 SMARTS 100 種 × コーパス分子での `GetSubstructMatches` 全数一致
(フィクスチャ生成は gen_substruct_fixture.py を拡張)。

### Phase R3: 分子 I/O

#### R3.1 MOL/SDF (V2000) リーダ
座標・結合ブロック・プロパティブロック。サニタイズは既存パイプラインへ合流。
ゲート: SDF サンプル集 (PubChem からダウンロード) を RDKit と突き合わせ
(原子数・結合次数・電荷・芳香族)。

#### R3.2 MOL/SDF ライタ + InChI 前段
V2000 書き出し (2D 座標は R6 まで 0 埋め)。
InChI は公式 C ライブラリの FFI が現実的なので、ここでは
「InChI 入力用の正規化 (H 処理・電荷)」までを整備し、FFI はオプション機能とする。

### Phase R4: 記述子

#### R4.1 組成系記述子
分子式 (Hill 順)、厳密質量・平均分子量 (同位体テーブル拡充)、
重原子数、環数、HBD/HBA (Lipinski 定義)、回転可能結合数、TPSA (Ertl)。
ゲート: RDKit `Descriptors` とコーパス全分子で一致
(TPSA/logP 系は小数を 1e-6 で比較)。

#### R4.2 Crippen logP / MR
Crippen 原子タイプ分類 (SMARTS ベース — R2 に依存) と寄与テーブル。
ゲート: RDKit `MolLogP`/`MolMR` 全数一致。

### Phase R5: フィンガープリントと類似度

#### R5.1 Morgan (ECFP) フィンガープリント
半径 r の反復環境ハッシュ、ビット折りたたみ (2048bit)、カウント版。
RDKit のハッシュ関数と完全一致させるか独自ハッシュにするかを最初に決める
(推奨: RDKit 一致。ゲートが機械化できるため)。
ゲート: `GetMorganFingerprintAsBitVect` とコーパス全分子でビット列一致。

#### R5.2 類似度 + RDKit FP/MACCS
Tanimoto/Dice。MACCS keys (166 鍵 — SMARTS 定義、R2 依存)。
ゲート: RDKit と鍵ビット全数一致。

### Phase R6: 2D 座標生成と描画

#### R6.1 2D 座標 (基本)
環系テンプレート + 鎖の逐次配置 (RDKit 互換は要求しない。
交差最小・結合長一定の品質メトリクスをゲートにする)。

#### R6.2 SVG 描画
結合・原子ラベル・電荷・立体ウェッジの SVG 出力。ゴールデンファイルテスト。

### Phase R7: 標準化ユーティリティ (MolStandardize サブセット)

塩の除去 (最大フラグメント選択)、電荷中和、正規化変換 (ニトロ基表記統一等)。
ゲート: RDKit `rdMolStandardize` と代表セットで一致。

## 3. 優先順位の考え方

- smiles2iupac の残り移植 (RUST_PORT_PLAN.md Phase 3〜6) が主目的なら、
  本計画は **R1.1 と R1.2 だけ**先行する価値がある (128 原子上限と立体出力は
  PubChem 全域を扱う S6.2 で顕在化し得る)。他は移植完了後で良い
- 汎用ライブラリとして独立させるのが目的なら R1 → R2 → R4 → R5 の順が
  費用対効果が高い (R2 の SMARTS が R4.2/R5.2 の前提になる点に注意)
- 各ステップの規模目安: R1.x = 小〜中 (1 セッション)、R2/R5.1 = 大 (2 セッション分割可)、
  R3/R4/R6/R7 = 中

## 4. リスク

| リスク | 対策 |
|---|---|
| RDKit のバージョン差 (フィクスチャの再現性) | フィクスチャ生成時に RDKit バージョンを JSON に記録。参照は Release_2023_09 系に固定 |
| Morgan FP のハッシュ完全一致が困難 | RDKit の boost::hash 連鎖を先に小さな単体で移植・検証してから全体に進む |
| SMARTS 論理式の意味論の細部 | R2.1 で「受理/拒否 + マッチ数」の大規模フィクスチャを先に作る |
| molrs 変更による smiles2iupac の退行 | 全ゲート (rdkit_compat / canon / substruct / corpus) を CI で常時実行 |
