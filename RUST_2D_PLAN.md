# 2D 構造式描画 (depict) 実装計画

> **完了 (2026-07-17)**: D0〜D13 全ステップ実装済み。実測結果:
> - コーパス 7,453 分子 + 補助立体セット 500 分子で性質ゲート
>   (`tests/depict_gate.rs`, DEPICT_GATE_FULL=1) 全項目 green —
>   レイアウト成功 100% / 非環結合長・30° 量子化 100% / 環結合長 99.4%
>   (橋かけ内部のみ逸脱) / 立体 round-trip 708/708 / SVG well-formed / 決定的
> - RDKit 外部オラクル (smiles2iupac tools/check_depict_stereo.py):
>   くさび付き 2D MOL → RDKit 再認識 → 異性体込み正準 SMILES 一致 **708/708**
> - 既知の制約: (1) かご型 (キュバン) は 2D で本質的に交差 — ゲート例外 1 分子。
>   (2) CIP ラベルの付かないパリティのみの立体中心 (橋頭など) はグラフが
>   CIP ラベルしか持たないため表現不能。(3) Wiley プリセットは ACS 同値の
>   暫定値 (公式一次情報が未確認)。
> - 目視レビュー: `... | cargo run --release --bin depict_gallery > gallery.html`

関連文書: [RUST_PORT_HANDOFF.md](RUST_PORT_HANDOFF.md) (開発サイクル)、
[RDKIT_RUST_PLAN.md](RDKIT_RUST_PLAN.md)、[RUST_3D_PLAN.md](RUST_3D_PLAN.md) (完了済み 3D 計画)

## Context

molrs (SMILES → 分子グラフ → 3D 配座まで実装済みの Rust ケモインフォマティクスライブラリ) に、
分子の 2D 構造式描画機能を追加する。要件は 2 つ:

1. **IUPAC 標準準拠**: Graphical Representation Standards for Chemical Structure
   Diagrams (IUPAC Recommendations 2008, GR-0〜GR-13) および立体配置の表現は
   Graphical Representation of Stereochemical Configuration (IUPAC Recommendations
   2006) に従ったレイアウト・描画を行う。
2. **雑誌書式対応**: 描画パラメータ (結合長・線幅・フォント等) をスタイルプリセットとして
   差し替え可能にし、ACS 1996 ドキュメント設定等の雑誌別プリセットを同梱する。

2D 描画は「**座標生成 (layout)**」と「**描画 (render)**」の 2 段からなり、
これは既存の 3D パイプライン (embed → optimize) と同型の構成。

## 調査済みの規格要件 (一次情報から抽出)

### IUPAC 2008 (iupac.qmul.ac.uk/drawing/drawing.html — 全文 193K 字を取得済み)

レイアウトに効く定量規則:
- **鎖**: ジグザグ、隣接結合は 120° (GR-4.1)。真に直線的な原子 (sp、アレン中心、
  2 個以上の二重結合を持つ原子) は 180°。4 結合は 90°×4 も可 (GR-4.1)
- **環**: 正多角形で描く (GR-3.3)。角度は 120° を基本に、五員環 108°、七員環 129°、
  八員環 135° を許容。五員環の「切頂六角形」形 (120°×3 + 90°×2) も許容
- **環置換基**: 環外結合は隣接環結合の角を二等分する方向 (GR-4.2)。環原子に置換基
  2 つなら 60° 開き (四員環のみ 90° 可)
- **配向**: 主鎖は水平、±30° の結合を最大化 (GR-3.2, GR-4.2.3)
- **原子ラベル** (GR-2): 無標識原子=C (適正数の H を暗黙に持つ)。ヘテロ原子は
  元素記号+H 数 (NH, OH, NH₂)。結合が左から付く場合は H を左置き (HO-, H₂N-)
- **二重結合** (GR-1.6, GR-1.10): 平行 2 本線。環内は内側に短い 2 本目
  (sidedness 規則)、末端・対称位置ではセンター振り分け
- **重なり回避** (GR-4.3): 置換基同士の衝突は角度調整・結合の伸長よりも
  対称的な角度再配分を優先
- E/Z は例外なく正しい幾何で描く (GR-1.10 系)。cis/trans を歪めてはならない
- 芳香環は原則ケクレ形 (内円は非推奨 — GR-6 で「曲線=非局在」の乱用を禁止)
- 電荷は原子の右上に付す (GR-5.1)、など

### IUPAC 2006 (立体、参照文献 [6]): 楔形結合
- 細端を立体中心側に置く solid wedge (手前) / hashed wedge (奥)
- 立体中心 2 つを結ぶ結合には楔を使わない、環内結合への楔は避ける等の配置規則
- 実装時は 2006 勧告 (Pure Appl. Chem. 78, 1897) の ST 規則を参照する

### ACS 書式 (pubs.acs.org graphics_prep — Wayback 経由で全設定値を取得済み)

ChemDraw「ACS Document 1996」相当の描画設定:

| 項目 | 値 |
|---|---|
| chain angle | 120° |
| fixed (bond) length | **14.4 pt (0.2 in)** |
| bond spacing (二重結合間隔) | **18% of length** |
| line width | **0.6 pt** |
| bold width | 2.0 pt |
| margin width (ラベル周囲の空白) | 1.6 pt |
| hash spacing | 2.5 pt |
| font / size | Helvetica (Arial) / **10 pt** |
| 図の最大幅 | 単段 3.25 in / 二段 7 in、高さ 9.5 in |

### 各誌プリセット値 (一次情報から確認済み)

| パラメータ | ACS 1996¹ | Nature² | RSC³ | Wiley⁴ |
|---|---|---|---|---|
| bond length | 14.4 pt (0.508 cm) | 10.8 pt (0.381 cm) | 12.2 pt (0.43 cm) | 14.4 pt |
| line width | 0.6 pt | 0.6 pt (0.021 cm) | 0.5 pt (0.016 cm) | 0.6 pt |
| bold/wedge width | 2.0 pt | 1.56 pt (0.055 cm) | 1.6 pt (0.056 cm) | 2.0 pt |
| double bond spacing | 18% of length | 18% | **20%** | 18% |
| hash spacing | 2.5 pt | 1.7 pt (0.06 cm) | 1.8 pt (0.062 cm) | 2.5 pt |
| margin width | 1.6 pt | 1.19 pt (0.042 cm) | (規定なし→1.6pt) | 1.6 pt |
| atom label font | Helvetica/Arial 10 pt | Arial/Helvetica **6 pt** | Arial/Helvetica **7 pt** | Helvetica 10 pt |
| chain angle | 120° | 120° | 120° | 120° |
| 図幅上限 | 3.25 in / 7 in | — | — | — |

¹ pubs.acs.org graphics_prep (Wayback 取得)。
² nature.com/documents/nr-chemical-structures-guide.pdf (公式 PDF、取得済み)。
  付随規則: 立体は楔のみ (太線結合不可)、楔の細端=立体中心 (遠近表現ではない)、
  H は立体指定に必要な場合のみ表示、ラベルに太字不可。
³ rsc.org Chemical Science 投稿ガイド (公式ページ取得)。stereo bond width 1.6 pt、
  TIFF 600 dpi の別規定あり。
⁴ Wiley は公式の数値一覧が確認できず。ChemDraw「Wiley Document」スタイル
  シート = ACS 同値との複数の二次情報に基づき ACS 値で実装し、要確認と注記する。

→ スタイル構造体はこの語彙 (bond_length, line_width, bold_width, bond_spacing_ratio,
hash_spacing, margin_width, font_family, font_size, max_width) を持ち、
プリセットは `Style::iupac_default() / acs_1996() / nature() / rsc() / wiley()` の
数値セット関数として提供する。レイアウトは結合長=1 の無次元座標で行い、
描画時にスタイルでスケールする (レイアウトとスタイルの直交性)。

## 確定済みの方針

- 出力: **SVG + 2D MOL** (V2000 くさびコード付き)。依存ゼロ維持 (文字列生成のみ)。
  PNG/EPS は外部変換に委ねる
- 置き場所: **`molrs::depict` モジュール** (conformer と同格)
- プリセット: **IUPAC 既定 / ACS 1996 / Nature / RSC / Wiley** の 5 種を初回から同梱

## 設計 (要点)

レイアウトは結合長 = 1.0 の無次元単位で計算し、描画時に Style がスケール。

1. **前処理**: 隠し H 決定 (H は明示ノードなので隣接 H を数えてラベルに畳む)、
   ケクレ化 (芳香環は Kekulé 描画)、`cip_ranks` 取得、環系グルーピング
   (SSSR を共有原子で連結成分化し fused/spiro/bridged を分類)
2. **環系**: 単環 = 正多角形。縮合 = 環隣接 BFS + 共有辺貼り付け。
   矛盾時 (橋かけ・高縮合)・大員環 (≥9) は既存の閉多角形ソルバ
   `planar_ring_coords` (conformer/bounds.rs:89、private → pub(crate) 化) で解く
3. **鎖**: ±120° ジグザグ、sp/累積は 180°、4 結合は 90°×4 許容。
   環置換基は外角二等分、同一環原子の 2 置換基は 60° (四員環 90°)
4. **E/Z**: stereo タグを cip_ranks で幾何側に翻訳しレイアウトで強制。
   以降の全変換は E/Z 保存操作のみ (GR: 例外なく正しい幾何)
5. **衝突解消** (GR-4.3): グリッド検出 → 部分木反転 → 対称角再配分 →
   部分木回転 → 最終手段の結合伸長 (記録)。決定的
6. **全体配向**: 2×2 慣性主軸の水平化 (`jacobi_eigen`) → 24 配向候補から
   「水平 ±30° の結合数」最大を選択。フラグメント (塩) は横並び
7. **くさび (IUPAC 2006)**: 細端 = 立体中心。候補結合を優先順 (再表示 H >
   末端非環単結合 > 非環単結合、立体中心間・環結合は回避) で選択。
   solid/hashed の向きは「+z/−z に持ち上げ → 符号付き体積 → CIP 再導出」が
   入力 chiral_tag (R/S) と一致する側を採用 (graph はパリティでなく CIP
   ラベルを持つため再導出一致方式が唯一の正解)

### モジュール構成と公開 API

```
molrs/src/depict/
  mod.rs          # 公開 API・DepictError・パイプライン統括
  point2.rs       # Point2 (演算・回転・30° スナップ)
  style.rs        # Style + 5 プリセット
  ring_layout.rs  # 環系グルーピング・多角形構築
  chain_layout.rs # 鎖・置換基角・E/Z 側決定
  place.rs        # 組み立て・配向・フラグメント
  collide.rs      # 衝突検出・解消
  stereo2d.rs     # くさび選択 + verify_stereo_2d
  label.rs        # ラベル・bbox (Helvetica 近似文字幅テーブル)
  svg.rs          # SVG 文字列生成
  molblock2d.rs   # V2000 2D + くさびコード
```

```rust
pub struct LayoutParams { pub seed: u64, /* 反復上限・伸長許可 */ }
pub enum WedgeDir { Up, Down }        // narrow 端 = bond.begin_idx
pub struct Coords2D { pub pos: Vec<Point2>, pub hidden: Vec<bool>,
                      pub wedge: Vec<Option<WedgeDir>> }
pub fn compute_coords_2d(g: &MoleculeGraph, p: &LayoutParams)
    -> Result<Coords2D, DepictError>;
pub fn to_svg(g: &MoleculeGraph, c: &Coords2D, s: &Style) -> String;
pub fn depict_svg(smiles: &str, s: &Style) -> Result<String, ChemError>;
pub fn to_mol_block_2d(g: &MoleculeGraph, c: &Coords2D, title: &str) -> String;
```

### 既存コードの再利用 (調査で検証済み、コピーはゼロ)

| 対象 | 現状 | 措置 |
|---|---|---|
| `conformer/bounds.rs:89 planar_ring_coords` | private | **pub(crate) 化** (D5) |
| `conformer/minimize.rs minimize_with` (汎用 CG) | pub(crate) | 変更不要 (z 勾配 0 で 2D 微調整に利用可) |
| `stereo.rs cip_ranks` | pub(crate) | 変更不要 |
| `aromaticity.rs kekulize` | pub(crate) だが内部型依存 | 新ヘルパ `kekulized_bond_orders(g) -> Vec<f64>` を追加 (D4) |
| `geometry.rs jacobi_eigen / SeededRng` | pub | そのまま使用 |
| `conformer/molblock.rs to_mol_block` | pub | 行フォーマットの雛形として踏襲 (新規実装) |

## 実装ステップ (D0–D13、各 1 Sonnet セッション)

- **D0 (S)**: 本計画を `RUST_2D_PLAN.md` として molrs リポジトリに保存
  (リポジトリの計画文書の慣行に従う)
- **D1 (S)**: 骨格 — `pub mod depict`、Point2、Style + **5 プリセット**
  (iupac_default/acs_1996/nature/rsc/wiley — 上表の確認済み数値、
  Wiley は ACS 同値である旨と要確認をコメントに明記)、
  LayoutParams/Coords2D/WedgeDir/DepictError。
  完了条件: 演算・プリセット値のユニットテスト、clippy clean
- **D2 (M)**: 前処理 + 鎖レイアウト (無環分子)。
  完了条件: 直鎖/分岐/アレン/E,Z-2-ブテンで結合長 1.0±1e-6、
  角度 {120,180,90}°、E/Z を座標再導出して入力一致
- **D3 (M)**: 最小 SVG レンダラ — **ここで end-to-end 成立**。
  完了条件: golden SVG ~5 分子 (バイト一致 = 決定性)、well-formed チェック
- **D4 (M)**: 単環 + ケクレ化ヘルパ。
  完了条件: ベンゼン/ピリジン/トルエン等 golden、コーパス単環で重なりなし
- **D5 (L)**: 縮合環系・スピロ (共有辺貼り付け、planar_ring_coords 可視性変更)。
  完了条件: ナフタレン/インドール/デカリン/ステロイド骨格/スピロ環 golden
- **D6 (M)**: 橋かけ環 (外周を planar_ring_coords + 橋を内側)・大員環 ≥9。
  完了条件: ノルボルナン/アダマンタン/クラウンエーテルで交差最小、
  コーパス全環系で座標生成 100%
- **D7 (M)**: 統合・全体配向・フラグメント横並び。
  完了条件: **コーパス 7,453 分子で成功 100%・決定性 (2 回実行一致)**
- **D8 (L)**: 衝突検出・解消。
  完了条件: 非結合原子間距離 ≥0.5L 違反 ≤0.5% (違反リストをテストに固定)
- **D9 (L)**: 立体くさび + `verify_stereo_2d`。
  完了条件: コーパス立体分子で round-trip 100%。コーパスは立体が薄い
  (@ 約 50 行・/ 約 158 行) ため **追加立体フィクスチャ ~500 分子を
  smiles2iupac/tools で生成**し同 100%
- **D10 (M)**: ラベル精緻化 — bbox クリップ (margin_width)、H 左置き (HO–)、
  下付き数字、電荷右上、二重結合の内側線 (sidedness)・短縮、くさび描画
  (hash_spacing 間隔)。完了条件: golden ~20 分子、ラベルと結合線の交差 0
- **D11 (S)**: 5 プリセット検証 (同一 Coords2D から全 preset の SVG 生成 =
  レイアウトとスタイルの直交性テスト) + `to_mol_block_2d` (V2000、2D フラグ、
  wedge codes 1/6)。完了条件: golden MOL + 全 preset well-formed
- **D12 (M)**: 性質ゲート `molrs/tests/depict_gate.rs` (conformer_gate 踏襲:
  決定的 1/10 サンプル + FULL 環境変数)。8 項目: 成功率 100% / 結合長 1.0±2%
  (例外リスト付き) / 無環結合の 30° 量子化率 ≥95% / 重なり違反 ≤0.5% /
  E/Z 再導出 100% / wedge→CIP round-trip 100% / SVG well-formed / 決定性
- **D13 (M)**: HTML ギャラリー dev bin (`src/bin/depict_gallery.rs`、目視用) +
  **RDKit 外部オラクル**: smiles2iupac/tools に 2D MOL → RDKit CIPLabeler で
  R/S・E/Z 再認識するフィクスチャ生成スクリプト → `rdkit_depict_compat.rs` で
  全数一致 (molrs の CIP 実装と独立な end-to-end 検証)

規模: S×3, M×8, L×3。D3 で最短の見える成果、D7 でコーパス全域、D12-13 が品質保証。

## リスク

| リスク | 対策 |
|---|---|
| 橋かけ多環の 2D 慣用配置に一般解がない | 外周+内部配置で「重なりなし」のみ保証。ゲートは審美でなく性質 |
| ケクレ化再実行が graph 構築時と食い違う | ビルダーと同経路で再構成、決定性テストで固定 |
| 混雑分子で衝突解消が非収束 | 反復上限 + 明示的結合伸長 + 例外リスト固定 (回帰検知) |
| フォントなしのラベル bbox 近似 | 文字幅テーブル + 安全マージン + ギャラリー目視 + D13 オラクル |
| golden SVG の脆さ | golden ≤20 分子に限定、主守りは性質ゲート |
| Wiley プリセットの公式値未確認 | ACS 同値で実装しコメントに要確認と出典を明記 |

## 検証 (エンドツーエンド)

1. `cargo test` — ユニット + golden + `depict_gate.rs` (コーパス性質ゲート 8 項目)
2. `cargo clippy --all-targets -- -D warnings`
3. `cargo run --bin depict_gallery > gallery.html` で目視レビュー
4. D13: `~/ghq/github.com/nobyt/smiles2iupac/tools/` (RDKit venv) で
   2D MOL round-trip フィクスチャ生成 → `rdkit_depict_compat.rs` 全数一致
5. 立体: verify_stereo_2d (内部) + RDKit CIPLabeler (外部) の二重検証
