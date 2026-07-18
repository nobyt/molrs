# SMILES → InChI / InChIKey (pure Rust) 実装計画

> **v1 実装状況 (2026-07-19)**: I0–I5 実装済み。実測 (コーパス 7,453 分子、
> RDKit オラクル):
> - 式層 (Hill): 99.13% 一致 (残差は電荷正規化・有機金属切断)
> - 正準番号 (AuxInfo /N:): **99.29%** 一致
> - **InChIKey 機構** (`inchi_key_from_string` vs RDKit `InchiToInchiKey`,
>   非立体 InChI 7,245 件): **7245/7245 = 完全一致**。SHA-256 (FIPS 180-4) +
>   base-26 (`ikey_base26.c` 移植、除外トリプレット対応) は公式とビット一致
> - フル InChI 文字列・InChIKey (v1 適用範囲 = 中性・単一成分・立体/同位体
>   なし): **74.17%** が RDKit と完全一致
> - CLI: `cargo run --bin inchi` (stdin SMILES → JSONL)
> - v1 の非対応 (→ `Unsupported`): 電荷 (q/p 層)・多成分・立体 (b/t/m/s)・
>   同位体 (i)・有機金属切断。残る文字列不一致の主因は可動 H 認識の未整備
>   (アミド/アミジン/ラクタム/環 N-H) と一部の縮合環 c 層直列化 — 逐次拡張。


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
