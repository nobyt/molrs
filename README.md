# molrs

Rust ケモインフォマティクスライブラリ。SMILES パース・分子グラフ・芳香族認識・
SSSR・CIP 立体化学・正準化・SMARTS サブセット・3D 配座生成 (距離幾何 + UFF +
ETKDGv2 実験トーション) を依存クレートなしで提供する `molrs` クレートと、
その上の IUPAC 2013 命名エンジン `smiles2iupac` クレート (移植中) からなる。

もとは [smiles2iupac](../smiles2iupac) (Python) リポジトリ内 `rust/` で開発され、
独立したもの。Python 実装 (参照コミット a01eccd) をオラクルとした差分テストで
検証している — フィクスチャは `corpus/` にあり、採取・codegen ツール
(RDKit/Python が必要) は smiles2iupac リポジトリの `tools/` にある。

- 開発サイクル・規約: [RUST_PORT_HANDOFF.md](RUST_PORT_HANDOFF.md)
- 汎用化ロードマップ: [RDKIT_RUST_PLAN.md](RDKIT_RUST_PLAN.md)
- 3D 配座生成 (完了): [RUST_3D_PLAN.md](RUST_3D_PLAN.md)

```bash
cargo test                             # 全ゲート (rdkit_compat 7453 分子 全数一致 ほか)
cargo test --release --test conformer_gate  # 3D 性質ゲート (release のみ)
cargo clippy --all-targets -- -D warnings
```
