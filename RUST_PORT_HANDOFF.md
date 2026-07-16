# Rust 移植 引き継ぎドキュメント

次の実装セッション (Sonnet クラス想定) がそのまま続きに着手できるようにまとめたもの。
全体計画は [RUST_PORT_PLAN.md](RUST_PORT_PLAN.md)。**着手前に本書と計画書の該当ステップ節を読むこと。**

## 1. 現在地 (2026-07-13 時点)

| 項目 | 状態 |
|---|---|
| 完了ステップ | S0.1〜S2.2 (計画 34 ステップ中 12) |
| ブランチ | `rust-port` (worktree: `.claude/worktrees/rust-port`) |
| 最終コミット | `4357733` S2.1+S2.2 |
| **Python 参照コミット** | **`a01eccd`** (ブランチ分岐点。コーパス・RDKit ダンプ・constants.rs は全てこのコミットの Python 実装から生成した) |
| コーパス合格 | **196 / 7,453** (保留名 + 直鎖アルカン) |
| テスト | cargo test 61 件、clippy `-D warnings` クリーン |

注意: 元リポジトリの `main` は分岐後も進んでいる (`fccb5b1` 時点確認)。
Python 側の変更を取り込む場合は rebase 後にコーパス抽出 (S0.1) と
RDKit ダンプ (下記) を再生成し、全ゲートを回し直すこと。

## 2. 何がどこにあるか

```
rust/
  molrs/            RDKit 代替層 (Phase 1 完了、7,453 分子で RDKit と全数一致検証済み)
    src/smiles/         SMILES パーサ (neighbor_order が立体化学の鍵)
    src/graph.rs        MoleculeGraph 構築 (build_molecule_graph が唯一の入口)
    src/aromaticity.rs  ケクレ化 + 芳香族認識 (RDKit デフォルトモデル互換)
    src/rings.rs        対称化 SSSR (RDKit FindRings.cpp の忠実移植)
    src/canon.rs        正規 SMILES (RDKit 非互換・自己完結。Strict/Lenient 2 モード)
    src/substructure.rs VF2 (mol クエリ + SMILES-as-SMARTS)
    src/stereo.rs       レガシー CIP (R/S, E/Z; RDKit Chirality.cpp の忠実移植)
  smiles2iupac/         命名エンジン (Phase 2 の骨格まで)
    src/constants.rs    ★自動生成 — 手で編集禁止。再生成: uv run python tools/gen_constants_rs.py
    src/retained.rs     保留名ルックアップ (複合キー: 正規SMILES + CIP/EZ署名)
    src/lib.rs          smiles_to_iupac パイプライン (骨格)
    tests/corpus.rs     コーパス回帰ハーネス
    tests/expected_pass.txt  合格リスト (★単調増加。1 件でも落ちたら fail)
corpus/
  corpus.jsonl              (smiles, expected, phases) 7,453 件 — tools/extract_corpus.py で生成
  rdkit_dump.jsonl.gz       Python 版 MoleculeGraph の全ダンプ — tools/dump_rdkit.py で生成
  canon_pairs.jsonl         同一分子の異表記グループ 498 件
  substruct_fixture.jsonl.gz  GetSubstructMatches 正解 2,995 ペア
tools/
  extract_corpus.py / dump_rdkit.py / gen_constants_rs.py /
  gen_canon_pairs.py / gen_substruct_fixture.py / diff_oracle.py
```

## 3. 開発サイクル (毎ステップ共通)

1. RUST_PORT_PLAN.md の該当ステップ節を読む
2. 参照元 Python ファイル (`src/smiles2iupac/*.py` @ a01eccd) を読み、**挙動を変えずに**移植する
3. 検証は必ず「Python をオラクルにした機械比較」で行う (目視・自作期待値は不可):
   - 中間関数の移植 → Python 側からスナップショットを JSON 抽出するツールを
     `tools/` に書き、フィクスチャを `corpus/` にコミットし、Rust 統合テストで突き合わせる
     (前例: `gen_substruct_fixture.py` + `tests/substruct_fixture.rs`)
   - 名前の一致 → コーパスハーネスに任せる
4. 合格リスト更新: `UPDATE_EXPECTED_PASS=1 cargo test -p smiles2iupac --test corpus`
   → 再実行してグリーン確認 → expected_pass.txt をコミットに含める
5. `cargo fmt` / `cargo clippy --all-targets -- -D warnings` / `cargo test` (rust/ ディレクトリで)
6. コミット (係: 1 ステップ = 1 コミット、メッセージにステップ番号)

**リグレッション規律: expected_pass.txt から 1 件でも脱落したら、原因を直すまで先に進まない。**

## 4. 確立済みの方式 (変えないこと)

- **Python が単一情報源**: データテーブルは手移植せず codegen (`gen_constants_rs.py` 方式) で生成
- **正規 SMILES は RDKit 非互換**で良い。SMILES キーのテーブルを追加するときは
  生 SMILES のまま持ち、初期化時に `retained.rs` の `canonical_key()` で再キー化する
- **複合キー**: 立体は正規 SMILES に書かず「正規出力位置 + CIP/EZ コード」の署名で区別する
- **molrs は凍結気味に扱う**: RDKit 互換ゲート (tests/rdkit_compat.rs) が
  7,453 分子で全数一致している。命名ロジック都合で molrs を変えたくなったら
  まず smiles2iupac 側で吸収できないか考える。変えた場合は必ず全ゲート再実行
- 数値比較の注意: Python の bond_order は float (1.0/1.5/2.0/3.0)。Rust も f64 で同じ値

## 5. ハマりどころ知識ベース (実測で確定済みの RDKit 挙動)

グラフ構築 (S1.2 で検証済み — 変更時はここを壊さない):
- `num_hs` は常に 0。H は明示原子ノード (AddHs 相当、重原子順に末尾付加)
- 素の `[H]` は隣接重原子へマージ。`[2H]`/`[H+]` は原子として残る
- 結合リスト順: 鎖結合出現順 → 環閉じ結合が末尾に環番号順 (同番号は閉じ順)
- 環閉じの向き: 開き側に次数記号 (`=1` 等) があれば (開,閉)、なければ (閉,開)

芳香族認識 (S1.3):
- 電子寄与はグローバル (環内二重結合 = 1、どの環かは不問)。孤立電子対は固定 2
- 環外二重結合: 相手が電気陰性 → 0 (キノン)、C=C → 候補外 (フルベン)
- ユニオンは縮合**ペアのみ**、周縁結合だけマーク (アズレンの縮合結合は単結合)

環 (S1.4): `AtomRings()` は**対称化** SSSR (ビシクロ[2.2.2]=3 環)。環内原子順は BFS 経路順

部分構造 (S1.6): mol クエリの電荷は「クエリ側が非ゼロのときだけ」一致要求。
SMARTS はサニタイズなし (ケクレ化不能パターンも有効)

立体 (S1.7): `[C@@H](Cl)(F)Br` (S) と `Cl[C@@H](F)Br` (R) は**別分子**
(先頭原子 + 明示 H の特例)。8 員未満の環内二重結合に E/Z は付かない

## 6. 次のステップ: S3.1 官能基検出

- 参照: `src/smiles2iupac/functional_group.py` (3,075 行) @ a01eccd
- 進め方の提案:
  1. `tools/dump_functional_groups.py` を書く: コーパス全分子について
     `detect_groups(graph)` / `principal_group(...)` の結果 (官能基種・原子インデックス) を
     JSON ダンプ (dump_rdkit.py の形式を踏襲、gzip で corpus/ へ)
  2. `smiles2iupac/src/functional_group.rs` に `FunctionalGroup` struct と
     `detect_groups` / 60 個超の `_is_*` 述語を移植 (述語は相互独立なので上から順に)
  3. 統合テストでダンプと全数突き合わせ (rdkit_compat.rs の形式を踏襲)
- Python 側の `graph.atoms[i].num_hs` は 0 である点に注意 (H は adjacency で数える)。
  述語が「H 数」を見る箇所は Python も明示 H ノードを数えているはずなので、そのまま移植する
- S3.2 (chain_finder) 以降も同じ「スナップショット駆動」で進める

## 7. 環境メモ

- ビルド/テストは `rust/` ディレクトリで実行 (ワークスペースルート)
- Python 側は `uv run python ...` (リポジトリルートで; .venv は uv が管理)
- RDKit ダンプ再生成: `uv run python tools/dump_rdkit.py` (Python 実装を変えた時のみ)
- CI: `.github/workflows/rust.yml` (fmt / clippy / test)
- コミットは日本語コミットメッセージ不可の縛りはないが、既存は英語 + ステップ番号
