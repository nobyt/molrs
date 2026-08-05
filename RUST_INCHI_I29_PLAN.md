# I29 — PubChem 実データ検証と残バグ一覧 (2026-08-06)

> **要旨**: リポジトリ内コーパス (7,453 分子) は 100% 一致だが、それは
> 中性・単一成分に偏ったコーパスへの過適合だった。PubChem から CID 空間
> 全体をサンプリングした **18,563 分子**で測ると **94.66%**。確認済みの
> molrs 側バグは **981 件 (5.28%)**。本書はその内訳と優先度。

## 検証データの作り方

`corpus/pubchem_inchi.jsonl.gz` (18,563 分子) は PubChem FTP の
`CURRENT-Full/SDF` から作った。

- PUG-REST は IP 単位で絞られていて 503 (`PUGREST.ServerBusy`) を返すため
  使えない。`https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/`
  は普通に引ける。
- 全 358 ファイル (CID 1〜1.79 億) から 16 ファイルを等間隔で選び、各
  ファイルの先頭 900KB だけを **HTTP Range** で取得。gzip ストリームは
  先頭からなら途中まで展開できるので、巨大ファイルを落とさずに済む。
  低 CID だけに偏らないよう CID 空間を横断してサンプリングしている。
- オラクルは PubChem が **IUPAC 公式 InChI ソフトで計算した**
  `PUBCHEM_IUPAC_INCHI` / `_INCHIKEY`。RDKit 経由より直接的。
- 取得スクリプトは `scratchpad` に置いた使い捨て
  (`fetch_pubchem.py`、本書末尾に要点を残す)。

## 「molrs のバグ」と「SMILES の情報落ち」の切り分け

PubChem の InChI は CTAB から計算されているので、SMILES に情報が落ちて
いれば molrs のせいでない不一致が出る。**同じ SMILES を RDKit (公式 InChI
ライブラリ) に通して切り分けた**:

| | 件数 |
|---|---|
| molrs 一致 | 17,146 |
| **molrs のバグ** (RDKit == PubChem なのに molrs だけ違う) | **1,407** |
| SMILES 由来 (RDKit も PubChem と違う) | 10 |

**不一致の 99.3% は molrs 側**。SMILES ラウンドトリップの問題ではない。

## 修正済み: 多成分の成分順序 (I29)

`formula.rs` の `component_sort_key` が
「炭素を含む成分が先 → **重原子数昇順** → H 数降順」だった。これは
リポジトリ内コーパスの多成分 32 例 (`2Na.H2O4S` 系が多数) から導出した
規則で、**PubChem 実データ 863 件では 33.7% しか再現できない**。

正しい規則は **炭素数降順 → 重原子数降順 → H 数降順 → 式の辞書順**
(実測 96.2%)。炭素数が主キーであることは 863 件で例外なく成立する。

```
got  InChI=1S/C2H4O2.C6H5.Hg/...     ← 重原子数昇順 (4 < 6)
want InChI=1S/C6H5.C2H4O2.Hg/...     ← 炭素数降順 (6 > 2)
```

**効果**: PubChem 92.37% → **94.66%** (+426 件)。

### 既知の残差 (無機塩、33/863 = 3.8%)

「単原子カチオン + 多原子アニオン」だけは実 InChI がカチオンを先に置く
(`2Na.H2O4S`、`5Na.H3O4P.H2O3S`、`Cu.N2O4.2NO3`)。ただし `FH.O3Si.2Zn` は
単原子の Zn が最後なので「単原子金属を先頭」では説明できない (電荷層との
関係も含めて未解明)。単原子金属フラグを足す案は実測でスコアが下がった。

この変更でリポジトリ内コーパスの `[Na+].[Na+].[O-]S(=O)(=O)[O-]` **1 件が
退行**する。PubChem で +426 件と引き換えの取引として受け入れ、
`inchi_gate.rs` に理由付きの明示的な例外として書いた
(`known_divergence_monoatomic_cation_salts` も参照)。

## 残バグ一覧 (981 件、優先度順)

### 1. 立体 `/t` — 未定義中心 `?` を出力していない (239 件) ★最優先

**最大かつ最も明確なバグ**。構造上は立体源性だが SMILES で配置が指定
されていない中心を、実 InChI は `9?` のように `/t` に列挙する。molrs は
**完全に省略している**。`/t` だけが違う 422 件のうち 239 件がこれで、
molrs が中心を**多く**出したケースは 0 件 (常に不足)。

```
smi   CCC(C(=O)O1)S(=O)(=O)c2ccc(C)cc2 …  (CC[C@@H]1CC(C(=O)O1)S(=O)(=O)c2ccc(C)cc2)
got   /t10-
want  /t10-,12?
```

`stereo.rs` の `/b` 側には未定義二重結合を `?` として併記する処理が既に
ある (`end_stereogenic`)。同じ考え方を四面体中心に用意すればよい。
`/m` `/s` の食い違い (計 ~145 件) の多くもこれに連動しているはず。

### 2. 立体 `/t` — パリティが反転する中心 (147 件)

中心の数は合っているが 1 個 (稀に複数) のパリティが逆。第四級炭素・
スピロ・S/B/N+ を含む置換基で目立つ (S 102 件、N+ 22、B 21、Si 11、P 10)。

```
smi   C[C@H]1C[C@@H]([C@](O1)(CO)[Si](C)(C)C2=CC=CC=C2)O
got   /t11-,13-,14+/m0/s1
want  /t11-,13-,14-/m0/s1
```

molrs は CIP の R/S を正準番号基準のパリティへ変換している
(`raw = '-' iff (rs_bit XOR perm(CIP昇順→正準昇順))`) ので、CIP ランクが
ずれると反転する。Si (原子番号 14) や S、ホウ素まわりの CIP 実装を疑う。

### 3. 電荷正規化 — 分子内塩 (zwitterion) (115 件)

正味電荷 0 の分子内塩を、molrs はカルボキシラートを protonate してから
`/p-1` で戻すため `/q+1/p-1` が付き、式も H が 1 個多くなる。

```
smi   CC(=O)OC(CC(=O)[O-])C[N+](C)(C)C     (アセチルカルニチン)
got   InChI=1S/C9H18NO4/…/q+1/p-1
want  InChI=1S/C9H17NO4/…              ← 層なし
```

正味 +1 の側も `/q+1` と `/p+1` を取り違える。`normalize.rs` が
「四級 N のような**恒久電荷**」と「プロトン化で動かせる電荷」を区別して
いないのが原因と見られる。

### 4. `/b` 層 (106 件)

余分に出す例が目立つ。可動 H 群にかかる C=N (アミジン型) で、実 InChI は
二重結合の位置自体が動くため `/b` を出さない。

```
smi   CC1=CC(=C(C(=C1S(=O)(=O)/N=C(\C2=CC=CC=C2)/N)C)C)OCCCC(=O)O
got   …(H2,21,22)(H,23,24)/b22-20+
want  …(H2,21,22)(H,23,24)          ← /b なし
```

`double_bond_layer` は未定義 (`?`) 側でのみ可動 H 群を除外している
(`tautomer_group_members`) が、**定義済みの `b.stereo` を持つ結合には
同じ除外を掛けていない**。

### 5. 128 原子超でパニック (87 件) ★堅牢性

`rings.rs:38` の `assert!(n_atoms <= 128, "molecule too large for ring
perception")`。ペプチド等でごく普通に踏む。**ライブラリ API が正当な入力で
パニックする**のは呼び出し側から回避できず、`catch_unwind` を強いる。
最低でも `Err(Unsupported)` を返すべき。環認識が 128 bit のビットセットに
依存しているなら動的ビットセットへの置き換えが要る。

### 6. 可動 H 群の分割 (`/h`) (25 件)

```
smi   C1=NC(=C(N1)C(=O)O)NC(=O)N
got   /h1H,(H,10,11)(H4,6,7,8,9,12)
want  /h1H,(H,7,8)(H,10,11)(H3,6,9,12)
```

I25〜I28 で詰めた群構築が、実データではまだ大きく取りすぎることがある。

### 7. `/i` 同位体層が未実装 (13 件)

```
smi   [3H]C([3H])(…)SC
want  …/t7-,8-,9-,10-/m1/s1/i5T2      ← got は /i 以降がない
```

`inchi_gate.rs` は `/i` を含む分子を除外しているので既存コーパスでは
見えていなかった。

### 8. SMILES パーサが超原子価ハロゲンを拒否 (8 件)

`valence 7 for atom N (Cl)` (過塩素酸型)、`valence 3 for atom N (I/Br)`
(超原子価ヨウ素) を `Invalid SMILES` で弾く。RDKit は受け付ける。

## 検証サイクル

```bash
rtk proxy cargo test --release --test inchi_gate   -- --nocapture --test-threads=1  # 7,453 件・完全一致
rtk proxy cargo test --release --test pubchem_gate -- --nocapture                   # 18,563 件・94.66%
```

`pubchem_gate.rs` は 128 原子超のパニックを `catch_unwind` で捕まえて
不一致として数えている (バグ 5 が直ったら外すこと)。

## 再取得スクリプトの要点

```python
BASE = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
# 358 ファイルから等間隔に 16 個選び、各先頭 900KB を Range で取得
req = Request(BASE + fn, headers={"User-Agent": ..., "Range": "bytes=0-900000"})
text = zlib.decompressobj(31).decompress(urlopen(req).read()).decode("utf8", "replace")
for rec in text.split("$$$$\n")[:-1]:          # 末尾は切れているので捨てる
    # PUBCHEM_SMILES / PUBCHEM_IUPAC_INCHI / PUBCHEM_IUPAC_INCHIKEY を拾う
```
