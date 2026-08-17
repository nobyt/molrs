//! PubChem 実データに対する InChI 差分ゲート (I29)。
//!
//! フィクスチャ `corpus/pubchem_inchi.jsonl.gz` は PubChem の
//! `CURRENT-Full/SDF` を CID 空間全体 (CID 1〜1.65 億) から 16 ファイル
//! 抽出してサンプリングした 18,563 分子。オラクルは PubChem が **IUPAC 公式
//! InChI ソフトで計算した** `PUBCHEM_IUPAC_INCHI` / `_INCHIKEY`。
//!
//! `inchi_gate.rs` の 7,453 分子コーパスは中性・単一成分に偏っており、
//! そこで 100% でもこちらでは 94.7% しか一致しない。塩・多成分・立体・
//! 電荷・同位体を含む実データでの現在地を測り、退行を防ぐのが目的。
//!
//! **不一致の 99.3% は molrs 側のバグ**であることを確認済み: 同じ SMILES を
//! RDKit (公式 InChI ライブラリ) に通すと PubChem と一致するため、SMILES の
//! 情報落ちではない (1,417 件中 1,407 件)。

use std::io::Read;
use std::path::PathBuf;

use flate2::read::GzDecoder;
use molrs::graph::build_molecule_graph;

struct Record {
    smiles: String,
    inchi: String,
}

fn load() -> Vec<Record> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../corpus/pubchem_inchi.jsonl.gz");
    let file =
        std::fs::File::open(&path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    let mut text = String::new();
    GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("gunzip");
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let v: serde_json::Value = serde_json::from_str(l).expect("json");
            Record {
                smiles: v["s"].as_str().unwrap().to_string(),
                inchi: v["inchi"].as_str().unwrap_or("").to_string(),
            }
        })
        .collect()
}

/// 生成できなければ None (超原子価ハロゲンの `InvalidSmiles` など) で、
/// 「不一致」として数える。
///
/// I29 の時点では 128 原子超が `rings.rs` の `assert!` でパニックしたため
/// ここで `catch_unwind` していた。I31 で `ChemError::Unsupported` に変え、
/// I33 で上限そのものを撤廃したので、どちらも不要になった。
fn try_inchi(smiles: &str) -> Option<String> {
    let g = build_molecule_graph(smiles).ok()?;
    molrs::inchi::to_inchi(&g).ok()
}

#[test]
fn pubchem_full_inchi() {
    let recs = load();
    let mut ok = 0usize;
    let mut n = 0usize;
    for r in &recs {
        if r.inchi.is_empty() {
            continue;
        }
        n += 1;
        if try_inchi(&r.smiles).as_deref() == Some(r.inchi.as_str()) {
            ok += 1;
        }
    }
    let acc = ok as f64 / n.max(1) as f64;
    println!("pubchem full InChI: {ok}/{n} exact ({:.2}%)", acc * 100.0);
    // I29 94.66% → I30 96.97% (未定義四面体中心 `?` + 第四級中心のパリティ)
    //            → I31 97.05% (成分単位の電荷中性化)
    //            → I32 97.80% (可動 H 群にかかる /b の除外、/m の `.` 連結)
    //            → I33 98.01% (環認識の 128 原子上限を撤廃)
    //            → I34 98.21% (孤立電子対の立体中心 = スルホキシド)
    //            → I35 98.28% (超原子価ハロゲンの原子価表)
    //            → I36 98.35% (脱プロトンを酸性 O-H に限定)
    //            → I37 98.44% (同位体層 /i)
    //            → I38 98.70% (可動電荷: 荷電可動 H 群の併合と酸点の拡大)
    //            → I39 98.71% (酸対に S を含む場合の N 端点、primary 限定を撤廃)
    //            → I40 98.77% (未定義四面体中心の判定に隣接済み立体中心を混ぜる)
    //            → I41 98.84% (プロトン化対象の絞り込み・金属ロック・アミジニウムリレー)
    //            → I42 98.92% (ビニロガスなアミジン-カルボニル互変異性の橋渡し・
    //               ホスホルアミド酸の二級 N 端点)
    //            → I43 99.04% (P(V) の未定義四面体中心・多成分 /m 層の連結方式)
    //            → I44 99.07% (式が同じ異性体成分どうしの順序タイブレーク)
    //            → I46 99.10% (/h の N* 圧縮に成分の重原子数一致を追加、
    //               IUPAC-InChI/InChI 公式 C ソース `ichiprt3.c` で確認)
    //            → I47 99.11% (橋頭ヘテロ原子の毒判定に真の縮環要件を追加)
    //            → I48 99.18% (成分順序タイブレークを接続層の文字列比較から
    //               整数配列比較に修正、`ichimake.c::CompINChI2` で確認)
    //            → I49 99.19% (Hückel 判定: π 電子 2 個での芳香族認定を
    //               3 員環限定に修正、4 員環ジオン類の誤った可動 H 検出を解消)
    //            → I50 99.21% (アミド型 N⁻ をプロトン化対象に追加、単純
    //               アルキルアミン型 N⁻ は引き続き対象外)。
    //            → I51 99.22% (チオホスホルアミド `R-NH-P(=S)` の P=S 硫黄を
    //               星型の受容体端点に追加、`n_double_o_ok` が O 限定
    //               だったのを P 中心に限り O/S どちらでも可に修正)。
    //            → I53 99.34% (成分並び順の主キーを `ichimake.c::
    //               CompareHillFormulasNoH` 準拠のトークン単位比較に全面
    //               書き直し。旧規則「炭素数降順」は例からの逆算による
    //               過適合だった)。
    //            → I54 99.35% (ペルオキシ酸アニオン `C(=O)O[O-]` の末端 O
    //               をプロトン化対象に追加。隣接原子が酸性中心でなくても、
    //               隣が単結合 O でその先に酸性中心があれば許可)。
    //            → I55 99.38% (`[H-]` のような未マージの孤立 H 原子が
    //               `build_molecule_graph` の出力で重原子より前に来てしまい、
    //               `neutralize` が前提とする「重原子が先頭に連続」という
    //               不変条件を壊して電荷正規化そのものが丸ごとスキップされて
    //               いた。孤立 H を重原子の後ろへ回すよう並び替え、電荷を
    //               保持したまま再構成し、/q 層にもその電荷を反映)。
    //            → I56 99.42% (成分順序タイブレークの接続表を「全隣接を
    //               フラット化」から実 InChI の `nConnTable` 相当 (各原子の
    //               正準番号自身 → 自分より小さい番号の隣接 (後退辺) のみを
    //               昇順) に修正。`ichicano.c::UpdateFullLinearCT` の
    //               `CT_ATOMID_IS_CURRANK` モードで確認)。
    //            → I57 99.43% (`is_acidic_oh` の「隣が炭素でなければ酸性」
    //               判定がホウ酸型 B-OH まで拾ってしまい、無関係な四級
    //               アンモニウムの +1 を中和するのにボロン酸からプロトンを
    //               奪っていたのを修正。B は非炭素の例外から除外)。
    //            → I58 99.44% (単純アルコキシド `CC[O-]` 等が孤立イオンでも
    //               プロトン化されていたのを修正。O⁻ の隣が sp3 アルキル
    //               炭素なら非酸性としてプロトン化対象から除外し、芳香族
    //               (フェノラート) と二重結合を持つ炭素 (エノラート・
    //               カルボキシラート) だけを引き続き対象とする。S/Se/Te は
    //               対象外 (チオラートは sp3 炭素隣接でもプロトン化される)。
    //            → I59 99.45% (スルホニル S(=O)(=O) のように実在二重結合を
    //               2 本同時に持つ過原子価中心が、ブロッサム法の `matched[]`
    //               に非対称に (後勝ちで) 登録され、環外酸素まで可動 H 群に
    //               誤併合していたのを修正。両端の実在二重結合数がそれぞれ
    //               ちょうど 1 本の結合だけを `matched[]` に登録する)。
    //            → I60 99.46% (未定義立体源性二重結合 `?` の判定で「環内は
    //               全部対象外」としていたのを、8 員未満の小環だけ対象外に
    //               縮小。マクロ環 (8 員以上) は実 InChI も E/Z を認識する)。
    //            → I61 99.47% (ヒドロキシルアミニウム型 `[NH3+]O` の N-O
    //               境界を陽イオン脱プロトン側にも適用。N の唯一の重原子
    //               隣接が O だけのときは脱プロトンせず `/q` のまま残す。
    //               第 2 パスの `is_acidic_oh` も、その O-H の隣が「それ
    //               自身が陽電荷を持つ N」なら対象外にする対称の修正)。
    //            → I62 99.49% (孤立した `[H+]` を、中性化で仮想的に足した
    //               プロトン (n_add) の実体として直接消費するよう修正。
    //               従来は `[H+].Cl -> C...` 型の入力で `[Cl-]` 側は
    //               `is_protonatable` 経由で正しく ClH になるが、明示的な
    //               `[H+]` 原子自体は別成分として残り `.H` + `/q+1/p-1` が
    //               余計に付いていた。消費した分だけ `n_add` を減らし、
    //               対応する孤立 H 原子を独立成分保持ループから除外)。
    // 残る不一致は RUST_INCHI_I29_PLAN.md を参照。
    assert!(acc >= 0.9948, "pubchem InChI accuracy {acc:.4} < 0.9948");
}
