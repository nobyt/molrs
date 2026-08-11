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
    //               ホスホルアミド酸の二級 N 端点)。
    // 残る不一致は RUST_INCHI_I29_PLAN.md を参照。
    assert!(acc >= 0.9891, "pubchem InChI accuracy {acc:.4} < 0.9891");
}
