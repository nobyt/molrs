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

/// 128 原子を超える分子は `build_molecule_graph` の環認識 (`rings.rs` の
/// `assert!`) でパニックするため (I29 で判明した既知のバグ)、パニックを
/// 捕まえて「不一致」として数える。SMILES パーサが弾く分 (超原子価ハロゲン)
/// も同様に不一致扱い。
fn try_inchi(smiles: &str) -> Option<String> {
    let s = smiles.to_string();
    std::panic::catch_unwind(|| {
        let g = build_molecule_graph(&s).ok()?;
        molrs::inchi::to_inchi(&g).ok()
    })
    .ok()
    .flatten()
}

#[test]
fn pubchem_full_inchi() {
    let recs = load();
    let prev = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {})); // 128 原子超のバックトレースを抑止
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
    std::panic::set_hook(prev);
    let acc = ok as f64 / n.max(1) as f64;
    println!("pubchem full InChI: {ok}/{n} exact ({:.2}%)", acc * 100.0);
    // I29 実測 94.66%。残る不一致の内訳は RUST_INCHI_I29_PLAN.md を参照
    // (立体 /t /m /s が最大、次いで電荷正規化・/b・128 原子超のパニック)。
    assert!(acc >= 0.945, "pubchem InChI accuracy {acc:.4} < 0.945");
}
