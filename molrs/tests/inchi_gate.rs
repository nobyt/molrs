//! InChI 差分ゲート (RUST_INCHI_PLAN.md §検証)。
//!
//! フィクスチャ corpus/inchi_dump.jsonl.gz (smiles2iupac tools/dump_inchi.py で
//! RDKit から採取) に対して層ごとに一致を検査する。実装の進行に応じて検査
//! 項目を増やす:
//! - I2: 式層 (Hill) の一致率
//! - I3: 正準番号 (AuxInfo /N:) の一致率
//! - I4/I5: フル InChI 文字列・InChIKey の一致率 (v1 適用範囲で 100%)

use std::io::Read;
use std::path::PathBuf;

use flate2::read::GzDecoder;
use molrs::graph::build_molecule_graph;

struct Record {
    smiles: String,
    formula: String,
    numbering: Vec<Vec<usize>>,
}

fn load_fixture() -> Vec<Record> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../corpus/inchi_dump.jsonl.gz");
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
            let numbering = v["n"]
                .as_array()
                .map(|comps| {
                    comps
                        .iter()
                        .map(|c| {
                            c.as_array()
                                .unwrap()
                                .iter()
                                .map(|x| x.as_u64().unwrap() as usize)
                                .collect()
                        })
                        .collect()
                })
                .unwrap_or_default();
            Record {
                smiles: v["s"].as_str().unwrap().to_string(),
                formula: v["formula"].as_str().unwrap_or("").to_string(),
                numbering,
            }
        })
        .collect()
}

#[test]
fn formula_layer_matches_rdkit() {
    let recs = load_fixture();
    let mut n = 0usize;
    let mut ok = 0usize;
    let mut mism: Vec<String> = Vec::new();
    for r in &recs {
        if r.formula.is_empty() {
            continue;
        }
        let Ok(g) = build_molecule_graph(&r.smiles) else {
            continue;
        };
        n += 1;
        let got = molrs::inchi::formula(&g);
        if got == r.formula {
            ok += 1;
        } else if mism.len() < 25 {
            mism.push(format!("{}: got {got}, want {}", r.smiles, r.formula));
        }
    }
    let rate = ok as f64 / n as f64;
    println!("formula layer: {ok}/{n} match ({:.2}%)", rate * 100.0);
    for m in &mism {
        println!("  MISMATCH {m}");
    }
    // v1: 電荷正規化 (プロトン移動) を要さない分子は一致するはず。
    // 正規化差の分子があるため閾値は段階的に引き上げる (I4 で normalize 実装後)。
    assert!(rate >= 0.90, "formula match rate {rate:.4} < 0.90");
}

#[test]
fn canonical_numbering_matches_auxinfo() {
    let recs = load_fixture();
    let mut n = 0usize;
    let mut ok = 0usize;
    let mut mism: Vec<String> = Vec::new();
    for r in &recs {
        if r.numbering.is_empty() {
            continue;
        }
        let Ok(g) = build_molecule_graph(&r.smiles) else {
            continue;
        };
        n += 1;
        let got: Vec<Vec<usize>> = molrs::inchi::number::canonical_numbering(&g)
            .iter()
            .map(|comp| comp.iter().map(|&i| i + 1).collect())
            .collect();
        if got == r.numbering {
            ok += 1;
        } else if mism.len() < 20 {
            mism.push(format!("{}: got {got:?}, want {:?}", r.smiles, r.numbering));
        }
    }
    let rate = ok as f64 / n as f64;
    println!("canonical numbering: {ok}/{n} match ({:.2}%)", rate * 100.0);
    for m in &mism {
        println!("  MISMATCH {m}");
    }
    // I3 段階の到達点 = 70.96%。残差は (a) InChI の正準タイブレーク基準
    // (辺シグネチャ最小化とは異なる) と (b) 可動 H 対称化の未実装 (I4)。
    // 閾値は退行検知用 (現状より少し下)。忠実度向上は I3→I4 の継続課題。
    assert!(rate >= 0.68, "numbering match rate {rate:.4} < 0.68");
}
