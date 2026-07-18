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
    inchi: String,
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
                inchi: v["inchi"].as_str().unwrap_or("").to_string(),
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
    // I3 到達点 = 99.29% (7400/7453)。残差の内訳:
    // - 有機金属 (金属-C 切断で多成分化) — I4 normalize で対応
    // - イソシアニド等の電荷正規化差 — I4
    // - 環内 N 互変異性 (ベンズイミダゾリン等) — v2
    // - 立体依存タイブレーク (E/Z 対称分子) — 骨格層の文字列には影響しない
    assert!(rate >= 0.985, "numbering match rate {rate:.4} < 0.985");
}

#[test]
fn full_inchi_matches_rdkit_where_produced() {
    // v1 が文字列を生成する範囲 (中性・単一成分) でフル InChI が一致すること。
    // カバレッジ (生成できた割合) も報告する。
    let recs = load_fixture();
    let mut produced = 0usize;
    let mut ok = 0usize;
    let mut total_neutral_single = 0usize;
    let mut mism: Vec<String> = Vec::new();
    for r in &recs {
        if r.inchi.is_empty() {
            continue;
        }
        let Ok(g) = build_molecule_graph(&r.smiles) else {
            continue;
        };
        // v1 適用範囲: 立体層 (/b /t /m /s) と同位体層 (/i) を含まない分子
        let has_stereo_or_iso = ["/b", "/t", "/m", "/s", "/i"]
            .iter()
            .any(|tag| r.inchi.contains(tag));
        let charged = g.atoms.iter().any(|a| a.formal_charge != 0);
        let multi = r.inchi[9..].split('/').next().unwrap().contains('.');
        if !charged && !multi && !has_stereo_or_iso {
            total_neutral_single += 1;
        }
        if has_stereo_or_iso {
            continue; // v2
        }
        if let Ok(got) = molrs::inchi::to_inchi(&g) {
            produced += 1;
            if got == r.inchi {
                ok += 1;
            } else if mism.len() < 30 {
                mism.push(format!("{} | got {got} | want {}", r.smiles, r.inchi));
            }
        }
    }
    let acc = ok as f64 / produced.max(1) as f64;
    let cov = produced as f64 / total_neutral_single.max(1) as f64;
    println!(
        "full InChI: produced {produced}, {ok} exact ({:.2}%); coverage of neutral-single {produced}/{total_neutral_single} ({:.2}%)",
        acc * 100.0,
        cov * 100.0
    );
    for m in &mism {
        println!("  MISMATCH {m}");
    }
    // I4 到達点 = 74.17% (v1 適用範囲 = 中性・単一成分・立体/同位体なし)。
    // 残差の主因: 可動 H 認識の未整備 (アミド/アミジン/ラクタム/環 N-H
    // 互変異性 — 逐次拡張) と一部の縮合環 c 層直列化。閾値は退行検知用。
    assert!(acc >= 0.72, "full InChI accuracy {acc:.4} < 0.72");
}
