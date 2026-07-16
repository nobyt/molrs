//! S1.6 完了条件ゲート: corpus/substruct_fixture.jsonl.gz (RDKit
//! GetSubstructMatches の正解データ) と全ペアで一致すること (順序無視)。

use molrs::substructure::{substruct_matches, substruct_matches_smarts};
use flate2::read::GzDecoder;
use std::io::Read;
use std::path::PathBuf;

#[test]
fn matches_rdkit_fixture() {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../corpus/substruct_fixture.jsonl.gz");
    let file = std::fs::File::open(&path)
        .unwrap_or_else(|e| panic!("cannot open {}: {e}", path.display()));
    let mut text = String::new();
    GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("gunzip");

    let mut n = 0usize;
    let mut failures: Vec<String> = Vec::new();
    // クエリグラフはキャッシュ (同じクエリを大量のターゲットに適用するため)
    let mut query_cache: std::collections::HashMap<String, molrs::graph::MoleculeGraph> =
        std::collections::HashMap::new();

    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec: serde_json::Value = serde_json::from_str(line).expect("json");
        let tsmi = rec["t"].as_str().expect("t");
        let qsmi = rec["q"].as_str().expect("q");
        let mode = rec["mode"].as_str().expect("mode");
        let expected: Vec<Vec<usize>> = rec["m"]
            .as_array()
            .expect("m")
            .iter()
            .map(|t| {
                t.as_array()
                    .expect("tuple")
                    .iter()
                    .map(|v| v.as_u64().expect("idx") as usize)
                    .collect()
            })
            .collect();
        n += 1;

        let target = match molrs::graph::build_molecule_graph(tsmi) {
            Ok(g) => g,
            Err(e) => {
                failures.push(format!("{tsmi}: target build error: {e}"));
                continue;
            }
        };
        let mut got = if mode == "smarts" {
            substruct_matches_smarts(&target, qsmi).expect("pattern parses")
        } else {
            let query = query_cache.entry(qsmi.to_string()).or_insert_with(|| {
                molrs::graph::build_molecule_graph(qsmi).expect("query builds")
            });
            substruct_matches(&target, query)
        };
        got.sort();
        if got != expected && failures.len() < 15 {
            failures.push(format!(
                "t={tsmi} q={qsmi} mode={mode}: got {} matches {:?}... != expected {} {:?}...",
                got.len(),
                got.first(),
                expected.len(),
                expected.first(),
            ));
        }
    }
    println!("substruct fixture: {n} pairs compared");
    assert!(
        failures.is_empty(),
        "{} mismatches:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
