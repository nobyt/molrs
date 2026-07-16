//! S1.5 完了条件ゲート:
//! 1. コーパス全分子で canon(smiles) == canon(reparse(canon(smiles))) (冪等性)
//! 2. 同一分子の異表記グループ (corpus/canon_pairs.jsonl) が同じ正規形に写ること

use std::path::PathBuf;

fn repo_path(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join(rel)
}

fn canon(smiles: &str) -> Result<String, String> {
    molrs::graph::build_molecule_graph(smiles)
        .map(|g| molrs::canon::to_canonical_smiles(&g))
        .map_err(|e| e.to_string())
}

#[test]
fn corpus_idempotency() {
    let text = std::fs::read_to_string(repo_path("corpus/corpus.jsonl")).expect("corpus");
    let mut n = 0usize;
    let mut failures: Vec<String> = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let row: serde_json::Value = serde_json::from_str(line).expect("json");
        let smiles = row["smiles"].as_str().expect("smiles");
        n += 1;
        let c1 = match canon(smiles) {
            Ok(c) => c,
            Err(e) => {
                if failures.len() < 20 {
                    failures.push(format!("{smiles}: canon error: {e}"));
                }
                continue;
            }
        };
        match canon(&c1) {
            Ok(c2) if c2 == c1 => {}
            Ok(c2) => {
                if failures.len() < 20 {
                    failures.push(format!("{smiles}: {c1} -> {c2}"));
                }
            }
            Err(e) => {
                if failures.len() < 20 {
                    failures.push(format!("{smiles}: reparse of {c1:?} failed: {e}"));
                }
            }
        }
    }
    println!("canon idempotency: {n} molecules");
    assert!(
        failures.is_empty(),
        "{} idempotency failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn variant_groups_agree() {
    let text = std::fs::read_to_string(repo_path("corpus/canon_pairs.jsonl")).expect("pairs");
    let mut n_groups = 0usize;
    let mut failures: Vec<String> = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let row: serde_json::Value = serde_json::from_str(line).expect("json");
        let variants: Vec<&str> = row["variants"]
            .as_array()
            .expect("variants")
            .iter()
            .map(|v| v.as_str().expect("str"))
            .collect();
        n_groups += 1;
        let canons: Vec<Result<String, String>> = variants.iter().map(|s| canon(s)).collect();
        let first = &canons[0];
        for (v, c) in variants.iter().zip(&canons) {
            if (c != first || c.is_err()) && failures.len() < 20 {
                failures.push(format!(
                    "group {:?}: {v} -> {c:?} != {first:?}",
                    variants[0]
                ));
            }
        }
    }
    println!("canon variant groups: {n_groups}");
    assert!(
        failures.is_empty(),
        "{} variant mismatches:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
