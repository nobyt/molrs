//! コーパス駆動回帰テストハーネス (RUST_PORT_PLAN.md S0.2)。
//!
//! `corpus/corpus.jsonl` の全 (smiles, expected) ペアを実行し、
//! `tests/expected_pass.txt` (合格リスト、1 行 1 SMILES) と突き合わせる:
//!
//! - 合格リスト掲載ケースが不合格になったら **リグレッションとして fail**
//! - 未掲載ケースが新たに合格したら件数を報告 (リストへの追記を促す)
//!
//! 合格リストは単調増加で管理する。更新は:
//!   UPDATE_EXPECTED_PASS=1 cargo test --test corpus

use serde::Deserialize;
use std::collections::BTreeSet;
use std::path::PathBuf;

#[derive(Deserialize)]
struct Case {
    smiles: String,
    expected: String,
    #[allow(dead_code)]
    phases: Vec<i64>,
    excluded: bool,
}

fn repo_root() -> PathBuf {
    // rust/smiles2iupac/ → リポジトリルート
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .canonicalize()
        .expect("repo root")
}

fn load_corpus() -> Vec<Case> {
    let path = repo_root().join("corpus/corpus.jsonl");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).expect("corpus line"))
        .collect()
}

fn expected_pass_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/expected_pass.txt")
}

fn load_expected_pass() -> BTreeSet<String> {
    match std::fs::read_to_string(expected_pass_path()) {
        Ok(text) => text
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
            .map(String::from)
            .collect(),
        Err(_) => BTreeSet::new(),
    }
}

#[test]
fn corpus_regression() {
    let corpus = load_corpus();
    let expected_pass = load_expected_pass();

    let mut passing = BTreeSet::new();
    let mut regressions: Vec<String> = Vec::new();

    for case in corpus.iter().filter(|c| !c.excluded) {
        let pass = matches!(
            smiles2iupac::smiles_to_iupac(&case.smiles),
            Ok(ref name) if *name == case.expected
        );
        if pass {
            passing.insert(case.smiles.clone());
        } else if expected_pass.contains(&case.smiles) {
            let got = match smiles2iupac::smiles_to_iupac(&case.smiles) {
                Ok(name) => format!("got {name:?}"),
                Err(e) => format!("error: {e}"),
            };
            regressions.push(format!(
                "  {} => want {:?}, {}",
                case.smiles, case.expected, got
            ));
        }
    }

    let new_passes: Vec<&String> = passing.difference(&expected_pass).collect();

    println!(
        "corpus: {} cases, {} passing ({} in expected list, {} new)",
        corpus.len(),
        passing.len(),
        expected_pass.len(),
        new_passes.len()
    );

    if std::env::var("UPDATE_EXPECTED_PASS").is_ok() {
        let mut merged: BTreeSet<String> = expected_pass.clone();
        merged.extend(passing.iter().cloned());
        let mut body = String::new();
        for s in &merged {
            body.push_str(s);
            body.push('\n');
        }
        std::fs::write(expected_pass_path(), body).expect("write expected_pass.txt");
        println!("expected_pass.txt updated: {} entries", merged.len());
    } else if !new_passes.is_empty() {
        println!(
            "note: {} newly passing cases not in expected_pass.txt — run with UPDATE_EXPECTED_PASS=1 to record them",
            new_passes.len()
        );
    }

    assert!(
        regressions.is_empty(),
        "{} regression(s) from expected_pass.txt:\n{}",
        regressions.len(),
        regressions.join("\n")
    );
}
