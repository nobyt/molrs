//! S1.1 完了条件: コーパス全 SMILES がパースエラーなく AST 化できること。

use std::path::PathBuf;

#[test]
fn all_corpus_smiles_parse() {
    let corpus = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../corpus/corpus.jsonl");
    let text = std::fs::read_to_string(&corpus)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", corpus.display()));

    let mut total = 0usize;
    let mut failures: Vec<String> = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let row: serde_json::Value = serde_json::from_str(line).expect("corpus line");
        let smiles = row["smiles"].as_str().expect("smiles field");
        total += 1;
        if let Err(e) = molrs::smiles::parse_smiles(smiles) {
            failures.push(format!("  {smiles}: {e}"));
        }
    }

    println!("corpus parse: {}/{} ok", total - failures.len(), total);
    assert!(
        failures.is_empty(),
        "{} corpus SMILES failed to parse:\n{}",
        failures.len(),
        failures[..failures.len().min(30)].join("\n")
    );
}
