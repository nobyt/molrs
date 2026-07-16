//! 開発用スモーク: コーパス全体での compute_coords_2d 成功率を測る。
//! (本ゲートは D12 の depict_gate.rs — これは進捗確認用)

use std::path::PathBuf;

use molrs::depict::{compute_coords_2d, DepictError, LayoutParams};

#[test]
#[ignore = "dev smoke: cargo test --test depict_smoke -- --ignored --nocapture"]
fn corpus_smoke() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../corpus/corpus.jsonl");
    let text = std::fs::read_to_string(&path).expect("corpus");
    let mut n = 0usize;
    let mut ok = 0usize;
    let mut unsupported = 0usize;
    let mut failed: Vec<(String, String)> = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec: serde_json::Value = serde_json::from_str(line).expect("json");
        let smiles = rec["smiles"].as_str().expect("smiles");
        let Ok(g) = molrs::graph::build_molecule_graph(smiles) else {
            continue;
        };
        n += 1;
        match compute_coords_2d(&g, &LayoutParams::default()) {
            Ok(_) => ok += 1,
            Err(DepictError::Unsupported(_)) => unsupported += 1,
            Err(e) => {
                if failed.len() < 15 {
                    failed.push((smiles.to_string(), e.to_string()));
                }
            }
        }
    }
    println!("depict smoke: {ok}/{n} ok, {unsupported} unsupported");
    for (s, e) in &failed {
        println!("  FAIL {s}: {e}");
    }
}
