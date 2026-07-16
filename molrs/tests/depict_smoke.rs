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
    let mut n_clash = 0usize;
    let mut clash_examples: Vec<String> = Vec::new();
    let mut failed: Vec<(String, String)> = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec: serde_json::Value = serde_json::from_str(line).expect("json");
        let smiles = rec["smiles"].as_str().expect("smiles");
        let Ok(g) = molrs::graph::build_molecule_graph(smiles) else {
            continue;
        };
        n += 1;
        match compute_coords_2d(&g, &LayoutParams::default()) {
            Ok(c) => {
                ok += 1;
                // 残存衝突の統計 (非結合可視原子対 < 0.5)
                let vis: Vec<usize> = (0..g.atoms.len()).filter(|&i| !c.hidden[i]).collect();
                let mut clash = false;
                for (k, &i) in vis.iter().enumerate() {
                    for &j in &vis[k + 1..] {
                        if !g.adjacency[i].contains(&j) && c.pos[i].distance(c.pos[j]) < 0.5 {
                            clash = true;
                        }
                    }
                }
                if clash {
                    n_clash += 1;
                    if clash_examples.len() < 8 {
                        clash_examples.push(smiles.to_string());
                    }
                }
            }
            Err(DepictError::Unsupported(_)) => unsupported += 1,
            Err(e) => {
                if failed.len() < 15 {
                    failed.push((smiles.to_string(), e.to_string()));
                }
            }
        }
    }
    println!("depict smoke: {ok}/{n} ok, {unsupported} unsupported, {n_clash} with residual clash");
    for s in &clash_examples {
        println!("  CLASH {s}");
    }
    for (s, e) in &failed {
        println!("  FAIL {s}: {e}");
    }
}
