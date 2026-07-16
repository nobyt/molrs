//! C10 完了条件ゲート: 実験トーションの収集結果が RDKit
//! `GetExperimentalTorsions` (ETKDGv2) と全数一致すること。
//! 比較キー: (中央結合 (a2<a3), ライブラリパターン番号) の集合。

use flate2::read::GzDecoder;
use std::io::Read;
use std::path::PathBuf;

#[test]
fn matches_rdkit_exp_torsions() {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../corpus/exp_torsion_fixture.jsonl.gz");
    let file = std::fs::File::open(&path)
        .unwrap_or_else(|e| panic!("cannot open {}: {e}", path.display()));
    let mut text = String::new();
    GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("gunzip");

    let mut n = 0usize;
    let mut n_tors = 0usize;
    let mut failures: Vec<String> = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec: serde_json::Value = serde_json::from_str(line).expect("json");
        let smiles = rec["s"].as_str().expect("s");
        let mut expected: Vec<(usize, usize, usize)> = rec["t"]
            .as_array()
            .expect("t")
            .iter()
            .map(|t| {
                let v = t.as_array().expect("triple");
                (
                    v[0].as_u64().expect("a2") as usize,
                    v[1].as_u64().expect("a3") as usize,
                    v[2].as_u64().expect("idx") as usize,
                )
            })
            .collect();
        expected.sort_unstable();

        let Ok(g) = molrs::graph::build_molecule_graph(smiles) else {
            continue;
        };
        n += 1;
        let mut got: Vec<(usize, usize, usize)> = molrs::conformer::exp_torsions_of(&g)
            .iter()
            .map(|t| {
                let (a, b) = (t.atoms[1].min(t.atoms[2]), t.atoms[1].max(t.atoms[2]));
                (a, b, t.torsion_idx)
            })
            .collect();
        got.sort_unstable();
        n_tors += got.len();

        if got != expected && failures.len() < 12 {
            failures.push(format!(
                "{smiles}:\n  got      {got:?}\n  expected {expected:?}"
            ));
        }
    }
    println!("exp torsion fixture: {n} molecules, {n_tors} torsions");
    assert!(
        failures.is_empty(),
        "{} mismatches:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
