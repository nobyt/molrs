//! S3.1 完了ゲート: detect_groups の出力が Python 版 (コミット a01eccd) と
//! コーパス全分子で全数一致すること。
//! 比較対象: (group_type, atom_indices, priority) の列 (順序込み)。
//! 正解データ: tools/dump_functional_groups.py で採取。

use flate2::read::GzDecoder;
use std::io::Read;
use std::path::PathBuf;

#[test]
fn matches_python_detect_groups() {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../corpus/functional_groups.jsonl.gz");
    let file = std::fs::File::open(&path)
        .unwrap_or_else(|e| panic!("cannot open {}: {e}", path.display()));
    let mut text = String::new();
    GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("gunzip");

    let mut n = 0usize;
    let mut n_groups = 0usize;
    let mut failures: Vec<String> = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec: serde_json::Value = serde_json::from_str(line).expect("json");
        let smiles = rec["s"].as_str().expect("s");
        let expected: Vec<(String, Vec<usize>, i64)> = rec["g"]
            .as_array()
            .expect("g")
            .iter()
            .map(|t| {
                let v = t.as_array().expect("triple");
                (
                    v[0].as_str().expect("type").to_string(),
                    v[1].as_array()
                        .expect("indices")
                        .iter()
                        .map(|x| x.as_u64().expect("idx") as usize)
                        .collect(),
                    v[2].as_i64().expect("priority"),
                )
            })
            .collect();

        let g = molrs::graph::build_molecule_graph(smiles)
            .unwrap_or_else(|e| panic!("graph build failed for {smiles}: {e}"));
        n += 1;
        let got: Vec<(String, Vec<usize>, i64)> = smiles2iupac::functional_group::detect_groups(&g)
            .into_iter()
            .map(|fg| {
                (
                    fg.group_type.to_string(),
                    fg.atom_indices,
                    fg.priority as i64,
                )
            })
            .collect();
        n_groups += got.len();

        if got != expected && failures.len() < 10 {
            failures.push(format!(
                "{smiles}:\n  got      {got:?}\n  expected {expected:?}"
            ));
        }
    }
    println!("functional group fixture: {n} molecules, {n_groups} groups");
    assert!(
        failures.is_empty(),
        "{} mismatches (first 10):\n{}",
        failures.len(),
        failures.join("\n")
    );
}
