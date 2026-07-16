//! S1.4 開発用: コーパス全分子の ring_atom_sets を JSONL で出力する。
//! tools/compare_rings.py が RDKit ダンプと突き合わせる。

use std::io::{BufRead, Write};

fn main() {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    for line in stdin.lock().lines() {
        let line = line.expect("stdin");
        let smiles = line.trim();
        if smiles.is_empty() {
            continue;
        }
        match molrs::graph::build_molecule_graph(smiles) {
            Ok(g) => {
                let rings: Vec<String> = g
                    .ring_atom_sets
                    .iter()
                    .map(|r| {
                        let items: Vec<String> = r.iter().map(|a| a.to_string()).collect();
                        format!("[{}]", items.join(","))
                    })
                    .collect();
                writeln!(out, "{{\"s\":{smiles:?},\"r\":[{}]}}", rings.join(",")).unwrap();
            }
            Err(e) => {
                writeln!(out, "{{\"s\":{smiles:?},\"error\":{:?}}}", e.to_string()).unwrap();
            }
        }
    }
}
