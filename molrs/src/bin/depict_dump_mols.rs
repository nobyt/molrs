//! D13 開発用: stdin の SMILES (1 行 1 分子) を 2D レイアウトし、
//! くさびコード付き 2D MOL ブロックを JSONL で stdout に出力する。
//! smiles2iupac リポジトリの tools/check_depict_stereo.py が RDKit で
//! 再認識し、立体の round-trip を検証する。

use std::io::{BufRead, Write};

use molrs::depict::{compute_coords_2d, to_mol_block_2d, LayoutParams};

fn main() {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    let params = LayoutParams::default();

    for line in stdin.lock().lines() {
        let line = line.expect("stdin");
        let smiles = line.trim();
        if smiles.is_empty() {
            continue;
        }
        match molrs::graph::build_molecule_graph(smiles) {
            Ok(g) => match compute_coords_2d(&g, &params) {
                Ok(c) => {
                    let mol = to_mol_block_2d(&g, &c, smiles);
                    writeln!(out, "{{\"s\":{smiles:?},\"mol\":{mol:?}}}").unwrap();
                }
                Err(e) => {
                    writeln!(out, "{{\"s\":{smiles:?},\"error\":{:?}}}", e.to_string()).unwrap();
                }
            },
            Err(e) => {
                writeln!(out, "{{\"s\":{smiles:?},\"error\":{:?}}}", e.to_string()).unwrap();
            }
        }
    }
}
