//! stdin の SMILES (1 行 1 分子) を InChI / InChIKey に変換し JSONL 出力 (I5)。
//!
//! 使い方:
//!   jq -r .smiles ../corpus/corpus.jsonl | cargo run --release --bin inchi
//!
//! v1 適用範囲外 (電荷・多成分・立体・同位体・有機金属) は "error" を返す。

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
        match molrs::inchi::inchi_of(smiles) {
            Ok(inchi) => {
                let key = molrs::inchi::inchi_key_from_string(&inchi);
                writeln!(
                    out,
                    "{{\"s\":{smiles:?},\"inchi\":{inchi:?},\"key\":{key:?}}}"
                )
                .unwrap();
            }
            Err(e) => {
                writeln!(out, "{{\"s\":{smiles:?},\"error\":{:?}}}", e.to_string()).unwrap();
            }
        }
    }
}
