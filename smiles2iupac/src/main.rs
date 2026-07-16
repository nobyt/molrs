//! CLI: Python 版 `_cli.py` 互換 (本実装は S6.1)。
//! 単一 SMILES 引数、または `-` で stdin から 1 行 1 SMILES。

use std::io::BufRead;
use std::process::ExitCode;

fn convert_and_print(smiles: &str) -> bool {
    match smiles2iupac::smiles_to_iupac(smiles) {
        Ok(name) => {
            println!("{name}");
            true
        }
        Err(e) => {
            eprintln!("error: {e}");
            false
        }
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [input] = args.as_slice() else {
        eprintln!("usage: smiles2iupac <SMILES | ->");
        return ExitCode::from(2);
    };

    let mut ok = true;
    if input == "-" {
        for line in std::io::stdin().lock().lines() {
            let line = line.expect("stdin read error");
            let s = line.trim();
            if !s.is_empty() {
                ok &= convert_and_print(s);
            }
        }
    } else {
        ok = convert_and_print(input);
    }
    if ok {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}
