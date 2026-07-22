//! 命名ゲート: molrs の smiles_to_iupac が Python 版 (names.jsonl.gz) と
//! 一致する分子の割合を測る。非環式・基本官能基の範囲で一致を目標とする。

use std::io::Read;
use std::path::PathBuf;

use flate2::read::GzDecoder;

#[test]
fn acyclic_names_match_python() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../corpus/names.jsonl.gz");
    let file = std::fs::File::open(&path)
        .unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    let mut text = String::new();
    GzDecoder::new(file).read_to_string(&mut text).expect("gunzip");

    let mut n = 0usize; // Python が名付けられた分子
    let mut attempted = 0usize; // molrs が名前を出した
    let mut ok = 0usize; // 一致
    let mut mism: Vec<String> = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let v: serde_json::Value = serde_json::from_str(line).expect("json");
        let smiles = v["s"].as_str().unwrap();
        let want = v["name"].as_str().unwrap_or("");
        if want.is_empty() {
            continue;
        }
        n += 1;
        if let Ok(got) = smiles2iupac::smiles_to_iupac(smiles) {
            attempted += 1;
            if got == want {
                ok += 1;
            } else if mism.len() < 25 {
                mism.push(format!("{smiles}: got {got} | want {want}"));
            }
        }
    }
    let cov = attempted as f64 / n as f64;
    let acc = ok as f64 / attempted.max(1) as f64;
    println!(
        "names: {ok}/{attempted} exact of {attempted} attempted ({:.1}% acc); coverage {attempted}/{n} ({:.1}%)",
        acc * 100.0,
        cov * 100.0
    );
    for m in &mism {
        println!("  MISMATCH {m}");
    }
    // 名前を出した分子は高精度で一致すべき (誤名は出さない方針)
    assert!(acc >= 0.99, "name accuracy {acc:.4} < 0.99");
    // カバレッジ (何割の分子を名付けられたか) — 段階的に上げる
    assert!(attempted >= 420, "too few molecules named: {attempted}");
}
