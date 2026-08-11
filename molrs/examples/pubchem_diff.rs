//! PubChem 不一致の層別診断 (I29 系の作業用、非ゲート)。
//!
//! 使い方: cargo run --release --example pubchem_diff [-- --layer h] [--limit 20]
//!
//! `corpus/pubchem_inchi.jsonl.gz` を読み、molrs の出力と PubChem InChI を
//! 層ごとに比較して、どの層が食い違っているかで分類する。

use std::collections::BTreeMap;
use std::io::Read;
use std::path::PathBuf;

use flate2::read::GzDecoder;
use molrs::graph::build_molecule_graph;

struct Record {
    smiles: String,
    inchi: String,
}

fn load() -> Vec<Record> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../corpus/pubchem_inchi.jsonl.gz");
    let file =
        std::fs::File::open(&path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    let mut text = String::new();
    GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("gunzip");
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let v: serde_json::Value = serde_json::from_str(l).expect("json");
            Record {
                smiles: v["s"].as_str().unwrap().to_string(),
                inchi: v["inchi"].as_str().unwrap_or("").to_string(),
            }
        })
        .collect()
}

/// InChI 文字列を層プレフィックス (先頭の `/X` の X) -> 内容 に分割する。
/// 式層は先頭 (プレフィックスなし) を "F" とする。
fn split_layers(inchi: &str) -> BTreeMap<String, String> {
    let mut m = BTreeMap::new();
    let body = inchi.strip_prefix("InChI=1S/").or_else(|| inchi.strip_prefix("InChI=1/"));
    let Some(body) = body else { return m };
    let mut cur_key = "F".to_string();
    let mut cur_val = String::new();
    for part in body.split('/') {
        if let Some(rest) = part.strip_prefix(|c: char| c.is_ascii_lowercase()) {
            let key = part.chars().next().unwrap().to_string();
            if cur_key != key || !cur_val.is_empty() {
                m.insert(cur_key.clone(), cur_val.clone());
            }
            cur_key = key;
            cur_val = rest.to_string();
        } else {
            if !cur_val.is_empty() {
                cur_val.push('/');
            }
            cur_val.push_str(part);
        }
    }
    m.insert(cur_key, cur_val);
    m
}

fn try_inchi(smiles: &str) -> Option<String> {
    let g = build_molecule_graph(smiles).ok()?;
    molrs::inchi::to_inchi(&g).ok()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut layer_filter: Option<String> = None;
    let mut limit = 20usize;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--layer" => {
                layer_filter = args.get(i + 1).cloned();
                i += 2;
            }
            "--limit" => {
                limit = args.get(i + 1).and_then(|s| s.parse().ok()).unwrap_or(20);
                i += 2;
            }
            _ => {
                i += 1;
            }
        }
    }

    let recs = load();
    let mut category_counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut examples: BTreeMap<String, Vec<(String, String, String)>> = BTreeMap::new();
    let mut total_mismatch = 0usize;

    for r in &recs {
        if r.inchi.is_empty() {
            continue;
        }
        let Some(got) = try_inchi(&r.smiles) else {
            let cat = "no-output".to_string();
            *category_counts.entry(cat.clone()).or_default() += 1;
            examples.entry(cat).or_default().push((
                r.smiles.clone(),
                "<none>".to_string(),
                r.inchi.clone(),
            ));
            total_mismatch += 1;
            continue;
        };
        if got == r.inchi {
            continue;
        }
        total_mismatch += 1;
        let want_layers = split_layers(&r.inchi);
        let got_layers = split_layers(&got);
        let mut diff_keys: Vec<String> = Vec::new();
        let all_keys: std::collections::BTreeSet<&String> =
            want_layers.keys().chain(got_layers.keys()).collect();
        for k in all_keys {
            if want_layers.get(k) != got_layers.get(k) {
                diff_keys.push(k.clone());
            }
        }
        diff_keys.sort();
        let cat = diff_keys.join("");
        *category_counts.entry(cat.clone()).or_default() += 1;
        if let Some(f) = &layer_filter {
            if &cat != f {
                continue;
            }
        }
        let ex = examples.entry(cat).or_default();
        if ex.len() < limit {
            ex.push((r.smiles.clone(), got, r.inchi.clone()));
        }
    }

    println!("total mismatches: {total_mismatch}");
    println!("--- category counts (diff layer set) ---");
    let mut sorted: Vec<_> = category_counts.iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(a.1));
    for (cat, n) in &sorted {
        println!("{n:4}  {cat}");
    }

    if let Some(f) = &layer_filter {
        println!("\n--- examples for category {f:?} ---");
        if let Some(ex) = examples.get(f) {
            for (smi, got, want) in ex {
                println!("SMI : {smi}");
                println!("got : {got}");
                println!("want: {want}");
                println!();
            }
        }
    }
}
