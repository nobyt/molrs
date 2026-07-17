//! D12: 2D 描画の決定的性質ゲート (RUST_2D_PLAN.md §検証)。
//!
//! コーパスの決定的 1/10 サンプル (DEPICT_GATE_FULL=1 で全数) に対して:
//! 1. レイアウト成功率 100%
//! 2. 非環結合の長さ = 1.0 ± 2% (100%)
//! 3. 環結合の長さ = 1.0 ± 2% (≥ 99% — 橋かけ系の内部結合のみ逸脱可)
//! 4. 非環結合の 30° 量子化率 ≥ 95%
//! 5. 重なり違反 (非結合可視原子対 < 0.5L) は既知例外リストのみ
//! 6. E/Z 幾何再導出 100% + くさび→CIP round-trip 100% (verify_stereo_2d)
//! 7. SVG well-formed (全サンプル)
//! 8. 決定性 (再計算で SVG バイト一致)
//!
//! 補助立体セット corpus/depict_stereo.jsonl があれば 6 をそこにも適用する。

use std::path::PathBuf;

use molrs::depict::{compute_coords_2d, to_svg, verify_stereo_2d, LayoutParams, Style};
use molrs::graph::{build_molecule_graph, MoleculeGraph};

/// 既知の重なり例外 (かご型は 2D で本質的に交差・接近する)。
const CLASH_EXCEPTIONS: &[&str] = &[
    "C12C3C4C1C5C4C3C25", // キュバン
];

fn corpus_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("../corpus/{name}"))
}

fn is_ring_bond(g: &MoleculeGraph, i: usize, j: usize) -> bool {
    g.ring_atom_sets.iter().any(|ring| {
        let n = ring.len();
        (0..n).any(|k| {
            let (a, b) = (ring[k], ring[(k + 1) % n]);
            (a == i && b == j) || (a == j && b == i)
        })
    })
}

fn assert_well_formed(svg: &str) -> Result<(), String> {
    let mut stack: Vec<String> = Vec::new();
    let bytes = svg.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] != b'<' {
            i += 1;
            continue;
        }
        let end = svg[i..].find('>').map(|e| i + e).ok_or("unclosed '<'")?;
        let tag = &svg[i + 1..end];
        if !tag.matches('"').count().is_multiple_of(2) {
            return Err(format!("unbalanced quotes in <{tag}>"));
        }
        if let Some(name) = tag.strip_prefix('/') {
            let open = stack.pop().ok_or_else(|| format!("stray </{name}>"))?;
            if open != name {
                return Err("mismatched close tag".into());
            }
        } else if !tag.ends_with('/') && !tag.starts_with('?') && !tag.starts_with('!') {
            let name: String = tag
                .chars()
                .take_while(|c| !c.is_whitespace() && *c != '>')
                .collect();
            stack.push(name);
        }
        i = end + 1;
    }
    if !stack.is_empty() {
        return Err(format!("unclosed tags: {stack:?}"));
    }
    Ok(())
}

#[test]
fn depict_property_gate() {
    let full = std::env::var("DEPICT_GATE_FULL").is_ok();
    let text = std::fs::read_to_string(corpus_path("corpus.jsonl")).expect("corpus");
    let mut smiles_list: Vec<String> = Vec::new();
    for (k, line) in text.lines().filter(|l| !l.trim().is_empty()).enumerate() {
        if !full && k % 10 != 0 {
            continue;
        }
        let rec: serde_json::Value = serde_json::from_str(line).expect("json");
        smiles_list.push(rec["smiles"].as_str().expect("smiles").to_string());
    }
    // 補助立体セット (あれば)
    if let Ok(extra) = std::fs::read_to_string(corpus_path("depict_stereo.jsonl")) {
        for line in extra.lines().filter(|l| !l.trim().is_empty()) {
            let rec: serde_json::Value = serde_json::from_str(line).expect("json");
            smiles_list.push(rec["smiles"].as_str().expect("smiles").to_string());
        }
    }

    let style = Style::acs_1996();
    let params = LayoutParams::default();

    let mut n = 0usize;
    let mut n_chain_bonds = 0usize;
    let mut n_chain_quantized = 0usize;
    let mut n_ring_bonds = 0usize;
    let mut n_ring_bond_ok = 0usize;
    let mut n_stereo_mols = 0usize;
    let mut errors: Vec<String> = Vec::new();
    let mut err = |msg: String, errors: &mut Vec<String>| {
        if errors.len() < 20 {
            errors.push(msg);
        }
    };

    for smiles in &smiles_list {
        let Ok(g) = build_molecule_graph(smiles) else {
            continue;
        };
        n += 1;
        // 1. 成功率
        let c = match compute_coords_2d(&g, &params) {
            Ok(c) => c,
            Err(e) => {
                err(format!("{smiles}: layout failed: {e}"), &mut errors);
                continue;
            }
        };
        let vis: Vec<usize> = (0..g.atoms.len()).filter(|&i| !c.hidden[i]).collect();

        // 2-4. 結合長・量子化
        for b in &g.bonds {
            let (i, j) = (b.begin_idx, b.end_idx);
            if c.hidden[i] || c.hidden[j] {
                continue;
            }
            let d = c.pos[i].distance(c.pos[j]);
            let ok_len = (d - 1.0).abs() <= 0.02;
            if is_ring_bond(&g, i, j) {
                n_ring_bonds += 1;
                if ok_len {
                    n_ring_bond_ok += 1;
                }
            } else {
                n_chain_bonds += 1;
                if !ok_len {
                    err(
                        format!("{smiles}: chain bond {i}-{j} length {d:.4}"),
                        &mut errors,
                    );
                }
                let ang = (c.pos[j] - c.pos[i]).angle();
                let steps = ang / (std::f64::consts::PI / 6.0);
                if (steps - steps.round()).abs() < 1e-6 {
                    n_chain_quantized += 1;
                }
            }
        }

        // 5. 重なり
        let excepted = CLASH_EXCEPTIONS.contains(&smiles.as_str());
        if !excepted {
            'outer: for (k, &i) in vis.iter().enumerate() {
                for &j in &vis[k + 1..] {
                    if g.adjacency[i].contains(&j) {
                        continue;
                    }
                    let d = c.pos[i].distance(c.pos[j]);
                    if d < 0.5 {
                        err(
                            format!("{smiles}: clash atoms {i},{j} d={d:.3}"),
                            &mut errors,
                        );
                        break 'outer;
                    }
                }
            }
        }

        // 6. 立体 round-trip
        let has_stereo = g.atoms.iter().any(|a| a.chiral_tag.is_some())
            || g.bonds.iter().any(|b| b.stereo.is_some());
        if has_stereo {
            n_stereo_mols += 1;
            let failures = verify_stereo_2d(&g, &c);
            if !failures.is_empty() {
                err(format!("{smiles}: stereo {failures:?}"), &mut errors);
            }
        }

        // 7-8. SVG well-formed + 決定性
        let svg = to_svg(&g, &c, &style);
        if let Err(e) = assert_well_formed(&svg) {
            err(format!("{smiles}: svg: {e}"), &mut errors);
        }
        let c2 = compute_coords_2d(&g, &params).expect("relayout");
        let svg2 = to_svg(&g, &c2, &style);
        if svg != svg2 {
            err(format!("{smiles}: non-deterministic"), &mut errors);
        }
    }

    let quant_rate = n_chain_quantized as f64 / n_chain_bonds.max(1) as f64;
    let ring_ok_rate = n_ring_bond_ok as f64 / n_ring_bonds.max(1) as f64;
    println!(
        "depict gate: {n} molecules ({} stereo), chain quantization {:.2}%, ring bond ok {:.2}%",
        n_stereo_mols,
        quant_rate * 100.0,
        ring_ok_rate * 100.0
    );
    assert!(
        errors.is_empty(),
        "{} gate violations (first 20):\n{}",
        errors.len(),
        errors.join("\n")
    );
    assert!(
        quant_rate >= 0.95,
        "chain bond 30° quantization {quant_rate:.4} < 0.95"
    );
    assert!(
        ring_ok_rate >= 0.99,
        "ring bond length ok rate {ring_ok_rate:.4} < 0.99"
    );
    assert!(n_stereo_mols > 0, "no stereo molecules sampled");
}
