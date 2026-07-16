//! RDKit 互換性ゲート: corpus/rdkit_dump.jsonl.gz (Python 版 MoleculeGraph の
//! ダンプ) と Rust 版 build_molecule_graph の出力を全分子で突き合わせる。
//!
//! 比較段階 (RUST_PORT_PLAN.md):
//! - S1.2: 原子数・元素記号・形式電荷・環メンバーシップ・結合トポロジー
//! - S1.3: 芳香族フラグと全結合次数 (ケクレ化 + 芳香族認識)
//! - S1.4: ring_atom_sets (対称化 SSSR、環順・環内原子順も一致)
//! - S1.7: CIP コード (R/S) と結合 E/Z (レガシー AssignStereochemistry)

use flate2::read::GzDecoder;
use std::io::Read;
use std::path::PathBuf;

fn load_dump() -> Vec<serde_json::Value> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../corpus/rdkit_dump.jsonl.gz");
    let file = std::fs::File::open(&path)
        .unwrap_or_else(|e| panic!("cannot open {}: {e}", path.display()));
    let mut text = String::new();
    GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("gunzip");
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).expect("dump line"))
        .collect()
}

#[test]
fn graph_matches_rdkit() {
    let dump = load_dump();
    assert!(dump.len() > 7000, "dump looks truncated: {}", dump.len());

    let mut n_checked = 0usize;
    let mut failures: Vec<String> = Vec::new();
    let fail = |msg: String, failures: &mut Vec<String>| {
        if failures.len() < 30 {
            failures.push(msg);
        }
    };

    for rec in &dump {
        let smiles = rec["s"].as_str().expect("s");
        let ref_atoms = rec["a"].as_array().expect("a");
        let ref_bonds = rec["b"].as_array().expect("b");
        let ref_cip: std::collections::HashMap<usize, char> = rec["cip"]
            .as_array()
            .expect("cip")
            .iter()
            .map(|e| {
                let pair = e.as_array().expect("cip pair");
                (
                    pair[0].as_u64().expect("idx") as usize,
                    pair[1]
                        .as_str()
                        .expect("code")
                        .chars()
                        .next()
                        .expect("char"),
                )
            })
            .collect();
        let ref_rings: Vec<Vec<usize>> = rec["r"]
            .as_array()
            .expect("r")
            .iter()
            .map(|r| {
                r.as_array()
                    .expect("ring")
                    .iter()
                    .map(|v| v.as_u64().expect("atom idx") as usize)
                    .collect()
            })
            .collect();

        let g = match molrs::graph::build_molecule_graph(smiles) {
            Ok(g) => g,
            Err(e) => {
                fail(format!("{smiles}: build error: {e}"), &mut failures);
                continue;
            }
        };
        n_checked += 1;

        if g.atoms.len() != ref_atoms.len() {
            fail(
                format!(
                    "{smiles}: atom count {} != rdkit {}",
                    g.atoms.len(),
                    ref_atoms.len()
                ),
                &mut failures,
            );
            continue;
        }

        for (i, (ours, r)) in g.atoms.iter().zip(ref_atoms).enumerate() {
            let (sym, arom, chg, in_ring) = (
                r[0].as_str().expect("sym"),
                r[1].as_i64().expect("arom") == 1,
                r[2].as_i64().expect("chg") as i8,
                r[3].as_i64().expect("in_ring") == 1,
            );
            if ours.symbol != sym {
                fail(
                    format!("{smiles}: atom {i} symbol {} != {sym}", ours.symbol),
                    &mut failures,
                );
            }
            if ours.is_aromatic != arom {
                fail(
                    format!("{smiles}: atom {i} aromatic {} != {arom}", ours.is_aromatic),
                    &mut failures,
                );
            }
            if ours.formal_charge != chg {
                fail(
                    format!("{smiles}: atom {i} charge {} != {chg}", ours.formal_charge),
                    &mut failures,
                );
            }
            if ours.in_ring != in_ring {
                fail(
                    format!("{smiles}: atom {i} in_ring {} != {in_ring}", ours.in_ring),
                    &mut failures,
                );
            }
            if ours.chiral_tag != ref_cip.get(&i).copied() {
                fail(
                    format!(
                        "{smiles}: atom {i} CIP {:?} != rdkit {:?}",
                        ours.chiral_tag,
                        ref_cip.get(&i)
                    ),
                    &mut failures,
                );
            }
        }

        if g.ring_atom_sets != ref_rings {
            fail(
                format!(
                    "{smiles}: rings {:?} != rdkit {:?}",
                    g.ring_atom_sets, ref_rings
                ),
                &mut failures,
            );
        }

        if g.bonds.len() != ref_bonds.len() {
            fail(
                format!(
                    "{smiles}: bond count {} != rdkit {}",
                    g.bonds.len(),
                    ref_bonds.len()
                ),
                &mut failures,
            );
            continue;
        }
        for (k, (ours, r)) in g.bonds.iter().zip(ref_bonds).enumerate() {
            let (bi, ei, order) = (
                r[0].as_u64().expect("begin") as usize,
                r[1].as_u64().expect("end") as usize,
                r[2].as_f64().expect("order"),
            );
            if (ours.begin_idx, ours.end_idx) != (bi, ei) {
                fail(
                    format!(
                        "{smiles}: bond {k} ({},{}) != rdkit ({bi},{ei})",
                        ours.begin_idx, ours.end_idx
                    ),
                    &mut failures,
                );
            }
            if ours.bond_order != order {
                fail(
                    format!(
                        "{smiles}: bond {k} order {} != rdkit {order}",
                        ours.bond_order
                    ),
                    &mut failures,
                );
            }
            let ref_stereo = match r[3].as_str().expect("stereo") {
                "" => None,
                st => st.chars().next(),
            };
            if ours.stereo != ref_stereo {
                fail(
                    format!(
                        "{smiles}: bond {k} stereo {:?} != rdkit {ref_stereo:?}",
                        ours.stereo
                    ),
                    &mut failures,
                );
            }
        }
    }

    println!(
        "rdkit compat (S1.2): {n_checked}/{} molecules compared",
        dump.len()
    );
    assert!(
        failures.is_empty(),
        "{} mismatches (showing up to 30):\n{}",
        failures.len(),
        failures.join("\n")
    );
}
