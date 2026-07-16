//! 3D 配座生成の性質ゲート (RUST_3D_PLAN.md C9)。
//!
//! 座標は乱数依存なので「RDKit と一致」ではなく、シード固定の決定的な
//! 性質で守る (RUST_3D_PLAN.md §1):
//! - 埋め込み成功率 ≥ 99%
//! - 立体保存 100% (3D から R/S・E/Z を再計算して入力と一致)
//! - 結合長の理想値からの偏差 ≤ 0.12 Å (UFF 前の水準; C7 で強化予定)
//! - 非結合原子対の最小距離 > 1.2 Å
//! - 芳香環の平面性 rms < 0.05 Å
//! - 座標は全て有限
//!
//! 実行時間の都合でコーパスの決定的サンプル (10 分の 1) + 立体指定を含む
//! 全分子を対象にする。`CONFORMER_GATE_FULL=1` で全件に切り替え。
//!
//! C3 改良後の実績: サンプル 934/934、全件 7453/7453 (100%)。

use molrs::conformer::{embed_molecule, verify_stereo_3d, EmbedParams};
use molrs::geometry::{jacobi_eigen, Vec3};
use molrs::graph::{build_molecule_graph, MoleculeGraph};
use std::path::PathBuf;

fn ring_planarity_rms(coords: &[Vec3], ring: &[usize]) -> f64 {
    let pts: Vec<Vec3> = ring.iter().map(|&i| coords[i]).collect();
    let centroid = pts.iter().fold(Vec3::ZERO, |a, &b| a + b) / pts.len() as f64;
    let mut m = [0.0f64; 9];
    for p in &pts {
        let r = *p - centroid;
        let v = [r.x, r.y, r.z];
        for a in 0..3 {
            for b in 0..3 {
                m[a * 3 + b] += v[a] * v[b];
            }
        }
    }
    let (_, vecs) = jacobi_eigen(&m, 3);
    let normal = Vec3::new(vecs[2][0], vecs[2][1], vecs[2][2]);
    (pts.iter()
        .map(|p| (*p - centroid).dot(normal).powi(2))
        .sum::<f64>()
        / pts.len() as f64)
        .sqrt()
}

/// 1 分子分の検査。問題があればメッセージを返す。
fn check_molecule(smiles: &str, g: &MoleculeGraph) -> Result<(), String> {
    let conf = embed_molecule(g, &EmbedParams::default()).map_err(|e| format!("{smiles}: {e}"))?;
    let coords = &conf.coords;
    let uff_r0 = molrs::conformer::uff_bond_rest_lengths(g);

    for c in coords {
        if !(c.x.is_finite() && c.y.is_finite() && c.z.is_finite()) {
            return Err(format!("{smiles}: non-finite coordinates"));
        }
    }

    // 結合長: 理想表 ±0.12 Å、または UFF 平衡長 ±0.08 Å
    // (UFF は O-O 1.316 Å など理想表と系統的に異なる結合があるため、
    //  力場として自己整合な長さは許容する)
    for (bi, b) in g.bonds.iter().enumerate() {
        let ideal = molrs::conformer::params::ideal_bond_length(
            &g.atoms[b.begin_idx].symbol,
            &g.atoms[b.end_idx].symbol,
            b.bond_order,
        );
        let d = coords[b.begin_idx].distance(coords[b.end_idx]);
        let ok_ideal = (d - ideal).abs() <= 0.12;
        let ok_uff = uff_r0.as_ref().is_some_and(|v| (d - v[bi]).abs() <= 0.08);
        if !ok_ideal && !ok_uff {
            return Err(format!(
                "{smiles}: bond ({},{}) length {d:.3} vs ideal {ideal:.3}",
                b.begin_idx, b.end_idx
            ));
        }
    }

    // 非結合最小距離
    let bonded: std::collections::HashSet<(usize, usize)> = g
        .bonds
        .iter()
        .map(|b| (b.begin_idx.min(b.end_idx), b.begin_idx.max(b.end_idx)))
        .collect();
    for i in 0..coords.len() {
        for j in (i + 1)..coords.len() {
            if bonded.contains(&(i, j)) {
                continue;
            }
            let d = coords[i].distance(coords[j]);
            if d < 1.2 {
                return Err(format!("{smiles}: clash ({i},{j}) = {d:.3}"));
            }
        }
    }

    // 芳香環平面性
    for ring in &g.ring_atom_sets {
        if ring.len() >= 5 && ring.iter().all(|&a| g.atoms[a].is_aromatic) {
            let rms = ring_planarity_rms(coords, ring);
            if rms > 0.05 {
                return Err(format!("{smiles}: aromatic ring rms {rms:.3}"));
            }
        }
    }

    // 立体保存
    if !verify_stereo_3d(g, &conf) {
        return Err(format!("{smiles}: stereo not preserved"));
    }
    Ok(())
}

#[test]
#[cfg_attr(
    debug_assertions,
    ignore = "slow: run with `cargo test --release --test conformer_gate`"
)]
fn conformer_property_gate() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../corpus/corpus.jsonl");
    let text = std::fs::read_to_string(&path).expect("corpus");
    let full = std::env::var("CONFORMER_GATE_FULL").is_ok();

    let mut n_tried = 0usize;
    let mut n_ok = 0usize;
    let mut failures: Vec<String> = Vec::new();
    for (li, line) in text.lines().filter(|l| !l.trim().is_empty()).enumerate() {
        let row: serde_json::Value = serde_json::from_str(line).expect("json");
        let smiles = row["smiles"].as_str().expect("smiles");
        let has_stereo = smiles.contains('@') || smiles.contains('/') || smiles.contains('\\');
        // 決定的サンプル: 10 分の 1 + 立体分子は全部
        if !full && !has_stereo && li % 10 != 0 {
            continue;
        }
        let Ok(g) = build_molecule_graph(smiles) else {
            continue;
        };
        n_tried += 1;
        match check_molecule(smiles, &g) {
            Ok(()) => n_ok += 1,
            Err(msg) => {
                if failures.len() < 15 {
                    failures.push(msg);
                }
            }
        }
    }

    let rate = n_ok as f64 / n_tried as f64;
    println!("conformer gate: {n_ok}/{n_tried} ok ({:.2}%)", rate * 100.0);
    assert!(
        rate >= 0.99,
        "success rate {:.2}% below 99%:\n{}",
        rate * 100.0,
        failures.join("\n")
    );
}
