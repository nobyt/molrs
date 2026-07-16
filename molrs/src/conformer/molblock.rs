//! V2000 MOL ブロック書き出し (RUST_3D_PLAN.md C8)。
//!
//! 生成した 3D 座標の可視化・RDKit 相互チェック (tools/compare_conformers.py)
//! に使う最小実装。プロパティは電荷 (`M  CHG`) のみ対応。

use crate::conformer::Conformer;
use crate::graph::MoleculeGraph;

fn bond_type_code(order: f64) -> u8 {
    if order == 2.0 {
        2
    } else if order == 3.0 {
        3
    } else if order == 1.5 {
        4 // aromatic
    } else {
        1
    }
}

/// V2000 MOL ブロックを生成する。
pub fn to_mol_block(g: &MoleculeGraph, conf: &Conformer, title: &str) -> String {
    let n_atoms = g.atoms.len();
    let n_bonds = g.bonds.len();
    assert_eq!(conf.coords.len(), n_atoms, "coordinate count mismatch");

    let mut s = String::new();
    s.push_str(title);
    s.push('\n');
    s.push_str("  molrs 3D\n");
    s.push('\n');
    s.push_str(&format!(
        "{n_atoms:3}{n_bonds:3}  0  0  0  0  0  0  0  0999 V2000\n"
    ));
    for (a, c) in g.atoms.iter().zip(&conf.coords) {
        s.push_str(&format!(
            "{:10.4}{:10.4}{:10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            c.x, c.y, c.z, a.symbol
        ));
    }
    for b in &g.bonds {
        s.push_str(&format!(
            "{:3}{:3}{:3}  0\n",
            b.begin_idx + 1,
            b.end_idx + 1,
            bond_type_code(b.bond_order)
        ));
    }
    // 形式電荷
    let charged: Vec<(usize, i8)> = g
        .atoms
        .iter()
        .filter(|a| a.formal_charge != 0)
        .map(|a| (a.idx + 1, a.formal_charge))
        .collect();
    for chunk in charged.chunks(8) {
        s.push_str(&format!("M  CHG{:3}", chunk.len()));
        for &(idx, chg) in chunk {
            s.push_str(&format!("{idx:4}{chg:4}"));
        }
        s.push('\n');
    }
    s.push_str("M  END\n");
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformer::{embed_molecule, EmbedParams};
    use crate::graph::build_molecule_graph;

    #[test]
    fn mol_block_structure() {
        let g = build_molecule_graph("CC(=O)[O-]").unwrap();
        let conf = embed_molecule(&g, &EmbedParams::default()).unwrap();
        let block = to_mol_block(&g, &conf, "acetate");
        let lines: Vec<&str> = block.lines().collect();
        assert_eq!(lines[0], "acetate");
        assert!(lines[3].contains("V2000"));
        // counts 行: 原子 7 (C,C,O,O + 3H) 結合 6
        assert!(lines[3].starts_with("  7  6"));
        // 電荷行がある
        assert!(block.contains("M  CHG  1"));
        assert!(block.trim_end().ends_with("M  END"));
        // 原子行の元素記号
        assert!(lines[4].contains(" C "));
        assert!(lines[6].contains(" O "));
    }

    #[test]
    fn benzene_block_has_aromatic_bonds() {
        let g = build_molecule_graph("c1ccccc1").unwrap();
        let conf = embed_molecule(&g, &EmbedParams::default()).unwrap();
        let block = to_mol_block(&g, &conf, "benzene");
        // 芳香族結合コード 4 が 6 本
        let n_arom = block
            .lines()
            .skip(4 + 12)
            .filter(|l| l.len() >= 12 && l[6..9].trim() == "4")
            .count();
        assert_eq!(n_arom, 6);
    }
}
