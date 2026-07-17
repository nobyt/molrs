//! 2D V2000 MOL ブロック書き出し (D11)。
//!
//! くさびコード付き (wedge up = 1, down = 6; 方向は第 1 原子 = 細端から)。
//! 隠し H は出力せず可視原子のみを再番号付けする (読み手が暗黙 H を補う)。
//! ケクレ次数で書く (芳香族 4 は使わない — 読み手依存を減らす)。
//! RDKit での再認識 (CIPLabeler round-trip, D13) の入出力に使う。

use crate::graph::MoleculeGraph;

use super::{Coords2D, WedgeDir};

/// 結合が環内か。
fn is_ring_bond(g: &MoleculeGraph, i: usize, j: usize) -> bool {
    g.ring_atom_sets.iter().any(|ring| {
        let n = ring.len();
        (0..n).any(|k| {
            let (a, b) = (ring[k], ring[(k + 1) % n]);
            (a == i && b == j) || (a == j && b == i)
        })
    })
}

/// 立体を持ちうるのに未指定の非環二重結合か (両端に他の重原子隣接がある)。
/// 2D 描画は見かけ上の幾何を持ってしまうため、MOL では crossed bond
/// (stereo code 3 = cis/trans either) で「未指定」を明示する。
fn is_unspecified_stereogenic_double(g: &MoleculeGraph, c: &Coords2D, bi: usize) -> bool {
    let b = &g.bonds[bi];
    if g.kekule_bond_orders[bi] != 2.0 || b.stereo.is_some() {
        return false;
    }
    let (i, j) = (b.begin_idx, b.end_idx);
    if is_ring_bond(g, i, j) {
        return false;
    }
    let _ = c;
    let has_other_heavy = |x: usize, other: usize| {
        g.adjacency[x]
            .iter()
            .any(|&nb| nb != other && g.atoms[nb].symbol != "H")
    };
    has_other_heavy(i, j) && has_other_heavy(j, i)
}

fn bond_type_code(order: f64) -> u8 {
    if order == 2.0 {
        2
    } else if order == 3.0 {
        3
    } else {
        1
    }
}

/// 2D MOL ブロックを生成する。座標はレイアウト単位 (結合長 = 1.0)。
pub fn to_mol_block_2d(g: &MoleculeGraph, c: &Coords2D, title: &str) -> String {
    // 可視原子の再番号付け (1 始まり)
    let mut new_idx = vec![0usize; g.atoms.len()];
    let mut visible: Vec<usize> = Vec::new();
    for (i, ni) in new_idx.iter_mut().enumerate() {
        if !c.hidden[i] {
            *ni = visible.len() + 1;
            visible.push(i);
        }
    }
    let bonds: Vec<usize> = (0..g.bonds.len())
        .filter(|&bi| {
            let b = &g.bonds[bi];
            !c.hidden[b.begin_idx] && !c.hidden[b.end_idx]
        })
        .collect();

    let mut s = String::new();
    s.push_str(title);
    s.push('\n');
    s.push_str("  molrs 2D\n");
    s.push('\n');
    s.push_str(&format!(
        "{:3}{:3}  0  0  0  0  0  0  0  0999 V2000\n",
        visible.len(),
        bonds.len()
    ));
    for &i in &visible {
        s.push_str(&format!(
            "{:10.4}{:10.4}{:10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            c.pos[i].x, c.pos[i].y, 0.0, g.atoms[i].symbol
        ));
    }
    for &bi in &bonds {
        let b = &g.bonds[bi];
        let (mut a1, mut a2) = (b.begin_idx, b.end_idx);
        let stereo_code = match &c.wedge[bi] {
            Some(w) => {
                // 細端 (立体中心) を第 1 原子に
                if w.narrow != a1 {
                    std::mem::swap(&mut a1, &mut a2);
                }
                match w.dir {
                    WedgeDir::Up => 1,
                    WedgeDir::Down => 6,
                }
            }
            None => {
                if is_unspecified_stereogenic_double(g, c, bi) {
                    3 // crossed bond: cis/trans either
                } else {
                    0
                }
            }
        };
        s.push_str(&format!(
            "{:3}{:3}{:3}{:3}\n",
            new_idx[a1],
            new_idx[a2],
            bond_type_code(g.kekule_bond_orders[bi]),
            stereo_code
        ));
    }
    // 形式電荷
    let charged: Vec<(usize, i8)> = visible
        .iter()
        .filter(|&&i| g.atoms[i].formal_charge != 0)
        .map(|&i| (new_idx[i], g.atoms[i].formal_charge))
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
    use crate::depict::{compute_coords_2d, LayoutParams};
    use crate::graph::build_molecule_graph;

    fn mol_of(smiles: &str) -> String {
        let g = build_molecule_graph(smiles).unwrap();
        let c = compute_coords_2d(&g, &LayoutParams::default()).unwrap();
        to_mol_block_2d(&g, &c, smiles)
    }

    #[test]
    fn header_and_counts() {
        let m = mol_of("CCO");
        let lines: Vec<&str> = m.lines().collect();
        assert_eq!(lines[1], "  molrs 2D");
        assert!(lines[3].contains("  3  2  0"));
        assert!(m.ends_with("M  END\n"));
    }

    #[test]
    fn wedge_codes_present_for_stereocenter() {
        let m = mol_of("C[C@H](O)CC");
        let has_wedge = m
            .lines()
            .skip(4)
            .any(|l| l.len() >= 12 && (l.ends_with("  1") || l.ends_with("  6")));
        assert!(has_wedge, "{m}");
    }

    #[test]
    fn charges_written() {
        let m = mol_of("C[N+](C)(C)C");
        assert!(m.contains("M  CHG  1"), "{m}");
    }

    #[test]
    fn aromatic_written_kekule() {
        let m = mol_of("c1ccccc1");
        // 結合タイプ 4 (芳香族) を使わない
        for l in m.lines().skip(10).take(6) {
            assert!(!l.trim_end().ends_with(" 4"), "{l}");
        }
        // 二重結合 3 本
        let n_double = m
            .lines()
            .skip(10)
            .filter(|l| l.split_whitespace().nth(2) == Some("2"))
            .count();
        assert_eq!(n_double, 3);
    }
}
