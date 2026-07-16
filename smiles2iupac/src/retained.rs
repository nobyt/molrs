//! 保留名テーブルのルックアップ (S2.2)。
//!
//! Python 版は RDKit 正規 SMILES (立体込み) をキーにする。Rust 版は
//! RDKit と独立に自己完結させるため、**複合キー** を使う:
//!
//! ```text
//! <立体なし正規SMILES> ["|" <出力位置><CIPコード>,...] ["|" <位置>-<位置><E/Z>,...]
//! ```
//!
//! CIP コード (R/S) と E/Z は入力表記に依らない不変量なので、
//! エナンチオマー対 (D/L-アラニン等、テーブルに 52 の立体キーがある) を
//! 正しく区別でき、立体指定のない入力が立体キーに誤ヒットすることもない。
//! テーブル側は生 SMILES (constants::RETAINED_NAMES_RAW) を初期化時に
//! 同じ関数で再キー化する。

use std::collections::HashMap;
use std::sync::LazyLock;

use molrs::canon::to_canonical_smiles_with_order;
use molrs::graph::MoleculeGraph;

use crate::constants::RETAINED_NAMES_RAW;

/// 分子の複合正規化キー。
pub fn canonical_key(g: &MoleculeGraph) -> String {
    let (smi, order) = to_canonical_smiles_with_order(g);
    // グラフ原子 idx → 出力位置
    let mut pos = HashMap::new();
    for (p, &gi) in order.iter().enumerate() {
        pos.insert(gi, p);
    }

    let mut atom_sig: Vec<(usize, char)> = g
        .atoms
        .iter()
        .filter_map(|a| a.chiral_tag.and_then(|c| pos.get(&a.idx).map(|&p| (p, c))))
        .collect();
    atom_sig.sort_unstable();

    let mut bond_sig: Vec<(usize, usize, char)> = g
        .bonds
        .iter()
        .filter_map(|b| {
            b.stereo.and_then(|c| {
                let (Some(&pa), Some(&pb)) = (pos.get(&b.begin_idx), pos.get(&b.end_idx)) else {
                    return None;
                };
                Some((pa.min(pb), pa.max(pb), c))
            })
        })
        .collect();
    bond_sig.sort_unstable();

    let mut key = smi;
    if !atom_sig.is_empty() {
        key.push('|');
        for (i, (p, c)) in atom_sig.iter().enumerate() {
            if i > 0 {
                key.push(',');
            }
            key.push_str(&format!("{p}{c}"));
        }
    }
    if !bond_sig.is_empty() {
        key.push('|');
        for (i, (pa, pb, c)) in bond_sig.iter().enumerate() {
            if i > 0 {
                key.push(',');
            }
            key.push_str(&format!("{pa}-{pb}{c}"));
        }
    }
    key
}

/// 再キー化済みの保留名テーブル。
static RETAINED_NAMES: LazyLock<HashMap<String, &'static str>> = LazyLock::new(|| {
    let mut map = HashMap::with_capacity(RETAINED_NAMES_RAW.len());
    for &(smiles, name) in RETAINED_NAMES_RAW {
        match molrs::graph::build_molecule_graph(smiles) {
            Ok(g) => {
                let key = canonical_key(&g);
                let prev = map.insert(key, name);
                debug_assert!(
                    prev.is_none(),
                    "retained-name key collision for {smiles} ({name} vs {})",
                    prev.unwrap()
                );
            }
            Err(e) => {
                debug_assert!(false, "retained-name SMILES fails to parse: {smiles}: {e}");
            }
        }
    }
    map
});

/// 保留名の照合 (Python `_try_retained_name` 相当)。
pub fn try_retained_name(g: &MoleculeGraph) -> Option<&'static str> {
    RETAINED_NAMES.get(&canonical_key(g)).copied()
}

#[cfg(test)]
mod tests {
    use super::*;
    use molrs::graph::build_molecule_graph;

    fn lookup(smiles: &str) -> Option<&'static str> {
        try_retained_name(&build_molecule_graph(smiles).expect("valid"))
    }

    #[test]
    fn table_builds_without_collisions() {
        // LazyLock 初期化を強制 (debug_assert が衝突を検出する)
        assert!(RETAINED_NAMES.len() > 190);
    }

    #[test]
    fn enantiomers_are_distinguished() {
        assert_eq!(
            lookup("C[C@H](N)C(=O)O"),
            Some("(2S)-2-aminopropanoic acid")
        );
        assert_eq!(
            lookup("C[C@@H](N)C(=O)O"),
            Some("(2R)-2-aminopropanoic acid")
        );
        // 立体指定なしは立体キーにヒットしない
        assert_ne!(lookup("CC(N)C(=O)O"), Some("(2S)-2-aminopropanoic acid"));
        assert_ne!(lookup("CC(N)C(=O)O"), Some("(2R)-2-aminopropanoic acid"));
    }

    #[test]
    fn notation_invariance() {
        // 同じ分子の別表記は同じ保留名に届く
        let a = lookup("N[C@@H](C)C(=O)O"); // L-アラニンの別表記
        assert_eq!(a, Some("(2S)-2-aminopropanoic acid"));
    }
}
