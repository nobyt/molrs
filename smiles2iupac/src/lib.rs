//! smiles2iupac: SMILES から IUPAC 2013 優先名を生成するライブラリ (Rust 移植版)。
//!
//! Python 実装 (`src/smiles2iupac/`) の移植。移植計画は RUST_PORT_PLAN.md を参照。
//!
//! 現在の対応範囲 (S2.2 パイプライン骨格):
//! - 保留名テーブル (立体署名付き複合キーで照合)
//! - 直鎖アルカン (C1〜C30)
//!
//! 未対応の構造は [`NameError::Unsupported`] を返す (移植の進行に応じて拡大)。

pub mod assemble;
pub mod chain;
pub mod constants;
pub mod functional_group;
pub mod retained;
pub mod substituent;

use molrs::graph::{build_molecule_graph, MoleculeGraph};
pub use molrs::ChemError;

/// 命名エラー。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NameError {
    /// 不正な SMILES
    InvalidSmiles(String),
    /// この構造はまだ移植されていない (移植完了までの暫定エラー)
    Unsupported(String),
}

impl std::fmt::Display for NameError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NameError::InvalidSmiles(s) => write!(f, "Invalid SMILES: {s}"),
            NameError::Unsupported(s) => write!(f, "Unsupported structure: {s}"),
        }
    }
}

impl std::error::Error for NameError {}

/// SMILES 文字列を IUPAC 2013 優先名に変換する。
///
/// Python 版 `smiles_to_iupac()` に対応。
pub fn smiles_to_iupac(smiles: &str) -> Result<String, NameError> {
    let smiles = smiles.trim();
    if smiles.contains('.') {
        // 多成分 (塩・混合物) は S5.8 で対応
        return Err(NameError::Unsupported(format!("multicomponent: {smiles}")));
    }
    let graph =
        build_molecule_graph(smiles).map_err(|e| NameError::InvalidSmiles(e.to_string()))?;

    // 保留名テーブル (Python _try_retained_name 相当)
    if let Some(name) = retained::try_retained_name(&graph) {
        return Ok(name.to_string());
    }

    // 直鎖アルカン (パイプライン骨格の最小ケース)
    if let Some(name) = try_linear_alkane(&graph) {
        return Ok(name);
    }

    // 非環式の系統名 (S3.2-S3.5)
    if let Some(name) = try_acyclic_name(&graph) {
        return Ok(name);
    }

    Err(NameError::Unsupported(smiles.to_string()))
}

/// 非環式分子の系統名を試みる。環・未対応官能基・命名不能置換基は None。
fn try_acyclic_name(g: &MoleculeGraph) -> Option<String> {
    use functional_group::{detect_groups, principal_group};

    // 環を含む分子は未対応
    if g.atoms.iter().any(|a| a.in_ring) {
        return None;
    }
    // 電荷は未対応
    if g.atoms.iter().any(|a| a.formal_charge != 0) {
        return None;
    }
    // C 以外の重原子が主鎖外に許されるのは対応済み置換基/官能基のみ。
    // まず官能基検出。
    let groups = detect_groups(g);
    let principal = principal_group(&groups);
    let gtype = principal.map(|p| p.group_type).unwrap_or("alkane");

    // 対応する主基の suffix マッピング
    let suffix = acyclic_suffix(gtype)?;

    let is_amine = matches!(gtype, "amine" | "diamine");

    // アミンは N のみを主基とし、最長鎖を N 隣接炭素が低ロカントになるよう選ぶ。
    // (N 上の他の炭素は N-置換基になる)
    let (chain, principal_atoms, mut suffix_locants);
    if is_amine {
        let n_atoms: Vec<usize> = principal
            .map(|p| {
                p.atom_indices
                    .iter()
                    .copied()
                    .filter(|&a| g.atoms[a].symbol == "N")
                    .collect()
            })
            .unwrap_or_default();
        let mut c = chain::find_principal_chain(g, None);
        if c.length() == 0 {
            return None;
        }
        // N 隣接炭素の最小ロカントで向きを決める
        let n_adj_locant = |ch: &chain::PrincipalChain| -> usize {
            ch.atom_indices
                .iter()
                .enumerate()
                .filter(|(_, &cc)| {
                    n_atoms
                        .iter()
                        .any(|&n| g.adjacency[cc].contains(&n))
                })
                .map(|(i, _)| i + 1)
                .min()
                .unwrap_or(usize::MAX)
        };
        let loc_fwd = n_adj_locant(&c);
        let rev: Vec<usize> = c.atom_indices.iter().rev().copied().collect();
        let crev = chain::PrincipalChain { atom_indices: rev };
        let loc_rev = n_adj_locant(&crev);
        if loc_rev < loc_fwd {
            c = crev;
        }
        // 主基 N に隣接する全鎖炭素のロカント (diamine 用)
        let mut locs: Vec<usize> = c
            .atom_indices
            .iter()
            .enumerate()
            .filter(|(_, &cc)| n_atoms.iter().any(|&n| g.adjacency[cc].contains(&n)))
            .map(|(i, _)| i + 1)
            .collect();
        locs.sort_unstable();
        locs.dedup();
        if locs.is_empty() {
            return None; // 鎖が N に隣接しない (異常)
        }
        chain = c;
        principal_atoms = n_atoms;
        suffix_locants = locs;
    } else {
        let c = chain::find_principal_chain(g, principal);
        if c.length() == 0 {
            return None;
        }
        let mut sl: Vec<usize> = Vec::new();
        let mut pa: Vec<usize> = Vec::new();
        if let Some(p) = principal {
            pa = p.atom_indices.clone();
            for &ai in &p.atom_indices {
                if g.atoms[ai].symbol == "C" {
                    if let Some(loc) = c.locant_of(ai) {
                        sl.push(loc);
                    }
                }
            }
            sl.sort_unstable();
            sl.dedup();
        }
        chain = c;
        principal_atoms = pa;
        suffix_locants = sl;
    }
    suffix_locants.sort_unstable();
    suffix_locants.dedup();

    let (ene, yne) = chain::multiple_bond_locants(g, &chain);

    // 置換基収集 (命名不能なら None)
    let (subs, n_subs) =
        substituent::collect_substituents(g, &chain.atom_indices, &principal_atoms)?;

    // 立体記述子 (E/Z 二重結合、R/S 不斉中心) を主鎖上のもののみ収集
    let stereo = collect_stereo(g, &chain);

    assemble::assemble_name(
        chain.length(),
        gtype,
        &ene,
        &yne,
        &subs,
        &n_subs,
        suffix,
        &suffix_locants,
        &stereo,
    )
}

/// 主鎖上の立体記述子を収集する。(ロカント, "E"/"Z"/"R"/"S")。
fn collect_stereo(g: &MoleculeGraph, chain: &chain::PrincipalChain) -> Vec<(usize, String)> {
    use molrs::graph::get_bond_order;
    let mut out: Vec<(usize, String)> = Vec::new();
    let path = &chain.atom_indices;
    // E/Z: 主鎖二重結合の bond.stereo
    for i in 0..path.len().saturating_sub(1) {
        let (a, b) = (path[i], path[i + 1]);
        if get_bond_order(g, a, b) != 2.0 {
            continue;
        }
        if let Some(bond) = g
            .bonds
            .iter()
            .find(|bd| (bd.begin_idx == a && bd.end_idx == b) || (bd.begin_idx == b && bd.end_idx == a))
        {
            if let Some(ez) = bond.stereo {
                out.push((i + 1, ez.to_string()));
            }
        }
    }
    // R/S: 主鎖炭素の chiral_tag
    for (i, &c) in path.iter().enumerate() {
        if let Some(rs) = g.atoms[c].chiral_tag {
            out.push((i + 1, rs.to_string()));
        }
    }
    out
}

/// 非環式で対応する主基 → assemble 用 suffix。未対応は None。
fn acyclic_suffix(gtype: &str) -> Option<&'static str> {
    Some(match gtype {
        "alkane" => "ane",
        "alkene" => "ene",
        "alkyne" => "yne",
        "alcohol" => "ol",
        "diol" => "diol",
        "triol" => "triol",
        "ketone" => "one",
        "dione" => "dione",
        "trione" => "trione",
        "aldehyde" => "al",
        "dial" => "dial",
        "carboxylic_acid" => "oic acid",
        "dioic_acid" => "dioic acid",
        "amine" => "amine",
        "diamine" => "diamine",
        "nitrile" => "nitrile",
        "dinitrile" => "dinitrile",
        _ => return None,
    })
}

/// 直鎖アルカン (メタン〜トリアコンタン) の命名。
/// 全原子が中性の炭素、環なし、分岐なし、単結合のみの場合のみ Some。
fn try_linear_alkane(g: &MoleculeGraph) -> Option<String> {
    let n_kept = g.parser_to_graph.iter().flatten().count();
    if n_kept == 0 || n_kept > 30 {
        return None;
    }
    let heavy_degree = |i: usize| g.adjacency[i].iter().filter(|&&x| x < n_kept).count();
    for i in 0..n_kept {
        let a = &g.atoms[i];
        if a.symbol != "C" || a.formal_charge != 0 || a.in_ring || a.is_aromatic {
            return None;
        }
        if heavy_degree(i) > 2 {
            return None;
        }
    }
    // 全結合が単結合
    for b in &g.bonds {
        if b.begin_idx < n_kept && b.end_idx < n_kept && b.bond_order != 1.0 {
            return None;
        }
    }
    Some(format!("{}ane", constants::CHAIN_PREFIX[n_kept]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_alkanes() {
        assert_eq!(smiles_to_iupac("C").unwrap(), "methane");
        assert_eq!(smiles_to_iupac("CC").unwrap(), "ethane");
        assert_eq!(smiles_to_iupac("CCC").unwrap(), "propane");
        assert_eq!(smiles_to_iupac("C(C)C").unwrap(), "propane"); // 表記違い
        assert_eq!(
            smiles_to_iupac("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC").unwrap(),
            "triacontane"
        );
    }

    #[test]
    fn acyclic_systematic() {
        // 非環式の系統名 (S3.2-S3.5)
        assert_eq!(smiles_to_iupac("CC(C)C").unwrap(), "2-methylpropane");
        assert_eq!(smiles_to_iupac("CCCCO").unwrap(), "butan-1-ol");
        assert_eq!(smiles_to_iupac("CC(=O)C").unwrap(), "propan-2-one");
        assert_eq!(smiles_to_iupac("OCCO").unwrap(), "ethane-1,2-diol");
        assert_eq!(smiles_to_iupac("C/C=C/C").unwrap(), "(2E)-but-2-ene");
        // 環はまだ未対応
        assert!(matches!(
            smiles_to_iupac("C1CC1"),
            Err(NameError::Unsupported(_))
        ));
    }

    #[test]
    fn invalid_smiles() {
        assert!(matches!(
            smiles_to_iupac("hello"),
            Err(NameError::InvalidSmiles(_))
        ));
    }

    #[test]
    fn constants_table_counts() {
        // 生成テーブルの件数が Python 側と一致すること
        assert_eq!(constants::FUNCTIONAL_GROUPS.len(), 165);
        assert_eq!(constants::CHAIN_PREFIX.len(), 31);
        assert_eq!(constants::RETAINED_NAMES_RAW.len(), 205);
        assert_eq!(constants::CHAIN_PREFIX[1], "meth");
        assert_eq!(constants::CHAIN_PREFIX[30], "triacont");
        assert_eq!(constants::MULTIPLIER[2], "di");
        assert_eq!(constants::halogen_prefix("Cl"), Some("chloro"));
        assert_eq!(constants::suffix_of("carboxylic_acid"), Some("oic acid"));
    }
}
