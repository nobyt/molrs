//! 金属結合の切断 (標準 InChI の disconnected-metal 表現、I20)。
//!
//! 標準 InChI は有機金属化合物を「金属を切り離した形」で表現する
//! (再結合形 `/r` 層は非標準オプション `RecMet` 専用なので標準には現れない)。
//! 例: `C[Hg]C` → `InChI=1S/2CH3.Hg`、`[BiH3]` → `InChI=1S/Bi.3H`。
//!
//! 切断時の電荷配分は結合相手の電気陰性度で決まる (コーパス実測):
//! - 金属-C / 金属-H → **等方開裂**。双方中性のまま (メチルラジカル等)。
//!   `C[Hg]C` に `/q` が付かないのはこのため。
//! - 金属-ヘテロ原子 (ハロゲン・O・N・S 等) → **異方開裂**。配位子が
//!   陰イオン、金属が同数の陽イオンになる。`C[Hg]Cl` が
//!   `/q;;+1/p-1` (Cl⁻ は後段の [`super::normalize`] で HCl + `/p-1` に
//!   中性化され、Hg⁺ が `/q` に残る) となるのはこのため。
//!
//! 金属-H の切断で生じた H は、どの重原子にも結合しない**独立した H 成分**に
//! なる (`[BiH3]` の式 `Bi.3H`)。この孤立 H 成分は c/h 層には何も寄与しない。

use crate::graph::{BondInfo, MoleculeGraph};
use std::collections::HashMap;

/// InChI が「金属」とみなす元素 (公式 InChI ソース `util.c` の金属表)。
///
/// 半金属のうち Sb・Sn・Bi・Po は金属側、B・Si・Ge・As・Te は非金属側に
/// 分類される点に注意 (コーパスでも `CC[AsH2]` や `Br[SiH2]C` は連結の
/// まま、`C[SbH2]` は切断されることを確認済み)。
pub(crate) fn is_metal(sym: &str) -> bool {
    matches!(
        sym,
        "Li" | "Be"
            | "Na"
            | "Mg"
            | "Al"
            | "K"
            | "Ca"
            | "Sc"
            | "Ti"
            | "V"
            | "Cr"
            | "Mn"
            | "Fe"
            | "Co"
            | "Ni"
            | "Cu"
            | "Zn"
            | "Ga"
            | "Rb"
            | "Sr"
            | "Y"
            | "Zr"
            | "Nb"
            | "Mo"
            | "Tc"
            | "Ru"
            | "Rh"
            | "Pd"
            | "Ag"
            | "Cd"
            | "In"
            | "Sn"
            | "Sb"
            | "Cs"
            | "Ba"
            | "La"
            | "Ce"
            | "Pr"
            | "Nd"
            | "Pm"
            | "Sm"
            | "Eu"
            | "Gd"
            | "Tb"
            | "Dy"
            | "Ho"
            | "Er"
            | "Tm"
            | "Yb"
            | "Lu"
            | "Hf"
            | "Ta"
            | "W"
            | "Re"
            | "Os"
            | "Ir"
            | "Pt"
            | "Au"
            | "Hg"
            | "Tl"
            | "Pb"
            | "Bi"
            | "Po"
            | "Fr"
            | "Ra"
            | "Ac"
            | "Th"
            | "Pa"
            | "U"
            | "Np"
            | "Pu"
            | "Am"
            | "Cm"
            | "Bk"
            | "Cf"
            | "Es"
            | "Fm"
            | "Md"
            | "No"
            | "Lr"
    )
}

/// 金属結合の切断結果。`metal_locked_ligands` は異方開裂で電荷を得た配位子
/// のうちハロゲン**以外** (O/N/S 等) の原子集合 ([`normalize`] 参照)。
pub(crate) struct Disconnected {
    pub(crate) graph: MoleculeGraph,
    pub(crate) metal_locked_ligands: std::collections::HashSet<usize>,
}

/// 金属原子に接続する全ての結合を切断したグラフを返す。
/// 金属結合がなければ元のグラフをそのまま複製して返す。
pub(crate) fn disconnect_metals(g: &MoleculeGraph) -> Disconnected {
    let is_m = |i: usize| is_metal(g.atoms[i].symbol.as_str());
    if !g.bonds.iter().any(|b| is_m(b.begin_idx) || is_m(b.end_idx)) {
        return Disconnected {
            graph: g.clone(),
            metal_locked_ligands: std::collections::HashSet::new(),
        };
    }

    let mut charges: Vec<i8> = g.atoms.iter().map(|a| a.formal_charge).collect();
    let mut bonds: Vec<BondInfo> = Vec::with_capacity(g.bonds.len());
    let mut kekule: Vec<f64> = Vec::with_capacity(g.bonds.len());
    let mut metal_locked_ligands = std::collections::HashSet::new();

    for (bi, b) in g.bonds.iter().enumerate() {
        let (i, j) = (b.begin_idx, b.end_idx);
        let (mi, mj) = (is_m(i), is_m(j));
        if !mi && !mj {
            bonds.push(b.clone());
            kekule.push(g.kekule_bond_orders[bi]);
            continue;
        }
        // 金属-金属結合は双方中性のまま切るだけ (電荷配分の非対称性がない)。
        if mi && mj {
            continue;
        }
        let (metal, ligand) = if mi { (i, j) } else { (j, i) };
        let lsym = g.atoms[ligand].symbol.as_str();
        let order = b.bond_order.round().max(1.0) as i8;
        // 多重結合 (order ≥ 2) の金属=配位子は**均等開裂**で中性に切る。
        // 実 InChI: `S=[Mo]` → `Mo.S` (電荷なし)、`O=[Mo]=O` → `Mo.2O`、
        // クロム酸 `O=[Cr](=O)([O-])[O-]` の `=O` も中性 O になる。電荷を
        // 移すのは order 1 の配位結合だけ。単結合カルコゲン化物 (`CO[Mo]`
        // → `/q-1;+1`) や `COCCO[Hg]` (I41) は従来どおり異方開裂する。
        if lsym != "C" && lsym != "H" && order == 1 {
            // 異方開裂: 結合次数分の電荷を配位子(−)と金属(+)に振る
            charges[ligand] -= order;
            charges[metal] += order;
            // ハロゲン化物イオンは通常どおりプロトン化 (`C[Hg]Cl` → HCl/p-1
            // が既に検証済み)。O/N/S 等は金属由来の電荷を**恒久**として残す
            // (`COCCO[Hg]` の実 InChI は `/q-1;+1` で O をプロトン化しない、
            // I41)。ただし配位子が「金属以外の重原子を持たない中性の
            // -OH/-SH」("aqua"/"sulfido-H" 型配位子、`inchi-1` で
            // `O[W](=O)(=O)O` → `2H2O.2O.W/q;;;;+2/p-2` などを確認済み)
            // の場合は例外: 金属に結合する前の原子価がちょうど H だけで
            // 満たされていた (=他の重原子置換基がゼロ、金属結合前は電荷
            // なし) ので、切断後は H2O/H2S に中性化させる (`is_protonatable`
            // に任せる) — アルコキシド `CO[Mo]` (炭素置換基を持つ) や
            // 既に電荷を持っていた配位子 `[O-][Fe]` (`inchi-1` で確認:
            // こちらは中性化されない) とはこの重原子置換基の有無と
            // 事前電荷の有無で区別する。
            let is_aqua_ligand = matches!(lsym, "O" | "S")
                && g.atoms[ligand].formal_charge == 0
                && g.adjacency[ligand]
                    .iter()
                    .filter(|&&nb| nb != metal)
                    .all(|&nb| g.atoms[nb].symbol == "H")
                && g.adjacency[ligand].iter().any(|&nb| nb != metal);
            if !matches!(lsym, "F" | "Cl" | "Br" | "I") && !is_aqua_ligand {
                metal_locked_ligands.insert(ligand);
            }
        }
    }

    let mut atoms = g.atoms.clone();
    for (a, &c) in atoms.iter_mut().zip(charges.iter()) {
        a.formal_charge = c;
    }

    let mut adjacency = vec![Vec::new(); atoms.len()];
    let mut bond_orders = HashMap::new();
    for b in &bonds {
        adjacency[b.begin_idx].push(b.end_idx);
        adjacency[b.end_idx].push(b.begin_idx);
        bond_orders.insert(
            (b.begin_idx.min(b.end_idx), b.begin_idx.max(b.end_idx)),
            b.bond_order,
        );
    }

    Disconnected {
        graph: MoleculeGraph {
            atoms,
            bonds,
            adjacency,
            bond_orders,
            ring_atom_sets: g.ring_atom_sets.clone(),
            kekule_bond_orders: kekule,
            parsed: g.parsed.clone(),
            parser_to_graph: g.parser_to_graph.clone(),
        },
        metal_locked_ligands,
    }
}

/// 重原子を一切含まない H だけの連結成分のサイズ一覧 (H 原子数)。
///
/// 金属水素化物の切断で生じた H は互いに結合していないのでサイズ 1 の成分が
/// 並ぶ (`[BiH3]` → `[1,1,1]` → 式 `3H`)。一方、水素分子 `[H][H]` は 2 個の H
/// が結合したサイズ 2 の成分 1 つで、実 InChI はこれを「骨格原子 1 個 +
/// その結合水素 1 個」として `InChI=1S/H2/h1H` と表現する。
/// H だけの各連結成分を (サイズ, 正味電荷) で走査する共通ヘルパ。
/// サイズ 1 かつ電荷が正の成分 (`[H+]` のように、他に結合先の重原子が
/// 一切ない裸のプロトン) は実 InChI では物質として数えない —
/// `[H+]` 単独が `InChI=1S/p+1` (式なし) になり、`[H+].[H-]` でも
/// H- 側だけが `H` の式として残り余った `[H+]` は `/p+1` に回るのが
/// その根拠 (I64)。これらは呼び出し側で `/p` に繰り込む。
fn hydrogen_components(g: &MoleculeGraph) -> Vec<(usize, i32)> {
    let n = g.atoms.len();
    let is_lone_h = |i: usize| {
        g.atoms[i].symbol == "H" && !g.adjacency[i].iter().any(|&nb| g.atoms[nb].symbol != "H")
    };
    let mut seen = vec![false; n];
    let mut out = Vec::new();
    for start in 0..n {
        if !is_lone_h(start) || seen[start] {
            continue;
        }
        let mut size = 0usize;
        let mut charge = 0i32;
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(v) = stack.pop() {
            size += 1;
            charge += g.atoms[v].formal_charge as i32;
            for &nb in &g.adjacency[v] {
                if is_lone_h(nb) && !seen[nb] {
                    seen[nb] = true;
                    stack.push(nb);
                }
            }
        }
        out.push((size, charge));
    }
    out
}

fn is_bare_proton(&(size, charge): &(usize, i32)) -> bool {
    size == 1 && charge > 0
}

pub(crate) fn hydrogen_component_sizes(g: &MoleculeGraph) -> Vec<usize> {
    hydrogen_components(g)
        .into_iter()
        .filter(|c| !is_bare_proton(c))
        .map(|(size, _)| size)
        .collect()
}

/// [`hydrogen_component_sizes`] と同じ走査順で、H だけの各連結成分の
/// 正味電荷を返す。ほとんどの場合 0 (金属水素化物切断由来・水素分子) だが、
/// `[H-]` のように入力 SMILES が孤立 H 自体に電荷を持たせるケースがある
/// (`[H-].C(=O)O[O-].[K+]` 等) — その場合だけ /q 層にこの成分の電荷を出す
/// 必要がある。
pub(crate) fn hydrogen_component_charges(g: &MoleculeGraph) -> Vec<i32> {
    hydrogen_components(g)
        .into_iter()
        .filter(|c| !is_bare_proton(c))
        .map(|(_, charge)| charge)
        .collect()
}

/// 裸のプロトン (他に結合先の重原子を持たない `[H+]`) の電荷合計。
/// `/p` 層に繰り込む (I64)。
pub(crate) fn bare_proton_charge(g: &MoleculeGraph) -> i32 {
    hydrogen_components(g)
        .into_iter()
        .filter(is_bare_proton)
        .map(|(_, charge)| charge)
        .sum()
}

/// H だけの成分の式 (サイズ 1 なら `H`、2 なら `H2`)。
pub(crate) fn hydrogen_component_formula(size: usize) -> String {
    if size > 1 {
        format!("H{size}")
    } else {
        "H".to_string()
    }
}

/// H だけの成分の h 層 (骨格原子 1 個に残り (size-1) 個がぶら下がる)。
/// サイズ 1 (孤立 H 原子) は h 層に何も出さない。
pub(crate) fn hydrogen_component_h_layer(size: usize) -> String {
    match size {
        0 | 1 => String::new(),
        2 => "1H".to_string(),
        k => format!("1H{}", k - 1),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    #[test]
    fn metal_carbon_bonds_break_without_charges() {
        let g = build_molecule_graph("C[Hg]C").unwrap();
        let d = disconnect_metals(&g).graph;
        assert!(d.atoms.iter().all(|a| a.formal_charge == 0));
        // Hg が孤立している
        let hg = d.atoms.iter().position(|a| a.symbol == "Hg").unwrap();
        assert!(d.adjacency[hg].is_empty());
    }

    #[test]
    fn metal_halogen_bond_breaks_heterolytically() {
        let g = build_molecule_graph("C[Hg]Cl").unwrap();
        let d = disconnect_metals(&g).graph;
        let hg = d.atoms.iter().position(|a| a.symbol == "Hg").unwrap();
        let cl = d.atoms.iter().position(|a| a.symbol == "Cl").unwrap();
        assert_eq!(d.atoms[hg].formal_charge, 1);
        assert_eq!(d.atoms[cl].formal_charge, -1);
    }

    #[test]
    fn metal_chalcogen_double_bond_breaks_neutral() {
        // 多重結合の金属=カルコゲンは中性で切る (I71)。
        // 実 InChI: `S=[Mo]` → `Mo.S` (電荷なし)、単結合の `CO[Mo]` は
        // 従来どおり異方開裂して `/q-1;+1`。
        let g = build_molecule_graph("S=[Mo]").unwrap();
        let d = disconnect_metals(&g).graph;
        assert!(d.atoms.iter().all(|a| a.formal_charge == 0));
        let s = d.atoms.iter().position(|a| a.symbol == "S").unwrap();
        assert!(d.adjacency[s].is_empty());

        // 単結合カルコゲン化物は引き続き異方開裂 (回帰防止)。
        let g = build_molecule_graph("CO[Mo]").unwrap();
        let d = disconnect_metals(&g);
        let o = d.graph.atoms.iter().position(|a| a.symbol == "O").unwrap();
        let mo = d.graph.atoms.iter().position(|a| a.symbol == "Mo").unwrap();
        assert_eq!(d.graph.atoms[o].formal_charge, -1);
        assert_eq!(d.graph.atoms[mo].formal_charge, 1);
        assert!(d.metal_locked_ligands.contains(&o));
    }

    #[test]
    fn metal_hydride_yields_lone_hydrogens() {
        let g = build_molecule_graph("[BiH3]").unwrap();
        let d = disconnect_metals(&g).graph;
        assert_eq!(hydrogen_component_sizes(&d), vec![1, 1, 1]);
        assert!(d.atoms.iter().all(|a| a.formal_charge == 0));
    }

    #[test]
    fn non_metals_are_left_connected() {
        for smi in ["Br[SiH2]C", "CC[AsH2]", "CB(O)O"] {
            let g = build_molecule_graph(smi).unwrap();
            let d = disconnect_metals(&g).graph;
            assert_eq!(d.bonds.len(), g.bonds.len(), "{smi}");
        }
    }
}
