//! トポロジー由来の理想幾何パラメータ (RUST_3D_PLAN.md C2)。
//!
//! - 混成の推定 (結合次数・芳香族から sp/sp2/sp3)
//! - 理想結合長: 代表的な結合はテーブル値 (実測文献値)、
//!   それ以外は Pyykkö 共有結合半径の和にフォールバック
//! - 理想結合角: 混成から (sp3 109.47° / sp2 120° / sp 180°)
//! - vdW 半径 (Bondi): 非結合下界に使う

use crate::graph::MoleculeGraph;

/// 混成状態。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Hybridization {
    Sp,
    Sp2,
    Sp3,
}

/// 全原子 (付加 H 含む) の混成を推定する。
/// H・ハロゲン等の末端原子も便宜上 Sp3 を返す (使われない)。
pub fn perceive_hybridization(g: &MoleculeGraph) -> Vec<Hybridization> {
    let n = g.atoms.len();
    let mut hyb = vec![Hybridization::Sp3; n];
    // 原子ごとの結合次数を集める
    let mut orders: Vec<Vec<f64>> = vec![Vec::new(); n];
    for b in &g.bonds {
        orders[b.begin_idx].push(b.bond_order);
        orders[b.end_idx].push(b.bond_order);
    }
    for i in 0..n {
        let a = &g.atoms[i];
        let deg = orders[i].len();
        let n_double = orders[i].iter().filter(|&&o| o == 2.0).count();
        let n_triple = orders[i].iter().filter(|&&o| o == 3.0).count();
        hyb[i] = if a.is_aromatic {
            Hybridization::Sp2
        } else if n_triple > 0 || (n_double >= 2 && deg <= 2) {
            // 三重結合、またはクムレン中心 (=C=) は直線。
            // スルホン S(=O)(=O) のような次数 4 の多重結合原子は含めない
            Hybridization::Sp
        } else if n_double >= 1 && deg <= 3 && matches!(a.symbol.as_str(), "C" | "N" | "O" | "B") {
            // 平面 sp2 は C/N/O/B のみ。S=O や P=O は錐体 (sp3 扱い) にする
            // (スルホキシドの S キラリティを平面化で壊さないため)
            Hybridization::Sp2
        } else {
            Hybridization::Sp3
        };
    }
    hyb
}

/// 混成から理想結合角 (ラジアン)。
pub fn ideal_angle(hyb: Hybridization) -> f64 {
    match hyb {
        Hybridization::Sp => std::f64::consts::PI,
        Hybridization::Sp2 => 120.0_f64.to_radians(),
        Hybridization::Sp3 => 109.471_f64.to_radians(),
    }
}

/// 元素・配位数補正つきの理想結合角 (ラジアン)。UFF θ0 に整合する値。
/// 低配位の第 16/15 族は四面体角より狭い (H2S 92°, ホスフィン 94°) が、
/// 超原子価 (スルホン S、リン酸 P など 4 配位) は四面体角に戻る。
pub fn ideal_angle_for(symbol: &str, aromatic: bool, hyb: Hybridization, degree: usize) -> f64 {
    let deg: f64 = match (symbol, aromatic, hyb) {
        ("S", true, _) => 92.2, // チオフェン型 (UFF S_R)
        ("Se" | "Te", true, _) => 90.0,
        ("O", true, _) => 110.0, // フラン型 (UFF O_R)
        ("S" | "Se" | "Te", false, Hybridization::Sp3) => match degree {
            0..=2 => 95.0, // スルフィド
            3 => 103.2,    // スルホキシド (UFF S_3+4)
            _ => 109.471,  // スルホン等の 4 配位
        },
        ("P" | "As" | "Sb" | "Bi", false, Hybridization::Sp3) => {
            if degree <= 3 {
                95.0 // ホスフィン型
            } else {
                109.471 // リン酸型 (UFF P_3+5)
            }
        }
        ("O", false, Hybridization::Sp3) => 104.51,
        ("N", false, Hybridization::Sp3) if degree <= 3 => 106.7,
        _ => return ideal_angle(hyb),
    };
    deg.to_radians()
}

/// 結合次数クラス (理想結合長の索引用)。
fn order_class(order: f64) -> u8 {
    if order == 1.5 {
        4 // 芳香族
    } else if order == 2.0 {
        2
    } else if order == 3.0 {
        3
    } else {
        1
    }
}

/// 代表的な結合の実測結合長 (Å)。キーは (記号昇順ペア, 次数クラス)。
/// 次数クラス: 1=単, 2=二重, 3=三重, 4=芳香族。
/// 出典: CRC Handbook / Allen et al. の典型値。
static BOND_LENGTHS: &[(&str, &str, u8, f64)] = &[
    ("C", "C", 1, 1.535),
    ("C", "C", 2, 1.339),
    ("C", "C", 3, 1.203),
    ("C", "C", 4, 1.394),
    ("C", "H", 1, 1.089),
    ("C", "N", 1, 1.469),
    ("C", "N", 2, 1.279),
    ("C", "N", 3, 1.158),
    ("C", "N", 4, 1.339),
    ("C", "O", 1, 1.426),
    ("C", "O", 2, 1.210),
    ("C", "O", 4, 1.370),
    ("C", "S", 1, 1.812),
    ("C", "S", 2, 1.600),
    ("C", "S", 4, 1.720),
    ("C", "F", 1, 1.350),
    ("C", "Cl", 1, 1.767),
    ("Br", "C", 1, 1.938),
    ("C", "I", 1, 2.139),
    ("C", "P", 1, 1.840),
    ("C", "Si", 1, 1.863),
    ("B", "C", 1, 1.560),
    ("C", "Se", 1, 1.970),
    ("H", "N", 1, 1.012),
    ("H", "O", 1, 0.962),
    ("H", "S", 1, 1.340),
    ("H", "Si", 1, 1.480),
    ("H", "P", 1, 1.420),
    ("N", "N", 1, 1.450),
    ("N", "N", 2, 1.252),
    ("N", "N", 3, 1.098),
    ("N", "N", 4, 1.350),
    ("N", "O", 1, 1.404),
    ("N", "O", 2, 1.212),
    ("O", "O", 1, 1.469),
    ("O", "P", 1, 1.630),
    ("O", "P", 2, 1.480),
    ("O", "S", 1, 1.658),
    ("O", "S", 2, 1.440),
    ("S", "S", 1, 2.048),
    ("N", "S", 1, 1.710),
    ("N", "P", 1, 1.700),
];

/// Pyykkö (2009) 共有結合半径 (Å)。フォールバック用。[単結合, 二重, 三重]
fn covalent_radii(symbol: &str) -> [f64; 3] {
    match symbol {
        "H" => [0.32, 0.32, 0.32],
        "B" => [0.85, 0.78, 0.73],
        "C" => [0.75, 0.67, 0.60],
        "N" => [0.71, 0.60, 0.54],
        "O" => [0.63, 0.57, 0.53],
        "F" => [0.64, 0.59, 0.53],
        "Na" => [1.55, 1.60, 1.60],
        "Mg" => [1.39, 1.32, 1.27],
        "Al" => [1.26, 1.13, 1.11],
        "Si" => [1.16, 1.07, 1.02],
        "P" => [1.11, 1.02, 0.94],
        "S" => [1.03, 0.94, 0.95],
        "Cl" => [0.99, 0.95, 0.93],
        "K" => [1.96, 1.93, 1.93],
        "Ca" => [1.71, 1.47, 1.33],
        "Fe" => [1.16, 1.09, 1.02],
        "Zn" => [1.18, 1.20, 1.20],
        "Ge" => [1.21, 1.11, 1.14],
        "As" => [1.21, 1.14, 1.06],
        "Se" => [1.16, 1.07, 1.07],
        "Br" => [1.14, 1.09, 1.10],
        "Sn" => [1.40, 1.30, 1.32],
        "Sb" => [1.40, 1.33, 1.27],
        "Te" => [1.36, 1.28, 1.21],
        "I" => [1.33, 1.29, 1.25],
        "Hg" => [1.32, 1.42, 1.42],
        "Pb" => [1.44, 1.35, 1.37],
        "Bi" => [1.51, 1.41, 1.35],
        _ => [1.4, 1.3, 1.3], // 未知元素の保守的なデフォルト
    }
}

/// 2 原子間の理想結合長 (Å)。
pub fn ideal_bond_length(sym_a: &str, sym_b: &str, order: f64) -> f64 {
    let oc = order_class(order);
    let (lo, hi) = if sym_a <= sym_b {
        (sym_a, sym_b)
    } else {
        (sym_b, sym_a)
    };
    for &(a, b, c, len) in BOND_LENGTHS {
        if a == lo && b == hi && c == oc {
            return len;
        }
    }
    // フォールバック: Pyykkö 半径の和 (芳香族は単・二重の平均)
    let ra = covalent_radii(sym_a);
    let rb = covalent_radii(sym_b);
    match oc {
        2 => ra[1] + rb[1],
        3 => ra[2] + rb[2],
        4 => (ra[0] + rb[0] + ra[1] + rb[1]) / 2.0,
        _ => ra[0] + rb[0],
    }
}

/// Bondi vdW 半径 (Å)。非結合原子対の下界に使う。
pub fn vdw_radius(symbol: &str) -> f64 {
    match symbol {
        "H" => 1.20,
        "B" => 1.92,
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "F" => 1.47,
        "Si" => 2.10,
        "P" => 1.80,
        "S" => 1.80,
        "Cl" => 1.75,
        "Ge" => 2.11,
        "As" => 1.85,
        "Se" => 1.90,
        "Br" => 1.85,
        "Sn" => 2.17,
        "Sb" => 2.06,
        "Te" => 2.06,
        "I" => 1.98,
        _ => 1.80,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    #[test]
    fn representative_bond_lengths_within_tolerance() {
        // 完了条件: 代表結合の理想長が文献値 ±0.03 Å
        let cases: &[(&str, &str, f64, f64)] = &[
            ("C", "C", 1.0, 1.54),
            ("C", "C", 2.0, 1.34),
            ("C", "C", 3.0, 1.20),
            ("C", "C", 1.5, 1.39),
            ("C", "H", 1.0, 1.09),
            ("C", "N", 1.0, 1.47),
            ("C", "N", 2.0, 1.28),
            ("C", "N", 3.0, 1.16),
            ("C", "N", 1.5, 1.34),
            ("C", "O", 1.0, 1.43),
            ("C", "O", 2.0, 1.21),
            ("C", "O", 1.5, 1.37),
            ("C", "S", 1.0, 1.81),
            ("C", "S", 2.0, 1.60),
            ("C", "F", 1.0, 1.35),
            ("C", "Cl", 1.0, 1.77),
            ("C", "Br", 1.0, 1.94),
            ("C", "I", 1.0, 2.14),
            ("C", "P", 1.0, 1.84),
            ("C", "Si", 1.0, 1.86),
            ("N", "H", 1.0, 1.01),
            ("O", "H", 1.0, 0.96),
            ("S", "H", 1.0, 1.34),
            ("N", "N", 1.0, 1.45),
            ("N", "N", 2.0, 1.25),
            ("N", "O", 1.0, 1.40),
            ("N", "O", 2.0, 1.21),
            ("O", "O", 1.0, 1.48),
            ("P", "O", 2.0, 1.48),
            ("S", "O", 2.0, 1.44),
            ("S", "S", 1.0, 2.05),
        ];
        for &(a, b, order, expect) in cases {
            let got = ideal_bond_length(a, b, order);
            assert!(
                (got - expect).abs() <= 0.03,
                "{a}-{b} (order {order}): got {got}, expected {expect}"
            );
            // 引数順に依らない
            assert_eq!(got, ideal_bond_length(b, a, order));
        }
    }

    #[test]
    fn fallback_is_sane() {
        // テーブルにないペアもフォールバックで妥当な範囲に入る
        let l = ideal_bond_length("Se", "Se", 1.0);
        assert!((2.0..2.7).contains(&l), "Se-Se = {l}");
        let l = ideal_bond_length("Sn", "H", 1.0);
        assert!((1.5..2.0).contains(&l), "Sn-H = {l}");
    }

    #[test]
    fn hybridization_perception() {
        // エタン: 全 sp3
        let g = build_molecule_graph("CC").unwrap();
        let h = perceive_hybridization(&g);
        assert_eq!(h[0], Hybridization::Sp3);
        // エチレン: sp2
        let g = build_molecule_graph("C=C").unwrap();
        let h = perceive_hybridization(&g);
        assert_eq!(h[0], Hybridization::Sp2);
        // アセチレン・ニトリル: sp
        let g = build_molecule_graph("C#C").unwrap();
        assert_eq!(perceive_hybridization(&g)[0], Hybridization::Sp);
        let g = build_molecule_graph("CC#N").unwrap();
        let h = perceive_hybridization(&g);
        assert_eq!(h[1], Hybridization::Sp);
        assert_eq!(h[2], Hybridization::Sp);
        // クムレン中心
        let g = build_molecule_graph("C=C=C").unwrap();
        assert_eq!(perceive_hybridization(&g)[1], Hybridization::Sp);
        // ベンゼン: 芳香族 sp2
        let g = build_molecule_graph("c1ccccc1").unwrap();
        assert!(perceive_hybridization(&g)[..6]
            .iter()
            .all(|&h| h == Hybridization::Sp2));
        // カルボニル: C も O も sp2
        let g = build_molecule_graph("CC(C)=O").unwrap();
        let h = perceive_hybridization(&g);
        assert_eq!(h[1], Hybridization::Sp2);
        assert_eq!(h[3], Hybridization::Sp2);
        // スルホン S は四面体 (sp3)、スルホキシド S も錐体 (sp3)
        let g = build_molecule_graph("CS(=O)(=O)C").unwrap();
        assert_eq!(perceive_hybridization(&g)[1], Hybridization::Sp3);
        let g = build_molecule_graph("CS(=O)C").unwrap();
        assert_eq!(perceive_hybridization(&g)[1], Hybridization::Sp3);
        // リン酸 P も sp3
        let g = build_molecule_graph("OP(=O)(O)O").unwrap();
        assert_eq!(perceive_hybridization(&g)[1], Hybridization::Sp3);
    }

    #[test]
    fn ideal_angles() {
        assert!((ideal_angle(Hybridization::Sp3).to_degrees() - 109.471).abs() < 1e-9);
        assert!((ideal_angle(Hybridization::Sp2).to_degrees() - 120.0).abs() < 1e-9);
        assert!((ideal_angle(Hybridization::Sp).to_degrees() - 180.0).abs() < 1e-9);
    }
}
