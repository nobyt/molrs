//! UFF 力場 (RUST_3D_PLAN.md C7)。Rappé et al. JACS 114, 10024 (1992)。
//!
//! 実装する項: 結合伸縮 (調和)、結合角 (Fourier / 直線形)、トーション、
//! vdW (LJ 12-6, 1-2/1-3 除外)。静電項は省略 (RDKit の UFF 既定と同じ)。
//! 面外反転 (Wilson 角) 項の代わりに、C6 の体積ペナルティ (sp2 平面・
//! キラル符号) をそのまま併用する — 実装が単純で立体保存が確実。
//!
//! パラメータは RDKit Params.cpp (= UFF 論文 Table 1) から転記し、
//! 転記値は代表結合の平衡長テストで固定する。エネルギーは全て
//! cosθ / cosφ の多項式として書き、解析勾配は数値微分テストで検証する。
//!
//! RDKit との既知の相違:
//! - チオカルボニル S を S_2 でなく S_3+2 と型付けする (S_2 だと C=S が
//!   1.43 Å と実測 1.60 Å から大きく外れるため。終端原子なので角度項に
//!   影響しない)
//! - アミド結合次数 1.41 の特例は未実装

use std::collections::HashSet;

use crate::conformer::exp_torsions::{collect_exp_torsions, ExpTorsion};
use crate::conformer::minimize::{minimize_with, signed_volume_of, VolumeConstraint};
use crate::conformer::params::{perceive_hybridization, Hybridization};
use crate::conformer::stereo3d::build_volume_constraints;
use crate::geometry::Vec3;
use crate::graph::MoleculeGraph;

/// UFF 原子パラメータ行: (label, r1, theta0[deg], x1, D1, Z1, V1, U1, chi)。
type UffRow = (&'static str, f64, f64, f64, f64, f64, f64, f64, f64);

static UFF_PARAMS: &[UffRow] = &[
    ("H_", 0.354, 180.0, 2.886, 0.044, 0.712, 0.0, 0.0, 4.528),
    ("Li", 1.336, 180.0, 2.451, 0.025, 1.026, 0.0, 2.0, 3.006),
    ("B_3", 0.838, 109.47, 4.083, 0.18, 1.755, 0.0, 2.0, 5.11),
    ("B_2", 0.828, 120.0, 4.083, 0.18, 1.755, 0.0, 2.0, 5.11),
    ("C_3", 0.757, 109.47, 3.851, 0.105, 1.912, 2.119, 2.0, 5.343),
    ("C_R", 0.729, 120.0, 3.851, 0.105, 1.912, 0.0, 2.0, 5.343),
    ("C_2", 0.732, 120.0, 3.851, 0.105, 1.912, 0.0, 2.0, 5.343),
    ("C_1", 0.706, 180.0, 3.851, 0.105, 1.912, 0.0, 2.0, 5.343),
    ("N_3", 0.7, 106.7, 3.66, 0.069, 2.544, 0.45, 2.0, 6.899),
    ("N_R", 0.699, 120.0, 3.66, 0.069, 2.544, 0.0, 2.0, 6.899),
    ("N_2", 0.685, 111.2, 3.66, 0.069, 2.544, 0.0, 2.0, 6.899),
    ("N_1", 0.656, 180.0, 3.66, 0.069, 2.544, 0.0, 2.0, 6.899),
    ("O_3", 0.658, 104.51, 3.5, 0.06, 2.3, 0.018, 2.0, 8.741),
    ("O_R", 0.68, 110.0, 3.5, 0.06, 2.3, 0.0, 2.0, 8.741),
    ("O_2", 0.634, 120.0, 3.5, 0.06, 2.3, 0.0, 2.0, 8.741),
    ("O_1", 0.639, 180.0, 3.5, 0.06, 2.3, 0.0, 2.0, 8.741),
    ("F_", 0.668, 180.0, 3.364, 0.05, 1.735, 0.0, 2.0, 10.874),
    ("Na", 1.539, 180.0, 2.983, 0.03, 1.081, 0.0, 1.25, 2.843),
    (
        "Mg3+2", 1.421, 109.47, 3.021, 0.111, 1.787, 0.0, 1.25, 3.951,
    ),
    ("Al3", 1.244, 109.47, 4.499, 0.505, 1.792, 0.0, 1.25, 4.06),
    (
        "Si3", 1.117, 109.47, 4.295, 0.402, 2.323, 1.225, 1.25, 4.168,
    ),
    ("P_3+3", 1.101, 93.8, 4.147, 0.305, 2.863, 2.4, 1.25, 5.463),
    (
        "P_3+5", 1.056, 109.47, 4.147, 0.305, 2.863, 2.4, 1.25, 5.463,
    ),
    (
        "S_3+2", 1.064, 92.1, 4.035, 0.274, 2.703, 0.484, 1.25, 6.928,
    ),
    (
        "S_3+4", 1.049, 103.2, 4.035, 0.274, 2.703, 0.484, 1.25, 6.928,
    ),
    (
        "S_3+6", 1.027, 109.47, 4.035, 0.274, 2.703, 0.484, 1.25, 6.928,
    ),
    ("S_R", 1.077, 92.2, 4.035, 0.274, 2.703, 0.0, 1.25, 6.928),
    ("Cl", 1.044, 180.0, 3.947, 0.227, 2.348, 0.0, 1.25, 8.564),
    ("K_", 1.953, 180.0, 3.812, 0.035, 1.165, 0.0, 0.7, 2.421),
    ("Ca6+2", 1.761, 90.0, 3.399, 0.238, 2.141, 0.0, 0.7, 3.231),
    ("Fe6+2", 1.335, 90.0, 2.912, 0.013, 2.43, 0.0, 0.7, 3.76),
    ("Zn3+2", 1.193, 109.47, 2.763, 0.124, 1.308, 0.0, 0.7, 5.106),
    ("Ge3", 1.197, 109.47, 4.28, 0.379, 2.789, 0.701, 0.7, 4.051),
    ("As3+3", 1.211, 92.1, 4.23, 0.309, 2.864, 1.5, 0.7, 5.188),
    ("Se3+2", 1.19, 90.6, 4.205, 0.291, 2.764, 0.335, 0.7, 6.428),
    ("Br", 1.192, 180.0, 4.189, 0.251, 2.519, 0.0, 0.7, 7.79),
    ("Sn3", 1.398, 109.47, 4.392, 0.567, 2.961, 0.199, 0.2, 3.987),
    ("Sb3+3", 1.407, 91.6, 4.42, 0.449, 2.704, 1.1, 0.2, 4.899),
    ("Te3+2", 1.386, 90.25, 4.47, 0.398, 2.882, 0.3, 0.2, 5.816),
    ("I_", 1.382, 180.0, 4.5, 0.339, 2.65, 0.0, 0.2, 6.822),
    ("Hg1+2", 1.34, 180.0, 2.705, 0.385, 1.75, 0.0, 0.1, 6.27),
    ("Pb3", 1.459, 109.47, 4.297, 0.663, 2.846, 0.1, 0.1, 3.9),
    ("Bi3+3", 1.512, 90.0, 4.37, 0.518, 2.47, 1.0, 0.1, 4.69),
];

#[derive(Clone, Copy)]
struct P {
    r1: f64,
    theta0: f64, // rad
    x1: f64,
    d1: f64,
    z1: f64,
    v1: f64,
    u1: f64,
    chi: f64,
}

fn param(label: &str) -> Option<P> {
    UFF_PARAMS
        .iter()
        .find(|t| t.0 == label)
        .map(|&(_, r1, th, x1, d1, z1, v1, u1, chi)| P {
            r1,
            theta0: th.to_radians(),
            x1,
            d1,
            z1,
            v1,
            u1,
            chi,
        })
}

/// 原子型の割当て。未対応元素があれば None (UFF をスキップ)。
fn assign_types(g: &MoleculeGraph, hyb: &[Hybridization]) -> Option<Vec<&'static str>> {
    let mut valence = vec![0.0f64; g.atoms.len()];
    for b in &g.bonds {
        valence[b.begin_idx] += b.bond_order;
        valence[b.end_idx] += b.bond_order;
    }
    g.atoms
        .iter()
        .map(|a| {
            let t: &'static str = match a.symbol.as_str() {
                "H" => "H_",
                "C" => match (a.is_aromatic, hyb[a.idx]) {
                    (true, _) => "C_R",
                    (_, Hybridization::Sp) => "C_1",
                    (_, Hybridization::Sp2) => "C_2",
                    _ => "C_3",
                },
                "N" => match (a.is_aromatic, hyb[a.idx]) {
                    (true, _) => "N_R",
                    (_, Hybridization::Sp) => "N_1",
                    (_, Hybridization::Sp2) => "N_2",
                    _ => "N_3",
                },
                "O" => match (a.is_aromatic, hyb[a.idx]) {
                    (true, _) => "O_R",
                    (_, Hybridization::Sp) => "O_1",
                    (_, Hybridization::Sp2) => "O_2",
                    _ => "O_3",
                },
                "F" => "F_",
                "Cl" => "Cl",
                "Br" => "Br",
                "I" => "I_",
                "B" => {
                    if a.is_aromatic || hyb[a.idx] == Hybridization::Sp2 {
                        "B_2"
                    } else {
                        "B_3"
                    }
                }
                "Si" => "Si3",
                "P" => {
                    if valence[a.idx] >= 4.5 {
                        "P_3+5"
                    } else {
                        "P_3+3"
                    }
                }
                "S" => {
                    if a.is_aromatic {
                        "S_R"
                    } else if valence[a.idx] >= 5.5 {
                        "S_3+6"
                    } else if valence[a.idx] >= 3.5 {
                        "S_3+4"
                    } else {
                        "S_3+2" // チオカルボニル S もこちら (冒頭コメント参照)
                    }
                }
                "Se" => "Se3+2",
                "Te" => "Te3+2",
                "As" => "As3+3",
                "Sb" => "Sb3+3",
                "Bi" => "Bi3+3",
                "Ge" => "Ge3",
                "Sn" => "Sn3",
                "Pb" => "Pb3",
                "Al" => "Al3",
                "Na" => "Na",
                "K" => "K_",
                "Li" => "Li",
                "Mg" => "Mg3+2",
                "Ca" => "Ca6+2",
                "Fe" => "Fe6+2",
                "Zn" => "Zn3+2",
                "Hg" => "Hg1+2",
                _ => return None,
            };
            Some(t)
        })
        .collect()
}

/// UFF 自然結合長 (RDKit calcBondRestLength 相当)。
fn rest_length(pi: &P, pj: &P, bond_order: f64) -> f64 {
    let r_bo = -0.1332 * (pi.r1 + pj.r1) * bond_order.ln();
    let dchi = pi.chi.sqrt() - pj.chi.sqrt();
    let r_en = pi.r1 * pj.r1 * dchi * dchi / (pi.chi * pi.r1 + pj.chi * pj.r1);
    pi.r1 + pj.r1 + r_bo - r_en
}

enum AngleForm {
    /// E = K (1 + cosθ) — 直線 (θ0 = 180°)
    Linear,
    /// E = K (C0 + C1 cosθ + C2 cos2θ)
    Fourier { c0: f64, c1: f64, c2: f64 },
}

struct AngleTerm {
    i: usize,
    j: usize,
    k: usize,
    ka: f64,
    form: AngleForm,
}

struct TorsionTerm {
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    /// V/2 (経路数で除算済み)
    v_half: f64,
    n: u8,
    cos_n_phi0: f64,
}

/// UFF エネルギー場。
pub(crate) struct UffField {
    bonds: Vec<(usize, usize, f64, f64)>, // (i, j, r0, kb)
    angles: Vec<AngleTerm>,
    torsions: Vec<TorsionTerm>,
    /// ETKDG 実験トーション (C10)。中央結合が一致する UFF 汎用トーションは
    /// 二重計上を避けるため生成しない
    exp_torsions: Vec<ExpTorsion>,
    vdw: Vec<(usize, usize, f64, f64)>, // (i, j, x_ij, d_ij)
    volumes: Vec<VolumeConstraint>,
}

/// 分子から UFF 場を構築する。未対応元素・孤立原子は None。
/// `use_exp_torsions` で ETKDG 実験トーション (C10) を有効化する。
pub(crate) fn build_uff(g: &MoleculeGraph, use_exp_torsions: bool) -> Option<UffField> {
    let exp = if use_exp_torsions {
        collect_exp_torsions(g)
    } else {
        Vec::new()
    };
    let exp_central: HashSet<(usize, usize)> = exp
        .iter()
        .map(|t| (t.atoms[1].min(t.atoms[2]), t.atoms[1].max(t.atoms[2])))
        .collect();
    let n = g.atoms.len();
    if n < 2 {
        return None;
    }
    let hyb = perceive_hybridization(g);
    let types = assign_types(g, &hyb)?;
    let p: Vec<P> = types
        .iter()
        .map(|t| param(t).expect("param exists"))
        .collect();

    let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for b in &g.bonds {
        adj[b.begin_idx].push((b.end_idx, b.bond_order));
        adj[b.end_idx].push((b.begin_idx, b.bond_order));
    }

    // ---- 結合伸縮 ----
    let mut bonds = Vec::with_capacity(g.bonds.len());
    for b in &g.bonds {
        let (i, j) = (b.begin_idx, b.end_idx);
        let r0 = rest_length(&p[i], &p[j], b.bond_order);
        let kb = 2.0 * 332.06 * p[i].z1 * p[j].z1 / (r0 * r0 * r0);
        bonds.push((i, j, r0, kb));
    }

    // ---- 結合角 ----
    let mut angles = Vec::new();
    for j in 0..n {
        for a in 0..adj[j].len() {
            for b in (a + 1)..adj[j].len() {
                let (i, o_ij) = adj[j][a];
                let (k, o_jk) = adj[j][b];
                let theta0 = p[j].theta0;
                let cos_t0 = theta0.cos();
                let r_ij = rest_length(&p[i], &p[j], o_ij);
                let r_jk = rest_length(&p[j], &p[k], o_jk);
                let r_ik = (r_ij * r_ij + r_jk * r_jk - 2.0 * r_ij * r_jk * cos_t0).sqrt();
                let beta = 2.0 * 332.06 / (r_ij * r_jk);
                let pre = beta * p[i].z1 * p[k].z1 / r_ik.powi(5);
                let r_term = r_ij * r_jk;
                let inner = 3.0 * r_term * (1.0 - cos_t0 * cos_t0) - r_ik * r_ik * cos_t0;
                let ka = pre * r_term * inner;
                let form = if theta0.to_degrees() > 175.0 {
                    AngleForm::Linear
                } else {
                    let sin2 = 1.0 - cos_t0 * cos_t0;
                    let c2 = 1.0 / (4.0 * sin2);
                    let c1 = -4.0 * c2 * cos_t0;
                    let c0 = c2 * (2.0 * cos_t0 * cos_t0 + 1.0);
                    AngleForm::Fourier { c0, c1, c2 }
                };
                angles.push(AngleTerm { i, j, k, ka, form });
            }
        }
    }

    // ---- トーション ----
    let group16 = |i: usize| matches!(g.atoms[i].symbol.as_str(), "O" | "S" | "Se" | "Te");
    let sp2ish = |i: usize| g.atoms[i].is_aromatic || hyb[i] == Hybridization::Sp2;
    let mut torsions = Vec::new();
    for b in &g.bonds {
        let (j, k) = (b.begin_idx, b.end_idx);
        if exp_central.contains(&(j.min(k), j.max(k))) {
            continue; // 実験トーションが受け持つ結合
        }
        if hyb[j] == Hybridization::Sp || hyb[k] == Hybridization::Sp {
            continue;
        }
        let nj = adj[j].len();
        let nk = adj[k].len();
        if nj < 2 || nk < 2 {
            continue;
        }
        let n_paths = ((nj - 1) * (nk - 1)) as f64;

        // 中央結合の (V, n, φ0)
        let (v, order, cos_n_phi0) = match (sp2ish(j), sp2ish(k)) {
            (false, false) => {
                if group16(j) && group16(k) {
                    // 両端が第 16 族 sp3: n=2, φ0=90°
                    let vj: f64 = if g.atoms[j].symbol == "O" { 2.0 } else { 6.8 };
                    let vk: f64 = if g.atoms[k].symbol == "O" { 2.0 } else { 6.8 };
                    ((vj * vk).sqrt(), 2u8, -1.0) // cos(2·90°) = −1
                } else {
                    // sp3-sp3: n=3, φ0=60° (cos 180° = −1)
                    ((p[j].v1 * p[k].v1).sqrt(), 3, -1.0)
                }
            }
            (true, true) => {
                // sp2-sp2: n=2, φ0=180° (cos 360° = +1)
                let v = 5.0 * (p[j].u1 * p[k].u1).sqrt() * (1.0 + 4.18 * b.bond_order.ln());
                (v, 2, 1.0)
            }
            (sj, _) => {
                // sp2-sp3 混合
                let sp2_atom = if sj { j } else { k };
                let another_sp2 = adj[sp2_atom]
                    .iter()
                    .any(|&(x, _)| x != j && x != k && sp2ish(x));
                if another_sp2 {
                    // プロペン型特例: n=3, φ0=180° (cos 540° = −1)
                    (2.0, 3, -1.0)
                } else {
                    // 一般: n=6, φ0=0° (cos 0 = 1)
                    (1.0, 6, 1.0)
                }
            }
        };
        if v.abs() < 1e-9 {
            continue;
        }
        let v_half = 0.5 * v / n_paths;
        for &(i, _) in &adj[j] {
            if i == k {
                continue;
            }
            for &(l, _) in &adj[k] {
                if l == j || l == i {
                    continue;
                }
                torsions.push(TorsionTerm {
                    i,
                    j,
                    k,
                    l,
                    v_half,
                    n: order,
                    cos_n_phi0,
                });
            }
        }
    }

    // ---- vdW (1-2, 1-3 除外) ----
    let mut excluded: HashSet<(usize, usize)> = HashSet::new();
    for b in &g.bonds {
        let (i, j) = (b.begin_idx.min(b.end_idx), b.begin_idx.max(b.end_idx));
        excluded.insert((i, j));
    }
    for nbrs in &adj {
        for a in 0..nbrs.len() {
            for b in (a + 1)..nbrs.len() {
                let (i, k) = (nbrs[a].0, nbrs[b].0);
                excluded.insert((i.min(k), i.max(k)));
            }
        }
    }
    let mut vdw = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if excluded.contains(&(i, j)) {
                continue;
            }
            vdw.push((i, j, (p[i].x1 * p[j].x1).sqrt(), (p[i].d1 * p[j].d1).sqrt()));
        }
    }

    Some(UffField {
        bonds,
        angles,
        torsions,
        exp_torsions: exp,
        vdw,
        // sp2 平面・芳香環・キラル符号は C6 の体積ペナルティで維持する
        volumes: build_volume_constraints(g),
    })
}

/// cos(nφ) とその dcos(nφ)/dcosφ (チェビシェフ多項式)。
fn cos_n(c: f64, n: u8) -> (f64, f64) {
    let c2 = c * c;
    match n {
        1 => (c, 1.0),
        2 => (2.0 * c2 - 1.0, 4.0 * c),
        3 => (4.0 * c2 * c - 3.0 * c, 12.0 * c2 - 3.0),
        4 => (8.0 * c2 * c2 - 8.0 * c2 + 1.0, 32.0 * c2 * c - 16.0 * c),
        5 => (
            16.0 * c2 * c2 * c - 20.0 * c2 * c + 5.0 * c,
            80.0 * c2 * c2 - 60.0 * c2 + 5.0,
        ),
        6 => (
            32.0 * c2 * c2 * c2 - 48.0 * c2 * c2 + 18.0 * c2 - 1.0,
            192.0 * c2 * c2 * c - 192.0 * c2 * c + 36.0 * c,
        ),
        _ => unreachable!("unsupported periodicity"),
    }
}

impl UffField {
    pub(crate) fn energy_and_grad(&self, coords: &[Vec3]) -> (f64, Vec<Vec3>) {
        let mut e = 0.0;
        let mut grad = vec![Vec3::ZERO; coords.len()];

        // 結合伸縮
        for &(i, j, r0, kb) in &self.bonds {
            let diff = coords[i] - coords[j];
            let d = diff.norm().max(1e-8);
            let delta = d - r0;
            e += 0.5 * kb * delta * delta;
            let gi = diff * (kb * delta / d);
            grad[i] = grad[i] + gi;
            grad[j] = grad[j] - gi;
        }

        // 結合角 (cosθ の多項式として評価)
        for t in &self.angles {
            let u = coords[t.i] - coords[t.j];
            let v = coords[t.k] - coords[t.j];
            let nu = u.norm().max(1e-8);
            let nv = v.norm().max(1e-8);
            let c = (u.dot(v) / (nu * nv)).clamp(-1.0, 1.0);
            let (energy, de_dc) = match t.form {
                AngleForm::Linear => (t.ka * (1.0 + c), t.ka),
                AngleForm::Fourier { c0, c1, c2 } => {
                    let cos2t = 2.0 * c * c - 1.0;
                    (
                        t.ka * (c0 + c1 * c + c2 * cos2t),
                        t.ka * (c1 + 4.0 * c2 * c),
                    )
                }
            };
            e += energy;
            // dc/du = v/(|u||v|) − c u/|u|²
            let dc_du = v / (nu * nv) - u * (c / (nu * nu));
            let dc_dv = u / (nu * nv) - v * (c / (nv * nv));
            grad[t.i] = grad[t.i] + dc_du * de_dc;
            grad[t.k] = grad[t.k] + dc_dv * de_dc;
            grad[t.j] = grad[t.j] - (dc_du + dc_dv) * de_dc;
        }

        // トーション (cosφ の多項式として評価)
        for t in &self.torsions {
            let b1 = coords[t.j] - coords[t.i];
            let b2 = coords[t.k] - coords[t.j];
            let b3 = coords[t.l] - coords[t.k];
            let n1 = b1.cross(b2);
            let n2 = b2.cross(b3);
            let m1 = n1.norm();
            let m2 = n2.norm();
            if m1 < 1e-6 || m2 < 1e-6 {
                continue; // 直線縮退
            }
            let c = (n1.dot(n2) / (m1 * m2)).clamp(-1.0, 1.0);
            let (cn, dcn_dc) = cos_n(c, t.n);
            e += t.v_half * (1.0 - t.cos_n_phi0 * cn);
            let de_dc = -t.v_half * t.cos_n_phi0 * dcn_dc;

            let dc_dn1 = n2 / (m1 * m2) - n1 * (c / (m1 * m1));
            let dc_dn2 = n1 / (m1 * m2) - n2 * (c / (m2 * m2));
            // n1 = b1×b2, n2 = b2×b3 より
            let dc_db1 = b2.cross(dc_dn1);
            let dc_db2 = dc_dn1.cross(b1) + b3.cross(dc_dn2);
            let dc_db3 = dc_dn2.cross(b2);
            grad[t.i] = grad[t.i] - dc_db1 * de_dc;
            grad[t.j] = grad[t.j] + (dc_db1 - dc_db2) * de_dc;
            grad[t.k] = grad[t.k] + (dc_db2 - dc_db3) * de_dc;
            grad[t.l] = grad[t.l] + dc_db3 * de_dc;
        }

        // ETKDG 実験トーション: E = Σ_{i=1..6} V_i (1 + s_i cos(iφ))
        for t in &self.exp_torsions {
            let b1 = coords[t.atoms[1]] - coords[t.atoms[0]];
            let b2 = coords[t.atoms[2]] - coords[t.atoms[1]];
            let b3 = coords[t.atoms[3]] - coords[t.atoms[2]];
            let n1 = b1.cross(b2);
            let n2 = b2.cross(b3);
            let m1 = n1.norm();
            let m2 = n2.norm();
            if m1 < 1e-6 || m2 < 1e-6 {
                continue;
            }
            let c = (n1.dot(n2) / (m1 * m2)).clamp(-1.0, 1.0);
            let mut de_dc = 0.0;
            for i in 0..6 {
                if t.v[i] == 0.0 {
                    continue;
                }
                let (cn, dcn) = cos_n(c, (i + 1) as u8);
                e += t.v[i] * (1.0 + t.signs[i] as f64 * cn);
                de_dc += t.v[i] * t.signs[i] as f64 * dcn;
            }
            let dc_dn1 = n2 / (m1 * m2) - n1 * (c / (m1 * m1));
            let dc_dn2 = n1 / (m1 * m2) - n2 * (c / (m2 * m2));
            let dc_db1 = b2.cross(dc_dn1);
            let dc_db2 = dc_dn1.cross(b1) + b3.cross(dc_dn2);
            let dc_db3 = dc_dn2.cross(b2);
            grad[t.atoms[0]] = grad[t.atoms[0]] - dc_db1 * de_dc;
            grad[t.atoms[1]] = grad[t.atoms[1]] + (dc_db1 - dc_db2) * de_dc;
            grad[t.atoms[2]] = grad[t.atoms[2]] + (dc_db2 - dc_db3) * de_dc;
            grad[t.atoms[3]] = grad[t.atoms[3]] + dc_db3 * de_dc;
        }

        // vdW (LJ 12-6)
        for &(i, j, x, d1) in &self.vdw {
            let diff = coords[i] - coords[j];
            let r = diff.norm().max(1e-8);
            if r > 12.0 {
                continue;
            }
            let q = x / r;
            let q6 = q.powi(6);
            let q12 = q6 * q6;
            e += d1 * (q12 - 2.0 * q6);
            // dE/dr = −12 D/r (q12 − q6)
            let de_dr = -12.0 * d1 * (q12 - q6) / r;
            let gi = diff * (de_dr / r);
            grad[i] = grad[i] + gi;
            grad[j] = grad[j] - gi;
        }

        // 体積ペナルティ (平面・キラル; kcal スケールに増幅)
        for vc in &self.volumes {
            let (v, gs) = signed_volume_of(coords, &vc.atoms);
            let viol = if v > vc.upper {
                v - vc.upper
            } else if v < vc.lower {
                v - vc.lower
            } else {
                continue;
            };
            let w = vc.weight * 10.0;
            e += w * viol * viol;
            let scale = 2.0 * w * viol;
            for (t, &ai) in vc.atoms.iter().enumerate() {
                grad[ai] = grad[ai] + gs[t] * scale;
            }
        }

        (e, grad)
    }
}

/// 各結合の UFF 平衡長 (g.bonds と同順)。
pub(crate) fn bond_rest_lengths(g: &MoleculeGraph) -> Option<Vec<f64>> {
    let hyb = perceive_hybridization(g);
    let types = assign_types(g, &hyb)?;
    let p: Vec<P> = types
        .iter()
        .map(|t| param(t).expect("param exists"))
        .collect();
    Some(
        g.bonds
            .iter()
            .map(|b| rest_length(&p[b.begin_idx], &p[b.end_idx], b.bond_order))
            .collect(),
    )
}

/// UFF 最適化。最終エネルギーを返す。
pub(crate) fn optimize(field: &UffField, coords: &mut [Vec3], max_iter: usize) -> f64 {
    minimize_with(&|c: &[Vec3]| field.energy_and_grad(c), coords, max_iter)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformer::{embed_molecule, EmbedParams};
    use crate::graph::build_molecule_graph;

    fn optimized(smiles: &str) -> (MoleculeGraph, Vec<Vec3>) {
        let g = build_molecule_graph(smiles).expect("valid");
        let conf = embed_molecule(&g, &EmbedParams::default()).expect("embeds");
        (g, conf.coords)
    }

    fn angle_deg(coords: &[Vec3], i: usize, j: usize, k: usize) -> f64 {
        let u = coords[i] - coords[j];
        let v = coords[k] - coords[j];
        (u.dot(v) / (u.norm() * v.norm())).acos().to_degrees()
    }

    /// 解析勾配 vs 数値微分 (中心差分)。転記・微分の全誤りを捕まえる。
    #[test]
    fn analytic_gradient_matches_numerical() {
        for smi in [
            "CCO",
            "c1ccccc1",
            "C/C=C/C",
            "CS(=O)(=O)C",
            "CC(=O)N",
            "COC",
        ] {
            let g = build_molecule_graph(smi).unwrap();
            let field = build_uff(&g, true).expect("typed");
            // DG の初期座標 (最適化前 = 勾配が大きい点) で比較
            let bm = crate::conformer::bounds::build_bounds(&g);
            let mut rng = crate::geometry::SeededRng::new(9);
            let coords = crate::conformer::embed::embed_from_bounds(&bm, &mut rng).unwrap();
            let (_, grad) = field.energy_and_grad(&coords);
            let h = 1e-6;
            for ai in 0..coords.len() {
                for dim in 0..3 {
                    let mut cp = coords.clone();
                    let mut cm = coords.clone();
                    match dim {
                        0 => {
                            cp[ai].x += h;
                            cm[ai].x -= h;
                        }
                        1 => {
                            cp[ai].y += h;
                            cm[ai].y -= h;
                        }
                        _ => {
                            cp[ai].z += h;
                            cm[ai].z -= h;
                        }
                    }
                    let num =
                        (field.energy_and_grad(&cp).0 - field.energy_and_grad(&cm).0) / (2.0 * h);
                    let ana = match dim {
                        0 => grad[ai].x,
                        1 => grad[ai].y,
                        _ => grad[ai].z,
                    };
                    let tol = 1e-4 + 1e-4 * num.abs();
                    assert!(
                        (num - ana).abs() < tol,
                        "{smi}: atom {ai} dim {dim}: numerical {num:.6} vs analytic {ana:.6}"
                    );
                }
            }
        }
    }

    #[test]
    fn ethane_geometry() {
        let (g, coords) = optimized("CC");
        // UFF の C_3-C_3 平衡長は 1.514 Å (rEN = 0)
        let d = coords[0].distance(coords[1]);
        assert!((d - 1.514).abs() < 0.02, "C-C = {d:.3}");
        // C-H
        for b in &g.bonds {
            if g.atoms[b.end_idx].symbol == "H" {
                let dh = coords[b.begin_idx].distance(coords[b.end_idx]);
                assert!((dh - 1.11).abs() < 0.03, "C-H = {dh:.3}");
            }
        }
        // H-C-C 角 ≈ 109.5°〜111.5°
        let a = angle_deg(&coords, 2, 0, 1);
        assert!((a - 110.0).abs() < 3.0, "H-C-C = {a:.1}");
    }

    #[test]
    fn benzene_geometry() {
        let (_, coords) = optimized("c1ccccc1");
        for i in 0..6 {
            let d = coords[i].distance(coords[(i + 1) % 6]);
            // UFF の平衡長 r0 = 1.379 Å だが、パラ位 1-4 の LJ 反発で
            // 環は 1.398 Å 前後まで膨張する (RDKit UFF と同じ挙動で、
            // 実測 1.394 Å にむしろ近い)
            assert!((d - 1.398).abs() < 0.015, "ring C-C = {d:.3}");
            let a = angle_deg(&coords, (i + 5) % 6, i, (i + 1) % 6);
            assert!((a - 120.0).abs() < 2.0, "ring angle = {a:.1}");
        }
    }

    #[test]
    fn cyclohexane_geometry() {
        let (_, coords) = optimized("C1CCCCC1");
        let mut sum = 0.0;
        for i in 0..6 {
            let d = coords[i].distance(coords[(i + 1) % 6]);
            assert!((d - 1.514).abs() < 0.03, "C-C = {d:.3}");
            sum += angle_deg(&coords, (i + 5) % 6, i, (i + 1) % 6);
        }
        let mean = sum / 6.0;
        // 椅子形の C-C-C ≈ 110°〜111.5°
        assert!((mean - 110.5).abs() < 2.5, "mean angle = {mean:.1}");
    }

    #[test]
    fn water_angle() {
        let (_, coords) = optimized("O");
        // H-O-H: UFF θ0 = 104.51°
        let a = angle_deg(&coords, 1, 0, 2);
        assert!((a - 104.51).abs() < 2.0, "H-O-H = {a:.1}");
    }

    /// C10: 実験トーションが二次アミドを trans (ω≈180°) に held する。
    #[test]
    fn amide_prefers_trans() {
        fn dihedral_deg(c: &[Vec3], i: usize, j: usize, k: usize, l: usize) -> f64 {
            let b1 = c[j] - c[i];
            let b2 = c[k] - c[j];
            let b3 = c[l] - c[k];
            let n1 = b1.cross(b2);
            let n2 = b2.cross(b3);
            let cos = (n1.dot(n2) / (n1.norm() * n2.norm())).clamp(-1.0, 1.0);
            cos.acos().to_degrees()
        }
        // N-メチルアセトアミド CC(=O)NC。trans アミドは
        // Cメチル-C(=O)-N-Cメチル ≈ 180° (ペプチドの ω)、O=C-N-C ≈ 0°
        let g = build_molecule_graph("CC(=O)NC").unwrap();
        for seed in [1u64, 7, 42] {
            let conf = embed_molecule(
                &g,
                &EmbedParams {
                    seed,
                    ..EmbedParams::default()
                },
            )
            .expect("embeds");
            let omega = dihedral_deg(&conf.coords, 0, 1, 3, 4); // C-C-N-C
            assert!(
                omega > 150.0,
                "seed {seed}: amide omega = {omega:.1} (want trans ~180)"
            );
            let oc = dihedral_deg(&conf.coords, 2, 1, 3, 4); // O=C-N-C
            assert!(oc < 30.0, "seed {seed}: O=C-N-C = {oc:.1} (want ~0)");
        }
    }

    #[test]
    fn stereo_survives_uff() {
        for smi in ["N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O", "C/C=C/C", "C/C=C\\C"] {
            let g = build_molecule_graph(smi).unwrap();
            let conf = embed_molecule(&g, &EmbedParams::default()).expect("embeds");
            assert!(
                crate::conformer::verify_stereo_3d(&g, &conf),
                "{smi}: stereo lost after UFF"
            );
        }
    }
}
