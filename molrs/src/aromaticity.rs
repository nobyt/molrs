//! ケクレ化と芳香族認識 (S1.3)。RDKit デフォルトモデル互換。
//!
//! RDKit の MolFromSmiles サニタイズと同じ 2 段階:
//! 1. **ケクレ化**: 芳香族表記の原子のうち「二重結合を 1 本必要とする」原子を
//!    環内芳香族結合上の完全マッチングで対にする (できなければ不正 SMILES)。
//! 2. **芳香族認識**: ケクレ構造に対し電子供与タイプを判定し、SSSR 環と
//!    縮合環ユニオンに Hückel 則 (π 電子数 ≡ 2 mod 4) を適用して
//!    原子・結合の芳香族フラグを付け直す。
//!
//! これにより `O=c1ccc(=O)cc1` (キノン) は芳香族表記でも非芳香族に、
//! `C1=CC=CC=C1` はケクレ表記でも芳香族になる (RDKit と同じ挙動)。

use crate::ChemError;

/// ケクレ化・認識対象の原子メタデータ (グラフ構築中の中間表現)。
pub(crate) struct AromAtom<'a> {
    pub symbol: &'a str,
    pub charge: i8,
    /// SMILES 上で芳香族表記だったか
    pub input_aromatic: bool,
    /// 総 H 数 (暗黙 + 明示 + マージ)
    pub num_hs: u8,
}

/// 結合の中間表現。order はケクレ化前は芳香族候補 = 1.0。
pub(crate) struct AromBond {
    pub a: usize,
    pub b: usize,
    pub order: f64,
    /// SMILES 上の芳香族候補 (小文字原子間の省略結合 or `:` 結合)
    pub aromatic_candidate: bool,
    /// 環内結合か (非橋エッジ)
    pub in_ring: bool,
}

/// 電荷調整済みのデフォルト原子価 (ケクレ化の「二重結合が必要か」判定用)。
/// `edit` モジュール (I45) も帯電原子の暗黙 H 再計算に再利用する。
pub(crate) fn adjusted_valence(symbol: &str, charge: i8) -> Option<i32> {
    let v = match symbol {
        "C" => 4 - charge.abs() as i32,
        "N" | "P" | "As" => 3 + charge as i32,
        "O" | "S" | "Se" | "Te" => 2 + charge as i32,
        "B" => 3 - charge as i32,
        _ => return None,
    };
    Some(v)
}

/// 最外殻電子数 (RDKit の getNouterElecs 相当; 芳香族関連元素のみ)。
fn outer_electrons(symbol: &str) -> Option<i32> {
    Some(match symbol {
        "B" => 3,
        "C" | "Si" => 4,
        "N" | "P" | "As" => 5,
        "O" | "S" | "Se" | "Te" => 6,
        "F" | "Cl" | "Br" | "I" => 7,
        "H" => 1,
        _ => return None,
    })
}

/// 芳香族になり得る元素 (RDKit デフォルトモデル)。
fn arom_allowed(symbol: &str) -> bool {
    matches!(
        symbol,
        "C" | "N" | "O" | "S" | "Se" | "Te" | "P" | "As" | "B"
    )
}

/// ケクレ化: 芳香族表記原子に環内二重結合を割り当てる。
/// `bonds` の order を書き換える (選ばれた結合 → 2.0)。
pub(crate) fn kekulize(atoms: &[AromAtom<'_>], bonds: &mut [AromBond]) -> Result<(), ChemError> {
    let n = atoms.len();

    // 芳香族原子が環内にあるかチェック
    let mut in_ring_atom = vec![false; n];
    for b in bonds.iter() {
        if b.in_ring {
            in_ring_atom[b.a] = true;
            in_ring_atom[b.b] = true;
        }
    }
    for (i, a) in atoms.iter().enumerate() {
        if a.input_aromatic && !in_ring_atom[i] {
            return Err(ChemError::InvalidSmiles(format!(
                "aromatic atom {i} is not in a ring"
            )));
        }
    }

    // 二重結合を必要とする原子 (needs == 1 以上)
    let mut needy = vec![false; n];
    for (i, a) in atoms.iter().enumerate() {
        if !a.input_aromatic {
            continue;
        }
        let Some(adj_val) = adjusted_valence(a.symbol, a.charge) else {
            continue;
        };
        let mut used = a.num_hs as i32;
        for b in bonds.iter() {
            if b.a == i || b.b == i {
                used += b.order as i32; // 芳香族候補は 1
            }
        }
        needy[i] = adj_val - used >= 1;
    }

    // マッチング対象エッジ: 環内の芳香族候補結合で両端が needy
    let cand_edges: Vec<usize> = (0..bonds.len())
        .filter(|&ei| {
            let b = &bonds[ei];
            b.aromatic_candidate && b.in_ring && needy[b.a] && needy[b.b]
        })
        .collect();
    let mut cand_adj: Vec<Vec<usize>> = vec![Vec::new(); n]; // atom -> cand edge idx
    for &ei in &cand_edges {
        cand_adj[bonds[ei].a].push(ei);
        cand_adj[bonds[ei].b].push(ei);
    }

    let needy_atoms: Vec<usize> = (0..n).filter(|&i| needy[i]).collect();
    let mut matched = vec![false; n];
    let mut chosen: Vec<usize> = Vec::new();
    if !backtrack_matching(&needy_atoms, &cand_adj, bonds, &mut matched, &mut chosen) {
        return Err(ChemError::InvalidSmiles(
            "cannot kekulize aromatic system".into(),
        ));
    }
    for ei in chosen {
        bonds[ei].order = 2.0;
    }
    Ok(())
}

/// needy 原子全部を被覆する完全マッチングをバックトラッキングで探す。
fn backtrack_matching(
    needy: &[usize],
    cand_adj: &[Vec<usize>],
    bonds: &[AromBond],
    matched: &mut [bool],
    chosen: &mut Vec<usize>,
) -> bool {
    // 未マッチで候補数最小の原子を選ぶ (枝刈り)
    let mut best: Option<(usize, usize)> = None; // (候補数, atom)
    for &u in needy {
        if matched[u] {
            continue;
        }
        let cnt = cand_adj[u]
            .iter()
            .filter(|&&ei| {
                let b = &bonds[ei];
                !matched[b.a] && !matched[b.b]
            })
            .count();
        if cnt == 0 {
            return false; // この原子はもうマッチできない
        }
        if best.is_none() || cnt < best.unwrap().0 {
            best = Some((cnt, u));
        }
    }
    let Some((_, u)) = best else {
        return true; // 全員マッチ済み
    };
    for &ei in &cand_adj[u] {
        let (a, b) = (bonds[ei].a, bonds[ei].b);
        if matched[a] || matched[b] {
            continue;
        }
        matched[a] = true;
        matched[b] = true;
        chosen.push(ei);
        if backtrack_matching(needy, cand_adj, bonds, matched, chosen) {
            return true;
        }
        chosen.pop();
        matched[a] = false;
        matched[b] = false;
    }
    false
}

/// 電子供与タイプ: 原子が環系に供出する π 電子数。None = 芳香族になれない。
///
/// RDKit 実測に基づく規則:
/// - 環内二重結合 (どの環かは問わずグラフ上の環に属する結合) → 1 電子。
///   アントラセン中央環の原子は二重結合が外環側にあっても 1 電子と数える。
/// - 環に属さない二重結合 → 相手がより電気陰性 (最外殻電子数が多い) なら
///   0 電子 (キノンの C=O)、そうでなければ候補外 (フルベンの C=C)。
/// - 多重結合なし → 孤立電子対から固定値 (ピロール N は常に 2 電子。
///   これにより 2H-キノリジンの sp2 環は 7 電子で不合格になる)。
fn electron_contribution(atoms: &[AromAtom<'_>], bonds: &[AromBond], i: usize) -> Option<u32> {
    let a = &atoms[i];
    if !arom_allowed(a.symbol) {
        return None;
    }
    // sp2 条件: 重原子次数 + H ≤ 3
    let heavy_deg = bonds.iter().filter(|b| b.a == i || b.b == i).count();
    if heavy_deg + a.num_hs as usize > 3 {
        return None;
    }

    // 多重結合の分類
    let mut ring_multiple = 0usize;
    let mut exo_partner: Option<usize> = None;
    let mut n_multiple = 0usize;
    for b in bonds.iter() {
        if b.a != i && b.b != i {
            continue;
        }
        if b.order >= 2.0 {
            n_multiple += 1;
            if b.order >= 3.0 {
                return None; // 三重結合
            }
            if b.in_ring {
                ring_multiple += 1;
            } else {
                exo_partner = Some(if b.a == i { b.b } else { b.a });
            }
        }
    }
    if n_multiple > 1 {
        return None; // 累積二重結合 (例: 芳香族表記のスルホン S)
    }
    if ring_multiple == 1 {
        return Some(1);
    }
    if let Some(p) = exo_partner {
        let oe_i = outer_electrons(a.symbol)?;
        let oe_p = outer_electrons(atoms[p].symbol)?;
        return if oe_p > oe_i { Some(0) } else { None };
    }
    // 多重結合なし: 孤立電子対 (RDKit countAtomElec 相当)
    let oe = outer_electrons(a.symbol)?;
    let nelec = oe - a.charge as i32 - (heavy_deg as i32 + a.num_hs as i32);
    match nelec {
        i32::MIN..=-1 => None,
        0 => Some(0), // 空軌道 (B, C+ など)
        1 => Some(1),
        _ => Some(2),
    }
}

/// 芳香族認識。返り値: (原子ごとの芳香族フラグ, 結合ごとの芳香族フラグ)。
/// `rings` は SSSR (原子列)。
pub(crate) fn perceive_aromaticity(
    atoms: &[AromAtom<'_>],
    bonds: &[AromBond],
    rings: &[Vec<usize>],
) -> (Vec<bool>, Vec<bool>) {
    let n = atoms.len();
    let mut atom_arom = vec![false; n];
    let mut bond_arom = vec![false; bonds.len()];
    if rings.is_empty() {
        return (atom_arom, bond_arom);
    }

    // 結合インデックスの逆引き
    let bond_index = |u: usize, v: usize| -> Option<usize> {
        bonds
            .iter()
            .position(|b| (b.a == u && b.b == v) || (b.a == v && b.b == u))
    };

    // 環ごとの結合集合 (SSSR 環のうち結合が引けるもの)
    struct RingInfo {
        atom_set: Vec<usize>,
        bond_set: Vec<usize>,
    }
    let mut ring_infos: Vec<RingInfo> = Vec::new();
    for ring in rings {
        let mut bond_set = Vec::new();
        let mut ok = true;
        for k in 0..ring.len() {
            let (u, v) = (ring[k], ring[(k + 1) % ring.len()]);
            match bond_index(u, v) {
                Some(bi) => bond_set.push(bi),
                None => {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            ring_infos.push(RingInfo {
                atom_set: ring.clone(),
                bond_set,
            });
        }
    }
    if ring_infos.is_empty() {
        return (atom_arom, bond_arom);
    }

    let contrib: Vec<Option<u32>> = (0..n)
        .map(|i| electron_contribution(atoms, bonds, i))
        .collect();

    // Hückel 判定: 全原子が候補で π 電子和 ≡ 2 mod 4。
    //
    // ただし π 電子総和が最小値の 2 で通るのは 3 員環 (シクロプロペニル
    // カチオン型) だけに限る (I49)。4 員環以上で電子数 2 は「環内の実二重
    // 結合 1 本だけがあり、残りの環メンバーは全員が環外二重結合で 0 電子
    // 寄与」というだけの局所的な状況を意味し、環全体で電子が非局在化して
    // いるわけではない — トロポン型 (7 員環、環外 C=O が 0 電子、残り 6 個の
    // 環メンバーが環内で完全に交互する 3 本の二重結合を持ち計 6 電子) の
    // ような正当な芳香族性とは違う。例: シクロブテン-1,2-ジオン環に
    // アミノ・ヒドロキシ置換基が付いた分子 (`CCNC1=C(O)C(=O)C1=O`、
    // スクエア酸誘導体) は環内二重結合が 1 本だけで残り 2 個の環員は環外
    // C=O が 0 電子寄与するだけなので電子数 2 で誤って芳香族と判定され、
    // 可動 H 検出のブリッジ探索が誤発火して本来固定のはずの環外 NH/OH を
    // 環内カルボニルと誤併合していた (実 InChI はこの環を非芳香族・非可動
    // として扱う)。
    let huckel = |atom_set: &[usize]| -> bool {
        let mut total = 0u32;
        for &i in atom_set {
            match contrib[i] {
                Some(c) => total += c,
                None => return false,
            }
        }
        if total == 2 {
            return atom_set.len() <= 3;
        }
        total >= 2 && (total - 2).is_multiple_of(4)
    };

    // 単環: 合格した環は全結合をマーク (ナフタレンの縮合結合はここで芳香族になる)
    let nr = ring_infos.len();
    for r in &ring_infos {
        if huckel(&r.atom_set) {
            for &i in &r.atom_set {
                atom_arom[i] = true;
            }
            for &bi in &r.bond_set {
                bond_arom[bi] = true;
            }
        }
    }

    // 縮合環ペアのユニオン (RDKit 実測に基づく規則):
    // - ペア (2 環) のみ試す。3 環以上のユニオンは行わない
    //   (ペリミジンのアミジン環がナフタレン全体と合算されて 14 電子で
    //    合格してしまうのを防ぐ。アズレン・インドリジンはペアで足りる)
    // - 合格したユニオンは原子と周縁結合 (ちょうど 1 つのメンバー環に属する
    //   結合) をマークする。共有結合は単結合のまま (アズレンの縮合結合)
    if nr >= 2 {
        for i in 0..nr {
            for j in i + 1..nr {
                let shares_bond = ring_infos[i]
                    .bond_set
                    .iter()
                    .any(|b| ring_infos[j].bond_set.contains(b));
                if !shares_bond {
                    continue;
                }
                let mut atom_set: Vec<usize> = ring_infos[i]
                    .atom_set
                    .iter()
                    .chain(ring_infos[j].atom_set.iter())
                    .copied()
                    .collect();
                atom_set.sort_unstable();
                atom_set.dedup();
                if !huckel(&atom_set) {
                    continue;
                }
                for &a in &atom_set {
                    atom_arom[a] = true;
                }
                // 周縁結合 = どちらか一方の環だけに属する結合
                for &bi in &ring_infos[i].bond_set {
                    if !ring_infos[j].bond_set.contains(&bi) {
                        bond_arom[bi] = true;
                    }
                }
                for &bi in &ring_infos[j].bond_set {
                    if !ring_infos[i].bond_set.contains(&bi) {
                        bond_arom[bi] = true;
                    }
                }
            }
        }
    }

    (atom_arom, bond_arom)
}
