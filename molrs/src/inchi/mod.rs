//! SMILES/分子グラフ → 標準 InChI (`InChI=1S/…`) と InChIKey。
//!
//! IUPAC 公式 InChI (C 実装) とビット完全一致を目標とする。計画:
//! RUST_INCHI_PLAN.md。
//!
//! 対象: 骨格層 (式・接続 `c`・水素 `h`・電荷 `q`/`p`)、立体 (`b`/`t`/`m`/`s`)、
//! 同位体 (`i`)、InChIKey、可動電荷 (`(H3-,…)` の荷電可動 H 群)。
//!
//! InChIKey は標準 InChI 文字列の SHA-256 ハッシュ (base-26 符号化) なので、
//! 依存クレートゼロを保つため SHA-256 を自前実装している ([`sha256`])。

pub mod base26;
pub(crate) mod blossom;
pub(crate) mod disconnect;
pub(crate) mod formula;
pub(crate) mod layers;
pub(crate) mod normalize;
pub mod number;
pub mod sha256;
pub(crate) mod stereo;

pub use base26::inchi_key_from_string;

use crate::graph::MoleculeGraph;

/// InChI 生成のエラー。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InchiError {
    /// 不正な SMILES ([`inchi_of`] 経由のみ)
    InvalidSmiles(String),
    /// v1 で未対応の構造クラス (実装の進行に応じて縮小)
    Unsupported(String),
}

impl std::fmt::Display for InchiError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            InchiError::InvalidSmiles(s) => write!(f, "invalid SMILES: {s}"),
            InchiError::Unsupported(s) => write!(f, "InChI unsupported: {s}"),
        }
    }
}

impl std::error::Error for InchiError {}

/// 分子グラフの Hill 式層を返す (I2)。式のみが必要な場合の公開 API。
pub fn formula(g: &MoleculeGraph) -> String {
    formula::formula_layer(g)
}

/// 成分ごとの層文字列を `;` で連結する (I20)。
///
/// 連続する同一の非空文字列は `N*` で圧縮する (`c2*1-2;` = 同一成分 2 つ +
/// 空の成分 1 つ)。空文字列は圧縮せず `;` を並べる。全成分が空なら空文字列を
/// 返し、呼び出し側が層自体を省略する (`InChI=1S/Bi.3H` に c/h 層がないのは
/// このため)。
fn join_components(parts: &[String]) -> String {
    join_components_keyed(parts, None)
}

/// [`join_components`] の本体。`keys` が `Some` なら、文字列が一致しても
/// `keys[i] != keys[j]` の隣接成分は圧縮しない (I46)。
///
/// 実 InChI (`ichiprt3.c`) の `/h` 層 (`str_H_atoms`) は圧縮前に
/// `pINChI_Prev->nNumberOfAtoms == pINChI->nNumberOfAtoms` を要求してから
/// `nNum_H` 配列を `memcmp` する — **レンダリング後の文字列同士**ではなく
/// **成分の重原子数 + 生の H 数配列**で同一性を判定している。molrs は
/// レンダリング後の文字列だけを比較していたため、式も重原子数も違う
/// 2 成分 (末端 CH3 が 2 個・内部 CH2 が 2 個という部分だけが偶然同じ文字列
/// になる、例: ペンタン-3-オン `CCC(=O)CC` とブタン `CCCC`) を誤って `2*…`
/// に圧縮してしまっていた。`/c` (`str_Connections`) は `lenConnTable` の
/// 一致 (= 重原子数と等価) を要求してから接続表を `memcmp` するので、
/// 実質的に同じ制約を最初から満たしており today's fix の対象外
/// (エタン+アセチレンのように式が違っても接続表がたまたま一致すれば
/// 圧縮されるのは正しい実 InChI の挙動、詳細は RUST_INCHI_I29_PLAN.md I46 節)。
fn join_components_keyed(parts: &[String], keys: Option<&[usize]>) -> String {
    if parts.iter().all(|s| s.is_empty()) {
        return String::new();
    }
    let mut out = String::new();
    let mut i = 0;
    while i < parts.len() {
        let mut j = i + 1;
        while j < parts.len()
            && parts[j] == parts[i]
            && keys.is_none_or(|k| k[j] == k[i])
        {
            j += 1;
        }
        let count = j - i;
        if i > 0 {
            out.push(';');
        }
        if parts[i].is_empty() {
            // 空成分の連続: 区切りだけを並べる
            for _ in 1..count {
                out.push(';');
            }
        } else {
            if count > 1 {
                out.push_str(&count.to_string());
                out.push('*');
            }
            out.push_str(&parts[i]);
        }
        i = j;
    }
    out
}

/// 標準 InChI (`InChI=1S/…`) を生成する。
///
/// 電荷は q/p 層で中性化して扱う。金属結合は標準 InChI の規約どおり切断し
/// (I20、[`disconnect`])、多成分は `;` 区切りで直列化する。多中心の環互変異性の
/// 一部は未対応。
pub fn to_inchi(g: &MoleculeGraph) -> Result<String, InchiError> {
    // 金属結合の切断 (標準 InChI は disconnected-metal 表現) → 電荷正規化
    let disconnected = disconnect::disconnect_metals(g);
    // 電荷正規化 (負の酸点を中性化・陽イオンを脱プロトン、残余 → q、移動 → p)
    let (ng, _q, p) =
        normalize::neutralize(&disconnected.graph, &disconnected.metal_locked_ligands);
    let g = &ng;

    let comps = layers::build_components(g);
    // 重原子を含まない H だけの成分 (常に末尾)。c/q/立体層には何も寄与せず、
    // h 層だけ水素分子 `[H][H]` が `1H` を出す。
    let h_sizes = disconnect::hydrogen_component_sizes(g);
    let pad = |mut v: Vec<String>| {
        v.resize(v.len() + h_sizes.len(), String::new());
        v
    };
    let pad_h = |mut v: Vec<String>| {
        v.extend(
            h_sizes
                .iter()
                .map(|&k| disconnect::hydrogen_component_h_layer(k)),
        );
        v
    };

    let formula = formula::formula_layer(g);
    let c = join_components(&pad(comps.iter().map(layers::connection_layer).collect()));
    // /h の N* 圧縮は文字列一致だけでなく成分の重原子数一致も要る (I46)。
    // 末尾の H のみ成分には実成分と衝突しない番兵として 0 を割り当てる
    // (どの実成分も重原子数 >= 1 を持つ)。
    let h_keys: Vec<usize> = comps
        .iter()
        .map(|c| c.inv.len())
        .chain(h_sizes.iter().map(|_| 0))
        .collect();
    let h = join_components_keyed(
        &pad_h(comps.iter().map(layers::hydrogen_layer).collect()),
        Some(&h_keys),
    );
    // /q は成分ごとの残余電荷 (中性成分は空欄)
    let q_parts: Vec<String> = comps
        .iter()
        .map(|comp| {
            let sum: i32 = comp
                .inv
                .iter()
                .map(|&a| g.atoms[a].formal_charge as i32)
                .sum();
            if sum == 0 {
                String::new()
            } else {
                format!("{sum:+}")
            }
        })
        .collect();
    let q = join_components(&pad(q_parts));
    let b = join_components(&pad(comps
        .iter()
        .map(|c| stereo::double_bond_layer(g, c))
        .collect()));
    let tms: Vec<(String, Option<char>, Option<char>)> = comps
        .iter()
        .map(|c| stereo::tetrahedral_layers(g, c))
        .collect();
    let t = join_components(&pad(tms.iter().map(|x| x.0.clone()).collect()));
    // /m は成分ごとの値を**区切り無しで連結**する (`;` でも、成分ごとに
    // `.` で区切るのでもない、I43)。立体中心の無い成分だけ `.` を
    // プレースホルダとして置く — 定義済みの値どうしは隣り合っても区切りを
    // 挟まない。例: 5 成分中 1〜5 番目が全て定義済みなら `/m01010`、末尾
    // 2 成分が未定義なら `/m11111..`、1 番目が未定義・2〜3 番目が定義済み
    // なら `/m.10` のようになる (以前は `join(".")` で定義済みどうしの間にも
    // 余計な `.` が挟まっていた)。
    let m_parts = pad(tms
        .iter()
        .map(|x| x.1.map(String::from).unwrap_or_default())
        .collect());
    let m = if m_parts.iter().all(|s| s.is_empty()) {
        String::new()
    } else {
        m_parts
            .iter()
            .map(|s| if s.is_empty() { "." } else { s.as_str() })
            .collect::<String>()
    };
    // /s は**構造全体で 1 個**であって成分ごとではない (コーパス 2,115 件すべて
    // `/s1`)。従来は成分ごとに並べて `/s1;` のようにしていた。
    let s_char = tms
        .iter()
        .find_map(|x| x.2)
        .map(String::from)
        .unwrap_or_default();

    let mut out = format!("InChI=1S/{formula}");
    if !c.is_empty() {
        out.push_str("/c");
        out.push_str(&c);
    }
    if !h.is_empty() {
        out.push_str("/h");
        out.push_str(&h);
    }
    // 電荷層 (h の後、立体の前): /q 残余電荷、/p プロトン化
    if !q.is_empty() {
        out.push_str("/q");
        out.push_str(&q);
    }
    if p != 0 {
        out.push_str(&format!("/p{p:+}"));
    }
    // 立体層 (順序: b, t, m, s)
    if !b.is_empty() {
        out.push_str("/b");
        out.push_str(&b);
    }
    if !t.is_empty() {
        out.push_str("/t");
        out.push_str(&t);
    }
    if !m.is_empty() {
        out.push_str("/m");
        out.push_str(&m);
    }
    if !s_char.is_empty() {
        out.push_str("/s");
        out.push_str(&s_char);
    }
    // 同位体層 (立体層の後、I37)
    let iso = join_components(&pad(comps
        .iter()
        .map(|c| layers::isotope_layer(g, c))
        .collect()));
    if !iso.is_empty() {
        out.push_str("/i");
        out.push_str(&iso);
    }
    Ok(out)
}

/// SMILES から標準 InChI を生成する便利関数。
pub fn inchi_of(smiles: &str) -> Result<String, InchiError> {
    let g = crate::graph::build_molecule_graph(smiles)
        .map_err(|e| InchiError::InvalidSmiles(e.to_string()))?;
    to_inchi(&g)
}

/// 分子グラフの InChIKey を生成する (I5、v1 範囲)。
pub fn to_inchi_key(g: &MoleculeGraph) -> Result<String, InchiError> {
    let inchi = to_inchi(g)?;
    Ok(inchi_key_from_string(&inchi))
}

/// SMILES から InChIKey を生成する便利関数。
pub fn inchi_key_of(smiles: &str) -> Result<String, InchiError> {
    let inchi = inchi_of(smiles)?;
    Ok(inchi_key_from_string(&inchi))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    #[test]
    fn formula_public_api() {
        let g = build_molecule_graph("CC(=O)O").unwrap();
        assert_eq!(formula(&g), "C2H4O2");
    }

    fn inchi(smiles: &str) -> String {
        inchi_of(smiles).unwrap()
    }

    /// 標準 InChI は金属結合を切断する (I20)。等方開裂 (金属-C/H) では電荷が
    /// 付かず、異方開裂 (金属-ハロゲン) では `/q`+`/p` が現れる。
    #[test]
    fn metal_bonds_are_disconnected() {
        assert_eq!(inchi("C[Hg]C"), "InChI=1S/2CH3.Hg/h2*1H3;");
        assert_eq!(inchi("C[Sn](C)(C)C"), "InChI=1S/4CH3.Sn/h4*1H3;");
        assert_eq!(inchi("CC[Hg]CC"), "InChI=1S/2C2H5.Hg/c2*1-2;/h2*1H2,2H3;");
        assert_eq!(inchi("C[Hg]Cl"), "InChI=1S/CH3.ClH.Hg/h1H3;1H;/q;;+1/p-1");
        assert_eq!(
            inchi("CC[Hg]Br"),
            "InChI=1S/C2H5.BrH.Hg/c1-2;;/h1H2,2H3;1H;/q;;+1/p-1"
        );
    }

    /// 金属水素化物は金属-H も切れ、H が独立成分になる (c/h 層は空欄のまま)。
    #[test]
    fn metal_hydrides_split_off_bare_hydrogens() {
        assert_eq!(inchi("[BiH3]"), "InChI=1S/Bi.3H");
        assert_eq!(inchi("[SnH4]"), "InChI=1S/Sn.4H");
        assert_eq!(inchi("C[PbH3]"), "InChI=1S/CH3.Pb.3H/h1H3;;;;");
        assert_eq!(
            inchi("CC[SnH2]CC"),
            "InChI=1S/2C2H5.Sn.2H/c2*1-2;;;/h2*1H2,2H3;;;"
        );
    }

    /// 水素分子は「骨格原子 1 個 + 結合 H 1 個」の単一成分 (孤立 H とは別扱い)。
    #[test]
    fn dihydrogen_is_a_single_component() {
        assert_eq!(inchi("[H][H]"), "InChI=1S/H2/h1H");
    }

    /// 塩は成分ごとに層が `;` で区切られ、`/q` は成分ごと・`/p` は全体で 1 つ。
    #[test]
    fn salts_serialize_per_component() {
        assert_eq!(
            inchi("[Na+].CC(=O)[O-]"),
            "InChI=1S/C2H4O2.Na/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1"
        );
        assert_eq!(
            inchi("[K+].[K+].[O-]C(=O)CCC(=O)[O-]"),
            "InChI=1S/C4H6O4.2K/c5-3(6)1-2-4(7)8;;/h1-2H2,(H,5,6)(H,7,8);;/q;2*+1/p-2"
        );
        assert_eq!(
            inchi("[NH4+].CC(=O)[O-]"),
            "InChI=1S/C2H4O2.H3N/c1-2(3)4;/h1H3,(H,3,4);1H3"
        );
    }

    /// 成分順序は「炭素数降順 → 重原子数降順 → H 数降順」(I29)。
    #[test]
    fn component_order_is_by_descending_carbon_then_size() {
        // 炭素数が主キー: C6H5 (6) が C2H4O2 (2) より先。重原子数は 6 < 4 では
        // ないので、旧規則 (重原子数昇順) では逆順になっていた。
        assert_eq!(
            inchi("CC(=O)O.C1=CC=C(C=C1)[Hg]"),
            "InChI=1S/C6H5.C2H4O2.Hg/c1-2-4-6-5-3-1;1-2(3)4;/h1-5H;1H3,(H,3,4);"
        );
        // 炭素数が並べば重原子数降順: C9H12O (10) → CO (2) → Fe (1)。
        assert!(
            inchi("C/C=C/C=C/C=C/C(=O)OC.[C-]#[O+].[C-]#[O+].[C-]#[O+].[Fe]")
                .starts_with("InChI=1S/C9H12O2.3CO.Fe/")
        );
    }

    /// **既知の残差** (I29): 「単原子カチオン + 多原子アニオン」の無機塩だけは
    /// 実 InChI がカチオンを先に置く。硫酸ナトリウムの実 InChI は
    /// `2Na.H2O4S` だが、炭素数降順→重原子数降順の規則では `H2O4S.2Na` になる。
    ///
    /// 単原子金属を先頭に特別扱いする規則も試したが、`FH.O3Si.2Zn` では
    /// 単原子の Zn が最後に来るため説明できず (電荷層との関係も含めて未解明)、
    /// PubChem 実データ 863 件中 33 件 (3.8%) がこの系統で残る。旧規則
    /// (重原子数昇順) はこの 1 件に合わせていたが、それは本リポジトリの
    /// 多成分 32 例への過適合で、PubChem では 33.7% しか再現できなかった。
    #[test]
    fn known_divergence_monoatomic_cation_salts() {
        assert_eq!(
            inchi("[Na+].[Na+].[O-]S(=O)(=O)[O-]"),
            // 実 InChI: "InChI=1S/2Na.H2O4S/c;;1-5(2,3)4/h;;(H2,1,2,3,4)/q2*+1;/p-2"
            "InChI=1S/H2O4S.2Na/c1-5(2,3)4;;/h(H2,1,2,3,4);;/q;2*+1/p-2"
        );
    }

    /// 非金属 (B/Si/As/Te) は切断されない — 金属表の境界の回帰テスト。
    #[test]
    fn metalloids_stay_connected() {
        assert_eq!(inchi("CC[AsH2]"), "InChI=1S/C2H7As/c1-2-3/h2-3H2,1H3");
        assert!(!inchi("Br[SiH2]C").contains('.'));
        assert!(!inchi("CB(O)O").contains('.'));
    }
}
