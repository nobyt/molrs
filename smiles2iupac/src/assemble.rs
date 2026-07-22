//! IUPAC 名の組み立て (S3.5, Python name_assembler.py の簡約移植)。
//!
//! 主鎖 stem + 不飽和 (ene/yne) + suffix + 置換基接頭辞を結合する。
//! 立体記述子・環・複雑接頭辞は未対応。

use crate::constants::{CHAIN_PREFIX, MULTIPLIER};

const BOND_MULT: [&str; 6] = ["", "", "di", "tri", "tetra", "penta"];

/// 置換基のロカント種別。
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Loc {
    /// 主鎖上の番号ロカント (None = ロカント省略)
    Num(Option<usize>),
    /// アミン N 上のロカント (0 = "N", 1 = "N'", ...)
    N(usize),
}

impl Loc {
    fn as_str(&self) -> String {
        match self {
            Loc::Num(Some(n)) => n.to_string(),
            Loc::Num(None) => String::new(),
            Loc::N(k) => format!("N{}", "'".repeat(*k)),
        }
    }
    fn is_n(&self) -> bool {
        matches!(self, Loc::N(_))
    }
    /// ソートキー (数字 < N)。
    fn sort_key(&self) -> (u8, usize) {
        match self {
            Loc::Num(Some(n)) => (0, *n),
            Loc::Num(None) => (0, 0),
            Loc::N(k) => (1, *k),
        }
    }
}

/// 置換基 (ロカント, 名前) → 接頭辞文字列。アルファベット順・倍数詞集約。
/// C 置換基と N 置換基は同名でも別グループ (Python 準拠)。
pub fn build_prefix(substituents: &[(Loc, String)]) -> String {
    if substituents.is_empty() {
        return String::new();
    }
    // (名前, N か) でグループ化
    let mut groups: Vec<(String, bool, Vec<Loc>)> = Vec::new();
    for (loc, name) in substituents {
        let is_n = loc.is_n();
        if let Some(e) = groups.iter_mut().find(|(nm, n, _)| nm == name && *n == is_n) {
            e.2.push(*loc);
        } else {
            groups.push((name.clone(), is_n, vec![*loc]));
        }
    }
    // アルファベット順、同名は C→N
    groups.sort_by(|a, b| alpha_key(&a.0).cmp(alpha_key(&b.0)).then(a.1.cmp(&b.1)));

    let all_loc_none = substituents
        .iter()
        .all(|(l, _)| matches!(l, Loc::Num(None)));

    let mut parts: Vec<String> = Vec::new();
    for (name, _is_n, locs) in &groups {
        let n = locs.len();
        let mut sorted = locs.clone();
        sorted.sort_by_key(|l| l.sort_key());
        let loc_strs: Vec<String> = sorted.iter().map(|l| l.as_str()).filter(|s| !s.is_empty()).collect();
        let complex = _needs_parens(name);
        if loc_strs.is_empty() {
            // ロカント省略
            if n > 1 {
                let mult = MULTIPLIER.get(n).copied().unwrap_or("");
                if complex {
                    parts.push(format!("{mult}({name})"));
                } else {
                    parts.push(format!("{mult}{name}"));
                }
            } else if complex {
                parts.push(format!("({name})"));
            } else {
                parts.push(name.clone());
            }
        } else {
            let loc_str = loc_strs.join(",");
            let mult = if n > 1 {
                MULTIPLIER.get(n).copied().unwrap_or("")
            } else {
                ""
            };
            if complex {
                // 複合置換基は bis/tris + 括弧
                let m = if n > 1 { BIS_MULT.get(n).copied().unwrap_or(mult) } else { "" };
                parts.push(format!("{loc_str}-{m}({name})"));
            } else {
                parts.push(format!("{loc_str}-{mult}{name}"));
            }
        }
    }
    let sep = if all_loc_none { "" } else { "-" };
    parts.join(sep)
}

const BIS_MULT: [&str; 5] = ["", "", "bis", "tris", "tetrakis"];

/// 複合置換基 (ロカント/括弧/置換されたヘテロ置換基) は括弧が必要。
fn _needs_parens(name: &str) -> bool {
    name.contains('-')
        || name.chars().any(|c| c.is_ascii_digit())
        || name.starts_with('(')
        // 置換された sulfanyl/amino など (例 methylsulfanyl) は囲む
        || (name.ends_with("sulfanyl") && name != "sulfanyl")
        || (name.ends_with("amino") && name != "amino")
}

fn alpha_key(name: &str) -> &str {
    for mult in ["di", "tri", "tetra", "penta", "hexa", "hepta", "octa"] {
        if let Some(rest) = name.strip_prefix(mult) {
            // "decyl" 等の誤集約を避けるため、既知倍数詞のみ
            return rest;
        }
    }
    name
}

/// ene/yne ロカント → suffix 直前の中間文字列。
fn format_multiple_bonds(ene: &[usize], yne: &[usize]) -> String {
    let needs_a = ene.len() > 1 || yne.len() > 1;
    let mut out = String::new();
    if !ene.is_empty() {
        let loc = ene.iter().map(|l| l.to_string()).collect::<Vec<_>>().join(",");
        let mult = BOND_MULT.get(ene.len()).copied().unwrap_or("");
        out.push_str(&format!("-{loc}-{mult}en"));
    }
    if !yne.is_empty() {
        let loc = yne.iter().map(|l| l.to_string()).collect::<Vec<_>>().join(",");
        let mult = BOND_MULT.get(yne.len()).copied().unwrap_or("");
        out.push_str(&format!("-{loc}-{mult}yn"));
    }
    if needs_a && out.starts_with('-') {
        out = format!("a{out}");
    }
    out
}

fn loc_list(locs: &[usize]) -> String {
    let mut v = locs.to_vec();
    v.sort_unstable();
    v.iter().map(|l| l.to_string()).collect::<Vec<_>>().join(",")
}

/// 語幹 + 不飽和 + suffix の本体。対応外の suffix は None。
pub fn build_name_body(
    stem: &str,
    suffix: &str,
    ene: &[usize],
    yne: &[usize],
    suffix_locants: &[usize],
    chain_length: usize,
) -> Option<String> {
    let has_mb = !ene.is_empty() || !yne.is_empty();
    let mb = || format_multiple_bonds(ene, yne);
    let sl = || suffix_locants.first().copied().unwrap_or(1);

    let body = match suffix {
        "ane" => format!("{stem}ane"),
        "ene" => {
            if !yne.is_empty() {
                format!("{stem}{}e", mb())
            } else if chain_length == 2 {
                format!("{stem}ene")
            } else if ene.len() == 1 {
                format!("{stem}-{}-ene", loc_list(ene))
            } else {
                let mult = ["", "", "di", "tri", "tetra", "penta"]
                    .get(ene.len())
                    .copied()
                    .unwrap_or("");
                format!("{stem}a-{}-{mult}ene", loc_list(ene))
            }
        }
        "yne" => {
            if chain_length == 2 {
                format!("{stem}yne")
            } else if yne.len() == 1 {
                format!("{stem}-{}-yne", loc_list(yne))
            } else {
                let mult = BOND_MULT.get(yne.len()).copied().unwrap_or("");
                format!("{stem}a-{}-{mult}yne", loc_list(yne))
            }
        }
        "al" => {
            if has_mb {
                format!("{stem}{}al", mb())
            } else {
                format!("{stem}anal")
            }
        }
        "dial" => {
            if has_mb {
                format!("{stem}{}edial", mb())
            } else {
                format!("{stem}anedial")
            }
        }
        "oic acid" => {
            if has_mb {
                format!("{stem}{}oic acid", mb())
            } else if chain_length == 1 {
                "formic acid".to_string()
            } else if chain_length == 2 {
                "acetic acid".to_string()
            } else {
                format!("{stem}anoic acid")
            }
        }
        "dioic acid" => {
            if has_mb {
                format!("{stem}{}edioic acid", mb())
            } else {
                format!("{stem}anedioic acid")
            }
        }
        "ol" => {
            let loc = sl();
            if has_mb {
                if chain_length == 2 && !ene.is_empty() && yne.is_empty() {
                    format!("{stem}enol")
                } else if chain_length == 2 && !yne.is_empty() && ene.is_empty() {
                    format!("{stem}ynol")
                } else {
                    format!("{stem}{}-{loc}-ol", mb())
                }
            } else if chain_length <= 2 && loc == 1 {
                format!("{stem}anol")
            } else {
                format!("{stem}an-{loc}-ol")
            }
        }
        "diol" | "triol" | "tetraol" => {
            let mult = suffix.strip_suffix("ol").unwrap_or("di");
            if has_mb {
                return None;
            }
            if chain_length == 1 {
                format!("{stem}ane{mult}ol")
            } else {
                format!("{stem}ane-{}-{mult}ol", loc_list(suffix_locants))
            }
        }
        "one" => {
            if chain_length == 1 {
                format!("{stem}anone")
            } else {
                let loc = suffix_locants.first().copied().unwrap_or(2);
                if has_mb {
                    if chain_length == 2 && ene == [1] && yne.is_empty() {
                        if loc == 1 {
                            format!("{stem}enone")
                        } else {
                            format!("{stem}en-{loc}-one")
                        }
                    } else {
                        format!("{stem}{}-{loc}-one", mb())
                    }
                } else {
                    format!("{stem}an-{loc}-one")
                }
            }
        }
        "dione" | "trione" | "tetraone" => {
            let mult = suffix.strip_suffix("one").unwrap_or("di");
            if has_mb {
                return None;
            }
            format!("{stem}ane-{}-{mult}one", loc_list(suffix_locants))
        }
        "amine" => {
            let loc = sl();
            if has_mb {
                if chain_length == 2 && ene == [1] && yne.is_empty() && loc == 1 {
                    format!("{stem}enamine")
                } else {
                    format!("{stem}{}-{loc}-amine", mb())
                }
            } else if chain_length <= 2 && loc == 1 {
                format!("{stem}anamine")
            } else {
                format!("{stem}an-{loc}-amine")
            }
        }
        "diamine" | "triamine" => {
            let mult = suffix.strip_suffix("amine").unwrap_or("di");
            if has_mb {
                return None;
            }
            if chain_length == 1 {
                format!("{stem}ane{mult}amine")
            } else {
                format!("{stem}ane-{}-{mult}amine", loc_list(suffix_locants))
            }
        }
        "nitrile" => {
            if has_mb {
                format!("{stem}{}enitrile", mb())
            } else {
                format!("{stem}anenitrile")
            }
        }
        "dinitrile" => {
            if has_mb {
                return None;
            }
            format!("{stem}anedinitrile")
        }
        _ => return None,
    };
    Some(body)
}

/// 主鎖長・主基・多重結合・置換基・suffix ロカント・立体記述子から
/// 完全な IUPAC 名を組み立てる。stereo = (ロカント, 記述子文字列) の列。
#[allow(clippy::too_many_arguments)]
pub fn assemble_name(
    chain_length: usize,
    principal_group_type: &str,
    ene: &[usize],
    yne: &[usize],
    substituents: &[(usize, String)],
    n_substituents: &[(usize, String)],
    suffix: &str,
    suffix_locants: &[usize],
    stereo: &[(usize, String)],
) -> Option<String> {
    let stem = CHAIN_PREFIX.get(chain_length).copied()?;

    // ロカント省略ルール (簡約): 1 炭素鎖、または特定の 2 炭素鎖単一置換基
    let drop_locant = chain_length == 1
        || (chain_length == 2
            && substituents.len() == 1
            && n_substituents.is_empty()
            && matches!(principal_group_type, "alkane" | "alkene"))
        || (chain_length == 2
            && principal_group_type == "alkane"
            && substituents.len() == 6
            && substituents.iter().map(|(_, n)| n).collect::<std::collections::HashSet<_>>().len()
                == 1);
    let mut eff: Vec<(Loc, String)> = substituents
        .iter()
        .map(|(l, n)| {
            let loc = if drop_locant { Loc::Num(None) } else { Loc::Num(Some(*l)) };
            (loc, n.clone())
        })
        .collect();
    for (ni, name) in n_substituents {
        eff.push((Loc::N(*ni), name.clone()));
    }

    let prefix = build_prefix(&eff);
    let body = build_name_body(stem, suffix, ene, yne, suffix_locants, chain_length)?;
    let result = format!("{prefix}{body}");

    if stereo.is_empty() {
        return Some(result);
    }
    // 立体記述子: ロカント昇順、(2E,3R) 形式で前置。
    // 単一の R/S 中心はロカントを省略する ((R)-…)。
    let mut st = stereo.to_vec();
    st.sort_by_key(|(l, _)| *l);
    let combined = if st.len() == 1 && matches!(st[0].1.as_str(), "R" | "S") && st[0].0 == 1 {
        st[0].1.clone()
    } else {
        st.iter()
            .map(|(l, d)| format!("{l}{d}"))
            .collect::<Vec<_>>()
            .join(",")
    };
    Some(format!("({combined})-{result}"))
}
