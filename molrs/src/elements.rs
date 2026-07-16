//! 周期表: 元素記号と原子番号の対応。

/// 原子番号順の元素記号 (index = Z - 1)。
pub const SYMBOLS: [&str; 118] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
    "Fl", "Mc", "Lv", "Ts", "Og",
];

/// 元素記号 (正しい大文字小文字) から原子番号を返す。未知なら None。
pub fn atomic_number(symbol: &str) -> Option<u8> {
    SYMBOLS
        .iter()
        .position(|&s| s == symbol)
        .map(|i| (i + 1) as u8)
}

/// SMILES で芳香族小文字表記が許される元素 (小文字表記) → 正規記号。
/// OpenSMILES: b, c, n, o, p, s, se, as (+ RDKit が受理する te)。
pub fn aromatic_symbol(lower: &str) -> Option<&'static str> {
    match lower {
        "b" => Some("B"),
        "c" => Some("C"),
        "n" => Some("N"),
        "o" => Some("O"),
        "p" => Some("P"),
        "s" => Some("S"),
        "se" => Some("Se"),
        "as" => Some("As"),
        "te" => Some("Te"),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atomic_numbers() {
        assert_eq!(atomic_number("H"), Some(1));
        assert_eq!(atomic_number("C"), Some(6));
        assert_eq!(atomic_number("Cl"), Some(17));
        assert_eq!(atomic_number("Se"), Some(34));
        assert_eq!(atomic_number("Og"), Some(118));
        assert_eq!(atomic_number("Xx"), None);
        assert_eq!(atomic_number("c"), None); // 小文字は aromatic_symbol の担当
    }

    #[test]
    fn aromatic_symbols() {
        assert_eq!(aromatic_symbol("c"), Some("C"));
        assert_eq!(aromatic_symbol("se"), Some("Se"));
        assert_eq!(aromatic_symbol("te"), Some("Te"));
        assert_eq!(aromatic_symbol("cl"), None);
        assert_eq!(aromatic_symbol("C"), None);
    }
}
