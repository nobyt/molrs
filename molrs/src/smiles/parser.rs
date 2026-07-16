//! SMILES パーサ本体。

use std::collections::HashMap;

use super::{AtomSpec, BondKind, BondSpec, Chirality, ParsedMolecule, RingClosure};
use crate::elements::{aromatic_symbol, atomic_number};
use crate::ChemError;

/// 環閉じ番号の未解決エントリ。
struct RingOpen {
    atom: usize,
    kind: BondKind,
    /// neighbor_order[atom] 内で環閉じ数字が現れた位置 (プレースホルダ)
    slot: usize,
}

struct Parser<'a> {
    input: &'a str,
    bytes: &'a [u8],
    pos: usize,
    atoms: Vec<AtomSpec>,
    bonds: Vec<BondSpec>,
    neighbor_order: Vec<Vec<usize>>,
    prev: Option<usize>,
    pending_bond: Option<BondKind>,
    branch_stack: Vec<usize>,
    ring_map: HashMap<u16, RingOpen>,
    /// 直前が `.` で、次に原子が必要
    after_dot: bool,
}

/// SMILES 文字列をパースする。
pub fn parse_smiles(smiles: &str) -> Result<ParsedMolecule, ChemError> {
    Parser::new(smiles).run()
}

impl<'a> Parser<'a> {
    fn new(input: &'a str) -> Self {
        Parser {
            input,
            bytes: input.as_bytes(),
            pos: 0,
            atoms: Vec::new(),
            bonds: Vec::new(),
            neighbor_order: Vec::new(),
            prev: None,
            pending_bond: None,
            branch_stack: Vec::new(),
            ring_map: HashMap::new(),
            after_dot: false,
        }
    }

    fn err(&self, msg: &str) -> ChemError {
        ChemError::InvalidSmiles(format!(
            "{msg} at position {} in {:?}",
            self.pos, self.input
        ))
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    fn bump(&mut self) -> Option<u8> {
        let c = self.peek();
        if c.is_some() {
            self.pos += 1;
        }
        c
    }

    fn run(mut self) -> Result<ParsedMolecule, ChemError> {
        if self.input.is_empty() {
            return Err(ChemError::InvalidSmiles("empty SMILES".into()));
        }
        while let Some(c) = self.peek() {
            match c {
                b'(' => {
                    self.pos += 1;
                    if self.pending_bond.is_some() {
                        return Err(self.err("bond symbol before '('"));
                    }
                    let p = self
                        .prev
                        .ok_or_else(|| self.err("branch without preceding atom"))?;
                    self.branch_stack.push(p);
                }
                b')' => {
                    self.pos += 1;
                    if self.pending_bond.is_some() {
                        return Err(self.err("dangling bond before ')'"));
                    }
                    let p = self
                        .branch_stack
                        .pop()
                        .ok_or_else(|| self.err("unmatched ')'"))?;
                    self.prev = Some(p);
                }
                b'.' => {
                    self.pos += 1;
                    if self.pending_bond.is_some() {
                        return Err(self.err("bond symbol before '.'"));
                    }
                    if self.atoms.is_empty() || self.after_dot {
                        return Err(self.err("misplaced '.'"));
                    }
                    self.prev = None;
                    self.after_dot = true;
                }
                b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\' => {
                    self.pos += 1;
                    if self.pending_bond.is_some() {
                        return Err(self.err("duplicate bond symbol"));
                    }
                    self.pending_bond = Some(match c {
                        b'-' => BondKind::Single,
                        b'=' => BondKind::Double,
                        b'#' => BondKind::Triple,
                        b'$' => BondKind::Quadruple,
                        b':' => BondKind::Aromatic,
                        b'/' => BondKind::Up,
                        _ => BondKind::Down,
                    });
                }
                b'0'..=b'9' => {
                    self.pos += 1;
                    let num = (c - b'0') as u16;
                    self.ring_bond(num)?;
                }
                b'%' => {
                    self.pos += 1;
                    let d1 = self.bump().filter(u8::is_ascii_digit);
                    let d2 = self.bump().filter(u8::is_ascii_digit);
                    let (Some(d1), Some(d2)) = (d1, d2) else {
                        return Err(self.err("'%' must be followed by two digits"));
                    };
                    let num = ((d1 - b'0') as u16) * 10 + (d2 - b'0') as u16;
                    self.ring_bond(num)?;
                }
                b'[' => {
                    self.pos += 1;
                    let atom = self.parse_bracket_atom()?;
                    self.add_atom(atom)?;
                }
                b'*' => {
                    self.pos += 1;
                    self.add_atom(AtomSpec {
                        symbol: "*".into(),
                        aromatic: false,
                        isotope: None,
                        charge: 0,
                        explicit_h: None,
                        chirality: None,
                        atom_class: None,
                        bracket: false,
                    })?;
                }
                _ => {
                    let atom = self.parse_organic_atom()?;
                    self.add_atom(atom)?;
                }
            }
        }
        if self.pending_bond.is_some() {
            return Err(self.err("dangling bond at end of SMILES"));
        }
        if !self.branch_stack.is_empty() {
            return Err(self.err("unclosed '('"));
        }
        if !self.ring_map.is_empty() {
            let mut nums: Vec<u16> = self.ring_map.keys().copied().collect();
            nums.sort_unstable();
            return Err(self.err(&format!("unclosed ring bond(s): {nums:?}")));
        }
        if self.after_dot {
            return Err(self.err("trailing '.'"));
        }
        debug_assert!(self
            .neighbor_order
            .iter()
            .flatten()
            .all(|&b| b != usize::MAX));
        Ok(ParsedMolecule {
            atoms: self.atoms,
            bonds: self.bonds,
            neighbor_order: self.neighbor_order,
        })
    }

    /// 有機サブセット原子 (角括弧なし): B, C, N, O, P, S, F, Cl, Br, I / b, c, n, o, p, s
    fn parse_organic_atom(&mut self) -> Result<AtomSpec, ChemError> {
        let rest = &self.input[self.pos..];
        let (symbol, aromatic, len) = if rest.starts_with("Cl") {
            ("Cl", false, 2)
        } else if rest.starts_with("Br") {
            ("Br", false, 2)
        } else {
            match rest.as_bytes()[0] {
                b'B' => ("B", false, 1),
                b'C' => ("C", false, 1),
                b'N' => ("N", false, 1),
                b'O' => ("O", false, 1),
                b'P' => ("P", false, 1),
                b'S' => ("S", false, 1),
                b'F' => ("F", false, 1),
                b'I' => ("I", false, 1),
                b'b' => ("B", true, 1),
                b'c' => ("C", true, 1),
                b'n' => ("N", true, 1),
                b'o' => ("O", true, 1),
                b'p' => ("P", true, 1),
                b's' => ("S", true, 1),
                c => {
                    return Err(self.err(&format!("unexpected character {:?}", c as char)));
                }
            }
        };
        self.pos += len;
        Ok(AtomSpec {
            symbol: symbol.into(),
            aromatic,
            isotope: None,
            charge: 0,
            explicit_h: None,
            chirality: None,
            atom_class: None,
            bracket: false,
        })
    }

    /// 角括弧原子: `[` isotope? symbol chiral? hcount? charge? class? `]`
    /// (先頭の `[` は消費済み)
    fn parse_bracket_atom(&mut self) -> Result<AtomSpec, ChemError> {
        // isotope
        let mut isotope: Option<u16> = None;
        while let Some(c @ b'0'..=b'9') = self.peek() {
            let v = isotope.unwrap_or(0) as u32 * 10 + (c - b'0') as u32;
            if v > 999 {
                return Err(self.err("isotope out of range"));
            }
            isotope = Some(v as u16);
            self.pos += 1;
        }

        // element symbol
        let (symbol, aromatic) = self.parse_bracket_symbol()?;

        // chirality
        let mut chirality = None;
        if self.peek() == Some(b'@') {
            self.pos += 1;
            if self.peek() == Some(b'@') {
                self.pos += 1;
                chirality = Some(Chirality::Clockwise);
            } else {
                chirality = Some(Chirality::Anticlockwise);
            }
            let rest = &self.input[self.pos..];
            for tag in ["TH", "AL", "SP", "TB", "OH"] {
                if rest.starts_with(tag) {
                    return Err(self.err(&format!("extended chirality @{tag} is not supported")));
                }
            }
        }

        // hcount
        let mut explicit_h: u8 = 0;
        if self.peek() == Some(b'H') {
            self.pos += 1;
            let mut n: Option<u32> = None;
            while let Some(c @ b'0'..=b'9') = self.peek() {
                n = Some(n.unwrap_or(0) * 10 + (c - b'0') as u32);
                self.pos += 1;
            }
            let n = n.unwrap_or(1);
            if n > 9 {
                return Err(self.err("H count out of range"));
            }
            explicit_h = n as u8;
        }

        // charge
        let mut charge: i8 = 0;
        if let Some(sign @ (b'+' | b'-')) = self.peek() {
            self.pos += 1;
            let mut repeats: i32 = 1;
            while self.peek() == Some(sign) {
                repeats += 1;
                self.pos += 1;
            }
            let mut magnitude = repeats;
            if repeats == 1 {
                let mut digits: Option<i32> = None;
                while let Some(c @ b'0'..=b'9') = self.peek() {
                    digits = Some(digits.unwrap_or(0) * 10 + (c - b'0') as i32);
                    self.pos += 1;
                }
                if let Some(d) = digits {
                    magnitude = d;
                }
            }
            if magnitude > 15 {
                return Err(self.err("charge out of range"));
            }
            charge = if sign == b'+' {
                magnitude as i8
            } else {
                -(magnitude as i8)
            };
        }

        // atom class
        let mut atom_class: Option<u32> = None;
        if self.peek() == Some(b':') {
            self.pos += 1;
            let mut n: Option<u32> = None;
            while let Some(c @ b'0'..=b'9') = self.peek() {
                n = Some(n.unwrap_or(0) * 10 + (c - b'0') as u32);
                self.pos += 1;
            }
            atom_class = Some(n.ok_or_else(|| self.err("':' in bracket requires digits"))?);
        }

        if self.bump() != Some(b']') {
            return Err(self.err("expected ']'"));
        }

        Ok(AtomSpec {
            symbol,
            aromatic,
            isotope,
            charge,
            explicit_h: Some(explicit_h),
            chirality,
            atom_class,
            bracket: true,
        })
    }

    /// 角括弧内の元素記号。最長一致 (2 文字優先)。
    fn parse_bracket_symbol(&mut self) -> Result<(String, bool), ChemError> {
        let rest = &self.bytes[self.pos..];
        match rest.first() {
            Some(b'*') => {
                self.pos += 1;
                Ok(("*".into(), false))
            }
            Some(c) if c.is_ascii_uppercase() => {
                // 2 文字記号を優先。ただし [CH4] の H は hcount なので、
                // 2 文字目が続く場合は有効な元素のときのみ採用する。
                if let Some(c2) = rest.get(1).filter(|c| c.is_ascii_lowercase()) {
                    let two = format!("{}{}", *c as char, *c2 as char);
                    if atomic_number(&two).is_some() {
                        self.pos += 2;
                        return Ok((two, false));
                    }
                }
                let one = (*c as char).to_string();
                if atomic_number(&one).is_none() {
                    return Err(self.err(&format!("unknown element {one:?}")));
                }
                self.pos += 1;
                Ok((one, false))
            }
            Some(c) if c.is_ascii_lowercase() => {
                if let Some(c2) = rest.get(1).filter(|c| c.is_ascii_lowercase()) {
                    let two = format!("{}{}", *c as char, *c2 as char);
                    if let Some(sym) = aromatic_symbol(&two) {
                        self.pos += 2;
                        return Ok((sym.into(), true));
                    }
                }
                let one = (*c as char).to_string();
                let sym = aromatic_symbol(&one)
                    .ok_or_else(|| self.err(&format!("invalid aromatic symbol {one:?}")))?;
                self.pos += 1;
                Ok((sym.into(), true))
            }
            _ => Err(self.err("expected element symbol in brackets")),
        }
    }

    /// 原子を追加し、直前原子との結合を張る。
    fn add_atom(&mut self, atom: AtomSpec) -> Result<(), ChemError> {
        let idx = self.atoms.len();
        self.atoms.push(atom);
        self.neighbor_order.push(Vec::new());
        if let Some(p) = self.prev {
            let kind = self.pending_bond.take().unwrap_or(BondKind::Elided);
            let bi = self.bonds.len();
            self.bonds.push(BondSpec {
                a: p,
                b: idx,
                kind,
                ring_closure: None,
            });
            self.neighbor_order[p].push(bi);
            self.neighbor_order[idx].push(bi);
        } else if self.pending_bond.is_some() {
            return Err(self.err("bond symbol without preceding atom"));
        }
        self.prev = Some(idx);
        self.after_dot = false;
        Ok(())
    }

    /// 環閉じ数字の処理 (開き / 閉じ)。
    fn ring_bond(&mut self, num: u16) -> Result<(), ChemError> {
        let cur = self
            .prev
            .ok_or_else(|| self.err("ring bond digit without preceding atom"))?;
        let kind = self.pending_bond.take().unwrap_or(BondKind::Elided);

        if let Some(open) = self.ring_map.remove(&num) {
            if open.atom == cur {
                return Err(self.err(&format!("ring bond {num} closes on the same atom")));
            }
            if self
                .bonds
                .iter()
                .any(|b| (b.a == open.atom && b.b == cur) || (b.a == cur && b.b == open.atom))
            {
                return Err(self.err(&format!("duplicate bond via ring closure {num}")));
            }
            let merged = merge_ring_bond_kind(open.kind, kind).ok_or_else(|| {
                self.err(&format!("conflicting bond symbols on ring closure {num}"))
            })?;
            let bi = self.bonds.len();
            let opened_with_order = matches!(
                open.kind,
                BondKind::Single
                    | BondKind::Double
                    | BondKind::Triple
                    | BondKind::Quadruple
                    | BondKind::Aromatic
            );
            self.bonds.push(BondSpec {
                a: open.atom,
                b: cur,
                kind: merged,
                ring_closure: Some(RingClosure {
                    num,
                    opened_with_order,
                }),
            });
            self.neighbor_order[open.atom][open.slot] = bi;
            self.neighbor_order[cur].push(bi);
        } else {
            let slot = self.neighbor_order[cur].len();
            self.neighbor_order[cur].push(usize::MAX); // 閉じ時に patch
            self.ring_map.insert(
                num,
                RingOpen {
                    atom: cur,
                    kind,
                    slot,
                },
            );
        }
        Ok(())
    }
}

/// 環閉じの両端で指定された結合記号をマージする。矛盾する場合 None。
fn merge_ring_bond_kind(a: BondKind, b: BondKind) -> Option<BondKind> {
    use BondKind::*;
    match (a, b) {
        (Elided, x) | (x, Elided) => Some(x),
        (x, y) if x == y => Some(x),
        // 方向付き結合は両端から見て相補 (/ と \) になる。開き側の向きを採用。
        (Up, Down) | (Down, Up) => Some(a),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::{BondKind::*, Chirality::*};

    fn parse(s: &str) -> ParsedMolecule {
        parse_smiles(s).unwrap_or_else(|e| panic!("{s}: {e}"))
    }

    fn symbols(m: &ParsedMolecule) -> Vec<&str> {
        m.atoms.iter().map(|a| a.symbol.as_str()).collect()
    }

    fn bond_kinds(m: &ParsedMolecule) -> Vec<BondKind> {
        m.bonds.iter().map(|b| b.kind).collect()
    }

    // ---- 正常系 ----

    #[test]
    fn simple_chains() {
        assert_eq!(symbols(&parse("C")), ["C"]);
        assert_eq!(symbols(&parse("CC")), ["C", "C"]);
        assert_eq!(symbols(&parse("CCO")), ["C", "C", "O"]);
        assert_eq!(symbols(&parse("NCC")), ["N", "C", "C"]);
        let m = parse("CCCC");
        assert_eq!(m.bonds.len(), 3);
        assert!(m
            .bonds
            .iter()
            .all(|b| b.kind == Elided && b.ring_closure.is_none()));
    }

    #[test]
    fn bond_symbols() {
        assert_eq!(bond_kinds(&parse("C=C")), [Double]);
        assert_eq!(bond_kinds(&parse("C#N")), [Triple]);
        assert_eq!(bond_kinds(&parse("C-C")), [Single]);
        assert_eq!(bond_kinds(&parse("S=C=S")), [Double, Double]);
        assert_eq!(bond_kinds(&parse("C=CC#N")), [Double, Elided, Triple]);
    }

    #[test]
    fn directional_bonds() {
        assert_eq!(bond_kinds(&parse("F/C=C/F")), [Up, Double, Up]);
        assert_eq!(bond_kinds(&parse(r"F/C=C\F")), [Up, Double, Down]);
    }

    #[test]
    fn branches() {
        let m = parse("CC(C)C");
        assert_eq!(m.atoms.len(), 4);
        // 中心原子 1 が 0, 2, 3 と結合
        let mut partners: Vec<usize> = m
            .bonds
            .iter()
            .map(|b| if b.a == 1 { b.b } else { b.a })
            .collect();
        partners.sort_unstable();
        assert_eq!(partners, [0, 2, 3]);

        let m = parse("CC(C)(C)C"); // ネオペンタン
        assert_eq!(m.atoms.len(), 5);
        assert_eq!(m.bonds.len(), 4);

        let m = parse("CC(C(C)C)C"); // ネスト分岐
        assert_eq!(m.atoms.len(), 6);
        assert_eq!(m.bonds.len(), 5);

        let m = parse("C(=O)O"); // 分岐先頭の結合記号
        assert_eq!(bond_kinds(&m), [Double, Elided]);
    }

    #[test]
    fn rings() {
        let m = parse("C1CC1");
        assert_eq!(m.bonds.len(), 3);
        assert_eq!(
            m.bonds.iter().filter(|b| b.ring_closure.is_some()).count(),
            1
        );
        let rc = m.bonds.iter().find(|b| b.ring_closure.is_some()).unwrap();
        assert_eq!((rc.a, rc.b), (0, 2));

        let m = parse("C1CCCCC1"); // シクロヘキサン
        assert_eq!(m.atoms.len(), 6);
        assert_eq!(m.bonds.len(), 6);

        // 番号の再利用
        let m = parse("C1CC1C1CC1");
        assert_eq!(
            m.bonds.iter().filter(|b| b.ring_closure.is_some()).count(),
            2
        );

        // トリフェニルホスフィン: 同じ番号を 3 回再利用
        let m = parse("P(c1ccccc1)(c1ccccc1)c1ccccc1");
        assert_eq!(m.atoms.len(), 19);
        assert_eq!(m.bonds.len(), 21);
    }

    #[test]
    fn ring_bond_symbols() {
        // 両側指定 (一致)
        assert!(bond_kinds(&parse("C=1CC=1")).contains(&Double));
        // 片側指定
        assert!(bond_kinds(&parse("C1CC=1")).contains(&Double));
        assert!(bond_kinds(&parse("C=1CC1")).contains(&Double));
        // %nn 表記
        let m = parse("C%12CC%12");
        assert_eq!(m.bonds.len(), 3);
        // 芳香族結合記号
        let m = parse("c:1:c:c:c:c:c:1");
        assert_eq!(m.atoms.len(), 6);
        assert!(m.bonds.iter().all(|b| b.kind == Aromatic));
    }

    #[test]
    fn aromatic_atoms() {
        let m = parse("c1ccccc1"); // ベンゼン
        assert_eq!(m.atoms.len(), 6);
        assert!(m.atoms.iter().all(|a| a.aromatic && a.symbol == "C"));
        assert!(m.bonds.iter().all(|b| b.kind == Elided));

        let m = parse("c1cc[nH]c1"); // ピロール
        assert_eq!(m.atoms[3].symbol, "N");
        assert!(m.atoms[3].aromatic);
        assert_eq!(m.atoms[3].explicit_h, Some(1));

        let m = parse("c1cc[se]c1"); // セレノフェン
        assert_eq!(m.atoms[3].symbol, "Se");
        assert!(m.atoms[3].aromatic);

        let m = parse("Clc1ccccc1"); // 2 文字有機サブセット + 芳香環
        assert_eq!(m.atoms[0].symbol, "Cl");
        assert_eq!(symbols(&parse("BrCCBr")), ["Br", "C", "C", "Br"]);
    }

    #[test]
    fn bracket_atoms() {
        let a = &parse("[CH4]").atoms[0];
        assert_eq!(
            (a.symbol.as_str(), a.explicit_h, a.bracket),
            ("C", Some(4), true)
        );

        let a = &parse("[H]").atoms[0];
        assert_eq!((a.symbol.as_str(), a.explicit_h), ("H", Some(0)));

        let a = &parse("[2H]").atoms[0]; // 重水素
        assert_eq!((a.symbol.as_str(), a.isotope), ("H", Some(2)));

        let a = &parse("[13CH4]").atoms[0];
        assert_eq!(a.isotope, Some(13));

        let a = &parse("[Si](C)(C)(C)C").atoms[0];
        assert_eq!((a.symbol.as_str(), a.explicit_h), ("Si", Some(0)));

        let a = &parse("[SiH3]C").atoms[0];
        assert_eq!(a.explicit_h, Some(3));

        for (smi, sym) in [
            ("[Se]", "Se"),
            ("[Te]", "Te"),
            ("[As]", "As"),
            ("[Hg]", "Hg"),
        ] {
            assert_eq!(parse(smi).atoms[0].symbol, sym);
        }
        // [te] は RDKit 拡張
        let a = &parse("[te]").atoms[0];
        assert!(a.aromatic);
        assert_eq!(a.symbol, "Te");
    }

    #[test]
    fn charges() {
        assert_eq!(parse("[O-]").atoms[0].charge, -1);
        assert_eq!(parse("[NH4+]").atoms[0].charge, 1);
        assert_eq!(parse("[NH4+]").atoms[0].explicit_h, Some(4));
        assert_eq!(parse("[Ca+2]").atoms[0].charge, 2);
        assert_eq!(parse("[Fe+3]").atoms[0].charge, 3);
        assert_eq!(parse("[O-2]").atoms[0].charge, -2);
        assert_eq!(parse("[O--]").atoms[0].charge, -2); // 重ね書き表記
        assert_eq!(parse("[N+](=O)([O-])C").atoms[0].charge, 1);
    }

    #[test]
    fn chirality() {
        let m = parse("N[C@H](C)C(=O)O"); // L-アラニン
        assert_eq!(m.atoms[1].chirality, Some(Anticlockwise));
        assert_eq!(m.atoms[1].explicit_h, Some(1));
        let m = parse("N[C@@H](C)C(=O)O");
        assert_eq!(m.atoms[1].chirality, Some(Clockwise));
        let m = parse("[C@](N)(C)(O)F");
        assert_eq!(m.atoms[0].chirality, Some(Anticlockwise));
    }

    #[test]
    fn atom_class_and_wildcard() {
        assert_eq!(parse("[CH4:2]").atoms[0].atom_class, Some(2));
        assert_eq!(parse("*").atoms[0].symbol, "*");
        assert_eq!(parse("[*]").atoms[0].symbol, "*");
        assert_eq!(parse("C*").bonds.len(), 1);
    }

    #[test]
    fn multicomponent() {
        let m = parse("[Na+].[Cl-]");
        assert_eq!(m.atoms.len(), 2);
        assert_eq!(m.bonds.len(), 0);

        let m = parse("C.C.C");
        assert_eq!(m.atoms.len(), 3);
        assert_eq!(m.bonds.len(), 0);

        let m = parse("CC(=O)[O-].[Na+]"); // 酢酸ナトリウム
        assert_eq!(m.atoms.len(), 5);
        assert_eq!(m.bonds.len(), 3);
    }

    #[test]
    fn neighbor_order_records_ring_digit_position() {
        // C1CC1: 原子 0 の近傍順は [環閉じ結合(数字位置), 鎖結合] の順
        let m = parse("C1CC1");
        let rc_bond = m
            .bonds
            .iter()
            .position(|b| b.ring_closure.is_some())
            .unwrap();
        assert_eq!(m.neighbor_order[0][0], rc_bond);
        // 原子 2 (閉じ側) では最後
        assert_eq!(*m.neighbor_order[2].last().unwrap(), rc_bond);
    }

    #[test]
    fn realistic_molecules() {
        for smi in [
            "CCS(=O)(=O)O",              // エタンスルホン酸
            "OC(=O)c1ccccc1",            // 安息香酸
            "N#Cc1ccccc1",               // ベンゾニトリル
            "Cn1c(=S)oc2ccccc21",        // N-メチルベンゾオキサゾールチオン
            "c1ccc2ocnc2c1",             // ベンゾオキサゾール
            "CC(C)(C)OC(=O)N",           // Boc-アミン
            "[O-][n+]1ccccc1",           // ピリジン N-オキシド
            "C/C=C/C=C/C",               // ジエン (E,E)
            "O=S(=O)(c1ccccc1)N1CCCCC1", // スルホンアミド
        ] {
            parse(smi);
        }
    }

    // ---- 異常系 ----

    #[test]
    fn errors() {
        let bad = [
            "",               // 空
            "C(",             // 括弧未閉じ
            "C)",             // 括弧過剰
            "(C)C",           // 先頭分岐
            "C1CC",           // 環未閉じ
            "C1CC2",          // 環 2 つ未閉じ
            "C=",             // 末尾結合
            "=C",             // 先頭結合
            "C==C",           // 結合記号重複
            "C=-C",           // 結合記号重複 (異種)
            "[C",             // 角括弧未閉じ
            "[]",             // 空角括弧
            "[Xx]",           // 未知元素
            "[cl]",           // 無効な芳香族記号
            "H",              // 裸の水素
            "hello",          // 非 SMILES
            "C11",            // 自己結合
            "C12CC12",        // 重複結合
            "C%1C",           // % の桁不足
            "C=1CC-1",        // 環閉じ結合記号の矛盾
            "C=(C)C",         // 分岐前の結合記号
            "C .C",           // 空白
            "C..C",           // 連続ドット
            "C.",             // 末尾ドット
            ".C",             // 先頭ドット
            "[C@TH1](N)(C)O", // 拡張キラリティ未対応
            "[CH+2H]",        // 電荷の後に H
            "1CC1",           // 原子より先の環数字
            "[C:]",           // クラス数字なし
            "[Ca+20]",        // 電荷範囲外
        ];
        for smi in bad {
            assert!(
                parse_smiles(smi).is_err(),
                "expected parse error for {smi:?}"
            );
        }
    }
}
