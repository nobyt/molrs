//! SMARTS サブセットのパーサとマッチャ (C10 / RDKIT_RUST_PLAN.md R2.1 の前倒し)。
//!
//! ETKDG トーションライブラリ (365 パターン) が使う機能を過不足なく実装する:
//!
//! - 原子プリミティブ: `*` `a` `A`、元素 (大文字=脂肪族, 小文字=芳香族)、
//!   `#n` (原子番号)、`X<n>` (総結合数)、`H<n>` (総 H 数)、`x<n>` (環結合数)、
//!   `r<n>` / `r` (n 員環に属する / 環内)、`R<n>` (所属環数)、`^<n>` (混成)、
//!   電荷 `+`/`-`/`+n`、再帰 SMARTS `$(...)` (入れ子可)、原子マップ `:n`
//! - 論理: `!` (否定) > `&`・隣接 (AND) > `,` (OR) > `;` (AND)
//! - 結合: `-` `=` `#` `:` `~` `@` と同じ論理演算子。無指定は「単結合または芳香族」
//! - 分岐 `()` と環閉じ数字
//!
//! マッチは `MoleculeGraph` の全原子 (明示 H 含む) に対して行う。
//! H 数・結合数は明示 H ノードを数えるので RDKit の AddHs 済み分子と同じ意味になる。

use std::collections::HashMap;

use crate::graph::MoleculeGraph;

// ---------------------------------------------------------------------------
// クエリ表現
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub enum AtomExpr {
    Any,               // *
    Aromatic,          // a
    Aliphatic,         // A
    AtomicNum(u8),     // #n (芳香族性を問わない)
    ElemAliphatic(u8), // C, N, ...
    ElemAromatic(u8),  // c, n, ...
    TotalConn(u8),     // X<n>
    TotalH(u8),        // H<n>
    RingBondCount(u8), // x<n>
    InRingSize(u8),    // r<n>
    InRing,            // r / R (数字なし)
    RingCount(u8),     // R<n>
    Hybridization(u8), // ^<n> (1=sp, 2=sp2, 3=sp3)
    Charge(i8),
    Recursive(usize), // recursives[idx] のパターンに根付きマッチ
    Not(Box<AtomExpr>),
    And(Vec<AtomExpr>),
    Or(Vec<AtomExpr>),
}

#[derive(Debug, Clone)]
pub enum BondExpr {
    Single,
    Double,
    Triple,
    AromaticB,
    Any,     // ~
    Ring,    // @
    Default, // 無指定: 単結合または芳香族
    Not(Box<BondExpr>),
    And(Vec<BondExpr>),
    Or(Vec<BondExpr>),
}

#[derive(Debug, Clone)]
pub struct SmartsPattern {
    pub atoms: Vec<AtomExpr>,
    /// 原子マップ番号 (`[C:2]` → Some(2))
    pub atom_maps: Vec<Option<u32>>,
    pub bonds: Vec<(usize, usize, BondExpr)>,
    /// 再帰 SMARTS の内部パターン
    recursives: Vec<SmartsPattern>,
    adj: Vec<Vec<(usize, usize)>>, // atom → [(相手, bond idx)]
}

// ---------------------------------------------------------------------------
// パーサ
// ---------------------------------------------------------------------------

struct Parser<'a> {
    b: &'a [u8],
    pos: usize,
    recursives: Vec<SmartsPattern>,
}

pub fn parse_smarts(s: &str) -> Result<SmartsPattern, String> {
    let mut p = Parser {
        b: s.as_bytes(),
        pos: 0,
        recursives: Vec::new(),
    };
    let mut pat = p.parse_chain()?;
    if p.pos != p.b.len() {
        return Err(format!("unexpected char at {} in {s:?}", p.pos));
    }
    pat.recursives = p.recursives;
    Ok(pat)
}

impl<'a> Parser<'a> {
    fn peek(&self) -> Option<u8> {
        self.b.get(self.pos).copied()
    }

    fn err(&self, msg: &str) -> String {
        format!(
            "{msg} at {} in {:?}",
            self.pos,
            std::str::from_utf8(self.b).unwrap_or("")
        )
    }

    fn num(&mut self) -> Option<u32> {
        let start = self.pos;
        while matches!(self.peek(), Some(b'0'..=b'9')) {
            self.pos += 1;
        }
        if self.pos == start {
            None
        } else {
            std::str::from_utf8(&self.b[start..self.pos])
                .ok()?
                .parse()
                .ok()
        }
    }

    /// 鎖 (分岐・環閉じ込み) をパースしてパターンを返す。
    fn parse_chain(&mut self) -> Result<SmartsPattern, String> {
        let mut atoms: Vec<AtomExpr> = Vec::new();
        let mut maps: Vec<Option<u32>> = Vec::new();
        let mut bonds: Vec<(usize, usize, BondExpr)> = Vec::new();
        let mut prev: Option<usize> = None;
        let mut stack: Vec<usize> = Vec::new();
        let mut pending: Option<BondExpr> = None;
        let mut ring_open: HashMap<u32, (usize, Option<BondExpr>)> = HashMap::new();

        loop {
            match self.peek() {
                None => break,
                Some(b')') if !stack.is_empty() => {
                    self.pos += 1;
                    prev = Some(stack.pop().expect("checked"));
                }
                Some(b')') => break, // 再帰 SMARTS の終端 (呼び出し側が消費)
                Some(b'(') => {
                    self.pos += 1;
                    let p = prev.ok_or_else(|| self.err("branch without atom"))?;
                    stack.push(p);
                }
                Some(c @ (b'-' | b'=' | b'#' | b':' | b'~' | b'@' | b'!' | b'/' | b'\\'))
                    if self.is_bond_context(c) =>
                {
                    if pending.is_some() {
                        return Err(self.err("duplicate bond expr"));
                    }
                    pending = Some(self.parse_bond_expr()?);
                }
                Some(b'0'..=b'9') => {
                    // 環閉じ
                    let d = (self.b[self.pos] - b'0') as u32;
                    self.pos += 1;
                    let cur = prev.ok_or_else(|| self.err("ring digit without atom"))?;
                    let bexpr = pending.take();
                    if let Some((open_atom, open_bond)) = ring_open.remove(&d) {
                        let be = bexpr.or(open_bond).unwrap_or(BondExpr::Default);
                        bonds.push((open_atom, cur, be));
                    } else {
                        ring_open.insert(d, (cur, bexpr));
                    }
                }
                _ => {
                    let (expr, map) = self.parse_atom()?;
                    let idx = atoms.len();
                    atoms.push(expr);
                    maps.push(map);
                    if let Some(p) = prev {
                        bonds.push((p, idx, pending.take().unwrap_or(BondExpr::Default)));
                    } else if pending.is_some() {
                        return Err(self.err("bond without preceding atom"));
                    }
                    prev = Some(idx);
                }
            }
        }
        if !ring_open.is_empty() {
            return Err(self.err("unclosed ring bond"));
        }
        if !stack.is_empty() {
            return Err(self.err("unclosed branch"));
        }
        if atoms.is_empty() {
            return Err(self.err("empty pattern"));
        }
        let mut adj = vec![Vec::new(); atoms.len()];
        for (ei, &(a, b, _)) in bonds.iter().enumerate() {
            adj[a].push((b, ei));
            adj[b].push((a, ei));
        }
        Ok(SmartsPattern {
            atoms,
            atom_maps: maps,
            bonds,
            recursives: Vec::new(), // トップレベルで設定
            adj,
        })
    }

    /// 現在位置の記号が結合表現か (原子の '!' 等と区別するための文脈判定)。
    /// 結合表現は原子の直後 (または閉じ括弧の直後) にだけ現れる。
    fn is_bond_context(&self, c: u8) -> bool {
        match c {
            b'-' | b'=' | b':' | b'~' | b'@' | b'/' | b'\\' => true,
            // '#' は結合 (三重) と #n (原子番号) の両方があるが、
            // ブラケット外の '#' は結合のみ
            b'#' => true,
            // '!' は結合否定 (!@ 等)。ブラケット外の ! は結合表現
            b'!' => true,
            _ => false,
        }
    }

    /// 結合表現 (論理込み)。
    fn parse_bond_expr(&mut self) -> Result<BondExpr, String> {
        // ; 低位 AND
        let mut semi: Vec<BondExpr> = vec![self.parse_bond_or()?];
        while self.peek() == Some(b';') {
            self.pos += 1;
            semi.push(self.parse_bond_or()?);
        }
        Ok(if semi.len() == 1 {
            semi.pop().expect("one")
        } else {
            BondExpr::And(semi)
        })
    }

    fn parse_bond_or(&mut self) -> Result<BondExpr, String> {
        let mut or: Vec<BondExpr> = vec![self.parse_bond_and()?];
        while self.peek() == Some(b',') {
            self.pos += 1;
            or.push(self.parse_bond_and()?);
        }
        Ok(if or.len() == 1 {
            or.pop().expect("one")
        } else {
            BondExpr::Or(or)
        })
    }

    fn parse_bond_and(&mut self) -> Result<BondExpr, String> {
        let mut and: Vec<BondExpr> = vec![self.parse_bond_factor()?];
        loop {
            match self.peek() {
                Some(b'&') => {
                    self.pos += 1;
                    and.push(self.parse_bond_factor()?);
                }
                // 隣接 = 暗黙 AND (例: `!@-`)
                Some(b'-' | b'=' | b'#' | b':' | b'~' | b'@' | b'!' | b'/' | b'\\') => {
                    and.push(self.parse_bond_factor()?);
                }
                _ => break,
            }
        }
        Ok(if and.len() == 1 {
            and.pop().expect("one")
        } else {
            BondExpr::And(and)
        })
    }

    fn parse_bond_factor(&mut self) -> Result<BondExpr, String> {
        match self.peek() {
            Some(b'!') => {
                self.pos += 1;
                Ok(BondExpr::Not(Box::new(self.parse_bond_factor()?)))
            }
            Some(b'-') => {
                self.pos += 1;
                Ok(BondExpr::Single)
            }
            Some(b'/') | Some(b'\\') => {
                self.pos += 1;
                Ok(BondExpr::Single) // 方向は無視して単結合扱い
            }
            Some(b'=') => {
                self.pos += 1;
                Ok(BondExpr::Double)
            }
            Some(b'#') => {
                self.pos += 1;
                Ok(BondExpr::Triple)
            }
            Some(b':') => {
                self.pos += 1;
                Ok(BondExpr::AromaticB)
            }
            Some(b'~') => {
                self.pos += 1;
                Ok(BondExpr::Any)
            }
            Some(b'@') => {
                self.pos += 1;
                Ok(BondExpr::Ring)
            }
            _ => Err(self.err("expected bond primitive")),
        }
    }

    /// 原子 1 つ ([...] または裸の有機サブセット/a/A/*)。
    fn parse_atom(&mut self) -> Result<(AtomExpr, Option<u32>), String> {
        match self.peek() {
            Some(b'[') => {
                self.pos += 1;
                let expr = self.parse_atom_expr()?;
                // 原子マップ
                let map = if self.peek() == Some(b':') {
                    self.pos += 1;
                    Some(self.num().ok_or_else(|| self.err("map number expected"))?)
                } else {
                    None
                };
                if self.peek() != Some(b']') {
                    return Err(self.err("expected ']'"));
                }
                self.pos += 1;
                Ok((expr, map))
            }
            Some(b'*') => {
                self.pos += 1;
                Ok((AtomExpr::Any, None))
            }
            Some(b'a') => {
                self.pos += 1;
                Ok((AtomExpr::Aromatic, None))
            }
            Some(b'A') => {
                self.pos += 1;
                Ok((AtomExpr::Aliphatic, None))
            }
            _ => {
                // 裸の有機サブセット元素
                let e = self.parse_bare_element()?;
                Ok((e, None))
            }
        }
    }

    fn parse_bare_element(&mut self) -> Result<AtomExpr, String> {
        let rest = &self.b[self.pos..];
        let two = |s: &mut Self, n: u8| {
            s.pos += 2;
            AtomExpr::ElemAliphatic(n)
        };
        if rest.starts_with(b"Cl") {
            return Ok(two(self, 17));
        }
        if rest.starts_with(b"Br") {
            return Ok(two(self, 35));
        }
        let e = match rest.first() {
            Some(b'B') => AtomExpr::ElemAliphatic(5),
            Some(b'C') => AtomExpr::ElemAliphatic(6),
            Some(b'N') => AtomExpr::ElemAliphatic(7),
            Some(b'O') => AtomExpr::ElemAliphatic(8),
            Some(b'P') => AtomExpr::ElemAliphatic(15),
            Some(b'S') => AtomExpr::ElemAliphatic(16),
            Some(b'F') => AtomExpr::ElemAliphatic(9),
            Some(b'I') => AtomExpr::ElemAliphatic(53),
            Some(b'b') => AtomExpr::ElemAromatic(5),
            Some(b'c') => AtomExpr::ElemAromatic(6),
            Some(b'n') => AtomExpr::ElemAromatic(7),
            Some(b'o') => AtomExpr::ElemAromatic(8),
            Some(b'p') => AtomExpr::ElemAromatic(15),
            Some(b's') => AtomExpr::ElemAromatic(16),
            _ => return Err(self.err("expected atom")),
        };
        self.pos += 1;
        Ok(e)
    }

    /// ブラケット内の原子式 (マップ・']' の手前まで)。
    fn parse_atom_expr(&mut self) -> Result<AtomExpr, String> {
        // RDKit 互換: ブラケット式の先頭が 'H' なら水素原子 (#1) を意味する
        // ([H:4] はトーションライブラリで明示 H を指す)。
        // それ以外の位置の H は「総 H 数」プリミティブ
        if self.peek() == Some(b'H')
            && matches!(
                self.b.get(self.pos + 1),
                None | Some(b']' | b':' | b';' | b',' | b'&')
            )
        {
            self.pos += 1;
            return Ok(AtomExpr::AtomicNum(1));
        }
        // ; 低位 AND
        let mut semi = vec![self.parse_atom_or()?];
        while self.peek() == Some(b';') {
            self.pos += 1;
            semi.push(self.parse_atom_or()?);
        }
        Ok(if semi.len() == 1 {
            semi.pop().expect("one")
        } else {
            AtomExpr::And(semi)
        })
    }

    fn parse_atom_or(&mut self) -> Result<AtomExpr, String> {
        let mut or = vec![self.parse_atom_and()?];
        while self.peek() == Some(b',') {
            self.pos += 1;
            or.push(self.parse_atom_and()?);
        }
        Ok(if or.len() == 1 {
            or.pop().expect("one")
        } else {
            AtomExpr::Or(or)
        })
    }

    fn parse_atom_and(&mut self) -> Result<AtomExpr, String> {
        let mut and = vec![self.parse_atom_factor()?];
        loop {
            match self.peek() {
                Some(b'&') => {
                    self.pos += 1;
                    and.push(self.parse_atom_factor()?);
                }
                // 隣接 = 暗黙 AND。式の終端記号以外なら続きのプリミティブ
                Some(c) if c != b';' && c != b',' && c != b']' && c != b':' => {
                    and.push(self.parse_atom_factor()?);
                }
                _ => break,
            }
        }
        Ok(if and.len() == 1 {
            and.pop().expect("one")
        } else {
            AtomExpr::And(and)
        })
    }

    fn parse_atom_factor(&mut self) -> Result<AtomExpr, String> {
        match self.peek() {
            Some(b'!') => {
                self.pos += 1;
                Ok(AtomExpr::Not(Box::new(self.parse_atom_factor()?)))
            }
            Some(b'$') => {
                self.pos += 1;
                if self.peek() != Some(b'(') {
                    return Err(self.err("expected '(' after '$'"));
                }
                self.pos += 1;
                let inner = self.parse_chain()?;
                if self.peek() != Some(b')') {
                    return Err(self.err("expected ')' after recursive SMARTS"));
                }
                self.pos += 1;
                // 入れ子の $() は同じパーサを通るため self.recursives に
                // 正しい ID で登録済み。ここでは inner 自身を登録するだけでよい
                let id = self.recursives.len();
                self.recursives.push(inner);
                Ok(AtomExpr::Recursive(id))
            }
            Some(b'*') => {
                self.pos += 1;
                Ok(AtomExpr::Any)
            }
            Some(b'#') => {
                self.pos += 1;
                let n = self
                    .num()
                    .ok_or_else(|| self.err("atomic number expected"))?;
                Ok(AtomExpr::AtomicNum(n as u8))
            }
            Some(b'X') => {
                self.pos += 1;
                let n = self.num().unwrap_or(1);
                Ok(AtomExpr::TotalConn(n as u8))
            }
            Some(b'H') => {
                // H は「総 H 数」プリミティブ。ただし他のプリミティブが
                // 先行しない位置では水素元素 ([H:4] 等)
                self.pos += 1;
                match self.num() {
                    Some(n) => Ok(AtomExpr::TotalH(n as u8)),
                    None => Ok(AtomExpr::TotalH(1)),
                }
            }
            Some(b'x') => {
                self.pos += 1;
                let n = self.num().ok_or_else(|| self.err("x needs a digit"))?;
                Ok(AtomExpr::RingBondCount(n as u8))
            }
            Some(b'r') => {
                self.pos += 1;
                match self.num() {
                    Some(n) => Ok(AtomExpr::InRingSize(n as u8)),
                    None => Ok(AtomExpr::InRing),
                }
            }
            Some(b'R') => {
                self.pos += 1;
                match self.num() {
                    Some(0) => Ok(AtomExpr::Not(Box::new(AtomExpr::InRing))),
                    Some(n) => Ok(AtomExpr::RingCount(n as u8)),
                    None => Ok(AtomExpr::InRing),
                }
            }
            Some(b'^') => {
                self.pos += 1;
                let n = self.num().ok_or_else(|| self.err("^ needs a digit"))?;
                Ok(AtomExpr::Hybridization(n as u8))
            }
            Some(b'+') => {
                self.pos += 1;
                let mut c = 1i8;
                while self.peek() == Some(b'+') {
                    self.pos += 1;
                    c += 1;
                }
                if let Some(n) = self.num() {
                    c = n as i8;
                }
                Ok(AtomExpr::Charge(c))
            }
            Some(b'-') => {
                self.pos += 1;
                let mut c = -1i8;
                while self.peek() == Some(b'-') {
                    self.pos += 1;
                    c -= 1;
                }
                if let Some(n) = self.num() {
                    c = -(n as i8);
                }
                Ok(AtomExpr::Charge(c))
            }
            Some(b'a') => {
                self.pos += 1;
                Ok(AtomExpr::Aromatic)
            }
            Some(b'A') => {
                self.pos += 1;
                Ok(AtomExpr::Aliphatic)
            }
            Some(c) if c.is_ascii_uppercase() => {
                // 元素 (2 文字対応)。H は上で処理済みなのでここは H 以外
                let rest = &self.b[self.pos..];
                let table: &[(&[u8], u8)] = &[
                    (b"Cl", 17),
                    (b"Br", 35),
                    (b"Si", 14),
                    (b"Se", 34),
                    (b"Sn", 50),
                    (b"Sb", 51),
                    (b"Te", 52),
                    (b"As", 33),
                    (b"B", 5),
                    (b"C", 6),
                    (b"N", 7),
                    (b"O", 8),
                    (b"F", 9),
                    (b"P", 15),
                    (b"S", 16),
                    (b"I", 53),
                ];
                for &(sym, num) in table {
                    if rest.starts_with(sym) {
                        self.pos += sym.len();
                        return Ok(AtomExpr::ElemAliphatic(num));
                    }
                }
                Err(self.err("unknown element"))
            }
            Some(c) if c.is_ascii_lowercase() => {
                let n = match c {
                    b'b' => 5,
                    b'c' => 6,
                    b'n' => 7,
                    b'o' => 8,
                    b'p' => 15,
                    b's' => 16,
                    _ => return Err(self.err("unknown aromatic element")),
                };
                self.pos += 1;
                Ok(AtomExpr::ElemAromatic(n))
            }
            _ => Err(self.err("expected atom primitive")),
        }
    }
}

// ---------------------------------------------------------------------------
// マッチャ
// ---------------------------------------------------------------------------

/// 分子側の前処理済みビュー (SMARTS 評価に必要な原子性質)。
pub struct MolView {
    pub n: usize,
    atomic_num: Vec<u8>,
    aromatic: Vec<bool>,
    charge: Vec<i8>,
    degree: Vec<u8>,
    total_h: Vec<u8>,
    ring_bond_count: Vec<u8>,
    ring_count: Vec<u8>,
    ring_sizes: Vec<Vec<u8>>,
    hyb: Vec<u8>,                     // 1/2/3
    adj: Vec<Vec<(usize, u8, bool)>>, // (相手, 次数クラス 0単/1芳香/2二重/3三重, 環内か)
}

/// RDKit 互換の混成 (SMARTS `^n` 用): 立体数 = 結合数 + 孤立電子対。
/// チオカルボニル S は sp2、スルホキシド/スルホン S は sp3 になる
/// (配座生成用の perceive_hybridization とは目的が異なる別実装)。
fn rdkit_like_hybridization(g: &MoleculeGraph, i: usize) -> u8 {
    let a = &g.atoms[i];
    if a.is_aromatic {
        return 2;
    }
    let outer: i32 = match a.atomic_num {
        1 => 1,
        5 => 3,
        6 | 14 => 4,
        7 | 15 | 33 | 51 | 83 => 5,
        8 | 16 | 34 | 52 => 6,
        9 | 17 | 35 | 53 => 7,
        _ => return 3,
    };
    let mut order_sum = 0.0f64;
    let mut degree = 0i32;
    for b in &g.bonds {
        if b.begin_idx == i || b.end_idx == i {
            order_sum += b.bond_order;
            degree += 1;
        }
    }
    let lone_pairs = ((outer - a.formal_charge as i32 - order_sum.round() as i32) / 2).max(0);
    match degree + lone_pairs {
        0..=2 => 1, // sp
        3 => 2,     // sp2
        _ => 3,     // sp3 以上
    }
}

fn order_class(order: f64) -> u8 {
    if order == 1.5 {
        1
    } else if order == 2.0 {
        2
    } else if order == 3.0 {
        3
    } else {
        0
    }
}

impl MolView {
    pub fn build(g: &MoleculeGraph) -> MolView {
        let n = g.atoms.len();
        // 環結合の判定 (対称化 SSSR の辺)
        let mut ring_bond: HashMap<(usize, usize), bool> = HashMap::new();
        let mut ring_count = vec![0u8; n];
        let mut ring_sizes: Vec<Vec<u8>> = vec![Vec::new(); n];
        for ring in &g.ring_atom_sets {
            for (t, &a) in ring.iter().enumerate() {
                ring_count[a] = ring_count[a].saturating_add(1);
                if !ring_sizes[a].contains(&(ring.len() as u8)) {
                    ring_sizes[a].push(ring.len() as u8);
                }
                let b = ring[(t + 1) % ring.len()];
                ring_bond.insert((a.min(b), a.max(b)), true);
            }
        }
        let mut adj: Vec<Vec<(usize, u8, bool)>> = vec![Vec::new(); n];
        let mut ring_bond_count = vec![0u8; n];
        for b in &g.bonds {
            let (i, j) = (b.begin_idx, b.end_idx);
            let in_ring = ring_bond.contains_key(&(i.min(j), i.max(j)));
            let oc = order_class(b.bond_order);
            adj[i].push((j, oc, in_ring));
            adj[j].push((i, oc, in_ring));
            if in_ring {
                ring_bond_count[i] += 1;
                ring_bond_count[j] += 1;
            }
        }

        let total_h: Vec<u8> = (0..n)
            .map(|i| {
                g.adjacency[i]
                    .iter()
                    .filter(|&&x| g.atoms[x].symbol == "H")
                    .count() as u8
            })
            .collect();
        MolView {
            n,
            atomic_num: (0..n).map(|i| g.atoms[i].atomic_num).collect(),
            aromatic: (0..n).map(|i| g.atoms[i].is_aromatic).collect(),
            charge: (0..n).map(|i| g.atoms[i].formal_charge).collect(),
            degree: (0..n).map(|i| g.adjacency[i].len() as u8).collect(),
            total_h,
            ring_bond_count,
            ring_count,
            ring_sizes,
            hyb: (0..n).map(|i| rdkit_like_hybridization(g, i)).collect(),
            adj,
        }
    }
}

fn atom_matches(e: &AtomExpr, view: &MolView, t: usize, recursives: &[SmartsPattern]) -> bool {
    match e {
        AtomExpr::Any => true,
        AtomExpr::Aromatic => view.aromatic[t],
        AtomExpr::Aliphatic => !view.aromatic[t],
        AtomExpr::AtomicNum(n) => view.atomic_num[t] == *n,
        AtomExpr::ElemAliphatic(n) => view.atomic_num[t] == *n && !view.aromatic[t],
        AtomExpr::ElemAromatic(n) => view.atomic_num[t] == *n && view.aromatic[t],
        AtomExpr::TotalConn(n) => view.degree[t] == *n,
        AtomExpr::TotalH(n) => view.total_h[t] == *n,
        AtomExpr::RingBondCount(n) => view.ring_bond_count[t] == *n,
        AtomExpr::InRingSize(n) => view.ring_sizes[t].contains(n),
        AtomExpr::InRing => view.ring_count[t] > 0,
        AtomExpr::RingCount(n) => view.ring_count[t] == *n,
        AtomExpr::Hybridization(n) => view.hyb[t] == *n,
        AtomExpr::Charge(c) => view.charge[t] == *c,
        AtomExpr::Recursive(id) => rooted_match(view, &recursives[*id], recursives, t),
        AtomExpr::Not(inner) => !atom_matches(inner, view, t, recursives),
        AtomExpr::And(v) => v.iter().all(|x| atom_matches(x, view, t, recursives)),
        AtomExpr::Or(v) => v.iter().any(|x| atom_matches(x, view, t, recursives)),
    }
}

fn bond_matches(e: &BondExpr, oc: u8, in_ring: bool) -> bool {
    match e {
        BondExpr::Single => oc == 0,
        BondExpr::Double => oc == 2,
        BondExpr::Triple => oc == 3,
        BondExpr::AromaticB => oc == 1,
        BondExpr::Any => true,
        BondExpr::Ring => in_ring,
        BondExpr::Default => oc == 0 || oc == 1,
        BondExpr::Not(inner) => !bond_matches(inner, oc, in_ring),
        BondExpr::And(v) => v.iter().all(|x| bond_matches(x, oc, in_ring)),
        BondExpr::Or(v) => v.iter().any(|x| bond_matches(x, oc, in_ring)),
    }
}

/// VF2 バックトラッキングの共通部。root が Some なら
/// クエリ原子 0 をその原子に固定する (再帰 SMARTS 用)。
fn backtrack_all(
    view: &MolView,
    pat: &SmartsPattern,
    recursives: &[SmartsPattern],
    root: Option<usize>,
    first_only: bool,
    results: &mut Vec<Vec<usize>>,
) {
    let nq = pat.atoms.len();
    // 探索順: 既出隣接優先
    let mut order = Vec::with_capacity(nq);
    let mut placed = vec![false; nq];
    order.push(0);
    placed[0] = true;
    while order.len() < nq {
        let next = (0..nq)
            .filter(|&i| !placed[i])
            .find(|&i| pat.adj[i].iter().any(|&(j, _)| placed[j]))
            .or_else(|| (0..nq).find(|&i| !placed[i]))
            .expect("unplaced exists");
        placed[next] = true;
        order.push(next);
    }

    const MAX_MATCHES: usize = 100_000;
    let mut mapping = vec![usize::MAX; nq];
    let mut used = vec![false; view.n];

    #[allow(clippy::too_many_arguments)]
    fn rec(
        depth: usize,
        order: &[usize],
        view: &MolView,
        pat: &SmartsPattern,
        recursives: &[SmartsPattern],
        root: Option<usize>,
        first_only: bool,
        mapping: &mut [usize],
        used: &mut [bool],
        results: &mut Vec<Vec<usize>>,
    ) -> bool {
        if results.len() >= MAX_MATCHES {
            return true;
        }
        if depth == order.len() {
            results.push(mapping.to_vec());
            return first_only;
        }
        let qi = order[depth];
        let candidates: Vec<usize> = if qi == 0 {
            match root {
                Some(r) => vec![r],
                None => (0..view.n).collect(),
            }
        } else {
            let anchor = pat.adj[qi]
                .iter()
                .find(|&&(j, _)| mapping[j] != usize::MAX)
                .map(|&(j, _)| mapping[j]);
            match anchor {
                Some(ta) => view.adj[ta].iter().map(|&(v, _, _)| v).collect(),
                None => (0..view.n).collect(),
            }
        };
        'cand: for ti in candidates {
            if used[ti] || !atom_matches(&pat.atoms[qi], view, ti, recursives) {
                continue;
            }
            for &(qj, ei) in &pat.adj[qi] {
                let tj = mapping[qj];
                if tj == usize::MAX {
                    continue;
                }
                let Some(&(_, oc, in_ring)) = view.adj[ti].iter().find(|&&(v, _, _)| v == tj)
                else {
                    continue 'cand;
                };
                if !bond_matches(&pat.bonds[ei].2, oc, in_ring) {
                    continue 'cand;
                }
            }
            mapping[qi] = ti;
            used[ti] = true;
            let stop = rec(
                depth + 1,
                order,
                view,
                pat,
                recursives,
                root,
                first_only,
                mapping,
                used,
                results,
            );
            mapping[qi] = usize::MAX;
            used[ti] = false;
            if stop {
                return true;
            }
        }
        false
    }
    rec(
        0,
        &order,
        view,
        pat,
        recursives,
        root,
        first_only,
        &mut mapping,
        &mut used,
        results,
    );
}

/// 再帰 SMARTS の根付き判定。
fn rooted_match(
    view: &MolView,
    pat: &SmartsPattern,
    recursives: &[SmartsPattern],
    root: usize,
) -> bool {
    let mut results = Vec::new();
    backtrack_all(view, pat, recursives, Some(root), true, &mut results);
    !results.is_empty()
}

/// 全マッチ列挙 (uniquify なし)。タプルはクエリ原子順。
pub fn smarts_matches(view: &MolView, pat: &SmartsPattern) -> Vec<Vec<usize>> {
    let mut results = Vec::new();
    backtrack_all(view, pat, &pat.recursives, None, false, &mut results);
    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::build_molecule_graph;

    fn matches(target: &str, smarts: &str) -> Vec<Vec<usize>> {
        let g = build_molecule_graph(target).expect("target");
        let view = MolView::build(&g);
        let pat = parse_smarts(smarts).expect("smarts");
        smarts_matches(&view, &pat)
    }

    #[test]
    fn primitives() {
        // X4H2: メチレン炭素
        assert_eq!(matches("CCC", "[CX4H2]").len(), 1);
        // H 数
        assert_eq!(matches("CC(C)C", "[CX4H1]").len(), 1);
        // 環サイズ r
        assert_eq!(matches("C1CC1CCC1CCCC1", "[C&r3]").len(), 3);
        assert_eq!(matches("C1CC1CCC1CCCC1", "[C&r5]").len(), 5);
        // 環結合数 x0
        assert_eq!(matches("c1ccccc1C", "[CX4x0]").len(), 1);
        // 混成 ^2
        assert_eq!(matches("CC=CC", "[C^2]").len(), 2);
        // 電荷
        assert_eq!(matches("CC(=O)[O-]", "[O-]").len(), 1);
        // 論理: OR と NOT
        assert_eq!(matches("CCO", "[C,O]").len(), 3);
        assert_eq!(matches("CCO", "[!#1;!O]").len(), 2);
        // 芳香族 a / 脂肪族 A
        assert_eq!(matches("c1ccccc1C", "[a]").len(), 6);
        assert_eq!(matches("c1ccccc1C", "[A]").len(), 1 + 5 + 3); // C + H×8
    }

    #[test]
    fn bonds_and_ring_bonds() {
        // !@;- : 環外単結合
        assert_eq!(matches("c1ccccc1c1ccccc1", "[c]!@;-[c]").len(), 2); // 両向き
        assert_eq!(matches("c1ccccc1", "[c]!@;-[c]").len(), 0);
        // ~ 任意
        assert_eq!(matches("C=C", "[C]~[C]").len(), 2);
        // @ 環内
        assert_eq!(matches("C1CC1C", "[C]@[C]").len(), 6);
    }

    #[test]
    fn recursive_smarts() {
        // カルボニル炭素に隣接する O
        let m = matches("CC(=O)OC", "[$(C=O)][O]");
        assert_eq!(m.len(), 1);
        // ネスト再帰
        let m = matches("NC(=O)NC", "[$([C](=O)[$([NX3H2])])]");
        assert_eq!(m.len(), 1);
        // 否定つき再帰
        let m = matches("CC(=O)C", "[CX4;!$(C=O)]");
        assert_eq!(m.len(), 2);
    }

    #[test]
    fn atom_maps_and_h_atom() {
        let pat = parse_smarts("[O:1]=[C:2]!@;-[NX3H1:3][!#1:4]").expect("parses");
        assert_eq!(pat.atom_maps[0], Some(1));
        assert_eq!(pat.atom_maps[3], Some(4));
        // N-メチルアセトアミド: O=C-N-C
        let m = matches("CC(=O)NC", "[O:1]=[C:2]!@;-[NX3H1:3][!#1:4]");
        assert_eq!(m.len(), 1);
        // [H:4] は明示 H にマッチ
        let m = matches("CC(=O)NC", "[O:1]=[C:2]!@;-[NX3H1:3][H:4]");
        assert_eq!(m.len(), 1);
    }
}
