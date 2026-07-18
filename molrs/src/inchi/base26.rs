//! InChIKey の base-26 符号化とキー組み立て (I5)。
//!
//! IUPAC 公式 InChI の `ikey_base26.c` を移植。標準 InChI 文字列を
//! 「major (式・c・h・q)」と「minor (立体・同位体)」に分け、それぞれの
//! SHA-256 から 14 文字・8 文字の base-26 ブロックを作る。プロトン化 (/p)
//! は末尾 1 文字に符号化する。
//!
//! t26 テーブル (14bit → 3 文字) は先頭が 'E' の 676 個と `TAA..TTV` の
//! 516 個を除いた 16384 個の有効トリプレット列。有効インデックス h から
//! プレーン base-26 インデックス p への写像は
//! `p = h + (h>=2704 ? 676 : 0) + (h>=12168 ? 516 : 0)`。

use super::sha256::sha256;

const A: &[u8; 26] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ";

/// 14bit 値 → 3 文字トリプレット (無効トリプレットを飛ばした列)。
fn t26(h: u32) -> [u8; 3] {
    let mut p = h as usize;
    if h >= 2704 {
        p += 676;
    }
    if h >= 12168 {
        p += 516;
    }
    [A[p / 676], A[(p / 26) % 26], A[p % 26]]
}

/// 9bit 値 → 2 文字ダブレット (プレーン base-26)。
fn d26(h: u32) -> [u8; 2] {
    let h = h as usize;
    [A[h / 26], A[h % 26]]
}

fn u(b: u8) -> u32 {
    b as u32
}

/// major ブロック用: 先頭 65bit から 14 文字。
fn block1(a: &[u8; 32]) -> String {
    let t1 = t26(u(a[0]) | (u(a[1]) & 0x3f) << 8);
    let t2 = t26(((u(a[1]) & 0xc0) | u(a[2]) << 8 | (u(a[3]) & 0x0f) << 16) >> 6);
    let t3 = t26(((u(a[3]) & 0xf0) | u(a[4]) << 8 | (u(a[5]) & 0x03) << 16) >> 4);
    let t4 = t26(((u(a[5]) & 0xfc) | u(a[6]) << 8) >> 2);
    let du = d26(u(a[7]) | (u(a[8]) & 0x01) << 8);
    let mut s = String::with_capacity(14);
    for part in [&t1[..], &t2, &t3, &t4, &du] {
        s.push_str(std::str::from_utf8(part).unwrap());
    }
    s
}

/// minor ブロック用: 先頭 37bit から 8 文字。
fn block2(a: &[u8; 32]) -> String {
    let t1 = t26(u(a[0]) | (u(a[1]) & 0x3f) << 8);
    let t2 = t26(((u(a[1]) & 0xc0) | u(a[2]) << 8 | (u(a[3]) & 0x0f) << 16) >> 6);
    let du = d26(((u(a[3]) & 0xf0) | (u(a[4]) & 0x1f) << 8) >> 4);
    let mut s = String::with_capacity(8);
    for part in [&t1[..], &t2, &du] {
        s.push_str(std::str::from_utf8(part).unwrap());
    }
    s
}

/// 標準 InChI 文字列 → InChIKey (`AAAAAAAAAAAAAA-BBBBBBBBSA-P`)。
///
/// 非標準 InChI (`InChI=1/`) や v2 の minor 文字列構成差には未対応
/// (v1 の生成 InChI は立体・同位体を持たないため minor は空 = `UHFFFAOY`)。
pub fn inchi_key_from_string(inchi: &str) -> String {
    // "InChI=1S/" を剥がす (標準)。非標準は "InChI=1/"。
    let body = inchi
        .strip_prefix("InChI=1S/")
        .or_else(|| inchi.strip_prefix("InChI=1/"))
        .unwrap_or(inchi);

    let mut formula = "";
    let mut major_layers: Vec<&str> = Vec::new();
    let mut minor_layers: Vec<&str> = Vec::new();
    let mut proton: i32 = 0;
    for (i, seg) in body.split('/').enumerate() {
        if i == 0 {
            formula = seg;
            continue;
        }
        let Some(tag) = seg.chars().next() else {
            continue;
        };
        match tag {
            'c' | 'h' | 'q' => major_layers.push(seg),
            'b' | 't' | 'm' | 's' | 'i' | 'f' => minor_layers.push(seg),
            'p' => {
                proton = seg[1..].parse().unwrap_or(0);
            }
            _ => {}
        }
    }

    // major 文字列 = 式 + 各層 (先頭スラッシュなし、'/' 区切り)
    let mut major = String::from(formula);
    for l in &major_layers {
        major.push('/');
        major.push_str(l);
    }
    let b1 = block1(&sha256(major.as_bytes()));

    // minor 文字列 = 立体・同位体層 (先頭 '/' 込み) を 2 重連結してハッシュ
    // (公式 ikey.c: strcpy(sminor+slen, sminor))。空なら sha256("") → "UHFFFAOY"。
    let b2 = if minor_layers.is_empty() {
        block2(&sha256(b""))
    } else {
        let mut minor = String::new();
        for l in &minor_layers {
            minor.push('/');
            minor.push_str(l);
        }
        let doubled = format!("{minor}{minor}");
        block2(&sha256(doubled.as_bytes()))
    };

    // フラグ: 標準 InChI = 'S'、バージョン = 'A'
    let flags = "SA";
    // プロトン化文字: 0 → 'N'、±k → 'N'±k
    let proton_char = proton_to_char(proton);

    format!("{b1}-{b2}{flags}-{proton_char}")
}

/// プロトン数 → 文字。0='N'、範囲外は '*'。
fn proton_to_char(p: i32) -> char {
    if p == 0 {
        return 'N';
    }
    let c = 'N' as i32 + p;
    if ('A' as i32..='Z' as i32).contains(&c) {
        (c as u8) as char
    } else {
        '*'
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_keys() {
        // RDKit InchiToInchiKey と一致すべき既知 InChI
        assert_eq!(
            inchi_key_from_string("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"),
            "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        );
        assert_eq!(
            inchi_key_from_string("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"),
            "UHOVQNZJYSORNB-UHFFFAOYSA-N"
        );
        assert_eq!(
            inchi_key_from_string("InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)"),
            "QTBSBXVTEAMEQO-UHFFFAOYSA-N"
        );
    }

    #[test]
    fn protonation() {
        // [NH4+] → /p+1 → 末尾 'O'
        assert_eq!(
            inchi_key_from_string("InChI=1S/H3N/h1H3/p+1"),
            "QGZKDVFQNNGYKY-UHFFFAOYSA-O"
        );
        // カルボキシラート → /p-1 → 末尾 'M'
        assert_eq!(
            inchi_key_from_string("InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)/p-1"),
            "QTBSBXVTEAMEQO-UHFFFAOYSA-M"
        );
    }

    #[test]
    fn empty_minor_is_uhfffaoy() {
        assert_eq!(block2(&sha256(b"")), "UHFFFAOY");
    }
}
