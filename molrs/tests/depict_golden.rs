//! D3 完了条件ゲート: 代表分子の SVG golden スナップショット (バイト一致 =
//! 決定性の担保) + well-formed チェック。
//!
//! golden の更新: `UPDATE_DEPICT_GOLDEN=1 cargo test --test depict_golden`

use std::path::PathBuf;

use molrs::depict::{depict_svg, Style};

const GOLDEN: &[(&str, &str)] = &[
    // 無環 (D3)
    ("ethanol", "CCO"),
    ("isobutane", "CC(C)C"),
    ("e2butene", "C/C=C/C"),
    ("acetonitrile", "CC#N"),
    ("acetic_acid", "CC(=O)O"),
    // 環系 (D4-D6)
    ("benzene", "c1ccccc1"),
    ("toluene", "Cc1ccccc1"),
    ("pyridine", "c1ccncc1"),
    ("naphthalene", "c1ccc2ccccc2c1"),
    ("indole", "c1ccc2[nH]ccc2c1"),
    ("norbornane", "C1CC2CCC1C2"),
    ("spiro45decane", "C1CCC2(CC1)CCCC2"),
];

fn golden_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

/// 依存なしの簡易 XML well-formed チェック: タグの開閉対応と属性引用。
fn assert_well_formed(svg: &str) {
    let mut stack: Vec<String> = Vec::new();
    let bytes = svg.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] != b'<' {
            i += 1;
            continue;
        }
        let end = svg[i..].find('>').map(|e| i + e).expect("unclosed '<'");
        let tag = &svg[i + 1..end];
        // 引用符の対応 (タグ内の " が偶数個)
        assert!(
            tag.matches('"').count().is_multiple_of(2),
            "unbalanced quotes in <{tag}>"
        );
        if let Some(name) = tag.strip_prefix('/') {
            let open = stack.pop().unwrap_or_else(|| panic!("stray </{name}>"));
            assert_eq!(open, name, "mismatched close tag");
        } else if !tag.ends_with('/') && !tag.starts_with('?') && !tag.starts_with('!') {
            let name: String = tag
                .chars()
                .take_while(|c| !c.is_whitespace() && *c != '>')
                .collect();
            stack.push(name);
        }
        i = end + 1;
    }
    assert!(stack.is_empty(), "unclosed tags: {stack:?}");
}

#[test]
fn golden_svg_snapshots() {
    let update = std::env::var("UPDATE_DEPICT_GOLDEN").is_ok();
    let dir = golden_dir();
    if update {
        std::fs::create_dir_all(&dir).expect("mkdir golden");
    }
    let style = Style::acs_1996();
    let mut failures = Vec::new();
    for (name, smiles) in GOLDEN {
        let svg = depict_svg(smiles, &style)
            .unwrap_or_else(|e| panic!("depict failed for {name} ({smiles}): {e}"));
        assert_well_formed(&svg);
        let path = dir.join(format!("{name}.svg"));
        if update {
            std::fs::write(&path, &svg).expect("write golden");
            continue;
        }
        let expected = std::fs::read_to_string(&path).unwrap_or_else(|e| {
            panic!(
                "missing golden {}: {e} (run with UPDATE_DEPICT_GOLDEN=1)",
                path.display()
            )
        });
        if svg != expected {
            failures.push(format!("{name} ({smiles}) differs from golden"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn well_formed_checker_rejects_bad_xml() {
    let result = std::panic::catch_unwind(|| assert_well_formed("<svg><line></svg>"));
    assert!(result.is_err());
}
