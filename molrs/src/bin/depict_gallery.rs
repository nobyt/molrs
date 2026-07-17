//! D13 開発用: stdin の SMILES (1 行 1 分子) を 2D 描画し、
//! 目視レビュー用の単一 HTML ギャラリーを stdout に出力する。
//!
//! 使い方:
//!   jq -r .smiles ../corpus/corpus.jsonl | awk 'NR%50==1' \
//!     | cargo run --release --bin depict_gallery > gallery.html

use std::io::{BufRead, Write};

use molrs::depict::{compute_coords_2d, to_svg, LayoutParams, Style};

fn main() {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    let style = Style::acs_1996();
    let params = LayoutParams::default();

    writeln!(
        out,
        "<!doctype html><meta charset=\"utf-8\"><title>molrs depict gallery</title>\
         <style>body{{font-family:sans-serif;display:flex;flex-wrap:wrap;gap:12px}}\
         figure{{margin:0;border:1px solid #ccc;padding:8px;max-width:280px}}\
         figcaption{{font-size:11px;word-break:break-all;color:#444}}</style>"
    )
    .unwrap();

    for line in stdin.lock().lines() {
        let line = line.expect("stdin");
        let smiles = line.trim();
        if smiles.is_empty() {
            continue;
        }
        let body = match molrs::graph::build_molecule_graph(smiles) {
            Ok(g) => match compute_coords_2d(&g, &params) {
                Ok(c) => to_svg(&g, &c, &style),
                Err(e) => format!("<p>layout error: {e}</p>"),
            },
            Err(e) => format!("<p>parse error: {e}</p>"),
        };
        let esc = smiles
            .replace('&', "&amp;")
            .replace('<', "&lt;")
            .replace('>', "&gt;");
        writeln!(out, "<figure>{body}<figcaption>{esc}</figcaption></figure>").unwrap();
    }
}
