//! InChI 差分ゲート (RUST_INCHI_PLAN.md §検証)。
//!
//! フィクスチャ corpus/inchi_dump.jsonl.gz (smiles2iupac tools/dump_inchi.py で
//! RDKit から採取) に対して層ごとに一致を検査する。実装の進行に応じて検査
//! 項目を増やす:
//! - I2: 式層 (Hill) の一致率
//! - I3: 正準番号 (AuxInfo /N:) の一致率
//! - I4/I5: フル InChI 文字列・InChIKey の一致率 (v1 適用範囲で 100%)

use std::io::Read;
use std::path::PathBuf;

use flate2::read::GzDecoder;
use molrs::graph::build_molecule_graph;

struct Record {
    smiles: String,
    formula: String,
    numbering: Vec<Vec<usize>>,
    inchi: String,
    key: String,
}

fn load_fixture() -> Vec<Record> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../corpus/inchi_dump.jsonl.gz");
    let file =
        std::fs::File::open(&path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    let mut text = String::new();
    GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("gunzip");
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let v: serde_json::Value = serde_json::from_str(l).expect("json");
            let numbering = v["n"]
                .as_array()
                .map(|comps| {
                    comps
                        .iter()
                        .map(|c| {
                            c.as_array()
                                .unwrap()
                                .iter()
                                .map(|x| x.as_u64().unwrap() as usize)
                                .collect()
                        })
                        .collect()
                })
                .unwrap_or_default();
            Record {
                smiles: v["s"].as_str().unwrap().to_string(),
                formula: v["formula"].as_str().unwrap_or("").to_string(),
                numbering,
                inchi: v["inchi"].as_str().unwrap_or("").to_string(),
                key: v["key"].as_str().unwrap_or("").to_string(),
            }
        })
        .collect()
}

#[test]
fn formula_layer_matches_rdkit() {
    let recs = load_fixture();
    let mut n = 0usize;
    let mut ok = 0usize;
    let mut mism: Vec<String> = Vec::new();
    for r in &recs {
        if r.formula.is_empty() {
            continue;
        }
        let Ok(g) = build_molecule_graph(&r.smiles) else {
            continue;
        };
        n += 1;
        let got = molrs::inchi::formula(&g);
        if got == r.formula {
            ok += 1;
        } else if mism.len() < 25 {
            mism.push(format!("{}: got {got}, want {}", r.smiles, r.formula));
        }
    }
    let rate = ok as f64 / n as f64;
    println!("formula layer: {ok}/{n} match ({:.2}%)", rate * 100.0);
    for m in &mism {
        println!("  MISMATCH {m}");
    }
    // v1: 電荷正規化 (プロトン移動) を要さない分子は一致するはず。
    // 正規化差の分子があるため閾値は段階的に引き上げる (I4 で normalize 実装後)。
    assert!(rate >= 0.90, "formula match rate {rate:.4} < 0.90");
}

#[test]
fn canonical_numbering_matches_auxinfo() {
    let recs = load_fixture();
    let mut n = 0usize;
    let mut ok = 0usize;
    let mut mism: Vec<String> = Vec::new();
    for r in &recs {
        if r.numbering.is_empty() {
            continue;
        }
        let Ok(g) = build_molecule_graph(&r.smiles) else {
            continue;
        };
        n += 1;
        let got: Vec<Vec<usize>> = molrs::inchi::number::canonical_numbering(&g)
            .iter()
            .map(|comp| comp.iter().map(|&i| i + 1).collect())
            .collect();
        if got == r.numbering {
            ok += 1;
        } else if mism.len() < 20 {
            mism.push(format!("{}: got {got:?}, want {:?}", r.smiles, r.numbering));
        }
    }
    let rate = ok as f64 / n as f64;
    println!("canonical numbering: {ok}/{n} match ({:.2}%)", rate * 100.0);
    for m in &mism {
        println!("  MISMATCH {m}");
    }
    // I3 到達点 = 99.29% (7400/7453)。残差の内訳:
    // - 有機金属 (金属-C 切断で多成分化) — I4 normalize で対応
    // - イソシアニド等の電荷正規化差 — I4
    // - 環内 N 互変異性 (ベンズイミダゾリン等) — v2
    // - 立体依存タイブレーク (E/Z 対称分子) — 骨格層の文字列には影響しない
    assert!(rate >= 0.985, "numbering match rate {rate:.4} < 0.985");
}

#[test]
fn full_inchi_matches_rdkit_where_produced() {
    // v1 が文字列を生成する範囲 (中性・単一成分) でフル InChI が一致すること。
    // カバレッジ (生成できた割合) も報告する。
    let recs = load_fixture();
    let mut produced = 0usize;
    let mut ok = 0usize;
    let mut total_neutral_single = 0usize;
    let mut mism: Vec<String> = Vec::new();
    for r in &recs {
        if r.inchi.is_empty() {
            continue;
        }
        let Ok(g) = build_molecule_graph(&r.smiles) else {
            continue;
        };
        // v1 適用範囲: 同位体層 (/i) を含まない中性・単一成分 (立体は対応済み)
        let has_iso = r.inchi.contains("/i");
        let charged = g.atoms.iter().any(|a| a.formal_charge != 0);
        let multi = r.inchi[9..].split('/').next().unwrap().contains('.');
        if !charged && !multi && !has_iso {
            total_neutral_single += 1;
        }
        if has_iso {
            continue; // v2
        }
        if let Ok(got) = molrs::inchi::to_inchi(&g) {
            produced += 1;
            if got == r.inchi {
                ok += 1;
            } else if mism.len() < 30 {
                mism.push(format!("{} | got {got} | want {}", r.smiles, r.inchi));
            }
        }
    }
    let acc = ok as f64 / produced.max(1) as f64;
    let cov = produced as f64 / total_neutral_single.max(1) as f64;
    println!(
        "full InChI: produced {produced}, {ok} exact ({:.2}%); coverage of neutral-single {produced}/{total_neutral_single} ({:.2}%)",
        acc * 100.0,
        cov * 100.0
    );
    for m in &mism {
        println!("  MISMATCH {m}");
    }
    // I14 到達点 = 93.33% (v1 適用範囲 = 中性・単一成分・立体/同位体なし)。
    // 89.33% からの追加修正 (`seed_groups`/`mobile_groups` の細部):
    // (1) 種形成時の「この中心に局所的な二重結合受容体端点を最低 1 つ」は
    //     維持しつつ、局所的な H 要件は撤廃 (真の供与体は後からブリッジで
    //     合流すればよい)。(2) 複数の中心がヘテロ端点を共有してもよい
    //     (`used` による早期排除を撤廃、union-find が正しく統合する)。
    //     (3) 三重結合や累積二重結合の中心 (中心自身が単結合の自由な端点を
    //     1 つも持たない場合) を誤って星型と判定しないガードを追加。
    //     (4) 縮環の共有原子は環ごとに実際の二重結合が「別の」環に割り当た
    //     ることがあるため、環内限定の仮想マッチ (同じ環で隣接する 2 つの
    //     「その環では二重結合を持たない」共有原子どうしのペア付け) で
    //     環の中心寄りのヘテロ原子も正しく拾う。(5) ブリッジ探索の辺は
    //     「環外側がヘテロ原子」の場合のみ採用 (無関係な環外炭素への
    //     誤った抜け道を防ぐ)。
    // I16 (2026-08-03): 実機 inchi-1 の BalancedNetworkSearch を計装調査し、
    //     (a) 原子ごとの価数スラック 0 (橋頭ヘテロ原子など) はどの辺にも
    //     参加できないという実 InChI の制約を `has_search_slack` として
    //     追加 (理論的に正しいが本コーパスでは中立、退行なし)、
    //     (b) 縮環越しブリッジの成立条件 (共有原子が両環のヘテロ端点に
    //     直接隣接) は既に正しく実装済みと実機で確認 (追加実装は不要)。
    //     既知の残る限界: 縮環系の橋頭ヘテロ原子 (例: インドリジン型の
    //     3 配位 N) 自身の辺を経由しない、単一環内だけの誤った橋渡し
    //     (例: 環外フェノール性 OH → 環内の遠い N) が 1 ケース種別として
    //     未解決 (2 回の実機計装調査でも根本メカニズム未特定、詳細は
    //     RUST_INCHI_PLAN.md I16 参照)。残る不一致の主因は一部の縮合環
    //     c 層直列化と立体絡みのケース。閾値は退行検知用。
    // I17 (2026-08-04): 電荷分離 (zwitterion) で固定された負電荷 (ニトロ
    //     N+(=O)O-・芳香族 N-オキシド n+ - O-) を可動プロトンとして数える
    //     誤りを修正。純粋なニトロは h 層に (H,...) 群を出さなくなった
    //     (2 つの O は正準番号付けの等価化のため群メンバーとしては残す)。
    //     加えてビアリール連結結合 (別々の芳香環を単結合でつなぐ、共有環
    //     なし) をブリッジ探索の辺から除外 (どの Kekule 構造でも単結合の
    //     ため互変異性経路にならない)。93.33% → 93.51%。
    // I18 (2026-08-04): 残差を層別に分類して以下を修正。(a) c 層の枝順序を
    //     「部分木の原子数 + 部分木内の環閉合数」の重みで決定 (実 InChI は
    //     環閉合も 1 項目として数える)。93.51% → 95.91%。(b) 立体源性だが
    //     未定義の C=N (イミン/ヒドラゾン) を /b 層に `?` 併記
    //     (定義済み立体が別にある分子のみ、可動 H 群メンバーは除く)。
    //     95.91% → 96.05%。(c) 縮環孤児のペアリングを 2 個限定から環に
    //     沿った貪欲ペアリングに一般化 (アクリドン型の中央環)。96.05% →
    //     96.12%。(d) ブリッジ探索の辺に「環内カルボニル型受容体 C」を
    //     追加 (イサチン/フタルイミド型が両カルボニル O に到達)。96.12% →
    //     96.18%。残る主因は縮環系の正準番号付け直列化と一部の可動 H。
    // I19 §3.1 (2026-08-04): リン中心 (P=O) の一級アミド端点を許可。硫黄の
    //     スルホンアミド判定 (S=O が 2 個必要) と異なり、リンは P=O が
    //     1 個でも一級 NH2 が可動になる (リン酸トリアミド等)。96.18% →
    //     96.20%。
    // I19 §3.2 (2026-08-04): O/S だけの酸系対に対する N 除外規則 (カルバミン
    //     酸のような置換 N を除外) が、末端の一級 NH2 (heavy_deg==1) まで
    //     過剰に除外していたのを修正。ジチオカルバミン酸・スルファミン酸の
    //     NH2 は酸対と一緒に可動になる。96.20% → 96.21%。
    // I19 §3.3 (2026-08-04): 縮環の共有原子 (頂点分割で分身を持つ) が、
    //     環ごとに独立した別々の探索から「たまたま両方とも通過する」
    //     ケースで、その原子自身がヘテロ原子でなくても union-find が
    //     原子 ID ベースで併合してしまい、本来無関係な 2 つの互変異性系
    //     (例: 縮環の一方の環の環外 OH↔環内 N と、もう一方の環内の
    //     独立した N-N-H) が誤って 1 群に統合されていた。到達原子が
    //     ヘテロのときだけ union するよう修正 (頂点分割は探索グラフでは
    //     正しく環を分離しているので、非ヘテロの経由点まで union する
    //     必要はない)。広範囲に影響するバグで大幅改善。96.21% → 97.41%
    //     (+89)。
    // I19 §3.2 追加修正 (2026-08-04): §3.2 の「末端一級 NH2 は酸対から
    //     除外しない」規則が広すぎ、無置換のカルバミン酸 `NC(=O)O`
    //     (NH2 も一級) まで誤って可動化していた (want は N 固定・O,O の
    //     みが酸対)。中心が炭素かつ酸対が純粋に酸素だけ (アミド性、
    //     ジチオカルバミン酸やスルファミン酸とは異なる) の場合に限り、
    //     一級 NH2 でも除外を維持するよう精緻化。97.41% → 97.42%。
    assert!(acc >= 0.9741, "full InChI accuracy {acc:.4} < 0.9741");
}

#[test]
fn inchi_key_from_string_matches_rdkit() {
    // キー機構の独立検証: RDKit の InChI 文字列 → キーが RDKit の
    // InchiToInchiKey と一致すること (立体/同位体を含む minor は v2 のため除外)。
    let recs = load_fixture();
    let mut n = 0usize;
    let mut ok = 0usize;
    let mut mism: Vec<String> = Vec::new();
    for r in &recs {
        if r.inchi.is_empty() || r.key.is_empty() {
            continue;
        }
        // minor 層 (立体/同位体) を含む InChI は minor 文字列構成が未対応 (v2)
        let has_minor = ["/b", "/t", "/m", "/s", "/i"]
            .iter()
            .any(|t| r.inchi.contains(t));
        if has_minor {
            continue;
        }
        n += 1;
        let got = molrs::inchi::inchi_key_from_string(&r.inchi);
        if got == r.key {
            ok += 1;
        } else if mism.len() < 20 {
            mism.push(format!("{} | got {got} | want {}", r.inchi, r.key));
        }
    }
    println!("inchi_key_from_string: {ok}/{n} match",);
    for m in &mism {
        println!("  MISMATCH {m}");
    }
    // 非立体 InChI ではキー機構は完全一致すべき
    assert_eq!(ok, n, "key machinery must be exact on non-stereo InChI");
}

#[test]
fn to_inchi_key_matches_rdkit_where_produced() {
    // エンドツーエンド: to_inchi_key が v1 範囲で RDKit MolToInchiKey と一致。
    let recs = load_fixture();
    let mut produced = 0usize;
    let mut ok = 0usize;
    for r in &recs {
        if r.key.is_empty() {
            continue;
        }
        let has_stereo_or_iso = ["/b", "/t", "/m", "/s", "/i"]
            .iter()
            .any(|t| r.inchi.contains(t));
        if has_stereo_or_iso {
            continue;
        }
        let Ok(g) = build_molecule_graph(&r.smiles) else {
            continue;
        };
        if let Ok(got) = molrs::inchi::to_inchi_key(&g) {
            produced += 1;
            if got == r.key {
                ok += 1;
            }
        }
    }
    let acc = ok as f64 / produced.max(1) as f64;
    println!(
        "to_inchi_key: produced {produced}, {ok} exact ({:.2}%)",
        acc * 100.0
    );
    // フル InChI 文字列一致と同率のはず (キー機構は非立体で完全)
    assert!(acc >= 0.9741, "to_inchi_key accuracy {acc:.4} < 0.9741");
    // 注: full InChI 97.42% / key 97.44% (I19 §3.2 精緻化時点)
}
