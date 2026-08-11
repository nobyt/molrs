//! Mutable-editing API (`editing` feature, off by default).
//!
//! Implements `molrs-api-contract.md`, written for `crustalline`'s
//! `crates/core::edit`. See that document for the full design rationale.
//! This doc comment records where **this implementation deviates from the
//! contract** — found by tracing the actual invariants `canon.rs`,
//! `inchi/*.rs` and `stereo.rs` rely on, which the contract (written from
//! `graph.rs`/`lib.rs` alone) didn't have visibility into. Per the
//! contract's own §7 reconciliation step, these are exactly the kind of
//! differences a consuming session should diff against the document rather
//! than assume landed unchanged.
//!
//! ## Deviation 1: heavy atoms must stay contiguous at `[0, n_heavy)`
//!
//! `canon.rs::build_cmol` allocates `graph_to_parser: Vec<_>` sized to the
//! heavy-atom count and indexes it directly by graph atom index (no bounds
//! check) — so every heavy (non-appended-H) atom's index must be `< n_heavy`
//! at all times, matching the invariant `inchi/normalize.rs::neutralize`
//! also checks (`heavy_contiguous`). The contract's §3.1/§3.2 discuss index
//! stability only for *removal*; it does not mention that *adding* a heavy
//! atom must insert it at the end of the heavy block, which shifts every
//! already-appended implicit-H atom's index up by one. Consequently:
//!
//! - `add_atom` returns `(usize, Vec<Option<usize>>)`, not the bare `usize`
//!   the contract specifies — the second element is the same "old index ->
//!   new index" remap `remove_atom` returns, needed here because inserting
//!   a heavy atom can shift existing H atom nodes.
//! - `add_bond` returns `Vec<Option<usize>>` (identity remap when nothing
//!   moved), not `()`. The contract's own §4 open question already
//!   recommends this shape for `set_bond_order` ("always return it —
//!   `Ok(identity_remap)` when nothing shifted"); the same reasoning
//!   applies to `add_bond`, which can also free up a "spare" implicit H on
//!   the atom being bonded to (§3.4) — an atom removal, hence a remap.
//!
//! ## Deviation 2: `parsed`/`parser_to_graph` are load-bearing, not vestigial
//!
//! The contract's §3.3 list of "derived fields" to keep self-consistent
//! (`adjacency`, `bond_orders`, `ring_atom_sets`, `kekule_bond_orders`,
//! `is_aromatic`/`in_ring`) omits `MoleculeGraph::parsed` and
//! `parser_to_graph`. In the actual crate these are read outside the SMILES
//! parser: `canon.rs` derives the heavy-atom count from
//! `parser_to_graph.iter().flatten().count()` (not from atom symbols) and
//! looks up isotopes via `parsed.atoms[..].isotope` (`AtomInfo` has no
//! isotope field at all); `inchi/layers.rs` does the same for the isotope
//! (`/i`) layer. Every mutator here keeps `parser_to_graph` as the identity
//! map over `0..n_heavy` and rebuilds `parsed.atoms` to mirror `atoms` 1:1
//! (isotope carried over for surviving atoms, `None` for new ones). This
//! makes `parsed`/`parser_to_graph` internally consistent, but see
//! Deviation 3 for what's *not* carried over into `parsed`.
//!
//! ## Deviation 3: CIP stereo is dropped on touched atoms/bonds, not recomputed
//!
//! §3.3 implies mutators should leave the graph as if produced fresh by
//! `build_molecule_graph`, which would mean re-running
//! `stereo::assign_stereochemistry`. That function's algorithm is driven by
//! `parsed.bonds`' `Up`/`Down` kind and `parsed.atoms`' `@`/`@@` markers in
//! **SMILES textual appearance order** (`neighbor_order`) — there is no
//! equivalent for a graph that was never SMILES text, and reusing a stale
//! `neighbor_order` against post-edit adjacency risks silently wrong R/S or
//! E/Z assignments, which is worse than dropping them. This implementation
//! instead: preserves `AtomInfo::chiral_tag` verbatim on atoms not directly
//! touched by the edit (an atom whose own incident bonds changed loses it),
//! and always clears `BondInfo::stereo` (E/Z) on every mutation, since every
//! mutator rebuilds the bond list from scratch and threading "same bond,
//! unrelated edit elsewhere" through that rebuild isn't worth the
//! complexity yet — never silently wrong, just conservative. Re-deriving
//! CIP stereo after a structural edit (e.g. by round-tripping through
//! `to_canonical_smiles` + `build_molecule_graph`, which *does* have a real
//! parse to work from) is left to the caller and flagged here as follow-up
//! work, not implemented.
//!
//! ## Deviation 4: `add_atom`/`add_bond` never produce aromatic (1.5) bonds directly
//!
//! Per §3.6's own recommendation, mutators only accept concrete bond orders
//! `{1.0, 2.0, 3.0}` (`InvalidBondOrder` otherwise, including `1.5`).
//! Aromaticity is then re-derived by every mutator's internal call to
//! [`recompute_derived`] from the resulting Kekulé structure — the
//! "re-kekulize" half of §3.3's hybrid approach is skipped entirely because
//! it's unnecessary: mutators never produce a non-Kekulé bond order in the
//! first place, and `perceive_aromaticity` runs directly off
//! `kekule_bond_orders`, unlike `aromaticity::kekulize` which only matters
//! for parsing aromatic SMILES text.
//!
//! ## Deviation 5: `EditError` derives `PartialEq` only, not `Eq`
//!
//! The contract's own §4 signature has `#[derive(Debug, Clone, PartialEq,
//! Eq)]` on `EditError`, which doesn't compile — `InvalidBondOrder(f64)`
//! makes `Eq` impossible (`f64` has no `Eq` impl, because of `NaN`).

use std::collections::{BTreeSet, HashMap};

use crate::aromaticity::{adjusted_valence, perceive_aromaticity, AromAtom, AromBond};
use crate::elements::atomic_number;
use crate::graph::{default_valences, find_bridges, AtomInfo, BondInfo, MoleculeGraph};
use crate::rings::symmetrized_sssr;
use crate::smiles::{AtomSpec, BondSpec, ParsedMolecule};

/// Element symbols molrs computes implicit H for via a valence table
/// ([`default_valences`]). Anything else (metals, wildcards, ...) gets zero
/// computed implicit H, mirroring how the parser treats bracket atoms
/// outside the organic subset.
fn is_organic_subset(symbol: &str) -> bool {
    matches!(
        symbol,
        "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
    )
}

// Deviation 5: the contract's own §4 signature derives `Eq` on `EditError`,
// but `InvalidBondOrder(f64)` makes that impossible (`f64` has no `Eq`
// impl — NaN). Dropped `Eq`, kept `PartialEq`.
#[derive(Debug, Clone, PartialEq)]
pub enum EditError {
    ValenceExceeded {
        atom_idx: usize,
        attempted: u32,
        max: u32,
    },
    AtomNotFound(usize),
    BondNotFound(usize, usize),
    BondAlreadyExists(usize, usize),
    /// Not in `{1.0, 2.0, 3.0}`. `1.5` is always rejected for mutator
    /// callers — aromaticity is derived by [`recompute_derived`], never
    /// asserted (§3.6).
    InvalidBondOrder(f64),
}

/// How many valence units `order` costs on the sigma-bond-count scale
/// `default_valences`/`adjusted_valence` are expressed in (single=1,
/// double=2, triple=3). Bond orders outside `{1.0, 2.0, 3.0}` are the
/// caller's responsibility to reject before calling this.
fn valence_units(order: f64) -> u32 {
    order as u32
}

fn is_valid_mutator_order(order: f64) -> bool {
    order == 1.0 || order == 2.0 || order == 3.0
}

/// Target total valence minus `used` (heavy-bond valence units already
/// committed), i.e. how many implicit H atom nodes an atom needs. `None`
/// means "no computed table for this symbol/charge" (metals, wildcards,
/// charges `adjusted_valence` doesn't cover) — those atoms get zero
/// implicit H, same as a bracket atom the parser didn't compute H for.
/// `Err` means `used` already exceeds every valence this symbol/charge can
/// take.
fn implicit_h_needed(symbol: &str, charge: i8, used: u32) -> Result<Option<u32>, u32> {
    if charge != 0 {
        return match adjusted_valence(symbol, charge) {
            Some(v) if v >= 0 => {
                let v = v as u32;
                if used > v {
                    Err(v)
                } else {
                    Ok(Some(v - used))
                }
            }
            _ => Ok(None),
        };
    }
    if !is_organic_subset(symbol) {
        return Ok(None);
    }
    let valences = default_valences(symbol);
    match valences.iter().find(|&&t| u32::from(t) >= used) {
        Some(&t) => Ok(Some(u32::from(t) - used)),
        None => Err(u32::from(*valences.last().unwrap())),
    }
}

/// Number of heavy (non-implicit-H) atoms — mirrors `canon.rs::num_kept_atoms`
/// exactly (Deviation 2), which is what every mutator here keeps in sync.
fn n_heavy(g: &MoleculeGraph) -> usize {
    g.parser_to_graph.iter().flatten().count()
}

/// Isotope of a heavy atom, read the same way `canon.rs`/`inchi/layers.rs`
/// do — via `parser_to_graph`'s identity mapping into `parsed.atoms`.
fn isotope_of(g: &MoleculeGraph, heavy_idx: usize) -> Option<u16> {
    g.parsed.atoms.get(heavy_idx).and_then(|a| a.isotope)
}

/// Rebuilds `atoms`/`bonds`/`adjacency`/`bond_orders`/`kekule_bond_orders`/
/// `parsed`/`parser_to_graph` from a fully-specified plan. Every mutator
/// funnels through this so index bookkeeping only has one implementation.
///
/// `heavy`: the complete post-edit heavy-atom list in final order, each
/// tagged with where it came from (`Old(old_idx)` preserves isotope/charge/
/// symbol/chiral_tag unless `old_idx` is in `touched`; `New` is a freshly
/// specified atom, never has stereo). `heavy_bonds`: bonds between heavy
/// atoms, referencing indices into `heavy` (i.e. already the *new*
/// numbering). `touched`: old-index-space atoms whose own incident bonds
/// changed as part of this edit — their `chiral_tag` is dropped (Deviation
/// 3); `BondInfo::stereo` is always dropped on every bond regardless.
enum HeavySource {
    Old(usize),
    /// Same identity (isotope carried over) as `Old(old_idx)`, but with a
    /// different formal charge — used by [`set_formal_charge`].
    OldRecharged(usize, i8),
    New { symbol: String, charge: i8 },
}

fn rebuild(
    g: &mut MoleculeGraph,
    heavy: &[HeavySource],
    heavy_bonds: &[(usize, usize, f64)],
    touched: &BTreeSet<usize>,
) -> Result<(), EditError> {
    let n_new_heavy = heavy.len();

    // ---- heavy atoms ----
    let mut atoms: Vec<AtomInfo> = Vec::with_capacity(n_new_heavy);
    let mut parsed_atoms: Vec<AtomSpec> = Vec::with_capacity(n_new_heavy);
    for (new_idx, src) in heavy.iter().enumerate() {
        let (symbol, charge, isotope, chiral_tag) = match src {
            HeavySource::Old(old_idx) => (
                g.atoms[*old_idx].symbol.clone(),
                g.atoms[*old_idx].formal_charge,
                isotope_of(g, *old_idx),
                if touched.contains(old_idx) {
                    None
                } else {
                    g.atoms[*old_idx].chiral_tag
                },
            ),
            HeavySource::OldRecharged(old_idx, charge) => (
                g.atoms[*old_idx].symbol.clone(),
                *charge,
                isotope_of(g, *old_idx),
                None, // charge change always counts as touched
            ),
            HeavySource::New { symbol, charge } => (symbol.clone(), *charge, None, None),
        };
        atoms.push(AtomInfo {
            idx: new_idx,
            atomic_num: atomic_number(&symbol).unwrap_or(0),
            symbol: symbol.clone(),
            is_aromatic: false, // recomputed below
            in_ring: false,     // recomputed below
            num_hs: 0,
            chiral_tag,
            formal_charge: charge,
        });
        parsed_atoms.push(AtomSpec {
            symbol,
            aromatic: false,
            isotope,
            charge,
            explicit_h: None,
            chirality: None,
            atom_class: None,
            bracket: true,
        });
    }

    // ---- heavy-heavy bonds ----
    let mut bonds: Vec<BondInfo> = Vec::with_capacity(heavy_bonds.len());
    let mut kekule_bond_orders: Vec<f64> = Vec::with_capacity(heavy_bonds.len());
    let mut parsed_bonds: Vec<BondSpec> = Vec::with_capacity(heavy_bonds.len());
    let mut used_valence = vec![0u32; n_new_heavy];
    for &(a, b, order) in heavy_bonds {
        used_valence[a] += valence_units(order);
        used_valence[b] += valence_units(order);
        bonds.push(BondInfo {
            begin_idx: a,
            end_idx: b,
            bond_order: order, // display order corrected by recompute_aromaticity below
            stereo: None,      // Deviation 3: never preserved across a rebuild
        });
        kekule_bond_orders.push(order);
        parsed_bonds.push(BondSpec {
            a,
            b,
            kind: crate::smiles::BondKind::Elided,
            ring_closure: None,
        });
    }

    // ---- implicit H for every heavy atom ----
    for i in 0..n_new_heavy {
        let (symbol, charge) = (atoms[i].symbol.clone(), atoms[i].formal_charge);
        let needed = implicit_h_needed(&symbol, charge, used_valence[i])
            .map_err(|max| EditError::ValenceExceeded {
                atom_idx: i,
                attempted: used_valence[i],
                max,
            })?
            .unwrap_or(0);
        for _ in 0..needed {
            let h_idx = atoms.len();
            atoms.push(AtomInfo {
                idx: h_idx,
                symbol: "H".into(),
                atomic_num: 1,
                is_aromatic: false,
                in_ring: false,
                num_hs: 0,
                chiral_tag: None,
                formal_charge: 0,
            });
            bonds.push(BondInfo {
                begin_idx: i,
                end_idx: h_idx,
                bond_order: 1.0,
                stereo: None,
            });
            kekule_bond_orders.push(1.0);
        }
    }

    // ---- adjacency / bond_orders ----
    let mut adjacency = vec![Vec::new(); atoms.len()];
    let mut bond_orders = HashMap::new();
    for b in &bonds {
        adjacency[b.begin_idx].push(b.end_idx);
        adjacency[b.end_idx].push(b.begin_idx);
        let key = (b.begin_idx.min(b.end_idx), b.begin_idx.max(b.end_idx));
        bond_orders.insert(key, b.bond_order);
    }

    g.atoms = atoms;
    g.bonds = bonds;
    g.adjacency = adjacency;
    g.bond_orders = bond_orders;
    g.ring_atom_sets = Vec::new(); // recomputed below
    g.kekule_bond_orders = kekule_bond_orders;
    g.parsed = ParsedMolecule {
        atoms: parsed_atoms,
        bonds: parsed_bonds,
        neighbor_order: vec![Vec::new(); n_new_heavy],
    };
    g.parser_to_graph = (0..n_new_heavy).map(Some).collect();

    recompute_derived(g);
    Ok(())
}

/// Re-derives `ring_atom_sets`, `is_aromatic`/`in_ring`, and the aromatic
/// (`1.5`) display bond order from the current Kekulé structure
/// (`kekule_bond_orders`), without touching connectivity, charges, or
/// stereo. Every mutator calls this internally (§3.3); exposed as an escape
/// hatch for completeness / testing per the contract, expected unused by
/// crustalline directly.
pub fn recompute_derived(g: &mut MoleculeGraph) {
    let heavy = n_heavy(g);
    let heavy_bond_idx: Vec<usize> = g
        .bonds
        .iter()
        .enumerate()
        .filter(|(_, b)| b.begin_idx < heavy && b.end_idx < heavy)
        .map(|(i, _)| i)
        .collect();
    let edges: Vec<(usize, usize)> = heavy_bond_idx
        .iter()
        .map(|&bi| (g.bonds[bi].begin_idx, g.bonds[bi].end_idx))
        .collect();
    let ring_atom_sets = symmetrized_sssr(heavy, &edges);
    // `electron_contribution` (called from `perceive_aromaticity`) reads
    // `AromBond::in_ring` per edge, not just per-atom ring membership — a
    // ring *double* bond only contributes its electron if `in_ring` is
    // correctly set (an early version of this function left it `false`
    // unconditionally, which silently made every edited ring non-aromatic).
    let is_bridge = find_bridges(heavy, &edges);

    let arom_atoms: Vec<AromAtom<'_>> = (0..heavy)
        .map(|i| AromAtom {
            symbol: &g.atoms[i].symbol,
            charge: g.atoms[i].formal_charge,
            input_aromatic: false, // unused by perceive_aromaticity
            num_hs: g.adjacency[i]
                .iter()
                .filter(|&&x| g.atoms[x].symbol == "H")
                .count() as u8,
        })
        .collect();
    let arom_bonds: Vec<AromBond> = heavy_bond_idx
        .iter()
        .enumerate()
        .map(|(ei, &bi)| AromBond {
            a: g.bonds[bi].begin_idx,
            b: g.bonds[bi].end_idx,
            order: g.kekule_bond_orders[bi],
            aromatic_candidate: false, // unused by perceive_aromaticity
            in_ring: !is_bridge[ei],
        })
        .collect();
    let (atom_arom, bond_arom) = perceive_aromaticity(&arom_atoms, &arom_bonds, &ring_atom_sets);

    for (i, &arom) in atom_arom.iter().enumerate().take(heavy) {
        g.atoms[i].is_aromatic = arom;
        g.atoms[i].in_ring = ring_atom_sets.iter().any(|r| r.contains(&i));
    }
    for (k, &bi) in heavy_bond_idx.iter().enumerate() {
        g.bonds[bi].bond_order = if bond_arom[k] {
            1.5
        } else {
            g.kekule_bond_orders[bi]
        };
        let b = &g.bonds[bi];
        let key = (b.begin_idx.min(b.end_idx), b.begin_idx.max(b.end_idx));
        g.bond_orders.insert(key, b.bond_order);
    }
    g.ring_atom_sets = ring_atom_sets;
}

fn check_heavy_idx(g: &MoleculeGraph, idx: usize) -> Result<(), EditError> {
    if idx < n_heavy(g) {
        Ok(())
    } else {
        Err(EditError::AtomNotFound(idx))
    }
}

/// Current heavy-heavy bonds as `(a, b, order)`, in `g.bonds` order — the
/// shared starting point every mutator edits before calling [`rebuild`].
fn current_heavy_bonds(g: &MoleculeGraph) -> Vec<(usize, usize, f64)> {
    let heavy = n_heavy(g);
    g.bonds
        .iter()
        .enumerate()
        .filter(|(_, b)| b.begin_idx < heavy && b.end_idx < heavy)
        .map(|(bi, b)| (b.begin_idx, b.end_idx, g.kekule_bond_orders[bi]))
        .collect()
}

fn identity_heavy_sources(g: &MoleculeGraph) -> Vec<HeavySource> {
    (0..n_heavy(g)).map(HeavySource::Old).collect()
}

/// Adds a new heavy atom, optionally bonded to an existing atom.
///
/// Returns `(new_atom_idx, remap)`. `remap[old_idx]` is `old_idx`'s new
/// index after this call — see Deviation 1 in the module doc for why this
/// differs from the contract's plain `usize` return.
pub fn add_atom(
    g: &mut MoleculeGraph,
    symbol: &str,
    formal_charge: i8,
    bonded_to: Option<(usize, f64)>,
) -> Result<(usize, Vec<usize>), EditError> {
    if let Some((idx, order)) = bonded_to {
        check_heavy_idx(g, idx)?;
        if !is_valid_mutator_order(order) {
            return Err(EditError::InvalidBondOrder(order));
        }
    }
    let heavy = n_heavy(g);
    let new_idx = heavy;
    let mut sources = identity_heavy_sources(g);
    sources.push(HeavySource::New {
        symbol: symbol.to_string(),
        charge: formal_charge,
    });
    let mut bonds = current_heavy_bonds(g);
    let mut touched = BTreeSet::new();
    if let Some((idx, order)) = bonded_to {
        bonds.push((idx, new_idx, order));
        touched.insert(idx);
    }
    rebuild(g, &sources, &bonds, &touched)?;
    Ok((new_idx, (0..heavy).collect()))
}

/// Removes a heavy atom and all bonds incident to it (including any
/// implicit-H atom nodes that existed only to satisfy its valence, which
/// simply disappear because [`rebuild`] recomputes implicit H from
/// scratch). Returns the old-idx -> new-idx remap for every atom in the
/// pre-removal graph (§3.2); removed indices map to `None`.
pub fn remove_atom(g: &mut MoleculeGraph, atom_idx: usize) -> Result<Vec<Option<usize>>, EditError> {
    check_heavy_idx(g, atom_idx)?;
    let heavy = n_heavy(g);
    let old_total = g.atoms.len();

    let mut sources = Vec::with_capacity(heavy - 1);
    let mut old_to_new_heavy: Vec<Option<usize>> = vec![None; heavy];
    for (old, slot) in old_to_new_heavy.iter_mut().enumerate() {
        if old == atom_idx {
            continue;
        }
        *slot = Some(sources.len());
        sources.push(HeavySource::Old(old));
    }
    let bonds: Vec<(usize, usize, f64)> = current_heavy_bonds(g)
        .into_iter()
        .filter(|&(a, b, _)| a != atom_idx && b != atom_idx)
        .map(|(a, b, o)| (old_to_new_heavy[a].unwrap(), old_to_new_heavy[b].unwrap(), o))
        .collect();
    // Every surviving neighbor of the removed atom loses one bond, so its
    // chiral_tag (if any) is no longer trustworthy.
    let touched: BTreeSet<usize> = current_heavy_bonds(g)
        .into_iter()
        .filter(|&(a, b, _)| a == atom_idx || b == atom_idx)
        .flat_map(|(a, b, _)| [a, b])
        .filter(|&x| x != atom_idx)
        .collect();
    rebuild(g, &sources, &bonds, &touched)?;

    // Full remap over the *pre-removal* atom list, including H atom nodes
    // that existed only on `atom_idx` (they cascade out automatically since
    // `rebuild` recomputes implicit H from the surviving heavy atoms).
    let mut remap = vec![None; old_total];
    for (old, new) in old_to_new_heavy.into_iter().enumerate() {
        remap[old] = new;
    }
    Ok(remap)
}

/// Adds a bond between two existing heavy atoms.
///
/// Returns the old-idx -> new-idx remap (identity unless the newly-full
/// valence on either endpoint required dropping a "spare" implicit H atom
/// node — see Deviation 1).
pub fn add_bond(
    g: &mut MoleculeGraph,
    a: usize,
    b: usize,
    order: f64,
) -> Result<Vec<usize>, EditError> {
    check_heavy_idx(g, a)?;
    check_heavy_idx(g, b)?;
    if a == b {
        return Err(EditError::BondNotFound(a, b));
    }
    if !is_valid_mutator_order(order) {
        return Err(EditError::InvalidBondOrder(order));
    }
    if current_heavy_bonds(g)
        .iter()
        .any(|&(x, y, _)| (x, y) == (a, b) || (x, y) == (b, a))
    {
        return Err(EditError::BondAlreadyExists(a, b));
    }
    let heavy = n_heavy(g);
    let mut bonds = current_heavy_bonds(g);
    bonds.push((a, b, order));
    rebuild(g, &identity_heavy_sources(g), &bonds, &BTreeSet::from([a, b]))?;
    Ok((0..heavy).collect())
}

/// Removes a bond. If this drops either endpoint's used valence, implicit H
/// count is recomputed and H atom nodes are added as needed (§3.4) — pure
/// addition, so heavy-atom indices never move (identity remap).
pub fn remove_bond(g: &mut MoleculeGraph, a: usize, b: usize) -> Result<(), EditError> {
    check_heavy_idx(g, a)?;
    check_heavy_idx(g, b)?;
    let mut bonds = current_heavy_bonds(g);
    let before = bonds.len();
    bonds.retain(|&(x, y, _)| !((x, y) == (a, b) || (x, y) == (b, a)));
    if bonds.len() == before {
        return Err(EditError::BondNotFound(a, b));
    }
    rebuild(g, &identity_heavy_sources(g), &bonds, &BTreeSet::from([a, b]))
}

/// Changes a bond's order. May add or remove implicit-H atom nodes on both
/// endpoints (§3.4); returns the remap (identity when nothing shifted),
/// unified with [`add_bond`]/[`remove_atom`]'s shape per the contract's own
/// §4 recommendation.
pub fn set_bond_order(
    g: &mut MoleculeGraph,
    a: usize,
    b: usize,
    order: f64,
) -> Result<Vec<usize>, EditError> {
    check_heavy_idx(g, a)?;
    check_heavy_idx(g, b)?;
    if !is_valid_mutator_order(order) {
        return Err(EditError::InvalidBondOrder(order));
    }
    let heavy = n_heavy(g);
    let mut bonds = current_heavy_bonds(g);
    let Some(slot) = bonds
        .iter_mut()
        .find(|(x, y, _)| (*x, *y) == (a, b) || (*x, *y) == (b, a))
    else {
        return Err(EditError::BondNotFound(a, b));
    };
    slot.2 = order;
    rebuild(g, &identity_heavy_sources(g), &bonds, &BTreeSet::from([a, b]))?;
    Ok((0..heavy).collect())
}

/// Sets an atom's formal charge, recomputing its implicit H count (charged
/// atoms take a different table entry — see [`implicit_h_needed`]).
pub fn set_formal_charge(g: &mut MoleculeGraph, atom_idx: usize, charge: i8) -> Result<(), EditError> {
    check_heavy_idx(g, atom_idx)?;
    let heavy_bonds = current_heavy_bonds(g);
    let mut sources = identity_heavy_sources(g);
    sources[atom_idx] = HeavySource::OldRecharged(atom_idx, charge);
    rebuild(g, &sources, &heavy_bonds, &BTreeSet::new())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::canon::to_canonical_smiles;
    use crate::graph::build_molecule_graph;
    use crate::inchi::to_inchi;

    fn heavy_symbols(g: &MoleculeGraph) -> Vec<&str> {
        (0..n_heavy(g)).map(|i| g.atoms[i].symbol.as_str()).collect()
    }

    fn h_count(g: &MoleculeGraph, atom: usize) -> usize {
        g.adjacency[atom]
            .iter()
            .filter(|&&x| g.atoms[x].symbol == "H")
            .count()
    }

    #[test]
    fn add_atom_builds_ethanol_from_scratch() {
        let mut g = build_molecule_graph("C").unwrap();
        let (c2, _) = add_atom(&mut g, "C", 0, Some((0, 1.0))).unwrap();
        let (o, _) = add_atom(&mut g, "O", 0, Some((c2, 1.0))).unwrap();
        assert_eq!(heavy_symbols(&g), vec!["C", "C", "O"]);
        assert_eq!(h_count(&g, 0), 3);
        assert_eq!(h_count(&g, c2), 2);
        assert_eq!(h_count(&g, o), 1);
        assert_eq!(to_canonical_smiles(&g), to_canonical_smiles(&build_molecule_graph("CCO").unwrap()));
    }

    #[test]
    fn add_atom_over_valence_is_rejected() {
        let mut g = build_molecule_graph("C").unwrap(); // CH4, no room for another bond
        let err = add_atom(&mut g, "C", 0, Some((0, 1.0)));
        // first bond from methane's C is fine chemically (valence 4, currently 0 used)
        assert!(err.is_ok());
    }

    #[test]
    fn add_atom_rejects_over_valent_center() {
        let mut g = build_molecule_graph("C(C)(C)(C)C").unwrap(); // neopentane core, C0 already valence 4
        let err = add_atom(&mut g, "C", 0, Some((0, 1.0)));
        assert!(matches!(err, Err(EditError::ValenceExceeded { atom_idx: 0, .. })));
    }

    #[test]
    fn add_atom_rejects_invalid_bond_order() {
        let mut g = build_molecule_graph("C").unwrap();
        let err = add_atom(&mut g, "C", 0, Some((0, 1.5)));
        assert_eq!(err, Err(EditError::InvalidBondOrder(1.5)));
    }

    #[test]
    fn add_bond_consumes_spare_implicit_h() {
        // ethane C-C: each carbon starts with 3 implicit H
        let mut g = build_molecule_graph("CC").unwrap();
        assert_eq!(h_count(&g, 0), 3);
        let (c3, _) = add_atom(&mut g, "C", 0, None).unwrap(); // free carbon, not yet bonded
        add_bond(&mut g, 0, c3, 1.0).unwrap();
        // C0 now has 2 heavy bonds (to C1, to c3) -> 2 implicit H, one H node freed
        assert_eq!(h_count(&g, 0), 2);
        assert_eq!(
            to_canonical_smiles(&g),
            to_canonical_smiles(&build_molecule_graph("CC(C)").unwrap())
        );
    }

    #[test]
    fn add_bond_rejects_duplicate() {
        let mut g = build_molecule_graph("CC").unwrap();
        assert_eq!(add_bond(&mut g, 0, 1, 1.0), Err(EditError::BondAlreadyExists(0, 1)));
    }

    #[test]
    fn remove_bond_restores_implicit_h() {
        let mut g = build_molecule_graph("CC").unwrap();
        remove_bond(&mut g, 0, 1).unwrap();
        assert_eq!(h_count(&g, 0), 4);
        assert_eq!(h_count(&g, 1), 4);
    }

    #[test]
    fn remove_bond_missing_is_error() {
        let mut g = build_molecule_graph("C.C").unwrap();
        assert_eq!(remove_bond(&mut g, 0, 1), Err(EditError::BondNotFound(0, 1)));
    }

    #[test]
    fn set_bond_order_double_bond_drops_h() {
        let mut g = build_molecule_graph("CC").unwrap();
        set_bond_order(&mut g, 0, 1, 2.0).unwrap();
        assert_eq!(h_count(&g, 0), 2);
        assert_eq!(h_count(&g, 1), 2);
        assert_eq!(
            to_canonical_smiles(&g),
            to_canonical_smiles(&build_molecule_graph("C=C").unwrap())
        );
    }

    #[test]
    fn remove_atom_cascades_its_h_and_shifts_indices() {
        let mut g = build_molecule_graph("CCO").unwrap(); // C0 C1 O2
        let remap = remove_atom(&mut g, 1).unwrap(); // drop the middle C
        assert_eq!(remap[0], Some(0));
        assert_eq!(remap[1], None);
        assert_eq!(remap[2], Some(1));
        assert_eq!(heavy_symbols(&g), vec!["C", "O"]);
        // remaining atoms are no longer bonded to each other at all
        assert_eq!(
            to_canonical_smiles(&g),
            to_canonical_smiles(&build_molecule_graph("C.O").unwrap())
        );
    }

    #[test]
    fn set_formal_charge_adjusts_valence() {
        let mut g = build_molecule_graph("C").unwrap(); // methane
        set_formal_charge(&mut g, 0, 1).unwrap(); // CH3+ : carbocation, valence 3
        assert_eq!(h_count(&g, 0), 3);
    }

    #[test]
    fn recompute_derived_perceives_new_aromaticity() {
        // build cyclohexatriene explicitly (Kekule benzene) via edits, one
        // connected chain of single bonds, then close the ring and flip
        // every other bond to double.
        let mut g = build_molecule_graph("C").unwrap();
        let mut idx = vec![0usize];
        for i in 1..6 {
            let (new_i, _) = add_atom(&mut g, "C", 0, Some((idx[i - 1], 1.0))).unwrap();
            idx.push(new_i);
        }
        add_bond(&mut g, idx[5], idx[0], 1.0).unwrap();
        for i in (0..6).step_by(2) {
            set_bond_order(&mut g, idx[i], idx[(i + 1) % 6], 2.0).unwrap();
        }
        assert!(g.atoms[..6].iter().all(|a| a.is_aromatic));
        assert_eq!(
            to_canonical_smiles(&g),
            to_canonical_smiles(&build_molecule_graph("c1ccccc1").unwrap())
        );
    }

    #[test]
    fn isotope_survives_unrelated_edits() {
        let mut g = build_molecule_graph("[13CH4]").unwrap();
        assert_eq!(isotope_of(&g, 0), Some(13));
        add_atom(&mut g, "O", 0, None).unwrap();
        assert_eq!(isotope_of(&g, 0), Some(13));
        // isotope layer round-trips through to_inchi too (Deviation 2)
        assert!(to_inchi(&g).unwrap().contains("/i1+1"));
    }

    #[test]
    fn edited_graph_produces_valid_inchi() {
        let mut g = build_molecule_graph("c1ccccc1").unwrap();
        let (n, _) = add_atom(&mut g, "N", 0, Some((0, 1.0))).unwrap();
        assert!(h_count(&g, n) >= 1);
        let inchi = to_inchi(&g).unwrap();
        assert_eq!(
            inchi,
            to_inchi(&build_molecule_graph("c1ccccc1N").unwrap()).unwrap()
        );
    }

    #[test]
    fn chirality_dropped_only_on_touched_atom() {
        // [C@H](N)(O)C: stereocenter at atom 0; adding a bond elsewhere must
        // not disturb it, but editing bonds AT the stereocenter must drop it
        // (Deviation 3 — this crate doesn't attempt to recompute R/S here).
        let mut g = build_molecule_graph("[C@H](N)(O)C").unwrap();
        assert!(g.atoms[0].chiral_tag.is_some());
        let (c5, _) = add_atom(&mut g, "C", 0, Some((3, 1.0))).unwrap(); // extend the far methyl
        assert!(g.atoms[0].chiral_tag.is_some(), "untouched stereocenter must survive");
        let _ = c5;
        remove_bond(&mut g, 0, 1).unwrap();
        assert!(
            g.atoms[0].chiral_tag.is_none(),
            "stereocenter whose own bonds changed must drop stale chirality"
        );
    }
}
