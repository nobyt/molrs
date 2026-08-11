# molrs mutable-editing API — contract spec

**Status: implemented in `molrs/src/edit.rs`** (`editing` Cargo feature, off
by default). **Before wiring `crates/core::edit` against this, read that
file's module doc comment** — it lists five concrete deviations from this
document, found by tracing invariants `canon.rs`/`inchi/*.rs`/`stereo.rs`
actually rely on that this spec (written from `graph.rs`/`lib.rs` alone)
didn't have visibility into. The two signature-shape ones matter most for
`core::edit`'s call sites: `add_atom` returns `(usize, Vec<usize>)` (not bare
`usize`) and `add_bond` returns `Vec<usize>` (not `()`) — both because
`canon.rs` requires heavy atoms stay contiguous at `[0, n_heavy)`, which the
original signatures can't preserve when a heavy atom is added or bonded onto
an under-substituted atom. This is exactly the §7 reconciliation step below,
already done for the molrs side.

molrs today (verified against `molrs/molrs/src/graph.rs`, `lib.rs`) is parse-in /
derive-out only: `graph::build_molecule_graph(smiles: &str) -> Result<MoleculeGraph, ChemError>`
is the sole way to construct a `MoleculeGraph`, and every other module
(`canon`, `conformer`, `depict`, `inchi`, `substructure`) only reads `&MoleculeGraph`.
There is no way to add/remove an atom or bond, or change a bond order, without
round-tripping through SMILES text. This spec adds that.

## 1. Why an in-molrs API (not a crustalline-side graph layer)

Decided with the user: editing must go through molrs, not a parallel mutable
graph crustalline maintains itself. Rationale: molrs already owns every
invariant a mutation needs to preserve — implicit-H bookkeeping, valence
tables (`default_valences`, `aromatic_takes_double_bond` in `graph.rs`),
SSSR (`rings::symmetrized_sssr`), aromaticity/kekulization
(`aromaticity::{kekulize, perceive_aromaticity}`), CIP stereo
(`stereo::assign_stereochemistry`). Reimplementing any of that in crustalline
to keep a second mutable graph consistent would duplicate RDKit-parity-gated
logic molrs already validates against a large corpus. crustalline should stay
a thin consumer.

## 2. Scope of this API

A new module, suggested path `molrs::edit`, exposing atom/bond mutation
functions over `&mut MoleculeGraph`. Suggested Cargo feature gate: `editing`
(default off), so the existing "frozen, RDKit-parity-gated" core in
`RUST_PORT_HANDOFF.md` stays untouched by default and this is treated as an
additive, separately-validated surface. This is a suggestion for the molrs
maintainer to accept or reject, not a requirement crustalline can enforce.

## 3. Design decisions this contract makes

### 3.1 Mutate in place, not return-a-new-graph

```rust
pub fn add_bond(g: &mut MoleculeGraph, a: usize, b: usize, order: f64) -> Result<(), EditError>;
```

not

```rust
pub fn add_bond(g: &MoleculeGraph, a: usize, b: usize, order: f64) -> Result<MoleculeGraph, EditError>;
```

`MoleculeGraph`'s existing fields (`adjacency: Vec<Vec<usize>>`,
`bond_orders: HashMap<(usize,usize), f64>`, `kekule_bond_orders: Vec<f64>`)
are already indexed by `atom_idx` / `(i,j)` pairs computed once at parse time
— in-place mutation lets molrs update only the touched indices instead of
rebuilding the whole graph. **crustalline's `core::edit` is responsible for
cloning the graph onto its own undo stack before calling a mutator** — molrs
does not snapshot history itself. (`MoleculeGraph` already derives `Clone`,
confirmed in `graph.rs`, so this is cheap for crustalline to do.)

### 3.2 Index stability on removal

Removing `atoms[i]` forces a choice: shift every later index down, or leave a
tombstoned/sparse graph. This spec chooses **shifting, with a returned remap
table**:

```rust
pub fn remove_atom(g: &mut MoleculeGraph, atom_idx: usize) -> Result<Vec<Option<usize>>, EditError>;
```

The returned `Vec<Option<usize>>` has one entry per *pre-removal* atom index:
`Some(new_idx)` if it survived (shifted), `None` if it was removed (the target
atom itself, plus any implicit-H atoms cascade-removed — see 3.4). crustalline
uses this to translate any selection/hover/undo-stack state it holds by atom
index. This was chosen over molrs inventing stable opaque atom IDs, which
would be a bigger change rippling through every existing molrs consumer
(`smiles2iupac`, `canon`, `depict`, tests) that currently assumes atom index
== array position.

### 3.3 Derived-data invalidation — mutators keep the graph self-consistent

`MoleculeGraph` carries several fields computed once in `build_molecule_graph`:
`adjacency`, `bond_orders`, `ring_atom_sets` (SSSR), `kekule_bond_orders`,
plus per-atom `is_aromatic`/`in_ring` in `AtomInfo`. Two options were
considered for keeping these correct after a mutation:

- (a) every mutator updates all derived fields itself before returning, or
- (b) mutators only touch what's cheap/local (`adjacency`, `bond_orders`) and
  leave ring/aromaticity fields stale, requiring callers to explicitly call a
  recompute function before any ring-dependent operation (aromaticity,
  canonicalization, 2D depiction).

**This spec picks a hybrid, but end state is (a) — the graph must be
self-consistent immediately after any mutator returns.** Concretely: mutators
update `adjacency`/`bond_orders` directly (cheap), then internally re-run the
existing ring/aromaticity perception routines (`rings::symmetrized_sssr`,
`aromaticity::perceive_aromaticity`) before returning — the same cost as a
fresh parse, but without re-parsing SMILES text. Rationale: crustalline calls
`to_canonical_smiles`/`compute_coords_2d`/etc. immediately after edits and
must not have to remember a separate invalidation step. An explicit escape
hatch is still provided for completeness / testing:

```rust
pub fn recompute_derived(g: &mut MoleculeGraph);
```

crustalline is expected to never call this directly — it exists in case a
mutator's internal recompute has a gap, or for molrs's own test harness.

### 3.4 Implicit hydrogens stay explicit atom nodes

molrs's invariant (per `graph.rs` doc comment) is that **all hydrogens are
explicit atom nodes** post-`build_molecule_graph`, `num_hs` is always 0 on
`AtomInfo`. Mutators must preserve this:

- `add_atom` computes the new atom's implicit H count the same way
  `build_molecule_graph` does (via `default_valences`/`aromatic_takes_double_bond`)
  and appends that many explicit `H` atom nodes, bonded with single bonds —
  not a `num_hs` field bump.
- `remove_bond` / `set_bond_order` (order decrease) must recompute the
  affected atoms' implicit-H count and add H atom nodes if valence now allows
  more; must also handle the case where an added bond's atom already had
  "spare" implicit Hs that should be removed (e.g. adding a bond that fills
  what used to be satisfied by an implicit H — remove one H atom node,
  re-run 3.2's index-shift machinery since that's itself an atom removal).
- `remove_atom` on a heavy atom must also cascade-remove any H atom nodes that
  existed only to satisfy *that* atom's valence, and the returned remap
  (3.2) must reflect all of them being removed, not just the target atom.

This is the trickiest correctness surface in the whole contract — it's the
same "H bookkeeping" logic `build_molecule_graph` already implements once at
parse time (see `graph.rs` §2–5), so the recommended implementation strategy
for the molrs-side session is: **factor the existing implicit-H computation
in `build_molecule_graph` into a reusable per-atom function, and call it from
both the parser path and every mutator**, rather than writing new H-counting
logic from scratch.

### 3.5 Conformers are a crustalline-side concern, not molrs's

`Conformer { coords: Vec<Vec3> }` (from `embed_molecule`) is a separate struct,
never stored on `MoleculeGraph`. Mutators must not attempt to update or
invalidate any conformer — they don't have one to touch. **`crates/core::derive`
is always responsible for re-running `embed_molecule` after a structural
edit** (reusing the same seed for view stability across incremental edits,
per crustalline's `MoleculeState::embed_seed`). This also means a batch of
programmatic edits (e.g. building a molecule atom-by-atom from a script)
should only pay the UFF-refinement cost once at the end, by calling mutators
directly in a loop and deferring `embed_molecule` until the caller is done —
already possible with this API shape, no batching API needed on molrs's side.

### 3.6 Validation is atomic

Every mutator validates **before** committing any change, and returns `Err`
without partially applying the edit:

```rust
pub enum EditError {
    ValenceExceeded { atom_idx: usize, attempted: u32, max: u32 },
    AtomNotFound(usize),
    BondNotFound(usize, usize),
    BondAlreadyExists(usize, usize),
    InvalidBondOrder(f64),   // not in {1.0, 1.5, 2.0, 3.0}, or aromatic (1.5) proposed on a non-ring bond
}
```

Rationale: crustalline's UI needs to reject a single edit and show one
specific error, not partially mutate the graph and require the user to
manually undo. (`1.5`/aromatic bond orders are almost certainly wrong for a
user-initiated `add_bond`/`set_bond_order` call since aromaticity is normally
*derived*, not asserted — flagged here as an open question for the molrs
session: should `InvalidBondOrder` reject `1.5` outright for interactively-set
bonds, only allowing `{1.0, 2.0, 3.0}` plus letting the aromaticity perception
step in 3.3 promote rings to 1.5 on its own? This spec's recommendation is
yes, reject `1.5` from mutator callers.)

## 4. Proposed signatures (module `molrs::edit`)

```rust
use crate::graph::MoleculeGraph;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EditError {
    ValenceExceeded { atom_idx: usize, attempted: u32, max: u32 },
    AtomNotFound(usize),
    BondNotFound(usize, usize),
    BondAlreadyExists(usize, usize),
    InvalidBondOrder(f64),
}

/// Adds a new heavy atom, optionally bonded to an existing atom. Implicit H
/// count for the new atom (and any adjustment to `bonded_to`'s implicit H
/// count) is computed the same way `build_molecule_graph` computes it for
/// parsed atoms (see §3.4). Returns the new atom's index.
pub fn add_atom(
    g: &mut MoleculeGraph,
    symbol: &str,
    formal_charge: i8,
    bonded_to: Option<(usize, f64)>, // (existing atom idx, bond order)
) -> Result<usize, EditError>;

/// Removes an atom and all incident bonds, cascading to any implicit-H atoms
/// that existed only to satisfy this atom's valence (§3.4). Returns the
/// old-idx -> new-idx remap for every atom in the pre-removal graph (§3.2).
pub fn remove_atom(g: &mut MoleculeGraph, atom_idx: usize) -> Result<Vec<Option<usize>>, EditError>;

pub fn add_bond(g: &mut MoleculeGraph, a: usize, b: usize, order: f64) -> Result<(), EditError>;

/// Removes a bond. If this drops either endpoint's used valence, implicit H
/// count is recomputed and H atom nodes are added as needed (§3.4) — this
/// itself does not require index-shifting since it only adds atoms.
pub fn remove_bond(g: &mut MoleculeGraph, a: usize, b: usize) -> Result<(), EditError>;

/// Changing order may add or remove implicit-H atom nodes on both endpoints
/// (§3.4); removing H atoms shifts indices, so this can also return a remap.
/// (Open question for the molrs session: unify this with remove_bond's
/// return type, i.e. should both return `Vec<Option<usize>>` for consistency
/// even when no removal occurred, or only when one did? Recommendation:
/// always return it — `Ok(identity_remap)` when nothing shifted — so
/// crustalline has one code path instead of branching on whether a remap
/// happened.)
pub fn set_bond_order(g: &mut MoleculeGraph, a: usize, b: usize, order: f64) -> Result<Vec<Option<usize>>, EditError>;

pub fn set_formal_charge(g: &mut MoleculeGraph, atom_idx: usize, charge: i8) -> Result<(), EditError>;

/// Escape hatch: forces re-derivation of ring_atom_sets / is_aromatic /
/// kekule_bond_orders without a structural change. Expected unused by
/// crustalline (§3.3) — documented for completeness / molrs's own testing.
pub fn recompute_derived(g: &mut MoleculeGraph);
```

## 5. Interaction with existing molrs functions

No signature changes needed anywhere else. `to_canonical_smiles(g)`,
`to_canonical_smiles_with_order(g)`, `conformer::embed_molecule(g, params)`,
`conformer::molblock::to_mol_block(g, conf, title)`,
`depict::compute_coords_2d(g, params)`, `depict::to_svg(g, coords, style)`,
`substructure::substruct_matches(g, ...)` all continue to take `&MoleculeGraph`
unchanged — per §3.3, a post-mutation graph is a valid input to all of them
with no new parameters, exactly as if it had been produced fresh by
`build_molecule_graph`.

## 6. How crustalline will consume this

`crates/core::edit::apply_edit` (crustalline repo) is the sole call site.
Each `EditCommand` variant (`AddAtom`, `RemoveAtom`, `AddBond`, `RemoveBond`,
`SetBondOrder`, `SetFormalCharge`) maps to exactly one call into this module.
Before calling a mutator, `core::edit` clones the current `MoleculeGraph` onto
`MoleculeState`'s undo stack (§3.1). After a successful mutation,
`core::derive` re-runs `embed_molecule` and `to_mol_block` (§3.5) and the
result is broadcast to the frontend. On `Err(EditError)`, no state changes and
the error is surfaced to the UI as-is (mapped to a `crustalline_ipc_types::IpcError`).

## 7. Reconciliation step when this lands

crustalline pins an exact molrs git rev/path; this spec is not binding on the
molrs implementation. Per the crustalline build-out plan, milestone M5's first
step is: diff the actual `molrs::edit` signatures against this document and
reconcile any differences in `crates/core::edit` before wiring Tauri commands
to it — do not assume the contract landed unchanged.
