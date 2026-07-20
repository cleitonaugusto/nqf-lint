# nqf-lint

**A pre-flight linter for quantum-chemistry cluster setups.**

A quantum-chemistry calculation can run for hours and return a number that *looks*
fine while the input was silently broken — an oxygen entered with no hydrogens, an
electron count made impossible by a forgotten pseudopotential core, a metal floating
with no coordination. The SCF either fails cryptically or, worse, converges to a
physically meaningless state that becomes a data point in someone's benchmark.

`nqf-lint` reads a declared cluster and runs cheap, deterministic checks that catch
these failure modes **before** the expensive calculation. It is written in Rust, has
no runtime dependencies, and exits with a nonzero code on error so it drops straight
into a `Makefile` or CI.

Every check corresponds to a real bug found in a production pipeline.

## Usage

```
nqf-lint geometry.xyz      # a standard XYZ geometry
nqf-lint cluster.json      # the tool's own cluster spec (charge/spin/ECP)
nqf-lint build_cluster.py  # Python source (conformer determinism)
```

**XYZ** is the universal format every quantum-chemistry package reads and writes,
so the linter runs on real inputs, not only its own JSON. XYZ carries no charge or
spin, so put them on the comment line — either as a bare `<charge> <mult>` pair or
as `charge=.. mult=..` — to enable the electron-count checks; without them, those
checks abstain and the geometry checks still run.

```
3
charge=0 mult=1
O   0.000   0.000   0.000
H   0.757   0.586   0.000
H  -0.757   0.586   0.000
```

The **JSON** spec adds what XYZ cannot express — the ECP list and the metal
oxidation state — for the ECP-aware parity and d¹⁰ spin checks:

```json
{
  "atoms": [
    { "element": "Hg", "x": 0.0, "y": 0.0, "z": 0.0 },
    { "element": "O",  "x": 2.30, "y": 0.0, "z": 0.0 },
    { "element": "H",  "x": 2.60, "y": 0.80, "z": 0.0 },
    { "element": "H",  "x": 2.60, "y": -0.80, "z": 0.0 }
  ],
  "charge": 2,
  "spin_multiplicity": 1,
  "ecp_elements": ["Hg"]
}
```

## What it checks

The geometry checks (bare heteroatom, metal coordination, overlapping atoms) run
on any input — `.xyz` or `.json`. The electron-count checks need charge/spin, and
the ECP/d¹⁰ checks need the extra fields only `.json` carries.

**Geometry + electron count:**

| check | catches |
|---|---|
| **electron parity** | `(N + M)` must be odd. Forgotten ECP core subtraction or a wrong formal charge silently flips parity — the SCF then fails or converges to nonsense. Counts electrons *with* the ECP core removed. |
| **bare heteroatom** | a C/N/O/P/S with no bonding partner in covalent range: waters entered as bare oxygens, dropped hydrogens, or a fragment sliced through a bond. |
| **metal coordination** | a metal center with zero neighbours in range is floating — the fragment was built wrong. |
| **metal spin state** | a closed-shell d¹⁰ cation (Zn²⁺, Hg²⁺, Cu⁺, …) declared open-shell. Only fires when the oxidation state is given and the ion is unambiguous; abstains for ligand-field-dependent ions like Ni²⁺. |
| **overlapping atoms** | two nuclei closer than 0.5 Å — below any real bond (H–H is 0.74 Å). A duplicated atom or a coordinate error that blows up the SCF or double-counts electrons. |

**Python source (`.py`):**

| check | catches |
|---|---|
| **conformer seed** | `EmbedMolecule` / `EmbedMultipleConfs` with no fixed `randomSeed`. The embedding is nondeterministic: the same SMILES gives a different geometry each run, so any descriptor computed on it is irreproducible. This is the bug that turned a published ρ=0.855 into a lottery. |

It **never guesses**. When it cannot verify something — an unknown ECP core, an
unknown element — it abstains with a warning rather than emit a false verdict,
because "cannot verify" is exactly where silent bugs live.

## Design principle

Every check must fail loudly on a known-bad input and pass a known-good one. The test
suite is built from real bugs: the mining cluster whose "waters" were bare oxygens
(`{Co:1, O:7, …}`, zero hydrogens), the Hg cluster whose charge was left at the
all-electron value after switching on an ECP, the metal built with no ligand nearby.
A check that cannot demonstrate both directions does not ship.

## Build

```
cargo test    # 23 tests, each a real bug
cargo build --release
./target/release/nqf-lint examples/bad_mining_cluster.json   # → 4 errors, exit 4
./target/release/nqf-lint examples/good_hg_cluster.json      # → clean, exit 0
./target/release/nqf-lint examples/bad_floating_metal.xyz    # → 2 errors, exit 2
./target/release/nqf-lint examples/water.xyz                 # → clean, exit 0
```

## Limitations

This tool is deliberately narrow and states what it cannot verify.

- **The electron-parity check assumes LANL2DZ core sizes.** With a different ECP
  (SDD, def2-ECP, …) the number of replaced core electrons differs, and the parity
  verdict would be wrong for a *known* element on a non-LANL2DZ ECP. For elements it
  is unsure about, the tool abstains (a warning) rather than guess. If you use another
  ECP, treat parity findings on those atoms as advisory and verify by hand.
- **Bonding is judged by distance, not a real bond-perception model.** A terminal
  metal-oxo/nitrido with an unusually long M=X bond (> 1.75 Å) could be flagged as a
  "bare heteroatom". Rare, but possible.
- **The conformer-seed check is a source-text heuristic, not a Python parser.** It
  recognises the common seeding patterns (`randomSeed=` kwarg, `params.randomSeed =`);
  an exotic way of setting the seed could be missed and reported as unseeded.

The guiding rule is the opposite of the bug it exists to catch: when it cannot be
sure, it says so, instead of emitting a confident wrong answer.

## License

MIT.