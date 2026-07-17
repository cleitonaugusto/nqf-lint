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
nqf-lint cluster.json
```

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

**Cluster spec (`.json`):**

| check | catches |
|---|---|
| **electron parity** | `(N + M)` must be odd. Forgotten ECP core subtraction or a wrong formal charge silently flips parity — the SCF then fails or converges to nonsense. Counts electrons *with* the ECP core removed. |
| **bare heteroatom** | a C/N/O/P/S with no bonding partner in covalent range: waters entered as bare oxygens, dropped hydrogens, or a fragment sliced through a bond. |
| **metal coordination** | a metal center with zero neighbours in range is floating — the fragment was built wrong. |
| **metal spin state** | a closed-shell d¹⁰ cation (Zn²⁺, Hg²⁺, Cu⁺, …) declared open-shell. Only fires when the oxidation state is given and the ion is unambiguous; abstains for ligand-field-dependent ions like Ni²⁺. |

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
cargo test    # 6 tests, each a real bug
cargo build --release
./target/release/nqf-lint examples/bad_mining_cluster.json   # → 4 errors, exit 4
./target/release/nqf-lint examples/good_hg_cluster.json      # → clean, exit 0
```

## License

MIT.