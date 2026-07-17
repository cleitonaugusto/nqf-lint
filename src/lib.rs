//! nqf-lint — pre-flight linter for quantum-chemistry cluster setups.
//!
//! A quantum-chemistry calculation can run for hours and return a number that
//! *looks* fine while the input was silently broken: an oxygen with no hydrogens,
//! an electron count made impossible by a forgotten ECP core, a metal floating
//! with no coordination. The SCF either fails cryptically or — worse — converges
//! to a physically meaningless state that becomes a data point in a benchmark.
//!
//! This crate reads a declared cluster (atoms + charge + spin + ECP list) and
//! runs cheap, deterministic checks that catch these failure modes *before* the
//! expensive calculation. Every check here corresponds to a real bug found in a
//! production quantum-chemistry pipeline.
//!
//! It never guesses. When it cannot verify something (e.g. an unknown ECP core),
//! it says so — because "cannot verify" is exactly where silent bugs live.

use serde::Deserialize;

pub mod source;
pub use source::lint_python_source;

#[derive(Debug, Clone, Deserialize)]
pub struct Atom {
    pub element: String,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Clone, Deserialize)]
pub struct ClusterSpec {
    pub atoms: Vec<Atom>,
    /// Formal total charge of the cluster.
    pub charge: i64,
    /// Spin multiplicity M = 2S+1 (1 = singlet, 2 = doublet, 3 = triplet, ...).
    pub spin_multiplicity: u32,
    /// Elements carrying an effective core potential (ECP), e.g. ["Hg"].
    /// Their core electrons are replaced by the pseudopotential and must be
    /// subtracted from the electron count. Omit for all-electron calculations.
    #[serde(default)]
    pub ecp_elements: Vec<String>,
    /// Optional formal oxidation state of the (single) metal center. When given,
    /// enables the spin-state check for closed-shell d¹⁰ cations. Omit to skip it.
    #[serde(default)]
    pub metal_oxidation_state: Option<i64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Severity {
    /// The setup is physically impossible or almost certainly wrong.
    Error,
    /// Suspicious; may be intentional. Cannot be auto-verified.
    Warning,
}

#[derive(Debug, Clone)]
pub struct Finding {
    pub severity: Severity,
    pub code: &'static str,
    pub message: String,
}

impl Finding {
    fn error(code: &'static str, message: String) -> Self {
        Finding {
            severity: Severity::Error,
            code,
            message,
        }
    }
    fn warn(code: &'static str, message: String) -> Self {
        Finding {
            severity: Severity::Warning,
            code,
            message,
        }
    }
}

// ── Reference data ───────────────────────────────────────────────────────────

/// Atomic number for the elements that appear in coordination clusters.
/// Deliberately partial: an unknown element makes the electron-count check
/// abstain (Warning) rather than guess.
fn atomic_number(el: &str) -> Option<i64> {
    Some(match el {
        "H" => 1,
        "He" => 2,
        "Li" => 3,
        "B" => 5,
        "C" => 6,
        "N" => 7,
        "O" => 8,
        "F" => 9,
        "Na" => 11,
        "Mg" => 12,
        "Al" => 13,
        "Si" => 14,
        "P" => 15,
        "S" => 16,
        "Cl" => 17,
        "K" => 19,
        "Ca" => 20,
        "Fe" => 26,
        "Co" => 27,
        "Ni" => 28,
        "Cu" => 29,
        "Zn" => 30,
        "Br" => 35,
        "Ag" => 47,
        "Cd" => 48,
        "I" => 53,
        "Au" => 79,
        "Hg" => 80,
        "Pb" => 82,
        _ => return None,
    })
}

/// Core electrons replaced by the LANL2DZ ECP, for the heavy elements where the
/// value is well established. Returns None when we are not confident — the caller
/// must then abstain from the parity check rather than produce a false verdict.
fn lanl2dz_ecp_core(el: &str) -> Option<i64> {
    Some(match el {
        // Third-row heavy elements: [Xe]4f14 = 60 core electrons.
        "Hg" | "Au" | "Pb" | "Tl" | "Bi" | "Pt" | "Ir" => 60,
        // Second-row transition / post-transition: [Ar]3d10 = 28.
        "Ag" | "Cd" | "I" => 28,
        _ => return None,
    })
}

const METALS: &[&str] = &[
    "Fe", "Co", "Ni", "Cu", "Zn", "Ag", "Cd", "Au", "Hg", "Pb", "Pt", "Ir", "Mn",
];

fn is_metal(el: &str) -> bool {
    METALS.contains(&el)
}

fn dist(a: &Atom, b: &Atom) -> f64 {
    ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt()
}

/// Nearest-neighbour distance to any *other* atom (∞ if the atom is alone).
fn nearest_neighbor(atoms: &[Atom], i: usize) -> f64 {
    atoms
        .iter()
        .enumerate()
        .filter(|(j, _)| *j != i)
        .map(|(_, other)| dist(&atoms[i], other))
        .fold(f64::INFINITY, f64::min)
}

// ── The checks ───────────────────────────────────────────────────────────────

/// C1 — Electron count vs spin multiplicity.
///
/// For N electrons and multiplicity M = 2S+1, the number of unpaired electrons is
/// M−1, so N and (M−1) must share parity ⇔ (N + M) is odd. An even (N + M) is a
/// physically impossible state. This is the classic ECP/charge footgun: forget to
/// subtract the pseudopotential core (or set the wrong formal charge) and the
/// electron count silently flips parity — the SCF then fails or converges to
/// nonsense. Here the count is computed *with* the ECP core removed.
fn check_electron_parity(spec: &ClusterSpec, out: &mut Vec<Finding>) {
    let mut z_sum: i64 = 0;
    let mut unknown: Vec<String> = Vec::new();
    for a in &spec.atoms {
        match atomic_number(&a.element) {
            Some(z) => z_sum += z,
            None => unknown.push(a.element.clone()),
        }
    }
    if !unknown.is_empty() {
        unknown.sort();
        unknown.dedup();
        out.push(Finding::warn(
            "electron-parity-abstain",
            format!(
                "electron count not verified: unknown atomic number for {}",
                unknown.join(", ")
            ),
        ));
        return;
    }

    // Subtract ECP core electrons per atom that carries an ECP.
    let ecp_set: std::collections::HashSet<&str> =
        spec.ecp_elements.iter().map(String::as_str).collect();
    let mut ecp_core_sum: i64 = 0;
    let mut unverifiable: Vec<String> = Vec::new();
    for a in &spec.atoms {
        if ecp_set.contains(a.element.as_str()) {
            match lanl2dz_ecp_core(&a.element) {
                Some(c) => ecp_core_sum += c,
                None => unverifiable.push(a.element.clone()),
            }
        }
    }
    if !unverifiable.is_empty() {
        unverifiable.sort();
        unverifiable.dedup();
        out.push(Finding::warn(
            "electron-parity-abstain",
            format!(
                "electron parity not verified: unknown ECP core for {} — declare it or verify by hand",
                unverifiable.join(", ")
            ),
        ));
        return;
    }

    let n_eff = z_sum - spec.charge - ecp_core_sum;
    if n_eff < 0 {
        out.push(Finding::error(
            "electron-count-negative",
            format!("electron count is negative ({n_eff}): charge or ECP setup is wrong"),
        ));
        return;
    }
    let m = spec.spin_multiplicity as i64;
    if m < 1 {
        out.push(Finding::error(
            "spin-multiplicity-invalid",
            format!("spin multiplicity must be ≥ 1, got {m}"),
        ));
        return;
    }
    if (n_eff + m) % 2 == 0 {
        out.push(Finding::error(
            "electron-parity",
            format!(
                "impossible state: {n_eff} electrons with multiplicity {m}. \
                 N and (M−1) must share parity — (N+M) must be odd, here it is even. \
                 Most common cause: forgotten ECP core subtraction or wrong formal charge."
            ),
        ));
    }
}

/// C2 — Bare heteroatom (missing hydrogens / floating atom).
///
/// A C, N, O, P or S atom with no bonding partner within a covalent distance is
/// almost always a bug: waters entered as bare oxygens, hydrogens dropped during
/// fragment cutting, or a fragment sliced through a bond leaving a dangling atom.
/// A metal–oxygen dative bond sits at ~2.0–2.5 Å, so a water oxygen whose *only*
/// near contact is the metal (i.e. its hydrogens are missing) is correctly caught.
fn check_bare_heteroatom(spec: &ClusterSpec, out: &mut Vec<Finding>) {
    for (i, a) in spec.atoms.iter().enumerate() {
        let threshold = match a.element.as_str() {
            "P" | "S" => 1.95, // longer single bonds to P/S
            "C" | "N" | "O" => 1.75,
            _ => continue, // metals, H, halides: skip
        };
        let nn = nearest_neighbor(&spec.atoms, i);
        if nn > threshold {
            out.push(Finding::error(
                "bare-heteroatom",
                format!(
                    "atom {i} ({}) has no bonding partner within {threshold:.2} Å \
                     (nearest is {nn:.2} Å). Likely a bare oxygen / missing hydrogens / \
                     a fragment cut through a bond.",
                    a.element
                ),
            ));
        }
    }
}

/// C3 — Metal coordination sanity.
///
/// The metal center must actually be coordinated. Zero neighbours within 2.7 Å
/// means the metal is floating (fragment built wrong); one neighbour is suspicious
/// for a coordination cluster and worth a look.
fn check_metal_coordination(spec: &ClusterSpec, out: &mut Vec<Finding>) {
    for (i, a) in spec.atoms.iter().enumerate() {
        if !is_metal(&a.element) {
            continue;
        }
        let coord = spec
            .atoms
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .filter(|(_, o)| dist(a, o) <= 2.7)
            .count();
        match coord {
            0 => out.push(Finding::error(
                "metal-uncoordinated",
                format!(
                    "metal {} (atom {i}) has no neighbour within 2.7 Å — it is floating",
                    a.element
                ),
            )),
            1 => out.push(Finding::warn(
                "metal-low-coordination",
                format!(
                    "metal {} (atom {i}) has coordination 1 — verify this is intended",
                    a.element
                ),
            )),
            _ => {}
        }
    }
}

/// A known closed-shell d¹⁰ cation must have a singlet ground state (M = 1).
/// Returns Some(true) only for ions where this is unambiguous. Everything else
/// returns None — the tool abstains rather than guess, because spin state for
/// open-shell ions (e.g. Ni²⁺ d⁸: singlet square-planar vs triplet octahedral)
/// depends on the ligand field and cannot be decided from the formula alone.
fn requires_singlet(element: &str, oxidation_state: i64) -> Option<bool> {
    match (element, oxidation_state) {
        ("Zn", 2) | ("Cd", 2) | ("Hg", 2) => Some(true), // d¹⁰
        ("Cu", 1) | ("Ag", 1) | ("Au", 1) => Some(true), // d¹⁰
        _ => None,
    }
}

/// C4 — Spin state of a closed-shell d¹⁰ metal.
///
/// A d¹⁰ cation (Zn²⁺, Hg²⁺, Cu⁺, …) has a filled d shell and a singlet ground
/// state; declaring it open-shell (multiplicity > 1) is a physical error and a
/// real class of bug — e.g. running a metal as a triplet when its d-count forbids
/// it. Only fires when the oxidation state is declared and the ion is unambiguous.
fn check_metal_spin_state(spec: &ClusterSpec, out: &mut Vec<Finding>) {
    let Some(ox) = spec.metal_oxidation_state else {
        return; // not declared → abstain
    };
    let metals: Vec<&Atom> = spec.atoms.iter().filter(|a| is_metal(&a.element)).collect();
    if metals.len() != 1 {
        return; // ambiguous which metal the oxidation state refers to → abstain
    }
    let metal = metals[0];
    if requires_singlet(&metal.element, ox) == Some(true) && spec.spin_multiplicity != 1 {
        out.push(Finding::error(
            "metal-spin-state",
            format!(
                "{}{:+} is a closed-shell d¹⁰ cation and must be a singlet (M=1), \
                 but multiplicity {} was declared.",
                metal.element, ox, spec.spin_multiplicity
            ),
        ));
    }
}

/// C5 — Overlapping atoms.
///
/// Two nuclei closer than 0.5 Å are not a chemical bond — the shortest real bond,
/// H–H in dihydrogen, is 0.74 Å. A sub-0.5-Å pair is a duplicated atom or a
/// coordinate error (e.g. a hydrogen written onto its parent heavy atom, or a
/// fragment merged twice). Left in, it blows up the SCF or double-counts electrons.
fn check_overlapping_atoms(spec: &ClusterSpec, out: &mut Vec<Finding>) {
    const OVERLAP: f64 = 0.5;
    let n = spec.atoms.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let d = dist(&spec.atoms[i], &spec.atoms[j]);
            if d < OVERLAP {
                out.push(Finding::error(
                    "overlapping-atoms",
                    format!(
                        "atoms {i} ({}) and {j} ({}) are {d:.3} Å apart — below any real \
                         bond (H–H is 0.74 Å). Likely a duplicated atom or a coordinate error.",
                        spec.atoms[i].element, spec.atoms[j].element
                    ),
                ));
            }
        }
    }
}

/// Run every check and return all findings, most severe first.
pub fn lint(spec: &ClusterSpec) -> Vec<Finding> {
    let mut out = Vec::new();
    check_electron_parity(spec, &mut out);
    check_bare_heteroatom(spec, &mut out);
    check_metal_coordination(spec, &mut out);
    check_metal_spin_state(spec, &mut out);
    check_overlapping_atoms(spec, &mut out);
    out.sort_by_key(|f| match f.severity {
        Severity::Error => 0,
        Severity::Warning => 1,
    });
    out
}

pub fn error_count(findings: &[Finding]) -> usize {
    findings
        .iter()
        .filter(|f| f.severity == Severity::Error)
        .count()
}

// ── Tests: every case is a REAL bug from a production pipeline ────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(el: &str, x: f64, y: f64, z: f64) -> Atom {
        Atom {
            element: el.into(),
            x,
            y,
            z,
        }
    }

    /// A well-formed Hg²⁺ cluster: Hg + two real waters (O with two H each),
    /// charge +2, singlet, Hg carrying the LANL2DZ ECP. Must pass clean.
    fn good_hg_cluster() -> ClusterSpec {
        ClusterSpec {
            atoms: vec![
                atom("Hg", 0.0, 0.0, 0.0),
                // water 1
                atom("O", 2.30, 0.0, 0.0),
                atom("H", 2.60, 0.80, 0.0),
                atom("H", 2.60, -0.80, 0.0),
                // water 2
                atom("O", -2.30, 0.0, 0.0),
                atom("H", -2.60, 0.80, 0.0),
                atom("H", -2.60, -0.80, 0.0),
            ],
            charge: 2,
            spin_multiplicity: 1,
            ecp_elements: vec!["Hg".into()],
            metal_oxidation_state: Some(2),
        }
    }

    #[test]
    fn good_cluster_is_clean() {
        let f = lint(&good_hg_cluster());
        assert_eq!(error_count(&f), 0, "well-formed Hg cluster flagged: {f:?}");
    }

    #[test]
    fn good_cluster_parity_holds() {
        // Hg 80 − ecp 60 = 20; two waters = 2*(8+1+1)=20; total N = 40 − charge 2 = 38.
        // M = 1. (38 + 1) = 39, odd ⇒ valid.
        let f = lint(&good_hg_cluster());
        assert!(!f.iter().any(|x| x.code == "electron-parity"));
    }

    /// The mining-cluster bug: waters entered as BARE oxygens, every hydrogen
    /// dropped. Real observed composition was {Co:1, O:7, P:1, C:3}, zero H.
    #[test]
    fn bare_oxygen_waters_are_caught() {
        let spec = ClusterSpec {
            atoms: vec![
                atom("Co", 0.0, 0.0, 0.0),
                atom("O", 2.10, 0.0, 0.0),  // bare: only the metal is near
                atom("O", -2.10, 0.0, 0.0), // bare
                atom("O", 0.0, 2.10, 0.0),  // bare
            ],
            charge: 2,
            spin_multiplicity: 4, // Co²⁺ d⁷ high spin
            ecp_elements: vec![],
            metal_oxidation_state: Some(2), // Co²⁺ is open-shell → spin check abstains
        };
        let f = lint(&spec);
        let bare = f.iter().filter(|x| x.code == "bare-heteroatom").count();
        assert_eq!(bare, 3, "should flag all three bare oxygens: {f:?}");
    }

    /// The ECP/charge parity bug: Hg cluster where the formal charge was left at
    /// the all-electron value and the ECP core was forgotten, flipping parity.
    /// Here we force an impossible (N+M) even.
    #[test]
    fn impossible_electron_parity_is_caught() {
        // Hg(H2O) : Hg 80−60=20, water 10 → N=30 at charge 0. M=2 (declared).
        // (30 + 2) = 32, even ⇒ impossible for a doublet.
        let spec = ClusterSpec {
            atoms: vec![
                atom("Hg", 0.0, 0.0, 0.0),
                atom("O", 2.30, 0.0, 0.0),
                atom("H", 2.60, 0.80, 0.0),
                atom("H", 2.60, -0.80, 0.0),
            ],
            charge: 0,
            spin_multiplicity: 2,
            ecp_elements: vec!["Hg".into()],
            metal_oxidation_state: None,
        };
        let f = lint(&spec);
        assert!(
            f.iter()
                .any(|x| x.code == "electron-parity" && x.severity == Severity::Error),
            "impossible parity not caught: {f:?}"
        );
    }

    /// A floating metal: the fragment was built without any ligand near the metal.
    #[test]
    fn floating_metal_is_caught() {
        let spec = ClusterSpec {
            atoms: vec![
                atom("Zn", 0.0, 0.0, 0.0),
                atom("O", 5.0, 0.0, 0.0), // far away — not coordinating
                atom("H", 5.3, 0.8, 0.0),
                atom("H", 5.3, -0.8, 0.0),
            ],
            charge: 2,
            spin_multiplicity: 1,
            ecp_elements: vec![],
            metal_oxidation_state: None,
        };
        let f = lint(&spec);
        assert!(
            f.iter().any(|x| x.code == "metal-uncoordinated"),
            "floating metal not caught: {f:?}"
        );
    }

    /// Unknown ECP core ⇒ the tool must ABSTAIN (warn), never emit a false verdict.
    #[test]
    fn unknown_ecp_abstains_rather_than_guessing() {
        let spec = ClusterSpec {
            atoms: vec![
                atom("Gd", 0.0, 0.0, 0.0),
                atom("O", 2.3, 0.0, 0.0),
                atom("H", 2.6, 0.8, 0.0),
                atom("H", 2.6, -0.8, 0.0),
            ],
            charge: 3,
            spin_multiplicity: 8,
            ecp_elements: vec!["Gd".into()],
            metal_oxidation_state: None,
        };
        let f = lint(&spec);
        // Gd's atomic number is unknown to our table → abstain on parity, no false Error.
        assert!(!f.iter().any(|x| x.code == "electron-parity"));
        assert!(f.iter().any(|x| x.code == "electron-parity-abstain"));
    }

    /// A d¹⁰ metal declared open-shell must be flagged.
    #[test]
    fn d10_metal_as_nonsinglet_is_caught() {
        let spec = ClusterSpec {
            atoms: vec![
                atom("Zn", 0.0, 0.0, 0.0),
                atom("O", 2.10, 0.0, 0.0),
                atom("H", 2.40, 0.80, 0.0),
                atom("H", 2.40, -0.80, 0.0),
            ],
            charge: 2,
            spin_multiplicity: 3, // wrong: Zn²⁺ d¹⁰ is a singlet
            ecp_elements: vec![],
            metal_oxidation_state: Some(2),
        };
        let f = lint(&spec);
        assert!(
            f.iter()
                .any(|x| x.code == "metal-spin-state" && x.severity == Severity::Error),
            "Zn²⁺ triplet not caught: {f:?}"
        );
    }

    /// Two atoms written to (nearly) the same coordinates — a duplicated atom.
    #[test]
    fn overlapping_atoms_are_caught() {
        let spec = ClusterSpec {
            atoms: vec![
                atom("Zn", 0.0, 0.0, 0.0),
                atom("O", 2.10, 0.0, 0.0),
                atom("O", 2.11, 0.0, 0.0), // 0.01 Å from the previous O — a duplicate
                atom("H", 2.40, 0.80, 0.0),
                atom("H", 2.40, -0.80, 0.0),
            ],
            charge: 2,
            spin_multiplicity: 1,
            ecp_elements: vec![],
            metal_oxidation_state: None,
        };
        let f = lint(&spec);
        assert!(
            f.iter()
                .any(|x| x.code == "overlapping-atoms" && x.severity == Severity::Error),
            "overlapping O atoms not caught: {f:?}"
        );
    }

    /// A normal O–H distance (0.96 Å) must NOT trip the overlap check.
    #[test]
    fn normal_bond_length_is_not_overlap() {
        let spec = ClusterSpec {
            atoms: vec![atom("O", 0.0, 0.0, 0.0), atom("H", 0.96, 0.0, 0.0)],
            charge: 0,
            spin_multiplicity: 1,
            ecp_elements: vec![],
            metal_oxidation_state: None,
        };
        let f = lint(&spec);
        assert!(!f.iter().any(|x| x.code == "overlapping-atoms"));
    }

    /// An open-shell / ambiguous ion (Ni²⁺ d⁸) must NOT be flagged on spin —
    /// its state depends on the ligand field.
    #[test]
    fn open_shell_metal_spin_abstains() {
        let spec = ClusterSpec {
            atoms: vec![
                atom("Ni", 0.0, 0.0, 0.0),
                atom("O", 2.05, 0.0, 0.0),
                atom("H", 2.35, 0.80, 0.0),
                atom("H", 2.35, -0.80, 0.0),
            ],
            charge: 2,
            spin_multiplicity: 1, // could be right (square-planar) — must not error
            ecp_elements: vec![],
            metal_oxidation_state: Some(2),
        };
        let f = lint(&spec);
        assert!(!f.iter().any(|x| x.code == "metal-spin-state"));
    }
}
