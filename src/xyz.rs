//! Parse standard XYZ geometry files into a `ClusterSpec`.
//!
//! XYZ is the universal interchange format for molecular geometries — every
//! quantum-chemistry package reads and writes it — so supporting it lets the
//! linter run on real inputs, not only the tool's own JSON.
//!
//! ```text
//! 3
//! water   (comment line — free text, or charge/spin, see below)
//! O   0.000  0.000  0.000
//! H   0.757  0.586  0.000
//! H  -0.757  0.586  0.000
//! ```
//!
//! XYZ carries no charge or spin. Two common comment-line conventions are
//! recognised so the electron-count checks can run when the information is there:
//!
//! - a bare integer pair `<charge> <multiplicity>` (e.g. `0 1`);
//! - `key=value` tokens `charge=0 mult=1` (also `multiplicity=1`).
//!
//! When neither is present, charge/spin stay unset and those checks abstain
//! rather than assume a value.

use crate::{Atom, ClusterSpec};

pub fn parse_xyz(text: &str) -> Result<ClusterSpec, String> {
    let mut lines = text.lines();

    let count_line = lines.next().ok_or("empty file")?;
    let declared: usize = count_line.trim().parse().map_err(|_| {
        format!(
            "first line must be the atom count, got {:?}",
            count_line.trim()
        )
    })?;

    let comment = lines.next().unwrap_or("");
    let (charge, spin_multiplicity) = parse_charge_spin(comment);

    let mut atoms: Vec<Atom> = Vec::with_capacity(declared);
    for (idx, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split_whitespace();
        let element = fields
            .next()
            .ok_or_else(|| format!("atom line {} has no element", idx + 1))?;
        let x = coord(fields.next(), idx, "x")?;
        let y = coord(fields.next(), idx, "y")?;
        let z = coord(fields.next(), idx, "z")?;
        atoms.push(Atom {
            element: normalize_element(element),
            x,
            y,
            z,
        });
        if atoms.len() == declared {
            break;
        }
    }

    if atoms.len() != declared {
        return Err(format!(
            "header declares {declared} atom(s) but {} were found",
            atoms.len()
        ));
    }

    Ok(ClusterSpec {
        atoms,
        charge,
        spin_multiplicity,
        ecp_elements: Vec::new(),
        metal_oxidation_state: None,
    })
}

fn coord(field: Option<&str>, line_idx: usize, axis: &str) -> Result<f64, String> {
    field
        .ok_or_else(|| {
            format!(
                "atom line {} is missing the {axis} coordinate",
                line_idx + 1
            )
        })?
        .parse()
        .map_err(|_| {
            format!(
                "atom line {}: {axis} coordinate is not a number",
                line_idx + 1
            )
        })
}

/// "hg" / "HG" / "Hg" → "Hg" so it matches the reference tables.
fn normalize_element(s: &str) -> String {
    let mut chars = s.chars();
    match chars.next() {
        Some(first) => {
            first.to_ascii_uppercase().to_string() + &chars.as_str().to_ascii_lowercase()
        }
        None => String::new(),
    }
}

fn parse_charge_spin(comment: &str) -> (Option<i64>, Option<u32>) {
    let lower = comment.to_lowercase();
    let charge = extract_kv(&lower, "charge=");
    let mult = extract_kv(&lower, "multiplicity=").or_else(|| extract_kv(&lower, "mult="));
    if charge.is_some() || mult.is_some() {
        return (charge, mult.and_then(|m| u32::try_from(m).ok()));
    }
    // Bare "<charge> <multiplicity>" convention (exactly two integers).
    let toks: Vec<&str> = comment.split_whitespace().collect();
    if toks.len() == 2 {
        if let (Ok(c), Ok(m)) = (toks[0].parse::<i64>(), toks[1].parse::<u32>()) {
            return (Some(c), Some(m));
        }
    }
    (None, None)
}

fn extract_kv(lower: &str, key: &str) -> Option<i64> {
    let pos = lower.find(key)?;
    let after = lower[pos + key.len()..].trim_start();
    let tok: String = after
        .chars()
        .take_while(|c| c.is_ascii_digit() || *c == '-' || *c == '+')
        .collect();
    tok.parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_plain_water() {
        let xyz = "3\nwater\nO 0.0 0.0 0.0\nH 0.757 0.586 0.0\nH -0.757 0.586 0.0\n";
        let spec = parse_xyz(xyz).unwrap();
        assert_eq!(spec.atoms.len(), 3);
        assert_eq!(spec.atoms[0].element, "O");
        assert_eq!(spec.charge, None); // no charge/spin in a bare comment
        assert_eq!(spec.spin_multiplicity, None);
    }

    #[test]
    fn reads_bare_charge_mult_pair() {
        let xyz = "1\n0 1\nNe 0.0 0.0 0.0\n";
        let spec = parse_xyz(xyz).unwrap();
        assert_eq!(spec.charge, Some(0));
        assert_eq!(spec.spin_multiplicity, Some(1));
    }

    #[test]
    fn reads_key_value_charge_mult() {
        let xyz = "1\ncharge=-2 multiplicity=3\nO 0.0 0.0 0.0\n";
        let spec = parse_xyz(xyz).unwrap();
        assert_eq!(spec.charge, Some(-2));
        assert_eq!(spec.spin_multiplicity, Some(3));
    }

    #[test]
    fn normalizes_element_case() {
        let spec = parse_xyz("1\n\nHG 0 0 0\n").unwrap();
        assert_eq!(spec.atoms[0].element, "Hg");
    }

    #[test]
    fn atom_count_mismatch_is_error() {
        assert!(parse_xyz("3\nc\nC 0 0 0\n").is_err());
    }

    #[test]
    fn non_numeric_count_is_error() {
        assert!(parse_xyz("water\n\nO 0 0 0\n").is_err());
    }
}
