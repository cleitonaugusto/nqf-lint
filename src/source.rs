//! Static checks on Python source that build molecular geometries.
//!
//! The single most expensive bug in this project's history was one line:
//!
//! ```python
//! AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())   # no randomSeed
//! ```
//!
//! Without a fixed `randomSeed`, RDKit's embedder is nondeterministic: the same
//! SMILES yields a different 3D geometry every run (2–6 Å RMSD apart). A descriptor
//! computed on that geometry then changes every run — a benchmark ρ that swings
//! from −0.10 to +0.86, and the published value is just the best draw. This check
//! finds unseeded embedder calls before they become an irreproducible number.
//!
//! It is a heuristic over source text, not a full Python parse. It errs toward
//! reporting (an unseeded call it cannot prove is seeded) rather than staying
//! silent — because a false "all clear" is exactly the failure it exists to prevent.

use crate::{Finding, Severity};

const EMBEDDERS: &[&str] = &["EmbedMolecule", "EmbedMultipleConfs"];

/// Byte offset → 1-based line number.
fn line_of(src: &str, idx: usize) -> usize {
    src[..idx].bytes().filter(|&b| b == b'\n').count() + 1
}

/// Given the index of an opening `(`, return the argument text up to the matching
/// `)`, respecting nesting.
fn balanced_args(src: &str, open: usize) -> &str {
    let bytes = src.as_bytes();
    let mut depth = 0i32;
    let mut i = open;
    while i < bytes.len() {
        match bytes[i] {
            b'(' => depth += 1,
            b')' => {
                depth -= 1;
                if depth == 0 {
                    return &src[open + 1..i];
                }
            }
            _ => {}
        }
        i += 1;
    }
    &src[open + 1..]
}

/// Split call arguments on top-level commas (ignoring commas inside nested parens).
fn top_level_split(args: &str) -> Vec<&str> {
    let mut parts = Vec::new();
    let bytes = args.as_bytes();
    let mut depth = 0i32;
    let mut start = 0usize;
    for (i, &b) in bytes.iter().enumerate() {
        match b {
            b'(' | b'[' | b'{' => depth += 1,
            b')' | b']' | b'}' => depth -= 1,
            b',' if depth == 0 => {
                parts.push(&args[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }
    if start <= args.len() {
        parts.push(&args[start..]);
    }
    parts
}

fn is_bare_identifier(s: &str) -> bool {
    let s = s.trim();
    !s.is_empty()
        && s.chars()
            .next()
            .is_some_and(|c| c.is_ascii_alphabetic() || c == '_')
        && s.chars().all(|c| c.is_ascii_alphanumeric() || c == '_')
}

/// From one argument, extract the bare identifier that could be a params object,
/// stripping an optional `name=` keyword prefix. Returns None for constructors,
/// literals, or dotted expressions.
fn param_identifier(arg: &str) -> Option<&str> {
    let arg = arg.trim();
    let value = match arg.split_once('=') {
        Some((_key, v)) => v, // keyword arg: params=<v>
        None => arg,          // positional
    };
    let value = value.trim();
    if is_bare_identifier(value) {
        Some(value)
    } else {
        None
    }
}

/// Parse the integer written immediately after the first `randomSeed` in `args`.
fn seed_value(args: &str) -> Option<i64> {
    let pos = args.find("randomSeed")?;
    let after = &args[pos + "randomSeed".len()..];
    // skip up to the '='
    let after = after.split_once('=').map(|(_, v)| v).unwrap_or(after);
    let mut chars = after.trim_start().chars().peekable();
    let mut num = String::new();
    if let Some(&c) = chars.peek() {
        if c == '-' || c == '+' {
            num.push(c);
            chars.next();
        }
    }
    while let Some(&c) = chars.peek() {
        if c.is_ascii_digit() {
            num.push(c);
            chars.next();
        } else {
            break;
        }
    }
    num.parse().ok()
}

/// Lint Python source for nondeterministic conformer embedding.
pub fn lint_python_source(src: &str) -> Vec<Finding> {
    let mut out = Vec::new();

    // Identifiers that get a `.randomSeed` assigned anywhere in the file are
    // considered seeded param objects.
    let seeded_var = |ident: &str| -> bool { src.contains(&format!("{ident}.randomSeed")) };

    for name in EMBEDDERS {
        let needle = format!("{name}(");
        let mut from = 0usize;
        while let Some(rel) = src[from..].find(&needle) {
            let call_start = from + rel;
            let open = call_start + name.len(); // index of '('
            let args = balanced_args(src, open);
            let line = line_of(src, call_start);

            let finding = classify_call(name, args, line, &seeded_var);
            if let Some(f) = finding {
                out.push(f);
            }
            from = open + 1;
        }
    }

    out.sort_by_key(|f| match f.severity {
        Severity::Error => 0,
        Severity::Warning => 1,
    });
    out
}

fn classify_call(
    name: &str,
    args: &str,
    line: usize,
    seeded_var: &dyn Fn(&str) -> bool,
) -> Option<Finding> {
    if args.contains("randomSeed") {
        return match seed_value(args) {
            Some(-1) => Some(Finding {
                severity: Severity::Warning,
                code: "conformer-explicit-random",
                message: format!(
                    "line {line}: {name}(...) sets randomSeed=-1, which is explicitly \
                     nondeterministic. Use a fixed nonnegative seed for reproducibility."
                ),
            }),
            _ => None, // a fixed seed is set in the call
        };
    }

    // No randomSeed kwarg: is a seeded params object passed by name?
    let args_seeded = top_level_split(args)
        .iter()
        .filter_map(|a| param_identifier(a))
        .any(seeded_var);
    if args_seeded {
        return None;
    }

    Some(Finding {
        severity: Severity::Error,
        code: "conformer-no-seed",
        message: format!(
            "line {line}: {name}(...) has no fixed randomSeed — the embedding is \
             nondeterministic and the 3D geometry (and every descriptor computed on it) \
             changes each run. Set params.randomSeed to a fixed value."
        ),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unseeded_inline_etkdg_is_caught() {
        let src = "    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())\n";
        let f = lint_python_source(src);
        assert_eq!(f.len(), 1);
        assert_eq!(f[0].code, "conformer-no-seed");
        assert_eq!(f[0].severity, Severity::Error);
    }

    #[test]
    fn seeded_params_object_passes() {
        // The real fix applied to 03_co_ni_benchmark_v3.py.
        let src = "\
    _params = AllChem.ETKDGv3()
    _params.randomSeed = 1337
    AllChem.EmbedMolecule(mol, _params)
";
        let f = lint_python_source(src);
        assert!(f.is_empty(), "seeded embed should pass: {f:?}");
    }

    #[test]
    fn seed_kwarg_in_call_passes() {
        let src = "AllChem.EmbedMolecule(mol, randomSeed=42)\n";
        assert!(lint_python_source(src).is_empty());
    }

    #[test]
    fn default_embedder_no_params_is_caught() {
        let src = "AllChem.EmbedMolecule(mol)\n";
        let f = lint_python_source(src);
        assert_eq!(f.len(), 1);
        assert_eq!(f[0].code, "conformer-no-seed");
    }

    #[test]
    fn explicit_minus_one_is_warned_not_errored() {
        let src = "AllChem.EmbedMolecule(mol, randomSeed=-1)\n";
        let f = lint_python_source(src);
        assert_eq!(f.len(), 1);
        assert_eq!(f[0].code, "conformer-explicit-random");
        assert_eq!(f[0].severity, Severity::Warning);
    }

    #[test]
    fn multiple_confs_unseeded_is_caught() {
        let src = "AllChem.EmbedMultipleConfs(mol, numConfs=10, params=AllChem.ETKDGv3())\n";
        let f = lint_python_source(src);
        assert_eq!(f.len(), 1);
        assert_eq!(f[0].code, "conformer-no-seed");
    }

    #[test]
    fn multiline_seeded_params_passes() {
        let src = "\
p = ETKDGv3()
p.randomSeed = 7
EmbedMultipleConfs(
    mol,
    numConfs=20,
    params=p,
)
";
        assert!(lint_python_source(src).is_empty());
    }
}
