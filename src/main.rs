//! nqf-lint CLI — lint a quantum-chemistry cluster setup before running it.
//!
//!     nqf-lint <cluster.json>
//!
//! Exit code = number of errors, so it drops straight into a Makefile or CI:
//! a broken cluster fails the build instead of quietly becoming a data point.

use std::process::exit;

use nqf_lint::{error_count, lint, lint_python_source, parse_xyz, ClusterSpec, Finding, Severity};

fn main() {
    let path = match std::env::args().nth(1) {
        Some(p) => p,
        None => {
            eprintln!("usage: nqf-lint <file>");
            eprintln!("  <cluster.json>  lints a quantum-chemistry cluster spec");
            eprintln!("  <geometry.xyz>  lints a standard XYZ geometry");
            eprintln!(
                "  <script.py>     lints Python source for nondeterministic conformer embedding"
            );
            exit(2);
        }
    };

    let text = match std::fs::read_to_string(&path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("nqf-lint: cannot read {path}: {e}");
            exit(2);
        }
    };

    // Dispatch by extension: .py → source checks, .xyz → geometry, else JSON cluster.
    let findings: Vec<Finding> = if path.ends_with(".py") {
        lint_python_source(&text)
    } else if path.ends_with(".xyz") {
        match parse_xyz(&text) {
            Ok(spec) => lint(&spec),
            Err(e) => {
                eprintln!("nqf-lint: {path} is not a valid XYZ file: {e}");
                exit(2);
            }
        }
    } else {
        let spec: ClusterSpec = match serde_json::from_str(&text) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("nqf-lint: {path} is not a valid cluster spec: {e}");
                exit(2);
            }
        };
        lint(&spec)
    };

    if findings.is_empty() {
        println!("✓ {path}: no issues found");
        exit(0);
    }

    for f in &findings {
        let tag = match f.severity {
            Severity::Error => "ERROR",
            Severity::Warning => "warn ",
        };
        println!("{tag} [{}] {}", f.code, f.message);
    }

    let errors = error_count(&findings);
    let warns = findings.len() - errors;
    println!("\n{path}: {errors} error(s), {warns} warning(s)");
    exit(errors.min(125) as i32);
}
