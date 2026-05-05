// recombination_analyzer.rs - DEPRECATED ALIAS for recombination_candidate_analyzer.
// Kept as a deprecation shim that forwards every argument to the new binary so existing
// scripts keep working.

use std::process::{exit, Command};

fn main() {
    eprintln!(
        "⚠️  Note: 'recombination_analyzer' is a deprecated alias for \
         'recombination_candidate_analyzer'. The legacy name is still accepted; \
         please update your scripts."
    );

    // Locate the renamed binary next to this one (same target/release directory).
    let mut new_path =
        std::env::current_exe().expect("could not determine path of current executable");
    new_path.set_file_name("recombination_candidate_analyzer");

    let args: Vec<String> = std::env::args().skip(1).collect();
    let status = Command::new(&new_path)
        .args(&args)
        .status()
        .unwrap_or_else(|e| {
            eprintln!("error: failed to invoke '{}': {}", new_path.display(), e);
            exit(127);
        });

    exit(status.code().unwrap_or(1));
}
