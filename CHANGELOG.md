# Changelog

All notable changes to cgDist are documented in this file. The format is based
on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html) with the
relaxed pre-1.0 convention (breaking changes are allowed in 0.x patch releases
until the API stabilizes).

## [0.1.2] — 2026-05-26

Maintenance release focused on installation, documentation, and
reproducibility — making cgDist easier to install on current toolchains and
its validation suite easier to reproduce from a fresh clone. The distance
algorithms and their numerical results are unchanged from 0.1.1.

### Changed

- Clarified the minimum supported Rust version and install guidance
  (latest stable Rust; `cargo install cgdist --locked`).
- Documented the recombination-candidate workflow through its supported
  path — an enriched cache (`--enrich-lengths`) analysed by
  `recombination_candidate_analyzer`; cache inspection via `cgdist --inspector`.
- Hardened CLI input validation, error messages, and detailed-alignment
  output (`--save-alignments`).

### Added

- A self-contained validation suite that runs in CI on every push — covering
  distance-mode correctness, cache consistency, the recombination-candidate
  workflow, filtering / missing-data / output-format behaviour, and a
  smoke-test over every CLI argument.

## [0.1.1] — 2026-05-05

### Added

- **CLI flag `--candidate-recombination-log`** (canonical name) for the
  per-locus mutation-density flagging log. The previous name
  `--recombination-log` is kept as a deprecated alias and prints a
  deprecation warning when used.
- **CLI flag `--candidate-recombination-threshold`** (canonical name).
  The previous name `--recombination-threshold` is kept as a deprecated
  alias with a deprecation warning.
- **Opt-in flag `--hamming-fallback`** to enable the +1 Hamming fallback
  in SNPs-only mode when only InDel differences exist between two
  alleles. The fallback is now disabled by default.
- **Binary `recombination_candidate_analyzer`** (canonical name)
  replacing `recombination_analyzer`, which is kept as a deprecation
  shim that forwards every argument to the new binary.
- **`examples/cgdist-config.toml`** — canonical TOML configuration
  example whose flat layout matches the parser.
- **MSRV declaration** `rust-version = "1.70"` in `Cargo.toml`, plus
  `readme`, `homepage`, `documentation`, `keywords`, `categories`, and
  `exclude` metadata in preparation for a future crates.io publication.
- **`validation_test/profiles/test_profiles_crc32.tsv`** committed to
  the repository so the validation suite is reproducible from a fresh
  clone.
- **GitHub Actions workflow** `.github/workflows/ci-and-docker.yml` that
  runs `cargo fmt --check`, `cargo clippy -D warnings`, `cargo test`,
  the four-mode validation smoke test, and (on master/`v*` tags) builds
  and pushes multi-arch Docker images to GHCR.
- **`CHANGELOG.md`** (this file).
- **Zenodo concept DOI**, **bioRxiv DOI**, and **MSRV** badges in
  `README.md`.

### Changed

- **Default `--threads` is now `1`** (previously implementation-defined,
  typically all available cores). Users who want parallelism must opt in
  explicitly with `--threads N`. This is a behaviour change for existing
  scripts.
- **`README.md`** rewritten:
  - The recombination section is reframed as
    "Recombination-Candidate Flagging" with an explicit disclaimer that
    cgDist is not a recombination detector and that confirmation
    requires phylogeny-aware tools (Gubbins, ClonalFrameML, fastGEAR).
  - The configuration file is documented as optional.
  - A new "CLI vs TOML precedence" subsection: the command-line value
    wins when both are provided.
  - The inline TOML example is rewritten to match the parser's flat
    layout.
  - The install path is updated to
    `cargo install --git ... --tag v0.1.1 cgdist`.
- **`examples/hamming-config.toml`**: the typo `hashery_type` is
  corrected to `hasher_type`.
- **CITATION.cff** version bumped to `0.1.1`, release date `2026-05-05`.
- A wide round of `cargo clippy` autofixes is applied across the
  codebase (idiomatic iterator usage, removal of unnecessary `unwrap`
  after `is_some`, `format!` argument inlining, etc.); no behaviour
  change.

### Removed

- The "Star clustering analysis" feature bullet from `README.md` — the
  feature was not implemented in code; the bullet was an aspirational
  overclaim.

### Deprecated

- `--recombination-log` → use `--candidate-recombination-log`.
- `--recombination-threshold` → use `--candidate-recombination-threshold`.
- Binary `recombination_analyzer` → use `recombination_candidate_analyzer`.

All three legacy names continue to work and print a deprecation notice
on use.

## [0.1.0] — 2025-12-23

Initial public release accompanying the bioRxiv preprint
(DOI: [10.1101/2025.10.16.682749](https://doi.org/10.1101/2025.10.16.682749)).

[0.1.1]: https://github.com/genpat-it/cgDist/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/genpat-it/cgDist/releases/tag/v0.1.0
