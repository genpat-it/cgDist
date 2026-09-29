# Changelog

All notable changes to cgDist are documented in this file. The format is based
on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html) with the
relaxed pre-1.0 convention (breaking changes are allowed in 0.x patch releases
until the API stabilizes).

## [Unreleased]

### Changed

- Alignments are much faster with bit-identical results. A pair is first
  aligned with a certified banded aligner (`src/core/banded.rs`): it fills
  only a diagonal band of the DP matrix, reproduces parasail's recurrences,
  tie-breaking and traceback, and returns a result only when a band
  certificate proves it equals the full-matrix result. Otherwise the pair is
  aligned by parasail as before, now with the scan kernel at 16-bit
  precision (escalating to 32/64-bit on saturation) and per-thread aligners
  instead of the slower striped kernel. The proof, the invariants and the
  verification (exhaustive over all short sequences, 739,554 real allele
  pairs, adversarial fuzzing: no differences) are in
  `docs/BANDED_ALIGNMENT_PROOF.md`. Distance matrices and cache statistics
  are identical to 0.1.4.

- A cache no longer goes through a separate enrichment pass after saving
  when every entry already has its allele lengths: newly aligned pairs record
  both lengths at alignment time, which avoids re-reading the whole schema.
  Cold runs with `--cache-file` (16 threads): 72 s -> 2.1 s
  (L. monocytogenes, 300 samples), 153 s -> 3.9 s (S. enterica, 120 samples);
  the alignment step itself is 120-140x faster per pair.

- `--save-alignments` uses the certified banded alignment too: it writes the
  same gapped strings as parasail, and the file content is byte-identical to
  0.1.4. Rows are now written as they are produced instead of being held in
  memory until the end: peak memory 2.9 GB -> 0.5 GB (L. monocytogenes, 300
  samples) and 5.5 GB -> 1.3 GB (S. enterica, 120 samples); run time
  72 s -> 2.4 s and 153 s -> 4.7 s.

### Fixed

- `--save-alignments` was silently ignored with `--cache-only`; the file is
  now written in that mode too.
- Cache enrichment looked up allele lengths in one CRC32 map for the whole
  schema. CRC32 values collide across loci of large schemas (906 colliding
  CRCs with different lengths in the S. enterica schema), so some entries
  received the length of another locus' allele, depending on file order.
  Lengths are now always taken per locus (215 affected entries in the 120
  sample S. enterica test set; distances were never affected, only
  recombination densities).

### Added

- `--save-cigar <file>`: one compact row per aligned pair (locus, hash1,
  hash2, CIGAR, snps, indel_events, indel_bases, alignment_score). The
  extended CIGAR (`=`, `X`, `I`, `D`) gives the position of every SNP and
  InDel, about 35x smaller than `--save-alignments`. It equals parasail's
  `get_cigar` (verified exhaustively on short sequences and on 739,554 real
  pairs).
- `cgdist-diff`: lists the SNPs, insertions and deletions between two
  alleles of a locus, given their hashes, with positions, bases and CIGAR.
  The alleles are aligned exactly as cgdist does, and the counts can be
  checked against a cache entry (`--cache-file`). Differences are annotated
  at protein level (new modules `core::protein` and `core::codon_tables`):
  any NCBI translation table (`--translation-table`, default 11; tables
  generated from and checked against Biopython), an optional first-codon Met
  rule (`--no-first-codon-as-met` turns it off), and Sequence Ontology effect
  terms (synonymous_variant, missense_variant, stop_gained, stop_lost,
  start_lost, stop/start_retained_variant, frameshift_variant,
  inframe_insertion/deletion) with codon and protein change. Checked against
  an independent Biopython-based classification (28,766 SNPs, table 11;
  5,706 differences each for tables 11, 4, 1 and 2, with and without the
  first-codon rule): no differences.
- Synonymous / nonsynonymous distances: `--mode nonsyn-snps` and
  `--mode custom --weights "key=w,..."`, a per-locus weighted sum over
  allele, snps, indel_events, indel_bases, syn, nonsyn and frame_disrupted
  (the built-in modes are particular weightings), with `--translation-table`,
  `--no-first-codon-as-met` and `--coding-stats`. The counts are optional
  cache fields tied to a `genetic_code` metadata entry, and remain readable
  by older cgdist. Checked on both test datasets:
  - existing modes identical to 0.1.4;
  - `custom` weightings identical to the built-in modes;
  - syn + nonsyn + frame_disrupted = snps on all 739,554 cache entries;
  - counts equal to cgdist-diff (Biopython-validated) on 600 pairs;
  - completing a 0.1.4 cache reproduces every cached SNP/InDel count.
- `cgdist-diff` also aligns the two proteins, with the same settings as the
  aa-* modes, and with `--protein-diffs` lists each amino-acid
  substitution and InDel. Its counts match the protein cache on 300 random
  pairs.
- Protein-level distances (new module `core::protein_distance`):
  - modes `aa-hamming`, `aa-substitutions`, `aa-substitutions-indel-events`
    and `aa-substitutions-indel-residues`, plus the `aa_*` keys for
    `--weights`;
  - alleles translated on the fly; proteins deduplicated by hash (a CRC32
    collision between distinct proteins of a locus stops the run);
  - distinct proteins aligned with parasail, with `--aa-matrix` (any
    BLOSUM/PAM or a file), `--aa-gap-open` and `--aa-gap-extend`;
  - a separate `--protein-cache-file` whose settings must match the run.

  Checked:
  - aa-hamming equals Biopython-translated proteins on 990 sample pairs;
  - the scan kernel equals parasail's striped kernel on 229,422 real
    protein pairs (BLOSUM62, PAM250);
  - custom weightings equal the built-in protein modes.
- `--verify-alignments <fraction>`: re-check a deterministic fraction of new
  alignments against parasail's original kernel and stop with an error on
  any difference (`1` = every pair).

## [0.1.4] — 2026-09-29

Bug-fix release. Two cache bugs, both present since 0.1.0, are fixed.
**Upgrading is recommended for anyone who reuses cache files.** Without a cache,
and with caches never used in Hamming mode, distance matrices are unchanged:
all four modes were checked byte-identical against 0.1.3 on real
*L. monocytogenes* (300 samples × 1,748 loci; also with `--min-loci`, all
output formats and `--emit-pairs`/`--report-ci`) and *S. enterica* (120
samples × 8,558 loci) datasets, and cached runs matched uncached runs in every
mode tested, including after a Hamming run on the same cache.

### ⚠️ Advisory — caches shared with Hamming mode

In cgdist ≤ 0.1.3, `--mode hamming` together with `--cache-file` stored a
placeholder result (1 SNP, 0 InDels), without aligning, for every allele pair
of that run not already in the cache. A later `snps`, `snps-indel-events` or
`snps-indel-bases` run on the **same cache file** read those placeholders as
real alignments:

- on a cache first written in Hamming mode, the SNP/InDel distances silently
  equalled the Hamming distances for the allele pairs of the Hamming run;
- on a cache that already held alignments, pairs involving alleles absent from
  the schema FASTA got distance 1 instead of 0 (with default settings).

You are affected only if one cache file was used both with `--mode hamming`
and with an SNP/InDel mode. Delete such caches, or rebuild them with
`--force-recompute`; SNP/InDel matrices computed from them should be
regenerated. cgdist now prints a warning when it loads a cache last written in
Hamming mode by a version ≤ 0.1.3; it cannot detect a cache that a later
SNP/InDel run saved again, so check your workflows rather than relying on the
warning alone. Runs using the `hamming` hasher, separate
caches per mode, or no cache were never affected. (The Hamming distances in the
cgDist article were computed with an external tool, cgmlst-dists, so the
article's results are not affected.)

### Fixed

- `--mode hamming` no longer reads or writes alignment statistics in the cache:
  Hamming distances are answered directly (different alleles = 1), and a
  Hamming run no longer adds entries to an existing cache file.
- Sequence lengths are no longer lost from enriched caches. With
  `--enrich-lengths`, saving a cache that gained new pairs wrote every entry
  without lengths, so the documented recombination workflow (fresh cache +
  `--enrich-lengths`) produced a cache with no lengths and
  `recombination_candidate_analyzer` silently flagged no loci. Newly aligned
  pairs now record their length immediately, and saving keeps known lengths.
- `recombination_candidate_analyzer` now refuses a cache without sequence
  lengths (instead of reporting zero candidates) and warns when only some
  entries have them.
- Allele pairs that cannot be aligned (an allele in the profiles has no
  sequence in the schema FASTA, e.g. profiles called with a newer or different
  schema) are now reported on every run with a warning giving their number,
  the affected loci and the effect: such pairs count as 0 (1 in snps mode with
  `--hamming-fallback`), so distances involving them are underestimated.
  Previously this happened silently and the log counted these pairs as
  computed alignments. Distances are unchanged.
- Cache compatibility compares the numeric alignment parameters only; a cache
  built with a preset and one built with identical custom parameters are now
  interchangeable.

### Changed

- Because new pairs carry lengths from the first run, `--report` and
  `--recomb-threshold` can use recombination signals without a separate
  enrichment pass.

## [0.1.3] — 2026-09-28

First release after publication of the cgDist article in *NAR Genomics and
Bioinformatics* (doi: [10.1093/nargab/lqag090](https://doi.org/10.1093/nargab/lqag090)).
It adds opt-in per-pair data-quality reporting and an HTML analyst dashboard.
**The distance algorithms and the default distance-matrix output are unchanged
from 0.1.2**: all four distance modes (with and without `--min-loci`) and the
TSV / CSV / PHYLIP / NEXUS outputs were checked byte-identical against 0.1.2
on a real *L. monocytogenes* dataset (300 samples × 1,748 loci), including runs
with the new flags enabled.

### Added

- `--emit-pairs <file>`: long-format per-pair table with data-quality columns
  (`sample_i, sample_j, distance, shared_loci, total_loci, missing_frac`),
  surfacing how many shared loci each distance actually rests on.
- `--report-ci` (with `--emit-pairs`, level via `--ci-level`, default 0.95):
  a missingness confidence interval per pair (`dist_norm, ci_low, ci_high,
  ci_reliable`) accounting for the loci a pair does not share. Uses a
  deterministic union of an exact Beta-Binomial count interval and a normal
  compound interval (self-contained ln_gamma / normal_ppf, no new dependency),
  validated against a Monte-Carlo posterior predictive. `ci_reliable` is `false`
  in the low-information regime (few differing loci + high missingness), where
  coverage is information-limited regardless of method.
- `--report <file.html>`: a self-contained HTML analyst dashboard (no external
  assets, works offline) with dataset summary, per-sample data quality, the
  pairwise-distance distribution, a clusters-vs-threshold curve, and an
  interactive single-linkage outbreak-clustering explorer that shows the
  missingness CI of each edge. Nord color theme.
- Dashboard recombination view (per-pair recombinant-loci load, top
  recombinant pairs, per-sample load), shown when an enriched cache provides
  sequence lengths (`--recomb-threshold`, percent, default 3).
- Library API: `calculate_pairs_table`, `calculate_pairs_recombination`,
  `calculate_sample_distance_detailed`, `PairRow`, `write_pairs_long`,
  `write_html_report`.

All of the above are opt-in and do not change the default distance-matrix output.

### Changed

- cgDist is now described as a distance calculator for core **and whole**
  genome MLST (cg/wgMLST) in the crate metadata, README and API docs.
- Citation: README and `CITATION.cff` (`preferred-citation`) now cite the
  published NAR Genomics and Bioinformatics article (8(3):lqag090) with the
  full author list; the bioRxiv preprint is still linked.

### Fixed

- `--ci-level` outside (0, 1) (e.g. `95` instead of `0.95`) and
  `--recomb-threshold` outside (0, 100] are now rejected with a clear error
  instead of silently producing infinite / meaningless intervals.
- `--report-ci` without `--emit-pairs` is now an error (it previously had no
  effect).

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

[Unreleased]: https://github.com/genpat-it/cgDist/compare/v0.1.4...HEAD
[0.1.4]: https://github.com/genpat-it/cgDist/compare/v0.1.3...v0.1.4
[0.1.3]: https://github.com/genpat-it/cgDist/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/genpat-it/cgDist/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/genpat-it/cgDist/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/genpat-it/cgDist/releases/tag/v0.1.0
