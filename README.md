# cgDist 🧬

[![Crates.io](https://img.shields.io/crates/v/cgdist.svg?logo=rust)](https://crates.io/crates/cgdist)
[![Rust](https://img.shields.io/badge/rust-1.88%2B-orange.svg?logo=rust)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/genpat-it/cgDist/actions/workflows/ci-and-docker.yml/badge.svg)](https://github.com/genpat-it/cgDist/actions/workflows/ci-and-docker.yml)
[![GitHub release](https://img.shields.io/github/v/release/genpat-it/cgDist?label=release&color=blue)](https://github.com/genpat-it/cgDist/releases/latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18025926.svg)](https://doi.org/10.5281/zenodo.18025926)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2025.10.16.682749-bd2024)](https://doi.org/10.1101/2025.10.16.682749)

**Ultra-fast SNP/indel-level distance calculator for core genome MLST analysis**

cgDist is a high-performance Rust implementation for calculating genetic distances in bacterial genomics, specifically designed for epidemiological outbreak investigations and phylogenetic analysis.

## 🚀 Features

- **⚡ Ultra-fast**: Parallel processing with optimized algorithms
- **🎯 Precision**: SNP/indel-level distance calculation
- **🔧 Flexible**: Multiple hashing algorithms (CRC32, MD5, SHA256)
- **📊 Comprehensive**: Built-in comparison tools and statistical analysis
- **🧬 Recombination-candidate flagging**: Per-locus mutation-density screen to flag loci as recombination candidates for downstream phylogenetic confirmation
- **💾 Efficient**: LZ4 compression for fast caching
- **📈 Scalable**: Memory-efficient processing of large datasets

## 📚 Table of Contents

- [Features](#-features)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Usage](#-usage)
- [Recombination-Candidate Flagging](#-recombination-candidate-flagging)
- [Cache Inspector](#-cache-inspector)
- [Custom Hashers Plugin System](#-custom-hashers-plugin-system)
- [API Documentation](#-api-documentation)
- [Citation](#-citation)
- [Support](#-support)
- [License](#-license)

## 🔧 Installation

### Prerequisites

- Rust **1.88 or later** (the minimum supported Rust version, MSRV, is also declared in `Cargo.toml`). The pinned dependency set requires a recent toolchain, so simply using the latest stable Rust is recommended. The easiest way to install or update Rust is via [rustup.rs](https://rustup.rs/):
  ```bash
  # Install rustup (skip if already installed)
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

  # If rustup is already installed but Rust is older than 1.88, update with:
  rustup update stable
  ```
- Python 3.8+ (only for the validation scripts in `validation_test/`)
- **System build dependencies for `parasail-rs`**: the alignment backend
  is built from C source via CMake, which requires a C compiler and
  `zlib` development headers. Install once per machine:
  ```bash
  # Debian / Ubuntu / WSL
  sudo apt install build-essential cmake zlib1g-dev

  # RHEL / AlmaLinux / Rocky / CentOS / Fedora
  sudo dnf install gcc gcc-c++ cmake zlib-devel

  # macOS (Homebrew; Xcode Command Line Tools provide compiler + zlib)
  xcode-select --install
  brew install cmake
  ```
  Windows users are encouraged to use the [Docker
  image](#-docker-install-pre-built-on-ghcr) or
  [WSL2](https://learn.microsoft.com/windows/wsl/install), which provide
  a ready-to-build Linux environment. Native Windows builds additionally
  require zlib via [vcpkg](https://vcpkg.io/) (`vcpkg install zlib`) or
  MSYS2.

### From Source

```bash
# Clone the repository
git clone https://github.com/genpat-it/cgDist.git
cd cgDist

# Build release version
cargo build --release

# The binary will be available at ./target/release/cgdist
```

### Install from crates.io (recommended)

```bash
cargo install cgdist
```

This fetches the latest published release from
[crates.io](https://crates.io/crates/cgdist), builds it locally with
your stable Rust toolchain, and installs the `cgdist`, `inspector`,
and `recombination_candidate_analyzer` binaries to `~/.cargo/bin/`
(which should already be on your `PATH` after a default rustup
install). The deprecated `recombination_analyzer` binary is also
installed and forwards every argument to
`recombination_candidate_analyzer` with a deprecation notice — existing
scripts continue to work.

To pin a specific published version:

```bash
cargo install cgdist --version 0.1.1
```

For a fully reproducible build that uses exactly the dependency
versions we tested against, add `--locked` (this reuses the published
`Cargo.lock` instead of re-resolving to the newest compatible
dependencies, which may otherwise require an even more recent Rust
toolchain):

```bash
cargo install cgdist --locked
```

### Install from GitHub (specific tag or unreleased commits)

To install directly from the GitHub repository — useful for installing
an unreleased commit or for fully self-contained reproducibility when
citing the manuscript:

```bash
# Specific release tag
cargo install --git https://github.com/genpat-it/cgDist --tag v0.1.2 cgdist

# Latest state on the default branch
cargo install --git https://github.com/genpat-it/cgDist cgdist
```

cgdist is a binary crate, so its `Cargo.lock` is committed to the
repository to guarantee reproducible builds — this is the convention
recommended in the
[official Cargo FAQ](https://doc.rust-lang.org/cargo/faq.html#why-have-cargolock-in-version-control)
for binary crates.

### Docker

A multi-arch (linux/amd64 + linux/arm64) image is published to GitHub
Container Registry on every release:

```bash
# Pull the public image (no authentication required)
docker pull ghcr.io/genpat-it/cgdist:0.1.2
# or pin to the minor / major series:
# docker pull ghcr.io/genpat-it/cgdist:0.1
# docker pull ghcr.io/genpat-it/cgdist:latest   # tracks master HEAD

# Run with the image (mount your working directory at /data).
# The image's ENTRYPOINT is `cgdist`, so flags are passed directly:
docker run --rm -v $(pwd):/data ghcr.io/genpat-it/cgdist:0.1.2 \
    --schema /data/schema_dir --profiles /data/profiles.tsv \
    --output /data/distances.tsv --mode snps-indel-bases
```

To build the image locally instead of pulling (useful for development):

```bash
docker build -t cgdist:dev .
docker run --rm -v $(pwd):/data cgdist:dev --help
```

## 🚀 Quick Start

### Basic Distance Calculation

```bash
# Calculate SNP distances from cgMLST profiles
cgdist --schema schema_dir/ --profiles profiles.tsv --output distances.tsv

# Use different distance mode
cgdist --schema schema_dir/ --profiles profiles.tsv --output distances.tsv --mode snps-indel-bases

# Use different hashing algorithm
cgdist --schema schema_dir/ --profiles profiles.tsv --output distances.tsv --hasher-type sha256

# Enable cache for faster recomputation
cgdist --schema schema_dir/ --profiles profiles.tsv --output distances.tsv --cache-file cache.lz4

# Specify number of threads
cgdist --schema schema_dir/ --profiles profiles.tsv --output distances.tsv --threads 16
```

### Validation / Smoke Test

A self-contained validation suite with a small embedded test dataset
(3 loci, 10 samples, ~3 KB) is provided in
[`validation_test/`](validation_test/). It verifies algorithmic
correctness across all four distance modes (Hamming, SNPs,
SNPs+InDel-events, SNPs+InDel-bases), checks the mathematical invariant
`cgDist ≥ Hamming`, and confirms Parasail alignment integration.

```bash
# After building cgDist (see Installation)
cd validation_test
../target/release/cgdist --schema schema_crc32 --profiles profiles/test_profiles_crc32.tsv --output results/crc32_hamming.tsv --mode hamming --hasher-type crc32
../target/release/cgdist --schema schema_crc32 --profiles profiles/test_profiles_crc32.tsv --output results/crc32_snps.tsv --mode snps --hasher-type crc32 --hamming-fallback
../target/release/cgdist --schema schema_crc32 --profiles profiles/test_profiles_crc32.tsv --output results/crc32_snps_indel_contiguous.tsv --mode snps-indel-contiguous --hasher-type crc32
../target/release/cgdist --schema schema_crc32 --profiles profiles/test_profiles_crc32.tsv --output results/crc32_snps_indel_bases.tsv --mode snps-indel-bases --hasher-type crc32

# Verify expected results
python3 run_validation.py
```

Expected output: `🎉 ALL VALIDATION TESTS PASSED!` See
[`validation_test/README.md`](validation_test/README.md) for details on
the test design, expected distances, and how to regenerate the fixture
from scratch.

The validation suite also runs automatically in CI on every push and
pull request (see `.github/workflows/ci-and-docker.yml`).

### Configuration File

A configuration file is **optional**: every parameter accepted by `cgdist` also has
a CLI flag. The configuration file simply lets you persist commonly-used
settings without retyping them. **CLI flags always override TOML values
when both are provided.**

A canonical example is shipped at
[`examples/cgdist-config.toml`](examples/cgdist-config.toml); a
Hamming-mode variant is at
[`examples/hamming-config.toml`](examples/hamming-config.toml). Both
files use the same flat key structure (no `[sections]`), and the same
key names as the corresponding CLI flags (the only normalization is
that CLI flag dashes become underscores in TOML — e.g. `--hasher-type`
becomes `hasher_type`).

You can also generate a fresh annotated sample with:

```bash
cgdist --generate-config > cgdist-config.toml
```

A minimal example (alignment-based mode):

```toml
profiles = "profiles.tsv"
schema = "schema/"
output = "distances.tsv"
hasher_type = "crc32"
mode = "snps"            # legacy alias snps-indel-events == snps-indel-contiguous (deprecated)
format = "tsv"
missing_char = "-"
threads = 1              # default; set to 0 for auto-detect
hamming_fallback = false # opt-in (see Hamming Fallback section below)
```

```bash
# Use a configuration file
cgdist --config cgdist-config.toml

# CLI overrides example: config says threads=1, but the CLI wins → 16 threads used
cgdist --config cgdist-config.toml --threads 16
```

#### CLI vs TOML precedence

When both the configuration file and the command line specify the same
parameter, the **command-line value wins**. Internally this is
implemented by loading the TOML first, then overlaying any CLI flag
that the user explicitly set. The same rule applies to switches: e.g.
if the TOML says `hamming_fallback = false` but you pass
`--hamming-fallback` on the command line, the fallback will be enabled
for that run.

## 📊 Usage

### Command Line Options

```bash
cgdist [OPTIONS]

MAIN OPTIONS:
    --schema <PATH>            Path to FASTA schema directory or schema file
    --profiles <PATH>          Path to allelic profile matrix (.tsv or .csv)
    --output <FILE>            Output distance matrix file
    --mode <MODE>              Distance mode [default: snps]
                               Options: snps, snps-indel-contiguous, snps-indel-bases, hamming
                               (legacy alias snps-indel-events == snps-indel-contiguous, deprecated)
    --format <FORMAT>          Output format [default: tsv]
                               Options: tsv, csv, phylip, nexus

FILTERING OPTIONS:
    --min-loci <N>             Minimum shared loci for distance calculation [default: 0]
    --sample-threshold <VAL>   Sample quality filter (0.0-1.0) [default: 0.0]
    --locus-threshold <VAL>    Locus quality filter (0.0-1.0) [default: 0.0]
    --include-samples <REGEX>  Include only samples matching regex pattern
    --exclude-samples <REGEX>  Exclude samples matching regex pattern
    --include-loci <REGEX>     Include only loci matching regex pattern
    --exclude-loci <REGEX>     Exclude loci matching regex pattern
    --include-samples-list <FILE>  Include samples from file (one per line)
    --exclude-samples-list <FILE>  Exclude samples from file (one per line)
    --include-loci-list <FILE>     Include loci from file (one per line)
    --exclude-loci-list <FILE>     Exclude loci from file (one per line)

ALIGNMENT OPTIONS:
    --alignment-mode <MODE>    Alignment mode [default: dna]
                               Options: dna, dna-strict, dna-permissive, custom
    --match-score <N>          Custom match score (enables custom mode)
    --mismatch-penalty <N>     Custom mismatch penalty (enables custom mode)
    --gap-open <N>             Custom gap open penalty (enables custom mode)
    --gap-extend <N>           Custom gap extend penalty (enables custom mode)
    --save-alignments <FILE>   Save detailed alignments to TSV file

PERFORMANCE OPTIONS:
    --threads <N>              Number of threads [default: 1; pass 0 for auto-detect]
    --cache-file <FILE>        Cache file path (.lz4 extension)
    --cache-note <TEXT>        Note to save with cache
    --cache-only               Build cache only without computing distance matrix
    --force-recompute          Force recomputation ignoring cache
    --hasher-type <TYPE>       Allele hasher type [default: crc32]
                               Options: crc32, sha256, md5, sequence, hamming

CACHE ENRICHMENT OPTIONS:
    --enrich-lengths           Enrich cache with nucleotide sequence lengths from schema
    --enrich-output <FILE>     Output file for enriched cache [default: overwrites input cache]

RECOMBINATION-CANDIDATE FLAGGING OPTIONS:
    --candidate-recombination-log <FILE>        Output flagging log (one row per flagged
                                                candidate locus)
    --candidate-recombination-threshold <N>     SNPs + InDel-bases threshold above which a
                                                locus is flagged as a recombination candidate
                                                [default: 20]
    (legacy aliases --recombination-log / --recombination-threshold are still accepted)

OTHER OPTIONS:
    --missing-char <CHAR>      Missing data character [default: -]
    --no-hamming-fallback      Disable Hamming fallback for SNPs mode
    --stats-only               Show matrix statistics only
    --benchmark                Measure alignment processing speed
    --benchmark-duration <N>   Benchmark duration in seconds [default: 15]
    --dry-run                  Validate inputs without computation
    --inspector <FILE>         Inspect cache file
    --config <FILE>            Path to TOML configuration file
    --generate-config          Generate sample configuration file
    --help                     Display usage information
```

### Supported Input Formats

**Schema** (FASTA directory):
- Individual FASTA files per locus
- Each file contains allele sequences
- File names correspond to locus names

**Profiles** (allelic profiles):
- **TSV**: Tab-separated values
- **CSV**: Comma-separated values
- Format: Sample name | Locus1 | Locus2 | ... | LocusN
- Missing data represented by configurable character (default: `-`)

**Cache files**:
- **LZ4**: Compressed cache files (.lz4 or .bin extension)
- Automatic compression/decompression

### Output Formats

- **TSV**: Tab-separated distance matrix (default)
- **CSV**: Comma-separated distance matrix
- **PHYLIP**: Phylogenetic analysis format
- **NEXUS**: Nexus format for phylogenetic tools

## 🧬 Recombination-Candidate Flagging

cgDist includes a companion screen that flags **candidate** recombinant loci based on per-locus mutation density. This is **not** a recombination detector: confirmation of recombination requires downstream phylogeny-aware tools (e.g. Gubbins, ClonalFrameML, fastGEAR). The flagging output identifies which loci warrant that follow-up.

### Features

- **Mutation Density Analysis**: Flags loci with high SNP/indel density per alignment as recombination candidates
- **Hamming Distance Filtering**: Focuses analysis on genetically related sample pairs
- **Pairwise Flagging Summary**: Per sample-pair count of flagged loci
- **EFSA Loci Support**: Compatible with standardized loci sets for food safety applications
- **Distance Matrix Correction**: Recomputes distances excluding flagged loci

### Tool 1: Built-in Flagging

The main `cgdist` binary can flag candidate recombinant loci during distance calculation:

```bash
# Flag candidate loci with default threshold (20 SNPs+InDel bases)
cgdist --schema schema_dir/ --profiles profiles.tsv --output distances.tsv \
    --candidate-recombination-log recombination_events.csv \
    --mode snps-indel-bases

# Custom threshold (e.g., 30 SNPs+InDel bases)
cgdist --schema schema_dir/ --profiles profiles.tsv --output distances.tsv \
    --candidate-recombination-log recombination_events.csv \
    --candidate-recombination-threshold 30 \
    --mode snps-indel-bases
```

> **Note** — `--recombination-log` and `--recombination-threshold` are kept as deprecated aliases of the canonical `--candidate-recombination-*` flags (a deprecation warning is printed when used). Existing scripts continue to work.

**Output**: CSV log with locus, sample pairs, divergence percentages, and sequence lengths

### Tool 2: Recombination-Candidate Analyzer (Post-processing)

For advanced flagging with Hamming filtering and EFSA loci support:

```bash
# Build the candidate analyzer
cargo build --release --bin recombination_candidate_analyzer

# Step 1: Create enriched cache with sequence lengths
cgdist --schema schema_dir/ --profiles profiles.tsv --output distances.tsv \
    --cache-file cache.bin --enrich-lengths --mode snps-indel-bases

# Step 2: Run the analyzer
./target/release/recombination_candidate_analyzer \
    --enriched-cache cache.bin \
    --profiles profiles.tsv \
    --distance-matrix distances.tsv \
    --output-matrix corrected_distances.tsv \
    --candidate-recombination-log recombination_events.tsv \
    --threshold 3.0

# Custom threshold (5% mutation density)
./target/release/recombination_candidate_analyzer \
    --enriched-cache cache.bin \
    --profiles profiles.tsv \
    --distance-matrix distances.tsv \
    --output-matrix corrected_distances.tsv \
    --candidate-recombination-log recombination_events.tsv \
    --threshold 5.0
```

> **Note** — the binary `recombination_analyzer` and the flag `--recombination-log` are kept as deprecated aliases for backward compatibility (a deprecation notice is printed when invoked). Existing scripts continue to work.

### Input Requirements

cgDist consumes standard cgMLST outputs. Profiles and schemas can be
generated, for example, with [ChewBBACA](https://github.com/B-UMMI/chewBBACA)
(Silva et al. 2018) or downloaded from [Chewie-NS](https://chewbbaca.online/)
(Mamede et al. 2020).

**For Tool 1 (Built-in Flagging)**:
1. **Schema**: FASTA directory with allele sequences (e.g. ChewBBACA schema directory)
2. **Profiles**: TSV/CSV file with sample-locus-allele matrix (e.g. ChewBBACA `results_alleles.tsv`)

**For Tool 2 (Recombination-Candidate Analyzer)**:
1. **Enriched Cache**: `.bin` file generated with `--enrich-lengths` option
2. **Allelic Profiles**: TSV file with sample-locus-allele matrix
3. **Distance Matrix**: Original distance matrix from cgdist
4. **EFSA Loci** (optional): TSV file listing loci of interest

### Output Files

**Tool 1 Output**: `recombination_events.csv`
- Locus name
- Sample pairs
- Divergence percentage
- Sequence lengths
- SNPs and InDel counts

**Tool 2 Outputs**:
1. **Corrected Distance Matrix**: Distance matrix with flagged candidate loci excluded
2. **Flagging Log**: Detailed list of flagged candidate loci with:
   - Sample pairs
   - Locus information
   - Mutation statistics (SNPs, InDels)
   - Density percentages
   - Sequence lengths

### Parameters

**Tool 1 (Built-in)**:
- `--candidate-recombination-threshold`: SNPs + InDel bases threshold (default: 20). Legacy alias `--recombination-threshold` is also accepted.
- `--candidate-recombination-log`: output flagging log path. Legacy alias `--recombination-log` is also accepted.

**Tool 2 (Recombination-Candidate Analyzer)**:
- `--threshold`: Mutation density percentage (default: 3.0%)
- `--candidate-recombination-log`: output flagging log path. Legacy alias `--recombination-log` is also accepted.

### Complete Workflow Example

```bash
# Option A: Quick flagging during distance calculation
cgdist --schema schema/ --profiles samples.tsv --output distances.tsv \
    --candidate-recombination-log events.csv --candidate-recombination-threshold 20 \
    --mode snps-indel-bases

# Option B: Advanced flagging with corrected distances
# Step 1: Create enriched cache
cgdist --schema schema/ --profiles samples.tsv --output distances.tsv \
    --cache-file cache.bin --enrich-lengths --mode snps-indel-bases

# Step 2: Flag and correct
./target/release/recombination_candidate_analyzer \
    --enriched-cache cache.bin \
    --profiles samples.tsv \
    --distance-matrix distances.tsv \
    --output-matrix corrected_distances.tsv \
    --candidate-recombination-log events.tsv \
    --threshold 3.0
```

### Interpretation Guidelines

- **High SNP Density**: > 3% flags a locus as a recombination candidate (confirm with phylogeny-aware tools)
- **High Indel Events**: May indicate mobile genetic elements; warrants downstream inspection
- **Pairwise Patterns**: Multiple flagged loci between the same sample pair suggests related strains
- **Hamming Filtering**: Ensures focus on epidemiologically relevant comparisons

### Performance Considerations

- **Memory Usage**: ~4-8GB for typical bacterial datasets (1000+ samples)
- **Processing Time**: 2-5 minutes for 21M cache entries on modern hardware
- **Scalability**: Linear with cache size, efficient for large epidemiological studies

### Scientific Applications

1. **Outbreak Investigation**: Flag candidate recombination loci in transmission chains for downstream confirmation
2. **Evolutionary Analysis**: Identify candidate horizontal gene transfer events
3. **Food Safety**: Screen for recombination signatures in foodborne pathogens
4. **Antimicrobial Resistance**: Flag candidate resistance gene transfer events
5. **Population Genomics**: Identify loci that may bias clonal-frame distance estimates

## 🔍 Cache Inspector

The `inspector` tool provides detailed analysis of cgDist cache files, including validation, statistics, and compatibility checks.

### Building the Inspector

```bash
cargo build --release --bin inspector
```

### Basic Usage

```bash
# Show cache summary
./target/release/inspector --cache cache.lz4

# Detailed information including all loci
./target/release/inspector --cache cache.lz4 --detailed

# Show entries for specific locus
./target/release/inspector --cache cache.lz4 --show-locus locus_name

# Validate cache integrity
./target/release/inspector --cache cache.lz4 --validate

# Export cache summary to TSV
./target/release/inspector --cache cache.lz4 --export-summary summary.tsv

# Check top N loci by entry count
./target/release/inspector --cache cache.lz4 --top-loci 20
```

### Advanced Features

```bash
# Detect alignment mode from parameters
./target/release/inspector --cache cache.lz4 --detect-mode

# Check compatibility with specific alignment parameters
./target/release/inspector --cache cache.lz4 \
    --check-compatibility "5,-4,-10,-1"  # match,mismatch,gap_open,gap_extend

# Quiet mode for scripting
./target/release/inspector --cache cache.lz4 --validate --quiet
```

### Use Cases

1. **Cache Validation**: Verify cache file integrity before reuse
2. **Troubleshooting**: Diagnose cache compatibility issues
3. **Statistics**: Understand cache size and loci distribution
4. **Auditing**: Track which alignment parameters were used
5. **Quality Control**: Ensure cache matches expected schema

## 🔌 Custom Hashers Plugin System

cgDist provides a powerful plugin architecture for implementing custom hashing algorithms. This is particularly useful for specialized applications or compatibility with other tools.

### Implementing a Custom Hasher

Create a new hasher by implementing the `AlleleHasher` trait:

```rust
use cgdist::hashers::{AlleleHasher, AlleleHash};

/// Example: Simple nucleotide composition hasher
#[derive(Debug)]
pub struct CompositionHasher;

impl AlleleHasher for CompositionHasher {
    fn hash_sequence(&self, sequence: &str) -> AlleleHash {
        // Count nucleotides: A, T, G, C
        let mut counts = [0u8; 4]; // A, T, G, C
        for nucleotide in sequence.chars() {
            match nucleotide.to_ascii_uppercase() {
                'A' => counts[0] += 1,
                'T' => counts[1] += 1,
                'G' => counts[2] += 1,
                'C' => counts[3] += 1,
                _ => {} // Ignore ambiguous bases
            }
        }
        
        // Create hash from composition: AAAAATTTTGGGGCCCC format
        let hash_string = format!("A{}T{}G{}C{}", 
            counts[0], counts[1], counts[2], counts[3]);
        AlleleHash::String(hash_string)
    }
    
    fn parse_allele(&self, allele_str: &str, missing_char: &str) -> Result<AlleleHash, String> {
        if allele_str == missing_char {
            Ok(AlleleHash::Missing)
        } else {
            // Parse composition string or return as-is
            Ok(AlleleHash::String(allele_str.to_string()))
        }
    }
    
    fn name(&self) -> &'static str {
        "composition"
    }
    
    fn description(&self) -> &'static str {
        "Nucleotide composition-based hasher (A/T/G/C counts)"
    }
    
    fn validate_sequence(&self, sequence: &str) -> Result<(), String> {
        // Only allow ATGC nucleotides
        for ch in sequence.chars() {
            match ch.to_ascii_uppercase() {
                'A' | 'T' | 'G' | 'C' | 'N' => {}
                _ => return Err(format!("Invalid nucleotide: {}", ch)),
            }
        }
        Ok(())
    }
}
```

### Registering Your Custom Hasher

```rust
use cgdist::hashers::HasherRegistry;

fn main() {
    let mut registry = HasherRegistry::new();
    
    // Register your custom hasher
    registry.register_hasher("composition", Box::new(CompositionHasher));
    
    // Use it like any built-in hasher
    let hasher = registry.get_hasher("composition").unwrap();
    let hash = hasher.hash_sequence("ATCGATCG");
    
    println!("Hash: {}", hash); // Output: A2T2G2C2
}
```

### Advanced Custom Hasher Examples

#### 1. K-mer Based Hasher
```rust
#[derive(Debug)]
pub struct KmerHasher {
    k: usize,
}

impl KmerHasher {
    pub fn new(k: usize) -> Self {
        Self { k }
    }
}

impl AlleleHasher for KmerHasher {
    fn hash_sequence(&self, sequence: &str) -> AlleleHash {
        let mut kmers = Vec::new();
        let seq_bytes = sequence.as_bytes();
        
        if seq_bytes.len() >= self.k {
            for i in 0..=(seq_bytes.len() - self.k) {
                let kmer = std::str::from_utf8(&seq_bytes[i..i + self.k])
                    .unwrap_or("")
                    .to_string();
                kmers.push(kmer);
            }
        }
        
        kmers.sort();
        let hash_string = kmers.join("|");
        AlleleHash::String(hash_string)
    }
    
    // ... implement other required methods
}
```

#### 2. Custom Numeric Hasher
```rust
#[derive(Debug)]
pub struct CustomNumericHasher;

impl AlleleHasher for CustomNumericHasher {
    fn hash_sequence(&self, sequence: &str) -> AlleleHash {
        // Convert sequence to custom numeric representation
        let mut hash_value = 0u32;
        for (i, nucleotide) in sequence.chars().enumerate() {
            let base_value = match nucleotide.to_ascii_uppercase() {
                'A' => 0,
                'T' => 1,
                'G' => 2,
                'C' => 3,
                _ => 0, // Default for ambiguous
            };
            // Simple polynomial rolling hash
            hash_value = hash_value.wrapping_mul(4).wrapping_add(base_value);
        }
        AlleleHash::Crc32(hash_value)
    }
    
    fn parse_allele(&self, allele_str: &str, missing_char: &str) -> Result<AlleleHash, String> {
        if allele_str == missing_char {
            Ok(AlleleHash::Missing)
        } else {
            match allele_str.parse::<u32>() {
                Ok(value) => Ok(AlleleHash::Crc32(value)),
                Err(_) => Err(format!("Invalid numeric allele: {}", allele_str)),
            }
        }
    }
    
    fn name(&self) -> &'static str {
        "custom-numeric"
    }
    
    fn description(&self) -> &'static str {
        "Custom polynomial rolling hash for sequences"
    }
}
```

### Integration with cgdist CLI

To use custom hashers with the cgdist command-line tool, you can:

1. **Fork and modify**: Add your hasher to the registry in `src/main.rs`
2. **Configuration file**: Load hashers from a configuration file
3. **Dynamic loading**: Use Rust's plugin system (advanced)

Example integration in `main.rs`:
```rust
fn create_registry() -> HasherRegistry {
    let mut registry = HasherRegistry::new();
    
    // Add your custom hashers here
    registry.register_hasher("composition", Box::new(CompositionHasher));
    registry.register_hasher("kmer3", Box::new(KmerHasher::new(3)));
    registry.register_hasher("custom-numeric", Box::new(CustomNumericHasher));
    
    registry
}
```

### Use Cases for Custom Hashers

1. **Legacy Compatibility**: Match existing tool formats
2. **Domain-Specific**: Specialized algorithms for specific organisms
3. **Research**: Experimental hashing strategies
4. **Performance**: Optimized for specific hardware or datasets
5. **Compliance**: Meet specific regulatory or institutional requirements

### Best Practices

1. **Deterministic**: Ensure same sequence always produces same hash
2. **Collision-Resistant**: Minimize hash collisions for your use case
3. **Performance**: Consider computational overhead
4. **Validation**: Implement robust input validation
5. **Documentation**: Provide clear usage examples and limitations

The plugin architecture makes cgDist highly extensible while maintaining backward compatibility with existing workflows.

### Running the Custom Hasher Example

See the complete working example:

```bash
# Run the custom hasher demonstration
cargo run --example custom_hasher

# Output shows different hashers applied to test sequences:
# 🔌 cgDist Custom Hasher Examples
# ===================================
# 
# 📊 Available Hashers:
#   • crc32: Fast CRC32 checksum (chewBBACA compatible)
#   • composition: Nucleotide composition-based hasher (A/T/G/C counts)
#   • kmer3: K-mer composition hasher (sorted k-mers)
#   • polynomial: Polynomial rolling hash for sequences
# 
# 🧬 Testing hasher: composition
#    Description: Nucleotide composition-based hasher (A/T/G/C counts)
#    ATCGATCGATCG → A3T3G3C3
#    AAATTTGGGCCC → A3T3G3C3
#    ATGCATGCATGC → A3T3G3C3
```

This example demonstrates practical implementation patterns for:
- **Composition-based hashing**: Count nucleotide frequencies
- **K-mer analysis**: Extract and sort sequence k-mers  
- **Polynomial hashing**: Mathematical sequence encoding
- **Error handling**: Validation and missing data management



## 📖 API Documentation

### Rust API

```rust
use cgdist::{DistanceCalculator, Config};

// Create calculator with custom config
let config = Config::new()
    .hasher("crc32")
    .threads(8)
    .cache_enabled(true);

let calculator = DistanceCalculator::new(config);

// Calculate distances
let distances = calculator.calculate_from_file("sequences.fasta")?;
```

### Python Integration

```python
import subprocess
import pandas as pd

# Run cgdist from Python
result = subprocess.run([
    'cgdist',
    '--schema', 'schema_dir/',
    '--profiles', 'profiles.tsv',
    '--output', 'distances.tsv',
    '--mode', 'snps-indel-bases'
], capture_output=True, text=True)

# Check for errors
if result.returncode != 0:
    print(f"Error: {result.stderr}")
else:
    # Load results
    distances = pd.read_csv('distances.tsv', sep='\t', index_col=0)
    print(f"Distance matrix shape: {distances.shape}")
    print(distances.head())
```



## 📜 Citation

If you use cgDist in your research, please cite our preprint:

**de Ruvo, A.; Castelli, P.; Bucciacchio, A.; Mangone, I.; Mixao, V.; Borges, V.; Radomski, N.; Di Pasquale, A.** (2025). *cgDist: An Enhanced Algorithm for Efficient Calculation of pairwise SNP and InDel differences from Core Genome Multilocus Sequence Typing*. bioRxiv. DOI: [10.1101/2025.10.16.682749](https://doi.org/10.1101/2025.10.16.682749)

```bibtex
@article{deruvo2025cgdist,
  title = {cgDist: An Enhanced Algorithm for Efficient Calculation of pairwise SNP and InDel differences from Core Genome Multilocus Sequence Typing},
  author = {de Ruvo, Andrea and Castelli, Pierluigi and Bucciacchio, Andrea and Mangone, Iolanda and Mixao, Verónica and Borges, Vítor and Radomski, Nicolas and Di Pasquale, Adriano},
  year = {2025},
  month = {October},
  doi = {10.1101/2025.10.16.682749},
  journal = {bioRxiv},
  note = {Preprint. Software: https://github.com/genpat-it/cgDist}
}
```

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/genpat-it/cgDist/issues)
- **Discussions**: [GitHub Discussions](https://github.com/genpat-it/cgDist/discussions)
- **Email**: a.deruvo@izs.it

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**Made with ❤️ for the bioinformatics community**