# cgDist 🧬

[![Rust](https://img.shields.io/badge/rust-%23000000.svg?style=for-the-badge&logo=rust&logoColor=white)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/genpat-it/cgDist/actions/workflows/ci-and-docker.yml/badge.svg)](https://github.com/genpat-it/cgDist/actions/workflows/ci-and-docker.yml)

**Ultra-fast SNP/indel-level distance calculator for core genome MLST analysis**

cgDist is a high-performance Rust implementation for calculating genetic distances in bacterial genomics, specifically designed for epidemiological outbreak investigations and phylogenetic analysis.

## 🚀 Features

- **⚡ Ultra-fast**: Parallel processing with optimized algorithms
- **🎯 Precision**: SNP/indel-level distance calculation
- **🔧 Flexible**: Multiple hashing algorithms (CRC32, MD5, SHA256)
- **📊 Comprehensive**: Built-in comparison tools and statistical analysis
- **🏥 Epidemiological**: Star clustering analysis for outbreak investigation
- **🧬 Recombination Analysis**: Advanced tools for detecting genetic recombination events
- **💾 Efficient**: LZ4 compression for fast caching
- **📈 Scalable**: Memory-efficient processing of large datasets

## 📚 Table of Contents

- [Features](#-features)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Usage](#-usage)
- [Recombination Analysis](#-recombination-analysis)
- [Cache Inspector](#-cache-inspector)
- [Custom Hashers Plugin System](#-custom-hashers-plugin-system)
- [API Documentation](#-api-documentation)
- [Citation](#-citation)
- [Related Repositories](#-related-repositories)
- [Support](#-support)
- [License](#-license)

## 🔧 Installation

### Prerequisites

- Rust 1.70+ (install from [rustup.rs](https://rustup.rs/))
- Python 3.8+ (for analysis scripts)

### From Source

```bash
# Clone the repository
git clone https://github.com/genpat-it/cgDist.git
cd cgDist

# Build release version
cargo build --release

# The binary will be available at ./target/release/cgdist
```

### Using Cargo

```bash
cargo install cgdist
```

### Docker

```bash
# Build Docker image
docker build -t cgdist-rs .

# Run with Docker
docker run -v $(pwd):/data cgdist-rs cgdist --help
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

### Configuration File

Create a `cgdist-config.toml` file:

```toml
[general]
hasher_type = "crc32"
threads = 8
cache_file = "./cache/cgdist_cache.lz4"
missing_char = "-"

[filtering]
min_loci = 0
sample_threshold = 0.0
locus_threshold = 0.0

[output]
format = "tsv"
mode = "snps"
```

```bash
# Use configuration file
cgdist --schema schema_dir/ --profiles profiles.tsv --config cgdist-config.toml --output distances.tsv
```

## 📊 Usage

### Command Line Options

```bash
cgdist [OPTIONS]

MAIN OPTIONS:
    --schema <PATH>            Path to FASTA schema directory or schema file
    --profiles <PATH>          Path to allelic profile matrix (.tsv or .csv)
    --output <FILE>            Output distance matrix file
    --mode <MODE>              Distance mode [default: snps]
                               Options: snps, snps-indel-events, snps-indel-bases, hamming
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
    --threads <N>              Number of threads [default: auto-detect]
    --cache-file <FILE>        Cache file path (.lz4 extension)
    --cache-note <TEXT>        Note to save with cache
    --cache-only               Build cache only without computing distance matrix
    --force-recompute          Force recomputation ignoring cache
    --hasher-type <TYPE>       Allele hasher type [default: crc32]
                               Options: crc32, sha256, md5, sequence, hamming

CACHE ENRICHMENT OPTIONS:
    --enrich-lengths           Enrich cache with nucleotide sequence lengths from schema
    --enrich-output <FILE>     Output file for enriched cache [default: overwrites input cache]

RECOMBINATION DETECTION OPTIONS:
    --recombination-log <FILE>      Output log of allele pairs exceeding threshold
    --recombination-threshold <N>   Threshold for recombination detection (SNPs + InDel bases) [default: 20]

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

## 🧬 Recombination Analysis

cgDist includes powerful tools for detecting and analyzing genetic recombination events in bacterial populations. This feature is particularly useful for epidemiological studies and understanding horizontal gene transfer.

### Features

- **Mutation Density Analysis**: Detects loci with high SNP/indel density indicating recombination
- **Hamming Distance Filtering**: Focuses analysis on genetically related samples
- **Pairwise Recombination Summary**: Comprehensive overview of recombination between sample pairs
- **EFSA Loci Support**: Compatible with standardized loci sets for food safety applications
- **Distance Matrix Correction**: Adjusts distances by excluding recombinant loci

### Tool 1: Built-in Recombination Detection

The main `cgdist` binary can detect potential recombination events during distance calculation:

```bash
# Detect recombination events with default threshold (20 SNPs+InDel bases)
cgdist --schema schema_dir/ --profiles profiles.tsv --output distances.tsv \
    --recombination-log recombination_events.csv \
    --mode snps-indel-bases

# Custom threshold (e.g., 30 SNPs+InDel bases)
cgdist --schema schema_dir/ --profiles profiles.tsv --output distances.tsv \
    --recombination-log recombination_events.csv \
    --recombination-threshold 30 \
    --mode snps-indel-bases
```

**Output**: CSV log with locus, sample pairs, divergence percentages, and sequence lengths

### Tool 2: Recombination Analyzer (Post-processing)

For advanced analysis with Hamming filtering and EFSA loci support:

```bash
# Build recombination analyzer
cargo build --release --bin recombination_analyzer

# Step 1: Create enriched cache with sequence lengths
cgdist --schema schema_dir/ --profiles profiles.tsv --output distances.tsv \
    --cache-file cache.bin --enrich-lengths --mode snps-indel-bases

# Step 2: Run recombination analyzer
./target/release/recombination_analyzer \
    --enriched-cache cache.bin \
    --profiles profiles.tsv \
    --distance-matrix distances.tsv \
    --output-matrix corrected_distances.tsv \
    --recombination-log recombination_events.tsv \
    --threshold 3.0

# Custom threshold (5% mutation density)
./target/release/recombination_analyzer \
    --enriched-cache cache.bin \
    --profiles profiles.tsv \
    --distance-matrix distances.tsv \
    --output-matrix corrected_distances.tsv \
    --recombination-log recombination_events.tsv \
    --threshold 5.0
```

### Input Requirements

**For Tool 1 (Built-in Detection)**:
1. **Schema**: FASTA directory with allele sequences
2. **Profiles**: TSV/CSV file with sample-locus-allele matrix

**For Tool 2 (Recombination Analyzer)**:
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
1. **Corrected Distance Matrix**: Distance matrix with recombinant loci excluded
2. **Recombination Events Log**: Detailed list of detected recombination events with:
   - Sample pairs
   - Locus information
   - Mutation statistics (SNPs, InDels)
   - Density percentages
   - Sequence lengths

### Parameters

**Tool 1 (Built-in)**:
- `--recombination-threshold`: SNPs + InDel bases threshold (default: 20)

**Tool 2 (Recombination Analyzer)**:
- `--threshold`: Mutation density percentage (default: 3.0%)

### Complete Workflow Example

```bash
# Option A: Quick detection during distance calculation
cgdist --schema schema/ --profiles samples.tsv --output distances.tsv \
    --recombination-log events.csv --recombination-threshold 20 \
    --mode snps-indel-bases

# Option B: Advanced analysis with corrected distances
# Step 1: Create enriched cache
cgdist --schema schema/ --profiles samples.tsv --output distances.tsv \
    --cache-file cache.bin --enrich-lengths --mode snps-indel-bases

# Step 2: Analyze and correct
./target/release/recombination_analyzer \
    --enriched-cache cache.bin \
    --profiles samples.tsv \
    --distance-matrix distances.tsv \
    --output-matrix corrected_distances.tsv \
    --recombination-log events.tsv \
    --threshold 3.0
```

### Interpretation Guidelines

- **High SNP Density**: > 3% suggests recombination vs point mutations
- **High Indel Events**: Indicates potential mobile genetic elements
- **Pairwise Patterns**: Multiple loci with recombination between same samples suggests related strains
- **Hamming Filtering**: Ensures focus on epidemiologically relevant comparisons

### Performance Considerations

- **Memory Usage**: ~4-8GB for typical bacterial datasets (1000+ samples)
- **Processing Time**: 2-5 minutes for 21M cache entries on modern hardware
- **Scalability**: Linear with cache size, efficient for large epidemiological studies

### Scientific Applications

1. **Outbreak Investigation**: Identify recombination hotspots in transmission chains
2. **Evolutionary Analysis**: Track horizontal gene transfer patterns
3. **Food Safety**: Monitor recombination in foodborne pathogens
4. **Antimicrobial Resistance**: Detect resistance gene transfer events
5. **Population Genomics**: Understand bacterial population structure

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

## 🔗 Related Repositories

- **📊 Comprehensive Study**: [cgDist-study](https://github.com/genpat-it/cgDist-study) - Complete epidemiological analysis and performance benchmarks

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/genpat-it/cgDist/issues)
- **Discussions**: [GitHub Discussions](https://github.com/genpat-it/cgDist/discussions)
- **Email**: a.deruvo@izs.it

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**Made with ❤️ for the bioinformatics community**