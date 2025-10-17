# cgDist 🧬

[![Rust](https://img.shields.io/badge/rust-%23000000.svg?style=for-the-badge&logo=rust&logoColor=white)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)](https://github.com/genpat-it/cgDist)

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
# Calculate distances from FASTA alignment
cgdist sequences.fasta --output distances.tsv

# Use different hashing algorithm
cgdist sequences.fasta --hasher sha256 --output distances.tsv

# Enable parallel processing
cgdist sequences.fasta --threads 8 --output distances.tsv
```

### Configuration File

Create a `cgdist-config.toml` file:

```toml
[general]
hasher = "crc32"
threads = 8
cache_dir = "./cache"

[output]
format = "tsv"
include_headers = true
```

```bash
# Use configuration file
cgdist sequences.fasta --config cgdist-config.toml
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
    --force-recompute          Force recomputation ignoring cache
    --hasher-type <TYPE>       Allele hasher type [default: crc32]
                               Options: crc32, sha256, md5, sequence, hamming

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

- **FASTA**: Multi-sequence alignment files
- **CSV/TSV**: Tabular sequence data
- **Compressed**: LZ4 compressed files (automatic detection)

### Output Formats

- **TSV**: Tab-separated values (default)
- **CSV**: Comma-separated values
- **JSON**: Structured JSON output
- **Matrix**: Full distance matrix

## 🧬 Recombination Analysis

cgDist includes powerful tools for detecting and analyzing genetic recombination events in bacterial populations. This feature is particularly useful for epidemiological studies and understanding horizontal gene transfer.

### Features

- **Mutation Density Analysis**: Detects loci with high SNP/indel density indicating recombination
- **Hamming Distance Filtering**: Focuses analysis on genetically related samples
- **Pairwise Recombination Summary**: Comprehensive overview of recombination between sample pairs
- **EFSA Loci Support**: Compatible with standardized loci sets for food safety applications

### Basic Usage

```bash
# Build recombination analyzer
cargo build --release --bin cache_recombination_clean

# Run analysis with default settings (3% threshold, Hamming ≤ 15)
./target/release/cache_recombination_clean \
    enriched_cache.bin \
    allelic_profiles.tsv \
    efsa_loci.tsv

# Custom thresholds
./target/release/cache_recombination_clean \
    enriched_cache.bin \
    allelic_profiles.tsv \
    efsa_loci.tsv \
    5.0 \    # 5% mutation density threshold
    - \      # Missing data character
    10       # Hamming distance threshold ≤ 10
```

### Input Requirements

1. **Enriched Cache**: `.bin` file generated by cgdist with `--cache-file` option
2. **Allelic Profiles**: TSV file with sample-locus-allele matrix
3. **EFSA Loci**: TSV file listing loci of interest (one per line)

### Output Files

#### 1. Detailed Recombination Events (`recombination_analysis_with_samples.tsv`)

Contains individual recombination events with detailed metrics:

| Column | Description |
|--------|-------------|
| Sample1, Sample2 | Sample pair showing recombination |
| Locus | Specific locus with recombination |
| Allele1, Allele2 | CRC hashes of different alleles |
| SNPs | Number of single nucleotide polymorphisms |
| IndelEvents | Number of insertion/deletion events |
| TotalMutations | SNPs + IndelEvents |
| AvgLength | Average sequence length |
| SNPDensity% | SNP density per sequence length |
| IndelDensity% | Indel density per sequence length |
| TotalDensity% | Total mutation density |

#### 2. Pairwise Summary (`pairwise_recombination_summary.tsv`)

Summarizes recombination between sample pairs:

| Column | Description |
|--------|-------------|
| Sample1, Sample2 | Sample pair |
| RecombiningLoci | Number of loci showing recombination |
| TotalEFSALoci | Total loci analyzed (from EFSA list) |
| RecombinationPercentage% | Percentage of loci with recombination |

### Parameters

- **Mutation Density Threshold**: Minimum mutation density (%) to classify as recombination (default: 3.0%)
- **Missing Data Character**: Character representing missing alleles (default: "-")
- **Hamming Threshold**: Maximum Hamming distance between samples to include in analysis (default: 15)

### Example Analysis Workflow

```bash
# 1. Generate enriched cache with cgdist
./target/release/cgdist \
    --profiles samples.tsv \
    --schema schema_dir/ \
    --cache-file enriched_cache.bin \
    --mode snps-indel-events

# 2. Run recombination analysis
./target/release/cache_recombination_clean \
    enriched_cache.bin \
    samples.tsv \
    efsa_loci.tsv \
    3.0 - 15

# 3. Analyze results
head -20 recombination_analysis_with_samples.tsv
head -10 pairwise_recombination_summary.tsv
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
    'cgdist', 'sequences.fasta', '--output', 'distances.tsv'
], capture_output=True, text=True)

# Load results
distances = pd.read_csv('distances.tsv', sep='\t')
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