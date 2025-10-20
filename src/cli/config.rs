// config.rs - Configuration file support

use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Config {
    // Input/Output
    pub profiles: Option<String>,
    pub schema: Option<String>,
    pub output: Option<String>,

    // Core settings
    pub hasher_type: Option<String>,
    pub mode: Option<String>,
    pub format: Option<String>,
    pub missing_char: Option<String>,

    // Performance
    pub threads: Option<usize>,
    pub cache_file: Option<String>,
    pub cache_note: Option<String>,

    // Quality filters
    pub sample_threshold: Option<f64>,
    pub locus_threshold: Option<f64>,
    pub min_loci: Option<usize>,

    // Sample/Loci filtering
    pub include_samples: Option<String>,
    pub exclude_samples: Option<String>,
    pub include_loci: Option<String>,
    pub exclude_loci: Option<String>,
    pub include_loci_list: Option<String>,
    pub exclude_loci_list: Option<String>,
    pub include_samples_list: Option<String>,
    pub exclude_samples_list: Option<String>,

    // Alignment settings (ignored for hamming)
    pub alignment_mode: Option<String>,
    pub match_score: Option<i32>,
    pub mismatch_penalty: Option<i32>,
    pub gap_open: Option<i32>,
    pub gap_extend: Option<i32>,

    // Flags
    pub no_hamming_fallback: Option<bool>,
    pub force_recompute: Option<bool>,
    pub dry_run: Option<bool>,
    pub save_alignments: Option<String>,
}

impl Config {
    /// Create a new empty configuration
    pub fn new() -> Self {
        Self {
            profiles: None,
            schema: None,
            output: None,
            hasher_type: None,
            mode: None,
            format: None,
            missing_char: None,
            threads: None,
            cache_file: None,
            cache_note: None,
            sample_threshold: None,
            locus_threshold: None,
            min_loci: None,
            include_samples: None,
            exclude_samples: None,
            include_loci: None,
            exclude_loci: None,
            include_loci_list: None,
            exclude_loci_list: None,
            include_samples_list: None,
            exclude_samples_list: None,
            alignment_mode: None,
            match_score: None,
            mismatch_penalty: None,
            gap_open: None,
            gap_extend: None,
            no_hamming_fallback: None,
            force_recompute: None,
            dry_run: None,
            save_alignments: None,
        }
    }

    /// Load configuration from TOML file
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, String> {
        let path = path.as_ref();
        let content = fs::read_to_string(path)
            .map_err(|e| format!("Failed to read config file '{}': {}", path.display(), e))?;

        let config: Config = toml::from_str(&content)
            .map_err(|e| format!("Failed to parse config file '{}': {}", path.display(), e))?;

        println!("📄 Loaded configuration from: {}", path.display());
        Ok(config)
    }

    /// Save configuration to TOML file
    pub fn to_file<P: AsRef<Path>>(&self, path: P) -> Result<(), String> {
        let path = path.as_ref();
        let content = toml::to_string_pretty(self)
            .map_err(|e| format!("Failed to serialize config: {}", e))?;

        fs::write(path, content)
            .map_err(|e| format!("Failed to write config file '{}': {}", path.display(), e))?;

        println!("📄 Saved configuration to: {}", path.display());
        Ok(())
    }

    /// Generate a sample configuration file with comments
    pub fn generate_sample() -> String {
        r#"# cgdist.toml - Configuration file for cgdist
# Command line arguments will override these settings

# =============================================================================
# INPUT/OUTPUT
# =============================================================================

# Path to allelic profile matrix (.tsv or .csv)
profiles = "/path/to/profiles.tsv"

# Path to FASTA schema directory or schema file (not needed for hamming hasher)
schema = "/path/to/schema"

# Output distance matrix file
output = "distances.tsv"

# =============================================================================
# CORE SETTINGS
# =============================================================================

# Allele hasher type: crc32, sha256, md5, sequence, hamming
hasher_type = "crc32"

# Distance mode: snps, snps-indel-events, snps-indel-bases, hamming
# Note: ignored when hasher_type = "hamming"
mode = "snps"

# Output format: tsv, csv, phylip, nexus
format = "tsv"

# Missing data character
missing_char = "-"

# =============================================================================
# PERFORMANCE
# =============================================================================

# Number of threads (omit for auto-detection)
threads = 128

# Cache file path for ultra-fast reuse (.lz4 extension)
# Note: not used with hamming hasher
cache_file = "cache.lz4"

# User note to save with the cache for future reference
cache_note = "My analysis run"

# =============================================================================
# QUALITY FILTERS
# =============================================================================

# Sample quality filter: minimum fraction of non-missing loci per sample (0.0-1.0)
sample_threshold = 0.8

# Locus quality filter: minimum fraction of non-missing samples per locus (0.0-1.0)
locus_threshold = 0.95

# Minimum number of shared loci required for distance calculation
min_loci = 100

# =============================================================================
# SAMPLE/LOCI FILTERING
# =============================================================================

# Include only samples matching regex pattern
# include_samples = "pattern.*"

# Exclude samples matching regex pattern
# exclude_samples = "control.*"

# Include only loci matching regex pattern
# include_loci = "INNUENDO.*"

# Exclude loci matching regex pattern
# exclude_loci = "deprecated.*"

# Include only loci listed in a file (one locus per line)
include_loci_list = "efsa_loci.tsv"

# Exclude loci listed in a file (one locus per line)
# exclude_loci_list = "blacklist.txt"

# Include only samples listed in a file (one sample per line)
# include_samples_list = "samples.txt"

# Exclude samples listed in a file (one sample per line)
# exclude_samples_list = "exclude.txt"

# =============================================================================
# ALIGNMENT SETTINGS (ignored for hamming hasher)
# =============================================================================

# Alignment mode: dna, dna-strict, dna-permissive, custom
alignment_mode = "dna"

# Custom alignment scores (overrides preset mode)
# match_score = 2
# mismatch_penalty = -1
# gap_open = 5
# gap_extend = 2

# =============================================================================
# FLAGS
# =============================================================================

# Disable Hamming fallback for SNPs-only mode
no_hamming_fallback = false

# Force recomputation ignoring cache compatibility
force_recompute = false

# Validate inputs without computation (dry run)
dry_run = false

# Save detailed alignments to file (TSV format)
# save_alignments = "alignments.tsv"
"#
        .to_string()
    }
}

impl Default for Config {
    fn default() -> Self {
        Self::new()
    }
}
