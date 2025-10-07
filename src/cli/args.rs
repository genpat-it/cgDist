// args.rs - Command line arguments definition

use argh::FromArgs;

#[derive(FromArgs)]
/// cgDist - High-performance distance matrix calculator
pub struct Args {
    /// path to FASTA schema directory or schema file
    #[argh(option)]
    pub schema: Option<String>,
    
    /// path to allelic profile matrix (.tsv or .csv)
    #[argh(option)]
    pub profiles: Option<String>,
    
    /// output distance matrix file
    #[argh(option)]
    pub output: Option<String>,
    
    /// distance mode: snps, snps-indel-events, snps-indel-bases, hamming (default: snps)
    #[argh(option, default = "String::from(\"snps\")")]
    pub mode: String,
    
    /// output format: tsv, csv, phylip, nexus (default: tsv)
    #[argh(option, default = "String::from(\"tsv\")")]
    pub format: String,
    
    /// missing data character (default: -)
    #[argh(option, default = "String::from(\"-\")")]
    pub missing_char: String,
    
    /// minimum number of shared loci required for distance calculation (default: 0)
    #[argh(option, default = "0")]
    pub min_loci: usize,
    
    /// number of threads (default: auto-detect)
    #[argh(option)]
    pub threads: Option<usize>,
    
    /// sample quality filter: minimum fraction of non-missing loci per sample (0.0-1.0, default: 0.0 = no filter)
    #[argh(option, default = "0.0")]
    pub sample_threshold: f64,
    
    /// locus quality filter: minimum fraction of non-missing samples per locus (0.0-1.0, default: 0.0 = no filter)
    #[argh(option, default = "0.0")]
    pub locus_threshold: f64,
    
    /// include only samples matching regex pattern
    #[argh(option)]
    pub include_samples: Option<String>,
    
    /// exclude samples matching regex pattern
    #[argh(option)]
    pub exclude_samples: Option<String>,
    
    /// include only loci matching regex pattern
    #[argh(option)]
    pub include_loci: Option<String>,
    
    /// exclude loci matching regex pattern
    #[argh(option)]
    pub exclude_loci: Option<String>,
    
    /// include only loci listed in a file (one locus per line)
    #[argh(option)]
    pub include_loci_list: Option<String>,
    
    /// exclude loci listed in a file (one locus per line)
    #[argh(option)]
    pub exclude_loci_list: Option<String>,
    
    /// include only samples listed in a file (one sample per line)
    #[argh(option)]
    pub include_samples_list: Option<String>,
    
    /// exclude samples listed in a file (one sample per line)
    #[argh(option)]
    pub exclude_samples_list: Option<String>,
    
    /// disable Hamming fallback for SNPs-only mode (default: enabled for SNPs mode, N/A for indel modes)
    #[argh(switch)]
    pub no_hamming_fallback: bool,
    
    /// cache file path for ultra-fast reuse (.lz4 extension)
    #[argh(option)]
    pub cache_file: Option<String>,
    
    /// user note to save with the cache for future reference
    #[argh(option)]
    pub cache_note: Option<String>,
    
    /// enrich cache with nucleotide sequence lengths from schema
    #[argh(switch)]
    pub enrich_lengths: bool,
    
    /// output file for enriched cache (default: overwrites input cache)
    #[argh(option)]
    pub enrich_output: Option<String>,
    
    /// save detailed alignments to file (TSV format)
    #[argh(option)]
    pub save_alignments: Option<String>,
    
    /// alignment mode: dna, dna-strict, dna-permissive, custom (default: dna)
    #[argh(option, default = "String::from(\"dna\")")]
    pub alignment_mode: String,
    
    /// custom match score (overrides preset mode, enables custom mode)
    #[argh(option)]
    pub match_score: Option<i32>,
    
    /// custom mismatch penalty (overrides preset mode, enables custom mode)
    #[argh(option)]
    pub mismatch_penalty: Option<i32>,
    
    /// custom gap open penalty (overrides preset mode, enables custom mode)
    #[argh(option)]
    pub gap_open: Option<i32>,
    
    /// custom gap extend penalty (overrides preset mode, enables custom mode)
    #[argh(option)]
    pub gap_extend: Option<i32>,
    
    /// force recomputation ignoring cache compatibility (start fresh)
    #[argh(switch)]
    pub force_recompute: bool,
    
    /// build cache only without computing distance matrix
    #[argh(switch)]
    pub cache_only: bool,
    
    /// detect potential recombination events: output log of allele pairs exceeding threshold
    #[argh(option)]
    pub recombination_log: Option<String>,
    
    /// threshold for recombination detection (SNPs + InDel bases, default: 20)
    #[argh(option, default = "20")]
    pub recombination_threshold: usize,
    
    
    /// show matrix statistics and diversity metrics only, then exit
    #[argh(switch)]
    pub stats_only: bool,
    
    /// benchmark mode: measure alignment processing speed (pairs/second) and exit
    #[argh(switch)]
    pub benchmark: bool,
    
    /// benchmark duration in seconds (default: 15)
    #[argh(option, default = "15")]
    pub benchmark_duration: u64,
    
    /// validate inputs without computation (dry run)
    #[argh(switch)]
    pub dry_run: bool,
    
    /// allele hasher type: crc32, sha256, md5, sequence, hamming (default: crc32)
    #[argh(option, default = "String::from(\"crc32\")")]
    pub hasher_type: String,
    
    /// inspect cache file instead of running distance calculation
    #[argh(option)]
    pub inspector: Option<String>,
    
    /// path to TOML configuration file
    #[argh(option)]
    pub config: Option<String>,
    
    /// generate sample configuration file and exit
    #[argh(switch)]
    pub generate_config: bool,
}