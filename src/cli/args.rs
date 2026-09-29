// args.rs - Command line arguments definition

use argh::FromArgs;

#[derive(FromArgs)]
/// cgDist - High-performance distance matrix calculator
pub struct Args {
    /// print the cgdist version and exit
    #[argh(switch, short = 'V')]
    pub version: bool,

    /// path to FASTA schema directory or schema file
    #[argh(option)]
    pub schema: Option<String>,

    /// path to allelic profile matrix (.tsv or .csv)
    #[argh(option)]
    pub profiles: Option<String>,

    /// output distance matrix file
    #[argh(option)]
    pub output: Option<String>,

    /// distance mode: snps, snps-indel-contiguous, snps-indel-bases, hamming (default: snps).
    /// Legacy alias: snps-indel-events == snps-indel-contiguous (deprecated).
    #[argh(option, default = "String::from(\"snps\")")]
    pub mode: String,

    /// output format: tsv, csv, phylip, nexus (default: tsv)
    #[argh(option, default = "String::from(\"tsv\")")]
    pub format: String,

    /// also emit a long-format per-pair table with data-quality columns
    /// (sample_i, sample_j, distance, shared_loci, total_loci, missing_frac) to this file
    #[argh(option)]
    pub emit_pairs: Option<String>,

    /// with --emit-pairs, add a missingness confidence interval per pair
    /// (dist_norm, ci_low, ci_high, ci_reliable) accounting for unobserved loci
    #[argh(switch)]
    pub report_ci: bool,

    /// confidence level for --report-ci (default: 0.95)
    #[argh(option, default = "0.95")]
    pub ci_level: f64,

    /// also write a self-contained HTML analyst dashboard (summary, per-sample
    /// quality, distance distribution, interactive outbreak-clustering explorer) to this file
    #[argh(option)]
    pub report: Option<String>,

    /// per-locus mutation-density threshold (percent) for the dashboard's
    /// recombination view; a locus above this density counts as recombinant.
    /// Requires an enriched cache (sequence lengths). Default: 3.0
    #[argh(option, default = "3.0")]
    pub recomb_threshold: f64,

    /// missing data character (default: -)
    #[argh(option, default = "String::from(\"-\")")]
    pub missing_char: String,

    /// minimum number of shared loci required for distance calculation (default: 0)
    #[argh(option, default = "0")]
    pub min_loci: usize,

    /// number of threads (default: 1; pass 0 for auto-detect = number of physical cores).
    /// On shared systems, auto-detect can interfere with other processes; users running on
    /// dedicated hardware should explicitly request the desired number of threads.
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

    /// enable Hamming fallback for SNPs-only mode (opt-in: when an allele pair has 0 SNPs but
    /// different hashes due to InDels, contribute +1 instead of 0; preserves cgDist >= Hamming
    /// ordering at the cost of counting non-SNP positions as "SNPs"). Default: disabled.
    #[argh(switch)]
    pub hamming_fallback: bool,

    /// [DEPRECATED] no longer needed: Hamming fallback is now opt-in via --hamming-fallback.
    /// This flag is accepted for backward compatibility and is a no-op (a warning is printed
    /// when supplied).
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

    /// weights for --mode custom: per-locus contribution = sum of
    /// weight*count, e.g. "nonsyn=1,frame_disrupted=1,indel_events=1"; keys:
    /// allele (1 per differing locus), snps, indel_events, indel_bases, syn,
    /// nonsyn, frame_disrupted; non-negative integers
    #[argh(option)]
    pub weights: Option<String>,

    /// NCBI translation table for synonymous/nonsynonymous classification
    /// (default: 11, Bacterial, Archaeal and Plant Plastid)
    #[argh(option, default = "11")]
    pub translation_table: u32,

    /// do not read an alternative start codon (e.g. GTG, TTG) as Met when it
    /// is the first codon
    #[argh(switch)]
    pub no_first_codon_as_met: bool,

    /// compute and store synonymous/nonsynonymous SNP counts in the cache
    /// even when the distance mode does not use them
    #[argh(switch)]
    pub coding_stats: bool,

    /// save one compact CIGAR row per aligned pair (TSV: locus, hash1, hash2,
    /// cigar, snps, indel_events, indel_bases, alignment_score); '=' match,
    /// 'X' SNP, 'I' query base vs reference gap, 'D' reference base vs query gap
    #[argh(option)]
    pub save_cigar: Option<String>,

    /// re-check this fraction (0-1) of new alignments against parasail's
    /// original kernel and stop with an error on any difference; the pairs
    /// are chosen deterministically (default: 0 = off, 1 = all)
    #[argh(option, default = "0.0")]
    pub verify_alignments: f64,

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
