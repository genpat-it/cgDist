// merge.rs - Merge configuration file with CLI arguments

use crate::cli::{Args, Config};

impl Args {
    /// Merge with configuration from file
    /// CLI arguments take precedence over config file values
    pub fn merge_with_config(mut self, config: Config) -> Self {
        // Input/Output
        if self.profiles.is_none() {
            self.profiles = config.profiles;
        }
        if self.schema.is_none() {
            self.schema = config.schema;
        }
        if self.output.is_none() {
            self.output = config.output;
        }
        
        // Core settings (only override defaults, not explicit CLI values)
        if self.hasher_type == "crc32" && config.hasher_type.is_some() {
            self.hasher_type = config.hasher_type.unwrap();
        }
        if self.mode == "snps" && config.mode.is_some() {
            self.mode = config.mode.unwrap();
        }
        if self.format == "tsv" && config.format.is_some() {
            self.format = config.format.unwrap();
        }
        if self.missing_char == "-" && config.missing_char.is_some() {
            self.missing_char = config.missing_char.unwrap();
        }
        
        // Performance
        if self.threads.is_none() {
            self.threads = config.threads;
        }
        if self.cache_file.is_none() {
            self.cache_file = config.cache_file;
        }
        if self.cache_note.is_none() {
            self.cache_note = config.cache_note;
        }
        
        // Quality filters (only override default 0.0)
        if self.sample_threshold == 0.0 && config.sample_threshold.is_some() {
            self.sample_threshold = config.sample_threshold.unwrap();
        }
        if self.locus_threshold == 0.0 && config.locus_threshold.is_some() {
            self.locus_threshold = config.locus_threshold.unwrap();
        }
        if self.min_loci == 0 && config.min_loci.is_some() {
            self.min_loci = config.min_loci.unwrap();
        }
        
        // Sample/Loci filtering
        if self.include_samples.is_none() {
            self.include_samples = config.include_samples;
        }
        if self.exclude_samples.is_none() {
            self.exclude_samples = config.exclude_samples;
        }
        if self.include_loci.is_none() {
            self.include_loci = config.include_loci;
        }
        if self.exclude_loci.is_none() {
            self.exclude_loci = config.exclude_loci;
        }
        if self.include_loci_list.is_none() {
            self.include_loci_list = config.include_loci_list;
        }
        if self.exclude_loci_list.is_none() {
            self.exclude_loci_list = config.exclude_loci_list;
        }
        if self.include_samples_list.is_none() {
            self.include_samples_list = config.include_samples_list;
        }
        if self.exclude_samples_list.is_none() {
            self.exclude_samples_list = config.exclude_samples_list;
        }
        
        // Alignment settings (only override default "dna")
        if self.alignment_mode == "dna" && config.alignment_mode.is_some() {
            self.alignment_mode = config.alignment_mode.unwrap();
        }
        if self.match_score.is_none() {
            self.match_score = config.match_score;
        }
        if self.mismatch_penalty.is_none() {
            self.mismatch_penalty = config.mismatch_penalty;
        }
        if self.gap_open.is_none() {
            self.gap_open = config.gap_open;
        }
        if self.gap_extend.is_none() {
            self.gap_extend = config.gap_extend;
        }
        
        // Flags (CLI flags take precedence, config only sets if not explicitly set)
        if !self.no_hamming_fallback && config.no_hamming_fallback.unwrap_or(false) {
            self.no_hamming_fallback = true;
        }
        if !self.force_recompute && config.force_recompute.unwrap_or(false) {
            self.force_recompute = true;
        }
        if !self.dry_run && config.dry_run.unwrap_or(false) {
            self.dry_run = true;
        }
        if self.save_alignments.is_none() {
            self.save_alignments = config.save_alignments;
        }
        
        self
    }
    
    /// Load configuration and merge with CLI args
    pub fn with_config_file(self, config_path: &str) -> Result<Self, String> {
        let config = Config::from_file(config_path)?;
        Ok(self.merge_with_config(config))
    }
}