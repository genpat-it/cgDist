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
        if let Some(val) = config.hasher_type {
            if self.hasher_type == "crc32" {
                self.hasher_type = val;
            }
        }
        if let Some(val) = config.mode {
            if self.mode == "snps" {
                self.mode = val;
            }
        }
        if let Some(val) = config.format {
            if self.format == "tsv" {
                self.format = val;
            }
        }
        if let Some(val) = config.missing_char {
            if self.missing_char == "-" {
                self.missing_char = val;
            }
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
        if let Some(val) = config.sample_threshold {
            if self.sample_threshold == 0.0 {
                self.sample_threshold = val;
            }
        }
        if let Some(val) = config.locus_threshold {
            if self.locus_threshold == 0.0 {
                self.locus_threshold = val;
            }
        }
        if let Some(val) = config.min_loci {
            if self.min_loci == 0 {
                self.min_loci = val;
            }
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
        if let Some(val) = config.alignment_mode {
            if self.alignment_mode == "dna" {
                self.alignment_mode = val;
            }
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

        // Flags (CLI flags take precedence, config only sets if not explicitly set).
        // hamming_fallback is opt-in; CLI flag wins, otherwise honour the config value
        // (default false).
        if !self.hamming_fallback && config.hamming_fallback.unwrap_or(false) {
            self.hamming_fallback = true;
        }
        // [DEPRECATED] config.no_hamming_fallback is accepted but no longer alters behaviour
        // (the new default is "no fallback"). It is read here only to suppress the unused-field
        // warning; a deprecation note is printed in main.rs when the legacy CLI flag is set.
        let _ = config.no_hamming_fallback;
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
