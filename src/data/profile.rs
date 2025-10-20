// profile.rs - Allelic profile and matrix data structures

use crate::hashers::AlleleHash;
use rayon::prelude::*;
use regex::Regex;
use std::collections::{HashMap, HashSet};

/// Genetic diversity metrics for the allelic matrix
#[derive(Debug, Clone)]
pub struct DiversityMetrics {
    pub avg_unique_alleles: f64,
    pub diversity_index: f64,
    pub diversity_category: &'static str,
    pub total_unique_pairs: usize,
}

/// Represents a single sample's allelic profile
#[derive(Debug)]
pub struct AllelicProfile {
    pub sample_id: String,
    pub loci_hashes: HashMap<String, AlleleHash>,
}

/// Collection of allelic profiles with associated metadata
#[derive(Debug)]
pub struct AllelicMatrix {
    pub samples: Vec<AllelicProfile>,
    pub loci_names: Vec<String>,
}

impl Default for AllelicMatrix {
    fn default() -> Self {
        Self::new()
    }
}

impl AllelicMatrix {
    /// Create a new empty matrix
    pub fn new() -> Self {
        Self {
            samples: Vec::new(),
            loci_names: Vec::new(),
        }
    }

    /// Load matrix from file with specified hasher
    #[allow(clippy::too_many_arguments)]
    pub fn from_file_with_hasher(
        file_path: &std::path::Path,
        missing_char: &str,
        hasher_type: &str,
        sample_threshold: f64,
        locus_threshold: f64,
        sample_include: Option<&Regex>,
        sample_exclude: Option<&Regex>,
        samples_include: Option<&HashSet<String>>,
        samples_exclude: Option<&HashSet<String>>,
        loci_include: Option<&Regex>,
        loci_exclude: Option<&Regex>,
        loci_include_set: Option<&HashSet<String>>,
        loci_exclude_set: Option<&HashSet<String>>,
    ) -> Result<Self, String> {
        println!(
            "📊 Loading allelic matrix with {} hasher: {}",
            hasher_type.to_uppercase(),
            file_path.display()
        );

        // Get the hasher from registry
        let registry = crate::hashers::HasherRegistry::new();
        let hasher = registry
            .get_hasher(hasher_type)
            .ok_or_else(|| format!("Unknown hasher type: {}", hasher_type))?;

        let extension = file_path
            .extension()
            .and_then(|s| s.to_str())
            .unwrap_or("tsv");

        // Load the matrix with the specific hasher
        let mut matrix = match extension {
            "csv" => Self::from_csv_with_hasher(file_path, missing_char, hasher)?,
            _ => Self::from_tsv_with_hasher(file_path, missing_char, hasher)?,
        };

        // Show initial statistics
        Self::print_matrix_statistics(&matrix, "INITIAL MATRIX");

        // Apply filtering
        matrix.apply_sample_filtering(
            sample_include,
            sample_exclude,
            samples_include,
            samples_exclude,
        )?;
        if sample_include.is_some()
            || sample_exclude.is_some()
            || samples_include.is_some()
            || samples_exclude.is_some()
        {
            Self::print_matrix_statistics(&matrix, "AFTER SAMPLE FILTERING");
        }

        matrix.apply_quality_filters(
            sample_threshold,
            locus_threshold,
            loci_include,
            loci_exclude,
            loci_include_set,
            loci_exclude_set,
        )?;
        Self::print_matrix_statistics(&matrix, "FINAL MATRIX");

        Ok(matrix)
    }

    /// Sample filtering function
    pub fn apply_sample_filtering(
        &mut self,
        sample_include: Option<&Regex>,
        sample_exclude: Option<&Regex>,
        samples_include: Option<&HashSet<String>>,
        samples_exclude: Option<&HashSet<String>>,
    ) -> Result<(), String> {
        let initial_samples = self.samples.len();

        // Apply sample filtering
        if sample_include.is_some()
            || sample_exclude.is_some()
            || samples_include.is_some()
            || samples_exclude.is_some()
        {
            self.samples.retain(|sample| {
                let sample_id = &sample.sample_id;

                // Include regex filter
                if let Some(regex) = sample_include {
                    if !regex.is_match(sample_id) {
                        return false;
                    }
                }

                // Exclude regex filter
                if let Some(regex) = sample_exclude {
                    if regex.is_match(sample_id) {
                        return false;
                    }
                }

                // Include set filter
                if let Some(set) = samples_include {
                    if !set.contains(sample_id) {
                        return false;
                    }
                }

                // Exclude set filter
                if let Some(set) = samples_exclude {
                    if set.contains(sample_id) {
                        return false;
                    }
                }

                true
            });

            let filtered_samples = self.samples.len();
            if initial_samples != filtered_samples {
                println!(
                    "Sample filters: kept {} samples (removed {})",
                    filtered_samples,
                    initial_samples - filtered_samples
                );
            }
        }

        Ok(())
    }

    /// Apply quality filters based on completeness thresholds
    pub fn apply_quality_filters(
        &mut self,
        sample_threshold: f64,
        locus_threshold: f64,
        loci_include: Option<&Regex>,
        loci_exclude: Option<&Regex>,
        loci_include_set: Option<&HashSet<String>>,
        loci_exclude_set: Option<&HashSet<String>>,
    ) -> Result<(), String> {
        println!("\n=== APPLYING QUALITY FILTERS ===");

        let original_samples = self.samples.len();
        let original_loci = self.loci_names.len();

        // Apply loci filtering first
        let mut loci_to_keep = self.loci_names.clone();

        // Regex filters
        if let Some(include_regex) = loci_include {
            loci_to_keep.retain(|locus| include_regex.is_match(locus));
        }
        if let Some(exclude_regex) = loci_exclude {
            loci_to_keep.retain(|locus| !exclude_regex.is_match(locus));
        }

        // Set filters
        if let Some(include_set) = loci_include_set {
            loci_to_keep.retain(|locus| include_set.contains(locus));
        }
        if let Some(exclude_set) = loci_exclude_set {
            loci_to_keep.retain(|locus| !exclude_set.contains(locus));
        }

        self.loci_names = loci_to_keep;
        let before = original_loci;
        let after = self.loci_names.len();
        if before != after {
            println!(
                "Loci filters: kept {} loci (removed {})",
                after,
                before - after
            );
        }

        // Update sample hashes in parallel
        let loci_names_set: HashSet<_> = self.loci_names.iter().cloned().collect();
        self.samples.par_iter_mut().for_each(|sample| {
            sample
                .loci_hashes
                .retain(|locus, _| loci_names_set.contains(locus));
        });

        // Apply sample quality threshold
        if sample_threshold > 0.0 {
            let removed_samples = self.samples.len();
            self.samples.retain(|sample| {
                let total_loci = self.loci_names.len();
                let present_loci = sample
                    .loci_hashes
                    .values()
                    .filter(|hash| !hash.is_missing())
                    .count();
                let completeness = present_loci as f64 / total_loci as f64;
                completeness >= sample_threshold
            });
            let removed_samples = removed_samples - self.samples.len();
            if removed_samples > 0 {
                println!(
                    "Sample completeness filter (threshold {:.1}%): removed {} samples",
                    sample_threshold * 100.0,
                    removed_samples
                );
            }
        }

        // Apply locus quality threshold
        if locus_threshold > 0.0 {
            let loci_to_keep: Vec<String> = self
                .loci_names
                .iter()
                .filter(|locus_name| {
                    let total_samples = self.samples.len();
                    let present_samples = self
                        .samples
                        .iter()
                        .filter(|sample| {
                            sample
                                .loci_hashes
                                .get(*locus_name)
                                .map(|hash| !hash.is_missing())
                                .unwrap_or(false)
                        })
                        .count();
                    let completeness = present_samples as f64 / total_samples as f64;
                    completeness >= locus_threshold
                })
                .cloned()
                .collect();

            if loci_to_keep.len() != self.loci_names.len() {
                let removed_loci = self.loci_names.len() - loci_to_keep.len();
                println!(
                    "Locus completeness filter (threshold {:.1}%): removed {} loci",
                    locus_threshold * 100.0,
                    removed_loci
                );
                self.loci_names = loci_to_keep;
                let loci_names_set: HashSet<_> = self.loci_names.iter().cloned().collect();
                self.samples.par_iter_mut().for_each(|sample| {
                    sample
                        .loci_hashes
                        .retain(|locus, _| loci_names_set.contains(locus));
                });
            }
        }

        let final_samples = self.samples.len();
        let final_loci = self.loci_names.len();
        println!("Filter summary:");
        println!(
            "  Samples: {} → {} (removed {})",
            original_samples,
            final_samples,
            original_samples - final_samples
        );
        println!(
            "  Loci: {} → {} (removed {})",
            original_loci,
            final_loci,
            original_loci - final_loci
        );

        if final_samples == 0 {
            return Err("No samples remain after filtering".to_string());
        }
        if final_loci == 0 {
            return Err("No loci remain after filtering".to_string());
        }

        Ok(())
    }

    /// Print matrix statistics
    pub fn print_matrix_statistics(&self, phase: &str) {
        println!("\n📊 === MATRIX STATISTICS ({}) ===", phase);
        let total_cells = self.samples.len() * self.loci_names.len();

        // Calculate missing data in parallel
        let sample_missing_counts: Vec<_> = self
            .samples
            .par_iter()
            .map(|sample| {
                let sample_missing = self
                    .loci_names
                    .iter()
                    .filter(|locus| {
                        sample
                            .loci_hashes
                            .get(*locus)
                            .map(|h| h.is_missing())
                            .unwrap_or(true)
                    })
                    .count();
                let missing_percent = 100.0 * sample_missing as f64 / self.loci_names.len() as f64;
                (sample.sample_id.clone(), missing_percent, sample_missing)
            })
            .collect();

        let locus_missing_counts: Vec<_> = self
            .loci_names
            .par_iter()
            .map(|locus| {
                let locus_missing = self
                    .samples
                    .iter()
                    .filter(|sample| {
                        sample
                            .loci_hashes
                            .get(locus)
                            .map(|h| h.is_missing())
                            .unwrap_or(true)
                    })
                    .count();
                let missing_percent = 100.0 * locus_missing as f64 / self.samples.len() as f64;
                (locus.clone(), missing_percent, locus_missing)
            })
            .collect();

        println!(
            "  📏 Dimensions: {} samples × {} loci = {} total cells",
            self.samples.len(),
            self.loci_names.len(),
            total_cells
        );

        // Global missing percentage
        let total_missing: usize = sample_missing_counts
            .iter()
            .map(|(_, _, missing)| missing)
            .sum();
        let global_missing_percent = 100.0 * total_missing as f64 / total_cells as f64;

        print!(
            "  📊 Missing data: {:.2}% ({} cells)",
            global_missing_percent, total_missing
        );
        if global_missing_percent <= 5.0 {
            println!("  🟢 EXCELLENT: Very low missing data");
        } else if global_missing_percent <= 15.0 {
            println!("  🟡 GOOD: Acceptable missing data");
        } else if global_missing_percent <= 30.0 {
            println!(
                "  🟠 FAIR: High missing data ({:.2}%) - consider quality filters",
                global_missing_percent
            );
        } else {
            println!(
                "  🔴 POOR: Very high missing data ({:.2}%) - quality filters recommended",
                global_missing_percent
            );
        }

        let complete_samples = sample_missing_counts
            .iter()
            .filter(|(_, missing, _)| *missing == 0.0)
            .count();
        let complete_loci = locus_missing_counts
            .iter()
            .filter(|(_, missing, _)| *missing == 0.0)
            .count();

        println!(
            "  ✅ Complete samples: {} ({:.1}%)",
            complete_samples,
            100.0 * complete_samples as f64 / self.samples.len() as f64
        );
        println!(
            "  ✅ Complete loci: {} ({:.1}%)",
            complete_loci,
            100.0 * complete_loci as f64 / self.loci_names.len() as f64
        );

        // Calculate genetic diversity metrics
        let diversity_stats = self.calculate_diversity_metrics();
        println!(
            "  🧬 Avg unique alleles per locus: {:.1}",
            diversity_stats.avg_unique_alleles
        );
        println!(
            "  📈 Allelic diversity index: {:.3} ({}/{})",
            diversity_stats.diversity_index,
            diversity_stats.diversity_category,
            if diversity_stats.diversity_index < 0.3 {
                "🟢 Low diversity - clonal population"
            } else if diversity_stats.diversity_index < 0.6 {
                "🟡 Moderate diversity"
            } else {
                "🔴 High diversity - very heterogeneous"
            }
        );
        println!(
            "  🔢 Total unique CRC pairs: {}",
            diversity_stats.total_unique_pairs
        );
    }

    /// Calculate genetic diversity metrics for the matrix
    pub fn calculate_diversity_metrics(&self) -> DiversityMetrics {
        // Calculate unique alleles per locus in parallel
        let locus_unique_counts: Vec<usize> = self
            .loci_names
            .par_iter()
            .map(|locus| {
                let mut unique_alleles = HashSet::new();

                // Collect all non-missing alleles for this locus
                for sample in &self.samples {
                    if let Some(hash) = sample.loci_hashes.get(locus) {
                        if !hash.is_missing() {
                            if let Some(crc) = hash.as_crc32() {
                                unique_alleles.insert(crc);
                            }
                        }
                    }
                }

                unique_alleles.len()
            })
            .collect();

        // Calculate statistics
        let total_unique_alleles: usize = locus_unique_counts.iter().sum();
        let avg_unique_alleles = total_unique_alleles as f64 / self.loci_names.len() as f64;

        // Calculate diversity index (normalized average unique alleles per locus)
        // Ranges from 0 (all identical) to 1 (maximum diversity)
        let max_possible_unique = self.samples.len();
        let diversity_index = avg_unique_alleles / max_possible_unique as f64;

        // Categorize diversity
        let diversity_category = if diversity_index < 0.3 {
            "Low"
        } else if diversity_index < 0.6 {
            "Moderate"
        } else {
            "High"
        };

        // Calculate total unique pairs (approximation)
        let total_unique_pairs: usize = locus_unique_counts
            .iter()
            .map(|&unique_count| {
                if unique_count > 1 {
                    unique_count * (unique_count - 1) / 2 // C(n,2)
                } else {
                    0
                }
            })
            .sum();

        DiversityMetrics {
            avg_unique_alleles,
            diversity_index,
            diversity_category,
            total_unique_pairs,
        }
    }
}
