// detector.rs - Recombination detection using enriched cache data

use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use serde::{Serialize, Deserialize};
use crate::cache::ModernCache;
use crate::hashers::AlleleHashPair;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RecombinationEvent {
    pub sample_i: String,
    pub sample_j: String,
    pub locus: String,
    pub cgdist_distance: f64,
    pub hamming_distance: f64,
    pub recombination_rate: f64,
    pub avg_sequence_length: f64,
    pub mutation_density: f64,
}

#[derive(Debug)]
pub struct RecombinationDetector {
    pub threshold_percent: f64,
    pub recombination_events: Vec<RecombinationEvent>,
}

impl RecombinationDetector {
    pub fn new(threshold_percent: f64) -> Self {
        Self {
            threshold_percent,
            recombination_events: Vec::new(),
        }
    }

    /// Analyze recombination from enriched cache and distance matrix
    pub fn detect_recombination(
        &mut self,
        enriched_cache: &ModernCache,
        distance_matrix_path: &str,
        sample_names: &[String],
    ) -> Result<(), String> {
        println!("🧬 Starting recombination detection with {}% threshold", self.threshold_percent);
        
        // Load distance matrix (assuming TSV format)
        let distance_matrix = self.load_distance_matrix(distance_matrix_path, sample_names)?;
        
        // Group cache entries by locus
        let mut locus_entries: HashMap<String, Vec<_>> = HashMap::new();
        for (key, entry) in &enriched_cache.entries {
            locus_entries.entry(key.locus.clone()).or_default().push((key, entry));
        }

        println!("📊 Processing {} loci for recombination analysis", locus_entries.len());
        
        let mut total_comparisons = 0;
        let mut recombination_events = 0;

        // Analyze each locus
        for (locus, entries) in locus_entries {
            for (cache_key, cache_entry) in entries {
                // Find corresponding samples in the distance matrix
                if let (Some(i), Some(j)) = (
                    self.find_sample_by_hash(&cache_key.hash_pair.hash1, &enriched_cache),
                    self.find_sample_by_hash(&cache_key.hash_pair.hash2, &enriched_cache)
                ) {
                    total_comparisons += 1;
                    
                    // Calculate average sequence length
                    let len1 = cache_entry.seq1_length.unwrap_or(0) as f64;
                    let len2 = cache_entry.seq2_length.unwrap_or(0) as f64;
                    
                    if len1 == 0.0 || len2 == 0.0 {
                        continue; // Skip if no length information
                    }
                    
                    let avg_length = (len1 + len2) / 2.0;
                    
                    // Calculate cgDist distance (SNPs + indels)
                    let cgdist_distance = cache_entry.snps as f64 + 
                                        cache_entry.indel_events as f64 + 
                                        cache_entry.indel_bases as f64;
                    
                    // Calculate Hamming distance (would be 1 if alleles differ, 0 if same)
                    let hamming_distance = if cache_key.hash_pair.hash1 != cache_key.hash_pair.hash2 { 1.0 } else { 0.0 };
                    
                    // Calculate mutation density (mutations per nucleotide)
                    let mutation_density = cgdist_distance / avg_length;
                    
                    // Check if exceeds threshold
                    let threshold_density = self.threshold_percent / 100.0;
                    
                    if mutation_density > threshold_density && hamming_distance > 0.0 {
                        let recombination_rate = (mutation_density - threshold_density) / threshold_density * 100.0;
                        
                        let event = RecombinationEvent {
                            sample_i: i.clone(),
                            sample_j: j.clone(),
                            locus: locus.clone(),
                            cgdist_distance,
                            hamming_distance,
                            recombination_rate,
                            avg_sequence_length: avg_length,
                            mutation_density: mutation_density * 100.0, // Convert to percentage
                        };
                        
                        self.recombination_events.push(event);
                        recombination_events += 1;
                    }
                }
            }
        }

        println!("✅ Recombination analysis complete:");
        println!("   Total comparisons: {}", total_comparisons);
        println!("   Recombination events detected: {}", recombination_events);
        println!("   Recombination rate: {:.2}%", 
               if total_comparisons > 0 { recombination_events as f64 / total_comparisons as f64 * 100.0 } else { 0.0 });

        Ok(())
    }

    /// Generate corrected distance matrix with hamming fallback for recombination events
    pub fn generate_corrected_matrix(
        &self,
        original_matrix_path: &str,
        output_matrix_path: &str,
        sample_names: &[String],
    ) -> Result<(), String> {
        println!("📝 Generating corrected distance matrix with hamming fallback");
        
        // Load original matrix
        let mut distance_matrix = self.load_distance_matrix(original_matrix_path, sample_names)?;
        
        // Apply hamming fallback to recombination events
        let mut corrections = 0;
        for event in &self.recombination_events {
            if let (Some(i), Some(j)) = (
                sample_names.iter().position(|s| s == &event.sample_i),
                sample_names.iter().position(|s| s == &event.sample_j)
            ) {
                distance_matrix[i][j] = event.hamming_distance;
                distance_matrix[j][i] = event.hamming_distance; // Symmetric matrix
                corrections += 1;
            }
        }

        // Write corrected matrix
        self.write_distance_matrix(&distance_matrix, output_matrix_path, sample_names)?;
        
        println!("✅ Applied {} corrections to distance matrix", corrections);
        
        Ok(())
    }

    /// Generate detailed recombination log
    pub fn write_recombination_log(&self, log_path: &str) -> Result<(), String> {
        let mut file = File::create(log_path)
            .map_err(|e| format!("Failed to create log file: {}", e))?;

        writeln!(file, "sample_i\tsample_j\tlocus\tcgdist_distance\thamming_distance\trecombination_rate_%\tavg_seq_length\tmutation_density_%")
            .map_err(|e| format!("Failed to write header: {}", e))?;

        for event in &self.recombination_events {
            writeln!(file, "{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.0}\t{:.4}",
                event.sample_i,
                event.sample_j,
                event.locus,
                event.cgdist_distance,
                event.hamming_distance,
                event.recombination_rate,
                event.avg_sequence_length,
                event.mutation_density
            ).map_err(|e| format!("Failed to write event: {}", e))?;
        }

        println!("📋 Wrote {} recombination events to {}", self.recombination_events.len(), log_path);
        Ok(())
    }

    // Helper functions
    fn load_distance_matrix(&self, path: &str, sample_names: &[String]) -> Result<Vec<Vec<f64>>, String> {
        // Implementation for loading TSV/CSV distance matrix
        // For now, return empty matrix as placeholder
        let n = sample_names.len();
        Ok(vec![vec![0.0; n]; n])
    }

    fn write_distance_matrix(&self, matrix: &[Vec<f64>], path: &str, sample_names: &[String]) -> Result<(), String> {
        let mut file = File::create(path)
            .map_err(|e| format!("Failed to create matrix file: {}", e))?;

        // Write header
        write!(file, "").map_err(|e| format!("Write error: {}", e))?;
        for sample in sample_names {
            write!(file, "\t{}", sample).map_err(|e| format!("Write error: {}", e))?;
        }
        writeln!(file).map_err(|e| format!("Write error: {}", e))?;

        // Write matrix
        for (i, row) in matrix.iter().enumerate() {
            write!(file, "{}", sample_names[i]).map_err(|e| format!("Write error: {}", e))?;
            for &value in row {
                write!(file, "\t{:.6}", value).map_err(|e| format!("Write error: {}", e))?;
            }
            writeln!(file).map_err(|e| format!("Write error: {}", e))?;
        }

        Ok(())
    }

    fn find_sample_by_hash(&self, hash: &crate::hashers::AlleleHash, cache: &ModernCache) -> Option<String> {
        // This is a placeholder - we would need to map CRC hashes back to sample names
        // This requires additional data structure or reverse lookup
        None
    }
}