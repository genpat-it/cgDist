// distance.rs - Core distance calculation engine

use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::time::Instant;
use rayon::prelude::*;
use rayon::iter::ParallelIterator;
use parasail_rs::{Aligner, Matrix};
use indicatif::{ProgressBar, ProgressStyle};
use crate::core::alignment::{DistanceMode, AlignmentConfig, compute_alignment_stats};
use crate::hashers::{HasherRegistry, AlleleHasher};
use crate::data::{AllelicProfile, SequenceDatabase};
use serde::{Serialize, Deserialize};
use chrono;

/// Internal cache entry storing all alignment statistics
#[derive(Debug, Clone, Copy)]
struct CacheEntry {
    snps: usize,
    indel_events: usize,
    indel_bases: usize,
}

/// Distance calculation engine
pub struct DistanceEngine {
    cache: HashMap<DistanceCacheKey, CacheEntry>,
    config: AlignmentConfig,
    sequence_db: Option<SequenceDatabase>,
    hasher_type: String,
    cache_note: Option<String>,
    has_new_entries: bool,
    // For saving detailed alignments
    save_alignments_path: Option<String>,
    alignment_details: Vec<String>, // TSV lines to write
}

/// Modern cache structure supporting any hasher type
#[derive(Debug, Serialize, Deserialize)]
pub struct ModernCache {
    pub data: HashMap<String, CacheValue>,  // Use string keys for JSON compatibility
    pub metadata: CacheMetadata,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CacheKey {
    pub locus: String,
    pub hash1: String,  // String to support any hasher (CRC32, SHA256, MD5, etc.)
    pub hash2: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CacheValue {
    pub snps: usize,
    pub indel_events: usize,
    pub indel_bases: usize,
    pub computed_at: String,
    // New fields for sequence lengths
    #[serde(skip_serializing_if = "Option::is_none")]
    pub seq1_length: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub seq2_length: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CacheMetadata {
    pub version: String,
    pub created: String,
    pub last_modified: String,
    pub alignment_config: AlignmentConfig,
    pub hasher_type: String,
    pub distance_mode: String,
    pub user_note: Option<String>,
    pub total_entries: usize,
    pub unique_loci: usize,
    pub format_version: u32,
}

// Legacy support for old inspector format
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
struct DistanceCacheKey {
    locus: String,
    crc1: u32,
    crc2: u32,
}

impl DistanceEngine {
    pub fn new(config: AlignmentConfig, hasher_type: String) -> Self {
        Self {
            cache: HashMap::new(),
            config,
            sequence_db: None,
            hasher_type,
            cache_note: None,
            has_new_entries: false,
            save_alignments_path: None,
            alignment_details: Vec::new(),
        }
    }
    
    pub fn with_sequences(config: AlignmentConfig, sequence_db: SequenceDatabase, hasher_type: String) -> Self {
        Self {
            cache: HashMap::new(),
            config,
            sequence_db: Some(sequence_db),
            hasher_type,
            cache_note: None,
            has_new_entries: false,
            save_alignments_path: None,
            alignment_details: Vec::new(),
        }
    }
    
    /// Set a user note for the cache
    pub fn set_cache_note(&mut self, note: String) {
        self.cache_note = Some(note);
    }
    
    /// Get distance between two CRCs for a specific locus (optimized)
    pub fn get_distance(
        &self,
        locus: &str,
        crc1: u32,
        crc2: u32,
        mode: DistanceMode,
        no_hamming_fallback: bool,
    ) -> usize {
        // Fast path for missing data
        if crc1 == u32::MAX || crc2 == u32::MAX {
            return 0; // Missing data contributes 0 to distance
        }
        
        // Fast path for identical alleles
        if crc1 == crc2 {
            return 0; // Identical alleles
        }
        
        // Fast path for hamming hasher
        if self.hasher_type == "hamming" {
            return 1; // Different CRCs = distance 1 for Hamming
        }
        
        // Optimized key lookup - temporarily create key for lookup only
        let (min_crc, max_crc) = if crc1 <= crc2 { (crc1, crc2) } else { (crc2, crc1) };
        
        // Use a temporary key for lookup - still allocates but only on cache miss
        let temp_key = DistanceCacheKey {
            locus: locus.to_string(),
            crc1: min_crc,
            crc2: max_crc,
        };
        
        if let Some(&entry) = self.cache.get(&temp_key) {
            let distance = match mode {
                DistanceMode::SnpsOnly => entry.snps,
                DistanceMode::SnpsAndIndelEvents => entry.snps + entry.indel_events,
                DistanceMode::SnpsAndIndelBases => entry.snps + entry.indel_bases,
                DistanceMode::Hamming => 1, // For hamming mode, different CRCs = 1
            };
            
            // Apply Hamming fallback ONLY for SNPs mode: if alignment found 0 differences 
            // but CRCs are different, return 1 to maintain consistency (different CRCs >= 1 difference)
            if distance == 0 && crc1 != crc2 && !no_hamming_fallback && mode == DistanceMode::SnpsOnly {
                return 1; // Hamming fallback: different CRCs = at least 1 difference
            }
            
            return distance;
        }
        
        // Cache miss - this should not happen if pre-computation worked
        // (Silent - individual misses not logged to reduce verbosity)
        
        // Cache miss - apply Hamming fallback ONLY for SNPs mode and when enabled
        if no_hamming_fallback {
            0 // Different CRCs with no alignment count as 0
        } else if mode == DistanceMode::SnpsOnly {
            1 // Hamming fallback for SNPs: different CRCs count as 1
        } else {
            // For indel modes, cache miss means we can't compute meaningful distance
            // This should be rare if sequences are available in schema
            0 // Conservative: assume no difference if we can't align
        }
    }
    
    /// Add distance to cache (optimized) - stores all alignment statistics
    pub fn cache_distance(&mut self, locus: &str, crc1: u32, crc2: u32, snps: usize, indel_events: usize, indel_bases: usize) {
        let (min_crc, max_crc) = if crc1 <= crc2 { (crc1, crc2) } else { (crc2, crc1) };
        let key = DistanceCacheKey {
            locus: locus.to_string(), // Only allocate when inserting into cache
            crc1: min_crc,
            crc2: max_crc,
        };
        let entry = CacheEntry { snps, indel_events, indel_bases };
        self.cache.insert(key, entry);
        self.has_new_entries = true; // Mark that cache has new entries
    }
    
    /// Get cache statistics
    pub fn cache_stats(&self) -> (usize, usize) {
        (self.cache.len(), self.cache.capacity())
    }
    
    /// Check if cache has new entries since last save/load
    pub fn has_new_entries(&self) -> bool {
        self.has_new_entries
    }
    
    /// Save cache to LZ4 compressed file (modern format)
    pub fn save_cache(&mut self, cache_path: &str, distance_mode: DistanceMode) -> Result<(), String> {
        println!("💾 Saving cache to {}...", cache_path);
        let start = Instant::now();
        
        // Convert internal cache to modern format
        let mut data = HashMap::new();
        let now = chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC").to_string();
        
        // Count unique loci efficiently without cloning
        let unique_loci: std::collections::HashSet<&str> = self.cache.keys()
            .map(|k| k.locus.as_str())
            .collect();
        
        for (key, entry) in &self.cache {
            // Create a string key with pre-allocated capacity to avoid reallocations
            let mut string_key = String::with_capacity(key.locus.len() + 24); // locus + ":4294967295:4294967295"
            string_key.push_str(&key.locus);
            string_key.push(':');
            string_key.push_str(&key.crc1.to_string());
            string_key.push(':');
            string_key.push_str(&key.crc2.to_string());
            
            let cache_value = CacheValue {
                snps: entry.snps,
                indel_events: entry.indel_events,
                indel_bases: entry.indel_bases,
                computed_at: now.clone(), // Keep this clone as now is reused
                seq1_length: None,  // Will be enriched later if requested
                seq2_length: None,  // Will be enriched later if requested
            };
            
            data.insert(string_key, cache_value);
        }
        
        let metadata = CacheMetadata {
            version: env!("CARGO_PKG_VERSION").to_string(),
            created: now.clone(),
            last_modified: now, // Move now instead of clone (last usage)
            alignment_config: self.config.clone(), // Keep clone - config is reused
            hasher_type: self.hasher_type.clone(), // Keep clone - hasher_type is reused
            distance_mode: match distance_mode {
                DistanceMode::SnpsOnly => "snps".to_string(),
                DistanceMode::SnpsAndIndelEvents => "snps-indel-events".to_string(),
                DistanceMode::SnpsAndIndelBases => "snps-indel-bases".to_string(),
                DistanceMode::Hamming => "hamming".to_string(),
            },
            user_note: self.cache_note.clone(),
            total_entries: self.cache.len(),
            unique_loci: unique_loci.len(),
            format_version: 2,  // Version 2 = modern format
        };
        
        let modern_cache = ModernCache { data, metadata };
        
        // Serialize cache with serde_json (more readable than bincode)
        let cache_data = serde_json::to_vec(&modern_cache)
            .map_err(|e| format!("Failed to serialize cache: {}", e))?;
        
        // Compress with LZ4
        let compressed = lz4_flex::compress_prepend_size(&cache_data);
        
        // Write to file
        std::fs::write(cache_path, &compressed)
            .map_err(|e| format!("Failed to write cache file: {}", e))?;
        
        let elapsed = start.elapsed();
        println!("✅ Cache saved in {:.2}s ({} entries, {} KB)", 
                 elapsed.as_secs_f64(), 
                 self.cache.len(),
                 compressed.len() / 1024);
        
        if let Some(note) = &self.cache_note {
            println!("📝 User note: {}", note);
        }
        
        // Reset the new entries flag since we just saved to disk
        self.has_new_entries = false;
        
        Ok(())
    }
    
    /// Quick compatibility check without loading full cache  
    pub fn check_cache_compatibility(&self, cache_path: &str, _distance_mode: DistanceMode) -> Result<(), String> {
        use std::fs::File;
        use std::io::Read;
        
        // Read only first 32KB which should contain metadata
        let mut file = File::open(cache_path)
            .map_err(|e| format!("Failed to open cache file: {}", e))?;
        
        let mut buffer = vec![0u8; 32768]; // 32KB should be enough for headers
        let bytes_read = file.read(&mut buffer)
            .map_err(|e| format!("Failed to read cache header: {}", e))?;
        buffer.truncate(bytes_read);
        
        // Try to find LZ4 size header (first 4 bytes = uncompressed size)
        if buffer.len() < 4 {
            return Err("Cache file too small".to_string());
        }
        
        let _uncompressed_size = u32::from_le_bytes([buffer[0], buffer[1], buffer[2], buffer[3]]);
        
        // Quick string search in compressed data for alignment config patterns
        let buffer_str = String::from_utf8_lossy(&buffer);
        
        // If metadata is not in the compressed header, skip quick check
        // This happens when JSON metadata is later in the file
        if !buffer_str.contains("alignment_config") && !buffer_str.contains("match_score") {
            // Can't do quick check - let the full load handle compatibility
            return Ok(());
        }
        
        // Check for specific mismatches between cache and current config
        if buffer_str.contains("match_score\":2") && self.config.match_score == 3 {
            return Err("Cache alignment config mismatch: cache uses match_score=2, current uses match_score=3".to_string());
        }
        if buffer_str.contains("match_score\":3") && self.config.match_score == 2 {
            return Err("Cache alignment config mismatch: cache uses match_score=3, current uses match_score=2".to_string());
        }
        
        // Check mismatch penalty
        if buffer_str.contains("mismatch_penalty\":-1") && self.config.mismatch_penalty == -2 {
            return Err("Cache alignment config mismatch: cache uses mismatch_penalty=-1, current uses mismatch_penalty=-2".to_string());
        }
        if buffer_str.contains("mismatch_penalty\":-2") && self.config.mismatch_penalty == -1 {
            return Err("Cache alignment config mismatch: cache uses mismatch_penalty=-2, current uses mismatch_penalty=-1".to_string());
        }
        
        // Check gap penalties
        if buffer_str.contains("gap_open\":5") && self.config.gap_open == 8 {
            return Err("Cache alignment config mismatch: cache uses gap_open=5, current uses gap_open=8".to_string());
        }
        if buffer_str.contains("gap_open\":8") && self.config.gap_open == 5 {
            return Err("Cache alignment config mismatch: cache uses gap_open=8, current uses gap_open=5".to_string());
        }
        
        Ok(())
    }

    /// Load cache from LZ4 compressed file (supports both modern and legacy formats)
    pub fn load_cache(&mut self, cache_path: &str, distance_mode: DistanceMode) -> Result<(), String> {
        println!("📂 Loading cache from {}...", cache_path);
        let start = Instant::now();
        
        // Read compressed file
        let compressed = std::fs::read(cache_path)
            .map_err(|e| format!("Failed to read cache file: {}", e))?;
        
        // Decompress with LZ4
        let decompressed = lz4_flex::decompress_size_prepended(&compressed)
            .map_err(|e| format!("Failed to decompress cache: {}", e))?;
        
        // Try modern format first (JSON)
        if let Ok(modern_cache) = serde_json::from_slice::<ModernCache>(&decompressed) {
            println!("🆕 Loading modern cache format (v{})", modern_cache.metadata.format_version);
            
            // Check alignment configuration compatibility
            if modern_cache.metadata.alignment_config != self.config {
                return Err(format!(
                    "Cache alignment config mismatch:\n  Cache: {:?}\n  Engine: {:?}",
                    modern_cache.metadata.alignment_config, self.config
                ));
            }
            
            // Check hasher compatibility
            if modern_cache.metadata.hasher_type != self.hasher_type {
                return Err(format!(
                    "Cache hasher type mismatch:\n  Cache: {}\n  Engine: {}",
                    modern_cache.metadata.hasher_type, self.hasher_type
                ));
            }
            
            // Check distance mode compatibility - TEMPORARILY DISABLED
            // The cache now contains all 3 values (SNPs, indel_events, indel_bases)
            // so it can be used with any distance mode safely
            let expected_mode = format!("{:?}", distance_mode);
            if modern_cache.metadata.distance_mode != expected_mode {
                println!("ℹ️  Cache was generated with mode '{}', using with mode '{}'", 
                         modern_cache.metadata.distance_mode, expected_mode);
                println!("   This is safe because cache contains all alignment statistics.");
            }
            
            // Convert modern format to internal format
            self.cache.clear();
            for (string_key, cache_value) in modern_cache.data {
                // Parse string key in format "locus:hash1:hash2"
                let parts: Vec<&str> = string_key.split(':').collect();
                if parts.len() == 3 {
                    let key = DistanceCacheKey {
                        locus: parts[0].to_string(),
                        crc1: parts[1].parse().unwrap_or(0),
                        crc2: parts[2].parse().unwrap_or(0),
                    };
                    let entry = CacheEntry {
                        snps: cache_value.snps,
                        indel_events: cache_value.indel_events,
                        indel_bases: cache_value.indel_bases,
                    };
                    self.cache.insert(key, entry);
                }
            }
            
            if let Some(note) = &modern_cache.metadata.user_note {
                println!("📝 Cache note: {}", note);
            }
            
            println!("✅ Modern cache loaded: {} loci, {} hasher, {} mode",
                     modern_cache.metadata.unique_loci,
                     modern_cache.metadata.hasher_type,
                     modern_cache.metadata.distance_mode);
            
        } else {
            // Fallback to legacy format (bincode) - for backward compatibility
            println!("🔄 Trying legacy cache format...");
            
            // Define legacy structure inline
            #[derive(Deserialize)]
            #[allow(clippy::type_complexity)]
            struct LegacyCache {
                data: HashMap<(String, u32, u32, u64), (usize, usize, usize)>,
                alignment_config: AlignmentConfig,
            }
            
            let legacy_cache: LegacyCache = bincode::deserialize(&decompressed)
                .map_err(|e| format!("Failed to deserialize cache (tried both modern and legacy formats): {}", e))?;
            
            // Check alignment configuration compatibility
            if legacy_cache.alignment_config != self.config {
                return Err(format!(
                    "Cache alignment config mismatch:\n  Cache: {:?}\n  Engine: {:?}",
                    legacy_cache.alignment_config, self.config
                ));
            }
            
            // Convert legacy format to internal format
            self.cache.clear();
            for ((locus, crc1, crc2, _config_hash), (snps, indel_events, indel_bases)) in legacy_cache.data {
                let key = DistanceCacheKey { locus, crc1, crc2 };
                let entry = CacheEntry { snps, indel_events, indel_bases };
                self.cache.insert(key, entry);
            }
            
            println!("⚠️  Loaded legacy cache format - consider regenerating with modern format");
        }
        
        let elapsed = start.elapsed();
        println!("✅ Cache loaded in {:.2}s ({} entries, {} KB)", 
                 elapsed.as_secs_f64(), 
                 self.cache.len(),
                 compressed.len() / 1024);
        
        // Reset the new entries flag since we just loaded from disk
        self.has_new_entries = false;
        
        Ok(())
    }
    
    
    /// Pre-compute all alignments for unique pairs in batch (cache-aware)
    pub fn precompute_alignments(&mut self, unique_pairs: &HashSet<(String, u32, u32)>, mode: DistanceMode) {
        let total_pairs = unique_pairs.len();
        
        // Filter out pairs already in cache (Strategy 1: Preventive filtering)
        let start_filter = Instant::now();
        let missing_pairs: Vec<_> = unique_pairs.iter()
            .filter(|(locus, crc1, crc2)| {
                // Pre-compute min/max to avoid multiple comparisons
                let (min_crc, max_crc) = if *crc1 <= *crc2 { (*crc1, *crc2) } else { (*crc2, *crc1) };
                let key = DistanceCacheKey {
                    locus: locus.clone(),
                    crc1: min_crc,
                    crc2: max_crc,
                };
                !self.cache.contains_key(&key)
            })
            .collect();
        
        let filter_elapsed = start_filter.elapsed();
        let cached_pairs = total_pairs - missing_pairs.len();
        let missing_count = missing_pairs.len();
        
        println!("🔍 Filtered unique pairs in {:.3}s:", filter_elapsed.as_secs_f64());
        println!("   📊 Total pairs needed: {}", total_pairs);
        println!("   ✅ Already in cache: {} ({:.1}%)", cached_pairs, 
                 (cached_pairs as f64 / total_pairs as f64) * 100.0);
        println!("   🔥 Missing pairs to compute: {} ({:.1}%)", missing_count,
                 (missing_count as f64 / total_pairs as f64) * 100.0);
        
        if missing_pairs.is_empty() {
            println!("🎯 All pairs already cached - no computation needed!");
            return;
        }
        
        // Setup progress bar for missing pairs only (update every 1% to reduce overhead)
        let pb = ProgressBar::new(missing_count as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({percent}%) {per_sec} ETA: {eta}")
                .unwrap()
                .progress_chars("#>-")
        );
        
        let start_compute = Instant::now();
        
        // Simple progress tracking - update every N completions
        let completed_count = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let update_frequency = 1000; // Update every 1000 completions
        
        // Compute alignments with periodic progress updates
        let results: Vec<_> = missing_pairs.into_par_iter()
            .filter_map(|(locus, crc1, crc2)| {
                let alignment_result = self.compute_single_alignment(locus, *crc1, *crc2, mode);
                
                // Increment and check if we should update progress
                let completed = completed_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                if completed % update_frequency == 0 {
                    pb.set_position(completed as u64);
                }
                
                // Only include pairs with successful alignment results
                alignment_result.map(|(snps, indel_events, indel_bases)| {
                    (locus, *crc1, *crc2, snps, indel_events, indel_bases)
                })
            })
            .collect();
        
        pb.finish_with_message("✅ Missing alignments computed!");
        
        // Store results in cache
        for (locus, crc1, crc2, snps, indel_events, indel_bases) in results {
            self.cache_distance(locus, crc1, crc2, snps, indel_events, indel_bases);
        }
        
        let compute_elapsed = start_compute.elapsed();
        
        println!("🚀 Cache-aware precompute completed:");
        println!("   ⚡ Computed {} new alignments in {:.2}s ({:.0} alignments/sec)", 
                 missing_count, compute_elapsed.as_secs_f64(), 
                 missing_count as f64 / compute_elapsed.as_secs_f64());
        println!("   📈 Total efficiency: {:.1}% time saved vs full recompute",
                 (cached_pairs as f64 / total_pairs as f64) * 100.0);
        println!("   💾 Cache now contains {} entries", self.cache.len());
    }
    
    /// Compute a single alignment (used in batch processing)  
    /// Returns None if sequences are not available (cache miss should apply Hamming fallback)
    fn compute_single_alignment(&self, locus: &str, crc1: u32, crc2: u32, mode: DistanceMode) -> Option<(usize, usize, usize)> {
        if crc1 == crc2 {
            return Some((0, 0, 0)); // Identical alleles
        }
        
        // For Hamming mode, different CRCs = distance 1 (no sequence analysis needed)
        if matches!(mode, DistanceMode::Hamming) {
            return Some((1, 0, 0)); // Different CRC = 1 Hamming distance unit
        }
        
        // Try to get sequences and align (for non-Hamming modes)
        if let Some(ref seq_db) = self.sequence_db {
            if let (Some(seq1), Some(seq2)) = (seq_db.get_sequence(locus, crc1), seq_db.get_sequence(locus, crc2)) {
                
                // Create traceback-enabled aligner for other modes
                let alphabet = b"ACGT";
                let matrix = match Matrix::create(alphabet, self.config.match_score, self.config.mismatch_penalty) {
                    Ok(m) => m,
                    Err(_) => return None, // Matrix creation failed, no sequences to align
                };
                
                let trace_aligner = Aligner::new()
                    .matrix(matrix)
                    .gap_open(self.config.gap_open)
                    .gap_extend(self.config.gap_extend)
                    .global()
                    .use_trace()  // Enable traceback
                    .build();
                
                match trace_aligner.align(Some(&seq1.sequence), &seq2.sequence) {
                    Ok(result) => {
                        // Get traceback strings with gaps
                        match result.get_traceback_strings(&seq1.sequence, &seq2.sequence) {
                            Ok(traceback) => {
                                // Use the proper alignment analysis with gaps
                                let (snps, indel_events, indel_bases) = compute_alignment_stats(&traceback.query, &traceback.reference);
                                
                                // TODO: Save alignment details if requested
                                // self.add_alignment_detail(locus, crc1, crc2, &seq1.sequence, &seq2.sequence,
                                //                          &traceback.query, &traceback.reference, 
                                //                          snps, indel_events, indel_bases, result.get_score() as f32);
                                
                                return Some((snps, indel_events, indel_bases));
                            }
                            Err(_) => {
                                // Traceback failed, use simple approach
                                let (snps, indel_events, indel_bases) = self.analyze_sequences(&seq1.sequence, &seq2.sequence);
                                return Some((snps, indel_events, indel_bases));
                            }
                        }
                    }
                    Err(_) => {
                        // Alignment failed, use simple comparison
                        let (snps, indel_events, indel_bases) = self.analyze_sequences(&seq1.sequence, &seq2.sequence);
                        return Some((snps, indel_events, indel_bases));
                    }
                }
            }
        }
        
        // No sequences available for alignment
        None
    }
    
    /// Analyze two sequences to count SNPs, indel events, and indel bases
    /// This implements a simplified alignment analysis that handles gaps (fallback only)
    fn analyze_sequences(&self, seq1: &[u8], seq2: &[u8]) -> (usize, usize, usize) {
        // For identical length sequences, do direct comparison
        if seq1.len() == seq2.len() {
            let mut snps = 0;
            for i in 0..seq1.len() {
                if seq1[i] != seq2[i] {
                    snps += 1;
                }
            }
            return (snps, 0, 0);
        }
        
        // For different lengths, implement a simple gap-aware analysis
        // This is a simplified approach that assumes optimal alignment
        let len_diff = seq1.len().abs_diff(seq2.len());
        let min_len = seq1.len().min(seq2.len());
        
        // Count SNPs in the overlapping region
        let mut snps = 0;
        for i in 0..min_len {
            if seq1[i] != seq2[i] {
                snps += 1;
            }
        }
        
        // Length difference represents indel events and bases
        let indel_events = if len_diff > 0 { 1 } else { 0 };
        let indel_bases = len_diff;
        
        (snps, indel_events, indel_bases)
    }
    
    /// Check if cache file exists and has sequence lengths enrichment
    pub fn cache_has_lengths(&self, cache_path: &str) -> Result<bool, String> {
        if !std::path::Path::new(cache_path).exists() {
            return Ok(false); // Cache doesn't exist = no lengths
        }
        
        // Try to load cache and check if it has sequence lengths
        let compressed = std::fs::read(cache_path)
            .map_err(|e| format!("Failed to read cache file: {}", e))?;
        
        // Decompress with LZ4
        let decompressed = lz4_flex::decompress_size_prepended(&compressed)
            .map_err(|e| format!("Failed to decompress cache: {}", e))?;
        
        // Load as ModernCache
        let modern_cache: ModernCache = serde_json::from_slice(&decompressed)
            .map_err(|e| format!("Failed to deserialize cache: {}", e))?;
        
        // Check if any entries have sequence lengths
        for (_, value) in modern_cache.data.iter().take(5) { // Check first 5 entries
            if value.seq1_length.is_some() || value.seq2_length.is_some() {
                return Ok(true);
            }
        }
        
        Ok(false)
    }

    /// Enrich existing cache with nucleotide sequence lengths from schema (with input/output paths)
    pub fn enrich_cache_with_lengths_from_input(&mut self, schema_path: &str, input_cache_path: &str, output_cache_path: &str) -> Result<usize, String> {
        // Load sequence lengths from schema with CRC mapping
        println!("🔍 Loading schema lengths with CRC mapping from {}...", schema_path);
        let (schema_lengths, crc_to_length) = self.load_schema_with_crc_mapping(schema_path)?;
        println!("📊 Loaded {} loci from schema with {} CRC mappings", schema_lengths.len(), crc_to_length.len());
        
        // Read from input cache file
        if !std::path::Path::new(input_cache_path).exists() {
            return Err("Input cache file does not exist".to_string());
        }
        
        println!("📂 Loading cache from {}...", input_cache_path);
        let compressed = std::fs::read(input_cache_path)
            .map_err(|e| format!("Failed to read cache file: {}", e))?;
        
        // Decompress with LZ4
        let decompressed = lz4_flex::decompress_size_prepended(&compressed)
            .map_err(|e| format!("Failed to decompress cache: {}", e))?;
        
        // Load as ModernCache
        let mut modern_cache: ModernCache = serde_json::from_slice(&decompressed)
            .map_err(|e| format!("Failed to deserialize cache: {}", e))?;
        
        println!("✅ Loaded cache with {} entries", modern_cache.data.len());
        
        // Enrich cache entries with sequence lengths
        println!("🔍 Enriching cache entries with sequence lengths...");
        let mut enriched_count = 0;
        let mut missing_entries = Vec::new();
        
        for (key, value) in &mut modern_cache.data {
            // Parse key format: "locus:crc1:crc2"
            let parts: Vec<&str> = key.split(':').collect();
            if parts.len() >= 3 {
                let locus = parts[0];
                let crc1_str = parts[1];
                let crc2_str = parts[2];
                
                let mut found_any = false;
                
                // Parse CRCs as u32
                if let (Ok(crc1), Ok(crc2)) = (crc1_str.parse::<u32>(), crc2_str.parse::<u32>()) {
                    // Look up lengths using CRC mapping
                    if let Some(&len1) = crc_to_length.get(&crc1) {
                        value.seq1_length = Some(len1);
                        found_any = true;
                    }
                    
                    if let Some(&len2) = crc_to_length.get(&crc2) {
                        value.seq2_length = Some(len2);
                        found_any = true;
                    }
                    
                    if found_any {
                        enriched_count += 1;
                    } else {
                        missing_entries.push(format!("{}:{}:{}", locus, crc1, crc2));
                    }
                } else {
                    missing_entries.push(format!("{}:{}:{} (invalid CRC format)", locus, crc1_str, crc2_str));
                }
            }
        }
        
        println!("✅ Enriched {} out of {} entries", enriched_count, modern_cache.data.len());
        
        if !missing_entries.is_empty() && missing_entries.len() <= 10 {
            println!("⚠️  Warning: {} entries with missing alleles:", missing_entries.len());
            for entry in &missing_entries {
                println!("   - {}", entry);
            }
        } else if !missing_entries.is_empty() {
            println!("⚠️  Warning: {} entries with missing alleles (showing first 10):", missing_entries.len());
            for entry in missing_entries.iter().take(10) {
                println!("   - {}", entry);
            }
        }
        
        // Update metadata
        modern_cache.metadata.last_modified = chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC").to_string();
        if let Some(ref mut note) = modern_cache.metadata.user_note {
            note.push_str(" [Enriched with sequence lengths]");
        } else {
            modern_cache.metadata.user_note = Some("Enriched with sequence lengths".to_string());
        }
        
        // Serialize and save
        println!("💾 Saving enriched cache to {}...", output_cache_path);
        let serialized = serde_json::to_vec(&modern_cache)
            .map_err(|e| format!("Failed to serialize cache: {}", e))?;
        
        // Compress with lz4_flex (matches cgDist format)
        let final_data = lz4_flex::compress_prepend_size(&serialized);
        
        std::fs::write(output_cache_path, final_data)
            .map_err(|e| format!("Failed to write cache file: {}", e))?;
        
        println!("✅ Enriched cache saved ({:.1}% success rate)", 
                (enriched_count as f64 / modern_cache.data.len() as f64) * 100.0);
        
        Ok(enriched_count)
    }
    
    /// Load schema with both lengths and CRC mapping for enrichment
    #[allow(clippy::type_complexity)]
    fn load_schema_with_crc_mapping(&self, schema_path: &str) -> Result<(HashMap<String, HashMap<String, usize>>, HashMap<u32, usize>), String> {
        use std::fs;
        use std::path::Path;
        
        let schema_dir = Path::new(schema_path);
        let mut all_lengths = HashMap::new();
        let mut crc_to_length = HashMap::new();
        
        // Get hasher for CRC calculation
        let registry = HasherRegistry::new();
        let hasher = registry.get_hasher(&self.hasher_type)
            .ok_or_else(|| format!("Unknown hasher type: {}", self.hasher_type))?;
        
        let entries = fs::read_dir(schema_dir)
            .map_err(|e| format!("Failed to read schema directory: {}", e))?;
        
        for entry in entries {
            let entry = entry.map_err(|e| format!("Failed to read directory entry: {}", e))?;
            let path = entry.path();
            
            if path.extension().and_then(|s| s.to_str()) == Some("fasta") {
                if let Some(filename) = path.file_stem().and_then(|s| s.to_str()) {
                    let locus_name = filename.to_string();
                    match self.load_fasta_with_crc_mapping(&path, hasher) {
                        Ok((lengths, crc_map)) => {
                            all_lengths.insert(locus_name, lengths);
                            crc_to_length.extend(crc_map);
                        }
                        Err(e) => {
                            eprintln!("⚠️  Warning: Failed to load {}: {}", filename, e);
                        }
                    }
                }
            }
        }
        
        Ok((all_lengths, crc_to_length))
    }
    
    /// Load FASTA file with both lengths and CRC mapping
    #[allow(clippy::type_complexity)]
    fn load_fasta_with_crc_mapping(&self, fasta_path: &Path, hasher: &dyn AlleleHasher) -> Result<(HashMap<String, usize>, HashMap<u32, usize>), String> {
        use std::fs;
        
        let content = fs::read_to_string(fasta_path)
            .map_err(|e| format!("Failed to read FASTA file: {}", e))?;
        
        let mut lengths = HashMap::new();
        let mut crc_to_length = HashMap::new();
        let mut current_id = String::new();
        let mut current_sequence = String::new();
        
        for line in content.lines() {
            if let Some(stripped) = line.strip_prefix('>') {
                // Save previous sequence if exists
                if !current_id.is_empty() && !current_sequence.is_empty() {
                    let seq_len = current_sequence.len();
                    lengths.insert(current_id.clone(), seq_len);

                    // Calculate hash for this sequence
                    let hash = hasher.hash_sequence(&current_sequence);
                    if let Some(crc) = hash.as_crc32() {
                        crc_to_length.insert(crc, seq_len);
                    }
                }

                // Parse new sequence ID (e.g., >INNUENDO_cgMLST-00031717_1)
                current_id = stripped.split_whitespace().next()
                    .unwrap_or("")
                    .to_string();

                // Extract just the allele number
                if let Some(underscore_pos) = current_id.rfind('_') {
                    current_id = current_id[underscore_pos + 1..].to_string();
                }

                current_sequence.clear();
            } else {
                // Add to current sequence
                current_sequence.push_str(line.trim());
            }
        }
        
        // Save last sequence
        if !current_id.is_empty() && !current_sequence.is_empty() {
            let seq_len = current_sequence.len();
            lengths.insert(current_id, seq_len);
            
            // Calculate hash for last sequence
            let hash = hasher.hash_sequence(&current_sequence);
            if let Some(crc) = hash.as_crc32() {
                crc_to_length.insert(crc, seq_len);
            }
        }
        
        Ok((lengths, crc_to_length))
    }
    
    /// Set path for saving detailed alignments
    pub fn set_save_alignments(&mut self, path: String) {
        self.save_alignments_path = Some(path);
        // Initialize with TSV header
        self.alignment_details.clear();
        self.alignment_details.push("locus\thash1\thash2\tseq1\tseq2\taligned_seq1\taligned_seq2\tsnps\tindel_events\tindel_bases\talignment_score".to_string());
    }
    
    /// Save alignment details to TSV file
    pub fn save_alignments(&self) -> Result<(), String> {
        if let Some(path) = &self.save_alignments_path {
            if !self.alignment_details.is_empty() {
                let content = self.alignment_details.join("\n") + "\n";
                std::fs::write(path, content)
                    .map_err(|e| format!("Failed to write alignments file: {}", e))?;
                println!("💾 Saved {} alignment details to: {}", self.alignment_details.len() - 1, path);
            }
        }
        Ok(())
    }
    
    // TODO: Implement add_alignment_detail for --save-alignments functionality
    // Currently disabled due to threading complexity with Rayon
}

/// Calculate distance between two samples
pub fn calculate_sample_distance(
    sample1: &AllelicProfile,
    sample2: &AllelicProfile,
    loci_names: &[String],
    engine: &DistanceEngine,
    mode: DistanceMode,
    min_loci: usize,
    no_hamming_fallback: bool,
) -> Option<usize> {
    let mut total_distance = 0;
    let mut shared_loci = 0;
    
    for locus in loci_names {
        let crc1 = sample1.loci_hashes.get(locus)
            .and_then(|h| h.as_crc32())
            .unwrap_or(u32::MAX);
        let crc2 = sample2.loci_hashes.get(locus)
            .and_then(|h| h.as_crc32())
            .unwrap_or(u32::MAX);
        
        if crc1 != u32::MAX && crc2 != u32::MAX {
            shared_loci += 1;
        }
        
        total_distance += engine.get_distance(locus, crc1, crc2, mode, no_hamming_fallback);
    }
    
    if shared_loci >= min_loci {
        Some(total_distance)
    } else {
        None
    }
}

/// Calculate full distance matrix
pub fn calculate_distance_matrix(
    samples: &[AllelicProfile],
    loci_names: &[String],
    engine: &DistanceEngine,
    mode: DistanceMode,
    min_loci: usize,
    no_hamming_fallback: bool,
) -> Vec<Vec<Option<usize>>> {
    let n_samples = samples.len();
    let mut matrix = vec![vec![None; n_samples]; n_samples];
    
    // Fill diagonal with zeros
    for (i, row) in matrix.iter_mut().enumerate().take(n_samples) {
        row[i] = Some(0);
    }
    
    // Calculate upper triangle in parallel with progress bar
    let start = Instant::now();
    let total_comparisons = n_samples * (n_samples - 1) / 2;
    println!("🔄 Computing distance matrix ({} × {} = {} comparisons)...", 
             n_samples, n_samples, total_comparisons);
    
    use indicatif::{ProgressBar, ProgressStyle};
    let pb = ProgressBar::new(total_comparisons as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({percent}%) {per_sec} ETA: {eta}")
            .unwrap()
            .progress_chars("#>-")
    );
    
    // Progress tracking with reduced contention
    let update_interval = std::cmp::max(1, total_comparisons / 100); // Update every 1%
    let progress_counter = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
    
    let upper_triangle: Vec<_> = (0..n_samples).into_par_iter()
        .flat_map(|i| {
            let progress_clone = progress_counter.clone();
            let pb_clone = pb.clone();
            (i + 1..n_samples).into_par_iter().map(move |j| {
                let distance = calculate_sample_distance(
                    &samples[i],
                    &samples[j],
                    loci_names,
                    engine,
                    mode,
                    min_loci,
                    no_hamming_fallback,
                );
                
                // Update progress periodically
                let count = progress_clone.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                if count % update_interval == 0 {
                    pb_clone.set_position(count as u64);
                }
                
                (i, j, distance)
            })
        })
        .collect();
        
    pb.finish_with_message("✅ Distance matrix computation completed!");
    
    // Fill matrix symmetrically
    for (i, j, distance) in upper_triangle {
        matrix[i][j] = distance;
        matrix[j][i] = distance;
    }
    
    let elapsed = start.elapsed();
    println!("✅ Distance matrix computed in {:.2}s", elapsed.as_secs_f64());
    
    matrix
}