// distance.rs - Core distance calculation engine

use crate::core::alignment::{
    cigar_from_aligned, compute_alignment_stats, AlignmentConfig, DistanceMode,
};
use crate::core::banded::{align_certified, align_certified_with_strings, Scoring};
use crate::data::{AllelicProfile, SequenceDatabase};
use crate::hashers::{AlleleHasher, HasherRegistry};
use chrono;
use indicatif::{ProgressBar, ProgressStyle};
use parasail_rs::{Aligner, Matrix};
use rayon::iter::ParallelIterator;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::time::Instant;

/// Internal cache entry storing all alignment statistics
#[derive(Debug, Clone, Copy)]
struct CacheEntry {
    snps: usize,
    indel_events: usize,
    indel_bases: usize,
    /// Nucleotide lengths of the two alleles (key order: smaller CRC first),
    /// when known. Enables per-locus mutation-density / recombination signals.
    lens: (Option<u32>, Option<u32>),
}

impl CacheEntry {
    /// Mean allele length, as cgdist has always derived it from the two
    /// stored lengths.
    fn mean_len(&self) -> Option<u32> {
        match self.lens {
            (Some(a), Some(b)) => Some((a + b) / 2),
            (Some(a), None) => Some(a),
            (None, Some(b)) => Some(b),
            (None, None) => None,
        }
    }
}

/// A TSV output written row by row as alignments are produced
/// (--save-alignments, --save-cigar). The file is created, with its header,
/// on the first row or at `finish`, so it exists even when no row is written.
struct RowWriter {
    path: String,
    header: &'static str,
    out: Option<std::io::BufWriter<std::fs::File>>,
    rows: usize,
    error: Option<String>,
}

impl RowWriter {
    fn new(path: String, header: &'static str) -> Self {
        Self {
            path,
            header,
            out: None,
            rows: 0,
            error: None,
        }
    }

    fn open(&mut self) {
        use std::io::Write;
        match std::fs::File::create(&self.path) {
            Ok(f) => {
                let mut w = std::io::BufWriter::with_capacity(1 << 20, f);
                match writeln!(w, "{}", self.header) {
                    Ok(()) => self.out = Some(w),
                    Err(e) => self.error = Some(format!("Failed to write {}: {e}", self.path)),
                }
            }
            Err(e) => self.error = Some(format!("Failed to write {}: {e}", self.path)),
        }
    }

    /// Append a row; a write error is kept, reported by `finish`, and later
    /// rows are skipped.
    fn write(&mut self, row: &str) {
        use std::io::Write;
        if self.error.is_some() {
            return;
        }
        if self.out.is_none() {
            self.open();
        }
        if let Some(w) = self.out.as_mut() {
            match writeln!(w, "{row}") {
                Ok(()) => self.rows += 1,
                Err(e) => self.error = Some(format!("Failed to write {}: {e}", self.path)),
            }
        }
    }

    /// Flush and close; returns the number of rows written.
    fn finish(&mut self) -> Result<usize, String> {
        use std::io::Write;
        if self.out.is_none() && self.error.is_none() {
            self.open();
        }
        if let Some(mut w) = self.out.take() {
            if let Err(e) = w.flush() {
                self.error
                    .get_or_insert(format!("Failed to write {}: {e}", self.path));
            }
        }
        match self.error.take() {
            Some(e) => Err(e),
            None => Ok(self.rows),
        }
    }
}

/// Per-pair output rows produced alongside the statistics.
#[derive(Default)]
struct PairDetail {
    /// --save-alignments row
    full: Option<String>,
    /// --save-cigar row
    cigar: Option<String>,
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
    // Fraction of new alignments re-checked against parasail's original kernel
    verify_fraction: f64,
    // --save-alignments / --save-cigar outputs, written as rows are produced
    alignments_out: Option<RowWriter>,
    cigar_out: Option<RowWriter>,
}

/// Modern cache structure supporting any hasher type
#[derive(Debug, Serialize, Deserialize)]
pub struct ModernCache {
    pub data: HashMap<String, CacheValue>, // Use string keys for JSON compatibility
    pub metadata: CacheMetadata,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CacheKey {
    pub locus: String,
    pub hash1: String, // String to support any hasher (CRC32, SHA256, MD5, etc.)
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
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
struct DistanceCacheKey {
    locus: String,
    crc1: u32,
    crc2: u32,
}

/// Borrowed view of a cache key, so lookups by (&str, u32, u32) need no
/// String allocation. Hash and equality are defined on the same parts for
/// the owned key and every view, as `HashMap`'s `Borrow` contract requires.
trait KeyView {
    fn parts(&self) -> (&str, u32, u32);
}

impl KeyView for DistanceCacheKey {
    fn parts(&self) -> (&str, u32, u32) {
        (&self.locus, self.crc1, self.crc2)
    }
}

impl KeyView for (&str, u32, u32) {
    fn parts(&self) -> (&str, u32, u32) {
        (self.0, self.1, self.2)
    }
}

impl<'a> std::borrow::Borrow<dyn KeyView + 'a> for DistanceCacheKey {
    fn borrow(&self) -> &(dyn KeyView + 'a) {
        self
    }
}

impl std::hash::Hash for dyn KeyView + '_ {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.parts().hash(state)
    }
}

impl PartialEq for dyn KeyView + '_ {
    fn eq(&self, other: &Self) -> bool {
        self.parts() == other.parts()
    }
}

impl Eq for dyn KeyView + '_ {}

impl std::hash::Hash for DistanceCacheKey {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.parts().hash(state)
    }
}

thread_local! {
    /// Per-thread parasail aligners keyed by (match, mismatch, gap_open,
    /// gap_extend, solution width). Building an aligner allocates a scoring
    /// matrix and looks up the kernel, so it is done once per thread.
    #[allow(clippy::type_complexity)]
    static ALIGNERS: std::cell::RefCell<Vec<((i32, i32, i32, i32, i32), Aligner)>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

/// Global (Needleman-Wunsch) alignment with traceback.
///
/// Uses parasail's scan kernel at 16-bit precision and redoes the pair at
/// 32-bit (then 64-bit) when parasail reports saturation. On 739,554 real allele pairs
/// (L. monocytogenes, S. enterica) this reproduces the previous striped
/// saturating kernel exactly (SNPs, InDel events, InDel bases) and is ~7x
/// faster: with traceback, parasail's striped kernel is much slower than scan.
/// Returns None only if the scoring matrix cannot be created.
fn align_global_trace(
    config: &AlignmentConfig,
    query: &[u8],
    reference: &[u8],
) -> Option<Result<parasail_rs::AlignResult, parasail_rs::AlignError>> {
    let run = |width: i32| -> Option<Result<parasail_rs::AlignResult, parasail_rs::AlignError>> {
        let key = (
            config.match_score,
            config.mismatch_penalty,
            config.gap_open,
            config.gap_extend,
            width,
        );
        ALIGNERS.with(|cell| {
            let mut aligners = cell.borrow_mut();
            let idx = match aligners.iter().position(|(k, _)| *k == key) {
                Some(i) => i,
                None => {
                    let matrix =
                        Matrix::create(b"ACGT", config.match_score, config.mismatch_penalty)
                            .ok()?;
                    let aligner = Aligner::new()
                        .matrix(matrix)
                        .gap_open(config.gap_open)
                        .gap_extend(config.gap_extend)
                        .global()
                        .use_trace()
                        .scan()
                        .solution_width(width)
                        .build();
                    aligners.push((key, aligner));
                    aligners.len() - 1
                }
            };
            Some(aligners[idx].1.align(Some(query), reference))
        })
    };
    // Escalate precision only for the pairs that need it. A saturated result
    // is never returned: parasail flags saturation conservatively (before an
    // actual overflow), and 64-bit cannot saturate for any real allele pair.
    for width in [16, 32] {
        match run(width)? {
            Ok(res) if res.is_saturated() => continue,
            other => return Some(other),
        }
    }
    match run(64)? {
        Ok(res) if res.is_saturated() => {
            panic!(
                "alignment saturated even at 64-bit precision (sequence lengths {} and {})",
                query.len(),
                reference.len()
            )
        }
        other => Some(other),
    }
}

/// A pairwise alignment with its gapped strings.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairAlignment {
    pub snps: usize,
    pub indel_events: usize,
    pub indel_bases: usize,
    pub score: i32,
    /// query (first allele) with '-' for gaps
    pub query: String,
    /// reference (second allele) with '-' for gaps
    pub reference: String,
}

/// Align two alleles exactly as cgdist does when it writes --save-alignments
/// rows: certified banded alignment when possible, else parasail; statistics
/// are computed from the gapped strings. Returns None when the pair cannot
/// be aligned (parasail refused it).
pub fn align_pair_with_strings(
    config: &AlignmentConfig,
    query: &[u8],
    reference: &[u8],
) -> Option<PairAlignment> {
    if !query.contains(&0) && !reference.contains(&0) && query.is_ascii() && reference.is_ascii() {
        let scoring = Scoring {
            match_score: config.match_score,
            mismatch: config.mismatch_penalty,
            gap_open: config.gap_open,
            gap_extend: config.gap_extend,
        };
        if let Some((b, st)) = align_certified_with_strings(query, reference, &scoring, 0.5) {
            let q = String::from_utf8_lossy(&st.query).into_owned();
            let r = String::from_utf8_lossy(&st.reference).into_owned();
            let (snps, indel_events, indel_bases) = compute_alignment_stats(&q, &r);
            return Some(PairAlignment {
                snps,
                indel_events,
                indel_bases,
                score: b.score,
                query: q,
                reference: r,
            });
        }
    }
    let res = align_global_trace(config, query, reference)?.ok()?;
    let tb = res.get_traceback_strings(query, reference).ok()?;
    let (snps, indel_events, indel_bases) = compute_alignment_stats(&tb.query, &tb.reference);
    Some(PairAlignment {
        snps,
        indel_events,
        indel_bases,
        score: res.get_score(),
        query: tb.query,
        reference: tb.reference,
    })
}

/// Deterministic choice of the pairs re-checked by --verify-alignments:
/// depends only on the two allele hashes, not on threads or run order.
fn verify_selected(crc1: u32, crc2: u32, fraction: f64) -> bool {
    if fraction >= 1.0 {
        return true;
    }
    let (a, b) = (crc1.min(crc2) as u64, crc1.max(crc2) as u64);
    let mut h = (a << 32 | b).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    h ^= h >> 29;
    h = h.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    h ^= h >> 32;
    ((h >> 11) as f64 / (1u64 << 53) as f64) < fraction
}

/// Statistics from parasail's original production kernel (striped,
/// saturating 8/16/32-bit), the reference for --verify-alignments.
fn reference_alignment_stats(
    config: &AlignmentConfig,
    query: &[u8],
    reference: &[u8],
) -> Option<(usize, usize, usize)> {
    let matrix = Matrix::create(b"ACGT", config.match_score, config.mismatch_penalty).ok()?;
    let aligner = Aligner::new()
        .matrix(matrix)
        .gap_open(config.gap_open)
        .gap_extend(config.gap_extend)
        .global()
        .use_trace()
        .build();
    let res = aligner.align(Some(query), reference).ok()?;
    let tb = res.get_traceback_strings(query, reference).ok()?;
    Some(compute_alignment_stats(&tb.query, &tb.reference))
}

/// True for cache files written by a cgdist older than 0.1.4, the first
/// release that no longer stores placeholder statistics in Hamming mode.
fn written_before_hamming_fix(version: &str) -> bool {
    let mut it = version.split('.').map(|p| {
        p.chars()
            .take_while(|c| c.is_ascii_digit())
            .collect::<String>()
            .parse::<u64>()
            .unwrap_or(0)
    });
    let v = (
        it.next().unwrap_or(0),
        it.next().unwrap_or(0),
        it.next().unwrap_or(0),
    );
    v < (0, 1, 4)
}

/// Alignment results are interchangeable iff the scoring parameters match;
/// the free-text description is irrelevant.
fn same_alignment_params(a: &AlignmentConfig, b: &AlignmentConfig) -> bool {
    a.match_score == b.match_score
        && a.mismatch_penalty == b.mismatch_penalty
        && a.gap_open == b.gap_open
        && a.gap_extend == b.gap_extend
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
            verify_fraction: 0.0,
            alignments_out: None,
            cigar_out: None,
        }
    }

    pub fn with_sequences(
        config: AlignmentConfig,
        sequence_db: SequenceDatabase,
        hasher_type: String,
    ) -> Self {
        Self {
            cache: HashMap::new(),
            config,
            sequence_db: Some(sequence_db),
            hasher_type,
            cache_note: None,
            has_new_entries: false,
            save_alignments_path: None,
            verify_fraction: 0.0,
            alignments_out: None,
            cigar_out: None,
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

        // Hamming distance never needs the alignment cache: different CRCs = 1.
        // (Answering from the cache here would make the result depend on
        // whatever happens to be cached.)
        if self.hasher_type == "hamming" || mode == DistanceMode::Hamming {
            return 1;
        }

        // Optimized key lookup - temporarily create key for lookup only
        let (min_crc, max_crc) = if crc1 <= crc2 {
            (crc1, crc2)
        } else {
            (crc2, crc1)
        };

        if let Some(&entry) = self.cache.get(&(locus, min_crc, max_crc) as &dyn KeyView) {
            let distance = match mode {
                DistanceMode::SnpsOnly => entry.snps,
                DistanceMode::SnpsAndIndelEvents => entry.snps + entry.indel_events,
                DistanceMode::SnpsAndIndelBases => entry.snps + entry.indel_bases,
                DistanceMode::Hamming => 1, // For hamming mode, different CRCs = 1
            };

            // Apply Hamming fallback ONLY for SNPs mode: if alignment found 0 differences
            // but CRCs are different, return 1 to maintain consistency (different CRCs >= 1 difference)
            if distance == 0
                && crc1 != crc2
                && !no_hamming_fallback
                && mode == DistanceMode::SnpsOnly
            {
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

    /// Whether per-locus mutation-density (recombination) signals are available,
    /// i.e. the loaded cache carries sequence lengths (enriched cache) and the
    /// hasher is sequence-based. Used to decide whether to include recombination
    /// in the dashboard.
    pub fn has_recomb_data(&self) -> bool {
        self.hasher_type != "hamming" && self.cache.values().any(|e| e.mean_len().is_some())
    }

    /// Per-locus recombination signal for an allele pair: `Some(true)` if the
    /// mutation density (SNPs + InDel bases) / allele length exceeds
    /// `thresh_frac` (e.g. 0.03 = 3%), `Some(false)` if below, `None` when the
    /// pair is missing/identical or length data is unavailable.
    pub fn locus_is_recombinant(
        &self,
        locus: &str,
        crc1: u32,
        crc2: u32,
        thresh_frac: f64,
    ) -> Option<bool> {
        if crc1 == u32::MAX || crc2 == u32::MAX || crc1 == crc2 || self.hasher_type == "hamming" {
            return None;
        }
        let (min_crc, max_crc) = if crc1 <= crc2 {
            (crc1, crc2)
        } else {
            (crc2, crc1)
        };
        let entry = self.cache.get(&(locus, min_crc, max_crc) as &dyn KeyView)?;
        let len = entry.mean_len()?;
        if len == 0 {
            return None;
        }
        let density = (entry.snps + entry.indel_bases) as f64 / len as f64;
        Some(density > thresh_frac)
    }

    /// Add distance to cache (optimized) - stores all alignment statistics
    pub fn cache_distance(
        &mut self,
        locus: &str,
        crc1: u32,
        crc2: u32,
        snps: usize,
        indel_events: usize,
        indel_bases: usize,
    ) {
        let (min_crc, max_crc) = if crc1 <= crc2 {
            (crc1, crc2)
        } else {
            (crc2, crc1)
        };
        let key = DistanceCacheKey {
            locus: locus.to_string(), // Only allocate when inserting into cache
            crc1: min_crc,
            crc2: max_crc,
        };
        let entry = CacheEntry {
            snps,
            indel_events,
            indel_bases,
            lens: (None, None), // set by precompute_alignments when sequences are known
        };
        self.cache.insert(key, entry);
        self.has_new_entries = true; // Mark that cache has new entries
    }

    /// True when every cache entry carries both allele lengths, i.e. a
    /// separate enrichment pass over the schema would change nothing.
    pub fn all_entries_have_lengths(&self) -> bool {
        self.cache
            .values()
            .all(|e| e.lens.0.is_some() && e.lens.1.is_some())
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
    pub fn save_cache(
        &mut self,
        cache_path: &str,
        distance_mode: DistanceMode,
    ) -> Result<(), String> {
        println!("💾 Saving cache to {cache_path}...");
        let start = Instant::now();

        // Convert internal cache to modern format
        let mut data = HashMap::new();
        let now = chrono::Utc::now()
            .format("%Y-%m-%d %H:%M:%S UTC")
            .to_string();

        // Count unique loci efficiently without cloning
        let unique_loci: std::collections::HashSet<&str> =
            self.cache.keys().map(|k| k.locus.as_str()).collect();

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
                // Keep the allele lengths already known (computed or enriched)
                seq1_length: entry.lens.0.map(|l| l as usize),
                seq2_length: entry.lens.1.map(|l| l as usize),
            };

            data.insert(string_key, cache_value);
        }

        // Same note the enrichment pass adds, so a cache whose lengths were
        // recorded at alignment time is labelled like an enriched one.
        let mut user_note = self.cache_note.clone();
        if !self.cache.is_empty() && self.all_entries_have_lengths() {
            match user_note {
                Some(ref mut note) if !note.contains("Enriched with sequence lengths") => {
                    note.push_str(" [Enriched with sequence lengths]")
                }
                None => user_note = Some("Enriched with sequence lengths".to_string()),
                _ => {}
            }
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
            user_note,
            total_entries: self.cache.len(),
            unique_loci: unique_loci.len(),
            format_version: 2, // Version 2 = modern format
        };

        let modern_cache = ModernCache { data, metadata };

        // Serialize cache with serde_json (more readable than bincode)
        let cache_data = serde_json::to_vec(&modern_cache)
            .map_err(|e| format!("Failed to serialize cache: {e}"))?;

        // Compress with LZ4
        let compressed = lz4_flex::compress_prepend_size(&cache_data);

        // Write to file
        std::fs::write(cache_path, &compressed)
            .map_err(|e| format!("Failed to write cache file: {e}"))?;

        let elapsed = start.elapsed();
        println!(
            "✅ Cache saved in {:.2}s ({} entries, {} KB)",
            elapsed.as_secs_f64(),
            self.cache.len(),
            compressed.len() / 1024
        );

        if let Some(note) = &self.cache_note {
            println!("📝 User note: {note}");
        }

        // Reset the new entries flag since we just saved to disk
        self.has_new_entries = false;

        Ok(())
    }

    /// Quick compatibility check without loading full cache  
    pub fn check_cache_compatibility(
        &self,
        cache_path: &str,
        _distance_mode: DistanceMode,
    ) -> Result<(), String> {
        use std::fs::File;
        use std::io::Read;

        // Read only first 32KB which should contain metadata
        let mut file =
            File::open(cache_path).map_err(|e| format!("Failed to open cache file: {e}"))?;

        let mut buffer = vec![0u8; 32768]; // 32KB should be enough for headers
        let bytes_read = file
            .read(&mut buffer)
            .map_err(|e| format!("Failed to read cache header: {e}"))?;
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
            return Err(
                "Cache alignment config mismatch: cache uses gap_open=5, current uses gap_open=8"
                    .to_string(),
            );
        }
        if buffer_str.contains("gap_open\":8") && self.config.gap_open == 5 {
            return Err(
                "Cache alignment config mismatch: cache uses gap_open=8, current uses gap_open=5"
                    .to_string(),
            );
        }

        Ok(())
    }

    /// Load cache from LZ4 compressed file (supports both modern and legacy formats)
    pub fn load_cache(
        &mut self,
        cache_path: &str,
        distance_mode: DistanceMode,
    ) -> Result<(), String> {
        println!("📂 Loading cache from {cache_path}...");
        let start = Instant::now();

        // Read compressed file
        let compressed =
            std::fs::read(cache_path).map_err(|e| format!("Failed to read cache file: {e}"))?;

        // Decompress with LZ4
        let decompressed = lz4_flex::decompress_size_prepended(&compressed)
            .map_err(|e| format!("Failed to decompress cache: {e}"))?;

        // Try modern format first (JSON)
        if let Ok(modern_cache) = serde_json::from_slice::<ModernCache>(&decompressed) {
            println!(
                "🆕 Loading modern cache format (v{})",
                modern_cache.metadata.format_version
            );

            // Check alignment configuration compatibility
            if !same_alignment_params(&modern_cache.metadata.alignment_config, &self.config) {
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

            // cgdist <= 0.1.3 stored placeholder statistics (1 SNP, 0 InDels)
            // for every pair it saw in Hamming mode. A cache last written in
            // Hamming mode may therefore hold values that are not alignments.
            if modern_cache.metadata.distance_mode == "hamming"
                && written_before_hamming_fix(&modern_cache.metadata.version)
            {
                eprintln!(
                    "⚠️  WARNING: this cache was last written by a Hamming-mode run. cgdist <= 0.1.3 \
                     stored placeholder values (1 SNP, 0 InDels) for pairs in Hamming mode, so \
                     SNP/InDel distances read from this cache may be wrong. Rebuild it (delete \
                     the file or use --force-recompute) unless it was written by cgdist >= 0.1.4."
                );
            }

            // Check distance mode compatibility - TEMPORARILY DISABLED
            // The cache now contains all 3 values (SNPs, indel_events, indel_bases)
            // so it can be used with any distance mode safely
            let expected_mode = format!("{distance_mode:?}");
            if modern_cache.metadata.distance_mode != expected_mode {
                println!(
                    "ℹ️  Cache was generated with mode '{}', using with mode '{}'",
                    modern_cache.metadata.distance_mode, expected_mode
                );
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
                        lens: (
                            cache_value.seq1_length.map(|l| l as u32),
                            cache_value.seq2_length.map(|l| l as u32),
                        ),
                    };
                    self.cache.insert(key, entry);
                }
            }

            if let Some(note) = &modern_cache.metadata.user_note {
                println!("📝 Cache note: {note}");
            }

            println!(
                "✅ Modern cache loaded: {} loci, {} hasher, {} mode",
                modern_cache.metadata.unique_loci,
                modern_cache.metadata.hasher_type,
                modern_cache.metadata.distance_mode
            );
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

            let legacy_cache: LegacyCache = bincode::deserialize(&decompressed).map_err(|e| {
                format!("Failed to deserialize cache (tried both modern and legacy formats): {e}")
            })?;

            // Check alignment configuration compatibility
            if !same_alignment_params(&legacy_cache.alignment_config, &self.config) {
                return Err(format!(
                    "Cache alignment config mismatch:\n  Cache: {:?}\n  Engine: {:?}",
                    legacy_cache.alignment_config, self.config
                ));
            }

            // Convert legacy format to internal format
            self.cache.clear();
            for ((locus, crc1, crc2, _config_hash), (snps, indel_events, indel_bases)) in
                legacy_cache.data
            {
                let key = DistanceCacheKey { locus, crc1, crc2 };
                let entry = CacheEntry {
                    snps,
                    indel_events,
                    indel_bases,
                    lens: (None, None), // legacy cache carries no sequence lengths
                };
                self.cache.insert(key, entry);
            }

            println!("⚠️  Loaded legacy cache format - consider regenerating with modern format");
        }

        let elapsed = start.elapsed();
        println!(
            "✅ Cache loaded in {:.2}s ({} entries, {} KB)",
            elapsed.as_secs_f64(),
            self.cache.len(),
            compressed.len() / 1024
        );

        // Reset the new entries flag since we just loaded from disk
        self.has_new_entries = false;

        Ok(())
    }

    /// Pre-compute all alignments for unique pairs in batch (cache-aware)
    #[allow(unknown_lints, clippy::manual_is_multiple_of)]
    pub fn precompute_alignments(
        &mut self,
        unique_pairs: &HashSet<(String, u32, u32)>,
        mode: DistanceMode,
    ) {
        let total_pairs = unique_pairs.len();

        // Hamming mode is answered without alignments (see get_distance). It
        // must never write placeholder statistics into the cache: a later
        // SNP/InDel run would read them back as real alignment results.
        if matches!(mode, DistanceMode::Hamming) {
            println!("🎯 Hamming mode: no alignments needed ({total_pairs} allele pairs)");
            return;
        }

        // Filter out pairs already in cache (Strategy 1: Preventive filtering)
        let start_filter = Instant::now();
        let missing_pairs: Vec<_> = unique_pairs
            .iter()
            .filter(|(locus, crc1, crc2)| {
                // Pre-compute min/max to avoid multiple comparisons
                let (min_crc, max_crc) = if *crc1 <= *crc2 {
                    (*crc1, *crc2)
                } else {
                    (*crc2, *crc1)
                };
                !self
                    .cache
                    .contains_key(&(locus.as_str(), min_crc, max_crc) as &dyn KeyView)
            })
            .collect();

        let filter_elapsed = start_filter.elapsed();
        let cached_pairs = total_pairs - missing_pairs.len();
        let missing_count = missing_pairs.len();

        println!(
            "🔍 Filtered unique pairs in {:.3}s:",
            filter_elapsed.as_secs_f64()
        );
        println!("   📊 Total pairs needed: {total_pairs}");
        println!(
            "   ✅ Already in cache: {} ({:.1}%)",
            cached_pairs,
            (cached_pairs as f64 / total_pairs as f64) * 100.0
        );
        println!(
            "   🔥 Missing pairs to compute: {} ({:.1}%)",
            missing_count,
            (missing_count as f64 / total_pairs as f64) * 100.0
        );

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

        // Compute alignments with periodic progress updates. Allele lengths
        // (known because both sequences were just aligned) are gathered in the
        // parallel part too, so the serial merge below is one insert per pair.
        // With --save-alignments the pairs are processed in chunks whose rows
        // are written out right away (same rows, same order), so memory does
        // not grow with the number of saved alignments.
        let chunk_size = if self.wants_details() {
            20_000
        } else {
            missing_pairs.len().max(1)
        };
        let mut unaligned: Vec<&(String, u32, u32)> = Vec::new();
        for chunk in missing_pairs.chunks(chunk_size) {
            #[allow(clippy::type_complexity)]
            let results: Vec<(
                &(String, u32, u32),
                Option<(
                    usize,
                    usize,
                    usize,
                    (Option<u32>, Option<u32>),
                    Option<PairDetail>,
                )>,
            )> = chunk
                .par_iter()
                .map(|&pair| {
                    let (locus, crc1, crc2) = pair;
                    let alignment_result = self.compute_single_alignment(locus, *crc1, *crc2);

                    // Increment and check if we should update progress
                    let completed =
                        completed_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                    if completed % update_frequency == 0 {
                        pb.set_position(completed as u64);
                    }

                    let with_lens =
                        alignment_result.map(|(snps, indel_events, indel_bases, detail)| {
                            let (lo, hi) = ((*crc1).min(*crc2), (*crc1).max(*crc2));
                            let lens = self.sequence_db.as_ref().map_or((None, None), |db| {
                                (
                                    db.get_sequence(locus, lo).map(|s| s.sequence.len() as u32),
                                    db.get_sequence(locus, hi).map(|s| s.sequence.len() as u32),
                                )
                            });
                            (snps, indel_events, indel_bases, lens, detail)
                        });
                    (pair, with_lens)
                })
                .collect();

            // Store results in cache (and write alignment detail rows for
            // --save-alignments). Pairs that could not be aligned are never
            // cached, so they are reported on every run until the schema is
            // fixed.
            for (pair, res) in results {
                let (locus, crc1, crc2) = pair;
                match res {
                    Some((snps, indel_events, indel_bases, lens, detail)) => {
                        self.cache.insert(
                            DistanceCacheKey {
                                locus: locus.clone(),
                                crc1: (*crc1).min(*crc2),
                                crc2: (*crc1).max(*crc2),
                            },
                            CacheEntry {
                                snps,
                                indel_events,
                                indel_bases,
                                lens,
                            },
                        );
                        self.has_new_entries = true;
                        if let Some(d) = detail {
                            if let (Some(row), Some(w)) = (d.full, self.alignments_out.as_mut()) {
                                w.write(&row);
                            }
                            if let (Some(row), Some(w)) = (d.cigar, self.cigar_out.as_mut()) {
                                w.write(&row);
                            }
                        }
                    }
                    None => unaligned.push(pair),
                }
            }
        }

        pb.finish_with_message("✅ Missing alignments computed!");

        let compute_elapsed = start_compute.elapsed();
        let aligned_count = missing_count - unaligned.len();

        println!("🚀 Cache-aware precompute completed:");
        println!(
            "   ⚡ Computed {} new alignments in {:.2}s ({:.0} alignments/sec)",
            aligned_count,
            compute_elapsed.as_secs_f64(),
            aligned_count as f64 / compute_elapsed.as_secs_f64()
        );
        if !unaligned.is_empty() {
            self.report_unaligned(&unaligned, total_pairs);
        }
        println!(
            "   📈 Total efficiency: {:.1}% time saved vs full recompute",
            (cached_pairs as f64 / total_pairs as f64) * 100.0
        );
        println!("   💾 Cache now contains {} entries", self.cache.len());
    }

    /// Warn about allele pairs that could not be aligned. They fall back to
    /// the cache-miss rule in `get_distance`, which undercounts differences.
    fn report_unaligned(&self, unaligned: &[&(String, u32, u32)], total_pairs: usize) {
        let mut absent: HashSet<(&str, u32)> = HashSet::new();
        let mut per_locus: HashMap<&str, usize> = HashMap::new();
        for (locus, c1, c2) in unaligned {
            *per_locus.entry(locus.as_str()).or_default() += 1;
            for c in [*c1, *c2] {
                let known = self
                    .sequence_db
                    .as_ref()
                    .is_some_and(|db| db.get_sequence(locus, c).is_some());
                if !known {
                    absent.insert((locus.as_str(), c));
                }
            }
        }
        let mut loci: Vec<(&str, usize)> = per_locus.into_iter().collect();
        loci.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(b.0)));
        let examples: Vec<String> = loci
            .iter()
            .take(5)
            .map(|(l, n)| format!("{l} ({n})"))
            .collect();

        eprintln!(
            "⚠️  WARNING: {} of {} allele pairs ({:.1}%) in {} loci could NOT be aligned.",
            unaligned.len(),
            total_pairs,
            unaligned.len() as f64 / total_pairs.max(1) as f64 * 100.0,
            loci.len()
        );
        if !absent.is_empty() {
            eprintln!(
                "   {} allele(s) found in the profiles have no sequence in the schema FASTA \
                 (profiles called with a newer or different schema?).",
                absent.len()
            );
        }
        eprintln!(
            "   These pairs count as 0 (1 in snps mode with --hamming-fallback), so distances \
             involving them are underestimated."
        );
        eprintln!("   Most affected loci: {}", examples.join(", "));
        eprintln!(
            "   Fix: use the schema the profiles were called with (including novel alleles)."
        );
    }

    /// Compute a single alignment (used in batch processing)  
    /// Returns None if sequences are not available (cache miss should apply Hamming fallback)
    fn compute_single_alignment(
        &self,
        locus: &str,
        crc1: u32,
        crc2: u32,
    ) -> Option<(usize, usize, usize, Option<PairDetail>)> {
        let result = self.compute_single_alignment_inner(locus, crc1, crc2);
        if self.verify_fraction > 0.0
            && crc1 != crc2
            && verify_selected(crc1, crc2, self.verify_fraction)
        {
            if let (Some((snps, ev, bases, _)), Some(db)) = (&result, &self.sequence_db) {
                if let (Some(s1), Some(s2)) =
                    (db.get_sequence(locus, crc1), db.get_sequence(locus, crc2))
                {
                    if let Some(want) =
                        reference_alignment_stats(&self.config, &s1.sequence, &s2.sequence)
                    {
                        if want != (*snps, *ev, *bases) {
                            panic!(
                                "alignment verification FAILED for locus {locus}, alleles {crc1}/{crc2}: \
                                 cgdist computed (snps, indel_events, indel_bases) = {:?}, parasail's \
                                 original kernel gives {:?}. Please report this at \
                                 https://github.com/genpat-it/cgDist/issues",
                                (*snps, *ev, *bases),
                                want
                            );
                        }
                    }
                }
            }
        }
        result
    }

    fn compute_single_alignment_inner(
        &self,
        locus: &str,
        crc1: u32,
        crc2: u32,
    ) -> Option<(usize, usize, usize, Option<PairDetail>)> {
        if crc1 == crc2 {
            return Some((0, 0, 0, None)); // Identical alleles
        }

        // Try to get sequences and align (for non-Hamming modes)
        if let Some(ref seq_db) = self.sequence_db {
            if let (Some(seq1), Some(seq2)) = (
                seq_db.get_sequence(locus, crc1),
                seq_db.get_sequence(locus, crc2),
            ) {
                // Fast path: certified banded alignment, bit-identical to
                // parasail (see core::banded). Sequences with a NUL byte are
                // left to parasail, whose C-string handling makes them fail
                // over to analyze_sequences.
                if !seq1.sequence.contains(&0) && !seq2.sequence.contains(&0) {
                    let scoring = Scoring {
                        match_score: self.config.match_score,
                        mismatch: self.config.mismatch_penalty,
                        gap_open: self.config.gap_open,
                        gap_extend: self.config.gap_extend,
                    };
                    if !self.wants_details() {
                        if let Some(b) =
                            align_certified(&seq1.sequence, &seq2.sequence, &scoring, 0.5)
                        {
                            return Some((b.snps, b.indel_events, b.indel_bases, None));
                        }
                    } else if seq1.sequence.is_ascii() && seq2.sequence.is_ascii() {
                        // --save-alignments / --save-cigar: the band's
                        // traceback yields the same gapped strings as
                        // parasail's (ASCII only: the parasail binding
                        // requires UTF-8). Statistics and rows are then built
                        // exactly as on the parasail path.
                        if let Some((b, st)) = align_certified_with_strings(
                            &seq1.sequence,
                            &seq2.sequence,
                            &scoring,
                            0.5,
                        ) {
                            let query = String::from_utf8_lossy(&st.query);
                            let reference = String::from_utf8_lossy(&st.reference);
                            let stats = compute_alignment_stats(&query, &reference);
                            let detail = self.pair_detail(
                                locus,
                                crc1,
                                crc2,
                                &seq1.sequence,
                                &seq2.sequence,
                                &query,
                                &reference,
                                stats,
                                b.score,
                            );
                            return Some((stats.0, stats.1, stats.2, Some(detail)));
                        }
                    }
                }

                // None only if the scoring matrix cannot be created
                let aligned = align_global_trace(&self.config, &seq1.sequence, &seq2.sequence)?;

                match aligned {
                    Ok(result) => {
                        // Get traceback strings with gaps
                        match result.get_traceback_strings(&seq1.sequence, &seq2.sequence) {
                            Ok(traceback) => {
                                // Use the proper alignment analysis with gaps
                                let (snps, indel_events, indel_bases) =
                                    compute_alignment_stats(&traceback.query, &traceback.reference);

                                // Build the detailed alignment row only when
                                // --save-alignments is active. Reading self here is
                                // fine in the parallel context; the row is returned and
                                // collected by the caller (which holds &mut self).
                                let detail = self.wants_details().then(|| {
                                    self.pair_detail(
                                        locus,
                                        crc1,
                                        crc2,
                                        &seq1.sequence,
                                        &seq2.sequence,
                                        &traceback.query,
                                        &traceback.reference,
                                        (snps, indel_events, indel_bases),
                                        result.get_score(),
                                    )
                                });

                                return Some((snps, indel_events, indel_bases, detail));
                            }
                            Err(_) => {
                                // Traceback failed, use simple approach
                                let (snps, indel_events, indel_bases) =
                                    self.analyze_sequences(&seq1.sequence, &seq2.sequence);
                                return Some((snps, indel_events, indel_bases, None));
                            }
                        }
                    }
                    Err(_) => {
                        // Alignment failed, use simple comparison
                        let (snps, indel_events, indel_bases) =
                            self.analyze_sequences(&seq1.sequence, &seq2.sequence);
                        return Some((snps, indel_events, indel_bases, None));
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
        let compressed =
            std::fs::read(cache_path).map_err(|e| format!("Failed to read cache file: {e}"))?;

        // Decompress with LZ4
        let decompressed = lz4_flex::decompress_size_prepended(&compressed)
            .map_err(|e| format!("Failed to decompress cache: {e}"))?;

        // Load as ModernCache
        let modern_cache: ModernCache = serde_json::from_slice(&decompressed)
            .map_err(|e| format!("Failed to deserialize cache: {e}"))?;

        // Check if any entries have sequence lengths
        for (_, value) in modern_cache.data.iter().take(5) {
            // Check first 5 entries
            if value.seq1_length.is_some() || value.seq2_length.is_some() {
                return Ok(true);
            }
        }

        Ok(false)
    }

    /// Enrich existing cache with nucleotide sequence lengths from schema (with input/output paths)
    pub fn enrich_cache_with_lengths_from_input(
        &mut self,
        schema_path: &str,
        input_cache_path: &str,
        output_cache_path: &str,
    ) -> Result<usize, String> {
        // Load sequence lengths from schema with CRC mapping
        println!("🔍 Loading schema lengths with CRC mapping from {schema_path}...");
        // Lengths are looked up per locus: CRC32 values collide across loci of
        // large schemas (e.g. ~900 in S. enterica), so a schema-wide CRC map
        // could assign the length of another locus' allele.
        let (schema_lengths, crc_by_locus) = self.load_schema_with_crc_mapping(schema_path)?;
        let crc_mappings: usize = crc_by_locus.values().map(|m| m.len()).sum();
        println!(
            "📊 Loaded {} loci from schema with {} CRC mappings",
            schema_lengths.len(),
            crc_mappings
        );

        // Read from input cache file
        if !std::path::Path::new(input_cache_path).exists() {
            return Err("Input cache file does not exist".to_string());
        }

        println!("📂 Loading cache from {input_cache_path}...");
        let compressed = std::fs::read(input_cache_path)
            .map_err(|e| format!("Failed to read cache file: {e}"))?;

        // Decompress with LZ4
        let decompressed = lz4_flex::decompress_size_prepended(&compressed)
            .map_err(|e| format!("Failed to decompress cache: {e}"))?;

        // Load as ModernCache
        let mut modern_cache: ModernCache = serde_json::from_slice(&decompressed)
            .map_err(|e| format!("Failed to deserialize cache: {e}"))?;

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
                    let locus_map = crc_by_locus.get(locus);
                    if let Some(&len1) = locus_map.and_then(|m| m.get(&crc1)) {
                        value.seq1_length = Some(len1);
                        found_any = true;
                    }

                    if let Some(&len2) = locus_map.and_then(|m| m.get(&crc2)) {
                        value.seq2_length = Some(len2);
                        found_any = true;
                    }

                    if found_any {
                        enriched_count += 1;
                    } else {
                        missing_entries.push(format!("{locus}:{crc1}:{crc2}"));
                    }
                } else {
                    missing_entries.push(format!(
                        "{locus}:{crc1_str}:{crc2_str} (invalid CRC format)"
                    ));
                }
            }
        }

        println!(
            "✅ Enriched {} out of {} entries",
            enriched_count,
            modern_cache.data.len()
        );

        if !missing_entries.is_empty() && missing_entries.len() <= 10 {
            println!(
                "⚠️  Warning: {} entries with missing alleles:",
                missing_entries.len()
            );
            for entry in &missing_entries {
                println!("   - {entry}");
            }
        } else if !missing_entries.is_empty() {
            println!(
                "⚠️  Warning: {} entries with missing alleles (showing first 10):",
                missing_entries.len()
            );
            for entry in missing_entries.iter().take(10) {
                println!("   - {entry}");
            }
        }

        // Update metadata
        modern_cache.metadata.last_modified = chrono::Utc::now()
            .format("%Y-%m-%d %H:%M:%S UTC")
            .to_string();
        if let Some(ref mut note) = modern_cache.metadata.user_note {
            note.push_str(" [Enriched with sequence lengths]");
        } else {
            modern_cache.metadata.user_note = Some("Enriched with sequence lengths".to_string());
        }

        // Serialize and save
        println!("💾 Saving enriched cache to {output_cache_path}...");
        let serialized = serde_json::to_vec(&modern_cache)
            .map_err(|e| format!("Failed to serialize cache: {e}"))?;

        // Compress with lz4_flex (matches cgDist format)
        let final_data = lz4_flex::compress_prepend_size(&serialized);

        std::fs::write(output_cache_path, final_data)
            .map_err(|e| format!("Failed to write cache file: {e}"))?;

        println!(
            "✅ Enriched cache saved ({:.1}% success rate)",
            (enriched_count as f64 / modern_cache.data.len() as f64) * 100.0
        );

        Ok(enriched_count)
    }

    /// Load schema with both lengths and CRC mapping for enrichment
    #[allow(clippy::type_complexity)]
    fn load_schema_with_crc_mapping(
        &self,
        schema_path: &str,
    ) -> Result<
        (
            HashMap<String, HashMap<String, usize>>,
            HashMap<String, HashMap<u32, usize>>,
        ),
        String,
    > {
        use std::fs;
        use std::path::Path;

        let schema_dir = Path::new(schema_path);
        let mut all_lengths = HashMap::new();
        let mut crc_by_locus = HashMap::new();

        // Get hasher for CRC calculation
        let registry = HasherRegistry::new();
        let hasher = registry
            .get_hasher(&self.hasher_type)
            .ok_or_else(|| format!("Unknown hasher type: {}", self.hasher_type))?;

        let entries = fs::read_dir(schema_dir)
            .map_err(|e| format!("Failed to read schema directory: {e}"))?;

        for entry in entries {
            let entry = entry.map_err(|e| format!("Failed to read directory entry: {e}"))?;
            let path = entry.path();

            if path.extension().and_then(|s| s.to_str()) == Some("fasta") {
                if let Some(filename) = path.file_stem().and_then(|s| s.to_str()) {
                    let locus_name = filename.to_string();
                    match self.load_fasta_with_crc_mapping(&path, hasher) {
                        Ok((lengths, crc_map)) => {
                            all_lengths.insert(locus_name.clone(), lengths);
                            crc_by_locus.insert(locus_name, crc_map);
                        }
                        Err(e) => {
                            eprintln!("⚠️  Warning: Failed to load {filename}: {e}");
                        }
                    }
                }
            }
        }

        Ok((all_lengths, crc_by_locus))
    }

    /// Load FASTA file with both lengths and CRC mapping
    #[allow(clippy::type_complexity)]
    fn load_fasta_with_crc_mapping(
        &self,
        fasta_path: &Path,
        hasher: &dyn AlleleHasher,
    ) -> Result<(HashMap<String, usize>, HashMap<u32, usize>), String> {
        use std::fs;

        let content = fs::read_to_string(fasta_path)
            .map_err(|e| format!("Failed to read FASTA file: {e}"))?;

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
                current_id = stripped.split_whitespace().next().unwrap_or("").to_string();

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
    /// Re-check a deterministic fraction of new alignments against parasail's
    /// original production kernel (striped, saturating); any difference
    /// aborts the run.
    pub fn set_verify_fraction(&mut self, fraction: f64) {
        self.verify_fraction = fraction.clamp(0.0, 1.0);
    }

    pub fn set_save_alignments(&mut self, path: String) {
        self.save_alignments_path = Some(path.clone());
        self.alignments_out = Some(RowWriter::new(path, ALIGNMENTS_HEADER));
    }

    /// Also write one CIGAR row per aligned pair (--save-cigar).
    pub fn set_save_cigar(&mut self, path: String) {
        self.cigar_out = Some(RowWriter::new(path, CIGAR_HEADER));
    }

    /// Whether per-pair output rows (gapped strings / CIGAR) are requested.
    fn wants_details(&self) -> bool {
        self.alignments_out.is_some() || self.cigar_out.is_some()
    }

    /// Build the requested output rows for one aligned pair.
    #[allow(clippy::too_many_arguments)]
    fn pair_detail(
        &self,
        locus: &str,
        crc1: u32,
        crc2: u32,
        seq1: &[u8],
        seq2: &[u8],
        query: &str,
        reference: &str,
        stats: (usize, usize, usize),
        score: i32,
    ) -> PairDetail {
        let (snps, indel_events, indel_bases) = stats;
        PairDetail {
            full: self.alignments_out.as_ref().map(|_| {
                format!(
                    "{locus}\t{crc1}\t{crc2}\t{}\t{}\t{}\t{}\t{snps}\t{indel_events}\t{indel_bases}\t{:.2}",
                    String::from_utf8_lossy(seq1),
                    String::from_utf8_lossy(seq2),
                    query,
                    reference,
                    score as f32,
                )
            }),
            cigar: self.cigar_out.as_ref().map(|_| {
                format!(
                    "{locus}\t{crc1}\t{crc2}\t{}\t{snps}\t{indel_events}\t{indel_bases}\t{:.2}",
                    cigar_from_aligned(query.as_bytes(), reference.as_bytes()),
                    score as f32,
                )
            }),
        }
    }

    /// Finish the per-pair output files (--save-alignments, --save-cigar).
    /// Each holds its header and one row per pair aligned in this run, and is
    /// written even when no pair was aligned.
    pub fn save_alignments(&mut self) -> Result<(), String> {
        let mut errors = Vec::new();
        for (out, what) in [
            (self.alignments_out.as_mut(), "alignment details"),
            (self.cigar_out.as_mut(), "CIGAR rows"),
        ] {
            if let Some(w) = out {
                match w.finish() {
                    Ok(n) => println!("💾 Saved {n} {what} to: {}", w.path),
                    Err(e) => errors.push(e),
                }
            }
        }
        if errors.is_empty() {
            Ok(())
        } else {
            Err(errors.join("; "))
        }
    }
}

const ALIGNMENTS_HEADER: &str = "locus\thash1\thash2\tseq1\tseq2\taligned_seq1\taligned_seq2\tsnps\tindel_events\tindel_bases\talignment_score";
const CIGAR_HEADER: &str =
    "locus\thash1\thash2\tcigar\tsnps\tindel_events\tindel_bases\talignment_score";

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
    calculate_sample_distance_detailed(
        sample1,
        sample2,
        loci_names,
        engine,
        mode,
        min_loci,
        no_hamming_fallback,
    )
    .0
}

/// Same computation as [`calculate_sample_distance`] but also returns the number
/// of shared (co-present in both samples) loci the distance was computed over.
/// Distance semantics are identical (`None` when `shared_loci < min_loci`); the
/// shared-loci count is always returned so callers can report per-pair data
/// quality (how many loci a distance actually rests on).
pub fn calculate_sample_distance_detailed(
    sample1: &AllelicProfile,
    sample2: &AllelicProfile,
    loci_names: &[String],
    engine: &DistanceEngine,
    mode: DistanceMode,
    min_loci: usize,
    no_hamming_fallback: bool,
) -> (Option<usize>, usize) {
    let d = calculate_sample_distance_full(
        sample1,
        sample2,
        loci_names,
        engine,
        mode,
        min_loci,
        no_hamming_fallback,
    );
    (d.0, d.1)
}

/// Full per-pair breakdown for reporting: `(distance, shared_loci, differing_loci,
/// sum_of_squared_contributions)`. `differing_loci` (h) and the sum of squared
/// per-locus contributions (q2) feed the missingness confidence interval; they
/// count only shared loci that actually contribute distance (> 0).
pub fn calculate_sample_distance_full(
    sample1: &AllelicProfile,
    sample2: &AllelicProfile,
    loci_names: &[String],
    engine: &DistanceEngine,
    mode: DistanceMode,
    min_loci: usize,
    no_hamming_fallback: bool,
) -> (Option<usize>, usize, usize, u64) {
    let mut total_distance = 0;
    let mut shared_loci = 0;
    let mut differing = 0usize;
    let mut sum_sq = 0u64;

    for locus in loci_names {
        let crc1 = sample1
            .loci_hashes
            .get(locus)
            .and_then(|h| h.as_crc32())
            .unwrap_or(u32::MAX);
        let crc2 = sample2
            .loci_hashes
            .get(locus)
            .and_then(|h| h.as_crc32())
            .unwrap_or(u32::MAX);

        if crc1 != u32::MAX && crc2 != u32::MAX {
            shared_loci += 1;
        }

        let dl = engine.get_distance(locus, crc1, crc2, mode, no_hamming_fallback);
        total_distance += dl;
        if dl > 0 {
            differing += 1;
            sum_sq += (dl as u64) * (dl as u64);
        }
    }

    let distance = if shared_loci >= min_loci {
        Some(total_distance)
    } else {
        None
    };
    (distance, shared_loci, differing, sum_sq)
}

/// One row of the long-format per-pair quality table (upper triangle, `i < j`).
/// `distance` is `None` when the pair fails the `min_loci` filter. `h` is the
/// number of differing shared loci and `q2` the sum of squared per-locus
/// contributions — both feed the optional missingness confidence interval.
#[derive(Debug, Clone, Copy)]
pub struct PairRow {
    pub i: usize,
    pub j: usize,
    pub distance: Option<usize>,
    pub shared: usize,
    pub h: usize,
    pub q2: u64,
}

/// Compute per-pair distances together with the number of shared loci, in long
/// format (upper triangle only). Data source for `--emit-pairs` and the HTML
/// dashboard. Parallelized like [`calculate_distance_matrix`]; distance lookups
/// hit the precomputed cache, so this is cheap despite re-deriving the pairs.
pub fn calculate_pairs_table(
    samples: &[AllelicProfile],
    loci_names: &[String],
    engine: &DistanceEngine,
    mode: DistanceMode,
    min_loci: usize,
    no_hamming_fallback: bool,
) -> Vec<PairRow> {
    let n_samples = samples.len();
    (0..n_samples)
        .into_par_iter()
        .flat_map(|i| {
            (i + 1..n_samples).into_par_iter().map(move |j| {
                let (distance, shared, h, q2) = calculate_sample_distance_full(
                    &samples[i],
                    &samples[j],
                    loci_names,
                    engine,
                    mode,
                    min_loci,
                    no_hamming_fallback,
                );
                PairRow {
                    i,
                    j,
                    distance,
                    shared,
                    h,
                    q2,
                }
            })
        })
        .collect()
}

/// Per-pair recombination load: number of loci whose mutation density exceeds
/// `thresh_frac`, in the same upper-triangle order as [`calculate_pairs_table`].
/// Returns `None` when the engine has no length data (non-enriched cache /
/// hamming hasher), so callers can hide the recombination view gracefully.
pub fn calculate_pairs_recombination(
    samples: &[AllelicProfile],
    loci_names: &[String],
    engine: &DistanceEngine,
    thresh_frac: f64,
) -> Option<Vec<u32>> {
    if !engine.has_recomb_data() {
        return None;
    }
    let n_samples = samples.len();
    let out = (0..n_samples)
        .into_par_iter()
        .flat_map(|i| {
            (i + 1..n_samples).into_par_iter().map(move |j| {
                let mut recomb = 0u32;
                for locus in loci_names {
                    let crc1 = samples[i]
                        .loci_hashes
                        .get(locus)
                        .and_then(|h| h.as_crc32())
                        .unwrap_or(u32::MAX);
                    let crc2 = samples[j]
                        .loci_hashes
                        .get(locus)
                        .and_then(|h| h.as_crc32())
                        .unwrap_or(u32::MAX);
                    if engine.locus_is_recombinant(locus, crc1, crc2, thresh_frac) == Some(true) {
                        recomb += 1;
                    }
                }
                recomb
            })
        })
        .collect();
    Some(out)
}

/// Calculate full distance matrix
#[allow(unknown_lints, clippy::manual_is_multiple_of)]
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
    println!(
        "🔄 Computing distance matrix ({n_samples} × {n_samples} = {total_comparisons} comparisons)..."
    );

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

    let upper_triangle: Vec<_> = (0..n_samples)
        .into_par_iter()
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
    println!(
        "✅ Distance matrix computed in {:.2}s",
        elapsed.as_secs_f64()
    );

    matrix
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::SequenceInfo;

    fn crc(s: &[u8]) -> u32 {
        let mut h = crc32fast::Hasher::new();
        h.update(s);
        h.finalize()
    }

    fn engine_with(seqs: &[&[u8]]) -> (DistanceEngine, Vec<u32>) {
        let mut db = SequenceDatabase::new();
        let mut crcs = Vec::new();
        for (i, s) in seqs.iter().enumerate() {
            let c = crc(s);
            crcs.push(c);
            db.add_sequence(
                "L1".to_string(),
                c,
                SequenceInfo {
                    sequence: s.to_vec(),
                    id: format!("L1_{i}"),
                },
            );
        }
        let e = DistanceEngine::with_sequences(AlignmentConfig::default(), db, "crc32".into());
        (e, crcs)
    }

    const A: &[u8] = b"ACGTACGTACGTACGTACGT";
    const B: &[u8] = b"ACGTACGAACGTACGTTCGT"; // 2 SNPs vs A

    fn pairs(c: &[u32]) -> HashSet<(String, u32, u32)> {
        [("L1".to_string(), c[0].min(c[1]), c[0].max(c[1]))].into()
    }

    #[test]
    fn hamming_mode_never_writes_placeholders_into_cache() {
        let (mut e, c) = engine_with(&[A, B]);
        e.precompute_alignments(&pairs(&c), DistanceMode::Hamming);
        assert_eq!(e.cache_stats().0, 0);
        assert!(!e.has_new_entries());
        assert_eq!(
            e.get_distance("L1", c[0], c[1], DistanceMode::Hamming, true),
            1
        );

        // A later SNP run must align instead of reading a Hamming placeholder.
        e.precompute_alignments(&pairs(&c), DistanceMode::SnpsOnly);
        assert_eq!(
            e.get_distance("L1", c[0], c[1], DistanceMode::SnpsOnly, true),
            2
        );
        assert_eq!(
            e.get_distance("L1", c[0], c[1], DistanceMode::Hamming, true),
            1
        );
    }

    #[test]
    fn computed_pairs_keep_lengths_through_save_and_load() {
        let (mut e, c) = engine_with(&[A, B]);
        e.precompute_alignments(&pairs(&c), DistanceMode::SnpsOnly);
        assert!(e.has_recomb_data());
        assert_eq!(e.locus_is_recombinant("L1", c[0], c[1], 0.05), Some(true));

        let path = std::env::temp_dir().join(format!("cgdist_len_test_{}.lz4", std::process::id()));
        let path = path.to_str().unwrap();
        e.save_cache(path, DistanceMode::SnpsOnly).unwrap();
        let (mut e2, _) = engine_with(&[A, B]);
        e2.load_cache(path, DistanceMode::SnpsOnly).unwrap();
        let _ = std::fs::remove_file(path);
        assert!(e2.has_recomb_data());
        assert_eq!(e2.locus_is_recombinant("L1", c[0], c[1], 0.05), Some(true));
        assert_eq!(e2.locus_is_recombinant("L1", c[0], c[1], 0.2), Some(false));
    }

    /// The previous production kernel (striped, saturating 8->16->32 bit).
    fn striped_reference(q: &[u8], r: &[u8]) -> (usize, usize, usize, i32) {
        let c = AlignmentConfig::default();
        let m = Matrix::create(b"ACGT", c.match_score, c.mismatch_penalty).unwrap();
        let a = Aligner::new()
            .matrix(m)
            .gap_open(c.gap_open)
            .gap_extend(c.gap_extend)
            .global()
            .use_trace()
            .build();
        let res = a.align(Some(q), r).unwrap();
        let tb = res.get_traceback_strings(q, r).unwrap();
        let (s, e, b) = compute_alignment_stats(&tb.query, &tb.reference);
        (s, e, b, res.get_score())
    }

    fn fast(q: &[u8], r: &[u8]) -> (usize, usize, usize, i32, bool) {
        let res = align_global_trace(&AlignmentConfig::default(), q, r)
            .unwrap()
            .unwrap();
        let tb = res.get_traceback_strings(q, r).unwrap();
        let (s, e, b) = compute_alignment_stats(&tb.query, &tb.reference);
        (s, e, b, res.get_score(), res.is_saturated())
    }

    fn lcg_seq(seed: &mut u64, n: usize) -> Vec<u8> {
        (0..n)
            .map(|_| {
                *seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(*seed >> 62) as usize]
            })
            .collect()
    }

    #[test]
    fn scan_kernel_matches_striped_including_32bit_fallback() {
        let mut seed = 42u64;
        let mut cases: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
        for &n in &[1usize, 7, 60, 300, 1200] {
            let a = lcg_seq(&mut seed, n);
            let mut b = a.clone();
            // sprinkle substitutions and an indel
            for k in (0..b.len()).step_by(37) {
                b[k] = if b[k] == b'A' { b'C' } else { b'A' };
            }
            if b.len() > 20 {
                b.drain(10..13);
            }
            cases.push((a.clone(), b));
            cases.push((a, lcg_seq(&mut seed, n + 5)));
        }
        // Long unrelated sequences: 16-bit must saturate and fall back to 32-bit.
        let long_a = lcg_seq(&mut seed, 20_000);
        let long_b = lcg_seq(&mut seed, 19_000);
        let fl = fast(&long_a, &long_b);
        assert!(!fl.4, "the returned result must come from the 32-bit rerun");
        cases.push((long_a, long_b));

        for (q, r) in &cases {
            let f = fast(q, r);
            assert_eq!(
                (f.0, f.1, f.2, f.3),
                striped_reference(q, r),
                "len {}",
                q.len()
            );
        }
    }

    #[test]
    fn long_pair_saturates_16bit() {
        let mut seed = 7u64;
        let a = lcg_seq(&mut seed, 20_000);
        let b = lcg_seq(&mut seed, 19_000);
        let m = Matrix::create(b"ACGT", 2, -1).unwrap();
        let a16 = Aligner::new()
            .matrix(m)
            .gap_open(5)
            .gap_extend(2)
            .global()
            .use_trace()
            .scan()
            .solution_width(16)
            .build();
        assert!(a16.align(Some(&a), &b).unwrap().is_saturated());

        // The 64-bit rung exists and agrees with 32-bit.
        let stats = |w: i32| {
            let m = Matrix::create(b"ACGT", 2, -1).unwrap();
            let al = Aligner::new()
                .matrix(m)
                .gap_open(5)
                .gap_extend(2)
                .global()
                .use_trace()
                .scan()
                .solution_width(w)
                .build();
            let res = al.align(Some(&a), &b).unwrap();
            assert!(!res.is_saturated());
            let tb = res.get_traceback_strings(&a, &b).unwrap();
            (
                compute_alignment_stats(&tb.query, &tb.reference),
                res.get_score(),
            )
        };
        assert_eq!(stats(32), stats(64));
    }

    #[test]
    fn verify_selection_is_deterministic_and_proportional() {
        let mut hits = 0;
        for i in 0..20_000u32 {
            let (a, b) = (i.wrapping_mul(2654435761), i ^ 0xDEADBEEF);
            let x = verify_selected(a, b, 0.1);
            assert_eq!(x, verify_selected(b, a, 0.1)); // order-independent
            assert_eq!(x, verify_selected(a, b, 0.1)); // repeatable
            hits += x as usize;
            assert!(verify_selected(a, b, 1.0));
            assert!(!verify_selected(a, b, 0.0));
        }
        assert!(
            (1600..2400).contains(&hits),
            "selected {hits} of 20000 at 10%"
        );
    }

    #[test]
    fn verification_accepts_correct_alignments() {
        let (mut e, c) = engine_with(&[A, B]);
        e.set_verify_fraction(1.0);
        e.precompute_alignments(&pairs(&c), DistanceMode::SnpsOnly);
        assert_eq!(
            e.get_distance("L1", c[0], c[1], DistanceMode::SnpsOnly, true),
            2
        );
        let want = reference_alignment_stats(&AlignmentConfig::default(), A, B);
        assert_eq!(want, Some((2, 0, 0)));
    }

    #[test]
    fn alignment_params_compared_numerically() {
        let dna = AlignmentConfig::from_mode("dna").unwrap();
        assert!(same_alignment_params(&dna, &AlignmentConfig::default()));
        let strict = AlignmentConfig::from_mode("dna-strict").unwrap();
        assert!(!same_alignment_params(&dna, &strict));
        assert!(written_before_hamming_fix("0.1.3"));
        assert!(written_before_hamming_fix("0.1.2-beta"));
        assert!(!written_before_hamming_fix("0.1.4"));
        assert!(!written_before_hamming_fix("0.2.0"));
    }
}
