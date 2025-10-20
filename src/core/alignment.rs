// alignment.rs - Alignment configuration and utilities

use serde::{Deserialize, Serialize};
use std::str::FromStr;

/// Configuration for sequence alignment
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AlignmentConfig {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub description: Option<String>,
}

impl Default for AlignmentConfig {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_penalty: -1,
            gap_open: 5,
            gap_extend: 2,
            description: Some("Default DNA alignment parameters".to_string()),
        }
    }
}

impl AlignmentConfig {
    /// Create configuration from mode string
    pub fn from_mode(mode: &str) -> Result<Self, String> {
        match mode {
            "dna" => Ok(Self {
                match_score: 2,
                mismatch_penalty: -1,
                gap_open: 5,
                gap_extend: 2,
                description: Some("Standard DNA alignment".to_string()),
            }),
            "dna-strict" => Ok(Self {
                match_score: 3,
                mismatch_penalty: -2,
                gap_open: 8,
                gap_extend: 3,
                description: Some("Strict DNA alignment (higher penalties)".to_string()),
            }),
            "dna-permissive" => Ok(Self {
                match_score: 1,
                mismatch_penalty: 0,
                gap_open: 3,
                gap_extend: 1,
                description: Some("Permissive DNA alignment (lower penalties)".to_string()),
            }),
            _ => Err(format!("Unknown alignment mode: {}", mode)),
        }
    }

    /// Create custom configuration
    pub fn custom(match_score: i32, mismatch_penalty: i32, gap_open: i32, gap_extend: i32) -> Self {
        Self {
            match_score,
            mismatch_penalty,
            gap_open,
            gap_extend,
            description: Some("Custom alignment parameters".to_string()),
        }
    }
}

/// Detailed alignment result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DetailedAlignment {
    pub locus: String,
    pub crc1: u32,
    pub crc2: u32,
    pub seq1_id: String,
    pub seq2_id: String,
    pub query_aligned: String,
    pub reference_aligned: String,
    pub alignment_score: i32,
    pub snps: usize,
    pub indel_events: usize,
    pub indel_bases: usize,
    pub alignment_length: usize,
    pub identity_percent: f64,
}

/// Distance calculation mode
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DistanceMode {
    SnpsOnly,
    SnpsAndIndelEvents,
    SnpsAndIndelBases,
    Hamming,
}

impl FromStr for DistanceMode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "snps" | "snps-only" => Ok(DistanceMode::SnpsOnly),
            "snps-indel-events" | "snps+indel-events" => Ok(DistanceMode::SnpsAndIndelEvents),
            "snps-indel-bases" | "snps+indel-bases" => Ok(DistanceMode::SnpsAndIndelBases),
            "hamming" => Ok(DistanceMode::Hamming),
            _ => Err(format!("Invalid distance mode: {}. Use: snps, snps-indel-events, snps-indel-bases, hamming", s))
        }
    }
}

impl DistanceMode {
    pub fn description(&self) -> &str {
        match self {
            DistanceMode::SnpsOnly => "SNPs only",
            DistanceMode::SnpsAndIndelEvents => "SNPs + indel events",
            DistanceMode::SnpsAndIndelBases => "SNPs + indel bases",
            DistanceMode::Hamming => "Hamming distance (all mismatches)",
        }
    }
}

/// Compute alignment statistics from aligned sequences
pub fn compute_alignment_stats(query: &str, reference: &str) -> (usize, usize, usize) {
    let query_bytes = query.as_bytes();
    let ref_bytes = reference.as_bytes();

    let mut snps = 0;
    let mut indel_events = 0;
    let mut indel_bases = 0;
    let mut in_gap = false;

    for i in 0..query_bytes.len().min(ref_bytes.len()) {
        let q = query_bytes[i];
        let r = ref_bytes[i];

        if q == b'-' || r == b'-' {
            if !in_gap {
                indel_events += 1;
                in_gap = true;
            }
            indel_bases += 1;
        } else {
            in_gap = false;
            if q != r {
                snps += 1;
            }
        }
    }

    (snps, indel_events, indel_bases)
}

/// Compute Hamming distance between two sequences
/// Counts all mismatches, treating gaps as mismatches
pub fn compute_hamming_distance(seq1: &[u8], seq2: &[u8]) -> usize {
    // Handle sequences of different lengths
    let min_len = seq1.len().min(seq2.len());
    let max_len = seq1.len().max(seq2.len());

    // Count mismatches in overlapping region
    let mismatches = seq1
        .iter()
        .take(min_len)
        .zip(seq2.iter().take(min_len))
        .filter(|(a, b)| a != b)
        .count();

    // Add length difference as additional mismatches
    mismatches + (max_len - min_len)
}
