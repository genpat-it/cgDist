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
            _ => Err(format!("Unknown alignment mode: {mode}")),
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
    /// nonsynonymous SNPs only (needs coding counts, see core::protein)
    NonsynSnps,
    /// user-defined weighted sum of per-pair counts (see DistanceWeights)
    Weighted,
    /// 1 per locus whose alleles encode different proteins
    AaHamming,
    /// amino-acid substitutions
    AaSubs,
    /// amino-acid substitutions + InDel events of the protein alignment
    AaSubsIndelEvents,
    /// amino-acid substitutions + inserted/deleted residues
    AaSubsIndelResidues,
}

impl DistanceMode {
    /// Needs protein-level results (core::protein_distance).
    pub fn is_protein(self) -> bool {
        matches!(
            self,
            DistanceMode::AaHamming
                | DistanceMode::AaSubs
                | DistanceMode::AaSubsIndelEvents
                | DistanceMode::AaSubsIndelResidues
        )
    }
}

/// Per-locus contribution of a pair of different alleles in `--mode custom`:
/// the sum of `weight * count` over the counts below. Integer weights keep
/// distances integer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct DistanceWeights {
    /// 1 per locus whose alleles differ (Hamming unit)
    pub allele: u32,
    pub snps: u32,
    pub indel_events: u32,
    pub indel_bases: u32,
    /// synonymous SNPs (synonymous, stop_retained, start_retained)
    pub syn: u32,
    /// nonsynonymous SNPs (missense, stop_gained, stop_lost, start_lost)
    pub nonsyn: u32,
    /// SNPs in codons shifted or split by an InDel
    pub frame_disrupted: u32,
    /// 1 per locus whose alleles encode different proteins
    pub aa_allele: u32,
    /// amino-acid substitutions (protein alignment)
    pub aa_subs: u32,
    /// InDel events of the protein alignment
    pub aa_indel_events: u32,
    /// inserted/deleted residues of the protein alignment
    pub aa_indel_residues: u32,
}

impl DistanceWeights {
    pub const KEYS: [&'static str; 11] = [
        "allele",
        "snps",
        "indel_events",
        "indel_bases",
        "syn",
        "nonsyn",
        "frame_disrupted",
        "aa_allele",
        "aa_subs",
        "aa_indel_events",
        "aa_indel_residues",
    ];

    /// Parse "key=weight,key=weight" (keys in `KEYS`, non-negative integers).
    pub fn parse(spec: &str) -> Result<Self, String> {
        let mut w = DistanceWeights::default();
        let mut any = false;
        for part in spec.split(',').map(str::trim).filter(|p| !p.is_empty()) {
            let (k, v) = part
                .split_once('=')
                .ok_or_else(|| format!("--weights: expected key=value, got '{part}'"))?;
            let v: u32 = v.trim().parse().map_err(|_| {
                format!(
                    "--weights: weight of '{}' must be a non-negative integer, got '{}'",
                    k.trim(),
                    v.trim()
                )
            })?;
            let slot = match k.trim() {
                "allele" => &mut w.allele,
                "snps" => &mut w.snps,
                "indel_events" => &mut w.indel_events,
                "indel_bases" => &mut w.indel_bases,
                "syn" => &mut w.syn,
                "nonsyn" => &mut w.nonsyn,
                "frame_disrupted" => &mut w.frame_disrupted,
                "aa_allele" => &mut w.aa_allele,
                "aa_subs" => &mut w.aa_subs,
                "aa_indel_events" => &mut w.aa_indel_events,
                "aa_indel_residues" => &mut w.aa_indel_residues,
                other => {
                    return Err(format!(
                        "--weights: unknown key '{other}' (use: {})",
                        Self::KEYS.join(", ")
                    ))
                }
            };
            *slot = v;
            any = true;
        }
        if !any {
            return Err("--weights: no weights given".to_string());
        }
        Ok(w)
    }

    /// Whether any weight needs the synonymous/nonsynonymous classification.
    pub fn needs_coding(&self) -> bool {
        self.syn > 0 || self.nonsyn > 0 || self.frame_disrupted > 0
    }

    /// Whether any weight needs a DNA alignment.
    pub fn needs_alignment(&self) -> bool {
        self.snps > 0 || self.indel_events > 0 || self.indel_bases > 0 || self.needs_coding()
    }

    /// Whether any weight needs protein-level results.
    pub fn needs_protein(&self) -> bool {
        self.aa_allele > 0
            || self.aa_subs > 0
            || self.aa_indel_events > 0
            || self.aa_indel_residues > 0
    }

    pub fn describe(&self) -> String {
        let vals = [
            self.allele,
            self.snps,
            self.indel_events,
            self.indel_bases,
            self.syn,
            self.nonsyn,
            self.frame_disrupted,
            self.aa_allele,
            self.aa_subs,
            self.aa_indel_events,
            self.aa_indel_residues,
        ];
        Self::KEYS
            .iter()
            .zip(vals)
            .filter(|(_, v)| *v > 0)
            .map(|(k, v)| format!("{v}*{k}"))
            .collect::<Vec<_>>()
            .join(" + ")
    }
}

impl FromStr for DistanceMode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "snps" | "snps-only" => Ok(DistanceMode::SnpsOnly),
            // "snps-indel-contiguous" is the primary name; "snps-indel-events" /
            // "snps+indel-events" are kept as backward-compat aliases.
            "snps-indel-contiguous" | "snps-indel-events" | "snps+indel-events" => {
                Ok(DistanceMode::SnpsAndIndelEvents)
            }
            "snps-indel-bases" | "snps+indel-bases" => Ok(DistanceMode::SnpsAndIndelBases),
            "hamming" => Ok(DistanceMode::Hamming),
            "nonsyn-snps" | "nonsynonymous-snps" => Ok(DistanceMode::NonsynSnps),
            "custom" => Ok(DistanceMode::Weighted),
            "aa-hamming" => Ok(DistanceMode::AaHamming),
            "aa-substitutions" => Ok(DistanceMode::AaSubs),
            "aa-substitutions-indel-events" => Ok(DistanceMode::AaSubsIndelEvents),
            "aa-substitutions-indel-residues" => Ok(DistanceMode::AaSubsIndelResidues),
            _ => Err(format!("Invalid distance mode: {s}. Use: snps, snps-indel-contiguous, snps-indel-bases, hamming, nonsyn-snps, aa-hamming, aa-substitutions, aa-substitutions-indel-events, aa-substitutions-indel-residues, custom"))
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
            DistanceMode::NonsynSnps => "nonsynonymous SNPs",
            DistanceMode::Weighted => "custom weighted counts",
            DistanceMode::AaHamming => "loci with different proteins",
            DistanceMode::AaSubs => "amino-acid substitutions",
            DistanceMode::AaSubsIndelEvents => "amino-acid substitutions + InDel events",
            DistanceMode::AaSubsIndelResidues => "amino-acid substitutions + InDel residues",
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

/// CIGAR of a global alignment given as two gapped strings, in parasail's
/// extended format (`parasail_result_get_cigar`): `=` match, `X` mismatch,
/// `I` a query base against a reference gap, `D` a reference base against a
/// query gap; runs are merged. As in parasail, bases are compared
/// case-insensitively (so `a`/`A` is `=`), unlike the byte-wise SNP count of
/// `compute_alignment_stats`.
pub fn cigar_from_aligned(query: &[u8], reference: &[u8]) -> String {
    let mut out = String::new();
    let mut run_op = 0u8;
    let mut run_len = 0usize;
    for (&q, &r) in query.iter().zip(reference) {
        let op = if q == b'-' {
            b'D'
        } else if r == b'-' {
            b'I'
        } else if q.eq_ignore_ascii_case(&r) {
            b'='
        } else {
            b'X'
        };
        if op == run_op {
            run_len += 1;
        } else {
            if run_len > 0 {
                out.push_str(&run_len.to_string());
                out.push(run_op as char);
            }
            run_op = op;
            run_len = 1;
        }
    }
    if run_len > 0 {
        out.push_str(&run_len.to_string());
        out.push(run_op as char);
    }
    out
}

#[cfg(test)]
mod cigar_tests {
    use super::cigar_from_aligned;

    #[test]
    fn cigar_runs_and_ops() {
        // SNP at column 3, query-gap (D) at 6, reference-gap (I) at 8
        assert_eq!(
            cigar_from_aligned(b"ATGCA-TTG", b"ATCCATT-G"),
            "2=1X2=1D1=1I1="
        );
        assert_eq!(cigar_from_aligned(b"ACGT", b"ACGT"), "4=");
        assert_eq!(cigar_from_aligned(b"--AC", b"GGAC"), "2D2=");
        assert_eq!(cigar_from_aligned(b"ACGG", b"AC--"), "2=2I");
        // parasail compares case-insensitively; N matches N
        assert_eq!(cigar_from_aligned(b"acNT", b"ACNA"), "3=1X");
        assert_eq!(cigar_from_aligned(b"", b""), "");
    }
}
