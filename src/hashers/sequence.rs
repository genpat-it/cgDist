// sequence.rs - Sequence hasher implementation

use super::traits::{AlleleHash, AlleleHasher};

/// Simple string hasher - uses the sequence itself as identifier
#[derive(Debug, Clone)]
pub struct SequenceHasher;

impl AlleleHasher for SequenceHasher {
    fn hash_sequence(&self, sequence: &str) -> AlleleHash {
        AlleleHash::String(sequence.to_string())
    }

    fn parse_allele(&self, allele_str: &str, missing_char: &str) -> Result<AlleleHash, String> {
        let cleaned = allele_str.trim();

        if cleaned.is_empty() || cleaned == "NA" || cleaned == missing_char {
            return Ok(AlleleHash::Missing);
        }

        // For sequence hasher, we use the allele string directly
        Ok(AlleleHash::String(cleaned.to_string()))
    }

    fn name(&self) -> &'static str {
        "SEQUENCE"
    }

    fn description(&self) -> &'static str {
        "Uses the actual sequence as identifier (no hashing)"
    }
}
