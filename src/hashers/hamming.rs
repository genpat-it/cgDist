// hamming.rs - Hamming distance hasher for CRC allelic comparison

use super::traits::{AlleleHasher, AlleleHash};

/// Hamming distance hasher - works at CRC allelic level only
/// For CRC alleles: different CRCs = distance 1, same CRCs = distance 0
/// Does not perform sequence alignment, only CRC comparison
#[derive(Debug, Clone)]
pub struct HammingHasher;

impl AlleleHasher for HammingHasher {
    fn hash_sequence(&self, sequence: &str) -> AlleleHash {
        // Use CRC32 for sequence hashing (maintains compatibility)
        use crc32fast::Hasher;
        let mut hasher = Hasher::new();
        hasher.update(sequence.as_bytes());
        let crc = hasher.finalize();
        AlleleHash::from_crc32(crc)
    }
    
    fn parse_allele(&self, allele_str: &str, missing_char: &str) -> Result<AlleleHash, String> {
        let cleaned = allele_str.trim();
        
        if cleaned.is_empty() || cleaned == "NA" || cleaned == missing_char {
            return Ok(AlleleHash::Missing);
        }
        
        let crc = cleaned.parse::<u32>()
            .map_err(|_| format!("Failed to parse '{}' as CRC32 for Hamming hasher", cleaned))?;
        
        Ok(AlleleHash::from_crc32(crc))
    }
    
    fn name(&self) -> &'static str {
        "Hamming"
    }
    
    fn description(&self) -> &'static str {
        "Hamming distance at CRC allelic level (different CRCs = 1, same CRCs = 0)"
    }
    
    fn validate_sequence(&self, _sequence: &str) -> Result<(), String> {
        // Hamming hasher accepts any sequence since it works at CRC level
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hamming_hasher_sequence() {
        let hasher = HammingHasher;
        let hash1 = hasher.hash_sequence("ATCG");
        let hash2 = hasher.hash_sequence("ATCG");
        let hash3 = hasher.hash_sequence("GCTA");
        
        assert_eq!(hash1, hash2);
        assert_ne!(hash1, hash3);
    }
    
    #[test]
    fn test_hamming_parse_allele() {
        let hasher = HammingHasher;
        
        // Test normal CRC parsing
        let result = hasher.parse_allele("12345", "-").unwrap();
        assert_eq!(result, AlleleHash::Crc32(12345));
        
        // Test missing allele
        let missing = hasher.parse_allele("-", "-").unwrap();
        assert!(missing.is_missing());
        
        let missing_na = hasher.parse_allele("NA", "-").unwrap();
        assert!(missing_na.is_missing());
        
        // Test invalid input
        let error = hasher.parse_allele("invalid", "-");
        assert!(error.is_err());
    }
    
    #[test]
    fn test_hamming_properties() {
        let hasher = HammingHasher;
        assert_eq!(hasher.name(), "Hamming");
        assert!(hasher.description().contains("CRC allelic level"));
        assert!(hasher.validate_sequence("ATCGATCG").is_ok());
    }
}