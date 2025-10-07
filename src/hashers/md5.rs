// md5.rs - MD5 hasher implementation

use super::traits::{AlleleHasher, AlleleHash};

/// MD5 hasher - legacy compatibility
#[derive(Debug, Clone)]
pub struct Md5Hasher;

impl AlleleHasher for Md5Hasher {
    fn hash_sequence(&self, sequence: &str) -> AlleleHash {
        let digest = md5::compute(sequence.as_bytes());
        let hash_str = format!("{:x}", digest);
        AlleleHash::String(hash_str)
    }
    
    fn parse_allele(&self, allele_str: &str, missing_char: &str) -> Result<AlleleHash, String> {
        let cleaned = allele_str.trim();
        
        if cleaned.is_empty() || cleaned == "NA" || cleaned == missing_char {
            return Ok(AlleleHash::Missing);
        }
        
        // For MD5, we expect the allele to be either a sequence or already a hash
        if cleaned.len() == 32 && cleaned.chars().all(|c| c.is_ascii_hexdigit()) {
            // Already an MD5 hash
            Ok(AlleleHash::String(cleaned.to_string()))
        } else {
            // Assume it's a sequence, compute the hash
            Ok(self.hash_sequence(cleaned))
        }
    }
    
    fn name(&self) -> &'static str {
        "MD5"
    }
    
    fn description(&self) -> &'static str {
        "MD5 hash for legacy compatibility"
    }
}