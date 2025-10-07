// sha256.rs - SHA256 hasher implementation

use super::traits::{AlleleHasher, AlleleHash};

/// SHA256 hasher - cryptographically secure alternative
#[derive(Debug, Clone)]
pub struct Sha256Hasher;

impl AlleleHasher for Sha256Hasher {
    fn hash_sequence(&self, sequence: &str) -> AlleleHash {
        use sha2::{Sha256, Digest};
        let mut hasher = Sha256::new();
        hasher.update(sequence.as_bytes());
        let hash_str = format!("{:x}", hasher.finalize());
        AlleleHash::String(hash_str)
    }
    
    fn parse_allele(&self, allele_str: &str, missing_char: &str) -> Result<AlleleHash, String> {
        let cleaned = allele_str.trim();
        
        if cleaned.is_empty() || cleaned == "NA" || cleaned == missing_char {
            return Ok(AlleleHash::Missing);
        }
        
        // For SHA256, we expect the allele to be either a sequence or already a hash
        if cleaned.len() == 64 && cleaned.chars().all(|c| c.is_ascii_hexdigit()) {
            // Already a SHA256 hash
            Ok(AlleleHash::String(cleaned.to_string()))
        } else {
            // Assume it's a sequence, compute the hash
            Ok(self.hash_sequence(cleaned))
        }
    }
    
    fn name(&self) -> &'static str {
        "SHA256"
    }
    
    fn description(&self) -> &'static str {
        "SHA256 hash for cryptographically secure allele identification"
    }
}