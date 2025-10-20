// crc32.rs - CRC32 hasher implementation

use super::traits::{AlleleHash, AlleleHasher};

/// CRC32 hasher - maintains compatibility with chewBACCA
#[derive(Debug, Clone)]
pub struct Crc32Hasher;

impl AlleleHasher for Crc32Hasher {
    fn hash_sequence(&self, sequence: &str) -> AlleleHash {
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

        let crc = cleaned
            .parse::<u32>()
            .map_err(|_| format!("Failed to parse '{}' as CRC32", cleaned))?;

        Ok(AlleleHash::from_crc32(crc))
    }

    fn name(&self) -> &'static str {
        "CRC32"
    }

    fn description(&self) -> &'static str {
        "CRC32 hash compatible with chewBACCA allele calling"
    }
}
