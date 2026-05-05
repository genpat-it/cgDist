// mod.rs - Hashers module root

pub mod crc32;
pub mod hamming;
pub mod md5;
pub mod registry;
pub mod sequence;
pub mod sha256;
pub mod traits;

// Re-export main types for convenience
pub use crc32::Crc32Hasher;
pub use hamming::HammingHasher;
pub use md5::Md5Hasher;
pub use registry::HasherRegistry;
pub use sequence::SequenceHasher;
pub use sha256::Sha256Hasher;
pub use traits::{AlleleHash, AlleleHashPair, AlleleHasher};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crc32_hasher() {
        let hasher = Crc32Hasher;
        let hash1 = hasher.hash_sequence("ATCG");
        let hash2 = hasher.hash_sequence("ATCG");
        let hash3 = hasher.hash_sequence("GCTA");

        assert_eq!(hash1, hash2);
        assert_ne!(hash1, hash3);
        assert_eq!(hasher.name(), "CRC32");
    }

    #[test]
    fn test_sha256_hasher() {
        let hasher = Sha256Hasher;
        let hash1 = hasher.hash_sequence("ATCG");
        let hash2 = hasher.hash_sequence("ATCG");
        let hash3 = hasher.hash_sequence("GCTA");

        assert_eq!(hash1, hash2);
        assert_ne!(hash1, hash3);
        assert_eq!(hasher.name(), "SHA256");
    }

    #[test]
    fn test_allele_hash_pair() {
        let hash1 = AlleleHash::Crc32(100);
        let hash2 = AlleleHash::Crc32(200);
        let pair1 = AlleleHashPair::new(hash1.clone(), hash2.clone());
        let pair2 = AlleleHashPair::new(hash2, hash1);

        assert_eq!(pair1, pair2); // Should be ordered
        assert_eq!(pair1.hash1, AlleleHash::Crc32(100));
        assert_eq!(pair1.hash2, AlleleHash::Crc32(200));
    }

    #[test]
    fn test_allele_hash_enum() {
        let crc_hash = AlleleHash::from_crc32(12345);
        let string_hash = AlleleHash::from_string("test".to_string(), "MISSING");
        let missing_hash = AlleleHash::Missing;

        assert_eq!(crc_hash.as_crc32(), Some(12345));
        assert_eq!(string_hash.as_string(), Some("test"));
        assert!(missing_hash.is_missing());

        // Test missing marker
        let missing_from_crc = AlleleHash::from_crc32(u32::MAX);
        let missing_from_string = AlleleHash::from_string("MISSING".to_string(), "MISSING");
        assert!(missing_from_crc.is_missing());
        assert!(missing_from_string.is_missing());
    }

    #[test]
    fn test_registry() {
        let registry = HasherRegistry::new();

        assert!(registry.has_hasher("crc32"));
        assert!(registry.has_hasher("sha256"));
        assert!(registry.has_hasher("md5"));
        assert!(registry.has_hasher("sequence"));
        assert!(registry.has_hasher("hamming"));
        assert!(!registry.has_hasher("nonexistent"));

        let hashers = registry.list_hashers();
        assert_eq!(hashers.len(), 5);

        let names = registry.get_hasher_names();
        assert!(names.contains(&"crc32"));
        assert!(names.contains(&"sha256"));
    }

    #[test]
    fn test_parse_allele() {
        let registry = HasherRegistry::new();

        // Test CRC32 hasher
        let crc_hasher = registry.get_hasher("crc32").unwrap();
        let result = crc_hasher.parse_allele("12345", "-").unwrap();
        assert_eq!(result, AlleleHash::Crc32(12345));

        let missing = crc_hasher.parse_allele("-", "-").unwrap();
        assert!(missing.is_missing());

        // Test SHA256 hasher
        let sha_hasher = registry.get_hasher("sha256").unwrap();
        let seq_result = sha_hasher.parse_allele("ATCG", "-").unwrap();
        match seq_result {
            AlleleHash::String(_) => {} // Should be a string hash
            _ => panic!("Expected string hash"),
        }
    }
}
