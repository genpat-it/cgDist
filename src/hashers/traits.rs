// traits.rs - Core traits and types for the hasher system

use serde::{Deserialize, Serialize};
use std::fmt::{Debug, Display};

/// Unified enum for all supported hash types - production ready
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum AlleleHash {
    /// CRC32 hash (u32) - compatible with chewBACCA
    Crc32(u32),
    /// String-based hash (SHA256, MD5, etc.)
    String(String),
    /// Missing allele marker
    Missing,
}

impl AlleleHash {
    pub fn is_missing(&self) -> bool {
        matches!(self, AlleleHash::Missing)
    }

    pub fn as_crc32(&self) -> Option<u32> {
        match self {
            AlleleHash::Crc32(val) => Some(*val),
            _ => None,
        }
    }

    pub fn as_string(&self) -> Option<&str> {
        match self {
            AlleleHash::String(s) => Some(s),
            _ => None,
        }
    }

    /// Create from CRC32 value
    pub fn from_crc32(val: u32) -> Self {
        if val == u32::MAX {
            AlleleHash::Missing
        } else {
            AlleleHash::Crc32(val)
        }
    }

    /// Create from string value
    pub fn from_string(val: String, missing_marker: &str) -> Self {
        if val == missing_marker {
            AlleleHash::Missing
        } else {
            AlleleHash::String(val)
        }
    }
}

impl Display for AlleleHash {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AlleleHash::Crc32(val) => write!(f, "{}", val),
            AlleleHash::String(s) => write!(f, "{}", s),
            AlleleHash::Missing => write!(f, "MISSING"),
        }
    }
}

/// Pair of allele hashes for distance calculations
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct AlleleHashPair {
    pub hash1: AlleleHash,
    pub hash2: AlleleHash,
}

impl AlleleHashPair {
    pub fn new(hash1: AlleleHash, hash2: AlleleHash) -> Self {
        if hash1 <= hash2 {
            Self { hash1, hash2 }
        } else {
            Self {
                hash1: hash2,
                hash2: hash1,
            }
        }
    }

    /// Check if either hash is missing
    pub fn has_missing(&self) -> bool {
        self.hash1.is_missing() || self.hash2.is_missing()
    }

    /// Convert to legacy CrcPair if both are CRC32 (used for backward compatibility)
    pub fn to_crc_pair(&self) -> Option<(u32, u32)> {
        match (&self.hash1, &self.hash2) {
            (AlleleHash::Crc32(c1), AlleleHash::Crc32(c2)) => Some((*c1, *c2)),
            _ => None,
        }
    }
}

/// Trait for allele identification strategies - Production Ready
/// This allows pluggable mechanisms for how alleles are identified/hashed
pub trait AlleleHasher: Send + Sync + Debug {
    /// Compute the hash/identifier for a nucleotide sequence
    fn hash_sequence(&self, sequence: &str) -> AlleleHash;

    /// Parse an allele string from a profile file
    fn parse_allele(&self, allele_str: &str, missing_char: &str) -> Result<AlleleHash, String>;

    /// Get a human-readable name for this hasher
    fn name(&self) -> &'static str;

    /// Get a description of this hasher
    fn description(&self) -> &'static str;

    /// Create a pair of hashes with ordering
    fn make_pair(&self, hash1: AlleleHash, hash2: AlleleHash) -> AlleleHashPair {
        AlleleHashPair::new(hash1, hash2)
    }

    /// Validate that a sequence is acceptable for this hasher
    fn validate_sequence(&self, _sequence: &str) -> Result<(), String> {
        // Default implementation accepts any sequence
        Ok(())
    }

    /// Get the missing allele marker
    fn missing_allele(&self) -> AlleleHash {
        AlleleHash::Missing
    }
}
