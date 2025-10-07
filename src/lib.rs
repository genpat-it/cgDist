// lib.rs - cgdist library root

//! # cgdist - High-performance SNP/indel-level distance calculator for core genome MLST analysis
//!
//! This library provides a high-performance implementation for calculating genetic distances
//! between bacterial samples using core genome MLST (cgMLST) data. It supports multiple
//! hashing algorithms and is compatible with chewBACCA allele calling.
//!
//! ## Features
//!
//! - **High performance**: Optimized parallel processing and caching
//! - **Plugin system**: Support for CRC32, SHA256, MD5, and custom hashers
//! - **Multiple formats**: TSV, CSV, PHYLIP, NEXUS output formats
//! - **Flexible filtering**: Sample and loci filtering with regex and file lists
//! - **Quality control**: Configurable thresholds for data completeness
//! - **chewBACCA compatible**: Full backward compatibility with existing workflows
//!
//! ## Basic Usage
//!
//! ```rust,no_run
//! use cgdist::prelude::*;
//!
//! // Load allelic profiles with CRC32 hasher (chewBACCA compatible)
//! let matrix = AllelicMatrix::from_file_with_hasher(
//!     std::path::Path::new("profiles.tsv"),
//!     "-",  // missing character
//!     "crc32",  // hasher type
//!     0.0,  // sample threshold
//!     0.0,  // locus threshold
//!     None, None, None, None,  // filters
//!     None, None, None, None,
//! )?;
//!
//! // Calculate distances
//! let engine = DistanceEngine::new(AlignmentConfig::default());
//! let distances = calculate_distance_matrix(
//!     &matrix.samples,
//!     &matrix.loci_names,
//!     &engine,
//!     DistanceMode::SnpsOnly,
//!     0,  // min loci
//!     false,  // hamming fallback
//! );
//! # Ok::<(), String>(())
//! ```

// Re-export all main modules
pub mod hashers;
pub mod data;
pub mod core;
pub mod cli;
pub mod output;

// Convenience prelude for common imports
pub mod prelude {
    pub use crate::hashers::{AlleleHash, AlleleHashPair, AlleleHasher, HasherRegistry};
    pub use crate::hashers::{Crc32Hasher, Sha256Hasher, Md5Hasher, SequenceHasher};
    pub use crate::data::{AllelicProfile, AllelicMatrix, SequenceInfo, SequenceDatabase};
    pub use crate::core::{AlignmentConfig, DistanceMode, DistanceEngine};
    pub use crate::core::{calculate_distance_matrix, calculate_sample_distance};
    pub use crate::cli::{Args, ValidationResult, validate_args};
    pub use crate::output::write_matrix;
}

// Re-export main types at the root level for convenience
pub use hashers::{AlleleHash, AlleleHashPair, AlleleHasher, HasherRegistry};
pub use data::{AllelicProfile, AllelicMatrix, SequenceInfo, SequenceDatabase};
pub use core::{AlignmentConfig, DistanceMode, DistanceEngine};
pub use cli::{Args, ValidationResult};

/// Library version
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Get library information
pub fn get_info() -> String {
    format!("cgdist v{} - High-performance distance calculator for cgMLST", VERSION)
}