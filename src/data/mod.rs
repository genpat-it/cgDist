// mod.rs - Data structures module

pub mod profile;
pub mod sequence;
pub mod loaders;

// Re-export main types for convenience
pub use profile::{AllelicProfile, AllelicMatrix};
pub use sequence::{SequenceInfo, SequenceDatabase};