// mod.rs - Data structures module

pub mod loaders;
pub mod profile;
pub mod sequence;

// Re-export main types for convenience
pub use profile::{AllelicMatrix, AllelicProfile};
pub use sequence::{SequenceDatabase, SequenceInfo};
