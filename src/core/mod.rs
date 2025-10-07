// mod.rs - Core logic module

pub mod alignment;
pub mod distance;
// pub mod recombination; // Disabled - pluggable system not needed for now

// Re-export main types for convenience
pub use alignment::{AlignmentConfig, DetailedAlignment, DistanceMode, compute_alignment_stats};
pub use distance::{DistanceEngine, calculate_sample_distance, calculate_distance_matrix};
// pub use recombination::{
//     RecombinationDetector, RecombinationResult, RecombinationDetectorConfig,
//     RecombinationDetectorFactory, ThresholdDetector, PhiTestDetector, RMRatioDetector
// };