// mod.rs - Core logic module

pub mod alignment;
pub mod distance;
// pub mod recombination; // Disabled - pluggable system not needed for now

// Re-export main types for convenience
pub use alignment::{compute_alignment_stats, AlignmentConfig, DetailedAlignment, DistanceMode};
pub use distance::{calculate_distance_matrix, calculate_sample_distance, DistanceEngine};
// pub use recombination::{
//     RecombinationDetector, RecombinationResult, RecombinationDetectorConfig,
//     RecombinationDetectorFactory, ThresholdDetector, PhiTestDetector, RMRatioDetector
// };
