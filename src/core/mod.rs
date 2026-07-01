// mod.rs - Core logic module

pub mod alignment;
pub mod distance;
// pub mod recombination; // Disabled - pluggable system not needed for now

// Re-export main types for convenience
pub use alignment::{compute_alignment_stats, AlignmentConfig, DetailedAlignment, DistanceMode};
pub use distance::{
    calculate_distance_matrix, calculate_pairs_recombination, calculate_pairs_table,
    calculate_sample_distance, calculate_sample_distance_detailed, DistanceEngine, PairRow,
};
// pub use recombination::{
//     RecombinationDetector, RecombinationResult, RecombinationDetectorConfig,
//     RecombinationDetectorFactory, ThresholdDetector, PhiTestDetector, RMRatioDetector
// };
