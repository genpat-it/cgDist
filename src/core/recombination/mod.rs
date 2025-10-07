// Pluggable recombination detection framework
// Similar architecture to hashers for extensibility

use std::collections::HashMap;
use serde::{Serialize, Deserialize};

/// Result of recombination detection analysis
#[derive(Debug, Clone)]
pub struct RecombinationResult {
    pub is_recombination: bool,
    pub confidence_score: f64,  // 0.0-1.0
    pub method_specific_data: HashMap<String, String>,
    pub statistical_significance: Option<f64>, // p-value if applicable
    pub distance_forward: Option<usize>,
    pub distance_reverse_complement: Option<usize>,
    pub best_distance: usize,
}

/// Configuration for specific recombination detector (for cache compatibility)
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct RecombinationDetectorConfig {
    pub detector_type: String,
    pub parameters: HashMap<String, String>,
    pub bidirectional_alignment: bool,
}

/// Trait for pluggable recombination detection methods
pub trait RecombinationDetector: Send + Sync {
    /// Get detector name
    fn name(&self) -> &'static str;
    
    /// Get detector description
    fn description(&self) -> &'static str;
    
    /// Detect recombination between two sequences
    fn detect(
        &self, 
        seq1: &[u8], 
        seq2: &[u8], 
        locus: &str,
        distance_forward: usize,
        distance_reverse_complement: Option<usize>
    ) -> RecombinationResult;
    
    /// Get detector-specific configuration for cache compatibility
    fn get_config(&self) -> RecombinationDetectorConfig;
    
    /// Check if this detector requires bidirectional alignment
    fn requires_bidirectional(&self) -> bool;
    
    /// Validate configuration parameters
    fn validate_config(config_str: &str) -> Result<(), String> where Self: Sized;
}

// Re-export detector implementations
pub mod threshold;
pub mod phi_test; 
pub mod rm_ratio;
pub mod factory;

pub use threshold::ThresholdDetector;
pub use phi_test::PhiTestDetector;
pub use rm_ratio::RMRatioDetector;
pub use factory::RecombinationDetectorFactory;