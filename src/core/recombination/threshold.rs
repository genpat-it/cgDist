// Legacy fixed threshold detector (current implementation)

use super::*;
use crate::core::recombination::factory::RecombinationDetectorFactory;

/// Fixed threshold recombination detector
/// Legacy method: distance > threshold = recombination
pub struct ThresholdDetector {
    threshold: usize,
    threshold_percent: Option<f64>,  // New: percentage-based threshold
    bidirectional: bool,
}

impl ThresholdDetector {
    pub fn new(config: &str) -> Result<Self, String> {
        let params = RecombinationDetectorFactory::parse_config(config);
        
        let threshold = params.get("threshold")
            .and_then(|s| s.parse().ok())
            .unwrap_or(5);
            
        let threshold_percent = params.get("threshold_percent")
            .and_then(|s| s.parse().ok());
            
        let bidirectional = params.get("bidirectional")
            .and_then(|s| s.parse().ok())
            .unwrap_or(true);
        
        // Validate that user doesn't specify both absolute and percentage
        if threshold != 5 && threshold_percent.is_some() {
            return Err("Cannot specify both 'threshold' and 'threshold_percent'. Use one or the other.".to_string());
        }
        
        // Validate percentage range
        if let Some(percent) = threshold_percent {
            if percent <= 0.0 || percent >= 100.0 {
                return Err("threshold_percent must be between 0.0 and 100.0".to_string());
            }
        }
        
        Ok(Self { threshold, threshold_percent, bidirectional })
    }
}

impl RecombinationDetector for ThresholdDetector {
    fn name(&self) -> &'static str {
        "threshold"
    }
    
    fn description(&self) -> &'static str {
        "Fixed SNP threshold detector (legacy method)"
    }
    
    fn detect(
        &self, 
        seq1: &[u8], 
        seq2: &[u8], 
        _locus: &str,
        distance_forward: usize,
        distance_reverse_complement: Option<usize>
    ) -> RecombinationResult {
        let best_distance = if self.bidirectional && distance_reverse_complement.is_some() {
            distance_forward.min(distance_reverse_complement.unwrap())
        } else {
            distance_forward
        };
        
        // Calculate average sequence length for percentage-based threshold
        let seq_len1 = seq1.len();
        let seq_len2 = seq2.len();
        let avg_length = if seq_len1 > 0 && seq_len2 > 0 {
            (seq_len1 + seq_len2) as f64 / 2.0
        } else {
            450.0 // Default MLST locus length
        };
        
        // Determine if recombination based on threshold type
        let (is_recombination, effective_threshold) = if let Some(percent_threshold) = self.threshold_percent {
            // Percentage-based threshold
            let absolute_threshold_from_percent = (avg_length * percent_threshold / 100.0).ceil() as usize;
            let is_recomb = best_distance > absolute_threshold_from_percent;
            (is_recomb, absolute_threshold_from_percent)
        } else {
            // Absolute threshold (legacy)
            let is_recomb = best_distance > self.threshold;
            (is_recomb, self.threshold)
        };
        
        let mut method_data = HashMap::new();
        if let Some(percent_threshold) = self.threshold_percent {
            method_data.insert("threshold_percent".to_string(), format!("{:.2}", percent_threshold));
            method_data.insert("effective_threshold".to_string(), effective_threshold.to_string());
            method_data.insert("avg_sequence_length".to_string(), format!("{:.1}", avg_length));
        } else {
            method_data.insert("threshold".to_string(), self.threshold.to_string());
        }
        method_data.insert("best_distance".to_string(), best_distance.to_string());
        
        if self.bidirectional {
            method_data.insert("distance_forward".to_string(), distance_forward.to_string());
            if let Some(rev_dist) = distance_reverse_complement {
                method_data.insert("distance_reverse_complement".to_string(), rev_dist.to_string());
            }
        }
        
        RecombinationResult {
            is_recombination,
            confidence_score: if is_recombination { 1.0 } else { 0.0 }, // Binary confidence
            method_specific_data: method_data,
            statistical_significance: None, // No p-value for threshold method
            distance_forward: Some(distance_forward),
            distance_reverse_complement,
            best_distance,
        }
    }
    
    fn get_config(&self) -> RecombinationDetectorConfig {
        let mut params = HashMap::new();
        if let Some(percent_threshold) = self.threshold_percent {
            params.insert("threshold_percent".to_string(), percent_threshold.to_string());
        } else {
            params.insert("threshold".to_string(), self.threshold.to_string());
        }
        params.insert("bidirectional".to_string(), self.bidirectional.to_string());
        
        RecombinationDetectorConfig {
            detector_type: "threshold".to_string(),
            parameters: params,
            bidirectional_alignment: self.bidirectional,
        }
    }
    
    fn requires_bidirectional(&self) -> bool {
        self.bidirectional
    }
    
    fn validate_config(config_str: &str) -> Result<(), String> {
        let params = RecombinationDetectorFactory::parse_config(config_str);
        
        let has_threshold = params.get("threshold").is_some();
        let has_threshold_percent = params.get("threshold_percent").is_some();
        
        // Check that user doesn't specify both
        if has_threshold && has_threshold_percent {
            return Err("Cannot specify both 'threshold' and 'threshold_percent'. Use one or the other.".to_string());
        }
        
        if let Some(threshold_str) = params.get("threshold") {
            if threshold_str.parse::<usize>().is_err() {
                return Err("Threshold must be a positive integer".to_string());
            }
        }
        
        if let Some(threshold_percent_str) = params.get("threshold_percent") {
            let percent: f64 = threshold_percent_str.parse()
                .map_err(|_| "threshold_percent must be a number")?;
            if percent <= 0.0 || percent >= 100.0 {
                return Err("threshold_percent must be between 0.0 and 100.0".to_string());
            }
        }
        
        if let Some(bidirectional_str) = params.get("bidirectional") {
            if bidirectional_str.parse::<bool>().is_err() {
                return Err("Bidirectional must be true or false".to_string());
            }
        }
        
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_threshold_detector_absolute() {
        let detector = ThresholdDetector::new("threshold=10,bidirectional=true").unwrap();
        assert_eq!(detector.threshold, 10);
        assert_eq!(detector.bidirectional, true);
        assert!(detector.threshold_percent.is_none());
        
        // Test detection
        let seq1 = b"ATCG";
        let seq2 = b"GTCA";  
        let result = detector.detect(seq1, seq2, "test_locus", 15, Some(8));
        
        // 15 > 10 threshold, but min(15,8)=8 < 10, so no recombination
        assert_eq!(result.is_recombination, false);
        assert_eq!(result.best_distance, 8);
    }
    
    #[test]
    fn test_threshold_detector_percentage() {
        let detector = ThresholdDetector::new("threshold_percent=5.0,bidirectional=true").unwrap();
        assert_eq!(detector.threshold, 5); // default when using percentage
        assert_eq!(detector.threshold_percent, Some(5.0));
        assert_eq!(detector.bidirectional, true);
        
        // Test detection with sequences of 100bp each (avg=100)
        // 5% of 100 = 5 events threshold
        let seq1 = &[b'A'; 100];
        let seq2 = &[b'T'; 100];
        let result = detector.detect(seq1, seq2, "test_locus", 6, Some(4)); // min(6,4) = 4
        
        // 4 events < 5 (5% of 100), so no recombination
        assert_eq!(result.is_recombination, false);
        assert_eq!(result.best_distance, 4);
        
        let result2 = detector.detect(seq1, seq2, "test_locus", 7, Some(6)); // min(7,6) = 6
        // 6 events > 5 (5% of 100), so recombination detected
        assert_eq!(result2.is_recombination, true);
        assert_eq!(result2.best_distance, 6);
    }
    
    #[test]
    fn test_threshold_detector_validation() {
        // Should fail if both threshold and threshold_percent specified
        assert!(ThresholdDetector::new("threshold=10,threshold_percent=5.0").is_err());
        
        // Should fail if threshold_percent out of range
        assert!(ThresholdDetector::new("threshold_percent=0.0").is_err());
        assert!(ThresholdDetector::new("threshold_percent=100.0").is_err());
        assert!(ThresholdDetector::new("threshold_percent=150.0").is_err());
        
        // Should succeed with valid values
        assert!(ThresholdDetector::new("threshold_percent=2.5").is_ok());
        assert!(ThresholdDetector::new("threshold=15").is_ok());
    }
}