// Factory for creating recombination detectors

use super::*;
use std::collections::HashMap;

pub struct RecombinationDetectorFactory;

impl RecombinationDetectorFactory {
    /// Create a recombination detector by name and configuration string
    pub fn create(detector_type: &str, config: &str) -> Result<Box<dyn RecombinationDetector>, String> {
        match detector_type {
            "threshold" => {
                ThresholdDetector::validate_config(config)?;
                Ok(Box::new(ThresholdDetector::new(config)?))
            },
            "phi-test" => {
                PhiTestDetector::validate_config(config)?;
                Ok(Box::new(PhiTestDetector::new(config)?))
            },
            "r-m-ratio" => {
                RMRatioDetector::validate_config(config)?;
                Ok(Box::new(RMRatioDetector::new(config)?))
            },
            _ => Err(format!("Unknown recombination detector: {}", detector_type))
        }
    }
    
    /// List all available detectors
    pub fn list_available() -> Vec<(&'static str, &'static str)> {
        vec![
            ("threshold", "Fixed SNP threshold (legacy method)"),
            ("phi-test", "EXPERIMENTAL: Phi statistical test (Bruen et al. 2006) - not validated for production use"),
            ("r-m-ratio", "EXPERIMENTAL: r/m ratio for bacterial populations (ClonalFrame approach) - not validated for production use"),
        ]
    }
    
    /// Get default detector (using threshold as most stable method)
    pub fn default() -> Result<Box<dyn RecombinationDetector>, String> {
        Self::create("threshold", "threshold=20")
    }
    
    /// Parse configuration string into HashMap
    pub fn parse_config(config_str: &str) -> HashMap<String, String> {
        let mut config = HashMap::new();
        
        if config_str.is_empty() {
            return config;
        }
        
        for pair in config_str.split(',') {
            let parts: Vec<&str> = pair.split('=').collect();
            if parts.len() == 2 {
                config.insert(parts[0].trim().to_string(), parts[1].trim().to_string());
            }
        }
        
        config
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_parse_config() {
        let config = RecombinationDetectorFactory::parse_config("threshold=5,bidirectional=true");
        assert_eq!(config.get("threshold"), Some(&"5".to_string()));
        assert_eq!(config.get("bidirectional"), Some(&"true".to_string()));
    }
    
    #[test]
    fn test_list_available() {
        let detectors = RecombinationDetectorFactory::list_available();
        assert!(detectors.len() >= 3);
        assert!(detectors.iter().any(|(name, _)| *name == "phi-test"));
    }
}