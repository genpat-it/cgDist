// Phi test for recombination detection (Bruen et al. 2006)
// "A simple and robust statistical test for detecting the presence of recombination" - Genetics

use super::*;
use crate::core::recombination::factory::RecombinationDetectorFactory;

/// Phi test recombination detector
/// Implements the Phi statistic as described in Bruen et al. (2006)
pub struct PhiTestDetector {
    significance_level: f64,  // p-value threshold (default 0.05)
    min_informative_sites: usize,  // Minimum sites needed for reliable test
    bidirectional: bool,
}

impl PhiTestDetector {
    pub fn new(config: &str) -> Result<Self, String> {
        let params = RecombinationDetectorFactory::parse_config(config);
        
        let significance_level = params.get("significance")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.05);
            
        let min_informative_sites = params.get("min_sites")
            .and_then(|s| s.parse().ok())  
            .unwrap_or(4);
            
        let bidirectional = params.get("bidirectional")
            .and_then(|s| s.parse().ok())
            .unwrap_or(true);
        
        if significance_level <= 0.0 || significance_level >= 1.0 {
            return Err("Significance level must be between 0 and 1".to_string());
        }
        
        Ok(Self { significance_level, min_informative_sites, bidirectional })
    }
    
    /// Calculate Phi statistic for two sequences
    /// Phi = (n11 * n00 - n10 * n01)^2 * n / ((n11 + n10)(n01 + n00)(n11 + n01)(n10 + n00))
    fn calculate_phi_statistic(&self, seq1: &[u8], seq2: &[u8]) -> (f64, usize) {
        let mut n00 = 0; // Both sequences have same nucleotide A
        let mut n01 = 0; // seq1 has A, seq2 has B  
        let mut n10 = 0; // seq1 has B, seq2 has A
        let mut n11 = 0; // Both sequences have same nucleotide B (but different from A)
        let mut informative_sites = 0;
        
        let min_len = seq1.len().min(seq2.len());
        
        // Collect all polymorphic sites first
        let mut polymorphic_sites = Vec::new();
        for i in 0..min_len {
            if seq1[i] != seq2[i] && seq1[i] != b'N' && seq2[i] != b'N' {
                polymorphic_sites.push((seq1[i], seq2[i]));
            }
        }
        
        if polymorphic_sites.len() < 2 {
            return (0.0, 0); // Need at least 2 polymorphic sites
        }
        
        // For each pair of polymorphic sites, calculate contingency table
        for i in 0..polymorphic_sites.len() {
            for j in (i+1)..polymorphic_sites.len() {
                let (a1, b1) = polymorphic_sites[i];
                let (a2, b2) = polymorphic_sites[j];
                
                // Create 2x2 contingency table
                // Site i: A/B, Site j: A/B
                match ((a1, a2), (b1, b2)) {
                    ((a, c), (b, d)) if a == b && c == d => n00 += 1, // Same-Same
                    ((a, c), (b, d)) if a == b && c != d => n01 += 1, // Same-Different  
                    ((a, c), (b, d)) if a != b && c == d => n10 += 1, // Different-Same
                    ((a, c), (b, d)) if a != b && c != d => n11 += 1, // Different-Different
                    _ => {}
                }
                informative_sites += 1;
            }
        }
        
        let n = (n00 + n01 + n10 + n11) as f64;
        if n == 0.0 {
            return (0.0, 0);
        }
        
        // Calculate Phi statistic
        let numerator = (n11 as f64 * n00 as f64 - n10 as f64 * n01 as f64).powi(2) * n;
        let denominator = (n11 + n10) as f64 * (n01 + n00) as f64 * (n11 + n01) as f64 * (n10 + n00) as f64;
        
        if denominator == 0.0 {
            return (0.0, informative_sites);
        }
        
        let phi = numerator / denominator;
        (phi, informative_sites)
    }
    
    /// Calculate p-value from Phi statistic using chi-square distribution
    /// Phi follows chi-square with 1 degree of freedom under null hypothesis
    fn phi_to_pvalue(&self, phi: f64) -> f64 {
        if phi <= 0.0 {
            return 1.0;
        }
        
        // Simplified chi-square p-value calculation for df=1
        // For more accuracy, would use proper gamma function
        // This is approximation for phi test
        let chi_square = phi;
        
        // Rough approximation for chi-square(1) p-value
        // For df=1: p ≈ 2 * (1 - Φ(√χ²)) where Φ is standard normal CDF
        let z = chi_square.sqrt();
        
        // Standard normal CDF approximation
        // p ≈ 2 * (1 - 0.5 * (1 + erf(z/√2)))
        let erf_approx = if z > 3.0 { 1.0 } else { z / (1.0 + 0.278393 * z + 0.230389 * z * z + 0.000972 * z * z * z + 0.078108 * z * z * z * z) };
        let cdf = 0.5 * (1.0 + erf_approx);
        let p_value = 2.0 * (1.0 - cdf);
        
        p_value.max(0.0).min(1.0)
    }
}

impl RecombinationDetector for PhiTestDetector {
    fn name(&self) -> &'static str {
        "phi-test"
    }
    
    fn description(&self) -> &'static str {
        "EXPERIMENTAL: Phi statistical test for recombination (Bruen et al. 2006) - not validated for production use"
    }
    
    fn detect(
        &self, 
        seq1: &[u8], 
        seq2: &[u8], 
        _locus: &str,
        distance_forward: usize,
        distance_reverse_complement: Option<usize>
    ) -> RecombinationResult {
        // Calculate Phi statistic for forward alignment
        let (phi_forward, sites_forward) = self.calculate_phi_statistic(seq1, seq2);
        let p_value_forward = self.phi_to_pvalue(phi_forward);
        
        let (best_phi, best_p_value, best_distance) = if self.bidirectional && distance_reverse_complement.is_some() {
            // For bidirectional, we need reverse-complement sequence
            let seq2_revcomp: Vec<u8> = seq2.iter()
                .rev()
                .map(|&c| match c.to_ascii_uppercase() {
                    b'A' => b'T',
                    b'T' => b'A', 
                    b'G' => b'C',
                    b'C' => b'G',
                    b'N' => b'N',
                    _ => c,
                })
                .collect();
                
            let (phi_reverse, _sites_reverse) = self.calculate_phi_statistic(seq1, &seq2_revcomp);
            let p_value_reverse = self.phi_to_pvalue(phi_reverse);
            
            // Use the more significant result (lower p-value)
            if p_value_forward <= p_value_reverse {
                (phi_forward, p_value_forward, distance_forward)
            } else {
                (phi_reverse, p_value_reverse, distance_reverse_complement.unwrap())
            }
        } else {
            (phi_forward, p_value_forward, distance_forward)
        };
        
        let is_recombination = sites_forward >= self.min_informative_sites && 
                              best_p_value < self.significance_level;
        
        let confidence_score = if sites_forward >= self.min_informative_sites {
            1.0 - best_p_value  // Higher confidence = lower p-value
        } else {
            0.0  // Not enough sites for reliable test
        };
        
        let mut method_data = HashMap::new();
        method_data.insert("phi_statistic".to_string(), format!("{:.6}", best_phi));
        method_data.insert("p_value".to_string(), format!("{:.6}", best_p_value));
        method_data.insert("significance_level".to_string(), self.significance_level.to_string());
        method_data.insert("informative_sites".to_string(), sites_forward.to_string());
        method_data.insert("min_required_sites".to_string(), self.min_informative_sites.to_string());
        
        RecombinationResult {
            is_recombination,
            confidence_score,
            method_specific_data: method_data,
            statistical_significance: Some(best_p_value),
            distance_forward: Some(distance_forward),
            distance_reverse_complement,
            best_distance,
        }
    }
    
    fn get_config(&self) -> RecombinationDetectorConfig {
        let mut params = HashMap::new();
        params.insert("significance".to_string(), self.significance_level.to_string());
        params.insert("min_sites".to_string(), self.min_informative_sites.to_string());
        params.insert("bidirectional".to_string(), self.bidirectional.to_string());
        
        RecombinationDetectorConfig {
            detector_type: "phi-test".to_string(),
            parameters: params,
            bidirectional_alignment: self.bidirectional,
        }
    }
    
    fn requires_bidirectional(&self) -> bool {
        self.bidirectional
    }
    
    fn validate_config(config_str: &str) -> Result<(), String> {
        let params = RecombinationDetectorFactory::parse_config(config_str);
        
        if let Some(sig_str) = params.get("significance") {
            let sig: f64 = sig_str.parse()
                .map_err(|_| "Significance level must be a number")?;
            if sig <= 0.0 || sig >= 1.0 {
                return Err("Significance level must be between 0 and 1".to_string());
            }
        }
        
        if let Some(sites_str) = params.get("min_sites") {
            let sites: usize = sites_str.parse()
                .map_err(|_| "Minimum sites must be a positive integer")?;
            if sites < 2 {
                return Err("Minimum sites must be at least 2".to_string());
            }
        }
        
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_phi_test_detector() {
        let detector = PhiTestDetector::new("significance=0.05,min_sites=4").unwrap();
        assert_eq!(detector.significance_level, 0.05);
        assert_eq!(detector.min_informative_sites, 4);
        
        // Test with sequences that should show recombination
        let seq1 = b"ATCGATCGATCG";  
        let seq2 = b"ATCGTTCGATCG"; // Different in middle
        let result = detector.detect(seq1, seq2, "test_locus", 2, None);
        
        assert!(result.statistical_significance.is_some());
        assert!(result.method_specific_data.contains_key("phi_statistic"));
        assert!(result.method_specific_data.contains_key("p_value"));
    }
}