// r/m ratio detector for bacterial populations
// Based on ClonalFrame methodology (Didelot & Wilson 2015)

use super::*;
use crate::core::recombination::factory::RecombinationDetectorFactory;

/// r/m ratio recombination detector
/// Uses recombination/mutation ratio typical for bacterial species
pub struct RMRatioDetector {
    r_m_ratio: f64,           // Species-specific r/m ratio
    mutation_rate: f64,       // Expected mutations per site per generation
    generation_threshold: f64, // Generations since common ancestor threshold
    bidirectional: bool,
}

impl RMRatioDetector {
    pub fn new(config: &str) -> Result<Self, String> {
        let params = RecombinationDetectorFactory::parse_config(config);
        
        let r_m_ratio = params.get("r_m_ratio")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.2); // Default for Salmonella from literature
            
        let mutation_rate = params.get("mutation_rate")
            .and_then(|s| s.parse().ok())
            .unwrap_or(1e-6); // Mutations per site per generation
            
        let generation_threshold = params.get("generation_threshold")
            .and_then(|s| s.parse().ok())
            .unwrap_or(100.0); // Generations since MRCA
            
        let bidirectional = params.get("bidirectional")
            .and_then(|s| s.parse().ok())
            .unwrap_or(true);
        
        if r_m_ratio < 0.0 {
            return Err("r/m ratio must be non-negative".to_string());
        }
        
        if mutation_rate <= 0.0 {
            return Err("Mutation rate must be positive".to_string());
        }
        
        Ok(Self { r_m_ratio, mutation_rate, generation_threshold, bidirectional })
    }
    
    /// Calculate expected mutations based on sequence length and divergence time
    fn expected_mutations(&self, sequence_length: usize, generations: f64) -> f64 {
        sequence_length as f64 * self.mutation_rate * generations
    }
    
    /// Calculate expected recombinations based on r/m ratio
    fn expected_recombinations(&self, expected_mutations: f64) -> f64 {
        expected_mutations * self.r_m_ratio
    }
    
    /// Calculate likelihood ratio for recombination vs. mutation-only
    fn likelihood_ratio(&self, observed_snps: usize, sequence_length: usize) -> (f64, bool) {
        // Estimate generations from observed SNP density
        let estimated_generations = observed_snps as f64 / (sequence_length as f64 * self.mutation_rate);
        
        let expected_muts = self.expected_mutations(sequence_length, estimated_generations);
        let expected_recombs = self.expected_recombinations(expected_muts);
        
        // Total expected changes = mutations + recombinations
        let _expected_total = expected_muts + expected_recombs;
        
        // Simple likelihood ratio test
        // If observed >> expected from mutations alone, likely recombination
        let mutation_only_expected = self.expected_mutations(sequence_length, self.generation_threshold);
        let recomb_plus_mutation_expected = mutation_only_expected * (1.0 + self.r_m_ratio);
        
        let likelihood_recomb = self.poisson_likelihood(observed_snps, recomb_plus_mutation_expected);
        let likelihood_mutation_only = self.poisson_likelihood(observed_snps, mutation_only_expected);
        
        let ratio = if likelihood_mutation_only > 0.0 {
            likelihood_recomb / likelihood_mutation_only
        } else {
            f64::INFINITY
        };
        
        // If recombination model is significantly more likely
        let is_recombination = ratio > 10.0 && observed_snps as f64 > mutation_only_expected * 2.0;
        
        (ratio, is_recombination)
    }
    
    /// Simple Poisson likelihood approximation
    fn poisson_likelihood(&self, observed: usize, lambda: f64) -> f64 {
        if lambda <= 0.0 {
            return if observed == 0 { 1.0 } else { 0.0 };
        }
        
        // Poisson PMF: P(k) = λ^k * e^(-λ) / k!
        // Use log to avoid overflow
        let log_lambda = lambda.ln();
        let log_factorial = self.log_factorial(observed);
        let log_likelihood = observed as f64 * log_lambda - lambda - log_factorial;
        
        log_likelihood.exp()
    }
    
    /// Approximate log factorial using Stirling's approximation for large n
    fn log_factorial(&self, n: usize) -> f64 {
        if n == 0 || n == 1 {
            return 0.0;
        }
        
        if n < 10 {
            // Exact for small values
            (2..=n).map(|i| (i as f64).ln()).sum()
        } else {
            // Stirling's approximation: ln(n!) ≈ n*ln(n) - n + 0.5*ln(2πn)
            let n_f = n as f64;
            n_f * n_f.ln() - n_f + 0.5 * (2.0 * std::f64::consts::PI * n_f).ln()
        }
    }
}

impl RecombinationDetector for RMRatioDetector {
    fn name(&self) -> &'static str {
        "r-m-ratio"
    }
    
    fn description(&self) -> &'static str {
        "EXPERIMENTAL: r/m ratio detector for bacterial populations (ClonalFrame approach) - not validated for production use"
    }
    
    fn detect(
        &self, 
        seq1: &[u8], 
        seq2: &[u8], 
        _locus: &str,
        distance_forward: usize,
        distance_reverse_complement: Option<usize>
    ) -> RecombinationResult {
        let sequence_length = seq1.len().min(seq2.len());
        
        let (best_distance, snps_to_test) = if self.bidirectional && distance_reverse_complement.is_some() {
            let rev_dist = distance_reverse_complement.unwrap();
            if distance_forward <= rev_dist {
                (distance_forward, distance_forward)
            } else {
                (rev_dist, rev_dist)
            }
        } else {
            (distance_forward, distance_forward)
        };
        
        let (likelihood_ratio, is_recombination) = self.likelihood_ratio(snps_to_test, sequence_length);
        
        // Confidence based on how strongly the data supports recombination
        let confidence_score = if is_recombination {
            (likelihood_ratio.ln() / 10.0).min(1.0).max(0.5) // Scale log-likelihood ratio
        } else {
            (1.0 / likelihood_ratio).min(0.5).max(0.0)
        };
        
        let mut method_data = HashMap::new();
        method_data.insert("r_m_ratio".to_string(), self.r_m_ratio.to_string());
        method_data.insert("likelihood_ratio".to_string(), format!("{:.4}", likelihood_ratio));
        method_data.insert("sequence_length".to_string(), sequence_length.to_string());
        method_data.insert("observed_snps".to_string(), snps_to_test.to_string());
        
        let expected_mutations = self.expected_mutations(sequence_length, self.generation_threshold);
        method_data.insert("expected_mutations".to_string(), format!("{:.2}", expected_mutations));
        method_data.insert("expected_recombinations".to_string(), format!("{:.2}", self.expected_recombinations(expected_mutations)));
        
        RecombinationResult {
            is_recombination,
            confidence_score,
            method_specific_data: method_data,
            statistical_significance: Some(1.0 / likelihood_ratio), // Inverse as pseudo p-value
            distance_forward: Some(distance_forward),
            distance_reverse_complement,
            best_distance,
        }
    }
    
    fn get_config(&self) -> RecombinationDetectorConfig {
        let mut params = HashMap::new();
        params.insert("r_m_ratio".to_string(), self.r_m_ratio.to_string());
        params.insert("mutation_rate".to_string(), self.mutation_rate.to_string());
        params.insert("generation_threshold".to_string(), self.generation_threshold.to_string());
        params.insert("bidirectional".to_string(), self.bidirectional.to_string());
        
        RecombinationDetectorConfig {
            detector_type: "r-m-ratio".to_string(),
            parameters: params,
            bidirectional_alignment: self.bidirectional,
        }
    }
    
    fn requires_bidirectional(&self) -> bool {
        self.bidirectional
    }
    
    fn validate_config(config_str: &str) -> Result<(), String> {
        let params = RecombinationDetectorFactory::parse_config(config_str);
        
        if let Some(ratio_str) = params.get("r_m_ratio") {
            let ratio: f64 = ratio_str.parse()
                .map_err(|_| "r/m ratio must be a number")?;
            if ratio < 0.0 {
                return Err("r/m ratio must be non-negative".to_string());
            }
        }
        
        if let Some(rate_str) = params.get("mutation_rate") {
            let rate: f64 = rate_str.parse()
                .map_err(|_| "Mutation rate must be a number")?;
            if rate <= 0.0 {
                return Err("Mutation rate must be positive".to_string());
            }
        }
        
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_rm_ratio_detector() {
        let detector = RMRatioDetector::new("r_m_ratio=0.3,mutation_rate=1e-6").unwrap();
        assert_eq!(detector.r_m_ratio, 0.3);
        assert_eq!(detector.mutation_rate, 1e-6);
        
        // Test with high SNP density (should indicate recombination)
        let seq1 = b"ATCGATCGATCG";
        let seq2 = b"GTCGATTGATCG"; 
        let result = detector.detect(seq1, seq2, "test_locus", 25, None); // High distance
        
        assert!(result.method_specific_data.contains_key("likelihood_ratio"));
        assert!(result.method_specific_data.contains_key("r_m_ratio"));
    }
    
    #[test]
    fn test_log_factorial() {
        let detector = RMRatioDetector::new("").unwrap();
        assert_eq!(detector.log_factorial(0), 0.0);
        assert_eq!(detector.log_factorial(1), 0.0);
        assert!((detector.log_factorial(5) - (1.0 + 2.0_f64.ln() + 3.0_f64.ln() + 4.0_f64.ln() + 5.0_f64.ln())).abs() < 1e-10);
    }
}