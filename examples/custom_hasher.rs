#!/usr/bin/env rust-script

//! # Custom Hasher Example
//!
//! This example demonstrates how to implement and use custom hashers in cgdist-rs.
//!
//! Usage:
//! ```bash
//! cargo run --example custom_hasher
//! ```

use cgdist::hashers::{AlleleHash, AlleleHasher, HasherRegistry};

/// Example 1: Simple nucleotide composition hasher
/// Counts A, T, G, C nucleotides and creates a hash from their counts
#[derive(Debug)]
pub struct CompositionHasher;

impl AlleleHasher for CompositionHasher {
    fn hash_sequence(&self, sequence: &str) -> AlleleHash {
        // Count nucleotides: A, T, G, C
        let mut counts = [0u32; 4]; // A, T, G, C
        for nucleotide in sequence.chars() {
            match nucleotide.to_ascii_uppercase() {
                'A' => counts[0] += 1,
                'T' => counts[1] += 1,
                'G' => counts[2] += 1,
                'C' => counts[3] += 1,
                _ => {} // Ignore ambiguous bases
            }
        }

        // Create hash from composition: A2T3G1C4 format
        let hash_string = format!("A{}T{}G{}C{}", counts[0], counts[1], counts[2], counts[3]);
        AlleleHash::String(hash_string)
    }

    fn parse_allele(&self, allele_str: &str, missing_char: &str) -> Result<AlleleHash, String> {
        if allele_str == missing_char {
            Ok(AlleleHash::Missing)
        } else {
            // Parse composition string or return as-is
            Ok(AlleleHash::String(allele_str.to_string()))
        }
    }

    fn name(&self) -> &'static str {
        "composition"
    }

    fn description(&self) -> &'static str {
        "Nucleotide composition-based hasher (A/T/G/C counts)"
    }

    fn validate_sequence(&self, sequence: &str) -> Result<(), String> {
        // Only allow ATGC nucleotides
        for ch in sequence.chars() {
            match ch.to_ascii_uppercase() {
                'A' | 'T' | 'G' | 'C' | 'N' => {}
                _ => return Err(format!("Invalid nucleotide: {}", ch)),
            }
        }
        Ok(())
    }
}

/// Example 2: K-mer based hasher
/// Creates hash from sorted k-mers in the sequence
#[derive(Debug)]
pub struct KmerHasher {
    k: usize,
}

impl KmerHasher {
    pub fn new(k: usize) -> Self {
        Self { k }
    }
}

impl AlleleHasher for KmerHasher {
    fn hash_sequence(&self, sequence: &str) -> AlleleHash {
        let mut kmers = Vec::new();
        let seq_bytes = sequence.as_bytes();

        if seq_bytes.len() >= self.k {
            for i in 0..=(seq_bytes.len() - self.k) {
                let kmer = std::str::from_utf8(&seq_bytes[i..i + self.k])
                    .unwrap_or("")
                    .to_string();
                kmers.push(kmer);
            }
        }

        kmers.sort();
        let hash_string = kmers.join("|");
        AlleleHash::String(hash_string)
    }

    fn parse_allele(&self, allele_str: &str, missing_char: &str) -> Result<AlleleHash, String> {
        if allele_str == missing_char {
            Ok(AlleleHash::Missing)
        } else {
            Ok(AlleleHash::String(allele_str.to_string()))
        }
    }

    fn name(&self) -> &'static str {
        "kmer"
    }

    fn description(&self) -> &'static str {
        "K-mer composition hasher (sorted k-mers)"
    }
}

/// Example 3: Custom numeric hasher using polynomial rolling hash
#[derive(Debug)]
pub struct PolynomialHasher;

impl AlleleHasher for PolynomialHasher {
    fn hash_sequence(&self, sequence: &str) -> AlleleHash {
        // Convert sequence to custom numeric representation
        let mut hash_value = 0u32;
        let base = 4u32;

        for nucleotide in sequence.chars() {
            let base_value = match nucleotide.to_ascii_uppercase() {
                'A' => 0,
                'T' => 1,
                'G' => 2,
                'C' => 3,
                _ => 0, // Default for ambiguous
            };
            // Polynomial rolling hash with overflow handling
            hash_value = hash_value.wrapping_mul(base).wrapping_add(base_value);
        }

        AlleleHash::Crc32(hash_value)
    }

    fn parse_allele(&self, allele_str: &str, missing_char: &str) -> Result<AlleleHash, String> {
        if allele_str == missing_char {
            Ok(AlleleHash::Missing)
        } else {
            match allele_str.parse::<u32>() {
                Ok(value) => Ok(AlleleHash::Crc32(value)),
                Err(_) => Err(format!("Invalid numeric allele: {}", allele_str)),
            }
        }
    }

    fn name(&self) -> &'static str {
        "polynomial"
    }

    fn description(&self) -> &'static str {
        "Polynomial rolling hash for sequences"
    }
}

fn main() {
    println!("🔌 cgdist-rs Custom Hasher Examples");
    println!("===================================\n");

    // Create registry with custom hashers
    let mut registry = HasherRegistry::new();

    // Register custom hashers
    registry.register_hasher("composition", Box::new(CompositionHasher));
    registry.register_hasher("kmer3", Box::new(KmerHasher::new(3)));
    registry.register_hasher("polynomial", Box::new(PolynomialHasher));

    // Test sequences
    let test_sequences = vec![
        "ATCGATCGATCG",
        "AAATTTGGGCCC",
        "ATGCATGCATGC",
        "NNNATCGATCGN",
    ];

    println!("📊 Available Hashers:");
    for (name, description) in registry.list_hashers() {
        println!("  • {}: {}", name, description);
    }
    println!();

    // Test each hasher with each sequence
    for (hasher_name, _) in registry.list_hashers() {
        if let Some(hasher) = registry.get_hasher(hasher_name) {
            println!("🧬 Testing hasher: {}", hasher.name());
            println!("   Description: {}", hasher.description());

            for sequence in &test_sequences {
                // Validate sequence
                match hasher.validate_sequence(sequence) {
                    Ok(()) => {
                        let hash = hasher.hash_sequence(sequence);
                        println!("   {} → {}", sequence, hash);
                    }
                    Err(e) => {
                        println!("   {} → ERROR: {}", sequence, e);
                    }
                }
            }
            println!();
        }
    }

    // Demonstrate hash consistency
    println!("🔄 Testing Hash Consistency:");
    let test_seq = "ATCGATCG";
    if let Some(hasher) = registry.get_hasher("composition") {
        let hash1 = hasher.hash_sequence(test_seq);
        let hash2 = hasher.hash_sequence(test_seq);
        println!("   Sequence: {}", test_seq);
        println!("   Hash 1: {}", hash1);
        println!("   Hash 2: {}", hash2);
        println!("   Consistent: {}", hash1 == hash2);
    }
    println!();

    // Demonstrate missing value handling
    println!("🚫 Testing Missing Value Handling:");
    if let Some(hasher) = registry.get_hasher("composition") {
        match hasher.parse_allele("-", "-") {
            Ok(AlleleHash::Missing) => println!("   Missing marker correctly parsed"),
            Ok(other) => println!("   Unexpected result: {}", other),
            Err(e) => println!("   Error: {}", e),
        }

        match hasher.parse_allele("A5T3G2C1", "-") {
            Ok(hash) => println!("   Valid allele parsed: {}", hash),
            Err(e) => println!("   Error: {}", e),
        }
    }

    println!("\n✅ Custom hasher examples completed!");
    println!("💡 Tip: Use these patterns to implement your own specialized hashers");
}
