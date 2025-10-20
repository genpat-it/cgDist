use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn load_matrix_samples(file_path: &str) -> Result<HashSet<String>, Box<dyn std::error::Error>> {
    println!("📋 Loading sample names from {}...", file_path);
    
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Skip comments
    for line in lines.by_ref() {
        let line = line?;
        if !line.starts_with('#') {
            // This is the header with sample names
            let samples: Vec<&str> = line.split('\t').collect();
            let sample_set: HashSet<String> = samples.iter()
                .skip(1) // Skip the first "Sample" column
                .map(|s| s.to_string())
                .collect();
            
            println!("✅ Found {} samples in matrix", sample_set.len());
            return Ok(sample_set);
        }
    }
    
    Err("No header found in matrix file".into())
}

fn load_recombination_samples(file_path: &str) -> Result<HashSet<String>, Box<dyn std::error::Error>> {
    println!("📋 Loading sample names from recombination file {}...", file_path);
    
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    
    // Skip header
    lines.next();
    
    let mut samples = HashSet::new();
    
    for line in lines.take(1000) { // Check only first 1000 lines for efficiency
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        
        if parts.len() >= 2 {
            samples.insert(parts[0].to_string()); // Sample1
            samples.insert(parts[1].to_string()); // Sample2
        }
    }
    
    println!("✅ Found {} unique samples in recombination file (from first 1000 lines)", samples.len());
    Ok(samples)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    
    if args.len() != 3 {
        eprintln!("Usage: {} <recombination.tsv> <matrix.tsv>", args[0]);
        std::process::exit(1);
    }
    
    let recombination_file = &args[1];
    let matrix_file = &args[2];
    
    println!("🔍 Checking sample overlap between recombination file and matrix");
    
    let matrix_samples = load_matrix_samples(matrix_file)?;
    let recomb_samples = load_recombination_samples(recombination_file)?;
    
    // Find intersection
    let common_samples: HashSet<_> = matrix_samples.intersection(&recomb_samples).collect();
    let missing_from_matrix: HashSet<_> = recomb_samples.difference(&matrix_samples).collect();
    
    println!("\n=== SAMPLE OVERLAP ANALYSIS ===");
    println!("📊 Samples in matrix: {}", matrix_samples.len());
    println!("📊 Samples in recombination file: {}", recomb_samples.len());
    println!("✅ Common samples: {}", common_samples.len());
    println!("❌ Missing from matrix: {}", missing_from_matrix.len());
    
    if !missing_from_matrix.is_empty() {
        println!("\n❓ First 20 samples missing from matrix:");
        let mut missing_list: Vec<_> = missing_from_matrix.iter().collect();
        missing_list.sort();
        for (i, sample) in missing_list.iter().take(20).enumerate() {
            println!("  {}: {}", i+1, sample);
        }
        
        if missing_from_matrix.len() > 20 {
            println!("  ... and {} more", missing_from_matrix.len() - 20);
        }
    }
    
    if !common_samples.is_empty() {
        println!("\n✅ First 10 common samples:");
        let mut common_list: Vec<_> = common_samples.iter().collect();
        common_list.sort();
        for (i, sample) in common_list.iter().take(10).enumerate() {
            println!("  {}: {}", i+1, sample);
        }
    }
    
    let coverage = (common_samples.len() as f64 / recomb_samples.len() as f64) * 100.0;
    println!("\n📈 Coverage: {:.1}% of recombination samples are in the matrix", coverage);
    
    if coverage < 50.0 {
        println!("⚠️  LOW COVERAGE: The matrix might have been generated with different sample filtering");
    } else {
        println!("✅ Good coverage for verification");
    }
    
    Ok(())
}