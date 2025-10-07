use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn load_distance_matrix(file_path: &str) -> Result<HashMap<(String, String), u32>, Box<dyn std::error::Error>> {
    println!("📋 Loading distance matrix from {}...", file_path);
    
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    
    // Skip comments
    while let Some(line) = lines.next() {
        let line = line?;
        if !line.starts_with('#') {
            break;
        }
    }
    
    // Read header with sample names
    let header_line = lines.next().ok_or("Missing header line")??;
    let samples: Vec<&str> = header_line.split('\t').collect();
    
    let mut matrix = HashMap::new();
    
    // Read matrix data
    for (row_idx, line) in lines.enumerate() {
        let line = line?;
        let distances: Vec<&str> = line.split('\t').collect();
        
        if distances.len() < 2 {
            continue; // Skip empty lines
        }
        
        let row_sample = distances[0];
        
        // Parse distances to other samples
        for (col_idx, &distance_str) in distances.iter().skip(1).enumerate() {
            if col_idx < samples.len() - 1 { // Skip the first column which is "Sample"
                let col_sample = samples[col_idx + 1]; // +1 because samples[0] is "Sample"
                
                if let Ok(distance) = distance_str.parse::<u32>() {
                    // Store both orientations
                    matrix.insert((row_sample.to_string(), col_sample.to_string()), distance);
                    matrix.insert((col_sample.to_string(), row_sample.to_string()), distance);
                }
            }
        }
    }
    
    println!("✅ Loaded {} distance pairs", matrix.len());
    Ok(matrix)
}

fn load_recombination_tsv(file_path: &str) -> Result<Vec<(String, String, String, u32, u32)>, Box<dyn std::error::Error>> {
    println!("📋 Loading recombination events from {}...", file_path);
    
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    
    // Skip header
    lines.next();
    
    let mut events = Vec::new();
    
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        
        if parts.len() >= 7 { // We need at least Sample1, Sample2, Locus, Allele1, Allele2, SNPs, IndelEvents
            let sample1 = parts[0].to_string();
            let sample2 = parts[1].to_string();
            let locus = parts[2].to_string();
            
            if let (Ok(snps), Ok(indel_events)) = (parts[5].parse::<u32>(), parts[6].parse::<u32>()) {
                events.push((sample1, sample2, locus, snps, indel_events));
            }
        }
    }
    
    println!("✅ Loaded {} recombination events", events.len());
    Ok(events)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    
    if args.len() != 4 {
        eprintln!("Usage: {} <recombination.tsv> <snps_matrix.tsv> <indel_matrix.tsv>", args[0]);
        std::process::exit(1);
    }
    
    let recombination_file = &args[1];
    let snps_matrix_file = &args[2];
    let indel_matrix_file = &args[3];
    
    println!("🔍 Verifying recombination distances against reference matrices");
    println!("📁 Recombination file: {}", recombination_file);
    println!("📁 SNPs matrix: {}", snps_matrix_file);
    println!("📁 Indel matrix: {}", indel_matrix_file);
    
    // Load reference matrices
    let snps_matrix = load_distance_matrix(snps_matrix_file)?;
    let indel_matrix = load_distance_matrix(indel_matrix_file)?;
    
    // Load recombination events
    let all_events = load_recombination_tsv(recombination_file)?;
    
    // Filter events to only include those where both samples exist in matrices
    let events: Vec<_> = all_events.iter()
        .filter(|(sample1, sample2, _, _, _)| {
            snps_matrix.contains_key(&(sample1.clone(), sample2.clone())) &&
            indel_matrix.contains_key(&(sample1.clone(), sample2.clone()))
        })
        .collect();
    
    println!("\n📊 Events analysis:");
    println!("   Total events: {}", all_events.len());
    println!("   Verifiable events (both samples in matrix): {}", events.len());
    println!("   Coverage: {:.1}%", (events.len() as f64 / all_events.len() as f64) * 100.0);
    
    if events.is_empty() {
        println!("❌ No events can be verified - no matching sample pairs in matrices!");
        return Ok(());
    }
    
    println!("\n🔍 Sampling and verifying distances...");
    
    let sample_size = std::cmp::min(50, events.len()); // Sample up to 50 events
    let mut verified_count = 0;
    let mut snp_mismatches = 0;
    let mut indel_mismatches = 0;
    let mut missing_in_matrix = 0;
    
    // Sample events evenly across the list
    let step = if events.len() > sample_size { events.len() / sample_size } else { 1 };
    
    println!("Sample1\t\tSample2\t\tLocus\t\t\tRec_SNPs\tMatrix_SNPs\tRec_Indels\tMatrix_Indels\tStatus");
    println!("--------------------------------------------------------------------------------------------------------");
    
    for i in 0..sample_size {
        let idx = i * step;
        if idx >= events.len() {
            break;
        }
        
        let (sample1, sample2, locus, rec_snps, rec_indels) = events[idx];
        
        // Look up distances in reference matrices (we know they exist due to filtering)
        let matrix_snps = *snps_matrix.get(&(sample1.clone(), sample2.clone())).unwrap();
        let matrix_indels = *indel_matrix.get(&(sample1.clone(), sample2.clone())).unwrap();
        
        let snp_match = *rec_snps == matrix_snps;
        let indel_match = *rec_indels == matrix_indels;
        
        let status = if snp_match && indel_match {
            "✅ MATCH"
        } else {
            if !snp_match { snp_mismatches += 1; }
            if !indel_match { indel_mismatches += 1; }
            "❌ MISMATCH"
        };
        
        if snp_match && indel_match {
            verified_count += 1;
        }
        
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            sample1, sample2, locus, rec_snps, matrix_snps, rec_indels, matrix_indels, status);
    }
    
    println!("\n=== VERIFICATION SUMMARY ===");
    println!("Total events sampled: {}", sample_size);
    println!("✅ Verified matches: {}", verified_count);
    println!("❌ SNP mismatches: {}", snp_mismatches);
    println!("❌ Indel mismatches: {}", indel_mismatches);
    println!("❓ Missing from matrices: {}", missing_in_matrix);
    println!("📊 Success rate: {:.1}%", (verified_count as f64 / sample_size as f64) * 100.0);
    
    if verified_count == sample_size - missing_in_matrix {
        println!("\n🎉 ALL VERIFIABLE DISTANCES MATCH! ✅");
    } else {
        println!("\n⚠️  Some distances don't match - investigation needed");
    }
    
    Ok(())
}