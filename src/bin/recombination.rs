use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{Write, BufRead, BufReader};
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicUsize, Ordering};
use rayon::prelude::*;
use argh::FromArgs;
use cgdist::core::distance::ModernCache;
use cgdist::hashers::{Crc32Hasher, AlleleHasher, AlleleHash};
use lz4_flex::decompress_size_prepended;

#[derive(FromArgs)]
/// cgDist recombination analyzer - Detect genetic recombination events from enriched cache data
struct Args {
    /// path to enriched cache file (.bin extension)
    #[argh(option)]
    cache_file: String,

    /// path to allelic profile matrix (.tsv or .csv)
    #[argh(option)]
    profiles: String,

    /// path to loci filter file (one locus per line) or "NONE" to use all loci
    #[argh(option)]
    include_loci_list: Option<String>,

    /// minimum mutation density threshold for recombination detection (default: 3.0%)
    #[argh(option, default = "3.0")]
    threshold: f64,

    /// missing data character (default: -)
    #[argh(option, default = "String::from(\"-\")")]
    missing_char: String,

    /// maximum Hamming distance between samples for analysis (default: 15)
    #[argh(option, default = "15")]
    hamming_threshold: u32,

    /// minimum locus completeness threshold percentage (default: 0.0 = no filter)
    #[argh(option, default = "0.0")]
    locus_threshold: f64,

    /// minimum sample completeness threshold percentage (default: 0.0 = no filter)
    #[argh(option, default = "0.0")]
    sample_threshold: f64,

    /// output directory for analysis results (default: current directory)
    #[argh(option)]
    output: Option<String>,
}

fn load_efsa_loci(efsa_loci_path: &str) -> Result<HashSet<String>, Box<dyn std::error::Error>> {
    println!("📋 Loading EFSA loci filter from {}...", efsa_loci_path);
    
    let file = File::open(efsa_loci_path)?;
    let reader = BufReader::new(file);
    let mut loci = HashSet::new();
    
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if !line.is_empty() && !line.starts_with('#') {
            // Take the first column (locus name)
            if let Some(locus) = line.split('\t').next() {
                loci.insert(locus.to_string());
            }
        }
    }
    
    println!("✅ Loaded {} EFSA loci for filtering", loci.len());
    Ok(loci)
}

fn filter_loci_by_completeness(sample_profiles: &HashMap<String, HashMap<String, String>>, efsa_loci: &HashSet<String>, missing_char: &str, completeness_threshold: f64) -> HashSet<String> {
    if completeness_threshold <= 0.0 {
        println!("📋 No completeness filtering applied");
        return efsa_loci.clone();
    }
    
    println!("🔍 Filtering loci by completeness >= {}%...", completeness_threshold);
    let total_samples = sample_profiles.len() as f64;
    let mut filtered_loci = HashSet::new();
    
    for locus in efsa_loci {
        let present_count = sample_profiles.values()
            .filter(|profile| {
                if let Some(allele) = profile.get(locus) {
                    allele != "0" && allele != missing_char && !allele.is_empty()
                } else {
                    false
                }
            })
            .count() as f64;
            
        let completeness = (present_count / total_samples) * 100.0;
        if completeness >= completeness_threshold {
            filtered_loci.insert(locus.clone());
        }
    }
    
    println!("✅ Filtered {} loci: {} → {} ({:.1}% retained)", 
        efsa_loci.len(), efsa_loci.len(), filtered_loci.len(), 
        (filtered_loci.len() as f64 / efsa_loci.len() as f64) * 100.0);
    
    filtered_loci
}

fn create_mapping_from_profiles(profiles_path: &str, efsa_loci: &HashSet<String>, missing_char: &str) -> Result<(HashMap<(String, u32), String>, HashMap<String, HashMap<String, String>>), Box<dyn std::error::Error>> {
    println!("📋 Creating CRC → Sample mapping from profiles...");
    
    let file = File::open(profiles_path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    
    // Read header to get loci names
    let header_line = lines.next().ok_or("Empty profiles file")??;
    let loci_names: Vec<&str> = header_line.split('\t').skip(1).collect(); // Skip sample column
    
    let mut profile_mapping: HashMap<(String, u32), String> = HashMap::new();
    let mut sample_profiles: HashMap<String, HashMap<String, String>> = HashMap::new();
    let hasher = Crc32Hasher;
    
    let mut processed_count = 0;
    
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        
        if parts.is_empty() {
            continue;
        }
        
        let sample_name = parts[0].to_string();
        let mut sample_profile = HashMap::new();
        
        // Process each locus for this sample
        for (locus_idx, &allele_str) in parts.iter().skip(1).enumerate() {
            if locus_idx < loci_names.len() {
                let locus_name = loci_names[locus_idx];
                
                // Store allele for Hamming distance calculation (all loci, not just EFSA-filtered)
                sample_profile.insert(locus_name.to_string(), allele_str.to_string());
                
                // Apply EFSA loci filter for CRC mapping
                if !efsa_loci.contains(locus_name) {
                    continue;
                }
                
                // Skip missing data
                if allele_str == "0" || allele_str == missing_char || allele_str.is_empty() {
                    continue;
                }
                
                // Parse the allele to get CRC
                match hasher.parse_allele(allele_str, missing_char)? {
                    AlleleHash::Crc32(crc) => {
                        // Store mapping: (locus_name, crc) → sample_name
                        profile_mapping.insert((locus_name.to_string(), crc), sample_name.clone());
                    }
                    AlleleHash::Missing => {
                        // Skip missing alleles
                        continue;
                    }
                    _ => {
                        // For non-CRC32 hashes, we can't use them directly
                        continue;
                    }
                }
            }
        }
        
        // Store the complete sample profile
        sample_profiles.insert(sample_name, sample_profile);
        
        processed_count += 1;
        if processed_count % 100 == 0 {
            println!("  Processed {} samples...", processed_count);
        }
    }
    
    println!("✅ Created profile mapping: {} sample entries", profile_mapping.len());
    println!("✅ Stored {} complete sample profiles", sample_profiles.len());
    Ok((profile_mapping, sample_profiles))
}

fn calculate_hamming_distance_matrix_parallel(
    sample_profiles: &HashMap<String, HashMap<String, String>>, 
    efsa_loci: &HashSet<String>, 
    missing_char: &str,
    threshold: u32
) -> HashSet<(String, String)> {
    println!("🔢 Computing Hamming distance matrix with parallel processing...");
    
    let samples: Vec<String> = sample_profiles.keys().cloned().collect();
    let total_pairs = samples.len() * (samples.len() - 1) / 2;
    
    println!("   {} samples, {} total pairs to compute", samples.len(), total_pairs);
    println!("   Using {} threads with threshold <= {}", rayon::current_num_threads(), threshold);
    
    let valid_pairs = Arc::new(Mutex::new(HashSet::new()));
    let processed_count = Arc::new(AtomicUsize::new(0));
    
    // Process pairs in parallel using rayon
    (0..samples.len()).into_par_iter().for_each(|i| {
        for j in (i + 1)..samples.len() {
            let sample1 = &samples[i];
            let sample2 = &samples[j];
            
            if let (Some(profile1), Some(profile2)) = 
                (sample_profiles.get(sample1), sample_profiles.get(sample2)) {
                
                let hamming_dist: u32 = efsa_loci.iter()
                    .map(|locus| {
                        let allele1 = profile1.get(locus).map(|s| s.as_str()).unwrap_or(missing_char);
                        let allele2 = profile2.get(locus).map(|s| s.as_str()).unwrap_or(missing_char);
                        
                        let allele1_present = allele1 != "0" && allele1 != missing_char && !allele1.is_empty();
                        let allele2_present = allele2 != "0" && allele2 != missing_char && !allele2.is_empty();
                        
                        if allele1_present && allele2_present {
                            if allele1 != allele2 { 1 } else { 0 }
                        } else {
                            0
                        }
                    })
                    .sum();
                
                if hamming_dist <= threshold {
                    valid_pairs.lock().unwrap().insert((sample1.clone(), sample2.clone()));
                }
            }
            
            let count = processed_count.fetch_add(1, Ordering::Relaxed) + 1;
            if count % 50000 == 0 {
                println!("  Processed {} / {} pairs ({:.1}%)...", count, total_pairs, (count as f64 / total_pairs as f64) * 100.0);
            }
        }
    });
    
    let result = valid_pairs.lock().unwrap().clone();
    println!("✅ Matrix computation completed: {} pairs within threshold", result.len());
    
    // Save the Hamming distance matrix for verification
    println!("💾 Saving Hamming distance matrix for verification...");
    if let Err(e) = save_hamming_matrix(&samples, sample_profiles, efsa_loci, missing_char, &result, ".") {
        eprintln!("Warning: Failed to save Hamming matrix: {}", e);
    }
    
    result
}

fn save_hamming_matrix(
    samples: &[String],
    sample_profiles: &HashMap<String, HashMap<String, String>>,
    efsa_loci: &HashSet<String>,
    missing_char: &str,
    valid_pairs: &HashSet<(String, String)>,
    output_dir: &str
) -> Result<(), Box<dyn std::error::Error>> {
    let matrix_file = format!("{}/hamming_distance_matrix.tsv", output_dir);
    let mut file = File::create(&matrix_file)?;
    
    // Write header
    write!(file, "Sample1\tSample2\tHammingDistance")?;
    for locus in efsa_loci.iter().take(10) { // Show first 10 loci for verification
        write!(file, "\t{}", locus)?;
    }
    writeln!(file)?;
    
    // Write matrix data for valid pairs only
    for (sample1, sample2) in valid_pairs {
        if let (Some(profile1), Some(profile2)) = 
            (sample_profiles.get(sample1), sample_profiles.get(sample2)) {
            
            let hamming_dist: u32 = efsa_loci.iter()
                .map(|locus| {
                    let allele1 = profile1.get(locus).map(|s| s.as_str()).unwrap_or(missing_char);
                    let allele2 = profile2.get(locus).map(|s| s.as_str()).unwrap_or(missing_char);
                    
                    let allele1_present = allele1 != "0" && allele1 != missing_char && !allele1.is_empty();
                    let allele2_present = allele2 != "0" && allele2 != missing_char && !allele2.is_empty();
                    
                    if allele1_present && allele2_present {
                        if allele1 != allele2 { 1 } else { 0 }
                    } else {
                        0
                    }
                })
                .sum();
            
            write!(file, "{}\t{}\t{}", sample1, sample2, hamming_dist)?;
            
            // Add first 10 loci values for verification
            for locus in efsa_loci.iter().take(10) {
                let allele1 = profile1.get(locus).map(|s| s.as_str()).unwrap_or(missing_char);
                let allele2 = profile2.get(locus).map(|s| s.as_str()).unwrap_or(missing_char);
                write!(file, "\t{}|{}", allele1, allele2)?;
            }
            writeln!(file)?;
        }
    }
    
    println!("✅ Hamming distance matrix saved to {}", matrix_file);
    Ok(())
}

fn load_cache(cache_path: &str) -> Result<ModernCache, Box<dyn std::error::Error>> {
    println!("📂 Loading cache...");
    let compressed = std::fs::read(cache_path)
        .map_err(|e| format!("Failed to read cache file: {}", e))?;
    
    println!("📦 Compressed size: {} MB", compressed.len() / 1_000_000);
    
    println!("🔓 Decompressing cache...");
    let decompressed = decompress_size_prepended(&compressed)
        .map_err(|e| format!("Failed to decompress cache: {}", e))?;
    
    println!("🔄 Parsing cache data efficiently...");
    // Use serde_json instead of bincode as in recombination_analyzer.rs
    let cache: ModernCache = serde_json::from_slice(&decompressed)
        .map_err(|e| format!("Failed to deserialize cache: {}", e))?;
    
    Ok(cache)
}

fn filter_samples_by_completeness(
    sample_profiles: HashMap<String, HashMap<String, String>>, 
    efsa_loci: &HashSet<String>, 
    missing_char: &str, 
    completeness_threshold: f64
) -> HashMap<String, HashMap<String, String>> {
    if completeness_threshold <= 0.0 {
        println!("📋 No sample completeness filtering applied");
        return sample_profiles;
    }
    
    println!("🔍 Filtering samples by completeness >= {}%...", completeness_threshold);
    let mut filtered_samples = HashMap::new();
    let initial_count = sample_profiles.len();
    
    for (sample_name, profile) in sample_profiles {
        let present_count = efsa_loci.iter()
            .filter(|locus| {
                if let Some(allele) = profile.get(*locus) {
                    allele != "0" && allele != missing_char && !allele.is_empty()
                } else {
                    false
                }
            })
            .count() as f64;
            
        let completeness = (present_count / efsa_loci.len() as f64) * 100.0;
        if completeness >= completeness_threshold {
            filtered_samples.insert(sample_name, profile);
        }
    }
    
    println!("✅ Filtered {} samples: {} → {} ({:.1}% retained)", 
        initial_count, initial_count, filtered_samples.len(), 
        (filtered_samples.len() as f64 / initial_count as f64) * 100.0);
    
    filtered_samples
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Args = argh::from_env();

    let cache_path = &args.cache_file;
    let profiles_path = &args.profiles;
    let include_loci_list = args.include_loci_list.as_deref().unwrap_or("NONE");
    let threshold_percent = args.threshold;
    let missing_char = &args.missing_char;
    let hamming_threshold = args.hamming_threshold;
    let locus_threshold = args.locus_threshold;
    let sample_threshold = args.sample_threshold;
    
    println!("🔍 Analyzing cache for recombination events");
    println!("📂 Cache file: {}", cache_path);
    println!("📋 Profiles file: {}", profiles_path);
    println!("🎯 Include loci list: {}", include_loci_list);
    println!("🎯 Threshold: {}%", threshold_percent);
    println!("🚫 Missing char: '{}'", missing_char);
    println!("🎯 Hamming threshold: <= {}", hamming_threshold);
    if locus_threshold > 0.0 {
        println!("🔍 Locus threshold: >= {}%", locus_threshold);
    }
    if sample_threshold > 0.0 {
        println!("🔍 Sample threshold: >= {}%", sample_threshold);
    }
    
    // Load loci filter or use ALL loci
    let efsa_loci = if include_loci_list == "NONE" {
        println!("📋 Using ALL loci (no filtering)...");
        HashSet::new() // Empty set - will be populated from profile headers
    } else {
        load_efsa_loci(include_loci_list)?
    };
    
    // Load cache
    let cache = load_cache(cache_path)?;
    
    println!("✅ Cache loaded successfully");
    println!("   Total entries: {}", cache.metadata.total_entries);
    println!("   Unique loci: {}", cache.metadata.unique_loci);
    
    // Create mapping from profiles: (locus, crc) → sample and get complete profiles
    let (profile_mapping, sample_profiles) = create_mapping_from_profiles(profiles_path, &efsa_loci, missing_char)?;
    
    // If no loci filter specified (NONE), use all loci from the profiles
    let final_efsa_loci = if efsa_loci.is_empty() {
        // Extract all loci from the first sample profile
        if let Some((_, first_profile)) = sample_profiles.iter().next() {
            let all_loci: HashSet<String> = first_profile.keys().cloned().collect();
            println!("✅ Using ALL {} loci from profiles", all_loci.len());
            all_loci
        } else {
            println!("⚠️ No profiles found, using empty loci set");
            HashSet::new()
        }
    } else {
        efsa_loci
    };
    
    // Apply locus completeness filtering in memory
    let filtered_efsa_loci = filter_loci_by_completeness(&sample_profiles, &final_efsa_loci, missing_char, locus_threshold);
    
    // Apply sample completeness filtering
    let filtered_sample_profiles = filter_samples_by_completeness(sample_profiles, &filtered_efsa_loci, missing_char, sample_threshold);
    
    // Calculate Hamming distance matrix and filter by threshold
    let valid_pairs = calculate_hamming_distance_matrix_parallel(&filtered_sample_profiles, &filtered_efsa_loci, missing_char, hamming_threshold);
    
    // Save the Hamming distance matrix for verification
    println!("💾 Saving Hamming distance matrix for verification...");
    let filtered_samples: Vec<String> = filtered_sample_profiles.keys().cloned().collect();
    let output_dir = args.output.as_deref().unwrap_or(".");
    if let Err(e) = save_hamming_matrix(&filtered_samples, &filtered_sample_profiles, &filtered_efsa_loci, missing_char, &valid_pairs, output_dir) {
        eprintln!("Warning: Failed to save Hamming matrix: {}", e);
    }
    
    println!("🎯 Proceeding with {} sample pairs within Hamming threshold", valid_pairs.len());
    
    // Process with efficient memory usage - using pre-filtered pairs
    let mut recombination_events = Vec::new();
    let mut total_pairs_with_lengths = 0;
    let mut locus_counts: HashMap<String, u32> = HashMap::new();
    let mut processed = 0;
    let mut processed_pairs: HashSet<(String, u32, u32)> = HashSet::new(); // Track (locus, min_crc, max_crc)
    let mut pairwise_recombination: HashMap<(String, String), HashSet<String>> = HashMap::new(); // Track (sample1, sample2) → set of recombining loci
    
    println!("🚀 Processing all {} entries with optimized algorithm...", cache.metadata.total_entries);
    
    for (key, entry) in cache.data.iter() {
        processed += 1;
        if processed % 500000 == 0 {
            println!("   Processed: {}/{} entries ({:.1}%) - {} recombination events found", 
                processed, 
                cache.metadata.total_entries, 
                (processed as f64 / cache.metadata.total_entries as f64) * 100.0,
                recombination_events.len());
        }
        
        // Parse key: locus:crc1:crc2
        let parts: Vec<&str> = key.split(':').collect();
        if parts.len() != 3 {
            continue;
        }
        
        let locus = parts[0];
        let crc1_str = parts[1];
        let crc2_str = parts[2];
        
        // Apply EFSA loci filter
        if !filtered_efsa_loci.contains(locus) {
            continue;
        }
        
        // Skip same allele pairs
        if crc1_str == crc2_str {
            continue;
        }
        
        // Parse CRC values
        let crc1: u32 = crc1_str.parse().unwrap_or(0);
        let crc2: u32 = crc2_str.parse().unwrap_or(0);
        
        // Normalize CRC pair: always use (min_crc, max_crc) to avoid duplicates
        let (min_crc, max_crc) = if crc1 <= crc2 { (crc1, crc2) } else { (crc2, crc1) };
        
        // Check if we already processed this pair
        let pair_key = (locus.to_string(), min_crc, max_crc);
        if processed_pairs.contains(&pair_key) {
            continue; // Skip duplicate
        }
        processed_pairs.insert(pair_key);
        
        // Find sample names using the profile mapping: (locus, crc) → sample
        let sample1_opt = profile_mapping.get(&(locus.to_string(), crc1));
        let sample2_opt = profile_mapping.get(&(locus.to_string(), crc2));
        
        // Skip if one or both CRCs correspond to missing data (not found in profiles)
        if sample1_opt.is_none() || sample2_opt.is_none() {
            continue; // Skip pairs involving missing alleles
        }
        
        let sample1 = sample1_opt.unwrap();
        let sample2 = sample2_opt.unwrap();
        
        // Check if this sample pair is within Hamming threshold (pre-computed)
        let normalized_pair = if sample1 <= sample2 {
            (sample1.clone(), sample2.clone())
        } else {
            (sample2.clone(), sample1.clone())
        };
        
        if !valid_pairs.contains(&normalized_pair) {
            continue; // Skip pairs not within Hamming threshold
        }
        
        // Check if we have sequence lengths for enriched cache
        if let (Some(len1), Some(len2)) = (entry.seq1_length, entry.seq2_length) {
            if len1 == 0 || len2 == 0 {
                continue;
            }
            
            total_pairs_with_lengths += 1;
            
            // Calculate average length
            let avg_length = (len1 + len2) as f64 / 2.0;
            
            // Calculate separate densities
            let total_mutations = entry.snps + entry.indel_events;
            let total_density = (total_mutations as f64 / avg_length) * 100.0;
            let snp_density = (entry.snps as f64 / avg_length) * 100.0;
            let indel_density = (entry.indel_events as f64 / avg_length) * 100.0;
            
            // Check if exceeds threshold (using total density)
            if total_density > threshold_percent {
                recombination_events.push((
                    sample1.to_string(),
                    sample2.to_string(),
                    locus.to_string(),
                    crc1_str.to_string(),
                    crc2_str.to_string(),
                    entry.snps,
                    entry.indel_events,
                    total_mutations,
                    avg_length,
                    snp_density,
                    indel_density,
                    total_density
                ));
                
                *locus_counts.entry(locus.to_string()).or_insert(0) += 1;
                
                // Track pairwise recombination - normalize sample pair order
                let sample_pair = if sample1 <= sample2 { 
                    (sample1.to_string(), sample2.to_string()) 
                } else { 
                    (sample2.to_string(), sample1.to_string()) 
                };
                pairwise_recombination.entry(sample_pair).or_insert_with(HashSet::new).insert(locus.to_string());
            }
        }
    }
    
    println!("\n=== RECOMBINATION ANALYSIS RESULTS ===");
    println!("Total pairs with length data: {}", total_pairs_with_lengths);
    println!("Recombination events detected: {}", recombination_events.len());
    println!("Percentage with recombination: {:.2}%", 
             (recombination_events.len() as f64 / total_pairs_with_lengths as f64) * 100.0);
    
    // Sort by mutation density (index 9 now)
    recombination_events.sort_by(|a, b| b.9.partial_cmp(&a.9).unwrap());
    
    // Show top loci
    let mut locus_vec: Vec<_> = locus_counts.iter().collect();
    locus_vec.sort_by(|a, b| b.1.cmp(a.1));
    
    println!("\n=== TOP 20 LOCI WITH RECOMBINATION ===");
    for (locus, count) in locus_vec.iter().take(20) {
        println!("  {}: {} allele pairs", locus, count);
    }
    
    // Show top events
    println!("\n=== TOP 30 RECOMBINATION EVENTS ===");
    println!("Sample1\t\tSample2\t\tLocus\t\t\tAllele1\t\tAllele2\t\tSNPs\tIndelEvents\tTotalMutations\tTotalDensity%");
    println!("--------------------------------------------------------------------------------------------------------------");
    
    for (sample1, sample2, locus, allele1, allele2, snps, indel_events, total_mutations, _avg_length, _snp_density, _indel_density, total_density) in recombination_events.iter().take(30) {
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{:.2}%", 
                 sample1, sample2, locus, allele1, allele2, snps, indel_events, total_mutations, total_density);
    }
    
    // Determine output directory
    let output_dir = args.output.as_deref().unwrap_or(".");

    // Write TSV output
    let output_file = format!("{}/recombination_analysis_with_samples.tsv", output_dir);
    let mut file = File::create(&output_file)?;
    
    writeln!(file, "Sample1\tSample2\tLocus\tAllele1\tAllele2\tSNPs\tIndelEvents\tTotalMutations\tAvgLength\tSNPDensity%\tIndelDensity%\tTotalDensity%")?;
    
    for (sample1, sample2, locus, allele1, allele2, snps, indel_events, total_mutations, avg_length, snp_density, indel_density, total_density) in &recombination_events {
        writeln!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}", 
                 sample1, sample2, locus, allele1, allele2, snps, indel_events, total_mutations, avg_length, snp_density, indel_density, total_density)?;
    }
    
    println!("\n✅ Results saved to {}", output_file);
    println!("   Total events: {}", recombination_events.len());
    
    // Write pairwise recombination summary
    let pairwise_file = format!("{}/pairwise_recombination_summary.tsv", output_dir);
    let mut pairwise_output = File::create(&pairwise_file)?;
    
    writeln!(pairwise_output, "Sample1\tSample2\tRecombiningLoci\tTotalEFSALoci\tRecombinationPercentage%")?;
    
    let total_efsa_loci = filtered_efsa_loci.len();
    let mut pairwise_list: Vec<_> = pairwise_recombination.iter().collect();
    pairwise_list.sort_by(|a, b| b.1.len().cmp(&a.1.len())); // Sort by number of recombining loci (descending)
    
    for ((sample1, sample2), recombining_loci) in &pairwise_list {
        let recombining_count = recombining_loci.len();
        let recombination_percentage = (recombining_count as f64 / total_efsa_loci as f64) * 100.0;
        
        writeln!(pairwise_output, "{}\t{}\t{}\t{}\t{:.2}", 
                 sample1, sample2, recombining_count, total_efsa_loci, recombination_percentage)?;
    }
    
    println!("✅ Pairwise recombination summary saved to {}", pairwise_file);
    println!("   Total sample pairs with recombination: {}", pairwise_list.len());
    
    Ok(())
}