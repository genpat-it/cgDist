// recombination_analyzer.rs - Standalone utility for recombination detection using enriched cache

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use clap::{Arg, Command};
use serde::{Serialize, Deserialize};

use cgdist::core::distance::ModernCache;

#[derive(Debug, Clone, Serialize, Deserialize)]
struct RecombinationEvent {
    sample_i: String,
    sample_j: String,
    locus: String,
    cgdist_snps: u32,
    cgdist_indel_events: u32,
    cgdist_indel_bases: u32,
    cgdist_total: f64,
    hamming_distance: u32,
    avg_sequence_length: f64,
    mutation_density_percent: f64,
    recombination_excess_percent: f64,
}

#[derive(Debug)]
struct DistanceMatrix {
    samples: Vec<String>,
    matrix: Vec<Vec<f64>>,
}

impl DistanceMatrix {
    fn load_from_tsv(path: &str) -> Result<Self, String> {
        let file = File::open(path)
            .map_err(|e| format!("Failed to open matrix file '{}': {}", path, e))?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        // Read header
        let header_line = lines.next()
            .ok_or("Empty matrix file")?
            .map_err(|e| format!("Failed to read header: {}", e))?;
        
        let samples: Vec<String> = header_line
            .split('\t')
            .skip(1) // Skip first empty cell
            .map(|s| s.trim().to_string())
            .collect();

        let n = samples.len();
        let mut matrix = vec![vec![0.0; n]; n];

        // Read matrix rows
        for (i, line) in lines.enumerate() {
            let line = line.map_err(|e| format!("Failed to read line {}: {}", i + 2, e))?;
            let values: Vec<&str> = line.split('\t').collect();
            
            if values.len() != n + 1 {
                return Err(format!("Row {} has {} columns, expected {}", i + 2, values.len(), n + 1));
            }

            for (j, value_str) in values.iter().skip(1).enumerate() {
                matrix[i][j] = value_str.trim().parse::<f64>()
                    .map_err(|e| format!("Failed to parse value at row {}, col {}: {}", i + 2, j + 1, e))?;
            }
        }

        Ok(DistanceMatrix { samples, matrix })
    }

    fn save_to_tsv(&self, path: &str) -> Result<(), String> {
        let mut file = File::create(path)
            .map_err(|e| format!("Failed to create output file '{}': {}", path, e))?;

        // Write header
        write!(file, "").map_err(|e| format!("Write error: {}", e))?;
        for sample in &self.samples {
            write!(file, "\t{}", sample).map_err(|e| format!("Write error: {}", e))?;
        }
        writeln!(file).map_err(|e| format!("Write error: {}", e))?;

        // Write matrix
        for (i, row) in self.matrix.iter().enumerate() {
            write!(file, "{}", self.samples[i]).map_err(|e| format!("Write error: {}", e))?;
            for &value in row {
                write!(file, "\t{:.6}", value).map_err(|e| format!("Write error: {}", e))?;
            }
            writeln!(file).map_err(|e| format!("Write error: {}", e))?;
        }

        println!("✅ Saved corrected distance matrix to: {}", path);
        Ok(())
    }

    fn get_sample_index(&self, sample: &str) -> Option<usize> {
        self.samples.iter().position(|s| s == sample)
    }

    fn set_distance(&mut self, sample_i: &str, sample_j: &str, distance: f64) -> bool {
        if let (Some(i), Some(j)) = (self.get_sample_index(sample_i), self.get_sample_index(sample_j)) {
            self.matrix[i][j] = distance;
            self.matrix[j][i] = distance; // Symmetric
            true
        } else {
            false
        }
    }
}

/// Load modern cache from file
fn load_modern_cache(cache_path: &str) -> Result<ModernCache, String> {
    let compressed = std::fs::read(cache_path)
        .map_err(|e| format!("Failed to read cache file: {}", e))?;
    
    // Decompress with LZ4
    let decompressed = lz4_flex::decompress_size_prepended(&compressed)
        .map_err(|e| format!("Failed to decompress cache: {}", e))?;
    
    // Load as ModernCache
    let modern_cache: ModernCache = serde_json::from_slice(&decompressed)
        .map_err(|e| format!("Failed to deserialize cache: {}", e))?;
    
    Ok(modern_cache)
}

struct RecombinationAnalyzer {
    threshold_percent: f64,
    enriched_cache: ModernCache,
    sample_to_crcs: HashMap<String, Vec<u32>>, // Sample -> list of CRCs for that sample
    crc_to_sample: HashMap<u32, String>, // CRC -> Sample name
}

impl RecombinationAnalyzer {
    fn new(threshold_percent: f64, enriched_cache: ModernCache) -> Self {
        Self {
            threshold_percent,
            enriched_cache,
            sample_to_crcs: HashMap::new(),
            crc_to_sample: HashMap::new(),
        }
    }

    fn build_crc_mappings(&mut self, profiles_path: &str) -> Result<(), String> {
        println!("📊 Building CRC to sample mappings from profiles...");
        
        let file = File::open(profiles_path)
            .map_err(|e| format!("Failed to open profiles file '{}': {}", profiles_path, e))?;
        
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        // Read header
        let header_line = lines.next()
            .ok_or("Empty profiles file")?
            .map_err(|e| format!("Failed to read header: {}", e))?;
        
        let headers: Vec<&str> = header_line.split('\t').collect();
        let loci_headers = &headers[1..]; // Skip sample name column

        // Read each sample
        for line in lines {
            let line = line.map_err(|e| format!("Failed to read line: {}", e))?;
            if line.trim().is_empty() { continue; }

            let values: Vec<&str> = line.split('\t').collect();
            if values.is_empty() { continue; }

            let sample_name = values[0].trim().to_string();
            let mut sample_crcs = Vec::new();

            // Process each allele value 
            for (locus_idx, &allele_value) in values.iter().skip(1).enumerate() {
                if locus_idx >= loci_headers.len() { break; }
                
                let allele_value = allele_value.trim();
                if allele_value.is_empty() || allele_value == "0" || allele_value == "INF-" {
                    continue; // Skip missing alleles
                }

                // Parse as CRC32 value
                if let Ok(crc) = allele_value.parse::<u32>() {
                    sample_crcs.push(crc);
                    self.crc_to_sample.insert(crc, sample_name.clone());
                }
            }

            if !sample_crcs.is_empty() {
                self.sample_to_crcs.insert(sample_name, sample_crcs);
            }
        }

        println!("✅ Built mappings for {} samples with {} total CRC mappings", 
                self.sample_to_crcs.len(), self.crc_to_sample.len());
        
        Ok(())
    }


    fn analyze_recombination(&self) -> Result<Vec<RecombinationEvent>, String> {
        println!("🧬 Analyzing recombination events with {}% threshold", self.threshold_percent);
        
        let mut events = Vec::new();
        let mut total_comparisons = 0;
        let threshold_density = self.threshold_percent / 100.0;

        // Parse cache entries (keys are strings like "locus:hash1:hash2")
        let mut locus_entries: HashMap<String, Vec<_>> = HashMap::new();
        for (key_str, entry) in &self.enriched_cache.data {
            let parts: Vec<&str> = key_str.split(':').collect();
            if parts.len() >= 3 {
                let locus = parts[0].to_string();
                locus_entries.entry(locus).or_default().push((key_str, entry));
            }
        }

        println!("📋 Processing {} loci", locus_entries.len());

        for (locus, entries) in locus_entries {
            for (key_str, cache_entry) in entries {
                // Parse key string to extract hash values
                let parts: Vec<&str> = key_str.split(':').collect();
                if parts.len() < 3 { continue; }
                
                // Parts: [locus, hash1, hash2]
                let hash1_str = parts[1];
                let hash2_str = parts[2];
                
                // Convert to u32 CRCs for sample lookup
                let hash1_crc = hash1_str.parse::<u32>().ok();
                let hash2_crc = hash2_str.parse::<u32>().ok();
                
                if let (Some(crc1), Some(crc2)) = (hash1_crc, hash2_crc) {
                    let sample_i = self.crc_to_sample.get(&crc1).cloned();
                    let sample_j = self.crc_to_sample.get(&crc2).cloned();

                    if let (Some(sample_i), Some(sample_j)) = (sample_i, sample_j) {
                    total_comparisons += 1;

                    // Get sequence lengths
                    let len1 = cache_entry.seq1_length.unwrap_or(0) as f64;
                    let len2 = cache_entry.seq2_length.unwrap_or(0) as f64;
                    
                    if len1 == 0.0 || len2 == 0.0 {
                        continue; // Skip if no length information
                    }

                    let avg_length = (len1 + len2) / 2.0;

                    // Calculate cgDist total distance
                    let cgdist_total = cache_entry.snps as f64 + 
                                     cache_entry.indel_events as f64 + 
                                     cache_entry.indel_bases as f64;

                        // Hamming distance (1 if different alleles, 0 if same)
                        let hamming_distance = if crc1 != crc2 { 1 } else { 0 };

                    // Calculate mutation density
                    let mutation_density = cgdist_total / avg_length;

                    // Check for recombination
                    if mutation_density > threshold_density && hamming_distance > 0 {
                        let recombination_excess = (mutation_density - threshold_density) / threshold_density * 100.0;

                        let event = RecombinationEvent {
                            sample_i,
                            sample_j,
                            locus: locus.clone(),
                            cgdist_snps: cache_entry.snps as u32,
                            cgdist_indel_events: cache_entry.indel_events as u32,
                            cgdist_indel_bases: cache_entry.indel_bases as u32,
                            cgdist_total,
                            hamming_distance,
                            avg_sequence_length: avg_length,
                            mutation_density_percent: mutation_density * 100.0,
                            recombination_excess_percent: recombination_excess,
                        };

                        events.push(event);
                        }
                    }
                }
            }
        }

        println!("✅ Analysis complete:");
        println!("   Total comparisons: {}", total_comparisons);
        println!("   Recombination events: {}", events.len());
        if total_comparisons > 0 {
            println!("   Recombination rate: {:.2}%", events.len() as f64 / total_comparisons as f64 * 100.0);
        }

        Ok(events)
    }

    fn generate_corrected_matrix(&self, events: &[RecombinationEvent], original_matrix_path: &str, output_matrix_path: &str) -> Result<(), String> {
        println!("📝 Generating corrected distance matrix with hamming fallback");
        
        let mut matrix = DistanceMatrix::load_from_tsv(original_matrix_path)?;
        let mut corrections = 0;

        for event in events {
            if matrix.set_distance(&event.sample_i, &event.sample_j, event.hamming_distance as f64) {
                corrections += 1;
            }
        }

        matrix.save_to_tsv(output_matrix_path)?;
        println!("✅ Applied {} corrections to distance matrix", corrections);
        
        Ok(())
    }

    fn write_recombination_log(&self, events: &[RecombinationEvent], log_path: &str) -> Result<(), String> {
        let mut file = File::create(log_path)
            .map_err(|e| format!("Failed to create log file '{}': {}", log_path, e))?;

        writeln!(file, "sample_i\tsample_j\tlocus\tcgdist_snps\tcgdist_indel_events\tcgdist_indel_bases\tcgdist_total\thamming_distance\tavg_seq_length\tmutation_density_%\trecombination_excess_%")
            .map_err(|e| format!("Failed to write header: {}", e))?;

        for event in events {
            writeln!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{:.0}\t{:.4}\t{:.2}",
                event.sample_i,
                event.sample_j,
                event.locus,
                event.cgdist_snps,
                event.cgdist_indel_events,
                event.cgdist_indel_bases,
                event.cgdist_total,
                event.hamming_distance,
                event.avg_sequence_length,
                event.mutation_density_percent,
                event.recombination_excess_percent,
            ).map_err(|e| format!("Failed to write event: {}", e))?;
        }

        println!("📋 Wrote {} recombination events to: {}", events.len(), log_path);
        Ok(())
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let matches = Command::new("Recombination Analyzer")
        .version("1.0.0")
        .about("Detects recombination events using enriched cgDist cache")
        .arg(Arg::new("enriched-cache")
            .long("enriched-cache")
            .value_name("FILE")
            .help("Path to enriched cache file (.bin)")
            .required(true))
        .arg(Arg::new("profiles")
            .long("profiles")
            .value_name("FILE")
            .help("Path to cgMLST profiles file (.tsv)")
            .required(true))
        .arg(Arg::new("distance-matrix")
            .long("distance-matrix")
            .value_name("FILE")
            .help("Path to original distance matrix (.tsv)")
            .required(true))
        .arg(Arg::new("threshold")
            .long("threshold")
            .value_name("PERCENT")
            .help("Recombination threshold as percentage (default: 3.0)")
            .default_value("3.0"))
        .arg(Arg::new("output-matrix")
            .long("output-matrix")
            .value_name("FILE")
            .help("Path to output corrected distance matrix")
            .required(true))
        .arg(Arg::new("recombination-log")
            .long("recombination-log")
            .value_name("FILE")
            .help("Path to output recombination events log")
            .required(true))
        .get_matches();

    let enriched_cache_path = matches.get_one::<String>("enriched-cache").unwrap();
    let profiles_path = matches.get_one::<String>("profiles").unwrap();
    let distance_matrix_path = matches.get_one::<String>("distance-matrix").unwrap();
    let threshold_str = matches.get_one::<String>("threshold").unwrap();
    let output_matrix_path = matches.get_one::<String>("output-matrix").unwrap();
    let recombination_log_path = matches.get_one::<String>("recombination-log").unwrap();

    let threshold = threshold_str.parse::<f64>()
        .map_err(|_| format!("Invalid threshold value: {}", threshold_str))?;

    println!("🔬 cgDist Recombination Analyzer");
    println!("================================");
    println!("📂 Enriched cache: {}", enriched_cache_path);
    println!("📊 Profiles: {}", profiles_path);
    println!("📈 Distance matrix: {}", distance_matrix_path);
    println!("🎯 Threshold: {}%", threshold);
    println!();

    // Load enriched cache
    println!("📂 Loading enriched cache...");
    let enriched_cache = load_modern_cache(enriched_cache_path)
        .map_err(|e| format!("Failed to load cache: {}", e))?;
    
    println!("✅ Loaded cache with {} entries", enriched_cache.data.len());

    // Initialize analyzer
    let mut analyzer = RecombinationAnalyzer::new(threshold, enriched_cache);

    // Build CRC mappings
    analyzer.build_crc_mappings(profiles_path)?;

    // Analyze recombination
    let events = analyzer.analyze_recombination()?;

    // Generate outputs
    analyzer.generate_corrected_matrix(&events, distance_matrix_path, output_matrix_path)?;
    analyzer.write_recombination_log(&events, recombination_log_path)?;

    println!();
    println!("🎉 Recombination analysis complete!");
    println!("📄 Corrected matrix: {}", output_matrix_path);
    println!("📋 Recombination log: {}", recombination_log_path);

    Ok(())
}