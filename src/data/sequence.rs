// sequence.rs - Sequence database for FASTA schema handling

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use bio::io::fasta;
use crc32fast::Hasher;

/// Information about a single sequence
#[derive(Debug, Clone)]
pub struct SequenceInfo {
    pub sequence: Vec<u8>,
    pub id: String,
}

/// Database of sequences organized by locus
#[derive(Debug)]
#[derive(Clone)]
pub struct SequenceDatabase {
    pub loci: HashMap<String, HashMap<u32, SequenceInfo>>,
    pub total_sequences: usize,
}

impl SequenceDatabase {
    pub fn new() -> Self {
        Self {
            loci: HashMap::new(),
            total_sequences: 0,
        }
    }
    
    pub fn empty() -> Self {
        Self::new()
    }
    
    /// Add a sequence to the database
    pub fn add_sequence(&mut self, locus: String, crc: u32, sequence_info: SequenceInfo) {
        self.loci.entry(locus)
            .or_default()
            .insert(crc, sequence_info);
        self.total_sequences += 1;
    }
    
    /// Get a sequence by locus and CRC
    pub fn get_sequence(&self, locus: &str, crc: u32) -> Option<&SequenceInfo> {
        self.loci.get(locus)?.get(&crc)
    }
    
    /// Get all sequences for a locus
    pub fn get_locus_sequences(&self, locus: &str) -> Option<&HashMap<u32, SequenceInfo>> {
        self.loci.get(locus)
    }
    
    /// Check if a sequence exists
    pub fn has_sequence(&self, locus: &str, crc: u32) -> bool {
        self.loci.get(locus)
            .map(|sequences| sequences.contains_key(&crc))
            .unwrap_or(false)
    }
    
    /// Get statistics about the database
    pub fn get_stats(&self) -> (usize, usize) {
        (self.loci.len(), self.total_sequences)
    }
}

impl Default for SequenceDatabase {
    fn default() -> Self {
        Self::new()
    }
}

/// Compute CRC32 hash of a sequence
fn compute_crc32(sequence: &[u8]) -> u32 {
    let mut hasher = Hasher::new();
    hasher.update(sequence);
    hasher.finalize()
}

impl SequenceDatabase {
    /// Load database from schema file or directory
    pub fn from_schema_file(schema_path: &Path) -> Result<Self, String> {
        let mut db = Self::new();
        
        if schema_path.is_dir() {
            println!("🧬 Loading FASTA schema from directory: {}", schema_path.display());
            db.load_from_directory(schema_path)?;
        } else if schema_path.is_file() {
            println!("🧬 Loading FASTA schema from file: {}", schema_path.display());
            db.load_from_schema_file(schema_path)?;
        } else {
            return Err(format!("Schema path does not exist: {}", schema_path.display()));
        }
        
        println!("✅ Schema loaded: {} loci, {} sequences", db.loci.len(), db.total_sequences);
        Ok(db)
    }
    
    /// Load database selectively for specific CRC pairs only
    pub fn from_schema_selective(
        schema_path: &Path, 
        unique_pairs: &HashSet<(String, u32, u32)>,
        _loci_names: &[String]
    ) -> Result<Self, String> {
        let mut db = Self::new();
        
        // Collect all required CRCs per locus
        let mut required_crcs: HashMap<String, HashSet<u32>> = HashMap::new();
        for (locus, crc1, crc2) in unique_pairs {
            let entry = required_crcs.entry(locus.clone()).or_default();
            entry.insert(*crc1);
            entry.insert(*crc2);
        }
        
        println!("🧬 Loading FASTA schema selectively from directory: {}", schema_path.display());
        let total_crcs = required_crcs.values().map(|s| s.len()).sum::<usize>();
        println!("📊 Will load {} loci with {} total required CRCs", 
                required_crcs.len(), total_crcs);
        
        use indicatif::{ProgressBar, ProgressStyle};
        let pb = ProgressBar::new(required_crcs.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} loci loaded")
                .unwrap()
        );
        
        if schema_path.is_dir() {
            db.load_from_directory_selective(schema_path, &required_crcs, &pb)?;
        } else {
            return Err("Selective loading only supports directory schemas".to_string());
        }
        
        pb.finish_with_message("✅ Schema loaded!");
        println!("✅ Schema loaded selectively: {} loci, {} sequences", db.loci.len(), db.total_sequences);
        Ok(db)
    }
    
    /// Load from schema directory selectively
    fn load_from_directory_selective(
        &mut self, 
        schema_dir: &Path, 
        required_crcs: &HashMap<String, HashSet<u32>>,
        pb: &indicatif::ProgressBar
    ) -> Result<(), String> {
        use rayon::prelude::*;
        use std::sync::atomic::{AtomicUsize, Ordering};
        
        let processed = AtomicUsize::new(0);
        
        // Collect all locus data in parallel
        let results: Vec<_> = required_crcs
            .par_iter()
            .map(|(locus_name, crcs_needed)| {
                let fasta_path = schema_dir.join(format!("{}.fasta", locus_name));
                if !fasta_path.exists() {
                    return Ok((locus_name.clone(), HashMap::new(), 0));
                }
                
                // Load sequences for this locus
                let file = File::open(&fasta_path)
                    .map_err(|e| format!("Failed to open FASTA file {}: {}", fasta_path.display(), e))?;
                
                let reader = fasta::Reader::new(BufReader::new(file));
                let mut locus_sequences = HashMap::new();
                let mut sequences_loaded = 0;
                
                for record_result in reader.records() {
                    let record = record_result
                        .map_err(|e| format!("Invalid FASTA record in {}: {}", fasta_path.display(), e))?;
                    
                    let sequence = record.seq().to_vec();
                    let crc32 = compute_crc32(&sequence);
                    
                    // Only load if we need this CRC
                    if crcs_needed.contains(&crc32) {
                        let id = record.id().to_string();
                        let seq_info = SequenceInfo { sequence, id };
                        locus_sequences.insert(crc32, seq_info);
                        sequences_loaded += 1;
                    }
                }
                
                let count = processed.fetch_add(1, Ordering::Relaxed) + 1;
                pb.set_position(count as u64);
                
                Ok((locus_name.clone(), locus_sequences, sequences_loaded))
            })
            .collect::<Result<Vec<_>, String>>()?;
        
        // Merge results into self
        for (locus_name, sequences, _count) in results {
            for (crc32, seq_info) in sequences {
                self.add_sequence(locus_name.clone(), crc32, seq_info);
            }
            // Don't add count here - add_sequence already increments total_sequences
        }
        
        Ok(())
    }
    
    
    /// Load from schema directory (each .fasta file is a locus)
    fn load_from_directory(&mut self, schema_dir: &Path) -> Result<(), String> {
        let entries = std::fs::read_dir(schema_dir)
            .map_err(|e| format!("Failed to read schema directory: {}", e))?;
        
        for entry in entries {
            let entry = entry.map_err(|e| format!("Failed to read directory entry: {}", e))?;
            let path = entry.path();
            
            if let Some(extension) = path.extension() {
                if extension == "fasta" || extension == "fa" {
                    if let Some(file_stem) = path.file_stem() {
                        let locus_name = file_stem.to_string_lossy().to_string();
                        self.load_fasta_into_locus(&path, &locus_name)?;
                    }
                }
            }
        }
        
        Ok(())
    }
    
    /// Load from schema file (tab-separated: locus_name<tab>fasta_path)
    fn load_from_schema_file(&mut self, schema_path: &Path) -> Result<(), String> {
        let content = std::fs::read_to_string(schema_path)
            .map_err(|e| format!("Failed to read schema file: {}", e))?;
        
        for (line_num, line) in content.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() != 2 {
                return Err(format!("Invalid schema line {}: expected locus_name<TAB>fasta_path", line_num + 1));
            }
            
            let locus_name = parts[0].to_string();
            let fasta_path = PathBuf::from(parts[1]);
            
            if !fasta_path.exists() {
                return Err(format!("FASTA file not found: {}", fasta_path.display()));
            }
            
            self.load_fasta_into_locus(&fasta_path, &locus_name)?;
        }
        
        Ok(())
    }
    
    /// Load FASTA file into a specific locus
    fn load_fasta_into_locus(&mut self, fasta_path: &Path, locus_name: &str) -> Result<(), String> {
        let file = File::open(fasta_path)
            .map_err(|e| format!("Failed to open FASTA file {}: {}", fasta_path.display(), e))?;
        
        let reader = fasta::Reader::new(BufReader::new(file));
        let mut sequences_added = 0;
        
        for record_result in reader.records() {
            let record = record_result
                .map_err(|e| format!("Invalid FASTA record in {}: {}", fasta_path.display(), e))?;
            
            let sequence = record.seq().to_vec();
            let id = record.id().to_string();
            let crc32 = compute_crc32(&sequence);
            
            // Check for CRC32 collision
            if let Some(existing_info) = self.get_sequence(locus_name, crc32) {
                if existing_info.sequence != sequence {
                    return Err(format!("CRC32 collision detected in {}: CRC32 {} maps to different sequences ({} vs {})", 
                                     fasta_path.display(), crc32, existing_info.id, id));
                }
            }
            
            let seq_info = SequenceInfo { sequence, id };
            self.add_sequence(locus_name.to_string(), crc32, seq_info);
            sequences_added += 1;
        }
        
        println!("  📄 {}: {} sequences loaded", locus_name, sequences_added);
        Ok(())
    }
}