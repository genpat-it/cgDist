// sequences.rs - Sequence database for FASTA schema loading

use std::path::{Path, PathBuf};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use bio::io::fasta;
use indicatif::{ProgressBar, ProgressStyle};
use crc32fast::Hasher;

/// Sequence information
#[derive(Debug, Clone)]
pub struct SequenceInfo {
    pub sequence: Vec<u8>,
    pub id: String,
}

/// Sequence database for storing schema sequences
#[derive(Debug)]
pub struct SequenceDatabase {
    pub loci: HashMap<String, HashMap<u32, SequenceInfo>>,
    pub total_sequences: usize,
}

/// Compute CRC32 hash of a sequence
fn compute_crc32(sequence: &[u8]) -> u32 {
    let mut hasher = Hasher::new();
    hasher.update(sequence);
    hasher.finalize()
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
    
    pub fn load_from_schema(&mut self, schema_path: &Path) -> Result<(), String> {
        if schema_path.is_dir() {
            self.load_from_directory(schema_path)
        } else {
            self.load_from_schema_file(schema_path)
        }
    }
    
    fn load_from_directory(&mut self, fasta_dir: &Path) -> Result<(), String> {
        println!("🧬 Loading FASTA schema from directory: {}", fasta_dir.display());
        
        let entries: Vec<_> = std::fs::read_dir(fasta_dir)
            .map_err(|e| format!("Failed to read directory: {}", e))?
            .filter_map(|entry| entry.ok())
            .filter(|entry| {
                let path = entry.path();
                path.is_file() && matches!(
                    path.extension().and_then(|s| s.to_str()).unwrap_or(""),
                    "fasta" | "fa" | "fas" | "fna"
                )
            })
            .collect();

        if entries.is_empty() {
            return Err("No FASTA files found in directory".to_string());
        }

        let pb = ProgressBar::new(entries.len() as u64);
        pb.set_style(ProgressStyle::with_template(
            "[{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({percent}%) {msg}"
        ).unwrap());

        for entry in entries {
            let path = entry.path();
            let locus_name = path.file_stem()
                .and_then(|s| s.to_str())
                .ok_or_else(|| format!("Invalid filename: {}", path.display()))?
                .to_string();

            let mut locus_sequences = HashMap::new();
            self.load_fasta_into_locus(&path, &mut locus_sequences)?;
            
            self.total_sequences += locus_sequences.len();
            self.loci.insert(locus_name, locus_sequences);
            
            pb.inc(1);
            pb.set_message(format!("Loaded: {} loci", self.loci.len()));
        }

        pb.finish_with_message(format!("✅ Schema loaded: {} loci, {} sequences", 
                                     self.loci.len(), self.total_sequences));
        Ok(())
    }
    
    fn load_from_schema_file(&mut self, schema_path: &Path) -> Result<(), String> {
        println!("🧬 Loading FASTA schema from file: {}", schema_path.display());
        
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
            
            let mut locus_sequences = HashMap::new();
            self.load_fasta_into_locus(&fasta_path, &mut locus_sequences)?;
            
            self.total_sequences += locus_sequences.len();
            self.loci.insert(locus_name, locus_sequences);
        }
        
        println!("✅ Schema loaded: {} loci, {} sequences", self.loci.len(), self.total_sequences);
        Ok(())
    }
    
    fn load_fasta_into_locus(&self, fasta_path: &Path, locus_sequences: &mut HashMap<u32, SequenceInfo>) -> Result<(), String> {
        let file = File::open(fasta_path)
            .map_err(|e| format!("Failed to open FASTA file {}: {}", fasta_path.display(), e))?;
        
        let reader = fasta::Reader::new(BufReader::new(file));
        
        for record_result in reader.records() {
            let record = record_result
                .map_err(|e| format!("Invalid FASTA record in {}: {}", fasta_path.display(), e))?;
            
            let sequence = record.seq().to_vec();
            let id = record.id().to_string();
            let crc32 = compute_crc32(&sequence);
            
            if let Some(existing_info) = locus_sequences.get(&crc32) {
                if existing_info.sequence != sequence {
                    return Err(format!("CRC32 collision detected in {}: CRC32 {} maps to different sequences ({} vs {})", 
                                     fasta_path.display(), crc32, existing_info.id, id));
                }
            }
            
            let seq_info = SequenceInfo { sequence, id };
            locus_sequences.insert(crc32, seq_info);
        }
        
        Ok(())
    }
    
    pub fn get_sequence(&self, locus: &str, crc32: u32) -> Option<&SequenceInfo> {
        self.loci.get(locus)?.get(&crc32)
    }
}