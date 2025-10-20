// tsv.rs - TSV file loader for allelic profiles

use crate::data::profile::{AllelicMatrix, AllelicProfile};
use crate::hashers::{AlleleHash, AlleleHasher};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Specialized allele parsing for CRC32 hasher
fn parse_allele_crc32(s: &str, missing_char: &str) -> Result<u32, String> {
    let cleaned = s.trim();

    if cleaned.is_empty() || cleaned == "NA" || cleaned == missing_char {
        return Ok(u32::MAX); // MISSING_ALLELE
    }

    let crc = cleaned
        .parse::<u32>()
        .map_err(|_| format!("Failed to parse '{}' as CRC32", cleaned))?;

    if crc == u32::MAX {
        return Err("CRC32 value cannot be u32::MAX (reserved for missing)".to_string());
    }

    Ok(crc)
}

impl AllelicMatrix {
    /// Load TSV file with hasher support
    pub fn from_tsv_with_hasher(
        file_path: &Path,
        missing_char: &str,
        hasher: &dyn AlleleHasher,
    ) -> Result<Self, String> {
        let file = File::open(file_path).map_err(|e| format!("Failed to open TSV file: {}", e))?;

        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        // Read header
        let header_line = lines
            .next()
            .ok_or("Empty TSV file")?
            .map_err(|e| format!("Failed to read header: {}", e))?;

        let header_parts: Vec<&str> = header_line.split('\t').collect();
        if header_parts.len() < 2 {
            return Err("TSV header must have at least 2 columns".to_string());
        }

        let loci_names: Vec<String> = header_parts[1..].iter().map(|s| s.to_string()).collect();
        let mut samples = Vec::new();

        for (line_num, line) in lines.enumerate() {
            let line = line.map_err(|e| format!("Failed to read line {}: {}", line_num + 2, e))?;
            let parts: Vec<&str> = line.split('\t').collect();

            if parts.len() != header_parts.len() {
                return Err(format!(
                    "Line {} has {} columns, expected {}",
                    line_num + 2,
                    parts.len(),
                    header_parts.len()
                ));
            }

            let sample_id = parts[0].to_string();
            let mut loci_hashes = HashMap::new();

            for (i, &allele_str) in parts[1..].iter().enumerate() {
                let hash = hasher.parse_allele(allele_str, missing_char).map_err(|e| {
                    format!(
                        "Invalid allele '{}' at line {} locus {}: {}",
                        allele_str,
                        line_num + 2,
                        loci_names[i],
                        e
                    )
                })?;

                loci_hashes.insert(loci_names[i].clone(), hash);
            }

            samples.push(AllelicProfile {
                sample_id,
                loci_hashes,
            });
        }

        println!(
            "✅ TSV loaded: {} samples, {} loci",
            samples.len(),
            loci_names.len()
        );
        Ok(Self {
            samples,
            loci_names,
        })
    }

    /// Legacy TSV loader for CRC32 compatibility
    pub fn from_tsv(file_path: &Path, missing_char: &str) -> Result<Self, String> {
        let file = File::open(file_path).map_err(|e| format!("Failed to open TSV file: {}", e))?;

        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        let header_line = lines
            .next()
            .ok_or("Empty TSV file")?
            .map_err(|e| format!("Failed to read header: {}", e))?;

        let header_parts: Vec<&str> = header_line.split('\t').collect();
        if header_parts.len() < 2 {
            return Err("Header must have at least sample_id and one locus".to_string());
        }

        let loci_names: Vec<String> = header_parts[1..].iter().map(|s| s.to_string()).collect();
        let mut samples = Vec::new();

        for (line_num, line_result) in lines.enumerate() {
            let line =
                line_result.map_err(|e| format!("Failed to read line {}: {}", line_num + 2, e))?;

            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() != header_parts.len() {
                return Err(format!(
                    "Line {} has {} fields, expected {}",
                    line_num + 2,
                    parts.len(),
                    header_parts.len()
                ));
            }

            let sample_id = parts[0].to_string();
            let mut loci_hashes = HashMap::new();

            for (i, &allele_str) in parts[1..].iter().enumerate() {
                let crc32 = parse_allele_crc32(allele_str, missing_char).map_err(|e| {
                    format!(
                        "Invalid allele '{}' at line {} locus {}: {}",
                        allele_str,
                        line_num + 2,
                        loci_names[i],
                        e
                    )
                })?;

                loci_hashes.insert(loci_names[i].clone(), AlleleHash::from_crc32(crc32));
            }

            samples.push(AllelicProfile {
                sample_id,
                loci_hashes,
            });
        }

        println!(
            "✅ Matrix loaded: {} samples, {} loci",
            samples.len(),
            loci_names.len()
        );
        Ok(AllelicMatrix {
            samples,
            loci_names,
        })
    }
}
