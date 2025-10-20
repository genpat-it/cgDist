// csv.rs - CSV file loader for allelic profiles

use crate::data::profile::{AllelicMatrix, AllelicProfile};
use crate::hashers::{AlleleHash, AlleleHasher};
use std::collections::HashMap;
use std::path::Path;

/// Parse allele string as CRC32
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
    /// Load CSV file with hasher support  
    pub fn from_csv_with_hasher(
        file_path: &Path,
        missing_char: &str,
        hasher: &dyn AlleleHasher,
    ) -> Result<Self, String> {
        let content = std::fs::read_to_string(file_path)
            .map_err(|e| format!("Failed to read CSV file: {}", e))?;

        let mut lines = content.lines();
        let header_line = lines.next().ok_or("Empty CSV file")?;
        let header_parts: Vec<&str> = header_line.split(',').collect();

        if header_parts.len() < 2 {
            return Err("CSV header must have at least 2 columns".to_string());
        }

        let loci_names: Vec<String> = header_parts[1..]
            .iter()
            .map(|s| s.trim().trim_matches('"').to_string())
            .collect();
        let mut samples = Vec::new();

        for (line_num, line) in lines.enumerate() {
            let parts: Vec<&str> = line.split(',').collect();

            if parts.len() != header_parts.len() {
                return Err(format!(
                    "CSV line {} has {} columns, expected {}",
                    line_num + 2,
                    parts.len(),
                    header_parts.len()
                ));
            }

            let sample_id = parts[0].trim().trim_matches('"').to_string();
            let mut loci_hashes = HashMap::new();

            for (i, &allele_str) in parts[1..].iter().enumerate() {
                let cleaned = allele_str.trim().trim_matches('"');
                let hash = hasher.parse_allele(cleaned, missing_char).map_err(|e| {
                    format!(
                        "Invalid allele '{}' at CSV line {} locus {}: {}",
                        cleaned,
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
            "✅ CSV loaded: {} samples, {} loci",
            samples.len(),
            loci_names.len()
        );
        Ok(Self {
            samples,
            loci_names,
        })
    }

    /// Legacy CSV loader for CRC32 compatibility  
    pub fn from_csv(file_path: &Path, missing_char: &str) -> Result<Self, String> {
        let content = std::fs::read_to_string(file_path)
            .map_err(|e| format!("Failed to read CSV file: {}", e))?;

        let mut lines = content.lines();
        let header_line = lines.next().ok_or("Empty CSV file")?;
        let header_parts: Vec<&str> = header_line.split(',').collect();

        if header_parts.len() < 2 {
            return Err("CSV header must have at least sample_id and one locus".to_string());
        }

        let loci_names: Vec<String> = header_parts[1..]
            .iter()
            .map(|s| s.trim().trim_matches('"').to_string())
            .collect();
        let mut samples = Vec::new();

        for (line_num, line) in lines.enumerate() {
            let parts: Vec<&str> = line.split(',').collect();
            if parts.len() != header_parts.len() {
                return Err(format!(
                    "CSV line {} has {} fields, expected {}",
                    line_num + 2,
                    parts.len(),
                    header_parts.len()
                ));
            }

            let sample_id = parts[0].trim().trim_matches('"').to_string();
            let mut loci_hashes = HashMap::new();

            for (i, &allele_str) in parts[1..].iter().enumerate() {
                let cleaned = allele_str.trim().trim_matches('"');
                let crc32 = parse_allele_crc32(cleaned, missing_char).map_err(|e| {
                    format!(
                        "Invalid allele '{}' at CSV line {} locus {}: {}",
                        cleaned,
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
            "✅ CSV matrix loaded: {} samples, {} loci",
            samples.len(),
            loci_names.len()
        );
        Ok(AllelicMatrix {
            samples,
            loci_names,
        })
    }
}
