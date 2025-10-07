// mod.rs - Output formatters module

use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::Path;
use crate::data::AllelicProfile;
use chrono;

/// Ensure parent directory exists before creating file
fn ensure_parent_dir(file_path: &str) -> Result<(), String> {
    if let Some(parent) = Path::new(file_path).parent() {
        create_dir_all(parent)
            .map_err(|e| format!("Failed to create parent directory '{}': {}", parent.display(), e))?;
    }
    Ok(())
}

/// Write distance matrix in TSV format
pub fn write_tsv(
    file_path: &str,
    samples: &[AllelicProfile],
    matrix: &[Vec<Option<usize>>],
    command_line: &str,
) -> Result<(), String> {
    ensure_parent_dir(file_path)?;
    let file = File::create(file_path)
        .map_err(|e| format!("Failed to create output file '{}': {}", file_path, e))?;
    let mut writer = BufWriter::new(file);
    
    // Write command header
    writeln!(writer, "# Command: {}", command_line).map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "# Generated: {}", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC")).map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "# cgDist v{}", env!("CARGO_PKG_VERSION")).map_err(|e| format!("Write error: {}", e))?;
    
    // Write header
    write!(writer, "Sample").map_err(|e| format!("Write error: {}", e))?;
    for sample in samples {
        write!(writer, "\t{}", sample.sample_id).map_err(|e| format!("Write error: {}", e))?;
    }
    writeln!(writer).map_err(|e| format!("Write error: {}", e))?;
    
    // Write matrix
    for (i, sample) in samples.iter().enumerate() {
        write!(writer, "{}", sample.sample_id).map_err(|e| format!("Write error: {}", e))?;
        for j in 0..samples.len() {
            let distance_str = match matrix[i][j] {
                Some(d) => d.to_string(),
                None => "NA".to_string(),
            };
            write!(writer, "\t{}", distance_str).map_err(|e| format!("Write error: {}", e))?;
        }
        writeln!(writer).map_err(|e| format!("Write error: {}", e))?;
    }
    
    writer.flush().map_err(|e| format!("Flush error: {}", e))?;
    println!("✅ Distance matrix written to: {}", file_path);
    Ok(())
}

/// Write distance matrix in CSV format
pub fn write_csv(
    file_path: &str,
    samples: &[AllelicProfile],
    matrix: &[Vec<Option<usize>>],
    command_line: &str,
) -> Result<(), String> {
    ensure_parent_dir(file_path)?;
    let file = File::create(file_path)
        .map_err(|e| format!("Failed to create output file '{}': {}", file_path, e))?;
    let mut writer = BufWriter::new(file);
    
    // Write command header
    writeln!(writer, "# Command: {}", command_line).map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "# Generated: {}", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC")).map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "# cgDist v{}", env!("CARGO_PKG_VERSION")).map_err(|e| format!("Write error: {}", e))?;
    
    // Write header
    write!(writer, "Sample").map_err(|e| format!("Write error: {}", e))?;
    for sample in samples {
        write!(writer, ",{}", sample.sample_id).map_err(|e| format!("Write error: {}", e))?;
    }
    writeln!(writer).map_err(|e| format!("Write error: {}", e))?;
    
    // Write matrix
    for (i, sample) in samples.iter().enumerate() {
        write!(writer, "{}", sample.sample_id).map_err(|e| format!("Write error: {}", e))?;
        for j in 0..samples.len() {
            let distance_str = match matrix[i][j] {
                Some(d) => d.to_string(),
                None => "NA".to_string(),
            };
            write!(writer, ",{}", distance_str).map_err(|e| format!("Write error: {}", e))?;
        }
        writeln!(writer).map_err(|e| format!("Write error: {}", e))?;
    }
    
    writer.flush().map_err(|e| format!("Flush error: {}", e))?;
    println!("✅ Distance matrix written to: {}", file_path);
    Ok(())
}

/// Write distance matrix in PHYLIP format
pub fn write_phylip(
    file_path: &str,
    samples: &[AllelicProfile],
    matrix: &[Vec<Option<usize>>],
    command_line: &str,
) -> Result<(), String> {
    ensure_parent_dir(file_path)?;
    let file = File::create(file_path)
        .map_err(|e| format!("Failed to create output file '{}': {}", file_path, e))?;
    let mut writer = BufWriter::new(file);
    
    // PHYLIP doesn't support comments, but we can add them as the first "sample"
    // Actually, let's put comments at the end to maintain PHYLIP compatibility
    
    // Write header
    writeln!(writer, "    {}", samples.len()).map_err(|e| format!("Write error: {}", e))?;
    
    // Write matrix (lower triangle)
    for (i, sample) in samples.iter().enumerate() {
        write!(writer, "{:<10}", sample.sample_id).map_err(|e| format!("Write error: {}", e))?;
        for j in 0..=i {
            let distance_str = match matrix[i][j] {
                Some(d) => format!("  {}", d),
                None => "  NA".to_string(),
            };
            write!(writer, "{}", distance_str).map_err(|e| format!("Write error: {}", e))?;
        }
        writeln!(writer).map_err(|e| format!("Write error: {}", e))?;
    }
    
    // Add command info as comments at the end (some PHYLIP parsers ignore trailing content)
    writeln!(writer).map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "# Command: {}", command_line).map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "# Generated: {}", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC")).map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "# cgDist v{}", env!("CARGO_PKG_VERSION")).map_err(|e| format!("Write error: {}", e))?;
    
    writer.flush().map_err(|e| format!("Flush error: {}", e))?;
    println!("✅ Distance matrix written to: {} (PHYLIP format)", file_path);
    Ok(())
}

/// Write distance matrix in NEXUS format
pub fn write_nexus(
    file_path: &str,
    samples: &[AllelicProfile],
    matrix: &[Vec<Option<usize>>],
    command_line: &str,
) -> Result<(), String> {
    ensure_parent_dir(file_path)?;
    let file = File::create(file_path)
        .map_err(|e| format!("Failed to create output file '{}': {}", file_path, e))?;
    let mut writer = BufWriter::new(file);
    
    // Write NEXUS header with command info
    writeln!(writer, "#NEXUS").map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "[Command: {}]", command_line).map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "[Generated: {}]", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC")).map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "[cgDist v{}]", env!("CARGO_PKG_VERSION")).map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "BEGIN DISTANCES;").map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "    DIMENSIONS NTAX={};", samples.len()).map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "    FORMAT LABELS LOWER DIAGONAL;").map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "    MATRIX").map_err(|e| format!("Write error: {}", e))?;
    
    // Write matrix (lower triangle)
    for (i, sample) in samples.iter().enumerate() {
        write!(writer, "        {}", sample.sample_id).map_err(|e| format!("Write error: {}", e))?;
        for j in 0..i {
            let distance_str = match matrix[i][j] {
                Some(d) => d.to_string(),
                None => "?".to_string(),
            };
            write!(writer, " {}", distance_str).map_err(|e| format!("Write error: {}", e))?;
        }
        writeln!(writer).map_err(|e| format!("Write error: {}", e))?;
    }
    
    writeln!(writer, "    ;").map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "END;").map_err(|e| format!("Write error: {}", e))?;
    
    writer.flush().map_err(|e| format!("Flush error: {}", e))?;
    println!("✅ Distance matrix written to: {} (NEXUS format)", file_path);
    Ok(())
}

/// Write distance matrix in the specified format
pub fn write_matrix(
    file_path: &str,
    format: &str,
    samples: &[AllelicProfile],
    matrix: &[Vec<Option<usize>>],
    command_line: &str,
) -> Result<(), String> {
    match format.to_lowercase().as_str() {
        "tsv" => write_tsv(file_path, samples, matrix, command_line),
        "csv" => write_csv(file_path, samples, matrix, command_line),
        "phylip" => write_phylip(file_path, samples, matrix, command_line),
        "nexus" => write_nexus(file_path, samples, matrix, command_line),
        _ => Err(format!("Unsupported output format: {}. Use: tsv, csv, phylip, nexus", format)),
    }
}