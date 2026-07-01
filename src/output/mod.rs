// mod.rs - Output formatters module

pub mod ci;
pub mod report;
pub use report::write_html_report;

use crate::core::PairRow;
use crate::data::AllelicProfile;
use chrono;
use std::fs::{create_dir_all, File};
use std::io::{BufWriter, Write};
use std::path::Path;

/// Ensure parent directory exists before creating file
fn ensure_parent_dir(file_path: &str) -> Result<(), String> {
    if let Some(parent) = Path::new(file_path).parent() {
        create_dir_all(parent).map_err(|e| {
            format!(
                "Failed to create parent directory '{}': {}",
                parent.display(),
                e
            )
        })?;
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
        .map_err(|e| format!("Failed to create output file '{file_path}': {e}"))?;
    let mut writer = BufWriter::new(file);

    // Write command header
    writeln!(writer, "# Command: {command_line}").map_err(|e| format!("Write error: {e}"))?;
    writeln!(
        writer,
        "# Generated: {}",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC")
    )
    .map_err(|e| format!("Write error: {e}"))?;
    writeln!(writer, "# cgDist v{}", env!("CARGO_PKG_VERSION"))
        .map_err(|e| format!("Write error: {e}"))?;

    // Write header
    write!(writer, "Sample").map_err(|e| format!("Write error: {e}"))?;
    for sample in samples {
        write!(writer, "\t{}", sample.sample_id).map_err(|e| format!("Write error: {e}"))?;
    }
    writeln!(writer).map_err(|e| format!("Write error: {e}"))?;

    // Write matrix
    for (i, sample) in samples.iter().enumerate() {
        write!(writer, "{}", sample.sample_id).map_err(|e| format!("Write error: {e}"))?;
        for cell in matrix[i].iter().take(samples.len()) {
            let distance_str = match cell {
                Some(d) => d.to_string(),
                None => "NA".to_string(),
            };
            write!(writer, "\t{distance_str}").map_err(|e| format!("Write error: {e}"))?;
        }
        writeln!(writer).map_err(|e| format!("Write error: {e}"))?;
    }

    writer.flush().map_err(|e| format!("Flush error: {e}"))?;
    println!("✅ Distance matrix written to: {file_path}");
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
        .map_err(|e| format!("Failed to create output file '{file_path}': {e}"))?;
    let mut writer = BufWriter::new(file);

    // Write command header
    writeln!(writer, "# Command: {command_line}").map_err(|e| format!("Write error: {e}"))?;
    writeln!(
        writer,
        "# Generated: {}",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC")
    )
    .map_err(|e| format!("Write error: {e}"))?;
    writeln!(writer, "# cgDist v{}", env!("CARGO_PKG_VERSION"))
        .map_err(|e| format!("Write error: {e}"))?;

    // Write header
    write!(writer, "Sample").map_err(|e| format!("Write error: {e}"))?;
    for sample in samples {
        write!(writer, ",{}", sample.sample_id).map_err(|e| format!("Write error: {e}"))?;
    }
    writeln!(writer).map_err(|e| format!("Write error: {e}"))?;

    // Write matrix
    for (i, sample) in samples.iter().enumerate() {
        write!(writer, "{}", sample.sample_id).map_err(|e| format!("Write error: {e}"))?;
        for cell in matrix[i].iter().take(samples.len()) {
            let distance_str = match cell {
                Some(d) => d.to_string(),
                None => "NA".to_string(),
            };
            write!(writer, ",{distance_str}").map_err(|e| format!("Write error: {e}"))?;
        }
        writeln!(writer).map_err(|e| format!("Write error: {e}"))?;
    }

    writer.flush().map_err(|e| format!("Flush error: {e}"))?;
    println!("✅ Distance matrix written to: {file_path}");
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
        .map_err(|e| format!("Failed to create output file '{file_path}': {e}"))?;
    let mut writer = BufWriter::new(file);

    // PHYLIP doesn't support comments, but we can add them as the first "sample"
    // Actually, let's put comments at the end to maintain PHYLIP compatibility

    // Write header
    writeln!(writer, "    {}", samples.len()).map_err(|e| format!("Write error: {e}"))?;

    // Write matrix (lower triangle)
    for (i, sample) in samples.iter().enumerate() {
        write!(writer, "{:<10}", sample.sample_id).map_err(|e| format!("Write error: {e}"))?;
        for cell in matrix[i].iter().take(i + 1) {
            let distance_str = match cell {
                Some(d) => format!("  {d}"),
                None => "  NA".to_string(),
            };
            write!(writer, "{distance_str}").map_err(|e| format!("Write error: {e}"))?;
        }
        writeln!(writer).map_err(|e| format!("Write error: {e}"))?;
    }

    // Add command info as comments at the end (some PHYLIP parsers ignore trailing content)
    writeln!(writer).map_err(|e| format!("Write error: {e}"))?;
    writeln!(writer, "# Command: {command_line}").map_err(|e| format!("Write error: {e}"))?;
    writeln!(
        writer,
        "# Generated: {}",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC")
    )
    .map_err(|e| format!("Write error: {e}"))?;
    writeln!(writer, "# cgDist v{}", env!("CARGO_PKG_VERSION"))
        .map_err(|e| format!("Write error: {e}"))?;

    writer.flush().map_err(|e| format!("Flush error: {e}"))?;
    println!("✅ Distance matrix written to: {file_path} (PHYLIP format)");
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
        .map_err(|e| format!("Failed to create output file '{file_path}': {e}"))?;
    let mut writer = BufWriter::new(file);

    // Write NEXUS header with command info
    writeln!(writer, "#NEXUS").map_err(|e| format!("Write error: {e}"))?;
    writeln!(writer, "[Command: {command_line}]").map_err(|e| format!("Write error: {e}"))?;
    writeln!(
        writer,
        "[Generated: {}]",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC")
    )
    .map_err(|e| format!("Write error: {e}"))?;
    writeln!(writer, "[cgDist v{}]", env!("CARGO_PKG_VERSION"))
        .map_err(|e| format!("Write error: {e}"))?;
    writeln!(writer, "BEGIN DISTANCES;").map_err(|e| format!("Write error: {e}"))?;
    writeln!(writer, "    DIMENSIONS NTAX={};", samples.len())
        .map_err(|e| format!("Write error: {e}"))?;
    writeln!(writer, "    FORMAT LABELS LOWER DIAGONAL;")
        .map_err(|e| format!("Write error: {e}"))?;
    writeln!(writer, "    MATRIX").map_err(|e| format!("Write error: {e}"))?;

    // Write matrix (lower triangle)
    for (i, sample) in samples.iter().enumerate() {
        write!(writer, "        {}", sample.sample_id).map_err(|e| format!("Write error: {e}"))?;
        for cell in matrix[i].iter().take(i) {
            let distance_str = match cell {
                Some(d) => d.to_string(),
                None => "?".to_string(),
            };
            write!(writer, " {distance_str}").map_err(|e| format!("Write error: {e}"))?;
        }
        writeln!(writer).map_err(|e| format!("Write error: {e}"))?;
    }

    writeln!(writer, "    ;").map_err(|e| format!("Write error: {e}"))?;
    writeln!(writer, "END;").map_err(|e| format!("Write error: {e}"))?;

    writer.flush().map_err(|e| format!("Flush error: {e}"))?;
    println!("✅ Distance matrix written to: {file_path} (NEXUS format)");
    Ok(())
}

/// Write a long-format per-pair table with data-quality columns.
///
/// Columns: `sample_i  sample_j  distance  shared_loci  total_loci  missing_frac`.
/// One row per upper-triangle pair. `distance` is `NA` when the pair failed the
/// `--min-loci` filter. `missing_frac = 1 - shared_loci / total_loci` is the
/// fraction of the schema unavailable for the pair — a distance resting on few
/// shared loci is less reliable. Opt-in via `--emit-pairs`; does not affect the
/// square-matrix output.
pub fn write_pairs_long(
    file_path: &str,
    samples: &[AllelicProfile],
    rows: &[PairRow],
    total_loci: usize,
    ci_level: Option<f64>,
    command_line: &str,
) -> Result<(), String> {
    ensure_parent_dir(file_path)?;
    let file = File::create(file_path)
        .map_err(|e| format!("Failed to create output file '{file_path}': {e}"))?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "# Command: {command_line}").map_err(|e| format!("Write error: {e}"))?;
    writeln!(
        writer,
        "# Generated: {}",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC")
    )
    .map_err(|e| format!("Write error: {e}"))?;
    writeln!(writer, "# cgDist v{}", env!("CARGO_PKG_VERSION"))
        .map_err(|e| format!("Write error: {e}"))?;
    let mut header =
        String::from("sample_i\tsample_j\tdistance\tshared_loci\ttotal_loci\tmissing_frac");
    if ci_level.is_some() {
        header.push_str("\tdist_norm\tci_low\tci_high\tci_reliable");
    }
    writeln!(writer, "{header}").map_err(|e| format!("Write error: {e}"))?;

    let denom = total_loci.max(1) as f64;
    for row in rows {
        let dist_str = match row.distance {
            Some(d) => d.to_string(),
            None => "NA".to_string(),
        };
        let missing_frac = 1.0 - (row.shared as f64) / denom;
        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{:.4}",
            samples[row.i].sample_id,
            samples[row.j].sample_id,
            dist_str,
            row.shared,
            total_loci,
            missing_frac
        )
        .map_err(|e| format!("Write error: {e}"))?;
        if let (Some(level), Some(d)) = (ci_level, row.distance) {
            let (dn, lo, hi, reliable) = ci::pair_ci(d, row.h, row.q2, row.shared, total_loci, level);
            write!(writer, "\t{dn:.2}\t{lo:.2}\t{hi:.2}\t{reliable}")
                .map_err(|e| format!("Write error: {e}"))?;
        } else if ci_level.is_some() {
            write!(writer, "\tNA\tNA\tNA\tNA").map_err(|e| format!("Write error: {e}"))?;
        }
        writeln!(writer).map_err(|e| format!("Write error: {e}"))?;
    }

    writer.flush().map_err(|e| format!("Flush error: {e}"))?;
    println!("✅ Per-pair quality table written to: {file_path}");
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
        _ => Err(format!(
            "Unsupported output format: {format}. Use: tsv, csv, phylip, nexus"
        )),
    }
}
