// Improved distance calculation with streaming recombination detection and bidirectional alignment

use std::io::Write;
use csv::Writer;
use std::fs::File;

/// Improved calculate_sample_distance with bidirectional alignment for recombination
pub fn calculate_sample_distance_with_bidirectional_recombination(
    sample1: &AllelicProfile,
    sample2: &AllelicProfile,
    loci_names: &[String],
    engine: &DistanceEngine,
    mode: DistanceMode,
    min_loci: usize,
    no_hamming_fallback: bool,
    recomb_writer: Option<&mut Writer<File>>, // Pass writer directly for streaming
) -> Option<usize> {
    let mut total_distance = 0;
    let mut shared_loci = 0;
    let mut recomb_count = 0;
    
    for locus in loci_names {
        let crc1 = sample1.loci_hashes.get(locus)
            .and_then(|h| h.as_crc32())
            .unwrap_or(u32::MAX);
        let crc2 = sample2.loci_hashes.get(locus)
            .and_then(|h| h.as_crc32())
            .unwrap_or(u32::MAX);
        
        if crc1 != u32::MAX && crc2 != u32::MAX {
            shared_loci += 1;
            
            // Check for potential recombination events with BIDIRECTIONAL alignment
            if let Some(threshold) = engine.recombination_threshold {
                if crc1 != crc2 { // Only check different alleles
                    // First try forward alignment
                    let distance_forward = engine.get_distance(
                        locus, crc1, crc2, DistanceMode::SnpsAndIndelBases, no_hamming_fallback
                    );
                    
                    let final_distance = if distance_forward > threshold {
                        // Try reverse alignment if forward exceeds threshold
                        let distance_reverse = engine.get_distance(
                            locus, crc2, crc1, DistanceMode::SnpsAndIndelBases, no_hamming_fallback
                        );
                        
                        if distance_reverse <= threshold {
                            // Use reverse distance if it's within threshold
                            distance_reverse
                        } else {
                            // Both directions exceed threshold = recombination event
                            if let Some(writer) = recomb_writer {
                                // STREAM WRITE immediately instead of accumulating
                                write_recombination_event(
                                    writer,
                                    locus,
                                    sample1,
                                    sample2,
                                    crc1,
                                    crc2,
                                    distance_forward.min(distance_reverse),
                                    threshold,
                                    engine,
                                )?;
                                recomb_count += 1;
                                
                                // Flush every 100 events to avoid buffering too much
                                if recomb_count % 100 == 0 {
                                    writer.flush()?;
                                }
                            }
                            distance_forward.min(distance_reverse) // Use minimum for distance calculation
                        }
                    } else {
                        // Forward distance is within threshold, use it
                        distance_forward
                    };
                    
                    // Use the final distance for total calculation
                    total_distance += final_distance;
                } else {
                    // Same alleles, distance is 0
                    total_distance += 0;
                }
            } else {
                // No recombination detection, use normal distance
                total_distance += engine.get_distance(locus, crc1, crc2, mode, no_hamming_fallback);
            }
        }
    }
    
    if shared_loci >= min_loci {
        Some(total_distance)
    } else {
        None
    }
}

/// Helper function to write recombination event immediately
fn write_recombination_event(
    writer: &mut Writer<File>,
    locus: &str,
    sample1: &AllelicProfile,
    sample2: &AllelicProfile,
    crc1: u32,
    crc2: u32,
    distance: usize,
    threshold: usize,
    engine: &DistanceEngine,
) -> Result<(), Box<dyn std::error::Error>> {
    // Get sequence lengths if available
    let (seq_len1, seq_len2) = if let Some(ref db) = engine.sequence_db {
        let len1 = db.get_sequence(locus, crc1).map(|s| s.sequence.len()).unwrap_or(0);
        let len2 = db.get_sequence(locus, crc2).map(|s| s.sequence.len()).unwrap_or(0);
        (len1, len2)
    } else {
        (0, 0)
    };
    
    let avg_length = if seq_len1 > 0 && seq_len2 > 0 {
        (seq_len1 + seq_len2) / 2
    } else {
        450 // Default MLST locus length
    };
    
    let divergence_percent = if avg_length > 0 {
        (distance as f64 / avg_length as f64) * 100.0
    } else {
        0.0
    };
    
    // Write directly to CSV
    writer.write_record(&[
        locus,
        &sample1.sample_id,
        &sample2.sample_id,
        &format!("{}", crc1),
        &format!("{}", crc2),
        &distance.to_string(),
        &threshold.to_string(),
        &seq_len1.to_string(),
        &seq_len2.to_string(),
        &format!("{:.2}", divergence_percent),
    ])?;
    
    Ok(())
}

/// New compute_distance_matrix that streams recombination events
pub fn compute_distance_matrix_streaming(
    matrix: &AllelicMatrix,
    engine: &DistanceEngine,
    mode: DistanceMode,
    min_loci: usize,
    no_hamming_fallback: bool,
    recomb_log_path: Option<&str>,
) -> Result<Vec<Vec<Option<usize>>>, Box<dyn std::error::Error>> {
    let n = matrix.samples.len();
    let mut distance_matrix = vec![vec![None; n]; n];
    
    // Create CSV writer if recombination log is requested
    let mut recomb_writer = if let Some(path) = recomb_log_path {
        let mut w = Writer::from_path(path)?;
        // Write header
        w.write_record(&[
            "locus", "sample1", "sample2", "allele1_hash", "allele2_hash", 
            "snps_indel_bases", "threshold", "seq_length1", "seq_length2", 
            "divergence_percent"
        ])?;
        Some(w)
    } else {
        None
    };
    
    // Process upper triangle with progress bar
    let total_comparisons = n * (n - 1) / 2;
    let pb = ProgressBar::new(total_comparisons as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({eta})")
        .progress_chars("█▓▒░ "));
    
    let mut comparisons_done = 0;
    let update_interval = (total_comparisons / 100).max(1);
    
    for i in 0..n {
        for j in i..n {
            if i == j {
                distance_matrix[i][j] = Some(0);
            } else {
                let distance = calculate_sample_distance_with_bidirectional_recombination(
                    &matrix.samples[i],
                    &matrix.samples[j],
                    &matrix.loci_names,
                    engine,
                    mode,
                    min_loci,
                    no_hamming_fallback,
                    recomb_writer.as_mut(),
                );
                
                distance_matrix[i][j] = distance;
                distance_matrix[j][i] = distance; // Symmetric
                
                comparisons_done += 1;
                if comparisons_done % update_interval == 0 {
                    pb.set_position(comparisons_done as u64);
                }
            }
        }
    }
    
    pb.finish();
    
    // Final flush and close writer
    if let Some(mut w) = recomb_writer {
        w.flush()?;
        println!("🔬 Recombination events written progressively to log");
    }
    
    Ok(distance_matrix)
}