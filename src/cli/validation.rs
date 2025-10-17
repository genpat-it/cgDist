// validation.rs - Input validation utilities

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufReader, BufRead};
use std::str::FromStr;
use regex::Regex;
use crate::cli::args::Args;
use crate::hashers::HasherRegistry;
use crate::core::{DistanceMode, AlignmentConfig};

pub struct ValidationResult {
    pub distance_mode: DistanceMode,
    pub alignment_config: AlignmentConfig,
    pub sample_include_regex: Option<Regex>,
    pub sample_exclude_regex: Option<Regex>,
    pub loci_include_regex: Option<Regex>,
    pub loci_exclude_regex: Option<Regex>,
    pub samples_include_set: Option<HashSet<String>>,
    pub samples_exclude_set: Option<HashSet<String>>,
    pub loci_include_set: Option<HashSet<String>>,
    pub loci_exclude_set: Option<HashSet<String>>,
}

/// Validate all command line arguments
pub fn validate_args(args: &Args) -> Result<ValidationResult, String> {
    // Validate hasher type
    let registry = HasherRegistry::new();
    if !registry.has_hasher(&args.hasher_type) {
        return Err(format!("Invalid hasher type '{}'. Available: {}", 
                          args.hasher_type,
                          registry.get_hasher_names().join(", ")));
    }

    // Validate hamming hasher incompatibilities
    if args.hasher_type == "hamming" {
        if args.schema.is_some() {
            return Err("--schema is not compatible with --hasher-type hamming (hamming hasher works at allelic level only)".to_string());
        }
        if args.cache_file.is_some() {
            return Err("--cache-file is not compatible with --hasher-type hamming (no alignment caching needed)".to_string());
        }
        if args.cache_only {
            return Err("--cache-only is not compatible with --hasher-type hamming (no alignment caching needed)".to_string());
        }
        if args.mode != "snps" {
            return Err(format!("Distance mode '{}' is not compatible with --hasher-type hamming (hamming works at allelic level)", args.mode));
        }
    }
    
    // Validate cache-only mode
    if args.cache_only && args.cache_file.is_none() {
        return Err("--cache-only requires --cache-file to specify where to save the cache".to_string());
    }

    // Validate distance mode
    let distance_mode = DistanceMode::from_str(&args.mode)?;

    // Validate and create alignment config (skip for hamming hasher)
    let alignment_config = if args.hasher_type == "hamming" {
        // Hamming hasher doesn't need alignment config
        AlignmentConfig::default()
    } else if args.match_score.is_some() || args.mismatch_penalty.is_some() ||
              args.gap_open.is_some() || args.gap_extend.is_some() {
        // Custom mode
        AlignmentConfig::custom(
            args.match_score.unwrap_or(2),
            args.mismatch_penalty.unwrap_or(-1),
            args.gap_open.unwrap_or(5),
            args.gap_extend.unwrap_or(2),
        )
    } else {
        // Preset mode
        AlignmentConfig::from_mode(&args.alignment_mode)?
    };

    // Early cache compatibility check (before loading any files) - using quick heuristic check
    if args.hasher_type != "hamming" {
        if let Some(ref cache_path) = args.cache_file {
            if std::path::Path::new(cache_path).exists() && !args.force_recompute {
                // Quick heuristic compatibility check without loading full data
                use crate::core::DistanceEngine;
                let temp_engine = DistanceEngine::new(alignment_config.clone(), args.hasher_type.clone());
                if let Err(e) = temp_engine.check_cache_compatibility(cache_path, distance_mode) {
                    eprintln!("❌ FATAL ERROR: Cache incompatible: {}", e);
                    eprintln!("💡 Solutions:");
                    eprintln!("   - Use --force-recompute to ignore existing cache");
                    eprintln!("   - Delete incompatible cache file");
                    eprintln!("   - Use consistent alignment parameters");
                    std::process::exit(1);
                }
            }
        }
    }

    // Validate thresholds
    if args.sample_threshold < 0.0 || args.sample_threshold > 1.0 {
        return Err("Sample threshold must be between 0.0 and 1.0".to_string());
    }
    if args.locus_threshold < 0.0 || args.locus_threshold > 1.0 {
        return Err("Locus threshold must be between 0.0 and 1.0".to_string());
    }

    // Compile regex patterns
    let sample_include_regex = if let Some(pattern) = &args.include_samples {
        Some(Regex::new(pattern).map_err(|e| format!("Invalid include_samples regex: {}", e))?)
    } else {
        None
    };

    let sample_exclude_regex = if let Some(pattern) = &args.exclude_samples {
        Some(Regex::new(pattern).map_err(|e| format!("Invalid exclude_samples regex: {}", e))?)
    } else {
        None
    };

    let loci_include_regex = if let Some(pattern) = &args.include_loci {
        Some(Regex::new(pattern).map_err(|e| format!("Invalid include_loci regex: {}", e))?)
    } else {
        None
    };

    let loci_exclude_regex = if let Some(pattern) = &args.exclude_loci {
        Some(Regex::new(pattern).map_err(|e| format!("Invalid exclude_loci regex: {}", e))?)
    } else {
        None
    };

    // Load filter sets from files
    let samples_include_set = if let Some(file_path) = &args.include_samples_list {
        Some(load_set_from_file(file_path)?)
    } else {
        None
    };

    let samples_exclude_set = if let Some(file_path) = &args.exclude_samples_list {
        Some(load_set_from_file(file_path)?)
    } else {
        None
    };

    let loci_include_set = if let Some(file_path) = &args.include_loci_list {
        Some(load_set_from_file(file_path)?)
    } else {
        None
    };

    let loci_exclude_set = if let Some(file_path) = &args.exclude_loci_list {
        Some(load_set_from_file(file_path)?)
    } else {
        None
    };

    Ok(ValidationResult {
        distance_mode,
        alignment_config,
        sample_include_regex,
        sample_exclude_regex,
        loci_include_regex,
        loci_exclude_regex,
        samples_include_set,
        samples_exclude_set,
        loci_include_set,
        loci_exclude_set,
    })
}

/// Load a set of strings from a file (one per line)
fn load_set_from_file(file_path: &str) -> Result<HashSet<String>, String> {
    let file = File::open(file_path)
        .map_err(|e| format!("Failed to open filter file '{}': {}", file_path, e))?;
    
    let reader = BufReader::new(file);
    let mut set = HashSet::new();
    
    for (line_num, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| format!("Failed to read line {} from '{}': {}", 
                                           line_num + 1, file_path, e))?;
        let trimmed = line.trim();
        if !trimmed.is_empty() {
            set.insert(trimmed.to_string());
        }
    }
    
    println!("📋 Loaded {} items from filter file '{}'", set.len(), file_path);
    Ok(set)
}