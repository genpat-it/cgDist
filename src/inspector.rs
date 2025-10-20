// inspector.rs - Cache-aware database inspector with alignment config analysis
// Features: LZ4 cache inspection, alignment mode detection, compatibility checks

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::Path;

use argh::FromArgs;
use lz4_flex::decompress_size_prepended;
use serde::{Deserialize, Serialize};

// ============================================================================
// CORE DATA STRUCTURES (must match cgdist.rs exactly)
// ============================================================================

// Type aliases for complex types to satisfy clippy
type CacheKey = (String, u32, u32, u64);
type CacheValue = (usize, usize, usize);

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CrcPair {
    pub crc1: u32,
    pub crc2: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct AlignmentConfig {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub description: Option<String>,
}

impl AlignmentConfig {
    fn detect_mode(&self) -> String {
        match (
            self.match_score,
            self.mismatch_penalty,
            self.gap_open,
            self.gap_extend,
        ) {
            (2, -1, 5, 2) => "dna".to_string(),
            (3, -2, 8, 3) => "dna-strict".to_string(),
            (2, -1, 3, 1) => "dna-permissive".to_string(),
            _ => "custom".to_string(),
        }
    }

    fn get_description(&self) -> String {
        if let Some(desc) = &self.description {
            desc.clone()
        } else {
            match self.detect_mode().as_str() {
                "dna" => "Balanced DNA alignment parameters (BWA-MEM style)".to_string(),
                "dna-strict" => "Strict DNA alignment parameters (HOXD style)".to_string(),
                "dna-permissive" => "Permissive DNA alignment parameters".to_string(),
                "custom" => format!(
                    "Custom parameters: match={}, mismatch={}, gap_open={}, gap_extend={}",
                    self.match_score, self.mismatch_penalty, self.gap_open, self.gap_extend
                ),
                _ => "Unknown alignment mode".to_string(),
            }
        }
    }
}

// Ultra-fast cache structure (matches cgdist.rs)
#[derive(Debug, Serialize, Deserialize)]
pub struct UltraFastCache {
    pub data: HashMap<CacheKey, CacheValue>,
    pub version: String,
    pub created: String,
    pub alignment_config: AlignmentConfig,
}

// ============================================================================
// CLI ARGUMENTS
// ============================================================================

#[derive(FromArgs)]
/// Inspect cgdist cache files with alignment config analysis
struct Args {
    /// path to the cache file (.lz4)
    #[argh(option)]
    cache: String,

    /// show detailed locus information
    #[argh(switch)]
    detailed: bool,

    /// show cache entries for specific locus
    #[argh(option)]
    show_locus: Option<String>,

    /// detect alignment mode from parameters
    #[argh(switch)]
    detect_mode: bool,

    /// check compatibility with given alignment parameters (format: match,mismatch,gap_open,gap_extend)
    #[argh(option)]
    check_compatibility: Option<String>,

    /// export cache summary to TSV file
    #[argh(option)]
    export_summary: Option<String>,

    /// validate cache integrity
    #[argh(switch)]
    validate: bool,

    /// show top N loci by entry count (default: 10)
    #[argh(option, default = "10")]
    top_loci: usize,

    /// quiet mode - minimal output
    #[argh(switch)]
    quiet: bool,
}

// ============================================================================
// CACHE LOADING
// ============================================================================

fn load_cache(
    cache_path: &Path,
    quiet: bool,
) -> Result<UltraFastCache, Box<dyn std::error::Error>> {
    if !cache_path.exists() {
        return Err(format!("Cache file does not exist: {}", cache_path.display()).into());
    }

    if !quiet {
        println!("Loading cache: {}", cache_path.display());
    }

    let compressed_data = std::fs::read(cache_path)?;
    let decompressed = decompress_size_prepended(&compressed_data)
        .map_err(|e| format!("LZ4 decompression failed: {}", e))?;

    let cache: UltraFastCache = bincode::deserialize(&decompressed)?;

    if !quiet {
        println!("✅ Cache loaded successfully");
        let compression_ratio = compressed_data.len() as f64 / decompressed.len() as f64;
        println!(
            "Compression ratio: {:.2}:1 ({:.1}% space saved)",
            1.0 / compression_ratio,
            (1.0 - compression_ratio) * 100.0
        );
    }

    Ok(cache)
}

// ============================================================================
// ANALYSIS FUNCTIONS
// ============================================================================

fn analyze_cache_overview(cache: &UltraFastCache, args: &Args) {
    if args.quiet {
        return;
    }

    println!("\n=== CACHE SUMMARY ===");
    println!("Version: {}", cache.version);
    println!("Created: {}", cache.created);
    println!("Total entries: {}", cache.data.len());

    // Count unique loci and CRCs
    let mut loci = HashSet::new();
    let mut all_crcs = HashSet::new();

    for (locus, crc1, crc2, _) in cache.data.keys() {
        loci.insert(locus.clone());
        all_crcs.insert(*crc1);
        all_crcs.insert(*crc2);
    }

    println!("Unique loci: {}", loci.len());
    println!("Unique CRCs: {}", all_crcs.len());

    // Overall statistics
    let mut total_snps = 0;
    let mut total_indel_events = 0;
    let mut total_indel_bases = 0;

    for (snps, indel_events, indel_bases) in cache.data.values() {
        total_snps += snps;
        total_indel_events += indel_events;
        total_indel_bases += indel_bases;
    }

    let count = cache.data.len();
    if count > 0 {
        println!("\n=== DISTANCE STATISTICS ===");
        println!(
            "Average SNPs per pair: {:.2}",
            total_snps as f64 / count as f64
        );
        println!(
            "Average indel events per pair: {:.2}",
            total_indel_events as f64 / count as f64
        );
        println!(
            "Average indel bases per pair: {:.2}",
            total_indel_bases as f64 / count as f64
        );
    }
}

fn analyze_alignment_config(cache: &UltraFastCache, args: &Args) {
    if args.quiet && !args.detect_mode {
        return;
    }

    println!("\n=== ALIGNMENT CONFIGURATION ===");
    println!("Match score: {}", cache.alignment_config.match_score);
    println!(
        "Mismatch penalty: {}",
        cache.alignment_config.mismatch_penalty
    );
    println!("Gap open penalty: {}", cache.alignment_config.gap_open);
    println!("Gap extend penalty: {}", cache.alignment_config.gap_extend);

    if args.detect_mode {
        let detected_mode = cache.alignment_config.detect_mode();
        println!("Detected mode: {}", detected_mode);
        println!("Description: {}", cache.alignment_config.get_description());

        println!("\n=== PARAMETER INTERPRETATION ===");
        match detected_mode.as_str() {
            "dna" => {
                println!("🎯 Standard DNA mode (BWA-MEM style)");
                println!("   - Balanced parameters for typical bacterial genomes");
                println!("   - Good for core genome MLST analysis");
            }
            "dna-strict" => {
                println!("🎯 Strict DNA mode (HOXD style)");
                println!("   - Higher penalties for mismatches and gaps");
                println!("   - More conservative alignment scoring");
                println!("   - Good for high-quality, closely related sequences");
            }
            "dna-permissive" => {
                println!("🎯 Permissive DNA mode");
                println!("   - Lower gap penalties");
                println!("   - More tolerant of sequence differences");
                println!("   - Good for divergent sequences (~75% identity)");
            }
            "custom" => {
                println!("🎯 Custom alignment mode");
                println!("   - User-defined parameters");
                println!("   - Check if parameters are appropriate for your data");
            }
            _ => {}
        }
    }
}

fn analyze_loci_overview(cache: &UltraFastCache, args: &Args) {
    if args.quiet && !args.detailed {
        return;
    }

    println!("\n=== LOCI OVERVIEW ===");

    // Count entries per locus
    let mut locus_distribution = HashMap::new();
    for (locus, _, _, _) in cache.data.keys() {
        *locus_distribution.entry(locus.clone()).or_insert(0) += 1;
    }

    // Sort loci by entry count
    let mut loci_info: Vec<_> = locus_distribution.iter().collect();
    loci_info.sort_by(|a, b| b.1.cmp(a.1));

    println!(
        "{:<25} {:>8} {:>10} {:>12} {:>14}",
        "Locus", "Entries", "Avg SNPs", "Avg Indels", "Avg IndelBases"
    );
    println!("{}", "=".repeat(75));

    let show_count = if args.detailed {
        loci_info.len()
    } else {
        args.top_loci.min(loci_info.len())
    };

    for (locus, &entry_count) in loci_info.iter().take(show_count) {
        // Calculate averages for this locus
        let mut locus_snps = 0;
        let mut locus_indel_events = 0;
        let mut locus_indel_bases = 0;
        let mut locus_count = 0;

        for ((cache_locus, _, _, _), (snps, indel_events, indel_bases)) in &cache.data {
            if cache_locus == *locus {
                locus_snps += snps;
                locus_indel_events += indel_events;
                locus_indel_bases += indel_bases;
                locus_count += 1;
            }
        }

        let avg_snps = if locus_count > 0 {
            locus_snps as f64 / locus_count as f64
        } else {
            0.0
        };
        let avg_indel_events = if locus_count > 0 {
            locus_indel_events as f64 / locus_count as f64
        } else {
            0.0
        };
        let avg_indel_bases = if locus_count > 0 {
            locus_indel_bases as f64 / locus_count as f64
        } else {
            0.0
        };

        println!(
            "{:<25} {:>8} {:>10.2} {:>12.2} {:>14.2}",
            locus, entry_count, avg_snps, avg_indel_events, avg_indel_bases
        );
    }

    if !args.detailed && loci_info.len() > args.top_loci {
        println!(
            "... and {} more loci (use --detailed to show all)",
            loci_info.len() - args.top_loci
        );
    }
}

fn analyze_specific_locus(cache: &UltraFastCache, locus_name: &str) {
    println!("\n=== LOCUS DETAILS: {} ===", locus_name);

    let mut locus_entries = Vec::new();
    let mut unique_crcs = HashSet::new();

    for ((cache_locus, crc1, crc2, config_hash), (snps, indel_events, indel_bases)) in &cache.data {
        if cache_locus == locus_name {
            locus_entries.push((
                (*crc1, *crc2, *config_hash),
                (*snps, *indel_events, *indel_bases),
            ));
            unique_crcs.insert(*crc1);
            unique_crcs.insert(*crc2);
        }
    }

    if locus_entries.is_empty() {
        println!("❌ Locus '{}' not found in cache", locus_name);

        let mut available_loci: Vec<String> = cache
            .data
            .keys()
            .map(|(locus, _, _, _)| locus.clone())
            .collect::<HashSet<_>>()
            .into_iter()
            .collect();
        available_loci.sort();

        println!("Available loci: {:?}", available_loci);
        return;
    }

    println!("Total entries: {}", locus_entries.len());
    println!("Unique CRCs: {}", unique_crcs.len());

    // Statistics
    let total_snps: usize = locus_entries.iter().map(|(_, (snps, _, _))| snps).sum();
    let total_indel_events: usize = locus_entries
        .iter()
        .map(|(_, (_, indel_events, _))| indel_events)
        .sum();
    let total_indel_bases: usize = locus_entries
        .iter()
        .map(|(_, (_, _, indel_bases))| indel_bases)
        .sum();

    let avg_snps = total_snps as f64 / locus_entries.len() as f64;
    let avg_indel_events = total_indel_events as f64 / locus_entries.len() as f64;
    let avg_indel_bases = total_indel_bases as f64 / locus_entries.len() as f64;

    println!("Average SNPs per pair: {:.2}", avg_snps);
    println!("Average indel events per pair: {:.2}", avg_indel_events);
    println!("Average indel bases per pair: {:.2}", avg_indel_bases);

    // Show some example entries
    println!("\nExample cache entries (showing up to 10):");
    println!(
        "{:<12} {:<12} {:<16} {:>6} {:>8} {:>10}",
        "CRC1", "CRC2", "ConfigHash", "SNPs", "Indels", "IndelBases"
    );
    println!("{}", "-".repeat(75));

    for (i, ((crc1, crc2, config_hash), (snps, indel_events, indel_bases))) in
        locus_entries.iter().enumerate()
    {
        if i >= 10 {
            break;
        }
        println!(
            "{:<12} {:<12} {:<16} {:>6} {:>8} {:>10}",
            crc1, crc2, config_hash, snps, indel_events, indel_bases
        );
    }

    if locus_entries.len() > 10 {
        println!("... and {} more entries", locus_entries.len() - 10);
    }
}

fn check_compatibility(cache: &UltraFastCache, params_str: &str) {
    println!("\n=== COMPATIBILITY CHECK ===");

    // Parse parameters string (format: "match,mismatch,gap_open,gap_extend")
    let parts: Vec<&str> = params_str.split(',').collect();
    if parts.len() != 4 {
        println!("❌ ERROR: Invalid parameter format. Use: match,mismatch,gap_open,gap_extend");
        println!("   Example: 2,-1,5,2");
        return;
    }

    let target_config = match (
        parts[0].parse::<i32>(),
        parts[1].parse::<i32>(),
        parts[2].parse::<i32>(),
        parts[3].parse::<i32>(),
    ) {
        (Ok(match_score), Ok(mismatch_penalty), Ok(gap_open), Ok(gap_extend)) => AlignmentConfig {
            match_score,
            mismatch_penalty,
            gap_open,
            gap_extend,
            description: None,
        },
        _ => {
            println!("❌ ERROR: Failed to parse parameters. All values must be integers.");
            return;
        }
    };

    println!("Cache alignment config:");
    println!(
        "  Match: {}, Mismatch: {}, Gap open: {}, Gap extend: {}",
        cache.alignment_config.match_score,
        cache.alignment_config.mismatch_penalty,
        cache.alignment_config.gap_open,
        cache.alignment_config.gap_extend
    );
    println!("  Mode: {}", cache.alignment_config.detect_mode());

    println!("\nTarget alignment config:");
    println!(
        "  Match: {}, Mismatch: {}, Gap open: {}, Gap extend: {}",
        target_config.match_score,
        target_config.mismatch_penalty,
        target_config.gap_open,
        target_config.gap_extend
    );
    println!("  Mode: {}", target_config.detect_mode());

    if cache.alignment_config == target_config {
        println!("\n✅ COMPATIBLE: Cache and target configurations match perfectly!");
        println!("   You can safely reuse this cache with your target parameters.");
    } else {
        println!("\n❌ INCOMPATIBLE: Cache and target configurations differ!");
        println!("   Using this cache with different parameters will produce incorrect results.");
        println!("   Recommendation: Use a different cache file or let cgdist create a new cache.");

        // Show specific differences
        println!("\nDifferences:");
        if cache.alignment_config.match_score != target_config.match_score {
            println!(
                "  • Match score: {} → {}",
                cache.alignment_config.match_score, target_config.match_score
            );
        }
        if cache.alignment_config.mismatch_penalty != target_config.mismatch_penalty {
            println!(
                "  • Mismatch penalty: {} → {}",
                cache.alignment_config.mismatch_penalty, target_config.mismatch_penalty
            );
        }
        if cache.alignment_config.gap_open != target_config.gap_open {
            println!(
                "  • Gap open: {} → {}",
                cache.alignment_config.gap_open, target_config.gap_open
            );
        }
        if cache.alignment_config.gap_extend != target_config.gap_extend {
            println!(
                "  • Gap extend: {} → {}",
                cache.alignment_config.gap_extend, target_config.gap_extend
            );
        }
    }
}

fn validate_cache_integrity(cache: &UltraFastCache) -> bool {
    println!("\n=== CACHE VALIDATION ===");

    let mut errors = 0;
    let mut warnings = 0;

    // Check that all entries use the same configuration hash
    let mut config_hashes = HashSet::new();
    for (_, _, _, config_hash) in cache.data.keys() {
        config_hashes.insert(*config_hash);
    }

    if config_hashes.len() > 1 {
        println!(
            "❌ ERROR: Multiple configuration hashes found ({} different hashes)",
            config_hashes.len()
        );
        println!("   This indicates cache corruption or mixed alignment configurations");
        errors += 1;
    }

    // Check CRC pair ordering
    let mut unordered_pairs = 0;
    for (locus, crc1, crc2, _) in cache.data.keys() {
        if crc1 > crc2 {
            if unordered_pairs < 5 {
                // Only show first 5 examples
                println!(
                    "⚠️  WARNING: Unordered CRC pair in locus '{}': {} > {}",
                    locus, crc1, crc2
                );
            }
            unordered_pairs += 1;
        }
    }

    if unordered_pairs > 0 {
        if unordered_pairs > 5 {
            println!(
                "⚠️  WARNING: {} total unordered CRC pairs found",
                unordered_pairs
            );
        }
        warnings += 1;
    }

    // Summary
    if errors == 0 && warnings == 0 {
        println!("✅ Cache validation passed - no issues found");
        true
    } else {
        println!(
            "⚠️  Cache validation completed: {} errors, {} warnings",
            errors, warnings
        );
        if errors > 0 {
            println!("❌ Cache has integrity issues that should be addressed");
            false
        } else {
            println!("✅ Cache is valid but has minor warnings");
            true
        }
    }
}

fn export_summary_to_tsv(
    cache: &UltraFastCache,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(output_path)?;

    // Write header with cache metadata
    writeln!(file, "# cgdist cache summary")?;
    writeln!(file, "# Version: {}", cache.version)?;
    writeln!(file, "# Created: {}", cache.created)?;
    writeln!(
        file,
        "# Alignment mode: {}",
        cache.alignment_config.detect_mode()
    )?;
    writeln!(
        file,
        "# Parameters: match={}, mismatch={}, gap_open={}, gap_extend={}",
        cache.alignment_config.match_score,
        cache.alignment_config.mismatch_penalty,
        cache.alignment_config.gap_open,
        cache.alignment_config.gap_extend
    )?;
    writeln!(file)?;

    writeln!(
        file,
        "locus\tentries\tavg_snps\tavg_indel_events\tavg_indel_bases"
    )?;

    // Count entries per locus and calculate averages
    let mut locus_data = HashMap::new();

    for ((locus, _, _, _), (snps, indel_events, indel_bases)) in &cache.data {
        let entry = locus_data.entry(locus.clone()).or_insert((0, 0, 0, 0));
        entry.0 += 1; // count
        entry.1 += snps; // total snps
        entry.2 += indel_events; // total indel events
        entry.3 += indel_bases; // total indel bases
    }

    // Sort loci by name for consistent output
    let mut loci: Vec<_> = locus_data.keys().collect();
    loci.sort();

    for locus in loci {
        let (count, total_snps, total_indel_events, total_indel_bases) = locus_data[locus];
        let avg_snps = total_snps as f64 / count as f64;
        let avg_indel_events = total_indel_events as f64 / count as f64;
        let avg_indel_bases = total_indel_bases as f64 / count as f64;

        writeln!(
            file,
            "{}\t{}\t{:.3}\t{:.3}\t{:.3}",
            locus, count, avg_snps, avg_indel_events, avg_indel_bases
        )?;
    }

    println!("✅ Summary exported to: {}", output_path);
    Ok(())
}

// ============================================================================
// MAIN FUNCTION
// ============================================================================

fn main() {
    let args: Args = argh::from_env();

    if !args.quiet {
        println!("🔍 CGDist Cache Inspector (Config-Aware)");
        println!("==========================================");
    }

    // Load cache
    let cache_path = Path::new(&args.cache);
    let cache = match load_cache(cache_path, args.quiet) {
        Ok(cache) => cache,
        Err(e) => {
            eprintln!("❌ ERROR loading cache: {}", e);
            std::process::exit(1);
        }
    };

    // Main analysis
    analyze_cache_overview(&cache, &args);
    analyze_alignment_config(&cache, &args);
    analyze_loci_overview(&cache, &args);

    // Show specific locus if requested
    if let Some(locus_name) = &args.show_locus {
        analyze_specific_locus(&cache, locus_name);
    }

    // Check compatibility with given parameters
    if let Some(params) = &args.check_compatibility {
        check_compatibility(&cache, params);
    }

    // Validation
    if args.validate {
        let is_valid = validate_cache_integrity(&cache);
        if !is_valid {
            std::process::exit(1);
        }
    }

    // Export summary
    if let Some(export_path) = &args.export_summary {
        if let Err(e) = export_summary_to_tsv(&cache, export_path) {
            eprintln!("❌ ERROR exporting summary: {}", e);
            std::process::exit(1);
        }
    }

    if !args.quiet {
        println!("\n✅ Cache inspection completed successfully");
        println!("\nUsage examples:");
        println!("  --detect-mode                        Detect alignment mode from parameters");
        println!("  --check-compatibility 2,-1,5,2      Check compatibility with parameters");
        println!("  --detailed                           Show all loci details");
        println!("  --show-locus locus1                  Show specific locus details");
        println!("  --validate                           Validate cache integrity");
        println!("  --export-summary out.tsv             Export summary to TSV");
        println!("  --quiet                              Minimal output mode");

        println!("\nAlignment mode detection:");
        println!("  dna (2,-1,5,2)                       Balanced DNA parameters (BWA-MEM style)");
        println!("  dna-strict (3,-2,8,3)                Strict DNA parameters (HOXD style)");
        println!("  dna-permissive (2,-1,3,1)            Permissive DNA parameters");
        println!("  custom                               Any other parameter combination");
    }
}
