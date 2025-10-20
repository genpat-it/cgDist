// main.rs - CLI entry point

#![allow(unknown_lints, clippy::manual_is_multiple_of)]

use std::collections::HashSet;
use std::time::Instant;
// Removed unused import
use cgdist::cli::Config;
use cgdist::core::distance::ModernCache;
use cgdist::data::SequenceDatabase;
use cgdist::prelude::*;

// Type alias for complex cache key/value types
type LegacyCacheData = std::collections::HashMap<(String, u32, u32, u64), (usize, usize, usize)>;

#[derive(serde::Deserialize)]
struct LegacyCache {
    #[allow(dead_code)]
    data: LegacyCacheData,
    #[allow(dead_code)]
    version: Option<String>,
    #[allow(dead_code)]
    created: Option<String>,
    #[allow(dead_code)]
    alignment_config: AlignmentConfig,
}

fn main() {
    if let Err(e) = run_main() {
        eprintln!("❌ ERROR: {}", e);
        std::process::exit(1);
    }
}

fn run_main() -> Result<(), String> {
    let mut args: Args = argh::from_env();
    let command_line = std::env::args().collect::<Vec<String>>().join(" ");

    // Handle generate config first
    if args.generate_config {
        let sample_config = Config::generate_sample();
        println!("{}", sample_config);
        println!("\n💡 Save this content to a .toml file and use --config /path/to/config.toml");
        return Ok(());
    }

    // Load configuration file if specified
    if let Some(config_path) = args.config.clone() {
        args = args.with_config_file(&config_path)?;
    }

    // Handle inspector mode
    if let Some(cache_path) = &args.inspector {
        run_inspector(cache_path);
        return Ok(());
    }

    // Validate required parameters when not using inspector
    let profiles = args.profiles.as_ref().ok_or("--profiles is required")?;

    let output = if args.stats_only || args.benchmark || args.cache_only {
        None
    } else {
        Some(args.output.as_ref().ok_or("--output is required")?)
    };

    // Schema is only required for non-hamming hashers and when not using stats-only or cache-only
    let schema = if args.hasher_type == "hamming" || args.stats_only {
        None
    } else {
        Some(
            args.schema
                .as_ref()
                .ok_or("--schema is required (not needed for hamming hasher or stats-only mode)")?,
        )
    };

    println!("🚀 cgDist v{}", env!("CARGO_PKG_VERSION"));
    println!("⚡ Strategy: Pre-compute unique pairs → Batch Parasail → Fast assembly");

    // Configure thread pool
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to configure thread pool");
        println!("🧵 Threads: {}", n);
    } else {
        let num_threads = rayon::current_num_threads();
        println!("🧵 Threads: {} (auto-detected)", num_threads);
    }

    // Initialize hasher registry and validate hasher type
    let registry = HasherRegistry::new();
    println!("🔨 Hasher: {}", args.hasher_type);

    if !registry.has_hasher(&args.hasher_type) {
        eprintln!("❌ Error: Invalid hasher type '{}'", args.hasher_type);
        eprintln!("Available hashers:");
        for (name, desc) in registry.list_hashers() {
            eprintln!("  - {}: {}", name, desc);
        }
        std::process::exit(1);
    }

    // Handle benchmark mode
    if args.benchmark {
        return run_benchmark(&args);
    }

    // Validate all arguments
    let validation_result = validate_args(&args)?;

    let total_start = Instant::now();

    // Load allelic profiles
    let matrix = match AllelicMatrix::from_file_with_hasher(
        std::path::Path::new(profiles), // Now always required
        &args.missing_char,
        &args.hasher_type,
        args.sample_threshold,
        args.locus_threshold,
        validation_result.sample_include_regex.as_ref(),
        validation_result.sample_exclude_regex.as_ref(),
        validation_result.samples_include_set.as_ref(),
        validation_result.samples_exclude_set.as_ref(),
        validation_result.loci_include_regex.as_ref(),
        validation_result.loci_exclude_regex.as_ref(),
        validation_result.loci_include_set.as_ref(),
        validation_result.loci_exclude_set.as_ref(),
    ) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("❌ ERROR loading profiles: {}", e);
            std::process::exit(1);
        }
    };

    if args.dry_run {
        println!("✅ Dry run completed successfully");
        println!(
            "📊 Final matrix: {} samples × {} loci",
            matrix.samples.len(),
            matrix.loci_names.len()
        );
        return Ok(());
    }

    // Handle stats-only mode
    if args.stats_only {
        println!("\n📈 === GENETIC DIVERSITY ANALYSIS ===");
        let diversity_stats = matrix.calculate_diversity_metrics();

        println!("🧬 Allelic Diversity Metrics:");
        println!(
            "  • Average unique alleles per locus: {:.1}",
            diversity_stats.avg_unique_alleles
        );
        println!(
            "  • Diversity index: {:.3} ({})",
            diversity_stats.diversity_index, diversity_stats.diversity_category
        );
        println!(
            "  • Total unique CRC pairs: {}",
            diversity_stats.total_unique_pairs
        );

        if diversity_stats.diversity_index < 0.3 {
            println!("  🟢 Population appears clonal with low genetic diversity");
        } else if diversity_stats.diversity_index < 0.6 {
            println!("  🟡 Population shows moderate genetic diversity");
        } else {
            println!("  🔴 Population is highly diverse and genetically heterogeneous");
        }

        println!("\n✅ Statistics analysis completed");
        return Ok(());
    }

    if args.hasher_type == "hamming" {
        println!("\n🎯 Distance calculation: Hamming hasher (allelic level)");
    } else {
        println!(
            "\n🎯 Distance calculation mode: {} ({})",
            validation_result.distance_mode.description(),
            if args.no_hamming_fallback {
                "no Hamming fallback"
            } else {
                "Hamming fallback enabled"
            }
        );
    }

    if args.min_loci > 0 {
        println!("📏 Minimum shared loci: {}", args.min_loci);
    }

    // Pre-compute unique CRC pairs from the matrix
    println!("🔍 Pre-computing unique allele pairs from matrix...");
    let unique_pairs = extract_unique_crc_pairs(&matrix);
    if args.hasher_type == "hamming" {
        println!(
            "✅ Found {} unique CRC pairs for Hamming distance",
            unique_pairs.len()
        );
    } else {
        println!("✅ Found {} unique CRC pairs to align", unique_pairs.len());
    }

    // Load only required sequences from schema (skip for hamming hasher)
    let sequence_db = if args.hasher_type == "hamming" {
        SequenceDatabase::empty()
    } else {
        match SequenceDatabase::from_schema_selective(
            std::path::Path::new(schema.unwrap()),
            &unique_pairs,
            &matrix.loci_names,
        ) {
            Ok(db) => db,
            Err(e) => {
                eprintln!("❌ ERROR loading schema: {}", e);
                std::process::exit(1);
            }
        }
    };

    // Create distance engine with hasher type
    let mut engine = DistanceEngine::with_sequences(
        validation_result.alignment_config,
        sequence_db,
        args.hasher_type.clone(),
    );

    // Set cache note if provided
    if let Some(note) = &args.cache_note {
        engine.set_cache_note(note.clone());
    }

    // Set alignment saving if requested
    if let Some(ref save_path) = args.save_alignments {
        engine.set_save_alignments(save_path.clone());
        println!("💾 Alignment details will be saved to: {}", save_path);
    }

    // Try to load existing cache first (skip for hamming hasher)
    if args.hasher_type != "hamming" {
        if let Some(ref cache_path) = args.cache_file {
            if std::path::Path::new(cache_path).exists() && !args.force_recompute {
                // Compatibility already checked early, proceed with full cache load
                match engine.load_cache(cache_path, validation_result.distance_mode) {
                    Ok(()) => {
                        let (loaded_entries, _) = engine.cache_stats();
                        println!(
                            "🎯 Cache loaded successfully with {} entries",
                            loaded_entries
                        );

                        // Check for automatic or explicit enrichment with sequence lengths
                        let should_enrich = if args.enrich_lengths {
                            // User explicitly requested enrichment
                            true
                        } else if loaded_entries == 0 {
                            // Empty cache - enrich automatically if schema is available
                            schema.is_some()
                        } else {
                            // Non-empty cache - no automatic enrichment unless explicitly requested
                            false
                        };

                        if should_enrich {
                            if let Some(schema_path) = schema {
                                let enrich_reason = if args.enrich_lengths {
                                    "explicitly requested"
                                } else {
                                    "empty cache (automatic)"
                                };
                                println!(
                                    "🔍 Enriching cache with sequence lengths from schema ({})...",
                                    enrich_reason
                                );
                                let output_cache_path =
                                    args.enrich_output.as_ref().unwrap_or(cache_path);
                                match engine.enrich_cache_with_lengths_from_input(
                                    schema_path,
                                    cache_path,
                                    output_cache_path,
                                ) {
                                    Ok(enriched_count) => {
                                        if args.enrich_output.is_some() {
                                            println!("✅ Enriched {} entries with sequence lengths, saved to: {}", enriched_count, output_cache_path);
                                        } else {
                                            println!("✅ Enriched {} entries with sequence lengths (cache updated in-place)", enriched_count);
                                        }
                                    }
                                    Err(e) => {
                                        println!("⚠️  Failed to enrich cache: {}", e);
                                    }
                                }
                            } else {
                                println!("⚠️  Cannot enrich cache: schema path required");
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!("❌ FATAL ERROR loading cache: {}", e);
                        eprintln!("💡 Solutions:");
                        eprintln!("   - Use --force-recompute to ignore existing cache");
                        eprintln!("   - Delete incompatible cache file");
                        eprintln!("   - Use consistent alignment parameters");
                        std::process::exit(1);
                    }
                }
            }

            // Always run precompute to handle missing pairs and show stats
            println!("🚀 Pre-computing missing alignments...");
            engine.precompute_alignments(&unique_pairs, validation_result.distance_mode);
        } else {
            println!("🚀 Pre-computing all alignments in batch...");
            engine.precompute_alignments(&unique_pairs, validation_result.distance_mode);
        }
    } else {
        println!("🔨 Using Hamming hasher - no alignment precomputation needed");
    }

    // Check if cache-only mode
    if args.cache_only {
        println!("\n✅ Cache-only mode: Alignments computed and cached successfully");
        if let Some(ref cache_path) = args.cache_file {
            let had_new_entries = engine.has_new_entries();
            if had_new_entries {
                if let Err(e) = engine.save_cache(cache_path, validation_result.distance_mode) {
                    eprintln!("⚠️  Warning: Failed to save cache: {}", e);
                } else {
                    println!("📁 Cache saved to: {}", cache_path);

                    // Check if we should automatically enrich the newly saved cache
                    let should_auto_enrich = !args.enrich_lengths && schema.is_some();

                    if should_auto_enrich {
                        if let Some(schema_path) = schema {
                            println!("🔍 Enriching cache with sequence lengths (automatic)...");
                            let output_cache_path =
                                args.enrich_output.as_ref().unwrap_or(cache_path);
                            match engine.enrich_cache_with_lengths_from_input(
                                schema_path,
                                cache_path,
                                output_cache_path,
                            ) {
                                Ok(enriched_count) => {
                                    if args.enrich_output.is_some() {
                                        println!("✅ Enriched {} entries with sequence lengths, saved to: {}", enriched_count, output_cache_path);
                                    } else {
                                        println!("✅ Enriched {} entries with sequence lengths (cache updated in-place)", enriched_count);
                                    }
                                }
                                Err(e) => {
                                    println!("⚠️  Warning: Failed to enrich cache: {}", e);
                                }
                            }
                        }
                    }
                }
            } else {
                println!("📌 Cache unchanged - skipping save");
            }
        } else {
            println!("⚠️  Warning: No cache file specified - computed alignments not saved");
        }
        let total_elapsed = total_start.elapsed();
        println!(
            "\n⏱️  Total execution time: {:.2}s",
            total_elapsed.as_secs_f64()
        );
        return Ok(());
    }

    // Set recombination detection if requested (TODO: implement)
    if let Some(ref _log_path) = args.recombination_log {
        // TODO: Implement recombination detection
        println!(
            "🧬 Recombination detection requested (threshold: {} bp) - not implemented yet",
            args.recombination_threshold
        );
    }

    // Calculate distance matrix
    println!("\n🔄 Computing distance matrix...");
    let distance_matrix = calculate_distance_matrix(
        &matrix.samples,
        &matrix.loci_names,
        &engine,
        validation_result.distance_mode,
        args.min_loci,
        args.no_hamming_fallback,
    );

    // Write output
    let output_path = output.unwrap(); // Safe because stats_only returns early
    if let Err(e) = write_matrix(
        output_path,
        &args.format,
        &matrix.samples,
        &distance_matrix,
        &command_line,
    ) {
        eprintln!("❌ ERROR writing output: {}", e);
        std::process::exit(1);
    }

    // Save cache if specified (skip for hamming hasher) and only if there are new entries
    if args.hasher_type != "hamming" {
        if let Some(ref cache_path) = args.cache_file {
            let had_new_entries = engine.has_new_entries();
            if had_new_entries {
                if let Err(e) = engine.save_cache(cache_path, validation_result.distance_mode) {
                    eprintln!("⚠️  Warning: Failed to save cache: {}", e);
                } else {
                    // Check if we should automatically enrich the newly saved cache
                    // This covers the case where cache was created for the first time
                    let should_auto_enrich = !args.enrich_lengths && schema.is_some();

                    if should_auto_enrich {
                        if let Some(schema_path) = schema {
                            println!("🔍 Enriching newly created cache with sequence lengths (automatic)...");
                            let output_cache_path =
                                args.enrich_output.as_ref().unwrap_or(cache_path);
                            match engine.enrich_cache_with_lengths_from_input(
                                schema_path,
                                cache_path,
                                output_cache_path,
                            ) {
                                Ok(enriched_count) => {
                                    if args.enrich_output.is_some() {
                                        println!("✅ Enriched {} entries with sequence lengths, saved to: {}", enriched_count, output_cache_path);
                                    } else {
                                        println!("✅ Enriched {} entries with sequence lengths (cache updated in-place)", enriched_count);
                                    }
                                }
                                Err(e) => {
                                    println!(
                                        "⚠️  Warning: Failed to enrich newly created cache: {}",
                                        e
                                    );
                                }
                            }
                        }
                    }
                }
            } else {
                println!("📌 Cache unchanged - skipping save");
            }
        }
    }

    // Save alignment details if requested
    if let Err(e) = engine.save_alignments() {
        eprintln!("⚠️  Warning: Failed to save alignment details: {}", e);
    }

    // Save recombination log if specified (TODO: implement)
    if let Some(ref _log_path) = args.recombination_log {
        println!("⚠️  Recombination log saving not implemented yet");
    }

    // Print summary
    let total_elapsed = total_start.elapsed();
    println!("\n🎉 === CGDIST COMPLETED SUCCESSFULLY ===");
    println!(
        "⏱️  Total execution time: {:.2}s",
        total_elapsed.as_secs_f64()
    );
    println!(
        "📊 Final matrix: {} samples × {} loci",
        matrix.samples.len(),
        matrix.loci_names.len()
    );
    println!("📁 Output written to: {}", output_path);
    println!("🔧 Command: {}", command_line);

    Ok(())
}

/// Run benchmark mode to measure alignment processing speed with real data
fn run_benchmark(args: &Args) -> Result<(), String> {
    println!("\n🏁 === BENCHMARK MODE ===");
    println!(
        "📊 Testing alignment speed with {} threads",
        rayon::current_num_threads()
    );

    // Get initial memory
    let _initial_memory = get_memory_usage();
    let _start_total = Instant::now();

    // Validate arguments
    let validation_result = validate_args(args)?;

    let profiles = args.profiles.as_ref().unwrap();

    // Quick setup - minimal data loading
    println!("\n🚀 Quick setup...");
    let setup_start = Instant::now();

    // Load just enough data to get unique pairs
    let matrix = match AllelicMatrix::from_file_with_hasher(
        std::path::Path::new(profiles),
        &args.missing_char,
        &args.hasher_type,
        args.sample_threshold,
        args.locus_threshold,
        validation_result.sample_include_regex.as_ref(),
        validation_result.sample_exclude_regex.as_ref(),
        validation_result.samples_include_set.as_ref(),
        validation_result.samples_exclude_set.as_ref(),
        validation_result.loci_include_regex.as_ref(),
        validation_result.loci_exclude_regex.as_ref(),
        validation_result.loci_include_set.as_ref(),
        validation_result.loci_exclude_set.as_ref(),
    ) {
        Ok(m) => m,
        Err(e) => return Err(format!("Failed to load profiles: {}", e)),
    };

    let unique_pairs = extract_unique_crc_pairs(&matrix);
    println!(
        "   ✅ Setup complete: {} samples, {} loci, {} unique pairs",
        matrix.samples.len(),
        matrix.loci_names.len(),
        unique_pairs.len()
    );

    // Load all necessary sequences for realistic benchmark
    let sequence_db = if args.hasher_type != "hamming" {
        let schema = args
            .schema
            .as_ref()
            .ok_or("--schema is required for benchmark mode")?;

        match SequenceDatabase::from_schema_selective(
            std::path::Path::new(schema),
            &unique_pairs,
            &matrix.loci_names,
        ) {
            Ok(db) => {
                println!("   📂 Loaded {} sequences from schema", db.total_sequences);
                Some(db)
            }
            Err(e) => {
                println!("   ⚠️  Failed to load schema: {}", e);
                None
            }
        }
    } else {
        None
    };

    let _engine = if let Some(ref seq_db) = sequence_db {
        DistanceEngine::with_sequences(
            validation_result.alignment_config.clone(),
            seq_db.clone(),
            args.hasher_type.clone(),
        )
    } else {
        DistanceEngine::new(
            validation_result.alignment_config.clone(),
            args.hasher_type.clone(),
        )
    };

    let setup_elapsed = setup_start.elapsed();
    println!("   ⏱️  Setup time: {:.2}s", setup_elapsed.as_secs_f64());

    // Create a single engine for the entire benchmark
    let _engine = if args.hasher_type != "hamming" {
        if let Some(seq_db) = sequence_db.as_ref() {
            DistanceEngine::with_sequences(
                validation_result.alignment_config.clone(),
                seq_db.clone(),
                args.hasher_type.clone(),
            )
        } else {
            DistanceEngine::new(
                validation_result.alignment_config.clone(),
                args.hasher_type.clone(),
            )
        }
    } else {
        DistanceEngine::new(
            validation_result.alignment_config.clone(),
            args.hasher_type.clone(),
        )
    };

    // Benchmark with configurable timeout
    println!(
        "\n⚡ Benchmarking alignment speed for {} seconds...",
        args.benchmark_duration
    );
    println!("📊 Mode: {:?}", validation_result.distance_mode);
    println!("🔧 Hasher: {}", args.hasher_type);

    // Analyze sequence lengths
    let (total_sequences, min_len, max_len, avg_len) = if let Some(ref seq_db) = sequence_db {
        let mut lengths = Vec::new();
        for locus_map in seq_db.loci.values() {
            for seq_info in locus_map.values() {
                lengths.push(seq_info.sequence.len());
            }
        }

        if !lengths.is_empty() {
            let min = *lengths.iter().min().unwrap();
            let max = *lengths.iter().max().unwrap();
            let avg = lengths.iter().sum::<usize>() as f64 / lengths.len() as f64;
            (lengths.len(), min, max, avg)
        } else {
            (0, 0, 0, 0.0)
        }
    } else {
        (0, 0, 0, 0.0)
    };

    println!("🧬 Loaded sequences: {}", total_sequences);
    if total_sequences > 0 {
        println!(
            "📏 Sequence lengths: min={}bp, max={}bp, avg={:.0}bp",
            min_len, max_len, avg_len
        );
    }

    let benchmark_start = Instant::now();
    let timeout = std::time::Duration::from_secs(args.benchmark_duration);
    let mut max_pairs_per_second = 0.0;
    let mut total_pairs_processed = 0;

    // Memory tracking
    let initial_memory = get_memory_usage();
    let peak_memory = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(initial_memory));
    let peak_memory_clone = peak_memory.clone();

    // Spawn memory monitor thread
    let monitor_handle = std::thread::spawn(move || {
        while peak_memory_clone.load(std::sync::atomic::Ordering::Relaxed) != 0 {
            let current_mem = get_memory_usage();
            let mut peak = peak_memory_clone.load(std::sync::atomic::Ordering::Relaxed);
            while current_mem > peak
                && peak_memory_clone
                    .compare_exchange_weak(
                        peak,
                        current_mem,
                        std::sync::atomic::Ordering::Relaxed,
                        std::sync::atomic::Ordering::Relaxed,
                    )
                    .is_err()
            {
                peak = peak_memory_clone.load(std::sync::atomic::Ordering::Relaxed);
            }
            std::thread::sleep(std::time::Duration::from_millis(100));
        }
    });

    // Convert to vector for chunked processing
    let pairs_vec: Vec<_> = unique_pairs.iter().cloned().collect();
    let chunk_size = 1000;
    let mut chunk_idx = 0;

    println!(
        "⏱️  Processing chunks of {} pairs until timeout...",
        chunk_size
    );

    while benchmark_start.elapsed() < timeout && chunk_idx * chunk_size < pairs_vec.len() {
        let start_idx = chunk_idx * chunk_size;
        let end_idx = std::cmp::min(start_idx + chunk_size, pairs_vec.len());
        let chunk_pairs: HashSet<_> = pairs_vec[start_idx..end_idx].iter().cloned().collect();

        // Create fresh engine for this chunk to avoid cache
        let mut fresh_engine = if let Some(ref seq_db) = sequence_db {
            DistanceEngine::with_sequences(
                validation_result.alignment_config.clone(),
                seq_db.clone(),
                args.hasher_type.clone(),
            )
        } else {
            DistanceEngine::new(
                validation_result.alignment_config.clone(),
                args.hasher_type.clone(),
            )
        };

        let chunk_start = Instant::now();

        // Debug: count how many pairs actually have sequences
        let mut pairs_with_sequences = 0;
        let mut pairs_without_sequences = 0;

        if let Some(ref seq_db) = sequence_db {
            for (locus, crc1, crc2) in &chunk_pairs {
                let has_seq1 = seq_db.get_sequence(locus, *crc1).is_some();
                let has_seq2 = seq_db.get_sequence(locus, *crc2).is_some();
                if has_seq1 && has_seq2 {
                    pairs_with_sequences += 1;
                } else {
                    pairs_without_sequences += 1;
                }
            }
        }

        fresh_engine.precompute_alignments(&chunk_pairs, validation_result.distance_mode);
        let chunk_elapsed = chunk_start.elapsed();

        // Show debug info for first chunk
        if chunk_idx == 0 {
            println!("   🔍 Debug chunk analysis:");
            println!("     • Pairs with both sequences: {}", pairs_with_sequences);
            println!(
                "     • Pairs missing sequences: {}",
                pairs_without_sequences
            );
            println!(
                "     • Alignment rate: {:.0} pairs/sec",
                chunk_pairs.len() as f64 / chunk_elapsed.as_secs_f64()
            );
        }

        let chunk_pairs_per_second = chunk_pairs.len() as f64 / chunk_elapsed.as_secs_f64();
        if chunk_pairs_per_second > max_pairs_per_second {
            max_pairs_per_second = chunk_pairs_per_second;
        }

        total_pairs_processed += chunk_pairs.len();
        chunk_idx += 1;

        // Progress update every 5 seconds
        if benchmark_start.elapsed().as_secs() % 5 == 0 && chunk_idx % 3 == 0 {
            let elapsed = benchmark_start.elapsed().as_secs();
            let remaining = args.benchmark_duration - elapsed;
            println!(
                "   ⏱️  {}s elapsed, {}s remaining, processed {} pairs, peak: {:.0} pairs/sec",
                elapsed, remaining, total_pairs_processed, max_pairs_per_second
            );
        }
    }

    // Stop memory monitor
    let final_peak_memory = peak_memory.load(std::sync::atomic::Ordering::Relaxed);
    peak_memory.store(0, std::sync::atomic::Ordering::Relaxed);
    let _ = monitor_handle.join();

    let benchmark_elapsed = benchmark_start.elapsed();
    let total_elapsed = benchmark_elapsed.as_secs_f64();

    // Calculate ETA for full dataset
    let total_pairs = unique_pairs.len();
    let eta_seconds = total_pairs as f64 / max_pairs_per_second;
    let eta_minutes = eta_seconds / 60.0;
    let eta_hours = eta_minutes / 60.0;

    // Calculate average performance
    let avg_pairs_per_second = total_pairs_processed as f64 / total_elapsed;

    // Calculate completion percentage
    let completion_percentage = (total_pairs_processed as f64 / total_pairs as f64) * 100.0;

    // Results
    println!("\n📈 === BENCHMARK RESULTS ===");
    println!("⏱️  Test duration: {:.1} seconds", total_elapsed);
    println!(
        "🔢 Pairs processed: {} / {} ({:.2}%)",
        total_pairs_processed, total_pairs, completion_percentage
    );
    println!(
        "⚡ Peak performance: {:.0} pairs/second",
        max_pairs_per_second
    );
    println!(
        "📊 Average performance: {:.0} pairs/second",
        avg_pairs_per_second
    );
    println!("🧵 Using {} threads", rayon::current_num_threads());
    println!(
        "📊 Per-thread peak: {:.0} pairs/sec/thread",
        max_pairs_per_second / rayon::current_num_threads() as f64
    );
    println!(
        "📊 Per-thread average: {:.0} pairs/sec/thread",
        avg_pairs_per_second / rayon::current_num_threads() as f64
    );

    println!("\n💾 Memory usage:");
    println!("   • Initial: {} MB", initial_memory / 1024 / 1024);
    println!(
        "   • Peak: {} MB (+{} MB)",
        final_peak_memory / 1024 / 1024,
        (final_peak_memory - initial_memory) / 1024 / 1024
    );

    // Estimate memory for full dataset using diversity-based approach
    let memory_per_pair = if total_pairs_processed > 0 {
        (final_peak_memory - initial_memory) as f64 / total_pairs_processed as f64
    } else {
        0.0
    };

    // Get diversity metrics for smarter estimation
    let diversity_stats = matrix.calculate_diversity_metrics();

    // Diversity-based memory estimation (empirical formula based on real-world data)
    // Base memory: setup + sequences
    // Cache memory: depends on diversity index and number of loci
    let base_memory_mb = initial_memory as f64 / 1024.0 / 1024.0;
    let cache_memory_mb = if diversity_stats.diversity_index < 0.05 {
        // Very low diversity (clonal) - small cache
        1000.0 + (matrix.loci_names.len() as f64 * 0.3)
    } else if diversity_stats.diversity_index < 0.15 {
        // Low-moderate diversity
        2000.0 + (matrix.loci_names.len() as f64 * 0.5)
    } else if diversity_stats.diversity_index < 0.3 {
        // Moderate diversity
        3000.0 + (matrix.loci_names.len() as f64 * 0.8)
    } else {
        // High diversity - larger cache needed
        4000.0 + (matrix.loci_names.len() as f64 * 1.2)
    };

    let diversity_based_estimate_mb = base_memory_mb + cache_memory_mb;
    let linear_estimate_mb =
        (initial_memory as f64 + (memory_per_pair * total_pairs as f64)) / 1024.0 / 1024.0;

    println!(
        "   • Memory per pair: {:.1} KB/pair",
        memory_per_pair / 1024.0
    );
    println!(
        "   • Linear estimate: {:.0} MB ({:.1} GB) ⚠️  Likely overestimated",
        linear_estimate_mb,
        linear_estimate_mb / 1024.0
    );
    println!(
        "   • Diversity-based estimate: {:.0} MB ({:.1} GB) 🎯 More realistic",
        diversity_based_estimate_mb,
        diversity_based_estimate_mb / 1024.0
    );
    println!(
        "     └─ Based on diversity index {:.3} ({})",
        diversity_stats.diversity_index, diversity_stats.diversity_category
    );
    println!(
        "   • Per sample estimate: {:.1} MB/sample",
        (final_peak_memory - initial_memory) as f64 / 1024.0 / 1024.0 / matrix.samples.len() as f64
    );
    println!("   📊 Diversity-based estimate accounts for cache stabilization patterns");

    println!("\n🎯 === FULL DATASET ETA ===");
    println!("🔢 Total unique pairs: {}", total_pairs);

    if eta_hours >= 1.0 {
        println!("⏱️  ETA at peak speed: {:.1} hours", eta_hours);
    } else if eta_minutes >= 1.0 {
        println!("⏱️  ETA at peak speed: {:.1} minutes", eta_minutes);
    } else {
        println!("⏱️  ETA at peak speed: {:.0} seconds", eta_seconds);
    }

    // Save benchmark results to JSON
    let timestamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_secs();
    let json_filename = format!("cgdist_benchmark_{}.json", timestamp);

    // Reconstruct the command line
    let command_line = std::env::args().collect::<Vec<String>>().join(" ");

    let benchmark_results = serde_json::json!({
        "timestamp": timestamp,
        "command_line": command_line,
        "test_info": {
            "duration_seconds": total_elapsed,
            "timeout_seconds": args.benchmark_duration,
            "threads": rayon::current_num_threads(),
            "hasher_type": args.hasher_type,
            "distance_mode": format!("{:?}", validation_result.distance_mode)
        },
        "dataset_info": {
            "samples": matrix.samples.len(),
            "loci": matrix.loci_names.len(),
            "total_unique_pairs": total_pairs,
            "pairs_processed": total_pairs_processed,
            "completion_percentage": completion_percentage
        },
        "performance": {
            "peak_pairs_per_second": max_pairs_per_second,
            "average_pairs_per_second": avg_pairs_per_second,
            "per_thread_peak": max_pairs_per_second / rayon::current_num_threads() as f64,
            "per_thread_average": avg_pairs_per_second / rayon::current_num_threads() as f64
        },
        "memory_usage": {
            "initial_mb": initial_memory / 1024 / 1024,
            "peak_mb": final_peak_memory / 1024 / 1024,
            "peak_increase_mb": (final_peak_memory - initial_memory) / 1024 / 1024,
            "memory_per_pair_kb": memory_per_pair / 1024.0,
            "linear_estimate_mb": linear_estimate_mb,
            "linear_estimate_gb": linear_estimate_mb / 1024.0,
            "diversity_estimate_mb": diversity_based_estimate_mb,
            "diversity_estimate_gb": diversity_based_estimate_mb / 1024.0,
            "diversity_index": diversity_stats.diversity_index,
            "diversity_category": diversity_stats.diversity_category,
            "per_sample_mb": (final_peak_memory - initial_memory) as f64 / 1024.0 / 1024.0 / matrix.samples.len() as f64
        },
        "eta": {
            "eta_seconds": eta_seconds,
            "eta_minutes": eta_minutes,
            "eta_hours": eta_hours,
            "eta_formatted": if eta_hours >= 1.0 {
                format!("{:.1} hours", eta_hours)
            } else if eta_minutes >= 1.0 {
                format!("{:.1} minutes", eta_minutes)
            } else {
                format!("{:.0} seconds", eta_seconds)
            }
        },
        "sequence_info": {
            "total_sequences": total_sequences,
            "min_length_bp": min_len,
            "max_length_bp": max_len,
            "avg_length_bp": avg_len
        },
        "notes": {
            "linear_estimate_warning": format!("Linear estimate extrapolated from {:.1}s test - likely overestimated", total_elapsed),
            "diversity_estimate_note": "Diversity-based estimate uses genetic diversity patterns for more realistic prediction",
            "recommendation": "Use diversity-based estimate for memory planning"
        }
    });

    if let Err(e) = std::fs::write(
        &json_filename,
        serde_json::to_string_pretty(&benchmark_results).unwrap(),
    ) {
        println!(
            "⚠️  Failed to save benchmark results to {}: {}",
            json_filename, e
        );
    } else {
        println!("💾 Benchmark results saved to: {}", json_filename);
    }

    Ok(())
}

/// Get current memory usage in bytes
fn get_memory_usage() -> u64 {
    // Read from /proc/self/status on Linux
    #[cfg(target_os = "linux")]
    {
        use std::fs::read_to_string;

        if let Ok(status) = read_to_string("/proc/self/status") {
            for line in status.lines() {
                if line.starts_with("VmRSS:") {
                    let parts: Vec<&str> = line.split_whitespace().collect();
                    if parts.len() >= 2 {
                        if let Ok(kb) = parts[1].parse::<u64>() {
                            return kb * 1024; // Convert KB to bytes
                        }
                    }
                }
            }
        }
    }

    // Fallback for non-Linux or if reading fails
    0
}

/// Extract unique CRC pairs from allelic matrix
fn extract_unique_crc_pairs(matrix: &AllelicMatrix) -> HashSet<(String, u32, u32)> {
    let mut unique_pairs = HashSet::new();

    for locus in &matrix.loci_names {
        let mut locus_crcs = HashSet::new();

        // Collect all CRCs for this locus
        for sample in &matrix.samples {
            if let Some(hash) = sample.loci_hashes.get(locus) {
                if let Some(crc) = hash.as_crc32() {
                    if crc != u32::MAX {
                        // Skip missing values
                        locus_crcs.insert(crc);
                    }
                }
            }
        }

        // Generate all unique pairs for this locus
        let crcs: Vec<u32> = locus_crcs.into_iter().collect();
        for i in 0..crcs.len() {
            for j in i + 1..crcs.len() {
                let crc1 = crcs[i];
                let crc2 = crcs[j];
                unique_pairs.insert((locus.clone(), crc1.min(crc2), crc1.max(crc2)));
            }
        }
    }

    unique_pairs
}

/// Run integrated inspector functionality
fn run_inspector(cache_path: &str) {
    println!("🔍 cgDist Cache Inspector");
    println!("==========================");

    // Read and decompress cache file
    let compressed = match std::fs::read(cache_path) {
        Ok(data) => data,
        Err(e) => {
            eprintln!("❌ ERROR reading cache file: {}", e);
            std::process::exit(1);
        }
    };

    let decompressed = match lz4_flex::decompress_size_prepended(&compressed) {
        Ok(data) => data,
        Err(e) => {
            eprintln!("❌ ERROR decompressing cache: {}", e);
            std::process::exit(1);
        }
    };

    // Try to load as modern format first
    if let Ok(modern_cache) = serde_json::from_slice::<ModernCache>(&decompressed) {
        print_modern_cache_info(&modern_cache, &compressed);
    } else {
        // Try legacy format
        println!("🔄 Attempting to load as legacy format...");

        match bincode::deserialize::<LegacyCache>(&decompressed) {
            Ok(legacy_cache) => {
                print_legacy_cache_info(&legacy_cache, &compressed);
            }
            Err(e) => {
                eprintln!("❌ ERROR: Unable to parse cache file (tried both modern and legacy formats): {}", e);
                std::process::exit(1);
            }
        }
    }
}

fn print_modern_cache_info(cache: &ModernCache, compressed: &[u8]) {
    let metadata = &cache.metadata;

    println!("✅ Modern cache format (v{})", metadata.format_version);
    println!();

    println!("=== METADATA ===");
    println!("Version: {}", metadata.version);
    println!("Created: {}", metadata.created);
    println!("Last modified: {}", metadata.last_modified);
    println!("Hasher type: {}", metadata.hasher_type);
    println!("Distance mode: {}", metadata.distance_mode);
    println!("Total entries: {}", metadata.total_entries);
    println!("Unique loci: {}", metadata.unique_loci);

    if let Some(note) = &metadata.user_note {
        println!("User note: {}", note);
    }

    println!();
    println!("=== ALIGNMENT CONFIG ===");
    println!("Match score: {}", metadata.alignment_config.match_score);
    println!(
        "Mismatch penalty: {}",
        metadata.alignment_config.mismatch_penalty
    );
    println!("Gap open: {}", metadata.alignment_config.gap_open);
    println!("Gap extend: {}", metadata.alignment_config.gap_extend);

    println!();
    println!("=== STORAGE INFO ===");
    println!("Compressed size: {} KB", compressed.len() / 1024);
    println!(
        "Compression ratio: {:.2}:1",
        cache.data.len() as f64 * 50.0 / compressed.len() as f64
    );

    // Show top loci by entry count
    let mut locus_counts = std::collections::HashMap::new();
    for string_key in cache.data.keys() {
        // Parse locus from string key "locus:hash1:hash2"
        if let Some(locus) = string_key.split(':').next() {
            *locus_counts.entry(locus.to_string()).or_insert(0) += 1;
        }
    }

    let mut loci_sorted: Vec<_> = locus_counts.iter().collect();
    loci_sorted.sort_by(|a, b| b.1.cmp(a.1));

    println!();
    println!("=== TOP LOCI BY ENTRY COUNT ===");
    for (locus, count) in loci_sorted.iter().take(10) {
        println!("{:<20} {:>8} entries", locus, count);
    }

    if loci_sorted.len() > 10 {
        println!("... and {} more loci", loci_sorted.len() - 10);
    }

    // Show sample cache values to verify content
    println!();
    println!("=== SAMPLE CACHE VALUES (first 5 entries) ===");
    for (count, (string_key, cache_value)) in cache.data.iter().enumerate() {
        if count >= 5 {
            break;
        }
        println!("Key: {}", string_key);
        println!(
            "  └─ SNPs: {}, indel_events: {}, indel_bases: {}, computed: {}",
            cache_value.snps,
            cache_value.indel_events,
            cache_value.indel_bases,
            cache_value.computed_at
        );
    }

    // Show statistics of the values
    println!();
    println!("=== CACHE VALUES ANALYSIS ===");
    let mut snps_sum = 0;
    let mut indel_events_sum = 0;
    let mut indel_bases_sum = 0;
    let mut non_zero_snps = 0;
    let mut non_zero_indel_events = 0;
    let mut non_zero_indel_bases = 0;

    for cache_value in cache.data.values() {
        snps_sum += cache_value.snps;
        indel_events_sum += cache_value.indel_events;
        indel_bases_sum += cache_value.indel_bases;

        if cache_value.snps > 0 {
            non_zero_snps += 1;
        }
        if cache_value.indel_events > 0 {
            non_zero_indel_events += 1;
        }
        if cache_value.indel_bases > 0 {
            non_zero_indel_bases += 1;
        }
    }

    let total_entries = cache.data.len();
    println!("📊 Values distribution:");
    println!(
        "  SNPs: sum={}, non-zero entries={}/{} ({:.1}%)",
        snps_sum,
        non_zero_snps,
        total_entries,
        (non_zero_snps as f64 / total_entries as f64) * 100.0
    );
    println!(
        "  Indel events: sum={}, non-zero entries={}/{} ({:.1}%)",
        indel_events_sum,
        non_zero_indel_events,
        total_entries,
        (non_zero_indel_events as f64 / total_entries as f64) * 100.0
    );
    println!(
        "  Indel bases: sum={}, non-zero entries={}/{} ({:.1}%)",
        indel_bases_sum,
        non_zero_indel_bases,
        total_entries,
        (non_zero_indel_bases as f64 / total_entries as f64) * 100.0
    );

    // Determine if this is a legacy cache with all zeros for indels
    if indel_events_sum == 0 && indel_bases_sum == 0 {
        println!();
        println!("⚠️  WARNING: This cache contains ONLY SNPs data!");
        println!("   All indel_events and indel_bases values are 0.");
        println!("   This cache was generated with old code or SNPs-only mode.");
        println!("   💡 To use with indel modes, regenerate cache with new code.");
    } else {
        println!();
        println!("✅ GOOD: This cache contains complete alignment data!");
        println!(
            "   Can be used with any distance mode (SNPs, SNPs+indel_events, SNPs+indel_bases)."
        );
    }
}

fn print_legacy_cache_info(_cache: &LegacyCache, compressed: &[u8]) {
    println!("⚠️  Legacy cache format detected");
    println!();

    // We can't easily introspect the legacy format without defining all structures
    // Just show basic info
    println!("=== STORAGE INFO ===");
    println!("Compressed size: {} KB", compressed.len() / 1024);
    println!();
    println!("💡 Recommendation: Regenerate cache with modern format for better features");
    println!("   - Hasher type support");
    println!("   - User notes");
    println!("   - Enhanced metadata");
    println!("   - Better inspection capabilities");
}
