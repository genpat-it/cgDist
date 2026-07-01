// report.rs - Self-contained HTML analyst dashboard
//
// Emits a single .html file (no external assets, works offline) summarising a
// cgDist run for a surveillance analyst: dataset summary, per-sample data
// quality, pairwise-distance distribution, a clusters-vs-threshold curve, and an
// INTERACTIVE single-linkage threshold explorer (the outbreak-clustering tool).
// Opt-in via --report; does not affect any other output.

use crate::core::PairRow;
use crate::data::AllelicProfile;
use serde_json::{json, Value};
use std::fs::{create_dir_all, File};
use std::io::{BufWriter, Write};
use std::path::Path;

const EMBED_BUDGET: usize = 150_000; // max edges embedded for the interactive explorer

fn ensure_parent_dir(file_path: &str) -> Result<(), String> {
    if let Some(parent) = Path::new(file_path).parent() {
        create_dir_all(parent)
            .map_err(|e| format!("Failed to create parent directory '{}': {e}", parent.display()))?;
    }
    Ok(())
}

/// Simple union-find for single-linkage clustering.
struct UnionFind {
    parent: Vec<usize>,
    size: Vec<usize>,
    n_components: usize,
}
impl UnionFind {
    fn new(n: usize) -> Self {
        Self { parent: (0..n).collect(), size: vec![1; n], n_components: n }
    }
    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        let (big, small) = if self.size[ra] >= self.size[rb] { (ra, rb) } else { (rb, ra) };
        self.parent[small] = big;
        self.size[big] += self.size[small];
        self.n_components -= 1;
    }
    fn largest(&self) -> usize {
        self.size.iter().copied().max().unwrap_or(0)
    }
}

/// Build the JSON data blob and render the dashboard.
#[allow(clippy::too_many_arguments)]
pub fn write_html_report(
    file_path: &str,
    samples: &[AllelicProfile],
    loci_names: &[String],
    pair_rows: &[PairRow],
    recomb: Option<&[u32]>,
    recomb_threshold: f64,
    mode: &str,
    hasher: &str,
    command_line: &str,
) -> Result<(), String> {
    ensure_parent_dir(file_path)?;
    let n = samples.len();
    let n_loci = loci_names.len();
    let denom = n_loci.max(1) as f64;

    // --- per-sample completeness ---
    let completeness: Vec<f64> = samples
        .iter()
        .map(|s| {
            let present = loci_names
                .iter()
                .filter(|l| s.loci_hashes.get(*l).map(|h| !h.is_missing()).unwrap_or(false))
                .count();
            present as f64 / denom
        })
        .collect();

    // --- per-sample nearest neighbour + mean shared loci (single pass) ---
    let mut nn_dist = vec![usize::MAX; n];
    let mut nn_idx = vec![usize::MAX; n];
    let mut shared_sum = vec![0usize; n];
    let mut shared_cnt = vec![0usize; n];
    let mut dists: Vec<usize> = Vec::with_capacity(pair_rows.len());
    let mut na_pairs = 0usize;
    let mut edges: Vec<(usize, usize, usize)> = Vec::new();
    for &(i, j, dist, shared) in pair_rows {
        shared_sum[i] += shared;
        shared_sum[j] += shared;
        shared_cnt[i] += 1;
        shared_cnt[j] += 1;
        match dist {
            Some(d) => {
                dists.push(d);
                edges.push((i, j, d));
                if d < nn_dist[i] {
                    nn_dist[i] = d;
                    nn_idx[i] = j;
                }
                if d < nn_dist[j] {
                    nn_dist[j] = d;
                    nn_idx[j] = i;
                }
            }
            None => na_pairs += 1,
        }
    }

    // --- pairwise distance summary + histogram ---
    dists.sort_unstable();
    let (dmin, dmax, dmed, dmean) = if dists.is_empty() {
        (0, 0, 0.0, 0.0)
    } else {
        let sum: usize = dists.iter().sum();
        (
            dists[0],
            dists[dists.len() - 1],
            dists[dists.len() / 2] as f64,
            sum as f64 / dists.len() as f64,
        )
    };
    let n_bins = 40usize;
    let mut counts = vec![0u64; n_bins];
    let span = (dmax as f64).max(1.0);
    for &d in &dists {
        let b = (((d as f64) / span) * (n_bins as f64 - 1.0)).round() as usize;
        counts[b.min(n_bins - 1)] += 1;
    }
    let bin_edges: Vec<f64> =
        (0..=n_bins).map(|b| (b as f64) * span / (n_bins as f64)).collect();

    // --- clusters-vs-threshold curve (incremental single-linkage over sorted edges) ---
    // Choose a threshold ceiling covering the epidemiologically interesting range.
    let t_ceiling = if dmax == 0 {
        1
    } else {
        // cap the curve/explorer to a sensible outbreak range
        (dmed.ceil() as usize).clamp(20, 300).min(dmax)
    };
    let mut sorted_edges = edges.clone();
    sorted_edges.sort_unstable_by_key(|e| e.2);
    let mut uf = UnionFind::new(n);
    let mut curve: Vec<Value> = Vec::new();
    let mut ei = 0usize;
    for t in 0..=t_ceiling {
        while ei < sorted_edges.len() && sorted_edges[ei].2 <= t {
            uf.union(sorted_edges[ei].0, sorted_edges[ei].1);
            ei += 1;
        }
        let mut singletons = 0usize;
        for x in 0..n {
            let r = uf.find(x);
            if uf.size[r] == 1 {
                singletons += 1;
            }
        }
        curve.push(json!({
            "t": t, "n_clusters": uf.n_components,
            "largest": uf.largest(), "singletons": singletons
        }));
    }

    // --- edges embedded for the interactive explorer (bounded) ---
    let mut embed_cap = t_ceiling;
    let count_le = |cap: usize| edges.iter().filter(|e| e.2 <= cap).count();
    while embed_cap > 0 && count_le(embed_cap) > EMBED_BUDGET {
        embed_cap = embed_cap.saturating_sub(1);
    }
    let embed_edges: Vec<Value> = edges
        .iter()
        .filter(|e| e.2 <= embed_cap)
        .map(|&(i, j, d)| json!([i, j, d]))
        .collect();
    let truncated = embed_cap < t_ceiling;

    // --- recombination (optional; only when enriched-cache length data present) ---
    let mut recomb_sample = vec![0f64; n];       // mean recombinant loci per sample's pairs
    let mut recomb_json = Value::Null;
    if let Some(rc) = recomb {
        let mut sum = vec![0u64; n];
        let mut per_pair: Vec<usize> = Vec::with_capacity(rc.len());
        let mut top: Vec<(usize, usize, u32)> = Vec::new();
        let mut pairs_with_recomb = 0usize;
        let mut max_rc = 0u32;
        for (&(i, j, _, _), &r) in pair_rows.iter().zip(rc.iter()) {
            sum[i] += r as u64;
            sum[j] += r as u64;
            per_pair.push(r as usize);
            if r > 0 {
                pairs_with_recomb += 1;
                top.push((i, j, r));
            }
            max_rc = max_rc.max(r);
        }
        for i in 0..n {
            recomb_sample[i] = if shared_cnt[i] > 0 {
                sum[i] as f64 / shared_cnt[i] as f64
            } else {
                0.0
            };
        }
        top.sort_by(|a, b| b.2.cmp(&a.2));
        let top_json: Vec<Value> = top
            .iter()
            .take(30)
            .map(|&(i, j, r)| {
                json!({"i": samples[i].sample_id, "j": samples[j].sample_id, "recomb": r})
            })
            .collect();
        // histogram of recombinant-loci-per-pair
        let rmax = (max_rc as usize).max(1);
        let rbins = 30usize.min(rmax + 1);
        let mut rcounts = vec![0u64; rbins];
        for &v in &per_pair {
            let b = ((v as f64 / rmax as f64) * (rbins as f64 - 1.0)).round() as usize;
            rcounts[b.min(rbins - 1)] += 1;
        }
        let mean_rc = per_pair.iter().sum::<usize>() as f64 / per_pair.len().max(1) as f64;
        recomb_json = json!({
            "threshold_pct": (recomb_threshold * 100.0 * 100.0).round() / 100.0,
            "mean_per_pair": (mean_rc * 100.0).round() / 100.0,
            "max_per_pair": max_rc,
            "pairs_with_recomb": pairs_with_recomb,
            "frac_pairs_with_recomb": ((pairs_with_recomb as f64 / per_pair.len().max(1) as f64) * 1e4).round() / 1e4,
            "hist_counts": rcounts,
            "hist_max": rmax,
            "top_pairs": top_json,
        });
    }

    // --- per-sample records ---
    let has_recomb = recomb.is_some();
    let sample_json: Vec<Value> = (0..n)
        .map(|i| {
            json!({
                "id": samples[i].sample_id,
                "completeness": (completeness[i] * 1e4).round() / 1e4,
                "nn_dist": if nn_dist[i] == usize::MAX { Value::Null } else { json!(nn_dist[i]) },
                "nn_id": if nn_idx[i] == usize::MAX { Value::Null }
                         else { json!(samples[nn_idx[i]].sample_id) },
                "mean_shared": if shared_cnt[i] > 0 {
                    (shared_sum[i] as f64 / shared_cnt[i] as f64).round()
                } else { 0.0 },
                "recomb_load": if has_recomb { json!((recomb_sample[i] * 100.0).round() / 100.0) } else { Value::Null },
            })
        })
        .collect();

    let mean_comp = completeness.iter().sum::<f64>() / (n.max(1) as f64);
    let mut comp_sorted = completeness.clone();
    comp_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med_comp = comp_sorted.get(n / 2).copied().unwrap_or(0.0);

    let data = json!({
        "meta": {
            "generated": chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC").to_string(),
            "version": env!("CARGO_PKG_VERSION"),
            "command": command_line,
            "mode": mode, "hasher": hasher,
            "n_samples": n, "n_loci": n_loci,
        },
        "summary": {
            "mean_completeness": (mean_comp * 1e4).round() / 1e4,
            "median_completeness": (med_comp * 1e4).round() / 1e4,
            "overall_missing_frac": ((1.0 - mean_comp) * 1e4).round() / 1e4,
            "n_pairs": pair_rows.len(),
            "na_pairs": na_pairs,
            "dist_min": dmin, "dist_max": dmax,
            "dist_median": dmed, "dist_mean": (dmean * 100.0).round() / 100.0,
        },
        "hist": { "bin_edges": bin_edges, "counts": counts },
        "cluster_curve": curve,
        "recombination": recomb_json,
        "samples": sample_json,
        "edges": embed_edges,
        "embed_cap": embed_cap,
        "embed_truncated": truncated,
        "t_ceiling": t_ceiling,
        "n": n,
    });

    let html = TEMPLATE.replace("/*__CGDIST_DATA__*/", &data.to_string());
    let file = File::create(file_path)
        .map_err(|e| format!("Failed to create report '{file_path}': {e}"))?;
    let mut w = BufWriter::new(file);
    w.write_all(html.as_bytes()).map_err(|e| format!("Write error: {e}"))?;
    w.flush().map_err(|e| format!("Flush error: {e}"))?;
    println!("✅ Analyst dashboard written to: {file_path}");
    if truncated {
        println!(
            "   ℹ️  interactive clustering explorable up to distance {embed_cap} \
             (edge budget); the clusters-vs-threshold curve covers up to {t_ceiling}."
        );
    }
    Ok(())
}

const TEMPLATE: &str = include_str!("report_template.html");
