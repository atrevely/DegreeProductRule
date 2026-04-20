use rayon::prelude::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

use crate::growth::{dpr_growth_process, achlioptas_growth_process};
use crate::percolation::{cluster_sizes_range, cluster_sizes_range_with_hist, max_diff};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProcessType {
    Dpr,
    Achlioptas,
}

/// Results from PercMaxJumpMeans.
pub struct MaxJumpResult {
    /// [forward_mean, backward_mean] maximum jump in order parameter.
    pub crit_jump: [f64; 2],
    /// [forward_std, backward_std].
    pub jump_std: [f64; 2],
    /// Per-run jump locations: (forward_loc, backward_loc), both offset back
    /// to the full-sequence timestep (i.e. + s_point_{f,b}).
    pub loc_max_jump: Vec<[usize; 2]>,
    /// Largest cluster at the mean-field critical point for each run.
    pub max_clust_mean_fld: Vec<f64>,
    /// Largest cluster at the self-identified critical point for each run.
    pub max_clust_self_crit: Vec<f64>,
    /// Fraction of runs that reached full connectivity.
    pub frac_fully_conn: f64,
}

/// Mirrors PercMaxJumpMeans.m.
///
/// `mean_crit`: critical edge count expressed as a *fraction of num_nodes*
/// (e.g. 0.5 means N/2 edges).  The MATLAB parameter `meancrit` carried the
/// same semantics.
#[allow(clippy::too_many_arguments)]
pub fn perc_max_jump_means(
    num_nodes: usize,
    num_choices: usize,
    s_point_f: usize,
    e_point_f: usize,
    s_point_b: usize,
    e_point_b: usize,
    num_runs: usize,
    len: usize,
    mean_crit: f64,
    process: ProcessType,
    alpha: f64,
) -> MaxJumpResult {
    // Generate seeds outside the parallel loop (mirrors the MATLAB passedSeed trick).
    let master_seed: u64 = rand::random();
    let seeds: Vec<u64> = (0..num_runs as u64)
        .map(|i| master_seed.wrapping_add(i * 6364136223846793005))
        .collect();

    let mean_crit_steps = (mean_crit * num_nodes as f64).round() as usize;

    let run_results: Vec<(f64, f64, usize, usize, f64, f64, f64)> = seeds
        .into_par_iter()
        .map(|seed| {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);

            // Forward percolation.
            let result = match process {
                ProcessType::Achlioptas => {
                    achlioptas_growth_process(num_nodes, num_choices, len, alpha, &[], &mut rng)
                }
                ProcessType::Dpr => {
                    dpr_growth_process(num_nodes, num_choices, len, alpha, &[], &mut rng)
                }
            };

            let frac_conn_f = cluster_sizes_range(
                &result.r_edge,
                &result.c_edge,
                num_nodes,
                s_point_f,
                e_point_f,
            );

            let (jump_f, loc_f) = max_diff(&frac_conn_f);

            // Clamp index to valid range before accessing.
            let mean_crit_idx = mean_crit_steps.saturating_sub(s_point_f);
            let mean_crit_idx = mean_crit_idx.min(frac_conn_f.len().saturating_sub(1));
            let max_clust_mf = frac_conn_f[mean_crit_idx];
            let max_clust_sc = frac_conn_f[(loc_f + 1).min(frac_conn_f.len() - 1)];
            let fcon_end = *frac_conn_f.last().unwrap_or(&0.0);

            // Reverse percolation: flip the edge sequence.
            let r_edge_b: Vec<u32> = result.r_edge.iter().rev().cloned().collect();
            let c_edge_b: Vec<u32> = result.c_edge.iter().rev().cloned().collect();

            let frac_conn_b = cluster_sizes_range(
                &r_edge_b,
                &c_edge_b,
                num_nodes,
                s_point_b,
                e_point_b,
            );
            let (jump_b, loc_b) = max_diff(&frac_conn_b);

            (
                jump_f,
                jump_b,
                loc_f + s_point_f,   // absolute forward location
                loc_b + s_point_b,   // absolute backward location
                max_clust_mf,
                max_clust_sc,
                fcon_end,
            )
        })
        .collect();

    let n = run_results.len() as f64;
    let jump_f_vals: Vec<f64> = run_results.iter().map(|r| r.0).collect();
    let jump_b_vals: Vec<f64> = run_results.iter().map(|r| r.1).collect();
    let loc_max_jump: Vec<[usize; 2]> =
        run_results.iter().map(|r| [r.2, r.3]).collect();
    let max_clust_mean_fld: Vec<f64> = run_results.iter().map(|r| r.4).collect();
    let max_clust_self_crit: Vec<f64> = run_results.iter().map(|r| r.5).collect();
    let fully_conn = run_results.iter().filter(|r| (r.6 - 1.0).abs() < 1e-12).count();

    let mean_f = jump_f_vals.iter().sum::<f64>() / n;
    let mean_b = jump_b_vals.iter().sum::<f64>() / n;
    let std_f = stddev(&jump_f_vals, mean_f);
    let std_b = stddev(&jump_b_vals, mean_b);

    MaxJumpResult {
        crit_jump: [mean_f, mean_b],
        jump_std: [std_f, std_b],
        loc_max_jump,
        max_clust_mean_fld,
        max_clust_self_crit,
        frac_fully_conn: fully_conn as f64 / n,
    }
}

/// Results from PercAveragesRange.
pub struct AveragesResult {
    /// Average order parameter at each timestep in [s_point, e_point).
    pub fcon_means: Vec<f64>,
    /// Combined cluster-size histograms across all runs at each dist_pt.
    /// Outer index = dist_pt position, inner = histogram counts by size.
    pub cluster_hists: Vec<Vec<u32>>,
}

/// Mirrors PercAveragesRange.m.
#[allow(clippy::too_many_arguments)]
pub fn perc_averages_range(
    num_nodes: usize,
    num_choices: usize,
    s_point: usize,
    e_point: usize,
    num_runs: usize,
    process: ProcessType,
    alpha: f64,
    dist_pts: &[usize],
) -> AveragesResult {
    let master_seed: u64 = rand::random();
    let seeds: Vec<u64> = (0..num_runs as u64)
        .map(|i| master_seed.wrapping_add(i * 6364136223846793005))
        .collect();

    let run_results: Vec<(Vec<f64>, Vec<Vec<u32>>)> = seeds
        .into_par_iter()
        .map(|seed| {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
            let result = match process {
                ProcessType::Achlioptas => {
                    achlioptas_growth_process(num_nodes, num_choices, e_point, alpha, &[], &mut rng)
                }
                ProcessType::Dpr => {
                    dpr_growth_process(num_nodes, num_choices, e_point, alpha, &[], &mut rng)
                }
            };
            cluster_sizes_range_with_hist(
                &result.r_edge,
                &result.c_edge,
                num_nodes,
                s_point,
                e_point,
                dist_pts,
            )
        })
        .collect();

    let window = e_point - s_point;
    let mut fcon_sum = vec![0.0f64; window];
    let num_hist_pts = dist_pts.len();
    let mut hists_combined: Vec<Vec<u32>> = vec![Vec::new(); num_hist_pts];

    for (fcon, hists) in &run_results {
        for (i, &v) in fcon.iter().enumerate() {
            fcon_sum[i] += v;
        }
        for (hi, hist) in hists.iter().enumerate() {
            let combined = &mut hists_combined[hi];
            if hist.len() > combined.len() {
                combined.resize(hist.len(), 0);
            }
            for (j, &count) in hist.iter().enumerate() {
                combined[j] += count;
            }
        }
    }

    let n = num_runs as f64;
    let fcon_means = fcon_sum.into_iter().map(|s| s / n).collect();

    AveragesResult { fcon_means, cluster_hists: hists_combined }
}

fn stddev(v: &[f64], mean: f64) -> f64 {
    if v.len() < 2 {
        return 0.0;
    }
    let variance = v.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (v.len() - 1) as f64;
    variance.sqrt()
}
