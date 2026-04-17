use rand::Rng;
use rand::seq::index::sample;
use std::collections::HashSet;

use crate::union_find::UnionFind;

pub struct GrowthResult {
    pub r_edge: Vec<u32>,
    pub c_edge: Vec<u32>,
    /// Degree snapshots at requested timesteps (one row per snapshot point).
    pub degree_snapshots: Vec<Vec<u32>>,
}

/// Encode an undirected edge (u < v) as a single u64 key.
#[inline(always)]
fn edge_key(u: usize, v: usize, n: usize) -> u64 {
    debug_assert!(u < v);
    u as u64 * n as u64 + v as u64
}

/// Compute the alpha-based switching probability.
/// Returns the probability of picking a random candidate instead of the minimum.
/// alpha=0 → 0.0 (pure DPR), alpha=+inf → 1.0 (Erdős-Rényi).
#[inline]
fn switch_prob(scores: &[f64], alpha: f64) -> f64 {
    if alpha == 0.0 {
        return 0.0;
    }
    if alpha.is_infinite() {
        return 1.0;
    }
    let min = scores.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = scores.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let kappa: f64 = scores.iter().sum();
    ((min - max) / (kappa * alpha)).exp()
}

/// Sample `num_choices` distinct undirected edges that are not already in `edge_set`.
/// Returns a Vec of (u, v) pairs with u < v.
fn sample_candidates(
    num_nodes: usize,
    num_choices: usize,
    remaining_edges: usize,
    edge_set: &HashSet<u64>,
    rng: &mut impl Rng,
) -> Vec<(usize, usize)> {
    let cap = num_choices.min(remaining_edges);
    let mut candidates = Vec::with_capacity(cap);
    let mut candidate_keys = HashSet::with_capacity(cap);
    while candidates.len() < cap {
        let idx = sample(rng, num_nodes, 2);
        let a = idx.index(0);
        let b = idx.index(1);
        let (u, v) = if a < b { (a, b) } else { (b, a) };
        let key = edge_key(u, v, num_nodes);
        if !edge_set.contains(&key) && candidate_keys.insert(key) {
            candidates.push((u, v));
        }
    }
    candidates
}

/// DPR Growth Process — faithfully mirrors DPRGrowthProcess.m.
///
/// * Degrees are initialised to 1 (the MATLAB convention).
/// * At each step, `num_choices` candidate edges compete; the pair with
///   the minimum degree product wins (with alpha-based randomisation).
/// * `deg_dist_pts`: 0-indexed timesteps at which to snapshot degrees.
///   Pass an empty slice to skip all snapshots.
pub fn dpr_growth_process(
    num_nodes: usize,
    num_choices: usize,
    len: usize,
    alpha: f64,
    deg_dist_pts: &[usize],
    rng: &mut impl Rng,
) -> GrowthResult {
    let mut degree = vec![1u32; num_nodes];
    let mut r_edge = Vec::with_capacity(len);
    let mut c_edge = Vec::with_capacity(len);
    let mut degree_snapshots = Vec::new();
    let mut edge_set: HashSet<u64> = HashSet::with_capacity(len + num_choices);

    // Maximum possible undirected edges (no self-loops).
    let max_edges = num_nodes * (num_nodes - 1) / 2;

    for step in 0..len {
        let remaining = max_edges.saturating_sub(step);
        let candidates = sample_candidates(num_nodes, num_choices, remaining, &edge_set, rng);
        let n_cands = candidates.len();

        // Score each candidate by degree product.
        let scores: Vec<f64> = candidates
            .iter()
            .map(|&(u, v)| (degree[u] as f64) * (degree[v] as f64))
            .collect();

        // Find minimum-score candidate.
        let mut winner = scores
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0);

        // Probabilistic switch to a random candidate.
        if n_cands > 1 && rng.gen::<f64>() < switch_prob(&scores, alpha) {
            winner = rng.gen_range(0..n_cands);
        }

        let (u, v) = candidates[winner];
        r_edge.push(u as u32);
        c_edge.push(v as u32);

        // Add winner to the permanent edge set; candidates not chosen are NOT kept.
        edge_set.insert(edge_key(u, v, num_nodes));

        degree[u] += 1;
        degree[v] += 1;

        if deg_dist_pts.contains(&step) {
            degree_snapshots.push(degree.clone());
        }
    }

    GrowthResult { r_edge, c_edge, degree_snapshots }
}

/// Achlioptas Growth Process — faithfully mirrors AchlioptasGrowthProcess.m.
///
/// Selection criterion: minimum product of *cluster sizes* (not degrees).
/// Uses Union-Find for O(α(n)) cluster queries instead of MATLAB's O(n) linear scan.
pub fn achlioptas_growth_process(
    num_nodes: usize,
    num_choices: usize,
    len: usize,
    alpha: f64,
    deg_dist_pts: &[usize],
    rng: &mut impl Rng,
) -> GrowthResult {
    let mut degree = vec![1u32; num_nodes];
    let mut r_edge = Vec::with_capacity(len);
    let mut c_edge = Vec::with_capacity(len);
    let mut degree_snapshots = Vec::new();
    let mut edge_set: HashSet<u64> = HashSet::with_capacity(len + num_choices);
    let mut uf = UnionFind::new(num_nodes);

    let max_edges = num_nodes * (num_nodes - 1) / 2;

    for step in 0..len {
        let remaining = max_edges.saturating_sub(step);
        let candidates = sample_candidates(num_nodes, num_choices, remaining, &edge_set, rng);
        let n_cands = candidates.len();

        // Score by cluster-size product.
        let scores: Vec<f64> = candidates
            .iter()
            .map(|&(u, v)| (uf.size_of(u) as f64) * (uf.size_of(v) as f64))
            .collect();

        let mut winner = scores
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0);

        if n_cands > 1 && rng.gen::<f64>() < switch_prob(&scores, alpha) {
            winner = rng.gen_range(0..n_cands);
        }

        let (u, v) = candidates[winner];
        r_edge.push(u as u32);
        c_edge.push(v as u32);

        edge_set.insert(edge_key(u, v, num_nodes));
        uf.union(u, v);

        degree[u] += 1;
        degree[v] += 1;

        if deg_dist_pts.contains(&step) {
            degree_snapshots.push(degree.clone());
        }
    }

    GrowthResult { r_edge, c_edge, degree_snapshots }
}
