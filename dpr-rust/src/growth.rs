use rand::Rng;
use rustc_hash::FxHashSet;

use crate::union_find::UnionFind;

pub struct GrowthResult {
    pub r_edge: Vec<u32>,
    pub c_edge: Vec<u32>,
    pub degree_snapshots: Vec<Vec<u32>>,
}

// ── Shared primitives ────────────────────────────────────────────────────────

/// Sample a uniform random pair (u, v) with u < v from [0, n).
/// Two RNG calls, zero heap allocation.
#[inline(always)]
fn sample_pair(n: usize, rng: &mut impl Rng) -> (u32, u32) {
    let u = rng.gen_range(0..n);
    let r = rng.gen_range(0..n - 1);
    let v = if r >= u { r + 1 } else { r };
    if u < v { (u as u32, v as u32) } else { (v as u32, u as u32) }
}

/// Encode an undirected edge (u < v) as a u64 for hash-set storage.
#[inline(always)]
fn edge_key(u: u32, v: u32, n: usize) -> u64 {
    u as u64 * n as u64 + v as u64
}

/// Boltzmann switching probability.
///   alpha = 0   → 0.0 (pure DPR / Achlioptas)
///   alpha = inf → 1.0 (Erdős-Rényi)
#[inline]
fn switch_prob(score_min: u64, score_max: u64, score_sum: u64, alpha: f64) -> f64 {
    if alpha == 0.0 {
        return 0.0;
    }
    if alpha.is_infinite() {
        return 1.0;
    }
    ((score_min as f64 - score_max as f64) / (score_sum as f64 * alpha)).exp()
}

/// Fill `buf` with `k` distinct candidate edges not already in the graph.
/// Checks edge_set first (O(1)), then linear-scan for intra-round duplicates
/// (cache-friendly for small k).
#[inline]
fn fill_candidates(
    buf: &mut Vec<(u32, u32)>,
    k: usize,
    n: usize,
    edge_set: &FxHashSet<u64>,
    rng: &mut impl Rng,
) {
    buf.clear();
    while buf.len() < k {
        let (u, v) = sample_pair(n, rng);
        let key = edge_key(u, v, n);
        if edge_set.contains(&key) {
            continue;
        }
        if buf.iter().any(|&(a, b)| edge_key(a, b, n) == key) {
            continue;
        }
        buf.push((u, v));
    }
}

/// Pick the winner index from `scores`, then apply Boltzmann switch if alpha ≠ 0.
#[inline]
fn pick_winner(scores: &[u64], alpha: f64, rng: &mut impl Rng) -> usize {
    let (min_idx, &min_s) = scores
        .iter()
        .enumerate()
        .min_by_key(|&(_, &s)| s)
        .unwrap();

    if alpha == 0.0 || scores.len() < 2 {
        return min_idx;
    }
    let max_s = *scores.iter().max().unwrap();
    let sum_s: u64 = scores.iter().sum();
    if rng.gen::<f64>() < switch_prob(min_s, max_s, sum_s, alpha) {
        rng.gen_range(0..scores.len())
    } else {
        min_idx
    }
}

// ── DPR Growth Process ───────────────────────────────────────────────────────

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
    let mut edge_set = FxHashSet::with_capacity_and_hasher(len, Default::default());

    // Both buffers pre-allocated once and reused across every step.
    let mut cands: Vec<(u32, u32)> = Vec::with_capacity(num_choices);
    let mut scores: Vec<u64> = Vec::with_capacity(num_choices);

    for step in 0..len {
        let (wu, wv) = if num_choices == 2 {
            // ── k=2 fast path ────────────────────────────────────────────────
            // Zero Vec allocations: two loops with O(1) hash lookups, then
            // one integer comparison.  k0 is cached to avoid recomputing it
            // inside the second loop.
            let (u0, v0) = loop {
                let p = sample_pair(num_nodes, rng);
                if !edge_set.contains(&edge_key(p.0, p.1, num_nodes)) {
                    break p;
                }
            };
            let k0 = edge_key(u0, v0, num_nodes);
            let (u1, v1) = loop {
                let p = sample_pair(num_nodes, rng);
                let k = edge_key(p.0, p.1, num_nodes);
                if k != k0 && !edge_set.contains(&k) {
                    break p;
                }
            };

            let s0 = degree[u0 as usize] as u64 * degree[v0 as usize] as u64;
            let s1 = degree[u1 as usize] as u64 * degree[v1 as usize] as u64;

            let use_first = if s0 <= s1 {
                alpha == 0.0 || rng.gen::<f64>() >= switch_prob(s0, s1, s0 + s1, alpha)
            } else {
                alpha != 0.0 && rng.gen::<f64>() < switch_prob(s1, s0, s0 + s1, alpha)
            };

            if use_first { (u0, v0) } else { (u1, v1) }
        } else {
            // ── General path (k > 2) ─────────────────────────────────────────
            fill_candidates(&mut cands, num_choices, num_nodes, &edge_set, rng);
            scores.clear();
            scores.extend(
                cands.iter().map(|&(u, v)| degree[u as usize] as u64 * degree[v as usize] as u64),
            );
            let w = pick_winner(&scores, alpha, rng);
            cands[w]
        };

        r_edge.push(wu);
        c_edge.push(wv);
        edge_set.insert(edge_key(wu, wv, num_nodes));
        degree[wu as usize] += 1;
        degree[wv as usize] += 1;

        if deg_dist_pts.contains(&step) {
            degree_snapshots.push(degree.clone());
        }
    }

    GrowthResult { r_edge, c_edge, degree_snapshots }
}

// ── Achlioptas Growth Process ────────────────────────────────────────────────

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
    let mut uf = UnionFind::new(num_nodes);
    let mut edge_set = FxHashSet::with_capacity_and_hasher(len, Default::default());

    let mut cands: Vec<(u32, u32)> = Vec::with_capacity(num_choices);
    let mut scores: Vec<u64> = Vec::with_capacity(num_choices);

    for step in 0..len {
        let (wu, wv) = if num_choices == 2 {
            let (u0, v0) = loop {
                let p = sample_pair(num_nodes, rng);
                if !edge_set.contains(&edge_key(p.0, p.1, num_nodes)) {
                    break p;
                }
            };
            let k0 = edge_key(u0, v0, num_nodes);
            let (u1, v1) = loop {
                let p = sample_pair(num_nodes, rng);
                let k = edge_key(p.0, p.1, num_nodes);
                if k != k0 && !edge_set.contains(&k) {
                    break p;
                }
            };

            let s0 = uf.size_of(u0 as usize) as u64 * uf.size_of(v0 as usize) as u64;
            let s1 = uf.size_of(u1 as usize) as u64 * uf.size_of(v1 as usize) as u64;

            let use_first = if s0 <= s1 {
                alpha == 0.0 || rng.gen::<f64>() >= switch_prob(s0, s1, s0 + s1, alpha)
            } else {
                alpha != 0.0 && rng.gen::<f64>() < switch_prob(s1, s0, s0 + s1, alpha)
            };

            if use_first { (u0, v0) } else { (u1, v1) }
        } else {
            fill_candidates(&mut cands, num_choices, num_nodes, &edge_set, rng);
            scores.clear();
            scores.extend(
                cands
                    .iter()
                    .map(|&(u, v)| uf.size_of(u as usize) as u64 * uf.size_of(v as usize) as u64),
            );
            let w = pick_winner(&scores, alpha, rng);
            cands[w]
        };

        r_edge.push(wu);
        c_edge.push(wv);
        edge_set.insert(edge_key(wu, wv, num_nodes));
        uf.union(wu as usize, wv as usize);
        degree[wu as usize] += 1;
        degree[wv as usize] += 1;

        if deg_dist_pts.contains(&step) {
            degree_snapshots.push(degree.clone());
        }
    }

    GrowthResult { r_edge, c_edge, degree_snapshots }
}
