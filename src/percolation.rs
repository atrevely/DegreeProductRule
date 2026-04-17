use crate::union_find::UnionFind;

/// Compute the order parameter (largest cluster / num_nodes) over the window
/// [s_point, e_point).
///
/// Returns a Vec of length `e_point - s_point`:
///   result[0]  = state after the first s_point edges have been added
///   result[k]  = state after s_point + k edges have been added
///   result[e_point - s_point - 1] = state after e_point - 1 edges
///
/// This mirrors ClusterSizesRange.m but uses Union-Find (O(α(n)) per edge)
/// instead of the O(n) cluster-label scan in the MATLAB version.
pub fn cluster_sizes_range(
    r_edge: &[u32],
    c_edge: &[u32],
    num_nodes: usize,
    s_point: usize,
    e_point: usize,
) -> Vec<f64> {
    debug_assert!(e_point > s_point);
    debug_assert!(e_point <= r_edge.len() + 1);

    let len = e_point - s_point;
    let mut result = Vec::with_capacity(len);
    let mut uf = UnionFind::new(num_nodes);

    for i in 0..s_point {
        uf.union(r_edge[i] as usize, c_edge[i] as usize);
    }
    result.push(uf.max_size as f64 / num_nodes as f64);

    for i in s_point..(e_point - 1) {
        uf.union(r_edge[i] as usize, c_edge[i] as usize);
        result.push(uf.max_size as f64 / num_nodes as f64);
    }

    result
}

/// Same as `cluster_sizes_range` but also records the full cluster-size
/// histogram at the timesteps in `dist_pts` (0-indexed relative to the full
/// edge sequence, not the window).
///
/// Returns (order_parameter_vec, histograms) where each histogram is a
/// Vec<u32> with histogram[k] = number of clusters of size k+1.
pub fn cluster_sizes_range_with_hist(
    r_edge: &[u32],
    c_edge: &[u32],
    num_nodes: usize,
    s_point: usize,
    e_point: usize,
    dist_pts: &[usize],
) -> (Vec<f64>, Vec<Vec<u32>>) {
    let len = e_point - s_point;
    let mut result = Vec::with_capacity(len);
    let mut histograms: Vec<Vec<u32>> = Vec::new();
    let mut uf = UnionFind::new(num_nodes);

    for i in 0..s_point {
        uf.union(r_edge[i] as usize, c_edge[i] as usize);
    }
    result.push(uf.max_size as f64 / num_nodes as f64);
    if dist_pts.contains(&s_point) {
        histograms.push(build_histogram(&mut uf, num_nodes));
    }

    for i in s_point..(e_point - 1) {
        uf.union(r_edge[i] as usize, c_edge[i] as usize);
        result.push(uf.max_size as f64 / num_nodes as f64);
        let step = i + 1;
        if dist_pts.contains(&step) {
            histograms.push(build_histogram(&mut uf, num_nodes));
        }
    }

    (result, histograms)
}

/// Build a histogram of cluster sizes from the current Union-Find state.
/// histogram[k] = number of clusters of size k+1.
fn build_histogram(uf: &mut UnionFind, num_nodes: usize) -> Vec<u32> {
    let mut size_counts: Vec<u32> = Vec::new();
    for node in 0..num_nodes {
        if uf.find(node) == node {
            let s = uf.size_of(node) as usize;
            if s > size_counts.len() {
                size_counts.resize(s, 0);
            }
            size_counts[s - 1] += 1;
        }
    }
    size_counts
}

/// Find the index and value of the maximum difference between consecutive
/// elements (equivalent to MATLAB's `max(diff(v))`).
/// Returns (max_diff, index_of_max_diff).
pub fn max_diff(v: &[f64]) -> (f64, usize) {
    v.windows(2)
        .enumerate()
        .map(|(i, w)| (w[1] - w[0], i))
        .max_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or((0.0, 0))
}
