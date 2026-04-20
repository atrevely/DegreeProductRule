# Degree Product Rule — Percolation Simulator

Simulation suite for the **Degree Product Rule (DPR)** and **Achlioptas** explosive-percolation processes on Erdős-Rényi-style random graphs, supporting finite-size scaling studies.

Based on the paper: *Rapid Explosive Percolation via the Degree Product Rule* (see `DegreeProductRulePRERapid.pdf`).

---

## Repository layout

```
dpr-original/   MATLAB reference implementation
dpr-rust/       High-performance Rust rewrite (this is what you should run)
```

---

## dpr-original — MATLAB reference

The original single-threaded MATLAB code. Requires the Parallel Computing Toolbox for `parfor` support.

| File | Description |
|------|-------------|
| `DPRGrowthProcess.m` | Grows a network under the Degree Product Rule: at each step, `k` candidate edges compete and the one with the smallest degree-product is added. An `alpha` parameter interpolates continuously to Erdős-Rényi. |
| `AchlioptasGrowthProcess.m` | Same process but scores edges by their connected-cluster sizes (product rule on component sizes). |
| `DPRGlobalChoice.m` | Analytical limit: considers all possible edges at each step rather than a random sample of `k`. |
| `ClusterSizesRange.m` | Computes the order parameter (largest cluster / N) over a window of timesteps from an edge list. |
| `PercAveragesRange.m` | Averages the order parameter over many `parfor`-parallelised realisations. |
| `PercMaxJumpMeans.m` | Locates the maximum jump in the order parameter for both forward and reverse percolation, averaged across runs. Seeds are generated outside `parfor` to avoid repeated RNG state. |

### Limitations of the MATLAB version

- O(N) cluster tracking via linear `nodeCluster` scans
- Single-node sparse-array edge deduplication is slow for large N
- `parfor` parallelism is coarse and requires a licence

---

## dpr-rust — Rust rewrite

A from-scratch Rust implementation that achieves **~2× higher throughput** than the MATLAB version while producing statistically identical results.

### Key algorithmic improvements

| Aspect | MATLAB | Rust |
|--------|--------|------|
| Cluster tracking | O(N) linear scan | O(α(N)) path-halving Union-Find |
| Edge deduplication | Cell array of sparse vectors | `FxHashSet<u64>` (rustc-hash) |
| Parallelism | `parfor` (Toolbox required) | Rayon work-stealing, no licence needed |
| RNG | `rand` (shared state) | Xoshiro256++ per thread, seeded outside loop |
| k=2 fast path | General loop | Inline, zero heap allocation per step |

### Build

Requires [Rust](https://rustup.rs/) (stable).

```bash
cd dpr-rust
cargo build --release
# binary: dpr-rust/target/release/dpr
```

### Subcommands

#### `perc-max` — maximum jump statistics (finite-size scaling)

Runs forward and reverse percolation across many realisations and returns the mean largest jump in the order parameter. This is the primary subcommand for scaling studies.

```
dpr perc-max -N <nodes> -k <choices> -r <runs> -l <len>
             --s-point-f <start> --e-point-f <end>
             --s-point-b <start> --e-point-b <end>
             --mean-crit <fraction> -p <dpr|achlioptas> -a <alpha>
             -o <output.csv>
```

Key flags:
- `-N` — number of nodes
- `-k` — candidate edges per step (2 = standard DPR/Achlioptas)
- `-r` — independent realisations (parallelised)
- `-l` — total edges to grow before reversing
- `--s-point-f/--e-point-f` — timestep window for forward jump search
- `--s-point-b/--e-point-b` — timestep window for reverse jump search
- `--mean-crit` — critical point as a fraction of N (DPR ≈ 0.763)
- `-p` — process: `dpr` or `achlioptas`
- `-a` — alpha: 0 = pure process, ∞ = Erdős-Rényi

Output CSV columns: `run, jump_forward, jump_backward, loc_forward, loc_backward, max_clust_mean_fld, max_clust_self_crit`

#### `perc-avg` — averaged order parameter curve

Computes the mean order parameter (largest cluster / N) at each timestep, averaged over many realisations. Optionally writes cluster-size histograms.

```
dpr perc-avg -N <nodes> -k <choices> -r <runs>
             --s-point <start> --e-point <end>
             -p <dpr|achlioptas> -a <alpha> -o <output.csv>
             [--dist-pts <t1,t2,...> --hist-prefix <prefix>]
```

#### `grow` — single growth realisation

Grows one network and writes the ordered edge list as CSV. Useful for inspecting individual realisations or debugging.

```
dpr grow -N <nodes> -k <choices> -l <len> -p <dpr|achlioptas> -o <edges.csv>
         [--deg-pts <t1,t2,...>]   # optional degree snapshots at given timesteps
```

#### `global` — global-choice limit

Runs the DPR process in the analytical (global-choice) limit where all possible edges are considered at each step.

```
dpr global -N <nodes> -r <runs> -o <output.csv>
```

---

## Bash scripts

Both scripts live in `dpr-rust/` and can be run from any directory — they locate the binary relative to their own path.

### `run_scaling.sh` — finite-size scaling sweep

Sweeps N from 100 to 10⁸ on a log scale (one point per order of magnitude and one halfway between, i.e. at √10 ≈ 3.162 intervals). At each N it runs `perc-max`, extracts S (raw order-parameter jump) and C/N (largest cluster fraction at the critical point), and writes `scaling.csv`.

Thread count is automatically capped based on available RAM (≈ 30 bytes × N per thread) to prevent out-of-memory errors at large N.

```bash
cd dpr-rust
./run_scaling.sh                    # uses ./target/release/dpr
./run_scaling.sh /path/to/dpr      # custom binary
```

Output: `scaling.csv` with columns `N, num_runs, threads, S, C_over_N, C_over_N_std, time_per_run_s`

### `run_benchmark.sh` — throughput benchmark

Sweeps N = 1 k → 10 M and reports wall-clock time per run, edge throughput (M edges/s), and estimated memory per thread.

```bash
cd dpr-rust
./run_benchmark.sh
```

Output: printed table + `benchmark.csv`

---

## Running on AWS

The Rust binary parallelises automatically across all available cores via Rayon. For large-N scaling studies an instance with high core count and sufficient RAM is recommended:

- **c7i.8xlarge** (32 vCPUs, 64 GB) — good for N up to ~10 M per thread
- **c7i.16xlarge** (64 vCPUs, 128 GB) — for N up to 10⁸

The thread count is capped automatically by `run_scaling.sh` based on `/proc/meminfo`, so there is no risk of OOM from over-threading. You can also set `RAYON_NUM_THREADS=<n>` manually before running the binary.
