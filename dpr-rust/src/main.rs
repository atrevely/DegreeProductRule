mod global_choice;
mod growth;
mod percolation;
mod statistics;
mod union_find;

use clap::{Parser, Subcommand, ValueEnum};
use std::io::Write as IoWrite;
use statistics::ProcessType;

// ---------------------------------------------------------------------------
// CLI definition
// ---------------------------------------------------------------------------

#[derive(Parser)]
#[command(
    name = "dpr",
    about = "Degree Product Rule / Achlioptas network percolation simulator",
    long_about = "\
Rust rewrite of the MATLAB DPR percolation suite published in Phys. Rev. E.\n\
Parallelises ensemble runs with Rayon; uses Union-Find for O(α(n)) cluster ops.\n\
\n\
Subcommands:\n\
  perc-avg    — PercAveragesRange: averaged order parameter + cluster histograms\n\
  perc-max    — PercMaxJumpMeans: max jump statistics (forward + reverse perc.)\n\
  global      — DPRGlobalChoice: DPR in the global-choice (analytical) limit\n\
  grow        — Single DPR or Achlioptas growth realisation, dump edge list"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Clone, Copy, ValueEnum)]
enum CliProcess {
    Dpr,
    Achlioptas,
}

impl From<CliProcess> for ProcessType {
    fn from(p: CliProcess) -> Self {
        match p {
            CliProcess::Dpr => ProcessType::Dpr,
            CliProcess::Achlioptas => ProcessType::Achlioptas,
        }
    }
}

#[derive(Subcommand)]
enum Commands {
    /// Average order parameter over many independent realisations.
    PercAvg {
        /// Number of nodes
        #[arg(short = 'N', long, default_value_t = 10_000)]
        num_nodes: usize,
        /// Candidate edges per step
        #[arg(short = 'k', long, default_value_t = 2)]
        num_choices: usize,
        /// Window start (timestep)
        #[arg(long, default_value_t = 4_000)]
        s_point: usize,
        /// Window end (timestep, exclusive)
        #[arg(long, default_value_t = 7_000)]
        e_point: usize,
        /// Independent realisations
        #[arg(short = 'r', long, default_value_t = 100)]
        num_runs: usize,
        /// Growth process
        #[arg(short = 'p', long, value_enum, default_value = "dpr")]
        process: CliProcess,
        /// Alpha (0 = pure DPR/Achlioptas, inf = Erdős-Rényi)
        #[arg(short = 'a', long, default_value_t = 0.0)]
        alpha: f64,
        /// Comma-separated timesteps at which to record cluster-size histograms
        #[arg(long, default_value = "")]
        dist_pts: String,
        /// Output CSV path for order parameter
        #[arg(short = 'o', long, default_value = "order_param.csv")]
        output: String,
        /// Output prefix for cluster histogram CSVs (omit to skip)
        #[arg(long, default_value = "")]
        hist_prefix: String,
    },

    /// Max jump in order parameter statistics (forward + reverse percolation).
    PercMax {
        #[arg(short = 'N', long, default_value_t = 10_000)]
        num_nodes: usize,
        #[arg(short = 'k', long, default_value_t = 2)]
        num_choices: usize,
        /// Forward window start
        #[arg(long, default_value_t = 4_000)]
        s_point_f: usize,
        /// Forward window end
        #[arg(long, default_value_t = 7_000)]
        e_point_f: usize,
        /// Backward window start
        #[arg(long, default_value_t = 3_000)]
        s_point_b: usize,
        /// Backward window end
        #[arg(long, default_value_t = 7_000)]
        e_point_b: usize,
        #[arg(short = 'r', long, default_value_t = 100)]
        num_runs: usize,
        /// Total simulation length (edges added)
        #[arg(short = 'l', long, default_value_t = 8_000)]
        len: usize,
        /// Mean-field critical point as a fraction of num_nodes (e.g. 0.5)
        #[arg(long, default_value_t = 0.5)]
        mean_crit: f64,
        #[arg(short = 'p', long, value_enum, default_value = "dpr")]
        process: CliProcess,
        #[arg(short = 'a', long, default_value_t = 0.0)]
        alpha: f64,
        #[arg(short = 'o', long, default_value = "max_jump.csv")]
        output: String,
    },

    /// DPR Global Choice (analytical limit, no random edge sampling).
    Global {
        #[arg(short = 'N', long, default_value_t = 10_000)]
        num_nodes: usize,
        #[arg(short = 'r', long, default_value_t = 1)]
        num_runs: usize,
        #[arg(short = 'o', long, default_value = "global_choice.csv")]
        output: String,
    },

    /// Single growth realisation — dump ordered edge list as CSV.
    Grow {
        #[arg(short = 'N', long, default_value_t = 10_000)]
        num_nodes: usize,
        #[arg(short = 'k', long, default_value_t = 2)]
        num_choices: usize,
        #[arg(short = 'l', long, default_value_t = 6_000)]
        len: usize,
        #[arg(short = 'p', long, value_enum, default_value = "dpr")]
        process: CliProcess,
        #[arg(short = 'a', long, default_value_t = 0.0)]
        alpha: f64,
        /// Comma-separated 0-indexed timesteps for degree snapshots
        #[arg(long, default_value = "")]
        deg_pts: String,
        #[arg(short = 'o', long, default_value = "edges.csv")]
        output: String,
    },
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn parse_usize_list(s: &str) -> Vec<usize> {
    if s.is_empty() {
        return Vec::new();
    }
    s.split(',')
        .filter_map(|t| t.trim().parse::<usize>().ok())
        .collect()
}

fn open_csv(path: &str) -> csv::Writer<std::fs::File> {
    csv::Writer::from_path(path)
        .unwrap_or_else(|e| panic!("Cannot open {path}: {e}"))
}

fn eprint_progress(msg: &str) {
    eprintln!("[dpr] {msg}");
    let _ = std::io::stderr().flush();
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

fn main() {
    let cli = Cli::parse();

    match cli.command {
        // ------------------------------------------------------------------ //
        Commands::PercAvg {
            num_nodes,
            num_choices,
            s_point,
            e_point,
            num_runs,
            process,
            alpha,
            dist_pts,
            output,
            hist_prefix,
        } => {
            let dist_pts = parse_usize_list(&dist_pts);
            eprint_progress(&format!(
                "perc-avg: N={num_nodes} k={num_choices} s={s_point} e={e_point} \
                 runs={num_runs} process={:?} alpha={alpha}",
                process as u8
            ));

            let result = statistics::perc_averages_range(
                num_nodes,
                num_choices,
                s_point,
                e_point,
                num_runs,
                process.into(),
                alpha,
                &dist_pts,
            );

            // Write order parameter CSV.
            let mut wtr = open_csv(&output);
            wtr.write_record(["timestep", "frac_largest_cluster"]).unwrap();
            for (i, &v) in result.fcon_means.iter().enumerate() {
                wtr.write_record(&[(s_point + i).to_string(), format!("{v:.8}")])
                    .unwrap();
            }
            wtr.flush().unwrap();
            eprintln!("[dpr] Order parameter written to {output}");

            // Write cluster histograms.
            if !hist_prefix.is_empty() {
                for (hi, (hist, &pt)) in
                    result.cluster_hists.iter().zip(dist_pts.iter()).enumerate()
                {
                    let path = format!("{hist_prefix}_t{pt}_{hi}.csv");
                    let mut wtr = open_csv(&path);
                    wtr.write_record(["cluster_size", "count"]).unwrap();
                    for (size_minus1, &count) in hist.iter().enumerate() {
                        wtr.write_record(&[
                            (size_minus1 + 1).to_string(),
                            count.to_string(),
                        ])
                        .unwrap();
                    }
                    wtr.flush().unwrap();
                    eprintln!("[dpr] Histogram written to {path}");
                }
            }
        }

        // ------------------------------------------------------------------ //
        Commands::PercMax {
            num_nodes,
            num_choices,
            s_point_f,
            e_point_f,
            s_point_b,
            e_point_b,
            num_runs,
            len,
            mean_crit,
            process,
            alpha,
            output,
        } => {
            eprint_progress(&format!(
                "perc-max: N={num_nodes} k={num_choices} runs={num_runs} \
                 process={:?} alpha={alpha}",
                process as u8
            ));

            let result = statistics::perc_max_jump_means(
                num_nodes,
                num_choices,
                s_point_f,
                e_point_f,
                s_point_b,
                e_point_b,
                num_runs,
                len,
                mean_crit,
                process.into(),
                alpha,
            );

            // Print summary.
            println!("=== PercMaxJumpMeans results ===");
            println!(
                "Forward  mean jump = {:.6}  std = {:.6}",
                result.crit_jump[0], result.jump_std[0]
            );
            println!(
                "Backward mean jump = {:.6}  std = {:.6}",
                result.crit_jump[1], result.jump_std[1]
            );
            println!("Fraction fully connected = {:.4}", result.frac_fully_conn);

            // Write per-run CSV.
            let mut wtr = open_csv(&output);
            wtr.write_record([
                "run",
                "jump_forward",
                "jump_backward",
                "loc_forward",
                "loc_backward",
                "max_clust_mean_fld",
                "max_clust_self_crit",
            ])
            .unwrap();
            for (i, locs) in result.loc_max_jump.iter().enumerate() {
                wtr.write_record(&[
                    i.to_string(),
                    format!("{:.8}", result.crit_jump[0]),
                    format!("{:.8}", result.crit_jump[1]),
                    locs[0].to_string(),
                    locs[1].to_string(),
                    format!("{:.8}", result.max_clust_mean_fld[i]),
                    format!("{:.8}", result.max_clust_self_crit[i]),
                ])
                .unwrap();
            }
            wtr.flush().unwrap();
            eprintln!("[dpr] Per-run stats written to {output}");
        }

        // ------------------------------------------------------------------ //
        Commands::Global { num_nodes, num_runs, output } => {
            eprint_progress(&format!(
                "global-choice: N={num_nodes} runs={num_runs}"
            ));

            let mut wtr = open_csv(&output);
            wtr.write_record(["run", "timestep", "max_clust", "max_jump", "no_finish"])
                .unwrap();

            for run in 0..num_runs {
                let mut rng = {
                    use rand::SeedableRng;
                    rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(rand::random())
                };
                let result = global_choice::dpr_global_choice_full(num_nodes, &mut rng);
                for (t, &v) in result.max_clust.iter().enumerate() {
                    wtr.write_record(&[
                        run.to_string(),
                        t.to_string(),
                        format!("{v:.8}"),
                        format!("{:.8}", result.max_jump),
                        (result.no_finish as u8).to_string(),
                    ])
                    .unwrap();
                }
                if (run + 1) % 10 == 0 || run + 1 == num_runs {
                    eprint_progress(&format!("  completed {}/{num_runs} runs", run + 1));
                }
            }
            wtr.flush().unwrap();
            eprintln!("[dpr] Global choice results written to {output}");
        }

        // ------------------------------------------------------------------ //
        Commands::Grow {
            num_nodes,
            num_choices,
            len,
            process,
            alpha,
            deg_pts,
            output,
        } => {
            let deg_pts = parse_usize_list(&deg_pts);
            eprint_progress(&format!(
                "grow: N={num_nodes} k={num_choices} len={len} process={:?} alpha={alpha}",
                process as u8
            ));

            let mut rng = {
                use rand::SeedableRng;
                rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(rand::random())
            };

            let result = match process.into() {
                ProcessType::Dpr => growth::dpr_growth_process(
                    num_nodes, num_choices, len, alpha, &deg_pts, &mut rng,
                ),
                ProcessType::Achlioptas => growth::achlioptas_growth_process(
                    num_nodes, num_choices, len, alpha, &deg_pts, &mut rng,
                ),
            };

            let mut wtr = open_csv(&output);
            wtr.write_record(["timestep", "node_u", "node_v"]).unwrap();
            for (t, (&u, &v)) in
                result.r_edge.iter().zip(result.c_edge.iter()).enumerate()
            {
                wtr.write_record(&[t.to_string(), u.to_string(), v.to_string()])
                    .unwrap();
            }
            wtr.flush().unwrap();
            eprintln!("[dpr] Edge list written to {output}");

            // Degree snapshots to stdout if requested.
            if !result.degree_snapshots.is_empty() {
                println!("degree_snapshot_index,node,degree");
                for (si, snap) in result.degree_snapshots.iter().enumerate() {
                    for (node, &deg) in snap.iter().enumerate() {
                        println!("{si},{node},{deg}");
                    }
                }
            }
        }
    }
}
