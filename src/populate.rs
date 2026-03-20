/// Batch-compute Ehrhart data for GT polytopes and store in MariaDB.
///
/// Iterates over all partitions λ with |λ| in [min_n, max_n],
/// all partition weights w with λ ≥ w in dominance order,
/// and optionally inner shapes μ (--skew) and flag vectors (--flag).
///
/// Results are stored in a MariaDB table; already-computed rows are skipped
/// on restart, enabling crash-resume.

use clap::Args;
use mysql::prelude::*;
use mysql::*;
use num_bigint::BigInt;
use num_traits::Zero;
use std::collections::HashSet;
use std::panic;
use std::time::Instant;

use crate::ehrhart::{compute_ehrhart, compute_hstar};
use crate::gt_dim::gt_polytope_dim_full;
use crate::kostka_dp::{flagged_skew_kostka, skew_kostka};
use crate::partition::Partition;

// ── CLI args ─────────────────────────────────────────────────────────────────

#[derive(Args)]
pub struct PopulateArgs {
    /// Maximum partition size |λ| to compute
    #[arg(long)]
    pub max_n: u32,

    /// Minimum partition size |λ| (default: 1)
    #[arg(long, default_value = "1")]
    pub min_n: u32,

    /// MariaDB connection URL, e.g. mysql://user:pass@localhost/dbname
    /// Can also be set via KOSTKA_DB_URL env var
    #[arg(long, env = "KOSTKA_DB_URL")]
    pub db_url: String,

    /// Also iterate over inner partitions μ ⊂ λ (skew shapes)
    #[arg(long)]
    pub skew: bool,

    /// Also iterate over weakly increasing upper flag vectors
    #[arg(long)]
    pub flag: bool,

    /// Halt and report if any Ehrhart polynomial has a negative coefficient
    #[arg(long)]
    pub kkt_counterexample: bool,

    /// Abort DP if any level exceeds N states (prevents OOM)
    #[arg(long)]
    pub max_states: Option<usize>,

    /// Skip entries where GT polytope dimension exceeds this value
    #[arg(long)]
    pub max_dim: Option<usize>,

    /// Skip entries where GT polytope dimension is below this value
    #[arg(long)]
    pub min_dim: Option<usize>,

    /// Dry run: enumerate work items without computing
    #[arg(long)]
    pub dry_run: bool,
}

// ── schema ───────────────────────────────────────────────────────────────────

const CREATE_TABLE: &str = r"
CREATE TABLE IF NOT EXISTS gt_ehrhart (
    id              BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    n               INT UNSIGNED NOT NULL,
    lambda          VARCHAR(255) NOT NULL,
    mu              VARCHAR(255) NOT NULL DEFAULT '',
    weight          VARCHAR(255) NOT NULL,
    upper_flags     VARCHAR(255) NOT NULL DEFAULT '',
    lower_flags     VARCHAR(255) NOT NULL DEFAULT '',
    kostka          TEXT NOT NULL,
    dimension       INT,
    polynomial      TEXT,
    hstar           TEXT,
    elapsed_ms      BIGINT UNSIGNED,
    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE KEY uq_shape (lambda, mu, weight, upper_flags, lower_flags)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
";

// ── helpers ──────────────────────────────────────────────────────────────────

fn parts_str(p: &Partition) -> String {
    p.parts().iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}

fn vec_str(v: &[u32]) -> String {
    v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}

fn hstar_json(h: &[BigInt]) -> String {
    let elems: Vec<String> = h.iter().map(|x| format!("\"{}\"", x)).collect();
    format!("[{}]", elems.join(","))
}

/// True if λ ≥ w in dominance order.
fn dominates(lam: &Partition, w: &[u32]) -> bool {
    let mut s_lam: u32 = 0;
    let mut s_w: u32 = 0;
    let max_len = lam.num_parts().max(w.len());
    for i in 0..max_len {
        s_lam += lam.part(i);
        s_w += if i < w.len() { w[i] } else { 0 };
        if s_lam < s_w {
            return false;
        }
    }
    true
}

/// Check that every row of λ/μ is non-empty and the diagram is connected.
fn valid_skew(lam: &Partition, mu: &Partition) -> bool {
    let n = lam.num_parts();
    // Every row non-empty: λ_i > μ_i for all i.
    for i in 0..n {
        if lam.part(i) <= mu.part(i) {
            return false;
        }
    }
    // Connected: column ranges of consecutive rows must overlap, i.e. μ_i < λ_{i+1}.
    for i in 0..n.saturating_sub(1) {
        if mu.part(i) >= lam.part(i + 1) {
            return false;
        }
    }
    true
}

/// All μ with μ ⊂ λ such that λ/μ has no empty rows and is connected.
/// Always includes the empty partition (straight shape).
fn inner_shapes(lambda: &Partition) -> Vec<Partition> {
    let mut result = vec![Partition::empty()];
    for s in 1..lambda.size() {
        for mu in Partition::all_of_size_bounded(
            s,
            lambda.num_parts(),
            lambda.parts()[0],
        ) {
            if mu.contained_in(lambda) && valid_skew(lambda, &mu) {
                result.push(mu);
            }
        }
    }
    result
}

/// All weakly increasing sequences (f_1, ..., f_k) with 1 ≤ f_i ≤ max_val.
fn weakly_increasing_sequences(k: usize, max_val: u32) -> Vec<Vec<u32>> {
    let mut result = Vec::new();
    let mut current = Vec::with_capacity(k);
    wi_recurse(k, max_val, 1, &mut current, &mut result);
    result
}

fn wi_recurse(k: usize, max_val: u32, min_val: u32, current: &mut Vec<u32>, result: &mut Vec<Vec<u32>>) {
    if current.len() == k {
        result.push(current.clone());
        return;
    }
    for v in min_val..=max_val {
        current.push(v);
        wi_recurse(k, max_val, v, current, result);
        current.pop();
    }
}

/// Batch-query: which (weight, upper_flags) pairs are already in the DB for a given (lambda, mu)?
fn already_computed(conn: &mut PooledConn, lambda_str: &str, mu_str: &str) -> HashSet<String> {
    let rows: Vec<(String, String)> = conn.exec(
        "SELECT weight, upper_flags FROM gt_ehrhart WHERE lambda = ? AND mu = ?",
        (lambda_str, mu_str),
    ).unwrap_or_default();
    rows.into_iter().map(|(w, uf)| format!("{}|{}", w, uf)).collect()
}

/// Batch-compute Ehrhart data and insert into MariaDB (crash-resumable).
pub fn run(args: PopulateArgs) {
    // Connect to DB (skip in dry-run mode).
    let pool = if !args.dry_run {
        let opts = Opts::from_url(&args.db_url).unwrap_or_else(|e| {
            eprintln!("Invalid database URL: {}", e);
            std::process::exit(1);
        });
        let p = Pool::new(opts).unwrap_or_else(|e| {
            eprintln!("Failed to connect to database: {}", e);
            std::process::exit(1);
        });
        let mut conn = p.get_conn().expect("Failed to get connection");
        conn.query_drop(CREATE_TABLE).expect("Failed to create table");
        eprintln!("Connected to database.");
        Some(p)
    } else {
        None
    };

    eprintln!("  max_n = {}, min_n = {}, skew = {}, flag = {}",
        args.max_n, args.min_n, args.skew, args.flag);
    if args.dry_run {
        eprintln!("  DRY RUN: enumerating work items only.");
    }
    eprintln!();

    let mut total_computed = 0u64;
    let mut total_skipped = 0u64;
    let mut total_errors = 0u64;

    for n in args.min_n..=args.max_n {
        let partitions = Partition::all_of_size(n);
        eprintln!("=== |λ| = {}: {} partitions ===", n, partitions.len());

        for lam in &partitions {
            let lambda_str = parts_str(lam);
            let inners: Vec<Partition> = if args.skew {
                inner_shapes(lam)
            } else {
                vec![Partition::empty()]
            };

            for mu in &inners {
                let mu_str = if mu.num_parts() == 0 { String::new() } else { parts_str(mu) };
                let skew_size = lam.size().saturating_sub(mu.size());
                if skew_size == 0 {
                    continue;
                }

                // Check containment for skew.
                if mu.num_parts() > 0 && !mu.contained_in(lam) {
                    continue;
                }

                // Partition weights of skew_size, filtered by dominance.
                let weights: Vec<Vec<u32>> = Partition::all_of_size(skew_size)
                    .into_iter()
                    .map(|p| p.parts().to_vec())
                    .filter(|w| dominates(lam, w))
                    .collect();

                // Batch resume: fetch already-computed keys for this (lambda, mu).
                let done = if let Some(ref p) = pool {
                    let mut conn = p.get_conn().expect("DB connection failed");
                    already_computed(&mut conn, &lambda_str, &mu_str)
                } else {
                    HashSet::new()
                };

                for w in &weights {
                    let w_str = vec_str(w);

                    // Flag vectors.
                    let flag_vecs: Vec<Option<Vec<u32>>> = if args.flag {
                        let n_rows = lam.num_parts();
                        let k = w.len();
                        if k > 0 && n_rows > 0 {
                            weakly_increasing_sequences(k, n_rows as u32)
                                .into_iter()
                                .map(Some)
                                .collect()
                        } else {
                            vec![None]
                        }
                    } else {
                        vec![None]
                    };

                    for uf in &flag_vecs {
                        let uf_str = uf.as_ref().map(|v| vec_str(v)).unwrap_or_default();
                        let key = format!("{}|{}", w_str, uf_str);

                        if done.contains(&key) {
                            total_skipped += 1;
                            continue;
                        }

                        if args.dry_run {
                            eprintln!("  [dry] λ=({}) μ=({}) w=({}) uf=({})",
                                lambda_str, mu_str, w_str, uf_str);
                            total_computed += 1;
                            continue;
                        }

                        // Check dimension before expensive computation.
                        if args.max_dim.is_some() || args.min_dim.is_some() {
                            let dim = gt_polytope_dim_full(lam.parts(), mu.parts(), w, uf.as_deref(), None);
                            if let Some(max) = args.max_dim {
                                if dim.map_or(false, |d| d > max) {
                                    total_skipped += 1;
                                    continue;
                                }
                            }
                            if let Some(min) = args.min_dim {
                                if dim.map_or(true, |d| d < min) {
                                    total_skipped += 1;
                                    continue;
                                }
                            }
                        }

                        // Compute.
                        let t0 = Instant::now();
                        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
                            compute_one(lam, mu, w, uf.as_deref(), args.max_states)
                        }));
                        let elapsed_ms = t0.elapsed().as_millis() as u64;

                        match result {
                            Ok(row) => {
                                // KKT check.
                                if args.kkt_counterexample {
                                    if let Some(ref poly) = row.polynomial {
                                        if poly.has_negative_coefficient() {
                                            eprintln!("!!! KKT COUNTEREXAMPLE FOUND !!!");
                                            eprintln!("  lambda = ({})", lambda_str);
                                            eprintln!("  mu     = ({})", mu_str);
                                            eprintln!("  weight = ({})", w_str);
                                            eprintln!("  flags  = ({})", uf_str);
                                            eprintln!("  degree = {:?}", row.dimension);
                                            eprintln!("  poly   = {}", poly.display());
                                            if let Some(ref h) = row.hstar {
                                                eprintln!("  h*     = {}", hstar_json(h));
                                            }
                                            std::process::exit(2);
                                        }
                                    }
                                }

                                // Insert.
                                let poly_str = row.polynomial.as_ref().map(|p| p.display());
                                let hstar_str = row.hstar.as_ref().map(|h| hstar_json(h));

                                if let Some(ref p) = pool {
                                    let mut conn = p.get_conn().expect("DB connection failed");
                                    conn.exec_drop(
                                        r"INSERT IGNORE INTO gt_ehrhart
                                          (n, lambda, mu, weight, upper_flags, lower_flags,
                                           kostka, dimension, polynomial, hstar, elapsed_ms)
                                          VALUES (?, ?, ?, ?, ?, '', ?, ?, ?, ?, ?)",
                                        (
                                            n,
                                            &lambda_str,
                                            &mu_str,
                                            &w_str,
                                            &uf_str,
                                            &row.kostka_str,
                                            row.dimension,
                                            &poly_str,
                                            &hstar_str,
                                            elapsed_ms,
                                        ),
                                    ).unwrap_or_else(|e| {
                                        eprintln!("  DB insert error: {}", e);
                                    });
                                }

                                total_computed += 1;
                            }
                            Err(_) => {
                                eprintln!("  SKIP λ=({}) w=({}) uf=({}): DP state limit exceeded",
                                    lambda_str, w_str, uf_str);
                                total_errors += 1;
                            }
                        }
                    }
                }

                let shape = if mu.num_parts() == 0 {
                    format!("({})", lambda_str)
                } else {
                    format!("({})/({})", lambda_str, mu_str)
                };
                let n_weights = weights.len();
                eprintln!("  {}: {} weights, {} computed, {} skipped",
                    shape, n_weights, total_computed, total_skipped);
            }
        }
    }

    eprintln!();
    eprintln!("Done: {} computed, {} skipped (already in DB), {} errors.",
        total_computed, total_skipped, total_errors);
}

// ── single-row computation ───────────────────────────────────────────────────

struct ComputedRow {
    kostka_str: String,
    dimension: Option<i32>,
    polynomial: Option<crate::ehrhart::EhrhartPoly>,
    hstar: Option<Vec<BigInt>>,
}

fn compute_one(
    lam: &Partition,
    mu: &Partition,
    w: &[u32],
    upper_flags: Option<&[u32]>,
    max_states: Option<usize>,
) -> ComputedRow {
    let has_flags = upper_flags.is_some();

    // Kostka value.
    let kostka_val = if has_flags {
        flagged_skew_kostka(lam, mu, w, upper_flags, None, max_states)
    } else {
        skew_kostka(lam, mu, w, max_states, true)
    };
    let kostka_str = kostka_val.to_string();

    // If K = 0, skip Ehrhart computation.
    if kostka_val.is_zero() {
        return ComputedRow {
            kostka_str,
            dimension: None,
            polynomial: None,
            hstar: None,
        };
    }

    // Dimension.
    let dim = gt_polytope_dim_full(lam.parts(), mu.parts(), w, upper_flags, None);

    // Ehrhart polynomial.
    let use_reciprocity = !has_flags;
    let poly = compute_ehrhart(lam, mu, w, upper_flags, None, false, max_states, use_reciprocity);
    let hstar = compute_hstar(&poly);

    ComputedRow {
        kostka_str,
        dimension: dim.map(|d| d as i32),
        polynomial: Some(poly),
        hstar: Some(hstar),
    }
}
