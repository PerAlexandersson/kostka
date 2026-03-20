//! CLI for computing Kostka coefficients, Ehrhart polynomials, and h*-vectors
//! of Gelfand-Tsetlin polytopes.  See `kostka --help` for usage.

mod partition;
mod kostka_dp;
mod gt_dim;
mod ehrhart;
mod syt;
mod table;
mod populate;
mod lr;

use clap::{Parser, Subcommand};
use partition::{parse_partition, parse_weight, Partition};
use kostka_dp::{kostka, skew_kostka, strict_kostka, strict_skew_kostka, flagged_skew_kostka};
use gt_dim::gt_polytope_dim_full;
use ehrhart::{compute_ehrhart, compute_hstar, is_palindromic, is_unimodal, verify_reciprocity};
use syt::{count_syt, hook_lengths};

// ── top-level CLI ───────────────────────────────────────────────────────────────

#[derive(Parser)]
#[command(
    name = "kostka",
    about = "Kostka coefficients, Ehrhart polynomials, and h*-vectors for GT polytopes",
    long_about = "\
Compute Kostka coefficients K(λ/μ, w) and related invariants of Gelfand-Tsetlin (GT) polytopes.

A Kostka coefficient K(λ, w) counts semistandard Young tableaux (SSYT) of shape λ
and content (weight) w.  The skew variant K(λ/μ, w) counts skew SSYT where entries
from the inner shape μ are removed.  Both equal the number of lattice points of the
GT polytope GT(λ/μ, w).

As the inputs are scaled uniformly (λ→n·λ, μ→n·μ, w→n·w), the count
    n ↦ K(n·λ / n·μ, n·w)
is a polynomial in n — the Ehrhart polynomial of GT(λ/μ, w).  Its h*-vector encodes
the distribution of lattice points across dilations via the Hilbert series.

Partitions are given as comma-separated weakly decreasing integers, e.g. '3,2,1'.
Weights are comma-separated non-negative integers summing to |λ|−|μ|, e.g. '2,2,2'.",
    after_help = "\
EXAMPLES:
  Single Kostka coefficient K(3,2,1 | 2,2,2):
    kostka kostka --lambda 3,2,1 -w 2,2,2

  Skew Kostka coefficient K(4,3,2 / 2,1 | 2,2,2):
    kostka skew --lambda 4,3,2 --mu 2,1 -w 2,2,2

  Ehrhart polynomial of GT(3,2,1 | 1,1,1,1,1,1):
    kostka ehrhart --lambda 3,2,1 -w 1,1,1,1,1,1

  h*-vector and palindromicity/unimodality:
    kostka hstar --lambda 3,2,1 -w 1,1,1,1,1,1

  Number of SYT of shape (5,3,2) with hook lengths:
    kostka syt --lambda 5,3,2 --hooks

  Full Kostka matrix for all partitions of 4:
    kostka table --n 4

  All Kostka coefficients for a fixed shape:
    kostka table --lambda 3,2,1 --all-weights"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

// ── subcommands ─────────────────────────────────────────────────────────────────

#[derive(Subcommand)]
enum Command {
    /// Compute K(λ, w): the Kostka coefficient for a straight shape
    ///
    /// K(λ, w) counts semistandard Young tableaux of shape λ and content w,
    /// equivalently the number of lattice points of the GT polytope GT(λ, w).
    /// The weight w must sum to |λ|; its order matters (different orderings
    /// generally give different values unless w is sorted).
    #[command(after_help = "\
EXAMPLES:
  K(3,2,1 | 3,2,1) — the diagonal entry:
    kostka kostka --lambda 3,2,1 -w 3,2,1

  K(4,2 | 2,2,1,1) in JSON format:
    kostka kostka --lambda 4,2 -w 2,2,1,1 --format json")]
    Kostka {
        /// Shape λ as comma-separated weakly decreasing parts, e.g. '3,2,1'
        #[arg(long)]
        lambda: String,

        /// Content vector w summing to |λ|, e.g. '2,2,2'
        #[arg(long, short = 'w')]
        weight: String,

        /// Output format: text (default), json, csv
        #[arg(long, default_value = "text")]
        format: String,

        /// Abort DP if any level exceeds N states (prevents OOM on huge inputs)
        #[arg(long)]
        max_states: Option<usize>,
    },

    /// Compute K(λ/μ, w): skew Kostka coefficient, with optional row flags
    ///
    /// Counts skew SSYT of shape λ/μ and content w, i.e. lattice points of
    /// GT(λ/μ, w).  Row flags restrict individual GT-pattern rows:
    ///   upper-flags u_i means row i is bounded above entry-wise by u_i,
    ///   lower-flags l_i means row i is bounded below entry-wise by l_i.
    /// Flags must have length equal to the number of rows of λ.
    /// With flags, Ehrhart reciprocity is not available.
    #[command(after_help = "\
EXAMPLES:
  K(4,3,2 / 2,1 | 2,2,2):
    kostka skew --lambda 4,3,2 --mu 2,1 -w 2,2,2

  Flagged skew Kostka with upper row bounds 4,4,4:
    kostka skew --lambda 4,3,2 --mu 2,1 -w 2,2,2 --upper-flags 4,4,4")]
    Skew {
        /// Outer shape λ
        #[arg(long)]
        lambda: String,

        /// Inner shape μ (must fit inside λ)
        #[arg(long)]
        mu: String,

        /// Content vector w summing to |λ|−|μ|
        #[arg(long, short = 'w')]
        weight: String,

        /// Upper row-flag bounds (comma-separated, length = rows of λ)
        #[arg(long)]
        upper_flags: Option<String>,

        /// Lower row-flag bounds (comma-separated, length = rows of λ)
        #[arg(long)]
        lower_flags: Option<String>,

        /// Output format: text (default), json, csv
        #[arg(long, default_value = "text")]
        format: String,

        /// Abort DP if any level exceeds N states
        #[arg(long)]
        max_states: Option<usize>,
    },

    /// Dimension of the GT polytope GT(λ/μ, w) (= degree of the Ehrhart polynomial)
    ///
    /// The GT polytope lives in a real vector space whose dimension equals the
    /// number of free (non-fixed) entries in a GT pattern for (λ/μ, w).
    /// Dimension 0 means the polytope is a single point: K(λ/μ, w) ≤ 1.
    #[command(after_help = "\
EXAMPLES:
  Dimension of GT(3,2,1 | 1,1,1,1,1,1):
    kostka degree --lambda 3,2,1 -w 1,1,1,1,1,1

  With an inner shape:
    kostka degree --lambda 4,3 --mu 2,1 -w 2,2")]
    Degree {
        /// Outer shape λ
        #[arg(long)]
        lambda: String,

        /// Inner shape μ (empty = straight shape)
        #[arg(long, default_value = "")]
        mu: String,

        /// Content vector w
        #[arg(long, short = 'w')]
        weight: String,

        /// Upper row-flag bounds
        #[arg(long)]
        upper_flags: Option<String>,

        /// Lower row-flag bounds
        #[arg(long)]
        lower_flags: Option<String>,

        #[arg(long)]
        verbose: bool,
    },

    /// Ehrhart polynomial P(n) = K(n·λ / n·μ, n·w) of the GT polytope
    ///
    /// By a theorem of Rassart (2004), the function n ↦ K(n·λ/n·μ, n·w) is
    /// a polynomial in n of degree = dim GT(λ/μ, w).  This command computes
    /// that polynomial by interpolation.
    ///
    /// By default, Ehrhart-Macdonald reciprocity is used to reduce the maximum
    /// dilation needed: P(0) = 1 is free, and P(−t) = (−1)^d · K_strict(t·λ/t·μ, t·w)
    /// where K_strict counts strict (interior) GT patterns.  The adaptive strategy
    /// greedily evaluates whichever side (positive or negative t) last returned the
    /// smaller count, collecting zeros on the negative side for free.
    ///
    /// Flags disable reciprocity and fall back to plain positive-dilation interpolation.
    #[command(after_help = "\
EXAMPLES:
  Ehrhart polynomial of GT(3,2,1 | 1,1,1,1,1,1):
    kostka ehrhart --lambda 3,2,1 -w 1,1,1,1,1,1

  Show values up to n=10 in JSON:
    kostka ehrhart --lambda 3,2,1 -w 1,1,1,1,1,1 --show-values 10 --format json

  Verify reciprocity for the first 5 dilations:
    kostka ehrhart --lambda 3,2,1 -w 1,1,1,1,1,1 --verify-reciprocity 5

  Skew shape (4,3/2,1) with weight 2,2,2:
    kostka ehrhart --lambda 4,3 --mu 2,1 -w 2,2,2")]
    Ehrhart {
        /// Outer shape λ
        #[arg(long)]
        lambda: String,

        /// Inner shape μ (empty = straight shape)
        #[arg(long, default_value = "")]
        mu: String,

        /// Content vector w
        #[arg(long, short = 'w')]
        weight: String,

        /// Upper row-flag bounds (disables reciprocity)
        #[arg(long)]
        upper_flags: Option<String>,

        /// Lower row-flag bounds (disables reciprocity)
        #[arg(long)]
        lower_flags: Option<String>,

        /// Print P(1), P(2), ..., P(N) below the polynomial
        #[arg(long, default_value = "5")]
        show_values: u64,

        /// Output format: text (default), json, csv
        #[arg(long, default_value = "text")]
        format: String,

        #[arg(long)]
        verbose: bool,

        /// Abort DP if any level exceeds N states
        #[arg(long)]
        max_states: Option<usize>,

        /// Verify Ehrhart-Macdonald reciprocity (−1)^d · K_strict(t·λ, t·w) = P(−t) for t=1..N
        #[arg(long, default_value = "0")]
        verify_reciprocity: u64,

        /// Use plain positive-dilation interpolation n=1,...,d+1 instead of reciprocity
        #[arg(long)]
        no_reciprocity: bool,
    },

    /// Interior lattice points of GT(λ/μ, w): strict GT pattern count
    ///
    /// Counts GT patterns with strictly increasing rows and columns (strict
    /// inequalities throughout).  This equals the number of lattice points in
    /// the relative interior of GT(λ/μ, w), denoted K_strict(λ/μ, w).
    ///
    /// By Ehrhart-Macdonald reciprocity:
    ///   K_strict(t·λ, t·μ, t·w) = (−1)^d · P(−t)
    /// where P is the Ehrhart polynomial and d = dim GT.
    #[command(after_help = "\
EXAMPLES:
  Interior points of GT(3,2,1 | 1,1,1,1,1,1):
    kostka strict-kostka --lambda 3,2,1 -w 1,1,1,1,1,1

  Skew interior count:
    kostka strict-kostka --lambda 4,3 --mu 2,1 -w 2,2,2")]
    StrictKostka {
        /// Outer shape λ
        #[arg(long)]
        lambda: String,

        /// Inner shape μ (empty = straight shape)
        #[arg(long, default_value = "")]
        mu: String,

        /// Content vector w
        #[arg(long, short = 'w')]
        weight: String,

        /// Abort DP if any level exceeds N states
        #[arg(long)]
        max_states: Option<usize>,
    },

    /// h*-vector of the GT polytope GT(λ/μ, w)
    ///
    /// The h*-vector (h*_0, ..., h*_d) encodes the Hilbert series of the
    /// Ehrhart ring:
    ///   Σ_{n≥0} P(n) t^n  =  (h*_0 + h*_1 t + ... + h*_d t^d) / (1−t)^{d+1}
    ///
    /// All h*_k are non-negative integers.  Palindromicity (h*_k = h*_{d−k})
    /// characterises Gorenstein polytopes; unimodality is expected but not
    /// always guaranteed.
    ///
    /// Computed from the Ehrhart polynomial using the formula:
    ///   h*_k = Σ_{j=0}^k (−1)^{k−j} C(d+1, k−j) P(j)
    #[command(after_help = "\
EXAMPLES:
  h*-vector of GT(3,2,1 | 1,1,1,1,1,1):
    kostka hstar --lambda 3,2,1 -w 1,1,1,1,1,1

  JSON output (useful for scripting):
    kostka hstar --lambda 3,2,1 -w 1,1,1,1,1,1 --format json

  Skew shape:
    kostka hstar --lambda 4,3 --mu 2,1 -w 2,2,2")]
    Hstar {
        /// Outer shape λ
        #[arg(long)]
        lambda: String,

        /// Inner shape μ (empty = straight shape)
        #[arg(long, default_value = "")]
        mu: String,

        /// Content vector w
        #[arg(long, short = 'w')]
        weight: String,

        /// Upper row-flag bounds (disables reciprocity)
        #[arg(long)]
        upper_flags: Option<String>,

        /// Lower row-flag bounds (disables reciprocity)
        #[arg(long)]
        lower_flags: Option<String>,

        /// Output format: text (default), json, csv
        #[arg(long, default_value = "text")]
        format: String,

        #[arg(long)]
        verbose: bool,

        /// Abort DP if any level exceeds N states
        #[arg(long)]
        max_states: Option<usize>,

        /// Use plain positive-dilation interpolation instead of reciprocity
        #[arg(long)]
        no_reciprocity: bool,
    },

    /// Number of standard Young tableaux of shape λ (hook-length formula)
    ///
    /// f^λ = |λ|! / Π_{(i,j)∈λ} hook(i,j)
    ///
    /// where hook(i,j) = (row length − j) + (column height − i) + 1.
    /// This equals K(λ, 1^|λ|): the Kostka coefficient with the all-ones weight.
    /// Computed exactly in O(|λ|) arithmetic operations.
    #[command(after_help = "\
EXAMPLES:
  f^(5,3,2):
    kostka syt --lambda 5,3,2

  f^(5,3,2) with hook lengths displayed:
    kostka syt --lambda 5,3,2 --hooks")]
    Syt {
        /// Shape λ as comma-separated weakly decreasing parts
        #[arg(long)]
        lambda: String,

        /// Also display hook lengths at each cell of the diagram
        #[arg(long)]
        hooks: bool,
    },

    /// Batch computation: all weights for a fixed shape, or full Kostka matrix
    ///
    /// With --lambda and --all-weights: prints K(λ, w) for every weight w of
    /// size |λ|−|μ| (optionally restricted to partitions with --weight-as-partition).
    ///
    /// With --n: prints the full Kostka matrix K(λ, w) for all partitions λ and w
    /// of size n, in the standard dominance order.
    ///
    /// Adding --ehrhart or --hstar computes the Ehrhart polynomial or h*-vector
    /// alongside each Kostka coefficient.  Filter flags (--palindromic, --unimodal,
    /// --min-degree, --max-degree) select interesting rows.
    #[command(after_help = "\
EXAMPLES:
  All Kostka coefficients for shape (3,2,1):
    kostka table --lambda 3,2,1 --all-weights

  Restrict to weight-partitions (deduplicate orderings):
    kostka table --lambda 3,2,1 --all-weights --weight-as-partition

  Full 4×4 Kostka matrix:
    kostka table --n 4

  Shapes with palindromic h* and degree ≥ 5:
    kostka table --lambda 4,3,2 --all-weights --hstar --palindromic --min-degree 5

  CSV output of all Ehrhart data for a shape:
    kostka table --lambda 3,2,1 --all-weights --ehrhart --format csv")]
    Table(table::TableArgs),

    /// Batch-compute Ehrhart data and store in MariaDB (crash-resumable)
    Populate(populate::PopulateArgs),

    /// Littlewood-Richardson coefficient c^λ_{μ,ν}
    ///
    /// Counts skew SSYT of shape λ/μ with content ν whose reverse reading word
    /// is a lattice word (Yamanouchi word).  Computed by two independent methods
    /// and cross-checked:
    ///
    /// 1. **GT DP**: augmented Gelfand-Tsetlin DP with Yamanouchi constraints
    ///    integrated directly into the horizontal-strip enumeration.
    ///
    /// 2. **Kostka inverse**: uses the identity K(λ/μ, α) = Σ_ν c_ν K(ν, α)
    ///    solved by back-substitution (the Kostka matrix is upper unitriangular
    ///    in dominance order).
    #[command(after_help = "\
EXAMPLES:
  c^(4,2,1)_{(2,1),(3,1,1)}:
    kostka lr --lambda 4,2,1 --mu 2,1 --nu 3,1,1

  JSON output:
    kostka lr --lambda 4,2,1 --mu 2,1 --nu 3,1,1 --format json")]
    Lr {
        /// Outer shape λ
        #[arg(long)]
        lambda: String,

        /// Inner shape μ (empty for straight shape, but c^λ_{∅,ν} = δ_{λ,ν})
        #[arg(long, default_value = "")]
        mu: String,

        /// Content partition ν (the Schur function index in s_{λ/μ} = Σ c·s_ν)
        #[arg(long)]
        nu: String,

        /// Output format: text (default), json
        #[arg(long, default_value = "text")]
        format: String,

        /// Abort DP if any level exceeds N states
        #[arg(long)]
        max_states: Option<usize>,
    },
}

fn main() {
    let cli = Cli::parse();
    match cli.command {
        Command::Kostka { lambda, weight, format, max_states } => {
            let lam = parse_part(&lambda);
            let w = parse_wt(&weight);
            let k = kostka(&lam, &w, max_states);
            match format.as_str() {
                "json" => println!("{{\"lambda\":{},\"weight\":{},\"value\":\"{}\"}}", json_ints(lam.parts()), json_ints(&w), k),
                "csv"  => println!("{},{},{}", lam, weight, k),
                _      => println!("K( {} | {} ) = {}", lam, format_weight(&w), k),
            }
        }
        Command::Skew { lambda, mu, weight, upper_flags, lower_flags, format, max_states } => {
            let lam = parse_part(&lambda);
            let inner = parse_part(&mu);
            let w = parse_wt(&weight);
            let uf = parse_opt_weight(upper_flags.as_deref());
            let lf = parse_opt_weight(lower_flags.as_deref());
            let k = if uf.is_some() || lf.is_some() {
                flagged_skew_kostka(&lam, &inner, &w, uf.as_deref(), lf.as_deref(), max_states)
            } else {
                skew_kostka(&lam, &inner, &w, max_states, true)
            };
            match format.as_str() {
                "json" => println!("{{\"lambda\":{},\"mu\":{},\"weight\":{},\"value\":\"{}\"}}", json_ints(lam.parts()), json_ints(inner.parts()), json_ints(&w), k),
                "csv"  => println!("{},{},{},{}", lam, inner, weight, k),
                _      => println!("K( {}/{} | {} ) = {}", lam, inner, format_weight(&w), k),
            }
        }
        Command::Degree { lambda, mu, weight, upper_flags, lower_flags, verbose: _ } => {
            let lam = parse_part(&lambda);
            let inner = parse_inner(&mu);
            let w = parse_wt(&weight);
            let uf = parse_opt_weight(upper_flags.as_deref());
            let lf = parse_opt_weight(lower_flags.as_deref());
            let deg = gt_polytope_dim_full(lam.parts(), inner.parts(), &w, uf.as_deref(), lf.as_deref());
            match deg {
                None    => println!("degree: empty polytope"),
                Some(d) => println!("degree: {}", d),
            }
        }
        Command::Ehrhart { lambda, mu, weight, upper_flags, lower_flags, show_values, format, verbose: _, max_states, verify_reciprocity: n_checks, no_reciprocity } => {
            let lam = parse_part(&lambda);
            let inner = parse_inner(&mu);
            let w = parse_wt(&weight);
            let uf = parse_opt_weight(upper_flags.as_deref());
            let lf = parse_opt_weight(lower_flags.as_deref());
            let poly = compute_ehrhart(&lam, &inner, &w, uf.as_deref(), lf.as_deref(), false, max_states, !no_reciprocity);
            let shape = skew_shape_str(&lam, &inner);
            match format.as_str() {
                "json" => {
                    let vals: Vec<String> = (1..=show_values).map(|n| poly.eval(n).to_string()).collect();
                    println!("{{\"lambda\":{},\"mu\":{},\"weight\":{},\"degree\":{},\"polynomial\":\"{}\",\"values\":[{}]}}",
                        json_ints(lam.parts()), json_ints(inner.parts()), json_ints(&w),
                        poly.degree, poly.display(), vals.join(","));
                }
                "csv" => println!("{},{},{},{},\"{}\"", lam, inner, weight, poly.degree, poly.display()),
                _ => {
                    println!("Ehrhart polynomial of GT( {} | {} ):", shape, format_weight(&w));
                    println!("  degree:     {}", poly.degree);
                    println!("  polynomial: {}", poly.display_factored());
                    let vals: Vec<String> = (1..=show_values).map(|n| format!("n={} -> {}", n, poly.eval(n))).collect();
                    println!("  values:     {}", vals.join(",  "));
                }
            }
            if n_checks > 0 {
                println!("Verifying Ehrhart-Macdonald reciprocity ({} checks):", n_checks);
                let ok = verify_reciprocity(&poly, &lam, &inner, &w, n_checks, max_states);
                println!("  result: {}", if ok { "ALL PASSED" } else { "FAILURES FOUND" });
            }
        }
        Command::StrictKostka { lambda, mu, weight, max_states } => {
            let lam = parse_part(&lambda);
            let inner = parse_inner(&mu);
            let w = parse_wt(&weight);
            let k = if inner.num_parts() == 0 {
                strict_kostka(&lam, &w, max_states)
            } else {
                strict_skew_kostka(&lam, &inner, &w, max_states, true)
            };
            let shape = skew_shape_str(&lam, &inner);
            println!("K_strict( {} | {} ) = {}", shape, format_weight(&w), k);
        }
        Command::Hstar { lambda, mu, weight, upper_flags, lower_flags, format, verbose: _, max_states, no_reciprocity } => {
            let lam = parse_part(&lambda);
            let inner = parse_inner(&mu);
            let w = parse_wt(&weight);
            let uf = parse_opt_weight(upper_flags.as_deref());
            let lf = parse_opt_weight(lower_flags.as_deref());
            let poly = compute_ehrhart(&lam, &inner, &w, uf.as_deref(), lf.as_deref(), false, max_states, !no_reciprocity);
            let hstar = compute_hstar(&poly);
            let pal = is_palindromic(&hstar);
            let uni = is_unimodal(&hstar);
            let hstar_str = format!("[{}]", hstar.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
            let shape = skew_shape_str(&lam, &inner);
            match format.as_str() {
                "json" => println!("{{\"lambda\":{},\"mu\":{},\"weight\":{},\"degree\":{},\"hstar\":[{}],\"palindromic\":{},\"unimodal\":{}}}",
                    json_ints(lam.parts()), json_ints(inner.parts()), json_ints(&w),
                    poly.degree, hstar.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","), pal, uni),
                "csv" => println!("{},{},{},{},\"{}\"", lam, inner, weight, poly.degree, hstar_str),
                _ => {
                    println!("h*-vector of GT( {} | {} ):", shape, format_weight(&w));
                    println!("  degree:      {}", poly.degree);
                    println!("  h* =         {}", hstar_str);
                    println!("  palindromic: {}", if pal { "yes" } else { "no" });
                    println!("  unimodal:    {}", if uni { "yes" } else { "no" });
                }
            }
        }
        Command::Syt { lambda, hooks } => {
            let lam = parse_part(&lambda);
            let count = count_syt(&lam);
            println!("f^({}) = {}", lam, count);
            if hooks {
                let hs = hook_lengths(&lam);
                let rows: Vec<String> = hs.iter().map(|row| row.iter().map(|h| h.to_string()).collect::<Vec<_>>().join(" ")).collect();
                println!("hook lengths: {}", rows.join(" / "));
            }
        }
        Command::Table(args) => {
            table::run(args);
        }
        Command::Populate(args) => {
            populate::run(args);
        }
        Command::Lr { lambda, mu, nu, format, max_states } => {
            let lam = parse_part(&lambda);
            let inner = parse_inner(&mu);
            let nu_part = parse_part(&nu);
            lr::run(&lam, &inner, &nu_part, &format, max_states);
        }
    }
}

// ── parsing helpers ────────────────────────────────────────────────────────────

fn parse_part(s: &str) -> Partition {
    parse_partition(s).unwrap_or_else(|e| { eprintln!("{}", e); std::process::exit(1); })
}

fn parse_wt(s: &str) -> Vec<u32> {
    parse_weight(s).unwrap_or_else(|e| { eprintln!("{}", e); std::process::exit(1); })
}

/// Parse an optional inner shape: empty string → empty partition.
fn parse_inner(s: &str) -> Partition {
    if s.is_empty() { Partition::empty() } else { parse_part(s) }
}

fn parse_opt_weight(s: Option<&str>) -> Option<Vec<u32>> {
    s.map(|x| parse_wt(x))
}

// ── formatting helpers ─────────────────────────────────────────────────────────

fn format_weight(w: &[u32]) -> String {
    w.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}

fn json_ints(v: &[u32]) -> String {
    format!("[{}]", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","))
}

fn skew_shape_str(lam: &Partition, inner: &Partition) -> String {
    if inner.num_parts() == 0 { format!("{}", lam) } else { format!("{}/{}", lam, inner) }
}
