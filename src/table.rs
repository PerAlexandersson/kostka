/// Batch / table mode.
use clap::Args;
use rayon::prelude::*;
use num_bigint::BigInt;

use crate::partition::{parse_partition, Partition};
use crate::kostka_dp::{kostka, skew_kostka};
use crate::ehrhart::{compute_ehrhart, compute_hstar, is_palindromic, is_unimodal};

#[derive(Args)]
pub struct TableArgs {
    /// Outer shape lambda (required unless --n is given)
    #[arg(long)]
    pub lambda: Option<String>,

    /// Inner shape mu for skew (optional)
    #[arg(long, default_value = "")]
    pub mu: String,

    /// Enumerate all weights of size |lambda|-|mu|
    #[arg(long)]
    pub all_weights: bool,

    /// Treat weight as partition (deduplicate orderings)
    #[arg(long)]
    pub weight_as_partition: bool,

    /// Compute full Kostka matrix for all lambda, w of size n
    #[arg(long)]
    pub n: Option<u32>,

    /// Max alphabet size for --all-weights (default: |lambda|-|mu|)
    #[arg(long)]
    pub alphabet: Option<usize>,

    /// Also compute and display Ehrhart polynomial
    #[arg(long)]
    pub ehrhart: bool,

    /// Also compute and display h*-vector
    #[arg(long)]
    pub hstar: bool,

    /// Only show rows with degree >= d
    #[arg(long)]
    pub min_degree: Option<usize>,

    /// Only show rows with degree <= d
    #[arg(long)]
    pub max_degree: Option<usize>,

    /// Only show palindromic h*-vectors
    #[arg(long)]
    pub palindromic: bool,

    /// Only show unimodal h*-vectors
    #[arg(long)]
    pub unimodal: bool,

    /// Maximum number of rows to display
    #[arg(long)]
    pub limit: Option<usize>,

    /// Abort if any DP level exceeds this many states (prevents OOM on large inputs)
    #[arg(long)]
    pub max_states: Option<usize>,

    /// Output format: text, json, csv
    #[arg(long, default_value = "text")]
    pub format: String,
}

pub fn run(args: TableArgs) {
    if let Some(n) = args.n {
        run_full_matrix(n, &args);
    } else if args.all_weights {
        run_all_weights(args);
    } else {
        eprintln!("table: specify --n for the full Kostka matrix, or --lambda + --all-weights.");
        std::process::exit(1);
    }
}

/// Full Kostka matrix: all lambda, w ⊢ n, in dominance order.
fn run_full_matrix(n: u32, args: &TableArgs) {
    let shapes = Partition::all_of_size(n);
    let weights = Partition::all_of_size(n);

    // Header
    if args.format == "text" {
        let header_shapes: Vec<String> = weights.iter().map(|p| format!("({})", p)).collect();
        println!("Kostka matrix for n={}\n", n);
        print!("{:>15}  ", "");
        for h in &header_shapes {
            print!("{:>10}", h);
        }
        println!();
    }

    for lam in &shapes {
        let max_states = args.max_states;
        let row: Vec<_> = weights.par_iter().map(|w| kostka(lam, w.parts(), max_states)).collect();
        match args.format.as_str() {
            "csv" => {
                let vals: Vec<String> = row.iter().map(|x| x.to_string()).collect();
                println!("{},{}", lam, vals.join(","));
            }
            "json" => {
                let vals: Vec<String> = row.iter().map(|x| x.to_string()).collect();
                println!("{{\"lambda\":[{}],\"values\":[{}]}}",
                    lam.parts().iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","),
                    vals.join(","));
            }
            _ => {
                print!("{:>15}  ", format!("({})", lam));
                for v in &row {
                    print!("{:>10}", v);
                }
                println!();
            }
        }
    }
}

/// All weights for a fixed shape (and optional skew), with optional Ehrhart/hstar.
fn run_all_weights(args: TableArgs) {
    let lambda_str = args.lambda.as_deref().unwrap_or_else(|| { eprintln!("--lambda required"); std::process::exit(1); });
    let lam = parse_partition(lambda_str).unwrap_or_else(|e| { eprintln!("{}", e); std::process::exit(1); });
    let inner = if args.mu.is_empty() { Partition::empty() } else { parse_partition(&args.mu).unwrap_or_else(|e| { eprintln!("{}", e); std::process::exit(1); }) };

    let skew_size = lam.size().saturating_sub(inner.size());
    let alphabet = args.alphabet.unwrap_or(skew_size as usize);

    // Generate all compositions (or partitions) of skew_size.
    let weights: Vec<Vec<u32>> = if args.weight_as_partition {
        Partition::all_of_size(skew_size)
            .into_iter()
            .map(|p| p.parts().to_vec())
            .collect()
    } else {
        // All compositions into at most `alphabet` positive parts.
        compositions_up_to(skew_size, alphabet)
    };

    let shape_str = if inner.num_parts() == 0 {
        format!("{}", lam)
    } else {
        format!("{}/{}", lam, inner)
    };

    if args.format == "text" {
        if args.ehrhart {
            println!("Ehrhart polynomials for GT( {} | w )  (|w| = {}):\n", shape_str, skew_size);
        } else {
            println!("K( {} | w ) for all weights w of size {}:\n", shape_str, skew_size);
        }
    }

    // Compute all rows in parallel.
    let max_states = args.max_states;
    let mut rows: Vec<TableRow> = weights.par_iter().map(|w| {
        let k = skew_kostka(&lam, &inner, w, max_states, true);
        let (poly, hv) = if args.ehrhart || args.hstar || args.min_degree.is_some() || args.max_degree.is_some() || args.palindromic || args.unimodal {
            let p = compute_ehrhart(&lam, &inner, w, None, None, false, max_states, true);
            let h = if args.hstar || args.palindromic || args.unimodal { Some(compute_hstar(&p)) } else { None };
            (Some(p), h)
        } else {
            (None, None)
        };
        TableRow { weight: w.clone(), kostka_val: k, poly, hstar: hv }
    }).collect();

    // Filter.
    rows.retain(|row| {
        if let Some(p) = &row.poly {
            if let Some(min_d) = args.min_degree { if p.degree < min_d { return false; } }
            if let Some(max_d) = args.max_degree { if p.degree > max_d { return false; } }
        }
        if let Some(h) = &row.hstar {
            if args.palindromic && !is_palindromic(h) { return false; }
            if args.unimodal && !is_unimodal(h) { return false; }
        }
        true
    });

    // Apply limit.
    if let Some(lim) = args.limit {
        rows.truncate(lim);
    }

    // Print.
    for row in &rows {
        print_row(row, &args);
    }
}

struct TableRow {
    weight: Vec<u32>,
    kostka_val: num_bigint::BigUint,
    poly: Option<crate::ehrhart::EhrhartPoly>,
    hstar: Option<Vec<BigInt>>,
}

fn print_row(row: &TableRow, args: &TableArgs) {
    let w_str = row.weight.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");

    match args.format.as_str() {
        "csv" => {
            let poly_str = row.poly.as_ref().map(|p| p.display()).unwrap_or_default();
            let deg_str = row.poly.as_ref().map(|p| p.degree.to_string()).unwrap_or_default();
            let hstar_str = row.hstar.as_ref().map(|h| {
                format!("\"[{}]\"", h.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","))
            }).unwrap_or_default();
            println!("{},{},{},{},\"{}\"", w_str, row.kostka_val, deg_str, poly_str, hstar_str);
        }
        "json" => {
            let poly_str = row.poly.as_ref().map(|p| format!("\"{}\"", p.display())).unwrap_or("null".into());
            let deg_str = row.poly.as_ref().map(|p| p.degree.to_string()).unwrap_or("null".into());
            let hstar_str = row.hstar.as_ref().map(|h| {
                format!("[{}]", h.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","))
            }).unwrap_or("null".into());
            println!("{{\"weight\":[{}],\"value\":\"{}\",\"degree\":{},\"polynomial\":{},\"hstar\":{}}}",
                w_str, row.kostka_val, deg_str, poly_str, hstar_str);
        }
        _ => {
            if args.ehrhart {
                let deg = row.poly.as_ref().map(|p| p.degree).unwrap_or(0);
                let poly_s = row.poly.as_ref().map(|p| p.display()).unwrap_or_default();
                if args.hstar {
                    let h_s = row.hstar.as_ref().map(|h| {
                        format!("[{}]", h.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "))
                    }).unwrap_or_default();
                    println!("  w = ({:20})  deg={}  poly= {:30}  h*={}", w_str, deg, poly_s, h_s);
                } else {
                    println!("  w = ({:20})  deg={}  poly= {}", w_str, deg, poly_s);
                }
            } else {
                println!("  w = ({:20})  ->  {}", w_str, row.kostka_val);
            }
        }
    }
}

/// All compositions of n into positive parts, with at most `max_parts` parts.
fn compositions_up_to(n: u32, max_parts: usize) -> Vec<Vec<u32>> {
    let mut result = Vec::new();
    compose(n, max_parts, &mut vec![], &mut result);
    result
}

fn compose(remaining: u32, parts_left: usize, current: &mut Vec<u32>, result: &mut Vec<Vec<u32>>) {
    if remaining == 0 {
        if !current.is_empty() {
            result.push(current.clone());
        }
        return;
    }
    if parts_left == 0 {
        return;
    }
    for p in 1..=remaining {
        current.push(p);
        compose(remaining - p, parts_left - 1, current, result);
        current.pop();
    }
}
