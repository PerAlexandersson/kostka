use clap::Parser;
use std::collections::BTreeMap;
use std::time::Instant;

use combinatoric_core::Partition;
use kostka::kostka_dp::skew_kostka_stats;
use kostka::partition::{parse_partition, parse_weight};

#[derive(Parser, Debug)]
#[command(about = "Explore how the order of weight entries affects Kostka DP performance")]
struct Args {
    #[arg(long)]
    lambda: Option<String>,

    #[arg(long, short = 'w')]
    weight: Option<String>,

    #[arg(long, default_value_t = 5)]
    top: usize,
}

#[derive(Debug, Clone)]
struct RunRecord {
    order: Vec<u32>,
    peak_states: usize,
    micros: u128,
}

fn unique_permutations(weight: &[u32]) -> Vec<Vec<u32>> {
    let mut counts = BTreeMap::<u32, usize>::new();
    for &x in weight {
        *counts.entry(x).or_insert(0) += 1;
    }

    let mut out = Vec::new();
    let mut current = Vec::with_capacity(weight.len());
    permute_rec(weight.len(), &mut counts, &mut current, &mut out);
    out
}

fn permute_rec(
    target_len: usize,
    counts: &mut BTreeMap<u32, usize>,
    current: &mut Vec<u32>,
    out: &mut Vec<Vec<u32>>,
) {
    if current.len() == target_len {
        out.push(current.clone());
        return;
    }

    let keys: Vec<u32> = counts.keys().copied().collect();
    for key in keys {
        let Some(count) = counts.get_mut(&key) else {
            continue;
        };
        if *count == 0 {
            continue;
        }
        *count -= 1;
        current.push(key);
        permute_rec(target_len, counts, current, out);
        current.pop();
        *counts.get_mut(&key).unwrap() += 1;
    }
}

fn descending_order(weight: &[u32]) -> Vec<u32> {
    let mut v = weight.to_vec();
    v.sort_unstable_by(|a, b| b.cmp(a));
    v
}

fn ascending_order(weight: &[u32]) -> Vec<u32> {
    let mut v = weight.to_vec();
    v.sort_unstable();
    v
}

fn middle_heavy_order(weight: &[u32]) -> Vec<u32> {
    let sorted = descending_order(weight);
    let n = sorted.len();
    let mut positions = Vec::with_capacity(n);
    let center = (n.saturating_sub(1)) / 2;
    positions.push(center);
    for step in 1..n {
        if center >= step {
            positions.push(center - step);
            if positions.len() == n {
                break;
            }
        }
        let right = center + step;
        if right < n {
            positions.push(right);
            if positions.len() == n {
                break;
            }
        }
    }

    let mut out = vec![0u32; n];
    for (value, pos) in sorted.into_iter().zip(positions.into_iter()) {
        out[pos] = value;
    }
    out
}

fn find_rank<F>(records: &[RunRecord], target: &[u32], key: F) -> Option<usize>
where
    F: Fn(&RunRecord) -> (usize, u128),
{
    let mut ranked = records.to_vec();
    ranked.sort_by_key(|r| key(r));
    ranked
        .iter()
        .position(|r| r.order == target)
        .map(|idx| idx + 1)
}

fn analyze_case(lambda: &Partition, weight: &[u32], top: usize) {
    let perms = unique_permutations(weight);
    let mut records = Vec::with_capacity(perms.len());
    let mut value_ref = None;

    for perm in perms {
        let start = Instant::now();
        let stats = skew_kostka_stats(lambda, &Partition::empty(), &perm, None, false);
        let micros = start.elapsed().as_micros();
        if let Some(ref value) = value_ref {
            assert_eq!(
                value, &stats.value,
                "Kostka value changed across permutations"
            );
        } else {
            value_ref = Some(stats.value.clone());
        }
        records.push(RunRecord {
            order: perm,
            peak_states: stats.peak_states,
            micros,
        });
    }

    let mut by_peak = records.clone();
    by_peak.sort_by_key(|r| (r.peak_states, r.micros));
    let mut by_time = records.clone();
    by_time.sort_by_key(|r| (r.micros, r.peak_states));

    let descending = descending_order(weight);
    let ascending = ascending_order(weight);
    let middle = middle_heavy_order(weight);

    println!();
    println!("lambda = ({})", lambda);
    println!(
        "weight multiset = ({})  |  unique permutations = {}  |  K = {}",
        weight
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(","),
        records.len(),
        value_ref.unwrap()
    );

    for (label, order) in [
        ("descending", descending.as_slice()),
        ("ascending", ascending.as_slice()),
        ("middle-heavy", middle.as_slice()),
    ] {
        let rec = records.iter().find(|r| r.order == order).unwrap();
        let peak_rank = find_rank(&records, order, |r| (r.peak_states, r.micros)).unwrap();
        let time_rank = find_rank(&records, order, |r| {
            (r.micros as usize, r.peak_states as u128)
        })
        .unwrap();
        println!(
            "  {:<12} order = ({:<15})  peak = {:<5} rank = {:<3}  time = {:<6}us rank = {}",
            label,
            order
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(","),
            rec.peak_states,
            peak_rank,
            rec.micros,
            time_rank
        );
    }

    println!("  best peak-state orders:");
    for rec in by_peak.iter().take(top) {
        println!(
            "    ({:<15})  peak = {:<5}  time = {}us",
            rec.order
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(","),
            rec.peak_states,
            rec.micros
        );
    }

    let largest = *weight.iter().max().unwrap();
    let mut peak_by_pos = vec![(0usize, 0u128, 0usize); weight.len()];
    for rec in &records {
        let pos = rec.order.iter().position(|&x| x == largest).unwrap();
        peak_by_pos[pos].0 += rec.peak_states;
        peak_by_pos[pos].1 += rec.micros;
        peak_by_pos[pos].2 += 1;
    }
    println!("  averages by position of the largest entry {}:", largest);
    for (idx, (peak_sum, micros_sum, count)) in peak_by_pos.into_iter().enumerate() {
        if count == 0 {
            continue;
        }
        println!(
            "    pos {:>2}: avg peak = {:>8.2}  avg time = {:>8.2}us",
            idx + 1,
            peak_sum as f64 / count as f64,
            micros_sum as f64 / count as f64
        );
    }
}

fn main() {
    let args = Args::parse();

    if let (Some(lambda), Some(weight)) = (args.lambda.as_deref(), args.weight.as_deref()) {
        let lambda = parse_partition(lambda).expect("invalid lambda");
        let weight = parse_weight(weight).expect("invalid weight");
        analyze_case(&lambda, &weight, args.top);
        return;
    }

    let presets = [
        ("7,5,3", "5,4,3,2,1"),
        ("5,4,3,2,1", "5,4,3,2,1"),
        ("8,4,2,1", "5,4,3,2,1"),
    ];

    for (lambda, weight) in presets {
        let lambda = parse_partition(lambda).expect("invalid preset lambda");
        let weight = parse_weight(weight).expect("invalid preset weight");
        analyze_case(&lambda, &weight, args.top);
    }
}
