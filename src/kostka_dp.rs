/// Level-by-level DP for computing skew Kostka coefficients K(lambda/mu, w).
///
/// A SSYT of shape lambda/mu and content w = (w_1,...,w_k) corresponds bijectively to a chain
///   mu = alpha^0 ⊂ alpha^1 ⊂ ... ⊂ alpha^k = lambda
/// where each alpha^i / alpha^{i-1} is a horizontal strip of size w_i.
///
/// DP state at level i: HashMap<Partition, BigUint> — number of paths from mu to each alpha.

use std::collections::HashMap;
use num_bigint::BigUint;
use num_traits::{Zero, One};
use crate::partition::Partition;

/// Enumerate all partitions beta such that beta/alpha is a horizontal strip of size `strip_size`,
/// with alpha ⊆ beta ⊆ lambda (containment constraint).
///
/// A horizontal strip means at most one box per column, i.e.:
///   c[r] = beta[r] - alpha[r] ≥ 0
///   c[r] ≤ lambda[r] - alpha[r]              (stay in lambda)
///   c[r+1] ≤ alpha[r] - alpha[r+1]           (horizontal strip condition)
///   sum(c) = strip_size
///
/// We enumerate valid increment vectors c by a simple bounded recursion.
pub fn horizontal_strip_extensions(
    alpha: &Partition,
    lambda: &Partition,
    strip_size: u32,
) -> Vec<Partition> {
    if strip_size == 0 {
        return vec![alpha.clone()];
    }
    let n_rows = lambda.num_parts();
    if n_rows == 0 {
        return vec![];
    }

    let mut results = Vec::new();
    let mut increments = vec![0u32; n_rows];
    enumerate_strips(alpha, lambda, strip_size, 0, n_rows, &mut increments, &mut results);
    results
}

fn enumerate_strips(
    alpha: &Partition,
    lambda: &Partition,
    remaining: u32,
    row: usize,
    n_rows: usize,
    increments: &mut Vec<u32>,
    results: &mut Vec<Partition>,
) {
    if row == n_rows {
        if remaining == 0 {
            // Build beta from alpha + increments.
            let mut parts = vec![0u32; n_rows];
            for r in 0..n_rows {
                parts[r] = alpha.part(r) + increments[r];
            }
            results.push(Partition::new(parts));
        }
        return;
    }

    // Max increment at this row:
    //   - can't exceed lambda[row] - alpha[row]  (stay in lambda)
    //   - horizontal strip: c[row] ≤ alpha[row-1] - alpha[row]  (for row > 0)
    //   - can't exceed remaining
    let max_from_lambda = lambda.part(row).saturating_sub(alpha.part(row));
    let max_from_strip = if row == 0 {
        remaining  // no upper-row constraint for the first row
    } else {
        alpha.part(row - 1).saturating_sub(alpha.part(row))
    };
    let max_c = remaining.min(max_from_lambda).min(max_from_strip);

    for c in 0..=max_c {
        increments[row] = c;
        enumerate_strips(alpha, lambda, remaining - c, row + 1, n_rows, increments, results);
    }
    increments[row] = 0;
}

/// Compute K(lambda/mu, w) using the level-by-level DP.
/// Returns 0 if lambda/mu/w are incompatible (sizes don't match, mu not contained in lambda, etc.).
///
/// `max_states`: if Some(limit), abort with an error message if any DP level exceeds `limit` states.
/// `sort_weight`: if true, sort w in decreasing order before the DP to minimise peak state count.
///   Only pass true when no flag bounds are active (sorting is always valid for K, but callers
///   that use w-ordering for degree/flag bookkeeping should pass false).
pub fn skew_kostka(
    lambda: &Partition,
    mu: &Partition,       // inner shape; use Partition::empty() for non-skew
    w: &[u32],            // weight composition
    max_states: Option<usize>,
    sort_weight: bool,
) -> BigUint {
    // Validate.
    let skew_size: u32 = lambda.size().saturating_sub(mu.size());
    let w_size: u32 = w.iter().sum();
    if skew_size != w_size {
        return BigUint::zero();
    }
    if !mu.contained_in(lambda) {
        return BigUint::zero();
    }

    // Optionally sort weight descending to reduce peak intermediate state count.
    let mut w_sorted;
    let w_eff: &[u32] = if sort_weight {
        w_sorted = w.to_vec();
        w_sorted.sort_unstable_by(|a, b| b.cmp(a));
        &w_sorted
    } else {
        w
    };

    // Initial DP state: the single partition mu with count 1.
    let mut dp: HashMap<Partition, BigUint> = HashMap::new();
    dp.insert(mu.clone(), BigUint::one());

    for &strip_size in w_eff {
        if strip_size == 0 {
            // No boxes added; state is unchanged.
            continue;
        }

        // Build new DP map by merging extensions directly, avoiding an intermediate Vec.
        let mut new_dp: HashMap<Partition, BigUint> = HashMap::new();
        for (alpha, count) in &dp {
            for beta in horizontal_strip_extensions(alpha, lambda, strip_size) {
                *new_dp.entry(beta).or_insert_with(BigUint::zero) += count;
            }
        }

        if let Some(limit) = max_states {
            if new_dp.len() > limit {
                eprintln!(
                    "Error: DP state count {} exceeds --max-states {}. \
                     Use a smaller input or raise the limit.",
                    new_dp.len(), limit
                );
                std::process::exit(1);
            }
        }

        dp = new_dp;
    }

    // The answer is the count at lambda.
    dp.remove(lambda).unwrap_or_else(BigUint::zero)
}

/// Convenience wrapper for non-skew K(lambda, w).
pub fn kostka(lambda: &Partition, w: &[u32], max_states: Option<usize>) -> BigUint {
    skew_kostka(lambda, &Partition::empty(), w, max_states, true)
}

// ── Flagged Kostka numbers ────────────────────────────────────────────────────
//
// A flagged SSYT restricts which rows each label may appear in:
//   upper_flags[i] = f  →  label i+1 only in rows 1..=f  (1-indexed)
//                      ⟺  strip i must have c[j] = 0 for j ≥ f  (0-indexed)
//   lower_flags[i] = g  →  label i+1 only in rows g..=n
//                      ⟺  strip i must have c[j] = 0 for j < g-1 (0-indexed)
//
// Both flags are per-step (length k = len(w)).  Missing entries mean no restriction.
// The weight cannot be sorted when flags are active.

/// Like `horizontal_strip_extensions` but restricted to rows `row_lo..row_hi` (0-indexed).
/// Rows outside this range are forced to zero increment.
fn horizontal_strip_extensions_restricted(
    alpha: &Partition,
    lambda: &Partition,
    strip_size: u32,
    row_lo: usize,
    row_hi: usize,
) -> Vec<Partition> {
    if strip_size == 0 {
        return vec![alpha.clone()];
    }
    let n_rows = lambda.num_parts();
    if n_rows == 0 {
        return vec![];
    }
    let mut results = Vec::new();
    let mut increments = vec![0u32; n_rows];
    enumerate_strips_restricted(
        alpha, lambda, strip_size, 0, n_rows,
        row_lo, row_hi, &mut increments, &mut results,
    );
    results
}

fn enumerate_strips_restricted(
    alpha: &Partition,
    lambda: &Partition,
    remaining: u32,
    row: usize,
    n_rows: usize,
    row_lo: usize,
    row_hi: usize,
    increments: &mut Vec<u32>,
    results: &mut Vec<Partition>,
) {
    if row == n_rows {
        if remaining == 0 {
            let mut parts = vec![0u32; n_rows];
            for r in 0..n_rows {
                parts[r] = alpha.part(r) + increments[r];
            }
            results.push(Partition::new(parts));
        }
        return;
    }

    // Rows outside the allowed range: forced zero.
    if row < row_lo || row >= row_hi {
        increments[row] = 0;
        enumerate_strips_restricted(
            alpha, lambda, remaining, row + 1, n_rows,
            row_lo, row_hi, increments, results,
        );
        return;
    }

    let max_from_lambda = lambda.part(row).saturating_sub(alpha.part(row));
    let max_from_strip = if row == 0 {
        remaining
    } else {
        alpha.part(row - 1).saturating_sub(alpha.part(row))
    };
    let max_c = remaining.min(max_from_lambda).min(max_from_strip);

    for c in 0..=max_c {
        increments[row] = c;
        enumerate_strips_restricted(
            alpha, lambda, remaining - c, row + 1, n_rows,
            row_lo, row_hi, increments, results,
        );
    }
    increments[row] = 0;
}

/// Compute the flagged skew Kostka coefficient K_flags(lambda/mu, w).
///
/// `upper_flags[i] = f`: label i+1 restricted to rows 1..=f (1-indexed).
/// `lower_flags[i] = g`: label i+1 restricted to rows g..=n (1-indexed).
/// Flags are per-step; lengths should equal len(w).  Missing entries = no restriction.
pub fn flagged_skew_kostka(
    lambda: &Partition,
    mu: &Partition,
    w: &[u32],
    upper_flags: Option<&[u32]>,
    lower_flags: Option<&[u32]>,
    max_states: Option<usize>,
) -> BigUint {
    let skew_size: u32 = lambda.size().saturating_sub(mu.size());
    let w_size: u32 = w.iter().sum();
    if skew_size != w_size {
        return BigUint::zero();
    }
    if !mu.contained_in(lambda) {
        return BigUint::zero();
    }

    let n = lambda.num_parts();
    let mut dp: HashMap<Partition, BigUint> = HashMap::new();
    dp.insert(mu.clone(), BigUint::one());

    for (i, &strip_size) in w.iter().enumerate() {
        if strip_size == 0 {
            continue;
        }

        // Allowed row range for this step (0-indexed, exclusive upper bound).
        let row_lo = lower_flags
            .and_then(|lf| lf.get(i))
            .map(|&g| (g as usize).saturating_sub(1).min(n))
            .unwrap_or(0);
        let row_hi = upper_flags
            .and_then(|uf| uf.get(i))
            .map(|&f| (f as usize).min(n))
            .unwrap_or(n);

        let mut new_dp: HashMap<Partition, BigUint> = HashMap::new();
        for (alpha, count) in &dp {
            for beta in horizontal_strip_extensions_restricted(
                alpha, lambda, strip_size, row_lo, row_hi,
            ) {
                *new_dp.entry(beta).or_insert_with(BigUint::zero) += count;
            }
        }

        if let Some(limit) = max_states {
            if new_dp.len() > limit {
                eprintln!(
                    "Error: DP state count {} exceeds --max-states {}.",
                    new_dp.len(), limit
                );
                std::process::exit(1);
            }
        }

        dp = new_dp;
    }

    dp.remove(lambda).unwrap_or_else(BigUint::zero)
}

// ── Strict (interior) counting for Ehrhart-Macdonald reciprocity ─────────────
//
// A chain μ=α⁰ ⊂ α¹ ⊂ … ⊂ αᵏ=λ is in the *relative interior* of the
// GT-polytope iff every interlacing constraint that is NOT in the affine hull
// of the polytope holds strictly.
//
// A constraint (between adjacent chain levels) is in the affine hull iff it is
// always an equality over all feasible chains.  After full propagation by the
// dimension algorithm (gt_polytope_bounds), this is equivalent to both
// endpoints being frozen at the same value:
//   globally_tight = (lb_src == ub_src) && (lb_dst == ub_dst) && (lb_src == lb_dst)
//
// The two constraints per step i (0-indexed), row j:
//   lower[i][j]:   α^{i+1}[j] >= α^i[j]
//   diagonal[i][j] (j≥1): α^{i+1}[j] <= α^i[j-1]
//
// For non-globally-tight lower:   require c[j] = β[j] - α[j] >= 1
// For non-globally-tight diagonal: require α[j-1] - β[j] >= 1
//
// By Ehrhart-Macdonald reciprocity:
//   (-1)^d * strict_skew_kostka(t·λ, t·μ, t·w) = P_Ehrhart(-t)

/// Enumerate horizontal strip extensions β of α with per-row strictness flags.
///
/// `strict_lower[j]`: require β[j] > α[j] (c[j] ≥ 1).
/// `strict_diag[j]`  (j≥1): require α[j-1] > β[j] (c[j] ≤ gap - 1).
fn strict_horizontal_strip_extensions_new(
    alpha: &Partition,
    lambda: &Partition,
    strip_size: u32,
    strict_lower: &[bool],
    strict_diag: &[bool],
) -> Vec<Partition> {
    let n_rows = lambda.num_parts();
    if n_rows == 0 {
        return vec![];
    }

    // Minimum boxes needed at each suffix for pruning.
    let mut min_needed_suffix = vec![0u32; n_rows + 1];
    for r in (0..n_rows).rev() {
        min_needed_suffix[r] = min_needed_suffix[r + 1]
            + if r < strict_lower.len() && strict_lower[r] { 1 } else { 0 };
    }

    if strip_size < min_needed_suffix[0] {
        return vec![];
    }

    let mut results = Vec::new();
    let mut increments = vec![0u32; n_rows];
    enumerate_strips_strict_new(
        alpha,
        lambda,
        strip_size,
        0,
        n_rows,
        strict_lower,
        strict_diag,
        &min_needed_suffix,
        &mut increments,
        &mut results,
    );
    results
}

#[allow(clippy::too_many_arguments)]
fn enumerate_strips_strict_new(
    alpha: &Partition,
    lambda: &Partition,
    remaining: u32,
    row: usize,
    n_rows: usize,
    strict_lower: &[bool],
    strict_diag: &[bool],
    min_needed_suffix: &[u32],
    increments: &mut Vec<u32>,
    results: &mut Vec<Partition>,
) {
    if row == n_rows {
        if remaining == 0 {
            let mut parts = vec![0u32; n_rows];
            for r in 0..n_rows {
                parts[r] = alpha.part(r) + increments[r];
            }
            results.push(Partition::new(parts));
        }
        return;
    }

    // Prune: not enough remaining to satisfy all strict-lower rows ahead.
    if remaining < min_needed_suffix[row] {
        return;
    }

    let need_strict_lower = row < strict_lower.len() && strict_lower[row];
    let need_strict_diag  = row > 0 && row < strict_diag.len() && strict_diag[row];

    let min_c: u32 = if need_strict_lower { 1 } else { 0 };

    let max_from_lambda = lambda.part(row).saturating_sub(alpha.part(row));

    // Regular horizontal-strip diagonal bound.
    let gap = if row == 0 {
        u32::MAX
    } else {
        alpha.part(row - 1).saturating_sub(alpha.part(row))
    };
    // Apply strict diagonal: reduce gap by 1.
    let max_from_diag = if need_strict_diag {
        gap.saturating_sub(1)
    } else {
        gap
    };

    // Can't use more than remaining minus what future strict-lower rows still need.
    let max_for_row = remaining.saturating_sub(min_needed_suffix[row + 1]);

    let max_c = remaining.min(max_from_lambda).min(max_from_diag).min(max_for_row);

    if max_c < min_c {
        return;
    }

    for c in min_c..=max_c {
        increments[row] = c;
        enumerate_strips_strict_new(
            alpha,
            lambda,
            remaining - c,
            row + 1,
            n_rows,
            strict_lower,
            strict_diag,
            min_needed_suffix,
            increments,
            results,
        );
    }
    increments[row] = 0;
}

/// Count interior lattice points of the GT-polytope GT(λ/μ, w).
///
/// Uses the lb/ub bounds from `gt_polytope_bounds` to identify which interlacing
/// constraints are globally tight (part of the affine hull) and only enforces
/// strictness for the remaining ones.  This correctly implements the relative
/// interior condition for all polytope dimensions.
///
/// `sort_weight` is ignored: sorting would break the lb/ub correspondence.
pub fn strict_skew_kostka(
    lambda: &Partition,
    mu: &Partition,
    w: &[u32],
    max_states: Option<usize>,
    _sort_weight: bool,
) -> BigUint {
    let skew_size: u32 = lambda.size().saturating_sub(mu.size());
    let w_size: u32 = w.iter().sum();
    if skew_size != w_size {
        return BigUint::zero();
    }
    if !mu.contained_in(lambda) {
        return BigUint::zero();
    }

    let n = lambda.num_parts();
    let k = w.len();

    // Obtain propagated bounds; None means the polytope is empty.
    let (_, lb, ub) = match crate::gt_dim::gt_polytope_bounds(lambda.parts(), mu.parts(), w) {
        None => return BigUint::zero(),
        Some(data) => data,
    };

    // Pad μ to length n.
    let mut mu_pad = mu.parts().to_vec();
    mu_pad.resize(n, 0);
    let lambda_parts = lambda.parts();

    // Tight bounds for α^i[j]:
    //   i = 0      → boundary μ, both lb and ub equal μ_j
    //   i = 1..k-1 → interior level ell = i-1: lb[i-1][j], ub[i-1][j]
    //   i = k      → boundary λ, both equal λ_j
    let src_lb = |i: usize, j: usize| -> u32 {
        if i == 0 { mu_pad[j] } else { lb[i - 1][j] }
    };
    let src_ub = |i: usize, j: usize| -> u32 {
        if i == 0 { mu_pad[j] } else { ub[i - 1][j] }
    };
    let dst_lb = |i: usize, j: usize| -> u32 {
        if i + 1 == k { lambda_parts[j] } else { lb[i][j] }
    };
    let dst_ub = |i: usize, j: usize| -> u32 {
        if i + 1 == k { lambda_parts[j] } else { ub[i][j] }
    };

    // Precompute strictness flags per step i and row j.
    // A constraint is globally tight iff both sides are frozen at the same value.
    let mut strict_lower = vec![vec![false; n]; k];
    let mut strict_diag  = vec![vec![false; n]; k];
    for i in 0..k {
        for j in 0..n {
            // Lower constraint: α^{i+1}[j] >= α^i[j]
            let (sl, su) = (src_lb(i, j), src_ub(i, j));
            let (dl, du) = (dst_lb(i, j), dst_ub(i, j));
            strict_lower[i][j] = !(sl == su && dl == du && sl == dl);

            // Diagonal constraint (j≥1): α^{i+1}[j] <= α^i[j-1]
            if j >= 1 {
                let (sl2, su2) = (src_lb(i, j - 1), src_ub(i, j - 1));
                // dst side is the same α^{i+1}[j]
                strict_diag[i][j] = !(sl2 == su2 && dl == du && sl2 == dl);
            }
        }
    }

    // Level-by-level DP with per-row strictness.
    let mut dp: HashMap<Partition, BigUint> = HashMap::new();
    dp.insert(mu.clone(), BigUint::one());

    for (i, &strip_size) in w.iter().enumerate() {
        if strip_size == 0 {
            // Zero-weight step: all its constraints are globally tight (both sides
            // frozen at the same value after propagation).  No strictness needed;
            // just skip as in the regular DP.
            continue;
        }

        let mut new_dp: HashMap<Partition, BigUint> = HashMap::new();
        for (alpha, count) in &dp {
            for beta in strict_horizontal_strip_extensions_new(
                alpha,
                lambda,
                strip_size,
                &strict_lower[i],
                &strict_diag[i],
            ) {
                *new_dp.entry(beta).or_insert_with(BigUint::zero) += count;
            }
        }

        if let Some(limit) = max_states {
            if new_dp.len() > limit {
                eprintln!(
                    "Error: DP state count {} exceeds --max-states {}.",
                    new_dp.len(),
                    limit
                );
                std::process::exit(1);
            }
        }

        dp = new_dp;
    }

    dp.remove(lambda).unwrap_or_else(BigUint::zero)
}

/// Convenience wrapper for non-skew strict K(lambda, w).
pub fn strict_kostka(lambda: &Partition, w: &[u32], max_states: Option<usize>) -> BigUint {
    strict_skew_kostka(lambda, &Partition::empty(), w, max_states, false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;

    fn p(parts: &[u32]) -> Partition { Partition::new(parts.to_vec()) }
    fn biguint(n: u64) -> BigUint { BigUint::from(n) }

    // ── Weak Kostka numbers ────────────────────────────────────────────────────

    #[test]
    fn kostka_321_222() {
        // K(3,2,1 | 2,2,2) = 2  (mentioned in session)
        assert_eq!(kostka(&p(&[3,2,1]), &[2,2,2], None), biguint(2));
    }

    #[test]
    fn kostka_311_11111_t1() {
        // K(3,1,1 | 1,1,1,1,1) = 6  [Ehrhart poly P(1)]
        assert_eq!(kostka(&p(&[3,1,1]), &[1,1,1,1,1], None), biguint(6));
    }

    #[test]
    fn kostka_311_11111_t2() {
        // K(6,2,2 | 2,2,2,2,2) = 20  [P(2)]
        assert_eq!(kostka(&p(&[6,2,2]), &[2,2,2,2,2], None), biguint(20));
    }

    #[test]
    fn kostka_311_11111_t3() {
        // K(9,3,3 | 3,3,3,3,3) = 50  [P(3)]
        assert_eq!(kostka(&p(&[9,3,3]), &[3,3,3,3,3], None), biguint(50));
    }

    #[test]
    fn kostka_311_11111_t4() {
        // K(12,4,4 | 4,4,4,4,4) = 105  [P(4)]
        assert_eq!(kostka(&p(&[12,4,4]), &[4,4,4,4,4], None), biguint(105));
    }

    // ── Strict Kostka numbers (Ehrhart–Macdonald reciprocity) ─────────────────
    //
    // For shape (3,1,1) / weight (1,1,1,1,1), Ehrhart polynomial:
    //   P(n) = (2n^4 + 16n^3 + 46n^2 + 56n + 24) / 24
    // Macdonald reciprocity: (-1)^4 * K_strict(t*λ, t*w) = P(-t)
    //   P(-1) = 0, P(-2) = 0, P(-3) = 0, P(-4) = 1, P(-5) = 6, P(-6) = 20

    #[test]
    fn strict_311_11111_t1() {
        assert_eq!(strict_kostka(&p(&[3,1,1]), &[1,1,1,1,1], None), biguint(0));
    }

    #[test]
    fn strict_311_11111_t2() {
        assert_eq!(strict_kostka(&p(&[6,2,2]), &[2,2,2,2,2], None), biguint(0));
    }

    #[test]
    fn strict_311_11111_t3() {
        assert_eq!(strict_kostka(&p(&[9,3,3]), &[3,3,3,3,3], None), biguint(0));
    }

    #[test]
    fn strict_311_11111_t4() {
        assert_eq!(strict_kostka(&p(&[12,4,4]), &[4,4,4,4,4], None), biguint(1));
    }

    #[test]
    fn strict_311_11111_t5() {
        assert_eq!(strict_kostka(&p(&[15,5,5]), &[5,5,5,5,5], None), biguint(6));
    }

    #[test]
    fn strict_311_11111_t6() {
        assert_eq!(strict_kostka(&p(&[18,6,6]), &[6,6,6,6,6], None), biguint(20));
    }
}
