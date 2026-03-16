/// Ehrhart polynomial computation for GT(lambda/mu, w).
///
/// By Rassart (2004), the function n ↦ K(n*lambda / n*mu, n*w) is a polynomial in n.
/// We:
///   1. Compute the degree d via gt_dim::gt_polytope_dim.
///   2. Collect d+1 sample points using Ehrhart-Macdonald reciprocity with an
///      adaptive strategy (or plain positive-dilation when flags are active).
///   3. Solve the resulting system over Q by Gaussian elimination.
///
/// The polynomial is stored as a Vec<BigRational> of length d+1,
/// where poly[k] is the coefficient of n^k (index 0 = constant term).

use num_bigint::{BigInt, BigUint, ToBigInt};
use num_rational::BigRational;
use num_traits::{Zero, One};
use rayon::prelude::*;
use crate::partition::Partition;
use crate::kostka_dp::{skew_kostka, strict_skew_kostka};
use crate::gt_dim::gt_polytope_dim;

pub struct EhrhartPoly {
    /// Coefficients of the polynomial in n: poly[k] = coeff of n^k.
    pub coeffs: Vec<BigRational>,
    /// Degree (index of highest non-zero coefficient).
    pub degree: usize,
}

impl EhrhartPoly {
    /// Evaluate the polynomial at a given n.
    pub fn eval(&self, n: u64) -> BigRational {
        let n_r = BigRational::from(BigInt::from(n));
        let mut result = BigRational::zero();
        let mut power = BigRational::one();
        for c in &self.coeffs {
            result += c * &power;
            power *= &n_r;
        }
        result
    }

    /// True if any coefficient of the polynomial is negative.
    pub fn has_negative_coefficient(&self) -> bool {
        self.coeffs.iter().any(|c| *c < BigRational::zero())
    }

    /// Display as a human-readable string, e.g. "(1/3)n^3 + n^2 + (5/3)n + 1".
    pub fn display(&self) -> String {
        let d = self.degree;
        if d == 0 {
            return format_rat(&self.coeffs[0]);
        }
        let mut terms: Vec<String> = Vec::new();
        for k in (0..=d).rev() {
            let c = &self.coeffs[k];
            if c.is_zero() { continue; }
            let c_str = format_rat(c);
            let term = match k {
                0 => c_str,
                1 => if c == &BigRational::one() { "n".into() } else { format!("{}n", c_str) },
                _ => if c == &BigRational::one() { format!("n^{}", k) } else { format!("{}n^{}", c_str, k) },
            };
            terms.push(term);
        }
        if terms.is_empty() { "0".into() } else { terms.join(" + ") }
    }

    /// Display as (1/d!) * (integer polynomial), where d = degree.
    /// Multiplying each coefficient by d! always yields integers for GT Ehrhart polynomials.
    pub fn display_factored(&self) -> String {
        let d = self.degree;
        if self.coeffs.iter().all(|c| c.is_zero()) {
            return "0".into();
        }

        // Compute d!
        let d_fact: BigInt = (1..=d as u64).fold(BigInt::one(), |acc, i| acc * BigInt::from(i));
        let d_fact_r = BigRational::from(d_fact.clone());

        // Integer coefficients: c_k * d!
        let int_coeffs: Vec<BigInt> = self.coeffs[..=d].iter()
            .map(|c| (c * &d_fact_r).to_integer())
            .collect();

        let poly_str = format_int_poly(&int_coeffs);

        if d <= 1 {
            poly_str
        } else {
            format!("(1/{}!) * ({})", d, poly_str)
        }
    }
}

fn format_rat(r: &BigRational) -> String {
    if r.denom() == &BigInt::one() {
        r.numer().to_string()
    } else {
        format!("({}/{})", r.numer(), r.denom())
    }
}

/// Format a polynomial with integer coefficients as "a_d n^d + ... + a_1 n + a_0".
fn format_int_poly(coeffs: &[BigInt]) -> String {
    let d = coeffs.len().saturating_sub(1);
    let mut terms: Vec<String> = Vec::new();
    for k in (0..=d).rev() {
        let c = &coeffs[k];
        if c.is_zero() { continue; }
        let c_abs = c.magnitude().clone();
        let sign: String = if terms.is_empty() {
            if *c < BigInt::zero() { "-".into() } else { "".into() }
        } else {
            if *c < BigInt::zero() { " - ".into() } else { " + ".into() }
        };
        let mag_str = if c_abs == num_bigint::BigUint::from(1u32) && k > 0 {
            "".into()
        } else {
            c_abs.to_string()
        };
        let var_str = match k {
            0 => "".into(),
            1 => "n".into(),
            _ => format!("n^{}", k),
        };
        terms.push(format!("{}{}{}", sign, mag_str, var_str));
    }
    if terms.is_empty() { "0".into() } else { terms.join("") }
}

/// Compute the Ehrhart polynomial of GT(lambda/mu, w).
///
/// When `use_reciprocity` is true (the default), the computation uses
/// Ehrhart-Macdonald reciprocity with a sequential adaptive strategy:
///   - P(0) = 1 is free.
///   - At each step, choose the side (positive or negative t) whose last
///     raw DP count was smaller, evaluating one point at a time until
///     d+1 points are collected.  Negative side wins ties.
///   - P(-t) = (-1)^d * K_strict(t·λ, t·μ, t·w)
///   - Strict Kostka is 0 for small t when w has many 1s, making those
///     evaluations free and heavily favouring the negative side early on.
///
/// Reciprocity requires no row flags; with flags the method silently
/// falls back to the plain positive-dilation scheme.
pub fn compute_ehrhart(
    lambda: &Partition,
    mu: &Partition,
    w: &[u32],
    upper_flags: Option<&[u32]>,
    lower_flags: Option<&[u32]>,
    _verbose: bool,
    max_states: Option<usize>,
    use_reciprocity: bool,
) -> EhrhartPoly {
    // Early exit: sizes must be compatible for a non-empty polytope.
    let skew_size = lambda.size().saturating_sub(mu.size());
    let w_size: u32 = w.iter().sum();
    if skew_size != w_size {
        return EhrhartPoly { coeffs: vec![BigRational::zero()], degree: 0 };
    }

    // Only reorder w when no flag bounds are active (flags are tied to w's row ordering).
    // Reciprocity also requires no flags (strict DP has no flagged variant).
    let sort_weight = upper_flags.is_none() && lower_flags.is_none();
    let use_recip = use_reciprocity && sort_weight;

    let d = match gt_polytope_dim(lambda.parts(), mu.parts(), w) {
        None => return EhrhartPoly { coeffs: vec![BigRational::zero()], degree: 0 },
        Some(d) => d,
    };

    if use_recip {
        // P(0) = 1 always (the trivial chain μ=μ is the unique point at dilation 0).
        if d == 0 {
            // Constant polynomial.  No DP evaluation needed.
            return EhrhartPoly { coeffs: vec![BigRational::one()], degree: 0 };
        }

        // Sequential adaptive reciprocity strategy:
        //
        // P(0) = 1 is free.  Then greedily evaluate one point at a time, choosing
        // the side (positive t or negative t) whose previous raw DP count was smaller.
        // Negative side wins ties and goes first (strict Kostka ≤ ordinary Kostka,
        // and equals 0 for small t when w has many 1s — those zeros are free).
        //
        // P(-t) = (-1)^d * K_strict(t·λ, t·μ, t·w)  (Ehrhart-Macdonald reciprocity)
        let sign_pos = d % 2 == 0; // (-1)^d is +1 iff d is even

        let mut points: Vec<(i64, BigRational)> = Vec::with_capacity(d + 1);
        points.push((0, BigRational::one())); // P(0) = 1 is free

        let mut pos_t: u64 = 0;
        let mut neg_t: u64 = 0;
        let mut last_pos_count = BigUint::one();
        let mut last_neg_count = BigUint::one(); // equal → negative wins the first tie

        while points.len() <= d {
            if last_neg_count <= last_pos_count {
                // Evaluate the next negative point via strict Kostka
                neg_t += 1;
                let tl = scale_partition(lambda, neg_t);
                let tm_p = scale_partition(mu, neg_t);
                let tw: Vec<u32> = w.iter().map(|&x| x * neg_t as u32).collect();
                let ks = strict_skew_kostka(&tl, &tm_p, &tw, max_states, false);
                let ks_r = BigRational::from(ks.to_bigint().unwrap());
                let p_val = if sign_pos { ks_r } else { -ks_r };
                points.push((-(neg_t as i64), p_val));
                last_neg_count = ks;
            } else {
                // Evaluate the next positive point via ordinary Kostka
                pos_t += 1;
                let tl = scale_partition(lambda, pos_t);
                let tm_p = scale_partition(mu, pos_t);
                let tw: Vec<u32> = w.iter().map(|&x| x * pos_t as u32).collect();
                let k = skew_kostka(&tl, &tm_p, &tw, max_states, true);
                let k_r = BigRational::from(k.to_bigint().unwrap());
                points.push((pos_t as i64, k_r));
                last_pos_count = k;
            }
        }

        let coeffs = poly_interpolate(&points);
        let true_degree = coeffs.iter().enumerate().rev()
            .find(|(_, c)| !c.is_zero())
            .map(|(i, _)| i)
            .unwrap_or(0);
        EhrhartPoly { coeffs, degree: true_degree }
    } else {
        // Plain method: evaluate at n = 1, ..., d+1.
        let values: Vec<BigRational> = (1..=(d + 1) as u64)
            .into_par_iter()
            .map(|n| {
                let nl = scale_partition(lambda, n);
                let nm = scale_partition(mu, n);
                let nw: Vec<u32> = w.iter().map(|&x| x * n as u32).collect();
                let k = skew_kostka(&nl, &nm, &nw, max_states, sort_weight);
                BigRational::from(k.to_bigint().unwrap())
            })
            .collect();

        let coeffs = vandermonde_solve(&values);
        let true_degree = coeffs.iter().enumerate().rev()
            .find(|(_, c)| !c.is_zero())
            .map(|(i, _)| i)
            .unwrap_or(0);
        EhrhartPoly { coeffs, degree: true_degree }
    }
}

fn scale_partition(p: &Partition, n: u64) -> Partition {
    Partition::new(p.parts().iter().map(|&x| x * n as u32).collect())
}

/// Interpolate a polynomial from d+1 arbitrary (x, y) sample points.
/// Returns coefficients [c_0, c_1, ..., c_d] for P(n) = c_0 + c_1·n + ... + c_d·n^d.
/// Uses Gaussian elimination over ℚ (d is small, typically ≤ 13).
fn poly_interpolate(points: &[(i64, BigRational)]) -> Vec<BigRational> {
    let d = points.len(); // d+1 unknowns, so this yields a degree-(d-1) polynomial
    // Build the augmented Vandermonde matrix: M[i][k] = x_i^k, last column = y_i.
    let mut mat: Vec<Vec<BigRational>> = points
        .iter()
        .map(|&(x, ref y)| {
            let xb = BigInt::from(x);
            let mut row: Vec<BigRational> = Vec::with_capacity(d + 1);
            let mut power = BigInt::one();
            for _ in 0..d {
                row.push(BigRational::from(power.clone()));
                power *= &xb;
            }
            row.push(y.clone());
            row
        })
        .collect();

    // Gauss-Jordan elimination.
    for col in 0..d {
        let pivot_row = (col..d)
            .find(|&r| !mat[r][col].is_zero())
            .expect("poly_interpolate: singular system (duplicate x-values?)");
        mat.swap(col, pivot_row);
        let pivot = mat[col][col].clone();
        for j in col..=d {
            let v = mat[col][j].clone() / &pivot;
            mat[col][j] = v;
        }
        for row in 0..d {
            if row == col { continue; }
            let factor = mat[row][col].clone();
            if factor.is_zero() { continue; }
            for j in col..=d {
                let sub = factor.clone() * &mat[col][j];
                mat[row][j] -= sub;
            }
        }
    }
    mat.iter().map(|row| row[d].clone()).collect()
}

/// Fit P to sample points (1, v_1), (2, v_2), ..., (d+1, v_{d+1}).
fn vandermonde_solve(values: &[BigRational]) -> Vec<BigRational> {
    let points: Vec<(i64, BigRational)> = values
        .iter()
        .enumerate()
        .map(|(i, v)| ((i + 1) as i64, v.clone()))
        .collect();
    poly_interpolate(&points)
}

/// Verify Ehrhart-Macdonald reciprocity for `n_checks` values of t.
///
/// For each t = 1 ..= n_checks, checks:
///   (-1)^d * strict_skew_kostka(t*λ, t*μ, t*w)  ==  P_Ehrhart(-t)
///
/// Prints a summary and returns true iff all checks pass.
pub fn verify_reciprocity(
    poly: &EhrhartPoly,
    lambda: &Partition,
    mu: &Partition,
    w: &[u32],
    n_checks: u64,
    max_states: Option<usize>,
) -> bool {
    let d = poly.degree;
    let sign_pos = d % 2 == 0; // (-1)^d is +1 iff d is even
    let sort_weight = true;

    // For degree-0 polytopes (a single point), the relative interior IS the point:
    // L_{P°}(t) = L_P(t) = constant.  Strict GT patterns always give 0, not 1, so
    // we cannot use the strict DP here — the check is trivially satisfied by P.
    if d == 0 {
        println!("  (degree 0: relative interior = polytope, check is trivial)");
        return true;
    }

    let mut all_ok = true;
    for t in 1..=n_checks {
        let tl = scale_partition(lambda, t);
        let tm = scale_partition(mu, t);
        let tw: Vec<u32> = w.iter().map(|&x| x * t as u32).collect();

        let interior = strict_skew_kostka(&tl, &tm, &tw, max_states, sort_weight);
        let interior_r = BigRational::from(interior.to_bigint().unwrap());

        // (-1)^d * interior
        let lhs = if sign_pos { interior_r.clone() } else { -interior_r.clone() };

        // P_Ehrhart(-t)
        let neg_t = BigRational::from(BigInt::from(-(t as i64)));
        let mut rhs = BigRational::zero();
        let mut power = BigRational::one();
        for c in &poly.coeffs {
            rhs += c * &power;
            power *= &neg_t;
        }

        let ok = lhs == rhs;
        if !ok {
            println!(
                "  FAIL t={}: (-1)^{} * L_{{P°}}({}) = {} but P(-{}) = {}",
                t, d, t, if sign_pos { interior_r.clone() } else { -interior_r.clone() }, t, rhs
            );
            all_ok = false;
        } else {
            println!("  OK   t={}: (-1)^{} * L_{{P°}}({}) = P(-{}) = {}", t, d, t, t, lhs);
        }
    }
    all_ok
}

/// Compute the h*-vector from the Ehrhart polynomial.
/// The h*-vector satisfies:
///   Σ_{n≥0} P(n) t^n  =  (Σ h*_k t^k) / (1-t)^{d+1}
/// Equivalently, express P(n) in the binomial basis C(n+d,d), C(n+d-1,d), ..., C(n,d):
///   P(n) = Σ_{k=0}^{d} h*_k * C(n+d-k, d)
/// The h*_k are the coefficients in this basis change.
pub fn compute_hstar(poly: &EhrhartPoly) -> Vec<BigInt> {
    let d = poly.degree;
    // Evaluate P(0), P(1), ..., P(d) (P(0) = h*_0 always).
    // Use the forward difference operator: h*_k = Δ^k P(0) / k! * ... (via finite differences).
    // Simpler: build the (d+1)×(d+1) basis-change matrix and solve.
    // For small d, just evaluate and use the known conversion:
    //   h*_k = Σ_{j=0}^{k} (-1)^{k-j} C(d+1, k-j) P(j)
    // This is the standard Ehrhart h*-vector formula.

    let d1 = d + 1;
    let mut hstar = vec![BigInt::zero(); d1];
    for k in 0..d1 {
        let mut val = BigInt::zero();
        for j in 0..=k {
            let p_j = poly.eval(j as u64);
            // p_j should be an integer (Kostka number at n=j).
            let p_j_int = p_j.numer().clone(); // denom should be 1
            let binom = binom_int(d1 as i64, (k - j) as i64);
            let sign = if (k - j) % 2 == 0 { BigInt::one() } else { -BigInt::one() };
            val += sign * binom * p_j_int;
        }
        hstar[k] = val;
    }
    hstar
}

fn binom_int(n: i64, k: i64) -> BigInt {
    if k < 0 || k > n { return BigInt::zero(); }
    let mut result = BigInt::one();
    for i in 0..k {
        result = result * BigInt::from(n - i) / BigInt::from(i + 1);
    }
    result
}

pub fn is_palindromic(v: &[BigInt]) -> bool {
    let v = trim_trailing_zeros(v);
    let n = v.len();
    (0..n / 2).all(|i| v[i] == v[n - 1 - i])
}

pub fn is_unimodal(v: &[BigInt]) -> bool {
    let v = trim_trailing_zeros(v);
    let n = v.len();
    if n <= 1 { return true; }
    let peak = v.iter().enumerate().max_by_key(|(_, x)| (*x).clone()).map(|(i, _)| i).unwrap_or(0);
    (0..peak).all(|i| v[i] <= v[i + 1]) && (peak..n - 1).all(|i| v[i] >= v[i + 1])
}

fn trim_trailing_zeros(v: &[BigInt]) -> &[BigInt] {
    let end = v.iter().rposition(|x| !x.is_zero()).map(|i| i + 1).unwrap_or(0);
    &v[..end]
}
