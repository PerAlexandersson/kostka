# Dimension of GT Polytopes via the Chain Model

## Overview

`src/gt_dim.rs` computes the dimension of the Gelfand-Tsetlin polytope
$\mathrm{GT}(\lambda/\mu, w)$.  This equals the degree of the Ehrhart polynomial
$n \mapsto K(n\lambda/n\mu,\, nw)$.

The function returns `Option<usize>`:
- `None` — the polytope is **empty** (infeasible constraints, $K = 0$ identically)
- `Some(0)` — a single lattice point (0-dimensional)
- `Some(d)` — a polytope of dimension $d \ge 1$

Both non-skew and skew shapes are supported, for arbitrary weight compositions
and optional row flags.


## The chain model

The polytope $\mathrm{GT}(\lambda/\mu, w)$ parameterizes chains of partitions

$$\mu = \alpha^0 \subset \alpha^1 \subset \cdots \subset \alpha^k = \lambda$$

where each successive difference $\alpha^i / \alpha^{i-1}$ is a **horizontal strip**
(at most one box per column) of size $w_i$, and $k = \mathrm{len}(w)$.  Each
intermediate partition $\alpha^i$ has at most $n = \mathrm{len}(\lambda)$ parts.

The **interior levels** are $\alpha^1, \ldots, \alpha^{k-1}$ (there are $k-1$ of them).
Each level has $n$ entries, giving a $(k-1) \times n$ grid of free variables indexed
by level $\ell = 0, \ldots, k-2$ and column $j = 0, \ldots, n-1$.


## Interlacing inequalities

The horizontal strip condition $\alpha^{i-1} \subset \alpha^i$ gives:

$$\alpha^i_0 \;\ge\; \alpha^{i-1}_0 \;\ge\; \alpha^i_1 \;\ge\; \alpha^{i-1}_1 \;\ge\; \cdots \;\ge\; \alpha^i_{n-1} \;\ge\; \alpha^{i-1}_{n-1} \;\ge\; 0$$

So each interior entry $\alpha^i_j$ is bounded by its four interlacing neighbors:

| direction       | bound                                          |
|-----------------|------------------------------------------------|
| level below (↓) | $\alpha^i_j \ge \alpha^{i-1}_j$               |
| level below (↙) | $\alpha^i_j \le \alpha^{i-1}_{j-1}$  ($j\ge 1$) |
| level above (↑) | $\alpha^i_j \le \alpha^{i+1}_j$               |
| level above (↗) | $\alpha^i_j \ge \alpha^{i+1}_{j+1}$           |


## Boundary bounds (initialization)

By chaining the interlacing inequalities all the way to the fixed boundary
levels $\alpha^0 = \mu$ and $\alpha^k = \lambda$, we obtain tight entry-wise bounds
that depend only on $\lambda$ and $\mu$:

$$\mathrm{lb}[\ell][j] = \max\!\big(\mu_j,\;\lambda_{j+k-(\ell+1)}\big)$$
$$\mathrm{ub}[\ell][j] = \min\!\big(\lambda_j,\;\mu_{j-(\ell+1)}\big)$$

with out-of-range indices contributing 0 to the lower bound and $+\infty$
(i.e., $\lambda_j$) to the upper bound.

**Derivation of the four chains:**
- **Vertical from bottom**: $\alpha^i_j \ge \alpha^{i-1}_j \ge \cdots \ge \mu_j$
- **Diagonal from top**: $\alpha^i_j \ge \alpha^{i+1}_{j+1} \ge \cdots \ge \lambda_{j+(k-i)}$
- **Vertical from top**: $\alpha^i_j \le \alpha^{i+1}_j \le \cdots \le \lambda_j$
- **Diagonal from bottom**: $\alpha^i_j \le \alpha^{i-1}_{j-1} \le \cdots \le \mu_{j-i}$

An entry is **frozen** when $\mathrm{lb} = \mathrm{ub}$ (its value is uniquely determined).
An infeasible entry ($\mathrm{lb} > \mathrm{ub}$) means the polytope is empty.


## Weight constraints

Each interior level $\ell$ has a fixed row sum (weight-sum target):
$$T_\ell = |\mu| + w_1 + \cdots + w_{\ell+1}$$

This constraint reduces the degrees of freedom: a level with $f$ free entries
contributes $\max(f-1, 0)$ to the dimension (one DOF is consumed by the
sum constraint).

Weight targets also **force** entries:
- If $T_\ell = \sum_j \mathrm{lb}[\ell][j]$: all free entries collapse to their lower bound.
- If $T_\ell = \sum_j \mathrm{ub}[\ell][j]$: all free entries collapse to their upper bound.
- If exactly one free entry remains: its value is determined by $T_\ell$ minus the
  sum of frozen entries.
- If $T_\ell < \sum \mathrm{lb}$ or $T_\ell > \sum \mathrm{ub}$: polytope is empty → return `None`.


## Flag constraints

Row flags restrict which rows of the SSYT may contain a given label, which
translates to equality constraints between adjacent interior levels:

- **upper_flags[ℓ] = f**: label $\ell+1$ appears only in rows $1\ldots f$ (1-indexed)
  $\Longleftrightarrow \alpha^{\ell+1}_j = \alpha^\ell_j$ for $j \ge f$ (0-indexed).

- **lower_flags[ℓ] = g**: label $\ell+1$ appears only in rows $g\ldots n$
  $\Longleftrightarrow \alpha^{\ell+1}_j = \alpha^\ell_j$ for $j < g-1$ (0-indexed).

These equalities are enforced by **bidirectional interval tightening**: for each
constrained column $j$ at the boundary between levels $\ell$ and $\ell+1$,

$$\mathrm{lb}[\ell][j] \leftarrow \max(\mathrm{lb}[\ell][j],\; \mathrm{lb}[\ell-1][j])$$
$$\mathrm{ub}[\ell][j] \leftarrow \min(\mathrm{ub}[\ell][j],\; \mathrm{ub}[\ell-1][j])$$

and the tightened bounds are propagated back to both levels.  Edge cases:
- $\ell = 0$: the lower boundary is $\mu$, so force $\mathrm{lb}[0][j] = \mathrm{ub}[0][j] = \mu_j$.
- $\ell = k-1$: the upper boundary is $\lambda$, so force $\mathrm{lb}[k-2][j] = \mathrm{ub}[k-2][j] = \lambda_j$.


## Algorithm (pseudocode)

```
function gt_polytope_dim(λ, μ, w, upper_flags, lower_flags):

    n ← len(λ),  k ← len(w)
    if n = 0 or k ≤ 1:  return Some(0)

    pad μ to length n with zeros
    compute weight-sum targets T[ℓ] = |μ| + Σ_{i≤ℓ+1} w_i   for ℓ = 0..k-2

    // Pass 0: initialize boundary bounds
    for ℓ = 0 to k-2, j = 0 to n-1:
        i ← ℓ + 1
        lb[ℓ][j] ← max(μ_j,  λ_{j+k-i}  if j+k-i < n  else 0)
        ub[ℓ][j] ← min(λ_j,  μ_{j-i}    if j ≥ i      else λ_j)
        if lb[ℓ][j] > ub[ℓ][j]:  return None   // empty

    // Main loop: propagate until stable
    repeat until no bounds change:

        // Pass 1: propagate frozen entries via interlacing
        for ℓ = 0 to k-2, j = 0 to n-1:
            if lb[ℓ][j] = ub[ℓ][j]:
                v ← lb[ℓ][j]
                if ℓ+1 < k-1:
                    lb[ℓ+1][j]   ← max(lb[ℓ+1][j],   v)    // α^{i+1}_j ≥ v
                    ub[ℓ+1][j+1] ← min(ub[ℓ+1][j+1], v)    // α^{i+1}_{j+1} ≤ v
                if ℓ > 0:
                    ub[ℓ-1][j]   ← min(ub[ℓ-1][j],   v)    // α^{i-1}_j ≤ v
                    lb[ℓ-1][j-1] ← max(lb[ℓ-1][j-1], v)    // α^{i-1}_{j-1} ≥ v

        check all lb[ℓ][j] ≤ ub[ℓ][j], else return None

        // Pass 2: weight forcing at each level
        for ℓ = 0 to k-2:
            free ← {j : lb[ℓ][j] < ub[ℓ][j]}
            min_sum ← Σ_j lb[ℓ][j],  max_sum ← Σ_j ub[ℓ][j]
            if T[ℓ] < min_sum or T[ℓ] > max_sum:  return None
            if T[ℓ] = min_sum:  ub[ℓ][j] ← lb[ℓ][j]  for all j in free
            if T[ℓ] = max_sum:  lb[ℓ][j] ← ub[ℓ][j]  for all j in free
            if |free| = 1:  force lb[ℓ][j] = ub[ℓ][j] ← T[ℓ] - Σ_{frozen} lb[ℓ][j]

        // Pass 3: flag constraints (bidirectional tightening)
        for each (ℓ, j) constrained by a flag:
            if ℓ = 0:
                lb[0][j] ← ub[0][j] ← μ_j
            else if ℓ < k-1:
                new_lb ← max(lb[ℓ][j], lb[ℓ-1][j])
                new_ub ← min(ub[ℓ][j], ub[ℓ-1][j])
                lb[ℓ][j] ← lb[ℓ-1][j] ← new_lb
                ub[ℓ][j] ← ub[ℓ-1][j] ← new_ub
            else:  // ℓ = k-1
                lb[k-2][j] ← ub[k-2][j] ← λ_j

        check all lb[ℓ][j] ≤ ub[ℓ][j], else return None

    // Compute dimension
    dim ← 0
    for ℓ = 0 to k-2:
        f ← |{j : lb[ℓ][j] < ub[ℓ][j]}|
        dim += max(f - 1, 0)
    return Some(dim)
```

**Complexity**: the propagation loop freezes at least one entry per iteration,
so it runs at most $(k-1) \times n$ times.  Each iteration scans $O(kn)$ entries.
Worst case $O((kn)^2)$, but in practice 1–3 iterations suffice because boundary
bounds already capture most frozen entries.


## The dimension formula

The dimension is:

$$\dim = \sum_{\ell=0}^{k-2} \max(f_\ell - 1,\; 0)$$

where $f_\ell$ is the number of entries at level $\ell$ with $\mathrm{lb} < \mathrm{ub}$
after propagation.

The "minus one" accounts for the weight constraint at each active level
(a level with all entries frozen contributes 0 regardless).


## Test coverage

The implementation is verified against empirical Ehrhart degrees (computed via
finite differences of actual Kostka numbers) for:

- All non-skew partitions with $|\lambda| \le 7$ and unit weight
- Various non-unit weights for shapes up to $|\lambda| = 5$
- All skew shapes with $|\lambda| \le 7$, $\mu \subset \lambda$, $|\mu| \le |\lambda| - 2$,
  over all partition weights — over 1600 cases
- 30 hand-picked larger cases with $8 \le |\lambda| \le 15$, both skew and non-skew
