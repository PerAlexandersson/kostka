# GT-Polytope Reference

## Source
De Loera & McAllister, "Vertices of Gelfand–Tsetlin Polytopes", Discrete Comput. Geom. 32:459–470 (2004).
Local copy: `gelfandtsetlin.pdf`

---

## GT-Pattern (Definition 1.1)

A **GT-pattern** is a triangular array $(x_{ij})_{1 \le i \le j \le n}$ with $x_{ij} \in \mathbb{R}$ satisfying:
- $x_{ij} \ge 0$ for all $1 \le i \le j \le n$
- $x_{i,j+1} \ge x_{ij} \ge x_{i+1,j+1}$ for all $1 \le i \le j \le n-1$ (interleaving)

**Layout** (rows counted from bottom, row $j$ has entries $x_{1j}, \ldots, x_{jj}$):
```
x_{1n}  x_{2n}  ...  x_{nn}      ← top row (row n), fixed to λ
  x_{1,n-1}  ...  x_{n-1,n-1}
      ...
    x_{12}  x_{22}
        x_{11}                    ← bottom row (row 1), fixed to μ₁
```

## GT-Polytope (Definition 1.2)

$GT(\lambda, \mu)$ is the set of GT-patterns where:
- **Top row fixed**: $x_{in} = \lambda_i$ for $1 \le i \le n$
- **Row sums fixed**: $\sum_{i=1}^{j} x_{ij} = \mu_j$ for $1 \le j \le n$
  (equivalently, $x_{11} = \mu_1$ and row $j$ sum minus row $j-1$ sum = $\mu_j - \mu_{j-1}$)

**Connection to Kostka numbers**: The integer lattice points of $GT(\lambda, \mu)$ are in bijection with $SSYT(\lambda, \mu)$. Hence $K(\lambda, \mu) = \#(GT(\lambda, \mu) \cap \mathbb{Z}^{n(n+1)/2})$.

**Ehrhart polynomial**: For integer $\lambda, \mu$, the function $m \mapsto \#(GT(m\lambda, m\mu) \cap \mathbb{Z}^{n(n+1)/2})$ is a polynomial in $m$ (Rassart 2004). This equals $K(m\lambda, m\mu)$.

---

## Kostka DP (chain-of-strips model)

$K(\lambda/\mu, w)$ counts chains $\mu = \alpha^0 \subset \alpha^1 \subset \cdots \subset \alpha^k = \lambda$ where each $\alpha^i / \alpha^{i-1}$ is a **horizontal strip** of size $w_i$ (at most one box per column). Implemented level-by-level in `src/kostka_dp.rs`.

---

## Row-Flagged GT-Polytopes

### Flags as row restrictions on labels

A **row flag** restricts which rows of the SSYT each label may occupy:

- **upper flag**: $\mathrm{uf}[\ell] = f$ means label $\ell+1$ (1-indexed) appears only in rows $1, \ldots, f$ of the SSYT.
- **lower flag**: $\mathrm{lf}[\ell] = g$ means label $\ell+1$ appears only in rows $g, \ldots, n$.

Both flags are indexed by label (length $k = \mathrm{len}(w)$). Missing entries mean no restriction.

### Translation to the chain model

In the chain $\mu = \alpha^0 \subset \alpha^1 \subset \cdots \subset \alpha^k = \lambda$, the strip $\alpha^\ell / \alpha^{\ell-1}$ contains exactly the boxes labelled $\ell$. A row flag on label $\ell$ therefore constrains which entries of $\alpha^\ell$ can differ from $\alpha^{\ell-1}$:

$$\mathrm{uf}[\ell-1] = f \;\Longrightarrow\; \alpha^\ell_j = \alpha^{\ell-1}_j \quad \text{for all } j \ge f \text{ (0-indexed)}$$

$$\mathrm{lf}[\ell-1] = g \;\Longrightarrow\; \alpha^\ell_j = \alpha^{\ell-1}_j \quad \text{for all } j < g-1 \text{ (0-indexed)}$$

In other words, **flags impose equality constraints between adjacent interior levels** of the chain.

### Effect on the GT-polytope dimension

The flagged polytope $\mathrm{GT}_{\mathrm{flags}}(\lambda/\mu, w)$ is the face of $\mathrm{GT}(\lambda/\mu, w)$ where these equalities hold. Flags can only reduce (or preserve) the dimension:

$$\dim \mathrm{GT}_{\mathrm{flags}} \;\le\; \dim \mathrm{GT}$$

Flags can also make the polytope **empty**: if a flag forces $\alpha^\ell_j$ to equal $\alpha^{\ell-1}_j$ for all $j$ at some level $\ell$ where the weight target $|\alpha^\ell| - |\alpha^{\ell-1}| = w_\ell > 0$, no feasible chain exists.

The dimension algorithm (Pass 3 in `src/gt_dim.rs`) enforces flag equalities by bidirectional bound tightening between adjacent levels, then continues with the same interlacing-and-weight-forcing propagation loop. The return type is `Option<usize>`: `None` = empty polytope, `Some(d)` = dimension $d$.

### Flagged Kostka DP

For the counting problem, flags translate directly into the strip enumeration: at step $\ell$, boxes can only be placed in the allowed rows. Concretely, the increment $c[j] = \alpha^\ell_j - \alpha^{\ell-1}_j$ is forced to zero outside the allowed row range. This is implemented in `flagged_skew_kostka` in `src/kostka_dp.rs`.

### Flagged vs. unflagged: a comparison

| Property | Unflagged $K(\lambda/\mu, w)$ | Flagged $K_{\mathrm{flags}}(\lambda/\mu, w)$ |
|---|---|---|
| Symmetric in $w$? | Yes ($s_{\lambda/\mu}$ is symmetric) | No (label order matters) |
| Weight sorting in DP? | Yes (reduces peak states) | No |
| Polytope | Full $\mathrm{GT}(\lambda/\mu, w)$ | Sub-polytope with equality constraints |
| Degree relation | $\dim \mathrm{GT}(\lambda/\mu, w)$ | $\le$ unflagged degree |

---

## Ehrhart–Macdonald Reciprocity

For a rational polytope $P$ of dimension $d$:
$$L_{P^\circ}(t) = (-1)^d L_P(-t)$$
where $L_P(t) = \#(tP \cap \mathbb{Z}^n)$ and $L_{P^\circ}$ counts **interior** (relative interior) lattice points.

**Application**: Instead of evaluating $K(m\lambda, m\mu)$ at $m = 1, 2, \ldots, d+1$ to interpolate the degree-$d$ polynomial, we can use negative evaluations $P(-t) = (-1)^d K_{\text{strict}}(t\lambda, t\mu)$ at $t = 1, \ldots, \lfloor d/2 \rfloor$ and positive evaluations at $m = 1, \ldots, \lceil d/2 \rceil + 1$, roughly halving the maximum dilation needed.

---

## Face Dimension via Tilings (Theorem 1.5 + Corollary 1.6)

### Tiling (Definition 1.3)
The **tiling** $\mathcal{P}$ of a GT-pattern $\mathbf{x}$ partitions $\{(i,j): 1 \le i \le j \le n\}$ into **tiles**: connected groups of entries that are equal and adjacent. Two entries $(i_1,j_1)$ and $(i_2,j_2)$ are adjacent if $(i_2,j_2) \in \{(i_1+1,j_1+1),(i_1,j_1+1),(i_1-1,j_1-1),(i_1,j_1-1)\}$.

### Free Tiles
A tile is **free** if it does not intersect:
- the bottom row ($j=1$, i.e., the single entry $x_{11}$), or
- the top row (entries $x_{in}$, $1 \le i \le n$), or
- any diagonal entry $(i,i)$ for $1 \le i \le n$.

(Free tiles = tiles not pinned by the boundary constraints.)

### Tiling Matrix
With free tiles $P_1, \ldots, P_s$, the **tiling matrix** $A_\mathcal{P}$ is the $(n-2) \times s$ matrix (rows indexed $j = 2, \ldots, n-1$):
$$a_{jk} = \#\{i : (i,j) \in P_k\}$$
i.e., $a_{jk}$ = number of entries in row $j$ of $\mathbf{x}$ that belong to free tile $P_k$.

### Theorem 1.5 (Key Result)
$$\dim(\ker A_\mathcal{P}) = \dim(\text{minimal face of } GT(\lambda,\mu) \text{ containing } \mathbf{x})$$

**Corollaries**:
- $\mathbf{x}$ is a **vertex** iff $\ker A_\mathcal{P} = \{0\}$ (trivial kernel), iff some $s \times s$ submatrix of $A_\mathcal{P}$ has nonzero determinant.
- $\mathbf{x}$ is in the **relative interior** iff $\dim(\ker A_\mathcal{P}) = d$ where $d = \dim GT(\lambda,\mu)$ is the polytope dimension (= Ehrhart degree).

### Algorithm for Interior Counting
To count interior lattice points of $GT(t\lambda, t\mu)$:
1. Enumerate all integer GT-patterns $\mathbf{x} \in GT(t\lambda, t\mu) \cap \mathbb{Z}^{\ldots}$
2. For each, compute the tiling $\mathcal{P}$, build $A_\mathcal{P}$, compute $\dim(\ker A_\mathcal{P})$
3. Count those with $\dim(\ker A_\mathcal{P}) = d$

**Current implementation status**: `strict_skew_kostka` in `src/kostka_dp.rs` uses an approximate "active rows" criterion that is **incorrect** for some cases (e.g., shape (3,1,1) with weight (1,1,1,1,1) fails at $t \ge 4$). The correct criterion requires computing $\dim(\ker A_\mathcal{P})$ per Theorem 1.5.

---

## Implementation Notes

- `src/kostka_dp.rs`: DP for $K(\lambda/\mu, w)$ (`skew_kostka`), flagged variant (`flagged_skew_kostka`), and `strict_skew_kostka` (approximate — see note on Theorem 1.5 above).
- `src/gt_dim.rs`: Chain-model dimension algorithm (`gt_polytope_dim`, `gt_polytope_dim_full`). Returns `Option<usize>`: `None` = empty polytope, `Some(d)` = dimension. Handles arbitrary $k$, skew shapes, and row flags.
- `src/ehrhart.rs`: Vandermonde interpolation to recover the Ehrhart polynomial; `verify_reciprocity` checks $(-1)^d \cdot K_{\text{strict}}(t\lambda, t\mu) = P(-t)$.
- `src/main.rs`: CLI with subcommands `kostka`, `skew` (supports `--upper-flags`, `--lower-flags`), `degree`, `ehrhart`, `hstar`, `strict-kostka`, `syt`, `table`.
