# kostka

A CLI utility for computing Kostka coefficients, skew Kostka coefficients, Ehrhart
polynomials of Gelfand-Tsetlin polytopes, and h*-vectors. All computations use exact
arithmetic over the integers or rationals.

---

## Notation

Throughout this tool, the following conventions apply:

- **λ** (`--lambda`): outer partition (the shape, or the outer shape in a skew pair)
- **μ** (`--mu`): inner partition (the hole in a skew shape λ/μ; omit for non-skew)
- **w** (`--weight`): the weight, a **composition** (ordered sequence of non-negative
  integers). The order of parts in w matters: different orderings give different GT
  polytopes and different results for row-flagged Schur functions. K(λ/μ, w) is
  invariant under permutations of w for standard (non-flagged) computations, but the
  tool retains the order throughout for consistency with the flagged setting.

Partitions and compositions are written as comma-separated integers, e.g. `3,2,1`.
Trailing zeros in w may be included to specify the alphabet size explicitly.

---

## Mathematical background

### Kostka coefficients  K(λ, w)

K(λ, w) counts the number of semistandard Young tableaux (SSYT) of shape λ and content w,
where λ is a partition and w is a composition of |λ|. Equivalently, K(λ, w) is the
multiplicity of the weight w in the GL(n) irreducible representation with highest weight λ.

### Skew Kostka coefficients  K(λ/μ, w)

K(λ/μ, w) counts SSYT of skew shape λ/μ and content w, where μ ⊆ λ are both partitions
and w is a composition of |λ| − |μ|.

### Gelfand-Tsetlin (GT) polytopes

For a skew shape λ/μ and weight w, the GT polytope GT(λ/μ, w) is the set of triangular
arrays (GT patterns) with:

- Top row fixed to λ
- Bottom row fixed to μ  (absent for the non-skew case)
- Interlacing inequalities between consecutive rows
- Row-sum constraints that fix the content to w

The lattice points of GT(λ/μ, w) are exactly the SSYT of shape λ/μ and content w:

    |GT(λ/μ, w) ∩ Z^m| = K(λ/μ, w)

### Ehrhart polynomial of GT(λ/μ, w)

The n-th dilate satisfies n · GT(λ/μ, w) = GT(nλ/nμ, nw), so the Ehrhart function is:

    i(GT(λ/μ, w), n) = K(nλ/nμ, nw)

By a theorem of Rassart (2004), this is a **polynomial** in n (not merely quasi-polynomial,
despite GT(λ/μ, w) being a rational polytope in general). The degree equals the dimension
of the GT polytope.

#### Degree of the Ehrhart polynomial

The dimension of GT(λ/μ, w) is computed via the **chain model**: interior levels
α^1, …, α^{k−1} each have n entries whose feasible ranges are tracked as intervals [lb, ub].
Forced entries (lb = ub) propagate through three mechanisms until stable:

1. **Boundary bounds**: lb[ℓ][j] = max(μ_j, λ_{j+k−ℓ−1}), ub[ℓ][j] = min(λ_j, μ_{j−ℓ−1}). Equal consecutive parts in λ or μ force entire triangular regions.
2. **Weight forcing**: the row-sum target at each level may force free entries to their bounds.
3. **Flag constraints** (if present): equality between adjacent levels freezes further entries.

If any lb > ub is ever detected the polytope is **empty** and `degree` reports this explicitly.
After propagation:

    dim = Σ_ℓ max(free_ℓ − 1, 0)

where free_ℓ is the number of entries at level ℓ with lb < ub.

### h*-vector

The h*-vector of GT(λ/μ, w) is extracted from the Ehrhart series:

    Σ_{n≥0} K(nλ/nμ, nw) t^n  =  h*(t) / (1 − t)^{d+1}

where d = deg. The coefficients h*_0, h*_1, ..., h*_d are non-negative integers. For GT
polytopes, the h*-vector is conjectured (and in many cases proved) to be unimodal.

### Row-flagged Kostka coefficients

A row flag restricts which rows of the SSYT each label may occupy:

- **upper_flags[i] = f**: label i+1 (1-indexed) may only appear in rows 1 through f.
- **lower_flags[i] = g**: label i+1 may only appear in rows g through n.

Both flags are per-label, so their length equals the number of parts of w.

#### Effect on the SSYT / chain model

In the chain model the SSYT corresponds to a chain
μ = α^0 ⊂ α^1 ⊂ ... ⊂ α^k = λ,
where α^i/α^{i−1} is a horizontal strip of size w_i (the boxes labelled i).
A row flag at step i restricts which rows gain new boxes:

| flag | rows that gain boxes | constraint on the chain |
|---|---|---|
| upper_flags[i] = f | rows 1..=f only | α^i_j = α^{i−1}_j for all j ≥ f (0-indexed) |
| lower_flags[i] = g | rows g..=n only | α^i_j = α^{i−1}_j for all j < g−1 (0-indexed) |

In other words, flags impose **equality constraints between adjacent levels** of the chain.

#### Effect on the GT polytope and its dimension

The GT polytope GT(λ/μ, w) with flags is the sub-polytope where the equality constraints
above hold. This can only reduce the dimension: flags freeze additional interior entries
(or create infeasibilities that make the polytope empty). The `degree` subcommand handles
both flagged and unflagged cases through the same constraint-propagation algorithm:

1. Boundary bounds are initialized from λ, μ, and interlacing.
2. Flag constraints tighten adjacent-level bounds bidirectionally (Pass 3 in the chain model).
3. Propagation repeats until stable; infeasibility → empty polytope → returns `None`.
4. Dimension = Σ_ℓ max(free_ℓ − 1, 0) over all interior levels.

### Standard Young tableaux (SYT)

For weight w = (1, 1, ..., 1), SSYT reduce to SYT. The count is given by the hook-length
formula:

    f^λ = n! / ∏_{(i,j) ∈ λ} hook(i, j)

This is computed in O(|λ|) without the GT machinery.

---

## Algorithms

### Skew Kostka DP

K(λ/μ, w) is computed by a level-by-level dynamic program over Young's lattice. A SSYT of
shape λ/μ and content w = (w_1, ..., w_k) corresponds bijectively to a chain:

    μ = α^0 ⊂ α^1 ⊂ ... ⊂ α^k = λ

where each α^i / α^{i−1} is a **horizontal strip** of size w_i. The horizontal strip
condition on the increments c[r] = α^i[r] − α^{i−1}[r] is:

    c[r] ≥ 0
    c[r] ≤ λ[r] − α^{i−1}[r]        (stay inside λ)
    c[r+1] ≤ α^{i−1}[r] − α^{i−1}[r+1]   (at most one box per column)
    Σ_r c[r] = w_i

DP state at level i: a map { partition α ↦ count of paths from μ to α }. When no flag
bounds are active, w is sorted in decreasing order before the DP to minimise the peak
number of intermediate states. Transitions are merged directly into the next-level map
to avoid peak intermediate allocations.

### Ehrhart polynomial computation

1. Compute the degree d by constraint propagation on GT(λ/μ, w).
2. Evaluate K(nλ/nμ, nw) for n = 1, ..., d+1 in parallel (independent DP calls).
3. Solve the (d+1) × (d+1) Vandermonde system over Q to recover polynomial coefficients.

### h*-vector

Express the Ehrhart polynomial in the binomial basis {C(n+d, d), C(n+d−1, d), ..., C(n, d)}.
The coefficients are h*_0, ..., h*_d by definition of the Ehrhart series.

---

## Installation

```
cd kostka
cargo build --release
```

Binary: `target/release/kostka`.

---

## CLI reference

### Subcommands

| Subcommand | Computes |
|---|---|
| `kostka`   | K(λ, w) — single Kostka coefficient |
| `skew`     | K(λ/μ, w) — single skew Kostka coefficient |
| `degree`   | Degree of the Ehrhart polynomial (= dim of GT polytope) |
| `ehrhart`  | Ehrhart polynomial n ↦ K(nλ/nμ, nw), with degree and values |
| `hstar`    | h*-vector of GT(λ/μ, w) |
| `syt`      | Number of SYT of shape λ via hook-length formula |
| `table`    | Batch computation; see below |

### Common flags

| Flag | Description |
|---|---|
| `--lambda <parts>` | Outer shape, e.g. `4,3,2,1` |
| `--mu <parts>` | Inner shape for skew (omit for non-skew) |
| `--weight <parts>` | Weight composition, e.g. `2,2,3,3` (alias: `-w`) |
| `--upper-flags <parts>` | Per-row upper bounds for flagged Schur (experimental) |
| `--lower-flags <parts>` | Per-row lower bounds for flagged Schur (experimental) |
| `--max-states <n>` | Abort if any DP level exceeds n states (prevents OOM) |
| `--format <fmt>` | Output format: `text` (default), `json`, `csv` |
| `--verbose` | Show propagation steps for `degree`; show DP level sizes elsewhere |

---

## Sample inputs and outputs

### Single Kostka coefficient

```
$ kostka kostka --lambda 3,2,1 --weight 2,2,2
K( 3,2,1 | 2,2,2 ) = 2
```

```
$ kostka kostka --lambda 4,2 --weight 3,2,1
K( 4,2 | 3,2,1 ) = 2
```

### Single skew Kostka coefficient

```
$ kostka skew --lambda 4,3,1 --mu 2,1 --weight 2,2,1
K( 4,3,1/2,1 | 2,2,1 ) = 8
```

```
$ kostka skew --lambda 3,2,1 --mu 1 --weight 2,2,1
K( 3,2,1/1 | 2,2,1 ) = 4
```

With row flags — label 1 restricted to row 1 only (upper_flags[0] = 1):

```
$ kostka skew --lambda 3,2,1 --mu 1 --weight 2,2,1 --upper-flags 1,3,3
K( 3,2,1/1 | 2,2,1 ) = 1
```

### Degree of the Ehrhart polynomial

```
$ kostka degree --lambda 3,2,1 --weight 2,2,2
degree: 1
```

```
$ kostka degree --lambda 4,3,2,1 --weight 2,2,3,3
degree: 3
```

Empty polytope (weight sum mismatches skew size, or constraints are infeasible):

```
$ kostka degree --lambda 3,2,1 --weight 3,3,3
degree: empty polytope
```

With flags (labels 1 and 2 restricted to row 1 only, reduces dimension from 3 to 2):

```
$ kostka degree --lambda 3,2 --weight 1,1,1,1,1 --upper-flags 1,1,2,2,2
degree: 2
```

### Ehrhart polynomial

```
$ kostka ehrhart --lambda 3,2,1 --weight 2,2,2
Ehrhart polynomial of GT( 3,2,1 | 2,2,2 ):
  degree:     1
  polynomial: n + 1
  values:     n=1 -> 2,  n=2 -> 3,  n=3 -> 4,  n=4 -> 5,  n=5 -> 6
```

```
$ kostka ehrhart --lambda 4,3,2,1 --weight 2,2,3,3
Ehrhart polynomial of GT( 4,3,2,1 | 2,2,3,3 ):
  degree:     3
  polynomial: (1/3!) * (n^3 + 6n^2 + 11n + 6)
  values:     n=1 -> 4,  n=2 -> 10,  n=3 -> 20,  n=4 -> 35,  n=5 -> 56
```

```
$ kostka ehrhart --lambda 4,3,2,1 --mu 2,1 --weight 2,2,2,1
Ehrhart polynomial of GT( 4,3,2,1/2,1 | 2,2,2,1 ):
  degree:     3
  polynomial: (1/3!) * (8200n^3 - 41616n^2 + 70016n - 36396)
  values:     n=1 -> 34,  n=2 -> 462,  n=3 -> 3418,  n=4 -> 17102,  n=5 -> 49714
```

### h*-vector

```
$ kostka hstar --lambda 6,4 --weight 2,2,2,2,2
h*-vector of GT( 6,4 | 2,2,2,2,2 ):
  degree:      3
  h* =         [1, 11, 11, 1]
  palindromic: yes
  unimodal:    yes
```

```
$ kostka hstar --lambda 3,2,1 --weight 2,2,2
h*-vector of GT( 3,2,1 | 2,2,2 ):
  degree:      1
  h* =         [1, 0]
  palindromic: no
  unimodal:    yes
```

### Standard Young tableaux

```
$ kostka syt --lambda 3,2,1 --hooks
f^(3,2,1) = 16
hook lengths: 5 3 1 / 3 1 / 1
```

```
$ kostka syt --lambda 5,4,3,2,1
f^(5,4,3,2,1) = 292864
```

---

## Batch / table mode

The `table` subcommand computes multiple values at once. The intent is always one of:

- Fix λ (and optionally μ), range over all weights w
- Fix w, range over all shapes λ (and optionally μ) of a given size
- Fix n, compute the full Kostka matrix for all λ, w ⊢ n

### Full Kostka matrix for a given size

Shapes λ and weights w (treated as partitions for the columns) are both listed in dominance
order. This is the classical Kostka matrix.

```
$ kostka table --n 4
Kostka matrix for n=4

                        (4)     (3,1)     (2,2)   (2,1,1) (1,1,1,1)
            (4)           1         1         1         1         1
          (3,1)           0         1         1         2         3
          (2,2)           0         0         1         1         2
        (2,1,1)           0         0         0         1         3
      (1,1,1,1)           0         0         0         0         1
```

### All weights for a fixed shape

Lists every composition w of |λ| (up to a maximum alphabet size), computing K(λ, w).
Compositions are grouped by their underlying partition (sorted w), then listed within
each group to make the w-symmetry visible.

```
$ kostka table --lambda 3,2 --all-weights
K( 3,2 | w ) for all weights w of size 5:

  w = (1,1,1,1,1           )  ->  5
  w = (1,1,1,2             )  ->  3
  w = (1,1,2,1             )  ->  3
  ...
  w = (5                   )  ->  0
```

Use `--weight-as-partition` to deduplicate (list only partitions w):

```
$ kostka table --lambda 3,2 --all-weights --weight-as-partition
K( 3,2 | w ) for all weights w of size 5:

  w = (5                   )  ->  0
  w = (4,1                 )  ->  0
  w = (3,2                 )  ->  1
  w = (3,1,1               )  ->  1
  w = (2,2,1               )  ->  2
  w = (2,1,1,1             )  ->  3
  w = (1,1,1,1,1           )  ->  5
```

### All Ehrhart polynomials for a fixed shape

Compute the polynomial n ↦ K(nλ, nw) for every weight w (as a partition), outputting one
polynomial per row.

```
$ kostka table --lambda 3,2 --all-weights --ehrhart --weight-as-partition
Ehrhart polynomials for GT( 3,2 | w )  (|w| = 5):

  w = (5                   )  deg=0  poly= 0
  w = (4,1                 )  deg=0  poly= 0
  w = (3,2                 )  deg=0  poly= 1
  w = (3,1,1               )  deg=0  poly= 1
  w = (2,2,1               )  deg=1  poly= n + 1
  w = (2,1,1,1             )  deg=2  poly= (1/2)n^2 + (3/2)n + 1
  w = (1,1,1,1,1           )  deg=3  poly= (1/2)n^3 + (3/2)n^2 + 2n + 1
```

Use `--hstar` to append the h*-vector to each row:

```
$ kostka table --lambda 3,2 --all-weights --ehrhart --hstar --weight-as-partition
  w = (2,2,1               )  deg=1  poly= n + 1                           h*=[1, 0]
  w = (2,1,1,1             )  deg=2  poly= (1/2)n^2 + (3/2)n + 1           h*=[1, 0, 0]
  w = (1,1,1,1,1           )  deg=3  poly= (1/2)n^3 + (3/2)n^2 + 2n + 1   h*=[1, 1, 1, 0]
  ...
```

### All Ehrhart polynomials for a fixed skew shape

```
$ kostka table --lambda 4,3,2,1 --mu 2,1 --all-weights --ehrhart --weight-as-partition
Ehrhart polynomials for GT( 4,3,2,1/2,1 | w )  (|w| = 7):

  w = (7                   )  deg=0  poly= 0
  w = (6,1                 )  deg=0  poly= 0
  w = (5,2                 )  deg=0  poly= 0
  w = (5,1,1               )  deg=0  poly= 0
  w = (4,3                 )  deg=0  poly= 1
  w = (4,2,1               )  deg=2  poly= (1/2)n^2 + (3/2)n + 1
  ...
```

### Filter and sort options for table mode

```
--min-degree <d>     Only show rows where deg ≥ d
--max-degree <d>     Only show rows where deg ≤ d
--palindromic        Only show rows where h* is palindromic
--unimodal           Only show rows where h* is unimodal
--sort-by <key>      Sort output by: weight (default), degree, leading-coeff, h*-lex
--limit <n>          Show at most n rows
--max-states <n>     Abort if any DP level exceeds n states (prevents OOM)
```

Example: find all weights w for shape (4,3,2,1) whose GT polytope has degree exactly 3:

```
$ kostka table --lambda 4,3,2,1 --all-weights --ehrhart --min-degree 3 --max-degree 3 --weight-as-partition
Ehrhart polynomials for GT( 4,3,2,1 | w )  (|w| = 10):

  w = (4,2,2,1,1           )  deg=3  poly= (1/6)n^3 + n^2 + (11/6)n + 1
  w = (3,3,2,2             )  deg=3  poly= (1/6)n^3 + n^2 + (11/6)n + 1
```

### Skew shape, all weights

```
$ kostka table --lambda 4,3,1 --mu 2,1 --all-weights --ehrhart --hstar --weight-as-partition
Ehrhart polynomials for GT( 4,3,1/2,1 | w )  (|w| = 5):

  w = (5                   )  deg=0  poly= 0                               h*=[0]
  w = (4,1                 )  deg=0  poly= 1                               h*=[1]
  w = (3,2                 )  deg=0  poly= 3                               h*=[3]
  w = (3,1,1               )  deg=0  poly= 5                               h*=[5]
  w = (2,2,1               )  deg=0  poly= 8                               h*=[8]
  w = (2,1,1,1             )  deg=0  poly= 14                              h*=[14]
  w = (1,1,1,1,1           )  deg=0  poly= 25                              h*=[25]
```

### Output as CSV (for downstream processing)

```
$ kostka table --lambda 3,2 --all-weights --ehrhart --hstar --weight-as-partition --format csv
5,0,0,0,""[0]""
4,1,0,0,0,""[0]""
3,2,1,0,1,""[1]""
3,1,1,1,0,1,""[1]""
2,2,1,2,1,n + 1,""[1,0]""
2,1,1,1,3,2,(1/2)n^2 + (3/2)n + 1,""[1,0,0]""
1,1,1,1,1,5,3,(1/2)n^3 + (3/2)n^2 + 2n + 1,""[1,1,1,0]""
```

### Output as JSON

```
$ kostka table --lambda 3,2 --all-weights --ehrhart --weight-as-partition --format json
{"weight":[5],"value":"0","degree":0,"polynomial":"0","hstar":null}
{"weight":[4,1],"value":"0","degree":0,"polynomial":"0","hstar":null}
{"weight":[3,2],"value":"1","degree":0,"polynomial":"1","hstar":null}
{"weight":[3,1,1],"value":"1","degree":0,"polynomial":"1","hstar":null}
{"weight":[2,2,1],"value":"2","degree":1,"polynomial":"n + 1","hstar":null}
{"weight":[2,1,1,1],"value":"3","degree":2,"polynomial":"(1/2)n^2 + (3/2)n + 1","hstar":null}
{"weight":[1,1,1,1,1],"value":"5","degree":3,"polynomial":"(1/2)n^3 + (3/2)n^2 + 2n + 1","hstar":null}
```

---

## Notes on weight ordering

K(λ/μ, w) is invariant under permutations of w for standard (non-flagged) computations,
because the Schur function s_{λ/μ} is symmetric. The tool normalizes w internally when
computing Kostka numbers or Ehrhart polynomials without flags.

However, the tool accepts w in any order and retains the ordering in all output, because:

- Row-flagged Schur functions break this symmetry: the i-th part of w indexes the i-th
  row of the GT pattern, and the flags f_i bound the entries of that row.
- Keeping the order explicit makes it straightforward to extend all subcommands to the
  flagged setting without changing the interface.

Zero parts in w specify explicit alphabet letters that contribute no boxes. They do not
affect K(λ/μ, w) but do affect the row structure for flagged computations.

---

---

## References

- Gelfand, I. M. and Tsetlin, M. L. (1950). Finite-dimensional representations of the group
  of unimodular matrices.
- Rassart, E. (2004). A polynomiality property for Littlewood-Richardson coefficients.
  *Journal of Combinatorial Theory, Series A*, 107(2), 161–179.
- Stanley, R. P. (1986). Two combinatorial applications of the Aleksandrov-Fenchel
  inequalities. *Journal of Combinatorial Theory, Series A*, 40(2), 183–208.
- Beck, M. and Robins, S. (2015). *Computing the Continuous Discretely*. Springer.
  (Ehrhart theory background.)
