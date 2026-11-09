# Independent validation against the source paper

Reviewed 2026-11-09. Mathematical reference: Pazit Haim-Kislev,
*On the symplectic size of convex polytopes*, arXiv:1712.03494v3 (2019),
the preprint of the GAFA paper cited by the upstream README.
The actual archived PDF was read via fresh `pdftotext -layout` extraction;
this review does not use the MATLAB implementation as its specification.

- Source PDF: `research/sources/1712.03494.pdf`.
- URL: <https://arxiv.org/pdf/1712.03494v3>.
- SHA-256: `37a977d8039d6fc7e772dba4dfb046f102d1693846a9f7f3bc8b318d68e5b5d5`.
- Formula (1), Theorem 1.1: printed page 2.
- Symmetric alternative, Corollary 1.2: printed page 2.
- Action definition and normalization: printed page 5.
- Formula proof and Remarks 3.8–3.9: printed pages 14–15.
- Symplectic invariance, Remark 3.10: printed page 16.

## Formula, coordinates, and witnesses

For unit outer facet normals n_i, heights h_i, and

    M = {beta >= 0 : sum beta_i n_i = 0, sum beta_i h_i = 1},
    Q(order,beta) = sum_{j<i} beta_order(i) beta_order(j)
                                  omega(n_order(i),n_order(j)),

the paper states `capacity = 1/(2 max Q)`. Here
`omega(a,b) = a_q.b_p - a_p.b_q`, in `(q_1,...,q_n,p_1,...,p_n)` order.
The symmetric lower-triangle matrix H in the implementation satisfies
`Q = beta^T H beta / 2`; the factors of two agree.

The port searches for `min Q` and reports `-1/(2 min Q)`. This is equivalent:
reversing the order changes the sign of every pair and preserves feasibility.
The output now exposes `maximum`, `maximizing_order`, and `maximizing_beta`,
in addition to the original minimum convention. The validation script evaluates
both printed witnesses directly with the paper's formula, checks nonnegativity,
closure and normalization, and verifies the inverse-objective capacity.
A maximizing witness refers to this finite-dimensional formula; it is not an
explicit closed characteristic on the polytope boundary.

Non-unit H-rows are handled by substitution, not an additional theorem:
if `a_i = s_i n_i`, `h_i = s_i h_unit_i`, `s_i>0`, then
`beta_unit_i = s_i beta_row_i`. Both constraints and Q are unchanged.
The implementation uses dimensionless variables `t_i = beta_row_i h_i`,
`P_i = L a_i/h_i`, with `L=min(h_i/||a_i||)`. Thus

    t_i >= 0, sum t_i = 1, sum t_i P_i = 0,
    Q_original = Q_P / L^2.

This algebra removes arbitrary row scaling and common geometric scale before
factorization. Weights in the output are converted back to the supplied rows.

## Why exhaustive support/stationary enumeration implements Formula (1)

After translating to an interior origin all h_i are positive, so M is compact.
For each order a quadratic function attains its extremum on M. Restrict it to
the positive support of an extremizer. Independent equality rows give the
stationary equations `H beta + A^T lambda = 0`, `A beta = b`; no assumption of
convexity of the objective is made. Boundary extrema occur on smaller supports,
which are also enumerated. Feasible saddle candidates do not invalidate the
search: they cannot improve past the true global extremum.

For an *exactly* singular stationary KKT block with independent equality rows,
a nonzero null vector has a nonzero beta component d, `A d=0`, and
`H d + A^T mu=0`. Hence `d^T H d=0`, and stationarity also gives
`beta^T H d=0`. Q is constant along the feasible line beta+t*d. Compactness
allows moving to a coordinate becoming zero and therefore a smaller support.
This justifies omission of exactly singular blocks in exact arithmetic. A
floating-point rank decision is not such a proof; see the limits below.

Anchoring the first facet of a support loses no values: moving a weighted normal
from first to last changes Q by its pairing with the sum of all other weighted
normals, which is its negative by closure, so that change is zero. Omitted
zero-weight facets can be inserted anywhere without changing Q. A single active
facet is impossible for nonzero normals, positive height normalization, and
closure. Thus supports of size 2 through m and anchored orders cover Formula (1).
No transition-graph or other research reduction is used in production.

## Inconsistencies found and corrected

| Issue | Before | Correction |
| --- | --- | --- |
| Violated quadratic homogeneity (Remark 3.9) | Square heights `1e-6`: no negative feasible score instead of `4e-12` | Dimensionless weights/geometry, objective and equality scaling, relative SVD solves instead of absolute determinant threshold |
| Depended on inequality row scaling | Multiplying all square rows/heights by `1e-4` caused failure instead of capacity 4 | Cancel row scaling algebraically before rank/solve decisions |
| Accepted input outside theorem domain | Unbounded full-rank H-input could be scored as capacity 4 (also reproduced in the legacy Python reference) | Validate nonzero normals and numerical positive spanning/boundedness before search |
| Reduced constraints alone were trusted | No check of the full original equalities after reduction | Check the full dimensionless equalities and nonnegative weights after every solve; only roundoff negatives are clamped before rechecking |
| Joggled facets were assumed to support original vertices | QJ simplices were accepted without containment check | Normalize the vertex cloud and check reconstructed planes against all original normalized vertices and interior origin |
| Paper uses maximum, output showed only minimum | Numerically correct but hard to compare signs directly | Preserve old fields and add reversed maximizing witness and residual fields |

The positive-spanning check is numerical: with positive heights the origin is
strictly interior, and boundedness is equivalent to full linear span plus a
strictly positive relation among the normals. Vertices of
`{t>=0: sum t_i P_i=0, sum t_i=1}` have support at most d+1. Finding feasible
vertices whose union includes every row establishes such a positive relation
by averaging. This input check adds no pruning to the capacity search.

## Independent executable evidence

`scripts/paper_oracle.py` is an exact rational implementation of Formula (1),
using RREF and free-coordinate differentiation of
`Q(beta_0+Bz)`. It does not build or solve the production KKT block. It uses
full permutations (including cyclic duplicates), exact consistency decisions,
and exact nonnegative constraints, and does not import MATLAB, the legacy
reference script, or the research probes.

Its second oracle implements the distinct Corollary 1.2 formula: choose one
normal per opposite pair, enumerate signs and orders, enforce
`sum |beta_i| h_i = 1/2`, and return `1/(4 max Q)`. This is validation-only:
the production solver still enumerates the full facets without symmetry dispatch.

`scripts/validate_paper.py` verifies the PDF hash, runs those exact oracles,
checks direct planar area, and tests the production binary. Results and exact
expected fractions are stored in `validation/paper-results.json`.

Coverage includes asymmetric planar polygons with independently derived edge
inequalities; rational simplices in R4 and R6; nonsymmetric truncated R4 simplices
with free weights; a box, the upstream README product, and a skew centrally
symmetric body via Corollary 1.2; scaling through `1e-9` and `1e9` in both vertex
and facet input; unequal facet row scalings; a rational R4 symplectic shear; and
rejection of unbounded and zero-normal inputs. The translated simplex
`x_i >= -1, sum x_i <= 1` has capacity `(d+1)^2/d`, not `1/d`: it is a scale
`d+1` copy of the standard simplex. Both tested simplex conventions agree.

Reproduce from the repository root:

```sh
cargo build --release --locked
cargo test --locked
python scripts/validate_paper.py
python scripts/validate.py
cargo fmt --check
cargo clippy --locked --all-targets -- -D warnings
```

## What is and is not established

The formula, sign, normalization, support coverage, and exact-singularity
argument agree with the cited PDF. The exact oracles establish the rational
expected answers for the tested inputs; the production results agree to the
stated numerical tolerances. These checks extend beyond MATLAB compatibility.

The production code remains floating point, with relative rank tolerance
`1e-12`, full-equality residual tolerance `1e-9`, and tolerance-based Qhull facet
reconstruction/merging. Nearly singular nonzero systems can still be skipped;
`singular_kkt` counts numerically rejected KKT solves, not proved exact nullities.
Input validation can reject ill-conditioned bounded bodies, and geometric
features below tolerance need not be resolved. These are not exact arithmetic,
interval bounds, or a certified exhaustive computation. Remark 3.8 makes a
truly feasible positive-Q witness an upper bound on capacity; a rounded witness
alone does not certify global optimality. No general certification claim is made.

No MATLAB/Octave execution or comparison with the published journal PDF was
performed. The preserved arXiv PDF is the explicitly identified mathematical
source for this validation. Original upstream/research artifacts remain unchanged.
