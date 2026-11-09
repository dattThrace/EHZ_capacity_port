# Rust port and source-paper review

Completed 2026-11-09. The Rust library and JSON CLI cover all six preserved
MATLAB files. Python supplies a Qhull simplex-index adapter, runner, fixture
generator, legacy floating reference and independent exact mathematical oracles.
Source: Pazit Haim-Kislev's *On the symplectic size of convex polytopes*, the
archived arXiv:1712.03494v3 PDF of the paper cited by upstream.

## Corrections after independent mathematical review

- Replaced absolute determinant/entry cutoffs with dimensionless weights, row
  scaling, relative SVD rank checks and full-constraint residual checks. The
  original solver failed on a square of height `1e-6` and on equivalent facet
  rows scaled by `1e-4`; these now return `4e-12` and `4`, respectively.
- Added numerical positive-spanning/boundedness validation. An unbounded input
  with full-rank normals and a closed loop in just one symplectic plane could
  previously return a finite number; it is now rejected.
- Normalize vertex geometry and check reconstructed supporting planes against
  all original vertices. The volume-weighted centroid is retained.
- Added maximizing order/weights, maximum objective, and witness residuals in the
  paper's convention, preserving the equivalent original minimum convention.

The sign and factor of two were already correct. Reversing the order changes
Q's sign, so `-1/(2 min Q) = 1/(2 max Q)`. No research optimization was introduced
in the production search. Symmetry and exact rational solvers are validation-only.

The full derivation, hypotheses, exact-singularity argument, numerical changes,
and limitations are recorded in [PAPER_VALIDATION.md](validation/PAPER_VALIDATION.md).
Upstream and pre-existing research artifacts remain unchanged.

## Validation

Passed release build, five Rust tests, formatting, strict Clippy, the 19-case
port regression suite, and 31 independent source-paper checks. Exact expected
fractions and numerical results are in
[paper-results.json](validation/paper-results.json).

The source-paper checks use exact rational Theorem 1.1 maximization by RREF and
free-coordinate quadratic differentiation, and a second oracle from Corollary
1.2. They do not use the MATLAB repository or the legacy Python translation as
reference. Representative results:

| Body | Exact expected capacity |
| --- | ---: |
| Original README cube × cross-polytope | 4 |
| Standard R4 simplex | 1/4 |
| Translated/scaled R4 simplex `x_i>=-1, sum x_i<=1` | 25/4 |
| Same construction in R6 | 49/6 |
| Truncated R4 simplex, added `x1+x2<=1` | 21/4 |
| Asymmetric planar pentagon | 17 |

## Remaining limits

The production solver is floating point, not a certified or exact capacity
algorithm. Nearly singular systems can still be skipped and geometry below the
stated tolerances can remain unresolved. Positive spanning is checked numerically.
Results are numerical estimates; a rounded witness does not prove global
optimality. Runtime remains factorial in facet count.

Vertex input requires Python/SciPy; facet input runs entirely in Rust. MATLAB
was not available, and no published-journal PDF comparison was performed; the
review identifies and hashes the actual arXiv PDF used. Toolchain and dependency
versions are pinned. Reproduction commands and usage are in README.md.
