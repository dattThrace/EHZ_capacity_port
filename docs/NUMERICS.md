# Numerics and MATLAB compatibility

This page describes the numerical backend, which reports `numerical-estimate`.
The [exact backend](EXACT.md) uses separate rational arithmetic, domain checks,
and certificate semantics.

| MATLAB | Port |
| --- | --- |
| `EHZ_perms.m`, including local routines | `src/lib.rs`: `ehz_perms`, `ehz_from_facets`, `omega` |
| `barycenter.m` | `barycenter`: volume-weighted Delaunay simplex centroid |
| `calcNormalsAndFaces.m` | `calc_normals_and_faces`: cofactor unit normals, outward orientation, face merging and incidence |
| `calcHeights.m` | `calc_heights` |
| `prepareAeqAndBeq.m` | `prepare_aeq_and_beq` |
| `KtimesT.m` | `ktimes_t` (K varies fastest) |

Rust owns centroid arithmetic, hull normal reconstruction, constraint reduction,
linear solves, symplectic pairings, support/permutation enumeration, and scoring.
`scripts/qhull_adapter.py` invokes SciPy's Qhull with `QJ Qt` and returns only
Delaunay/hull simplex indices. A simplex with exactly d+1 vertices is already its
own triangulation. This avoids a second numerical kernel implementation while
retaining the Qhull geometry adapter proposed in the research structure.
`EHZ_QHULL_ADAPTER` can override the adapter path; the default is the source-tree
path compiled into the library, so keep the source checkout when moving binaries
or set this variable to the installed adapter.

The search visits every support of size 2 through m and every cyclic order with
the first facet anchored. It streams candidates and runs serially instead of
`parfor`; no mathematical pruning or special-family solver is used. The paper's
`1/(2 max Q)` and the original `-1/(2 min Q)` agree by order reversal.

Dimensionless weights `t_i=beta_i h_i` and scaled polar rows remove common length
and inequality-row scaling before solving. Equality rows are scaled and selected
by pivoted, reorthogonalized Gram-Schmidt; KKT systems use relative SVD rank tests
(`1e-12`) and residual checks (`1e-9`), replacing the original absolute determinant
and entry-zeroing cutoffs. Every candidate is checked against the full equalities.
Roundoff negatives up to `1e-12` are clamped and constraints checked again.
`singular_kkt` counts numerical KKT rejections, not certified exact nullities.

Geometry retains the volume centroid, QJ triangulations, and cofactor unit
normals. Vertex clouds are translated and scaled before geometry arithmetic;
reconstructed planes are checked against all original normalized vertices.
Facet merging uses componentwise `1e-8` tolerance rather than MATLAB's scaled
`uniquetol`; normal orientation uses a side comparison. Qhull versions, rank
decisions and tie ordering can differ from MATLAB, especially near degeneracy.
Exhaustive floating-point enumeration is not a rigorous certificate.

For subsequent optimizations, follow the [optimization guide](OPTIMIZATION.md)
and the [rational solver contract](RATIONAL_SOLVER.md). They specify dispatch
hypotheses and the evidence required for pruning and global bounds.

## Mathematical reference and review

The algorithm implements Theorem 1.1 of Pazit Haim-Kislev's
[On the Symplectic Size of Convex Polytopes](https://doi.org/10.1007/s00039-019-00486-4).
The [source-paper review](../validation/PAPER_VALIDATION.md) derives the formula,
explains support enumeration and singular systems, and records corrections to
the original port's numerical behavior.

The original source is preserved in `upstream/EHZ-capacity/` at commit
`22dd473cc8b9b2397b44cd13acbceb763dd978e9`;
[upstream provenance](../upstream/README.md) records its repository and history.
