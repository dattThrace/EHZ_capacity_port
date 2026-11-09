# EHZ capacity: Rust port

A Rust library and JSON command-line port of Pazit Haim-Kislev's MATLAB
implementation, with Python geometry, reference, and validation scripts.
Input vertices use `(q_1,...,q_n,p_1,...,p_n)` coordinates.
The original source is preserved in `upstream/EHZ-capacity/` at commit
`22dd473cc8b9b2397b44cd13acbceb763dd978e9`.

The solver retains the exhaustive algorithm. A source-paper review corrected
scale-dependent numerical cutoffs and added theorem-domain checks. It does not implement the optimization proposals in
[RESEARCH.md](research/RESEARCH.md). Results are numerical estimates, not exact
or certified capacities. Runtime grows factorially with the number of facets.

## Build and run

Requires Rust (pinned to 1.85.1 by `rust-toolchain.toml`) and Python 3.12 with
NumPy/SciPy. Use the committed Cargo lockfile.

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
cargo build --release --locked
python scripts/run.py examples/square.json
python scripts/run.py examples/readme_product.json
```

The square and original README product both have capacity 4. The wrapper sets
`EHZ_PYTHON` to its Python interpreter. Alternatively run the binary directly:

```sh
EHZ_PYTHON=.venv/bin/python target/release/ehz-capacity examples/square.json
target/release/ehz-capacity --geometry examples/square.json
printf '%s\n' '{"normals":[[1,0],[-1,0],[0,1],[0,-1]],"heights":[1,1,1,1]}' | target/release/ehz-capacity
```

Input is either `{"vertices": [[...], ...]}` or
`{"normals": [[...], ...], "heights": [...]}`. A filename or `-`/stdin is
accepted. Facet input describes `normal_i.x <= height_i`, with a strictly
interior origin and positive heights. It bypasses geometry and requires no
Python installation. Supply a full-dimensional bounded polytope. The facet
interface checks dimensions, finite values, positive heights, nonzero normals,
and numerical positive spanning (boundedness). Ill-conditioned inputs can be
rejected; this is not a certified geometry check. Redundant rows are not removed.

Output includes capacity, signed minimum and minimizing weights/order, the
positive maximum and reversed maximizing weights/order in the paper's convention,
closure/normalization residuals, search counts, numerical status, and (for vertex
input) geometry and centroid.
All output indices are zero based; beta is in the returned order. `--geometry`
returns only geometry. Errors go to stderr with a nonzero exit code. Enumeration
is serial and has no GUI wait bar; interrupt the process to cancel a large job.

## Source mapping and compatibility

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

## Validation and scripts

```sh
cargo test --locked
python scripts/generate_fixtures.py
python scripts/validate.py
python scripts/validate_paper.py
python scripts/reference.py examples/triangle.json
cargo fmt --check
cargo clippy --locked --all-targets -- -D warnings
```

The independent source-paper suite uses exact rational maximization of
Theorem 1.1 and the separate symmetric formula in Corollary 1.2 from the archived
[PDF](research/sources/1712.03494.pdf). It does not import the MATLAB-style
reference. Its 31 checks cover planar area, R4/R6 simplices, nonsymmetric truncated
simplices with free weights, symmetric examples, extreme common/row scalings,
symplectic invariance, and invalid theorem-domain inputs. See
[the mathematical review and correction log](validation/PAPER_VALIDATION.md) and
[recorded exact expected results](validation/paper-results.json).

The additional 19-case port suite uses `scripts/reference.py`, a retained
MATLAB-style NumPy/SciPy reference, and checks original fixtures and invariants.
That legacy script intentionally retains the old cutoffs; it is not the source
of truth for the corrected solver. Five Rust integration tests also pass.
There is no MATLAB/Octave runtime comparison or general numerical certification.

The earlier research and probes remain unchanged in `research/`. They describe
later algorithm options and can be run with `research/requirements.txt`; they are
not used to accelerate this implementation. See [PORT_REPORT.md](PORT_REPORT.md)
for the recorded validation outcome and limitations.

## Attribution

The mathematical algorithm and original MATLAB implementation are by Pazit
Haim-Kislev; see [upstream provenance](upstream/README.md) and the original README.
The preserved upstream has no blanket license notice; no new license is assigned
to upstream code, papers, or these collected materials.
