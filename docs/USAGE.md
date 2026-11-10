# Installation and CLI usage

Run the commands below from the repository root.

Requires Rust (pinned to 1.85.1 by `rust-toolchain.toml`). Numerical vertex input
also requires Python 3.12 with NumPy/SciPy; facet input needs no Python runtime.
The independent exact verifier uses Python's standard library. Use the committed
Cargo lockfile.

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
These legacy input forms use the numerical backend, also available explicitly
with `--mode numerical`.

For exact rational H-polytopes, supply `"schema": "ehz-exact-h-v1"` and encode
each scalar as a string or numerator/denominator pair of integer strings. The
schema selects the exact backend automatically; `--mode exact` is also accepted.
Exact mode checks boundedness and removes redundancy with rational certificates.
It supports `--max-steps N` and `--timeout-ms N`, preserving an incumbent as an
upper bound when the search stops. Add `--certify-interval` to also compute a
conservative certified interval for a partial positive witness. The numerical
backend reports `numerical-estimate`. See [the exact contract](EXACT.md) for input
examples, status meanings, verification, and limitations.

Output includes capacity, signed minimum and minimizing weights/order, the
positive maximum and reversed maximizing weights/order in the paper's convention,
closure/normalization residuals, search counts, numerical status, and (for vertex
input) geometry and centroid.
All output indices are zero based; beta is in the returned order. `--geometry`
returns only geometry. Errors go to stderr with a nonzero exit code. Enumeration
is serial and has no GUI wait bar. Exact invalid inputs also produce a JSON
`invalid-input` result on stdout. Use exact-mode limits to obtain a partial
certificate; forcibly interrupting the process can prevent JSON output.

## Qhull adapter configuration

The Python runner sets `EHZ_PYTHON` to its own interpreter. When invoking the
binary directly, set `EHZ_PYTHON` to the Python executable containing NumPy/SciPy.
`EHZ_QHULL_ADAPTER` can override the adapter script path. Its default is the
source-tree path compiled into the library; keep the source checkout when moving
binaries or set the variable to the installed adapter.

See [numerics and compatibility](NUMERICS.md) for the algorithm, tolerances,
MATLAB function mapping, and geometry adapter details.
