# EHZ capacity

A Rust library and command-line tool for computing numerical estimates of the
Ekeland–Hofer–Zehnder capacity of convex polytopes. It implements the exhaustive
combinatorial formula from Pazit Haim-Kislev's
[On the Symplectic Size of Convex Polytopes](https://doi.org/10.1007/s00039-019-00486-4)
([paper PDF](https://arxiv.org/pdf/1712.03494v3)), building on the original
[MATLAB reference implementation](https://github.com/pazithaimkislev/EHZ-capacity).

Polytopes can be supplied as vertices or facet inequalities in
`(q_1,...,q_n,p_1,...,p_n)` coordinates. Rust handles the numerical and
combinatorial calculations; Python provides a Qhull geometry adapter and
validation tools. The solver enumerates facet supports and cyclic orders, with
factorial runtime in the number of facets. Results use floating-point arithmetic
and are not certified capacities.

## Quick start

Requires Rust, with toolchain 1.85.1 selected by `rust-toolchain.toml`. Vertex
input also requires Python 3.12 with NumPy/SciPy.

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
cargo build --release --locked
python scripts/run.py examples/square.json
```

The square example has capacity 4. The runner invokes the Rust CLI and uses its
Python interpreter for geometry. For facet input, the Rust binary can run
without Python:

```sh
printf '%s\n' '{"normals":[[1,0],[-1,0],[0,1],[0,-1]],"heights":[1,1,1,1]}' | target/release/ehz-capacity
```

The CLI reads JSON from a file or stdin and writes JSON containing the capacity,
objective values, facet orders, weights, and numerical diagnostics. Vertex input
also returns the reconstructed geometry. See [usage](docs/USAGE.md) for input
requirements, output conventions, geometry-only mode, and adapter configuration.

## Codebase

| Location | Purpose |
| --- | --- |
| [src/lib.rs](src/lib.rs) | Reusable library: centroids, facet normals and heights, constraint reduction, symplectic pairings, and exhaustive capacity search |
| [src/main.rs](src/main.rs) | JSON command-line interface for vertex and facet inputs |
| [scripts/qhull_adapter.py](scripts/qhull_adapter.py) | SciPy/Qhull adapter returning triangulation and hull simplex indices to Rust |
| [scripts/run.py](scripts/run.py) | Build-and-run wrapper using the active Python interpreter |
| [scripts/](scripts/) | Fixture generation, independent exact mathematical oracles, and a legacy MATLAB-style reference |
| [tests/](tests/) | Rust integration tests |
| [examples/](examples/) | JSON inputs, including the original repository's cube × cross-polytope example |
| [validation/](validation/) | Source-paper derivation, correction log, and recorded exact expected results |
| [research/](research/) | Mathematical research, experimental probes, and archived sources |
| [upstream/](upstream/) | Pinned original MATLAB source, provenance, and Git history bundle |

The library exposes `ehz_perms` for vertices and `ehz_from_facets` for facet
inequalities. Geometry and helper operations are also available separately.
The Python geometry adapter supplies simplex indices; centroid arithmetic,
normal reconstruction, linear solves, and the capacity search run in Rust.

## Documentation

- [Usage](docs/USAGE.md): installation, CLI examples, input/output, and Python configuration.
- [Numerics and compatibility](docs/NUMERICS.md): algorithm, tolerances, and mapping from the MATLAB functions.
- [Development and validation](docs/DEVELOPMENT.md): checks, fixtures, reference calculations, and research probes.
- [Mathematical review](validation/PAPER_VALIDATION.md): correspondence with the paper and numerical limitations.
- [Port report](PORT_REPORT.md): implementation history, corrections, and recorded validation outcomes.
- [Research memo](research/RESEARCH.md): mathematical background and proposed algorithms beyond the exhaustive solver.

## Attribution

The mathematical algorithm and original MATLAB implementation are by Pazit
Haim-Kislev. See [upstream provenance](upstream/README.md) for the preserved source
revision and history. The preserved upstream has no blanket license notice; no
new license is assigned to upstream code, papers, or these collected materials.
