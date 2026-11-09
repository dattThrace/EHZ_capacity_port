# Development and validation

Run these commands from the repository root after installing the dependencies
and building the release binary as described in [usage](USAGE.md).

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
[PDF](../research/sources/1712.03494.pdf). It does not import the MATLAB-style
reference. Its 31 checks cover planar area, R4/R6 simplices, nonsymmetric truncated
simplices with free weights, symmetric examples, extreme common/row scalings,
symplectic invariance, and invalid theorem-domain inputs. See
[the mathematical review and correction log](../validation/PAPER_VALIDATION.md) and
[recorded exact expected results](../validation/paper-results.json).

The additional 19-case port suite uses `scripts/reference.py`, a retained
MATLAB-style NumPy/SciPy reference, and checks original fixtures and invariants.
That legacy script intentionally retains the old cutoffs; it is not the source
of truth for the corrected solver. Five Rust integration tests also pass.
There is no MATLAB/Octave runtime comparison or general numerical certification.

The earlier research and probes remain unchanged in `research/`. They describe
later algorithm options and can be run with `research/requirements.txt`; they are
not used to accelerate this implementation. See [PORT_REPORT.md](../PORT_REPORT.md)
for the recorded validation outcome and limitations.

## Archived research probes

The research probes are separate from the production validation suites. To run
them with their recorded dependencies:

```sh
python -m pip install -r research/requirements.txt
python research/probes/check_ordering_dp.py
python research/probes/check_symmetric_product.py
python research/probes/check_itinerary_lp.py
```

The ordering and symmetric-product probes use only the Python standard library;
the itinerary LP probe also requires NumPy/SciPy. Run without `-O` or
`PYTHONOPTIMIZE`, because the probes use assertions. The
[research archive](../research/README.md) describes the original files and results.
