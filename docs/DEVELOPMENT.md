# Development and validation

Run these commands from the repository root after installing the dependencies
and building the release binary as described in [usage](USAGE.md).

```sh
cargo test --locked
cargo build --release --locked
python3 scripts/validate_exact.py --benchmark
python3 scripts/test_verify_exact.py
python scripts/test_optimization.py --binary target/release/ehz-capacity
python -O scripts/test_optimization.py --math-only
python scripts/generate_fixtures.py
python scripts/validate.py
python scripts/validate_paper.py
python scripts/reference.py examples/triangle.json
cargo fmt --check
cargo clippy --locked --all-targets -- -D warnings
```

The exact backend's acceptance suite compares canonical rational values with
`fractions.Fraction`; it never converts an expected fraction to binary64.
`scripts/validate_exact.py` checks independently derived planar areas, seeded
nonsymmetric rational free-weight examples, R4 simplex truncations, symplectic
shears, row duplicates and certified redundancies, exact decimal input,
equivalent large row scales, and perturbations smaller than binary64 resolution.
It also checks original-coordinate witnesses, cyclic/reversal conventions,
independent certificate verification, interruption semantics, and malformed or
geometrically invalid inputs. The deterministic input subset is in
[`examples/exact/`](../examples/exact/); the parameter families and seed are
defined in the validation script.

With `--benchmark`, the script writes a
[baseline report](../validation/RATIONAL_BENCHMARK.md) and
[machine-readable results](../validation/rational-results.json), recording
dimension, original and retained row counts, input and intermediate rational
bit lengths, support/order counts, elapsed and CPU time, per-process peak RSS,
and completion status. Timings are measurements for the exact exhaustive
baseline, without an acceleration threshold. Rust integration tests cover
the public exact API and structured CLI rejection. Independent verifier tests
deliberately corrupt witnesses, row mappings, and completion claims.

The linked guide branch's `scripts/test_optimization.py` retains its independent
ordering, physical-itinerary LP, and numerical production checks, extended with
exact production comparisons. Run `--math-only` without Rust or `--binary` for
all backends; a missing binary is an error. Checks remain active under Python
`-O`, and each binary invocation has a timeout. These tests do not certify the
numerical backend. The workflow `.github/workflows/validation.yml` runs the
Rust checks and all validation suites on Python 3.12 and Rust 1.85.1.

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
of truth for the corrected solver. Rust integration suites also cover the
numerical helpers and exact public API behavior.
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
