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
python scripts/test_optimization.py --binary target/release/ehz-capacity
python -O scripts/test_optimization.py --math-only
```

## Optimization quality and algorithm validity

The [optimization guide](OPTIMIZATION.md) defines the implementation sequence
and acceptance criteria for new backends. `scripts/test_optimization.py` adds
10 mathematical/checker/physical-LP test methods and 5 production integration
test methods, with parameterized cases inside each method. These counts are
test methods, not a benchmark instance count or a count of formal proofs.

Use `python scripts/test_optimization.py --math-only` to run the mathematical
and checker tests without Rust. Use `--binary target/release/ehz-capacity` for
the full suite; a missing binary is an error, never a silent skip. Each binary
invocation has a 120-second timeout. Checks use unittest assertions or explicit
exceptions, so they remain active under Python `-O`.

The suite tests exact ordering against full permutations, closure-dependent
cyclic reduction, a suboptimal but feasible square witness, the failure of a
d+1 support cutoff, exact duplicate-row plateaus, planar area, free-weight R4
fixtures, and physical itinerary LP solutions with independently checked
segments. A deliberately invalid itinerary demonstrates why the experimental
LP's normalized-closure precondition must be enforced by a future backend.
Production checks cover exact quality, original-coordinate witnesses, facet
relabeling, row scaling, translation, symplectic and antisymplectic transforms,
redundancy, symmetric products, tied optima, conditioning, and domain rejection.

These extend rather than replace the existing suites. They do not certify the
floating solver or test production backends that have not been implemented.
The workflow `.github/workflows/validation.yml` runs all suites on Python 3.12
and the repository's pinned Rust toolchain. Compare performance separately on
controlled hardware; CI timeouts catch hangs, not small timing regressions.

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

