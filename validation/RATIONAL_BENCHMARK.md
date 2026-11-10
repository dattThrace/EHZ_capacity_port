# Rational solver baseline benchmark

Measured at `2026-11-10T23:10:49.497965+00:00` on `Linux x86_64` with Python `3.12.14`.
Binary path, SHA256, and the pinned toolchain specification are recorded in
[`rational-results.json`](rational-results.json).

Reproduce with `cargo build --release --locked` followed by
`python3 scripts/validate_exact.py --benchmark`. Every expected value is compared
as an exact fraction. The independent paper oracle uses nullspace coordinates and
full permutations; Rust uses KKT systems and anchored cyclic orders. Planar
shoelace areas and the symmetric corollary supply additional independent checks.

These are correctness baseline measurements, with no pruning or speedup target.
Runtime and memory include one fresh Rust CLI process per case. Launch RSS comes
from Linux `wait4`, independently of the solver, and includes the pre-exec
process-launch footprint. Solver RSS is its `/proc/self/status` VmHWM near
search completion; these two measurements have different scopes.
Intermediate bit lengths are
the solver's observed elimination/search arithmetic statistic; see the machine
report for complete counters, CPU times, fractions, and input bit lengths.
The arithmetic statistic excludes transient BigInt implementation internals
during rational normalization. All emitted certificates are also independently
verified, including input identity, row transformations, domain, and global
coverage re-enumeration. The binary uses the pinned Rust 1.85.1 release profile.

| Case | d | Input/retained rows | Input bits | Intermediate bits | Supports | Orders | Time (s) | Launch/solver RSS (KiB) | Status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| square | 2 | 4/4 | 1 | 5 | 11 | 16 | 0.034104 | 12672/1900 | exact-complete |
| decimal_square | 2 | 4/4 | 4 | 8 | 11 | 16 | 0.007968 | 12672/1856 | exact-complete |
| asymmetric_quadrilateral | 2 | 4/4 | 4 | 22 | 11 | 10 | 0.004589 | 12672/1872 | exact-complete |
| seeded_rational_pentagon_0 | 2 | 5/5 | 12 | 103 | 26 | 64 | 0.038385 | 12672/1844 | exact-complete |
| seeded_rational_pentagon_1 | 2 | 5/5 | 11 | 114 | 26 | 64 | 0.108536 | 12800/1864 | exact-complete |
| seeded_rational_pentagon_2 | 2 | 5/5 | 11 | 119 | 26 | 64 | 0.047602 | 12800/1848 | exact-complete |
| seeded_rational_pentagon_3 | 2 | 5/5 | 11 | 151 | 26 | 64 | 0.048009 | 12800/1864 | exact-complete |
| simplex_r4 | 4 | 5/5 | 1 | 5 | 26 | 24 | 0.006472 | 12800/1840 | exact-complete |
| truncated_simplex_r4_1 | 4 | 6/6 | 1 | 13 | 57 | 236 | 0.213321 | 12800/1880 | exact-complete |
| truncated_simplex_r4_3_2 | 4 | 6/6 | 2 | 17 | 57 | 236 | 0.226141 | 12800/1872 | exact-complete |
| symplectic_shear_r4 | 4 | 5/5 | 8 | 20 | 26 | 24 | 0.011280 | 12800/1836 | exact-complete |
| duplicate_square | 2 | 6/4 | 2 | 6 | 11 | 16 | 0.005810 | 12800/1844 | exact-complete |
| redundant_square | 2 | 9/4 | 2 | 6 | 11 | 16 | 0.007727 | 12800/1872 | exact-complete |
| large_row_scales | 2 | 4/4 | 181 | 302 | 11 | 16 | 0.004882 | 12800/1828 | exact-complete |
| perturb_square_k0_-1 | 2 | 4/4 | 2 | 5 | 11 | 16 | 0.009327 | 12800/1828 | exact-complete |
| perturb_square_k0_+0 | 2 | 4/4 | 1 | 5 | 11 | 16 | 0.004699 | 12800/1816 | exact-complete |
| perturb_square_k0_+1 | 2 | 4/4 | 2 | 5 | 11 | 16 | 0.003349 | 12800/1828 | exact-complete |
| perturb_square_k20_-1 | 2 | 4/4 | 22 | 44 | 11 | 16 | 0.006476 | 12800/1828 | exact-complete |
| perturb_square_k20_+0 | 2 | 4/4 | 1 | 5 | 11 | 16 | 0.005990 | 12800/1836 | exact-complete |
| perturb_square_k20_+1 | 2 | 4/4 | 22 | 45 | 11 | 16 | 0.004451 | 12800/1824 | exact-complete |
| perturb_square_k60_-1 | 2 | 4/4 | 62 | 124 | 11 | 16 | 0.004619 | 12800/1836 | exact-complete |
| perturb_square_k60_+0 | 2 | 4/4 | 1 | 5 | 11 | 16 | 0.003732 | 12928/1820 | exact-complete |
| perturb_square_k60_+1 | 2 | 4/4 | 62 | 125 | 11 | 16 | 0.006416 | 12928/1832 | exact-complete |
| perturb_square_k120_-1 | 2 | 4/4 | 122 | 244 | 11 | 16 | 0.009941 | 12928/1876 | exact-complete |
| perturb_square_k120_+0 | 2 | 4/4 | 1 | 5 | 11 | 16 | 0.007102 | 12928/1840 | exact-complete |
| perturb_square_k120_+1 | 2 | 4/4 | 122 | 245 | 11 | 16 | 0.004185 | 12928/1828 | exact-complete |
| symmetric_product | 2 | 4/4 | 2 | 8 | 11 | 16 | 0.006231 | 12928/1816 | exact-complete |

The perturbation family holds facet count and dimension fixed while increasing
denominator bit length. All heights remain strictly positive, including the
negative perturbation cases. Values at 60 and 120 bits differ even when binary64
would identify the inputs. Large row scales preserve exactly the same square.

Budget exhaustion, zero-timeout behavior, preservation of a partial exact witness,
malformed fractions, zero denominators, JSON floats/numbers, and invalid geometric
domains are checked by the same command. Certificate mutation and independent
coverage re-enumeration are checked by `python3 scripts/test_verify_exact.py`.

Lifecycle checks also record incomplete domain checks, interrupted feasible
witnesses, and conservative certified intervals in the machine report. Partial
witnesses remain capacity upper bounds and never produce exact-complete status.
