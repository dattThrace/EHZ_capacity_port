# Exact rational H-polytopes

The exact backend implements the exhaustive Haim-Kislev formula using
arbitrary-precision integers and canonical rational numbers. It accepts bounded,
full-dimensional polytopes in a positive even dimension, with rational nonzero
facet normals and strictly positive heights about the supplied interior origin.
It does not recenter input or convert vertices to exact facets.

## Input and execution

Build with the pinned toolchain and lockfile:

```sh
cargo build --release --locked
```

This input is the square of capacity 4:

```json
{
  "schema": "ehz-exact-h-v1",
  "normals": [["1", "0"], ["-1", "0"], ["0", "1"], ["0", "-1"]],
  "heights": ["1", "1", "1", "1"]
}
```

Save it as `square-exact.json`, then run:

```sh
target/release/ehz-capacity square-exact.json > square-certificate.json
python scripts/verify_exact.py square-certificate.json --input square-exact.json
```

The schema selects exact mode automatically; `--mode exact` is also accepted.
The schema is required even with that flag. Legacy input without a schema keeps
the numerical backend, which can also be selected with `--mode numerical`.
An exact schema combined with `--mode numerical` is an error. `--geometry` is
available only for numerical vertex input.

Every exact scalar must be a string or a two-element pair of integer strings.
For example, `"3"`, `"-2/7"`, `"0.125"`, and `["1", "8"]` are valid.
Decimal strings denote their exact decimal value, so `"0.1"` means `1/10`.
JSON numbers, malformed fractions, and zero denominators are rejected. Exact
input is never converted through `f64`. To express an existing machine float,
the caller must supply its exact binary rational value explicitly; the parser
does not infer an intended fraction. Output scalars are canonical strings:
reduced fractions with a positive denominator, or integers when it is one.
Coordinates use `(q_1,...,q_n,p_1,...,p_n)` and indices are zero based.

## Preprocessing and exhaustive search

Facet rows retain their rational meaning. Positive rational scales clear
denominators and divide common factors; Euclidean unit normalization is never
used. The certificate records the scale and representative of every original
row. Positive multiples of identical inequalities are deduplicated. Other
redundant inequalities are removed only with nonnegative rational multiplier
certificates, retaining original row identities throughout.

Boundedness is established by exact full column rank and a strictly positive
rational vector in the kernel of the transposed original normal matrix. The
search enumerates all supports and cyclic orders anchored at one row. It caches
the equality reduction for each support, solves unique weights or nonsingular
stationary systems exactly, checks feasibility against the original equalities,
and independently rescores each candidate. Exactly singular stationary systems
are omitted because their flat directions reach a smaller support, all of which
are included in the enumeration. There is no dimension-plus-one support cutoff.
The quadratic can be indefinite; no convexity or floating eigenvalue test is
used.

## Results and limits

JSON output has a `result` object with schema `ehz-exact-result-v1`. It includes
the canonical original input, row transformations and domain certificates, a
feasible maximizing witness when available, completion information, and search
statistics. Witness order indices identify original input rows; witness weights
are in the original row coordinates. For a witness with positive objective
`Q`, `capacity_upper_bound` is exactly `1/(2Q)`.

| Status | Meaning |
| --- | --- |
| `exact-complete` | The exact search completed all justified enumeration; `capacity` is the resulting exact fraction. |
| `exact-witness/incomplete` | A positive exact witness was retained, but the search stopped before completion. Only its capacity upper bound is established. |
| `certified-interval` | A partial positive witness and a conservative objective envelope establish an exact capacity interval; enumeration remains incomplete. |
| `incomplete` | A limit stopped domain checking or search before a positive witness was available. This does not establish invalidity. |
| `invalid-input` | The input violates the parsing or geometric domain requirements. The CLI also exits nonzero. |

The numerical backend reports `numerical-estimate`. With independently verified
bounds `0 < L <= Q_* <= U`, the capacity lies in `[1/(2U), 1/(2L)]`.
The optional `--certify-interval` flag adds a conservative bound to a partial
positive witness. It computes
`M = max_ij |omega(a_i/h_i, a_j/h_j)|` in original input coordinates. Writing
`t_i = h_i beta_i`, nonnegative feasible weights have `sum t_i = 1`, so
`Q <= M sum_{j<i} t_i t_j <= M/2`. The resulting interval is
`[1/M, 1/(2Q_witness)]`. The result records `interval.lower`, `interval.upper`,
`interval.objective_upper_bound`, and method `pairing-envelope-v1`, while
`capacity` stays null and `completion.complete` stays false. This bound can be
wide. The flag leaves a complete result or a result without a witness unchanged;
the default partial-result status remains `exact-witness/incomplete`.

Use deterministic work limits or an elapsed-time limit for larger inputs:

```sh
target/release/ehz-capacity --max-steps 10000 square-exact.json
target/release/ehz-capacity --timeout-ms 5000 square-exact.json
target/release/ehz-capacity --max-steps 10000 --certify-interval square-exact.json
```

Both limits include domain preprocessing. Zero is a valid limit. A stopped
domain check returns `incomplete`, not `invalid-input`; a stopped search preserves
any feasible incumbent and never labels it exact-complete. Limits are checked
between exact operations, so a single large rational elimination can exceed a
requested elapsed-time limit. Use these cooperative limits to obtain a partial
certificate; forcibly terminating the process can prevent JSON output entirely.

## Independent verification

`scripts/verify_exact.py` uses Python's `fractions.Fraction` and the independent
paper oracle's free-coordinate stationary equations. It checks input identity,
row scales, domain and redundancy witnesses, nonnegative weights, closure,
normalization, the ordered objective, and the reciprocal. For a certified interval
it recomputes the objective envelope and both endpoints. For an exact-complete
claim it independently enumerates supports and full permutations, checking the
global optimum instead of trusting the optimizer's search counters.

```sh
python scripts/verify_exact.py square-certificate.json --input square-exact.json
python scripts/verify_exact.py square-certificate.json --witness-only
```

The second command checks the witness without certifying global optimality,
even if the optimizer claimed completion. Verification is exponential and may
take longer than the Rust solver. An exact feasible witness alone establishes
an upper bound on capacity; search logs or counts alone never certify coverage.

The library API is `exact::solve(&ExactInput, &SearchOptions)`; call
`exact::certify_interval(&mut result)` to add the optional partial-result bound.
See the
[implementation guide](RATIONAL_SOLVER.md) for the mathematical contract and
acceptance gates. Ordering DP, symmetric-product dispatch, numerical filters,
fraction-free elimination, modular reconstruction, and general exact V-to-H
conversion are subsequent work.
