# Optimization implementation guide

Reviewed against main at `f21d409d53fe742c41e82abf11948f2dc524df8e`.
This is the implementation sequence and acceptance contract for subsequent
optimization work. It updates the priorities in the archived
[research memo](../research/RESEARCH.md), without treating its probes as shipped
backends. It assumes familiarity with the mathematics and numerical optimization;
it specifies algorithms, invariants, and evidence rather than source templates.

**Sequencing update:** [RATIONAL_SOLVER.md](RATIONAL_SOLVER.md) specifies the
next milestone and takes precedence on arithmetic, verification, and delivery
order. This document retains the broader optimization roadmap.

## 1. Baseline and priorities

The current production path in `src/lib.rs` enumerates every support and anchored
cyclic order, with unique equality solutions or a stationary KKT solve. It has
dimensionless weights, relative rank tests, reorthogonalized equality reduction,
full-constraint residual checks, and a numerical positive-spanning input check.
These are already implemented; do not schedule them again as new optimizations.
Neither an exact production backend, ordering DP, transition graph, physical
itinerary LP, nor symmetric-product dispatch is implemented. H-input currently
means binary64 JSON numbers, not rational input with exact provenance.

Preserve this path as an explicitly selectable reference backend. Change one
reduction at a time and compare both paths on the same normalized input. The
next milestone is the complete rational H-input reference solver and independent
verification specified in [RATIONAL_SOLVER.md](RATIONAL_SOLVER.md), before
verified symmetric-product dispatch and unique-weight DP. General graph pruning
and itinerary LPs follow only after their completeness and numerical contracts
are implemented. Parallelism comes after deterministic serial validation.

The bottleneck is combinatorial coverage, not just the speed of a dense solve.
For m facets, the unfiltered cyclic-support count is

\[
N(m)=\sum_{s=2}^{m}\binom ms(s-1)!.
\]

Production's `orders` counts only orders entered after equality preprocessing;
it is not N(m). `supports` counts all size-at-least-two subsets. Do not report a
reduction against N(m) as measured speedup against the existing implementation.
Boundedness validation itself enumerates small supports and deserves a separate
timer: on larger problems it can dominate before the capacity search starts.

## 2. Mathematical and numerical result contracts

Fix coordinates `(q_1,...,q_n,p_1,...,p_n)` and
\(\omega(a,b)=a_q\cdot b_p-a_p\cdot b_q\). Let

\[
M=\{\beta\ge0:A^T\beta=0,\ h^T\beta=1\},\qquad
Q_\sigma=\sum_{j<i}\beta_{\sigma(i)}\beta_{\sigma(j)}
\omega(a_{\sigma(i)},a_{\sigma(j)}).
\]

The target is \(Q_* = \max_{\sigma,\beta\in M}Q_\sigma\), with
\(c=1/(2Q_*)\). The production minimum convention is equivalent under complete
order reversal. Keep weights aligned with the returned order under reversal and
relabeling. A formal optimizing weight/order witness is not yet a reconstructed
closed characteristic on the boundary.

In the current dimensionless coordinates, \(t_i=\beta_i h_i\) and
\(P_i=\ell a_i/h_i\), with \(\ell=\min_i h_i/\|a_i\|\), so
\(\sum t_i=1\), \(\sum t_iP_i=0\), and \(Q=Q_P/\ell^2\).
Every backend must document whether objectives, bounds, residuals,
and tolerances refer to original or normalized coordinates. Transport both ends
of a bound interval through the same scaling. Common scale and anisotropic
conditioning are different test axes; normalization addresses only part of the
latter.

Introduce distinct outcomes before interruptible or heuristic search:

| Outcome | Required evidence |
| --- | --- |
| Numerical estimate | Feasible to declared tolerance; coverage and numerical rejection counts disclosed |
| Feasible candidate | Independently checked weight/order or physical orbit; no global claim |
| Certified capacity interval | Valid global objective bounds and feasible incumbent, including certified pruning decisions |
| Exact capacity | Exact arithmetic or exact reconstruction plus exhaustive justified coverage |
| Incomplete / numerically unresolved | Budget exhaustion or unresolved decisions; retain incumbent and any valid bound |
| Invalid input | Verified domain failure, distinguished from inability to establish boundedness at tolerance |

A feasible positive objective L is a lower bound on Q*, hence an **upper** bound
on capacity. If \(0<L\le Q_*\le U\), report
\([1/(2U),1/(2L)]\). If L=0 there is no finite incumbent capacity bound. Without
U there is no global gap. Rounded weights that merely pass a residual tolerance
are not a rigorous feasible witness. For normalized objective error,
\((c_{candidate}-c_*)/c_*=(Q_*-Q_{candidate})/Q_{candidate}\); denominator inversion
can magnify optimization error. Track the objective as well as capacity.

## 3. Input, geometry, and exactness boundary

Keep H-input as the mathematical validation interface. A future rational schema
must preserve integer numerators and denominators without a detour through
binary64. Reject malformed dimension, nonfinite entries, zero rows, and invalid
interior normalization before dispatch. General H-input not centered at an
interior origin requires a separate interior-point stage; the present positive
height interface must remain documented until that stage exists.

Replace the expensive boundedness routine only with an equivalent verified
test. With full row span, maximizing delta under \(P^Tt=0\),
\(\sum t_i=1\), \(t_i\ge\delta\) is one suitable LP formulation. A positive
verified delta proves positive spanning. A floating solver's tiny or negative
delta does not prove unboundedness; retain an unresolved category and obtain a
recession witness or exact fallback where needed. Redundant inequalities can
force tiny positive coefficients, so absolute coefficient thresholds are not
geometric certificates.

For V-input, separate hull time from optimization time. Replacing the volume
centroid by the vertex mean is a legitimate optional simplification for a
validated full-dimensional cloud, but changes the MATLAB-compatible geometry
path and must be independently tested. QJ joggling, tolerance merging, and
containment checks alone do not establish an exact original hull: certification
also needs completeness of its facet description. Deduplicate exact positive
row multiples; remove a generally redundant inequality only with a redundancy
certificate. Preserve maps to original facet indices and beta units.

Do not accept approximate symmetry, an approximately symplectic transformation,
or almost-product sparsity as an exact dispatch condition. They may suggest a
candidate or a nearby comparison body. Any certified capacity conclusion for
the original body must use containment bounds or exact predicates.

## 4. Improve the exhaustive kernel without changing its search space

Cache support-level equality reductions and feasible affine parameterizations.
Keep matrices in fixed support coordinates and change the ordering-dependent
quadratic coefficients; this avoids repeatedly permuting equalities. Reuse
buffers and deterministic enumeration state. Precompute pairings once. Profile
allocations, equality reduction, and KKT factorizations separately before choosing
a faster factorization over the current SVD.

For \(t=t_0+Zy\), stationary equations are
\((Z^TH_\sigma Z)y=-Z^TH_\sigma t_0\). Use a rank-revealing equality factorization
and an appropriate symmetric-indefinite or rank-revealing reduced solve. This can
reduce matrix size substantially, but benchmark against the KKT formulation:
the quality of Z, residual propagation, and ill-conditioning may outweigh the
dimension saving. Re-evaluate candidates against the original full equalities,
nonnegativity, and the directly summed ordered objective after any clipping.

Separate inconsistent equalities, infeasible weights, zero objective, exact
singularity, numerical rank ambiguity, failed residual checks, and nonfinite
arithmetic in diagnostics. The existing `singular_kkt` combines numerical solve
rejections and cannot serve as a proof of safe omission. Equality/unique-solve
rejections are not currently covered by that counter.

An exactly singular stationary block can be omitted because its flat feasible
direction reaches a smaller support with unchanged objective; see the proof in
[PAPER_VALIDATION.md](../validation/PAPER_VALIDATION.md). The same inference is
invalid for a merely small singular value. Exact rational fallback on small
ambiguous supports can reuse the rational reference backend described in the
new guide; that complete backend is now the first milestone. Higher precision
with stability across precisions is useful evidence, but is not automatically a
certificate. Retain unresolved supports in any global upper bound.

Do not require a positive-definite Hessian: the fixed-order quadratic need not
be convex. A feasible saddle can be scored harmlessly in exhaustive search, but
stationarity and small residuals alone never establish optimality. Do not impose
a d+1 support cutoff: the planar square requires four directions. The d+1 bound
for vertices of the closure polytope does not bound support of a quadratic
extremum.

## 5. First specialized backends

### Verified symmetric Lagrangian products

Accept explicit factors, or verify that H-rows split into the designated q and p
blocks and that each bounded factor is centrally symmetric about a known center.
For \(K=\{|u_iq|\le h_i\}\), \(T=\{|v_jp|\le g_j\}\), use

\[
c(K\times T)=\frac4{\max_{i,j}|u_i\cdot v_j|/(h_i g_j)}.
\]

Return the maximizing facet pair and the four-segment characteristic described
in the archived symmetric-product probe. Verify containment, facet contact,
closure, characteristic orientation, and action independently. The inradius
theorem supplies the lower bound; the orbit proves attainment. Rational data
permit exact comparisons and an exact answer without normalizing normals.

Require both factors to be symmetric; general central symmetry in the ambient
space is insufficient. A symplectically transformed product may be accepted only
with a verified transformation and inverse. Benchmark detection and witness
verification too. Test boxes with tied minimizing planes, diamond products,
skew factor inequalities, redundant rows, and a nonsymmetric near-miss that must
fall back. General central symmetry via Corollary 1.2 is a separate later backend,
not an extension of this constant-pair formula.

### Unique-weight ordering DP

Dispatch only when equality constraints uniquely fix feasible weights on the
support under consideration. Compute
\(w_{ji}=\beta_j\beta_i\omega(a_j,a_i)\) and

\[
D(S)=\max_{j\in S}\left[D(S\setminus\{j\})+
\sum_{i\in S\setminus\{j\}}w_{ji}\right],\qquad D(\varnothing)=0.
\]

Store predecessor choices, reconstruct an order, and rescore it directly.
Tie-breaking must be deterministic but tests must allow different optimal orders.
The full recurrence works for arbitrary weights. Anchoring and deleting one
facet works here because closure makes all row sums of w zero; skew symmetry
alone does not justify it. An anchor contributes zero in total under closure.

The straightforward cost is O(s² 2^s) arithmetic and O(2^s) score/predecessor
storage. Precomputed subset sums can reduce arithmetic at O(s 2^s) storage.
Check bit-width, allocation bounds, and rational integer growth before allocating.
An explicit memory budget must select a fallback or return incomplete status,
never truncate the DP silently. Exact integer-scaled rational weights are a
natural initial backend; floating tie decisions need intervals or fallback for
certification.

For a whole simplex this replaces ordering enumeration. For a general polytope
it accelerates only eligible supports. Do not freeze one arbitrary feasible beta
and call the resulting DP an exact capacity solver. The test suite includes a
square with a fully solved ordering problem at suboptimal fixed weights.

## 6. General pruning and the physical-itinerary alternative

Start with feasibility pruning of a support's closure/normalization equations.
Reject only a certified infeasible support in certified mode. LP status, solver
timeouts, and small negative slacks must not all collapse into infeasibility.

For graph reductions use the positive characteristic transition graph from
Remark 3.11, not ridge adjacency or signs of individual pairings. Certify absent
edges by the reachability LP; uncertain edges stay. Remove vertices outside
cyclic strongly connected components, enumerate simple directed cycles with
canonical anchors, and retain lower-dimensional contact possibilities.

This graph concerns **physical characteristics**. Its use to restrict the formal
Q search needs the existence theorem that at least one global minimizer has an
eligible simple itinerary and the matching orientation convention. It does not
say every formal feasible order has a physical realization. Write down that
completeness argument in the backend documentation, including removal of zero
segments; do not infer correctness solely from agreeing on random polytopes.

For a physical itinerary, the archived LP uses positions
\(x_r=z+\sum_{t<r}\lambda_tJa_{i_t}\), \(\lambda\ge0\), boundary contact,
containment, closure, and linear action \(\tfrac12\sum\lambda_th_{i_t}\).
Its J convention must be matched to the Q convention using an oriented square.
The present probe is a low-level routine with an **unchecked precondition**:
the selected normals admit nonnegative normalized closure. Without that test,
adjacent facets sharing a point admit a spurious zero-action constant loop.
Enforce the precondition in the production wrapper; adding an arbitrary epsilon
lower bound on action or each duration is not a valid replacement.

Verify returned positions and segments against every original inequality,
including the final endpoint, independently of the LP's success flag. Classify
infeasible, unbounded, unresolved, and nonconstant feasible outcomes separately.
Use primal/dual rational certificates before making rigorous completeness or
pruning claims. Compare minimum action with the exact formula on full small
searches. This LP is an alternative backend, not yet a demonstrated speedup.
Do not compare its objective for each formal order with that order's stationary
Q value: most formal orders are not geometrically realizable.

## 7. Bounds, anytime behavior, and parallel search

Candidate generation and proof search should share independently verified
incumbents. In maximization coordinates an incumbent improves L monotonically;
only a valid relaxation reduces U. For normalized t on a support of size s,
one elementary but weak upper bound is
\(Q_P\le\tfrac12 C(1-1/s)\), where
\(C=\max_{i,j}|\omega(P_i,P_j)|\). It follows from the simplex mass constraint
and bounds every ordering; ignoring closure makes it conservative. Restore the
factor \(\ell^{-2}\) before comparing to original-coordinate incumbents. Implement this
bound first as a sign/scale/pruning test, not as a promise of useful pruning.

Tighter candidates include certified LP/McCormick relaxations of pair products,
with rigorously bounded t ranges, or other explicitly proved quadratic
relaxations. A local optimizer, stationary value, heuristic ordering, or sampled
orbit is not an upper bound on Q*. Partial-order bounds must cover every
completion and every feasible weight; a fixed-beta ordering bound does not.
Reject a search region only when its verified U cannot beat the incumbent.

On interruption, aggregate upper bounds for all unresolved regions, retain the
incumbent, and report the coverage ledger. Worker failures must leave their
regions unresolved. Deterministic task identifiers, canonical witness comparison,
and a deterministic reduction of equal scores are needed before claiming
reproducibility. Checkpoint input hashes, algorithm settings, completed regions,
and the incumbent/bound provenance. Never use the disproved general Viterbo
volume inequality as a bound.

## 8. Acceptance tests and quality measurements

The new executable contract is [scripts/test_optimization.py](../scripts/test_optimization.py).
Its mathematical mode checks exact ordering, singular/redundant examples,
counterexamples to unsafe reductions, independent planar areas, and physical LP
witnesses. Binary mode additionally checks production answers, both printed
witnesses, transformed free-weight cases, product identities, multiple optima,
conditioning, and invalid inputs. All checks remain active under Python -O.
The existing source-paper and legacy suites remain complementary checks.

| Layer | Required evidence before enabling an optimization |
| --- | --- |
| Formula and convention | Ordered score, reversal sign, mass/closure, factors of two, planar area |
| Weight optimization | Exact free-weight nonsymmetric fixtures; boundary optima; singular plateaus; demonstrably suboptimal feasible candidates rejected by quality checks |
| Ordering | Full permutations versus DP; predecessor reconstruction; zero/tied weights; closure-required anchoring counterexample |
| Dispatch | Verified family predicates, exact known values, near-misses falling back, same result with dispatch disabled |
| Pruning | Every retained optimum matches unpruned exact search; each removed region's certificate checks; uncertain decisions retained |
| Physical LP | Closure precondition, no constant loops, segment containment/contact, action, reversed orientation, independent capacity agreement |
| Numerics | Common scales, independent row scales, symplectic stretching/shears, near-rank cases, explicit unresolved outcomes |
| Search lifecycle | Budget interruption, resume equivalence, failed-worker handling, monotone incumbent and valid bounds |

Test currently unimplemented contracts as each backend lands; do not add skipped
placeholders and describe them as coverage. In particular, the present suite
does not verify production pruning, dispatch, DP reconstruction, interruption,
or certification because those facilities do not yet exist.

Build benchmark panels independently of oracle test data. Vary dimension, facet
count, closure-polytope dimension, symmetries, product structure, graph density,
rank conditioning, and optimum multiplicity. Include deterministic seed lists
and adversarial perturbation families. Keep small exact fixtures for correctness
CI, medium cases for scheduled regression, and hard instances for budgeted search.
The regular-pentagon product counterexample belongs in the medium panel with its
published value and clearly identified irrational approximation error; do not
force a factorial ten-facet exhaustive run into every pull request.

Record input/source hashes, backend, build profile, toolchain, numeric precision,
geometry/domain-validation/search/certificate times, peak memory, candidate
counts, rejection reasons, incumbent Q, independent reference Q or valid U,
capacity error, and witness residuals. Compare paired instances and fixed
budgets. Report failures, timeouts, unresolved cases, medians and tails; do not
discard difficult seeds or average only successful runs. Wall-time thresholds
belong on controlled benchmark hardware, not as brittle unit-test assertions.

For small exact cases, acceptance requires relative objective and capacity error
at most 2e-8 and independent witness checks at the documented tolerance; exact
backends require equality. This is a regression tolerance, not a rigorous bound
or a guarantee for arbitrary ill-conditioned inputs. Beyond the oracle regime,
report certified gaps where available, otherwise improvement and agreement as
empirical evidence. No speedup is accepted at the cost of silently weaker
coverage, relaxed accuracy, or omitted input validation.

## References and scope

The archived [research memo](../research/RESEARCH.md) identifies primary sources
and the scope of each theorem: Haim-Kislev Theorem 1.1/1.5 and Remark 3.11;
Artstein-Avidan–Karasev–Ostrover's symmetric-product inradius identity;
Leipold–Vallentin's ordering complexity result; Rudolf's strong-billiard
characterization. The DP and proposed bounds above are algorithmic derivations,
not claims of new capacity theorems. The exact-singularity proof and numerical
baseline are documented in [the source-paper review](../validation/PAPER_VALIDATION.md).
