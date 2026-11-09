# EHZ port: research and implementation recommendation

Research date: 2026-11-09. Upstream inspected:
`pazithaimkislev/EHZ-capacity`, commit
`22dd473cc8b9b2397b44cd13acbceb763dd978e9` (2024-05-25).

This is a design memo and three small validation probes. It is not the production
port, a MATLAB compatibility validation, or a performance benchmark. Both project
checkouts were left unchanged. No environment configuration was saved or published.

## Recommendation

Use Rust for a command-line tool and reusable numerical/combinatorial library.
Use Python for independent small reference calculations, fixture generation, and
research notebooks. Implement a transparent exhaustive reference algorithm first;
add individually validated exact reductions and specialized algorithms afterward.

Do not select C merely to accelerate the existing loops: the dominant problem is
the number of orders, not MATLAB scalar arithmetic. C is entirely viable, especially
for a project primarily integrating Qhull, LAPACK, and GMP. Rust provides convenient
dependency/toolchain pinning, safe combinatorial state management, and parallelism
without a relevant inherent disadvantage for this workload. Its computational
geometry ecosystem is a tradeoff: an optional C geometry adapter may be sensible.

Tentative components: nalgebra for small dense floating-point linear algebra,
num-bigint/num-rational for a small exact reference path, serde for versioned data,
and Rayon for independent jobs. Verify versions and backend behavior during
implementation; none of these packages was installed or benchmarked for this memo.
Choose an LP backend with accessible solutions and dual certificates; avoid making
a proprietary optimizer a requirement. A C ABI or Python binding can follow a
working CLI and stable input format.

## What the existing code actually computes

The input is a vertex matrix in `(q_1,...,q_n,p_1,...,p_n)` coordinates. It computes
a volume centroid using a Delaunay triangulation, obtains hull simplices with
Qhull's `QJ` option, recovers normals, approximately merges matching normals, and
builds their symplectic pairing matrix. `EHZ_perms.m` then enumerates every support
of facets and every cyclic order of that support. It removes dependent equality
constraints and solves a KKT linear system for each order, or directly scores the
order when the equality constraints already determine the weights.

With `m` facets, the number of possible cyclic support orders before filtering is

`sum_{s=2}^m binomial(m,s) (s-1)!`.

This is 16,064 at m=8; 1,112,073 at m=10; 119,481,284 at m=12; and
349,096,664,728,623,316 at m=20. The original already fixes cyclic rotation by
anchoring each support's first facet. That reduction should not be counted as a
new improvement. It still materializes combination/permutation arrays and arrays
of scores, which creates avoidable memory costs.

Using an H-representation `K = {x: a_i.x <= h_i}`, let
`Omega_ij = omega(a_i,a_j)` and

`M = { beta >= 0 : sum_i beta_i a_i = 0, sum_i beta_i h_i = 1 }`.

The theorem gives `c_EHZ(K) = 1/(2 max Q_sigma(beta))`, where
`Q_sigma(beta) = sum_{j<i} beta_sigma(i) beta_sigma(j) Omega_sigma(i),sigma(j)`.
The code minimizes this signed quadratic expression and returns `-1/(2 min Q)`.
Reversing the full ordering changes its sign, explaining the equivalence. Flow
orientation, coordinate ordering, and factors of two need explicit regression tests.

Rows need not be unit normals: positive rescaling of an inequality can be absorbed
in beta. The arbitrary-row version is also stated in [2]. Thus rational facet
inequalities allow exact arithmetic without square roots from normalizing vectors.

## Improvements justified by existing mathematics

### 1. Transition graph and cycle enumeration (general polytopes)

[1, Remark 3.11] is more precise than a generic permutation heuristic. It defines
a directed graph on facets: an edge i->j exists when a positive segment in the
characteristic direction of facet i can start on facet i and end on facet j.
The search may be restricted to cycles of this graph. Together with [1, Theorem
1.5], at least one minimizer uses each active facet direction at most once.

For a chosen direction `v_i`, pairwise reachability can be checked with a linear
program in x and t:

* `A x <= h`, `A(x+t v_i) <= h`;
* `a_i.x = h_i`, `a_j.(x+t v_i) = h_j`;
* maximize `t >= 0`.

The direction and whether the formula uses max Q or min Q must be consistent.
In the probes `v_i = (a_i,p, -a_i,q)`, the MATLAB pairing-matrix convention.
Each facet direction is only specified up to positive scaling for reachability.
Compactness bounds t, and a positive optimum establishes a possible transition.
When numerical classification is uncertain, retain the edge or certify the LP.

Stream simple cycles and minimum scores/witnesses, use canonical anchors, eliminate
graph regions that cannot lie on a cycle, and test support feasibility before
expensive work. Pairwise reachability is necessary; it does not establish a single
geometrically consistent orbit through an entire cycle. Do not replace this graph
with ordinary ridge adjacency or discard zero symplectic pairings without a proof:
segments can lie in lower-dimensional intersections, including Lagrangian faces.
Worst-case cycle counts remain large.

### 2. Central symmetry (already in the original paper)

[1, Corollary 1.2] gives a separate formula for centrally symmetric polytopes.
It uses one representative per opposite facet pair, signed coefficients, and a
weighted absolute-value constraint instead of the original closure equations.
This potentially reduces a problem on m facet directions to m/2 paired directions,
with the sign choices accounted for. It is more than just dividing a count by two.

The MATLAB implementation does not use this corollary. It is especially relevant
to cubes, cross-polytopes, and the README example. Verify symmetry about a common
center before dispatching to it; approximate symmetry is not an exact reduction.
Other automorphisms must preserve the relevant heights and symplectic pairings.

### 3. Weighted ordering and subset DP (unique weights, particularly simplices)

[2] proves NP-hardness of general EHZ computation, even for simplices, via feedback
arc sets in bipartite tournaments. Its rational decision problem is NP-complete.
Dimension varies in this hardness result; it does not settle complexity in fixed
dimension four. It rules out expecting a uniform polynomial exact solver unless
P=NP, not substantial practical improvements.

For a full-dimensional simplex the closure and normalization equations determine
beta uniquely. For any support with uniquely determined feasible beta, define
`w_ji = beta_j beta_i omega(a_j,a_i)`. The remaining optimization is the weighted
linear-ordering problem. A standard last-element subset recurrence is

`D(S) = max_{j in S} [ D(S\{j}) + sum_{i in S\{j}} w_ji ]`, with `D(empty)=0`.

The straightforward implementation takes O(s^2 2^s) arithmetic operations and
O(2^s) stored scores, instead of enumerating (s-1)! cyclic orders. Precomputing
all subset sums can trade more memory for fewer operations. Closure implies each
row of w sums to zero; anchoring one facet then permits a DP on the remaining
s-1 facets. This is a design application of a standard ordering recurrence, not
a claimed new theorem or an algorithm asserted to appear in [2].

When beta is not unique, one scalar score per subset is insufficient: the best
partial order depends on beta. Applying the same DP to arbitrary polytopes would
be unjustified. A weighted-ordering MILP with transitivity constraints is another
candidate for this special case, to benchmark against DP rather than require.

### 4. Lagrangian products: a different exact formulation

**Strongest immediate special case:** if both K and T are centrally symmetric
about the origin, [9, Remark 4.2 / equation (4.2)] gives

`c_EHZ(K x T) = 4 inrad_(T^polar)(K)`.

For symmetric H-representations `|u_i.q| <= h_i` and `|v_j.p| <= g_j`, with
positive heights, this becomes

`c_EHZ(K x T) = 4 / max_(i,j) [ |u_i.v_j| / (h_i g_j) ]`.

Indeed `h_(T^polar)(u_i)` is the gauge of T at u_i, equal to
`max_j |v_j.u_i|/g_j`, and containment of `r T^polar` in K gives the formula.
After the H-representations are available, this costs O(m_K m_T n) arithmetic
operations. Rational inputs give an exact rational answer. This bypasses all
permutation search for this class, in any dimension, and is the first specialized
fast path to implement. Both factors must have verified central symmetry; a
general centrally symmetric body in R^(2n) need not be a Lagrangian product.
Translations of each factor can first be removed.

For nonsymmetric products, the following broader billiard route remains useful.

[3, Theorem 1.3] proves the Minkowski billiard characterization for arbitrary
convex bodies K,T, including nonsmooth polytopes. The capacity of `K x T` can be
realized by a shortest **strong** Minkowski billiard with at most n+1 bouncing
points in K. In R^4, K and T are planar, so only two- and three-bounce cases are
needed. This bound is about projected billiard bounces, not the number of facets
in a general 2n-dimensional polytope itinerary.

[4, Section 7.3] develops a planar polytopal algorithm based on face selections
and LPs, with existing Python code [10]. It quotes O(v_K^2 v_T^2) iterations for
two-bounce cases and O(v_K^3) outer iterations for three-bounce cases, each
requiring further LP/geometry work. Its treatment of degenerate normal cones
needs scrutiny before claiming a fully certified implementation; a discussed
normal-direction sampling alternative is not automatically complete.

For fixed face selections, normal-cone membership gives linear relations between
q/p positions and cone coefficients. [4, Proposition 7.1] explains why the length
is constant among trajectories with the same relative-interior face data. This is
a strong basis for a dedicated product backend. Retain full normal cones and
lower-dimensional faces when completeness requires them.

Weak billiards are not an acceptable substitute: [3] explicitly gives situations
where their shortest length is smaller than that of strong billiards. The README
example is `Q x Q^polar`, for centrally symmetric Q, whose capacity is exactly 4
by [9]. Its special structure is worth accepting explicitly in the input rather
than losing it in a Cartesian-product vertex cloud.

### 5. Four-dimensional Reeb dynamics

[5] supplies a combinatorial Reeb-flow framework with action/rotation bounds and
an EHZ characterization (Corollary 1.15). Its symplectic-polytope hypothesis
excludes Lagrangian 2-faces; Type 2 orbits also matter for completeness. These
conditions prevent treating it as a universal replacement for the general solver.
It is a possible later backend for suitable four-dimensional inputs. Symmetry
and product polytopes can be precisely the degenerate examples requiring care.

## Additional formulation to validate: an LP per physical itinerary

This follows from piecewise-facet characteristics and [1, Theorem 1.5]; it is a
design derivation, not a claim of novel prior art. It differs from optimizing the
quadratic formula over beta for an arbitrary formal ordering.

For a fixed itinerary `(i_0,...,i_{s-1})`, use starting point z and segment
coefficients `lambda_r >= 0`, with

`x_r = z + sum_{t<r} lambda_t J a_i(t)`.

Impose `A x_r <= h`, `a_i(r).x_r = h_i(r)`, and
`sum_r lambda_r J a_i(r) = 0`. Segment endpoints lie in the polytope; convexity
places the whole segment there. Since `a_i.J a_i = 0`, the segment lies in its
named facet. Its velocity is in that facet's characteristic direction. With
consistent positive-action orientation the action is linear:

`action = (1/2) sum_r lambda_r h_i(r)`.

Minimize this action with an LP. Selected normals must admit normalized
nonnegative closure (the support-feasibility test above), excluding constant-loop
solutions: otherwise zero durations and a shared point of all named facets would
give a spurious zero answer. To see why the test suffices, any common facet point
would imply `sum beta_i h_i = (sum beta_i a_i).z = 0`, contradicting normalization
to 1. Zero-duration individual segments may be removed.

Every feasible nonconstant solution is a closed generalized characteristic, hence
its action is at least the capacity. Conversely, Theorem 1.5 guarantees a minimum
represented by an itinerary with no repeated active direction. Enumerating the
eligible itineraries therefore recovers the capacity. A production implementation
must preserve those hypotheses, validate degeneracies, and certify decisions.

This formulation may avoid repeated indefinite quadratic optimization and can
test geometric consistency that the pair graph cannot see. It has more inequalities
than a small KKT system, so it is **not yet established to be faster**. Compare
it with the faithful reference before choosing a default general backend. LP
solutions and infeasibility decisions would need rational primal/dual certificates
for certified search/pruning.

## What newer literature changes, and what it does not

* [2] provides the complexity limit and a route to mature ordering techniques.
* [6] (first preprint 2024, revised 2025) disproves Viterbo's general
  volume-capacity conjecture. A product of a regular unit-circumradius pentagon
  and its quarter-turn has capacity `2 cos(pi/10) (1+cos(pi/5))` and systolic
  ratio greater than 1. Include this regression case. Do not use the conjectured
  volume bound, or equality with Gromov width, as an unconditional pruning rule.
* [7] (2026, v2) proves new Viterbo inequalities for products with quadrilaterals
  and affinely regular hexagons and studies equality cases. This can help with
  specialized checks/bounds, but is not a general exact capacity algorithm.
* [8] (September 2026 preprint) reports exact formulas for certain rotated
  Lp-ellipsoids (with specified p and block-linear structure), and a Fourier-loop
  numerical optimizer. The optimizer is explicitly local/multistart; small
  residuals do not certify the globally minimal action. Potential special-family
  dispatch and candidate generation deserve testing, not substitution for exact
  polytope search. The p->infinity parallelepiped-product specialization is already
  covered by the older symmetric-product inradius identity above.
* [11] supplies constructed symplectically self-polar polytopes with exact
  capacities `2+1/n`, useful as additional structured validation families.
* [12] (2025 preprint) studies nonsmooth action minimizers and Zoll/cut-additivity
  properties under explicit hypotheses. It does not justify adding capacities
  of arbitrary simplex pieces.

I did not find a general shortcut that makes simplicial polytopes easy. Simplices
are already hard when dimension varies. Triangulating a polytope and summing
capacities is invalid in general. Likewise, there is no general 2n+1 bound on
active facet directions: a square in R^2 requires four edge directions to trace
its capacity-realizing boundary, while 2n+1=3. The billiard bounce bound has a
different domain of applicability.

## Numerical changes required for rigor and reproducibility

The centroid's exact value is unnecessary: translation invariance means any
certified interior point works. For a full-dimensional finite vertex cloud the
arithmetic mean is strictly interior to its convex hull. This can remove the
entire volume-triangulation calculation. For H-input, find/certify an interior
point directly. Detect and recenter actual central symmetry separately.

Offer rational H-input before coupling correctness to a hull implementation.
Validate even dimension, boundedness, full dimension, redundant/scaled facet
descriptions, and strict interior placement as appropriate. A V-input adapter
must reconstruct true facets and incidence information. Qhull `QJ` deliberately
joggles input; a joggled triangulation and tolerance-merging step cannot silently
serve as an exact hull certificate. For floating inputs, explicitly report what
geometry has been validated and retain uncertain cases.

Replace the absolute determinant cutoff `abs(det(KKT)) < 1e-10` with suitable
factorizations, rank/residual checks and exact fallback. Determinants are badly
scaled rank tests. In the exact exhaustive reference, dropping a singular KKT
system needs the argument that an optimal point can move in a flat feasible
direction to a smaller support with the same objective; simply calling a matrix
numerically singular is not that proof. Considering feasible saddle stationary
points is not inherently a mathematical error in the existing exhaustive search:
they cannot improve on the true global extremum over the feasible set.

Distinguish an estimated number, a certified interval, and an exactly verified
capacity. A feasible beta and positive Q gives an upper bound on capacity, not a
proof of optimality. If the global maximum Q is bounded by `0<L<=Qmax<=U`, then
`1/(2U) <= capacity <= 1/(2L)`. A witness alone is not a certificate of exhaustive
coverage. Certified mode needs justified pruning and global coverage or a bound
gap, not just an accurate final linear solve.

Pin toolchains, dependency versions and source commits; store input hashes,
fixture provenance, conventions, solver settings, numerical status, certificates,
and search counts. Keep deterministic task ordering/tie-breaking and explicitly
state the level of reproducibility; pinned Rust alone does not guarantee bitwise
identical floating-point behavior across CPUs and solver builds.

## Validation and implementation order

1. Freeze the upstream commit and develop the small exact formula oracle.
   Obtain MATLAB R2021b comparison outputs when available and preserve the input,
   facet data, and options that generated them. No MATLAB/Octave execution or
   output comparison was performed in this research turn.
2. Implement the general reference Rust CLI/library, rational H-input, witnesses,
   and a documented V-input geometry adapter. Stream enumeration. Verify the
   planar area identity, boxes, simplices, and the README example.
3. Add translation, quadratic scaling, facet-row rescaling, and rational
   symplectic-shear invariance tests, plus redundant vertices, nonsimplicial facets,
   near-degeneracies, and multiple optimal orbits. Cross-check exact and numerical
   paths on small random instances.
4. Validate transition-graph filtering and each other pruning rule against the
   unpruned reference. Compare the itinerary LP backend independently. Add the
   unique-weight ordering DP and verified central-symmetry formula.
5. Implement the symmetric-product inradius/pairing formula first, and the
   specialized planar nonsymmetric Lagrangian-product backend if these inputs
   are important, cross-checking with the general solver and published examples.
6. Add parallelism and bounds after correctness. Benchmark by dimension, true
   facet count, symmetry, graph density, numerical conditioning, and backend.
   Report candidates, rejected candidates, KKT/LP counts, peak memory, wall time,
   solver versions, seeds, and exactness/gap status. Make no speedup claim until
   measured against the same workload and accuracy requirements.

The original repository and [10] expose no license in the inspected repository
metadata/root listings. Confirm distribution terms for any literal source
translation or reuse; the mathematical specification and attribution should be
explicit in an independent implementation. This does not block the design work.

## Validation probes completed

`check_ordering_dp.py` uses only standard-library exact fractions. It compared
the subset DP with full permutations on 18 seeded rational skew-weight instances,
and with all anchored cyclic orders for four rational simplices in dimensions
2,4,6,8 (up to 40,320 cyclic orders). All 22 comparisons passed. This verifies
these instances and recurrence implementation, not a general solver.

`check_itinerary_lp.py` uses the environment's existing NumPy/SciPy/HiGHS. It
enumerated all anchored orders for three rational simplices in dimensions 2,4,6
and compared their minimum LP actions with exact ordering-DP capacities. It also
recovered area 4 for `[-1,1]^2` and an action-4 orbit on four selected facets of
the README product. The opposite README itinerary was rejected for the chosen
flow orientation. These LP checks use floating-point tolerances, not formal
certificates. The product check is not an exhaustive eight-facet search.

`check_symmetric_product.py` checks the direct symmetric-product formula and
constructs an exact four-segment orbit attaining it. It checks the README product,
diamond times diamond, interval times interval, and 20 seeded rational products
with factor dimensions 1 through 4. All containment, closure and action checks
passed with exact fractions. The lower-bound justification is the published
inradius theorem, not the finite sample checks.

Commands, from any directory with the noted dependencies:

```sh
python /workspace/scratch/ehz-research/check_ordering_dp.py
python /workspace/scratch/ehz-research/check_itinerary_lp.py
python /workspace/scratch/ehz-research/check_symmetric_product.py
```

Logs: `ordering-results.txt`, `itinerary-results.txt`,
`symmetric-product-results.txt`. Source PDF hashes and
environment versions: `sources.json`, `probe-environment.json`.

## Primary references

1. Pazit Haim-Kislev, *On the Symplectic Size of Convex Polytopes*, GAFA 29
   (2019), 440–463. Theorems 1.1,1.5; Corollary 1.2; Remarks 3.8–3.11.
   <https://arxiv.org/abs/1712.03494v3>;
   <https://doi.org/10.1007/s00039-019-00486-4>.
2. Karla Leipold and Frank Vallentin, *Computing the EHZ capacity is NP-hard*,
   Proc. Amer. Math. Soc. Ser. B 11 (2024), 603–611.
   <https://arxiv.org/abs/2402.09914v3>;
   <https://doi.org/10.1090/bproc/250>.
3. Daniel Rudolf, *The Minkowski Billiard Characterization of the EHZ-capacity
   of Convex Lagrangian Products*, v2 (2022).
   <https://arxiv.org/abs/2203.01718v2>.
4. Daniel Rudolf and Stefan Krupp, *Shortest Minkowski billiard trajectories on
   convex bodies*, v1 (2022), especially Theorem 1.3 and Section 7.3.
   <https://arxiv.org/abs/2203.01802v1>.
5. Julian Chaidez and Michael Hutchings, *Computing Reeb dynamics on 4d convex
   polytopes*, v2 (2021). Definition 1.4 and Corollary 1.15.
   <https://arxiv.org/abs/2008.10111v2>.
6. Pazit Haim-Kislev and Yaron Ostrover, *A Counterexample to Viterbo's
   Conjecture*, first posted 2024; inspected v3 (2025), Proposition 1.4.
   <https://arxiv.org/abs/2405.16513v3>.
7. Alexey Balitskiy, Ivan Mitrofanov and Alexander Polyanskii, *Triangle covering
   problems and the Viterbo inequality in the plane*, v2 (2026).
   <https://arxiv.org/abs/2603.12495v2>.
8. Vardan Oganesyan, *The Ekeland–Hofer–Zehnder Capacity of rotated Lp-Ellipsoids*,
   v1 (2026). Theorem 2 and numerical Appendix 6.
   <https://arxiv.org/abs/2609.37794v1>.
9. Shiri Artstein-Avidan, Roman Karasev and Yaron Ostrover, *From Symplectic
   Measurements to the Mahler Conjecture*, Duke Math. J. 163 (2014), 2003–2022.
   Theorem 1.7 and Remark 4.2 (general symmetric products).
   <https://arxiv.org/abs/1303.4197v2>.
10. Stefan Krupp, code accompanying *Calculating the EHZ Capacity of Polytopes*.
    Inspected repository metadata/root listing; code is a ZIP archive and was not
    executed or exhaustively audited in this turn.
    <https://github.com/S-Krupp/EHZ-Capacity-of-polytopes>.
11. Mark Berezovik, *Symplectically self-polar polytopes of minimal capacity*,
    inspected v2 (2025), Theorem 1.2.
    <https://arxiv.org/abs/2310.14998v2>.
12. Pazit Haim-Kislev, *Dynamical extensions of Zoll to nonsmooth convex bodies*,
    v1 (2025). Inspected abstract and introduction/results in arXiv HTML.
    <https://arxiv.org/abs/2511.16644v1>.

The broad arXiv searches and inspected source downloads are retained alongside
this memo. Search coverage is not a proof that no other useful result exists.
No claimed running-time improvement in this memo depends on an unverified
volume-capacity conjecture or on a withdrawn preprint.
