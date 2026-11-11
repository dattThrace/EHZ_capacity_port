# Rational solver: implementation guide

**Reference baseline implemented:** a complete exact solver for small rational
H-polytopes, with independent verification. This guide takes precedence over
the sequencing in [OPTIMIZATION.md](OPTIMIZATION.md). The original contract was
reviewed against `88a5e15` on 2026-11-09. See [EXACT.md](EXACT.md) for the shipped
input/output contract and [the baseline report](../validation/RATIONAL_BENCHMARK.md)
for acceptance evidence. Section 5 and fraction-free/modular elimination remain
subsequent work.

## 1. Mathematical contract

For bounded, full-dimensional $K=\{x:Ax\le h\}$ in standard symplectic
coordinates, use rational, unnormalized facet rows and a rational interior
origin, so $h>0$. Haim-Kislev's formula, written for arbitrary facet-row scales
in Leipold–Vallentin, is[^hk][^lv]

$$
c_{\rm EHZ}(K)=\frac1{2Q_*},\qquad
Q_* = \max_{\sigma,\,\beta\ge0,\,A^T\beta=0,\,h^T\beta=1}
\sum_{j<i}\beta_{\sigma(i)}\beta_{\sigma(j)}
\omega(a_{\sigma(i)},a_{\sigma(j)}).
$$

**Derived implementation fact:** rational input admits rational optimizing
weights and rational capacity. For a fixed order, choose a maximizer of minimal
positive support. With independent equality rows $E$, its stationary system
has rational coefficients:

$$
\begin{pmatrix}H_\sigma&E^T\\E&0\end{pmatrix}
\binom{\beta}{\lambda}=\binom0b,
\qquad Q_\sigma=\tfrac12\beta^TH_\sigma\beta.
$$

If singular, a kernel vector has a nonzero weight component $u$, with $Eu=0$.
Stationarity gives $\beta^TH_\sigma u=0$ and the kernel equations give
$u^TH_\sigma u=0$. Thus the objective is constant along $\beta+tu$.
Compactness allows reaching a smaller support, a contradiction. The minimal
system is therefore nonsingular and its solution rational. This argument
justifies exact singular-system omission **only with exhaustive smaller-support
coverage**; see also [the existing derivation](../validation/PAPER_VALIDATION.md).

The motivation is exact rank, feasibility, and ordering decisions. Classical
Diophantine small-divisor behavior is not an established obstruction here.
Exact arithmetic does not remove combinatorial hardness: the rational problem
is NP-hard even on simplices when dimension varies.[^lv]

## 2. Input and arithmetic foundation

- Add a versioned exact H-input mode with integer strings or numerator/denominator
  pairs; define decimal strings as exact decimal rationals. Reject zero
  denominators. Never pass exact input through `f64`. Conversion of a machine
  float must explicitly mean its exact binary value, not a recovered intended
  fraction.
- Use arbitrary-precision integers and canonical rationals. Rust's
  `num_rational::BigRational` is an appropriate initial representation.[^rust]
  Pin direct dependencies compatible with the existing toolchain; avoid a
  general symbolic-expression engine at this stage.
- Preserve unnormalized rows. Clear denominators and remove common factors
  using positive row scales only, recording their witness mapping. Euclidean
  unit normals and the current square-root-based length normalization do not
  belong in this backend. Rational rescaling is optional.
- Initially require positive heights and nonzero normals. Verify boundedness
  exactly: full column rank of $A$, together with a strictly positive rational
  vector $w$ satisfying $A^Tw=0$, is sufficient. The current small-support
  positive-spanning enumeration can be ported without a new LP dependency.
  Domain-check budget exhaustion means incomplete, not invalid.
- Deduplicate positive multiples of the same inequality. Accept other redundant
  inequalities through certified removal, or a separately justified extension
  of the formula to that representation. Preserve original row identities.
  General recentering and exact V-to-H conversion are subsequent work; the
  current QJ adapter cannot establish an exact result for the original vertices.

## 3. Complete reference algorithm

Retain the existing numerical backend under an explicit mode. Implement a
separate exact search over all supports and anchored cyclic orders, using
closure to justify cyclic reduction. Preserve zero-weight boundary cases;
there is no general dimension-plus-one support cutoff. The square already
requires four facet directions in dimension two.

Use exact row reduction to distinguish inconsistent, uniquely determined, and
free-weight equality systems. Solve each unique or nonsingular stationary
system exactly, check all original equalities and weight signs, and rescore
the ordered objective independently. The quadratic is generally indefinite;
neither a convex-QP assumption nor an eigenvalue tolerance is appropriate.
Classify inconsistent and exactly singular stationary systems separately.

Start with transparent rational elimination. After correctness is established,
benchmark fraction-free elimination on denominator-cleared integer systems;
Bareiss supplies the classical construction.[^bareiss] Handle zero pivots,
rank deficiency, and variable permutations explicitly. Cache support-level
equality reductions. Track peak integer bit lengths and memory as well as time.
Defer modular reconstruction until every reconstructed solution is checked
against the original rational equations.

## 4. Verification and result semantics

Keep the verifier independent of the optimizer's elimination routine. Serialize
canonical fractions, input identity, row transformations, support/order, weights,
objective, and completion status. Recheck $\beta\ge0$, $A^T\beta=0$,
$h^T\beta=1$, the objective, and its reciprocal using original input coordinates.

An exact positive witness gives an **upper bound on capacity**, not global
optimality. With verified $0<L\le Q_*\le U$, the capacity interval is
$[1/(2U),1/(2L)]$. An exact answer requires complete justified enumeration or a
verified proof covering all excluded regions. For small instances, independent
re-enumeration is an acceptable first implementation. Search counts and logs
alone are not a coverage certificate.

Expose distinct exact-complete, exact-witness/incomplete, certified-interval,
numerical-estimate, and invalid-input outcomes. On timeout, preserve any verified
incumbent and valid bound; never promote a partial search to an exact answer.

## 5. Acceleration after the exact baseline

Implement rational unique-weight ordering DP and verified symmetric-product
dispatch next, following their hypotheses in [OPTIMIZATION.md](OPTIMIZATION.md).
Compare each with the exact exhaustive backend. Exactness of individual solves
does not justify graph pruning or freezing arbitrary feasible weights.

Introduce numerical acceleration through certified filters: accept a comparison
only when an error enclosure determines its sign; otherwise recompute exactly.
CGAL documents this separation between predicates, constructions, and filtered
evaluation.[^cgal] Approximate search may suggest orders or bases, but must not
discard an unresolved region.

For future feasibility, redundancy, transition, and itinerary LPs, evaluate an
established rational solver such as SoPlex before building a simplex engine.
Its exact mode combines numerical work with rational checks; exact parsing and
output verification must be enabled explicitly.[^soplex] Check rational primal/dual
solutions or infeasibility certificates independently. Iterative-refinement
literature supports this architecture; additional precision alone is not a
certificate.[^refinement] SoPlex addresses these LP subproblems, not the general
indefinite fixed-order quadratic optimization.

## 6. Acceptance gates

| Gate | Required evidence |
| --- | --- |
| Arithmetic | Exact parsing/serialization round trips, large integers, equivalent row scales, and rejection of malformed fractions; no hidden float conversion |
| Linear algebra | Unique, inconsistent, singular, and nearly singular rational systems; exact residuals and kernel checks; pivot permutations tested |
| Global correctness | Exact equality with `scripts/paper_oracle.py` on small seeded rational cases, including nonsymmetric free-weight bodies; independent planar-area and symmetric-product checks |
| Structural edge cases | Square full support, duplicate-row plateaus, certified redundancy, tied optimizers, exact symplectic transformations, and cyclic/reversal conventions |
| Arithmetic sensitivity | Rational perturbations $0,\pm2^{-k}$ for growing $k$, including distinctions lost in binary64; vary input bit length independently of facet count |
| Certificates and lifecycle | Mutated weights, objectives, row mappings, or coverage rejected; interruption cannot produce exact-complete status; independent verification survives optimizer changes |

Extend the existing optimization tests with rational comparisons: converting the
expected fraction to float is insufficient. Keep tolerances only for comparisons
with the explicitly numerical backend. Parameterized families must remain in
the declared geometric domain, or deliberately test its boundary.

Deliver the exact input/output contract, reference solver, verifier, deterministic
fixtures, and a benchmark report before enabling new pruning. Report enumeration
counts, dimension, facet count, input/intermediate bit lengths, memory, time, and
completion status. No speedup threshold is required for this milestone. General
symbolic parameter families and algebraic-number inputs remain separate work.

## Sources

[^hk]: P. Haim-Kislev, *On the Symplectic Size of Convex Polytopes*, GAFA 29
    (2019), 440–463, Theorem 1.1 and §3.
    [Paper, v3](https://arxiv.org/html/1712.03494v3).

[^lv]: K. Leipold and F. Vallentin, *Computing the EHZ capacity is NP-hard*,
    Proc. AMS Series B 11 (2024), 603–611, equation (1) and Theorem 1.1.
    [Paper, v3](https://arxiv.org/html/2402.09914v3).
    The rationality argument above is a derivation from the formulation,
    not a claim that their NP-hardness theorem proves it.

[^rust]: `num-rational` 0.4.2, official API documentation:
    [BigRational](https://docs.rs/num-rational/0.4.2/num_rational/type.BigRational.html).

[^bareiss]: E. H. Bareiss, *Sylvester's Identity and Multistep Integer-Preserving
    Gaussian Elimination*, Math. Comp. 22 (1968), 565–578.
    [DOI](https://doi.org/10.1090/S0025-5718-1968-0226829-0).

[^cgal]: O. Devillers and S. Schirra, CGAL Manual, *Robustness Issues*,
    especially “The role of predicates and constructions.”
    [Official documentation](https://doc.cgal.org/latest/Manual/devman_robustness.html).

[^soplex]: SoPlex 8.0.2, *How to use SoPlex as an exact LP solver*.
    [Official documentation](https://soplex.zib.de/doc-8.0.2/html/EXACT.php).

[^refinement]: A. M. Gleixner, D. E. Steffy and K. Wolter, *Iterative Refinement
    for Linear Programming*, INFORMS Journal on Computing 28 (2016), 449–464.
    [DOI](https://doi.org/10.1287/ijoc.2016.0692).
