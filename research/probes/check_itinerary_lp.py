"""Research probe for a fixed-itinerary LP, not a general capacity solver.

Run: python /workspace/scratch/ehz-research/check_itinerary_lp.py
Requires the environment's existing NumPy and SciPy.
Exact ordering values come from check_ordering_dp.py; LPs are floating point.
"""

from itertools import permutations
from random import Random
from fractions import Fraction as F

import numpy as np
from scipy.optimize import linprog

from check_ordering_dp import omega, ordering_dp


def itinerary_lp(a, h, order):
    """Requires the selected normals to admit normalized nonnegative closure.

    This condition excludes the constant loops that would otherwise give zero
    objective when all selected facets share a point.
    """
    a, h = np.asarray(a, float), np.asarray(h, float)
    m, d = a.shape
    n, s = d // 2, len(order)
    normals = a[list(order)]
    velocities = np.column_stack((normals[:, n:], -normals[:, :n]))
    # Variables: starting point x0 and nonnegative segment coefficients lambda.
    # x_r = x0 + sum_{t<r} lambda_t * J a_order[t].
    ub, ub_rhs, eq, eq_rhs = [], [], [], []
    for r, facet in enumerate(order):
        prefix = np.zeros((d, s))
        prefix[:, :r] = velocities[:r].T
        point = np.column_stack((np.eye(d), prefix))
        ub.extend(a @ point)
        ub_rhs.extend(h)
        eq.append(a[facet] @ point)
        eq_rhs.append(h[facet])
    eq.extend(np.column_stack((np.zeros((d, d)), velocities.T)))
    eq_rhs.extend(np.zeros(d))
    c = np.r_[np.zeros(d), h[list(order)] / 2]
    return linprog(c, A_ub=ub, b_ub=ub_rhs, A_eq=eq, b_eq=eq_rhs,
                   bounds=[(None, None)] * d + [(0, None)] * s,
                   method="highs")


def minimum_for_support(a, h):
    feasible, best = 0, float("inf")
    s = len(a)
    for rest in permutations(range(1, s)):
        result = itinerary_lp(a, h, (0,) + rest)
        if result.success:
            assert result.fun > 0
            best = min(best, result.fun)
            feasible += 1
        else:
            assert result.status == 2, result.message  # solver reports infeasible
    return best, feasible


def main():
    rng = Random(20261109)
    for d in [2, 4, 6]:
        a = [[(1 if i == j else rng.randint(-2, 2) if j < i else 0)
              for j in range(d)] for i in range(d)]
        a.append([-sum(row[j] for row in a) for j in range(d)])
        s = d + 1
        w = [[F(omega(x, y), s*s) for y in a[1:]] for x in a[1:]]
        expected = 1 / (2 * ordering_dp(w))
        best, feasible = minimum_for_support(a, np.ones(s))
        assert np.isclose(best, float(expected), rtol=1e-9, atol=1e-9)
        print(f"PASS: R^{d} simplex; itinerary LP minimum={best:.12g}; "
              f"exact ordering capacity={expected}; feasible cyclic orders={feasible}.")

    a = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]])
    best, feasible = minimum_for_support(a, np.ones(4))
    assert np.isclose(best, 4)
    print(f"PASS: [-1,1]^2, capacity={best:g}; feasible cyclic orders={feasible}.")

    # README's K x T = [-1,1]^2 x {p: |p1|+|p2| <= 1}.
    # Four selected facets support an orbit of action 4; the known polar-product
    # identity supplies the independent lower bound 4. This is not an exhaustive
    # search of the eight-facet product or a general proof of the LP algorithm.
    a = np.array([[1,0,0,0], [-1,0,0,0], [0,1,0,0], [0,-1,0,0],
                  [0,0,1,1], [0,0,-1,-1], [0,0,1,-1], [0,0,-1,1]])
    order = (0, 5, 1, 4)
    result = itinerary_lp(a, np.ones(8), order)
    assert result.success and np.isclose(result.fun, 4)
    print(f"PASS: README product, selected four-facet itinerary has action={result.fun:g}.")
    opposite = itinerary_lp(a, np.ones(8), (0, 4, 1, 5))
    assert opposite.status == 2
    print("PASS: the reversed README itinerary is infeasible for the fixed flow orientation.")


if __name__ == "__main__":
    main()
