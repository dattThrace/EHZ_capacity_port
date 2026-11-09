"""Research probe, not an EHZ solver or a MATLAB port.

Verify an ordering DP using exact rational arithmetic and exhaustive permutations.
Run: python /workspace/scratch/ehz-research/check_ordering_dp.py
Only the Python standard library is needed.
"""

from fractions import Fraction as F
from itertools import permutations
from math import comb, factorial
from random import Random


def score(w, order):
    return sum((w[order[i]][order[j]] for i in range(len(order))
                for j in range(i)), F(0))


def ordering_dp(w):
    """O(s^2 2^s) arithmetic operations and O(2^s) stored values."""
    s = len(w)
    values = [None] * (1 << s)
    values[0] = F(0)
    for mask in range(1, 1 << s):
        options = []
        for j in range(s):
            if mask & (1 << j):
                rest = mask ^ (1 << j)
                options.append(values[rest] + sum(
                    (w[j][i] for i in range(s) if rest & (1 << i)), F(0)))
        values[mask] = max(options)
    return values[-1]


def omega(a, b):
    n = len(a) // 2
    return sum(a[i] * b[n+i] - a[n+i] * b[i] for i in range(n))


def main():
    rng = Random(20261109)
    arbitrary = 0
    for s in range(2, 8):
        for _ in range(3):
            w = [[F(0) for _ in range(s)] for _ in range(s)]
            for i in range(s):
                for j in range(i):
                    w[i][j] = F(rng.randint(-7, 7), rng.randint(1, 9))
                    w[j][i] = -w[i][j]
            expected = max(score(w, order) for order in permutations(range(s)))
            assert ordering_dp(w) == expected
            arbitrary += 1
    print(f"PASS: {arbitrary} rational skew-weight instances; DP equals exhaustive ordering.")

    # Invertible lower-triangular B, plus negative row sum, gives a bounded
    # full-dimensional simplex {x : A x <= 1}. Its unique beta is 1/(d+1).
    for d in [2, 4, 6, 8]:
        a = [[(1 if i == j else rng.randint(-2, 2) if j < i else 0)
              for j in range(d)] for i in range(d)]
        a.append([-sum(row[j] for row in a) for j in range(d)])
        s = d + 1
        w = [[F(omega(x, y), s*s) for y in a] for x in a]
        assert all(sum(row) == 0 for row in w)
        expected = max(score(w, (0,) + order)
                       for order in permutations(range(1, s)))
        reduced = [[w[i][j] for j in range(1, s)] for i in range(1, s)]
        assert ordering_dp(reduced) == expected
        capacity = 1 / (2 * expected)
        print(f"PASS: R^{d} simplex, {s} facets; anchored DP agrees with "
              f"{factorial(s-1)} cyclic orders; capacity={capacity}.")

    for m in [8, 10, 12, 20]:
        cycles = sum(comb(m, s) * factorial(s-1) for s in range(2, m+1))
        print(f"Complete facet graph, m={m}: {cycles:,} cyclic support orders.")


if __name__ == "__main__":
    main()
