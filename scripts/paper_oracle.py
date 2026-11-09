#!/usr/bin/env python3
"""Small exact oracles derived from Haim-Kislev, arXiv:1712.03494v3.

Theorem 1.1: MAXIMIZE the ordered quadratic sum on M(K).
Corollary 1.2: independent signed absolute-weight formula for symmetric bodies.
No Rust code, MATLAB code, floating rank decisions, or research probe is used.
Exponential reference only; all arithmetic is fractions.Fraction.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product


def affine_solution(a, b, variables):
    """Exact RREF: return a particular solution and a nullspace basis."""
    rows = [[F(v) for v in row] + [F(rhs)] for row, rhs in zip(a, b)]
    pivots = []
    for col in range(variables):
        pivot = next((i for i in range(len(pivots), len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        r = len(pivots)
        rows[r], rows[pivot] = rows[pivot], rows[r]
        scale = rows[r][col]
        rows[r] = [v / scale for v in rows[r]]
        for i in range(len(rows)):
            if i != r:
                scale = rows[i][col]
                rows[i] = [x - scale*y for x, y in zip(rows[i], rows[r])]
        pivots.append(col)
    if any(not any(row[:variables]) and row[-1] for row in rows):
        return None
    point = [F(0)] * variables
    for r, col in enumerate(pivots):
        point[col] = rows[r][-1]
    basis = []
    for col in range(variables):
        if col not in pivots:
            v = [F(0)] * variables
            v[col] = F(1)
            for r, pivot in enumerate(pivots):
                v[pivot] = -rows[r][col]
            basis.append(v)
    return point, basis


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


def matvec(a, x):
    return [dot(row, x) for row in a]


def omega(a, b):
    n = len(a)//2
    return sum((a[i]*b[n+i] - a[n+i]*b[i] for i in range(n)), F(0))


def stationary(h, point, basis):
    """Differentiate Q(point + B z); solve in free coordinates, not a KKT block."""
    if not basis:
        return point
    hb = [matvec(h, v) for v in basis]
    hp = matvec(h, point)
    reduced = [[dot(u, hv) for hv in hb] for u in basis]
    rhs = [-dot(u, hp) for u in basis]
    solution = affine_solution(reduced, rhs, len(basis))
    # At a singular stationary extremum, a flat null direction reaches a boundary
    # with the same value, which is covered by a smaller support.
    if solution is None or solution[1]:
        return None
    z = solution[0]
    return [point[i] + sum((z[j]*basis[j][i] for j in range(len(basis))), F(0))
            for i in range(len(point))]


def maximum_on_orders(normals, heights, closure=True, mass=F(1)):
    m, d = len(normals), len(normals[0])
    best = F(0)
    orders = 0
    for size in range(1, m+1):
        for support in combinations(range(m), size):
            a = [[normals[i][j] for i in support] for j in range(d)] if closure else []
            a.append([heights[i] for i in support])
            solution = affine_solution(a, ([F(0)]*d if closure else []) + [mass], size)
            if solution is None:
                continue
            # Full permutations, deliberately without the production cyclic anchor.
            for route in permutations(range(size)):
                orders += 1
                h = [[F(0) for _ in support] for _ in support]
                for i in range(size):
                    for j in range(i):
                        u, v = route[i], route[j]
                        h[u][v] = h[v][u] = omega(normals[support[u]], normals[support[v]])
                beta = stationary(h, *solution)
                if beta is not None and min(beta) >= 0:
                    best = max(best, dot(beta, matvec(h, beta))/2)
    return best, orders


def theorem_capacity(normals, heights):
    normals = [[F(v) for v in row] for row in normals]
    maximum, orders = maximum_on_orders(normals, list(map(F, heights)))
    if maximum <= 0:
        raise ValueError("no positive exact maximum")
    return F(1)/(2*maximum), orders


def symmetric_capacity(representatives, heights):
    """Corollary 1.2: sum |beta_i| h_i = 1/2, c = 1/(4 max Q)."""
    representatives = [[F(v) for v in row] for row in representatives]
    best, orders = F(0), 0
    for signs in product((-1, 1), repeat=len(representatives)):
        signed = [[sign*v for v in row] for sign, row in zip(signs, representatives)]
        maximum, count = maximum_on_orders(signed, list(map(F, heights)), closure=False, mass=F(1, 2))
        best = max(best, maximum)
        orders += count
    if best <= 0:
        raise ValueError("no positive exact symmetric maximum")
    return F(1)/(4*best), orders
