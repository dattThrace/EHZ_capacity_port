"""Exact probe of the H-input symmetric-product formula and its orbit witness.

The lower bound is the published inradius theorem, not established by these tests.
The constructed four-segment orbit verifies attainment with exact fractions.
Uses only the Python standard library. Run this file directly.
"""

from fractions import Fraction as F
from random import Random


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


def formula_and_witness(u, h, v, g):
    # K={q: |u_i.q|<=h_i}, T={p: |v_j.p|<=g_j}.
    m, i, j = max((abs(dot(a,b))/(ha*gb), i, j)
                  for i,(a,ha) in enumerate(zip(u,h))
                  for j,(b,gb) in enumerate(zip(v,g)))
    assert m > 0
    r = 1/m
    sign = 1 if dot(u[i],v[j]) > 0 else -1
    q0 = [r*sign*x/g[j] for x in v[j]]
    p0 = [r*x/h[i] for x in u[i]]
    assert all(abs(dot(a,q0)) <= ha for a,ha in zip(u,h))
    assert all(abs(dot(b,p0)) <= gb for b,gb in zip(v,g))
    assert dot(u[i],q0) == h[i]
    assert dot([sign*x for x in v[j]],p0) == g[j]
    # The orbit follows q-facet +u_i, p-facet -sign*v_j, then opposites.
    # From (q0,p0) to (q0,-p0), (-q0,-p0), (-q0,p0), (q0,p0).
    lambda_q, lambda_p = 2*r/h[i], 2*r/g[j]
    assert [lambda_q*x for x in u[i]] == [2*x for x in p0]
    assert [lambda_p*sign*x for x in v[j]] == [2*x for x in q0]
    action = (2*lambda_q*h[i] + 2*lambda_p*g[j])/2
    assert action == 4*r
    return action


def main():
    basis = [[F(1), F(0)], [F(0), F(1)]]
    diamond = [[F(1), F(1)], [F(1), F(-1)]]
    assert formula_and_witness(basis,[F(1)]*2,diamond,[F(1)]*2) == 4
    assert formula_and_witness(diamond,[F(1)]*2,diamond,[F(1)]*2) == 2
    assert formula_and_witness([[F(1)]],[F(2)],[[F(1)]],[F(3)]) == 24
    print('PASS: README cube x polar, diamond x diamond, and interval x interval.')
    rng = Random(20261109)
    for n in [1,2,3,4]:
        for _ in range(5):
            u = [[F(int(i==j)) for j in range(n)] for i in range(n)]
            v = [[F(int(i==j)) for j in range(n)] for i in range(n)]
            for _ in range(3):
                # All facets have positive heights and the basis rows keep
                # both symmetric polytopes bounded with nonempty interior.
                u.append([F(rng.randint(-3,3)) for _ in range(n)])
                v.append([F(rng.randint(-3,3)) for _ in range(n)])
            h = [F(rng.randint(1,4)) for _ in u]
            g = [F(rng.randint(1,4)) for _ in v]
            formula_and_witness(u,h,v,g)
    print('PASS: 20 seeded rational symmetric products; exact containment, closure, and action witnesses.')


if __name__ == '__main__':
    main()
