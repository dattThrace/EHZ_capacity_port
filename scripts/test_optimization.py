#!/usr/bin/env python3
"""Independent optimization checks; run --help for the explicit test modes.

No check disappears under python -O. Binary calls have hard timeouts.
Exact formula checks deliberately use the existing independent Fraction oracle;
area, product, and physical-itinerary checks use different characterizations.
"""
import argparse
import copy
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations
import json
import math
import os
from pathlib import Path
from random import Random
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "research/probes"))
from check_ordering_dp import ordering_dp, score
from paper_oracle import affine_solution, omega, theorem_capacity

BINARY = None
SQUARE = {"normals": [[1, 0], [0, 1], [-1, 0], [0, -1]], "heights": [1]*4}


def near(actual, expected, rtol=2e-8, atol=0.):
    actual, expected = float(actual), float(expected)
    if not math.isfinite(actual) or not math.isfinite(expected):
        raise ValueError("nonfinite comparison")
    if abs(actual - expected) > atol + rtol * abs(expected):
        raise ValueError(f"{actual} != {expected}")


def validate_witness(data, result):
    """Check printed data in original coordinates, independently of diagnostics.

    This establishes numerical feasibility and self-consistency, NOT optimality.
    Fraction arithmetic avoids introducing summation cancellation in the checker.
    """
    a = [[F(v) for v in row] for row in data["normals"]]
    h = list(map(F, data["heights"]))
    m, d = len(a), len(a[0])
    for prefix, field in (("", "minimum"), ("maximizing_", "maximum")):
        order, weights = result[prefix + "order"], result[prefix + "beta"]
        if (len(order) != len(weights) or len(order) < 2
                or len(set(order)) != len(order)
                or any(type(i) is not int or not 0 <= i < m for i in order)):
            raise ValueError("invalid witness indices")
        if any(not math.isfinite(b) or b < 0 for b in weights):
            raise ValueError("invalid witness weights")
        beta = list(map(F, weights))
        near(sum(b*h[i] for i, b in zip(order, beta)), 1)
        terms = [[b*a[i][j] for i, b in zip(order, beta)] for j in range(d)]
        # Normwise backward error, with no unit-dependent absolute floor.
        # Componentwise relative error wrongly rejects harmless positive
        # roundoff weights in coordinates inactive at a boundary optimum.
        residual = max(abs(sum(row)) for row in terms)
        scale = max(sum(map(abs, row)) for row in terms)
        if residual > F("2e-8") * scale:
            raise ValueError(f"closure failure: residual={float(residual)}, scale={float(scale)}")
        q = sum((beta[i]*beta[j]*omega(a[order[i]], a[order[j]])
                 for i in range(len(order)) for j in range(i)), F(0))
        if (field == "minimum" and q >= 0) or (field == "maximum" and q <= 0):
            raise ValueError("wrong objective orientation")
        near(q, result[field])
        near(result["capacity"], 1/(2*abs(q)))
    near(result["minimum"], -result["maximum"])


def square_witness(weights=(.25, .25, .25, .25)):
    # CCW facet order gives the negative-Q convention used by production.
    order = list(range(4))
    q = sum(weights[i]*weights[j]*omega(SQUARE["normals"][i], SQUARE["normals"][j])
            for i in range(4) for j in range(i))
    return {"capacity": -1/(2*q), "minimum": q, "maximum": -q,
            "order": order, "beta": list(weights),
            "maximizing_order": order[::-1], "maximizing_beta": list(weights)[::-1]}


def simplex(d=4, seed=0):
    rng = Random(seed)
    a = [[int(i == j) if i <= j else rng.randint(-2, 2)
          for j in range(d)] for i in range(d)]
    a.append([-sum(row[j] for row in a) for j in range(d)])
    return {"normals": a, "heights": [1]*(d+1)}


def truncated(height):
    return {"normals": [[-1,0,0,0], [0,-1,0,0], [0,0,-1,0],
                        [0,0,0,-1], [1,1,1,1], [1,1,0,0]],
            "heights": [1,1,1,1,1,height]}


@lru_cache(None)
def exact_capacity(serialized):
    return theorem_capacity(**json.loads(serialized))[0]


def expected(data):
    return exact_capacity(json.dumps(data, sort_keys=True))


def polygon(vertices):
    a, h = [], []
    for u, v in zip(vertices, vertices[1:] + vertices[:1]):
        row = [v[1]-u[1], u[0]-v[0]]
        a.append(row)
        h.append(sum(x*y for x, y in zip(row, u)))
    area = abs(sum(F(u[0]*v[1]-u[1]*v[0], 2)
                   for u, v in zip(vertices, vertices[1:] + vertices[:1])))
    return {"normals": a, "heights": h}, area


def transform(data, kind):
    out = copy.deepcopy(data)
    if kind == "relabel":
        indices = list(range(len(out["normals"])))
        Random(7301).shuffle(indices)
        out = {k: [v[i] for i in indices] for k, v in out.items()}
    elif kind == "row-scale":
        scales = [2.**(8*i-20) for i in range(len(out["normals"]))]
        out["normals"] = [[s*v for v in row] for s, row in zip(scales, out["normals"])]
        out["heights"] = [s*v for s, v in zip(scales, out["heights"])]
    elif kind == "translate":
        # K -> K+t; all new heights remain strictly positive in these fixtures.
        t = [F(1, 32), F(-1, 16), F(1, 32), F(1, 16)]
        out["heights"] = [float(F(h)+sum(F(x)*v for x,v in zip(row,t)))
                          for row,h in zip(out["normals"], out["heights"])]
    elif kind == "shear":
        # x'=Sx, S=(I B;0 I), B=(2 1;1 3); rows transform by S^-1.
        out["normals"] = [[a,b,c-2*a-b,d-a-3*b] for a,b,c,d in out["normals"]]
    elif kind == "antisymplectic":
        out["normals"] = [row[:2] + [-v for v in row[2:]] for row in out["normals"]]
    else:
        raise ValueError(kind)
    return out


class MathematicalTests(unittest.TestCase):
    def test_checker_rejects_corruption(self):
        good = square_witness()
        validate_witness(SQUARE, good)
        for field, value in (("capacity", 2), ("maximum", float("nan")),
                             ("order", [0,1,1,3]), ("order", [0,1,2,4]),
                             ("beta", [.3,.25,.25,.25]),
                             ("maximizing_beta", [-.25,.25,.25,.25])):
            with self.subTest(field=field, value=value):
                bad = copy.deepcopy(good)
                bad[field] = value
                with self.assertRaises(ValueError):
                    validate_witness(SQUARE, bad)

        # Embed the square witness in a box and add an almost-zero weight on
        # an otherwise inactive q2 facet. The checker must resolve roundoff,
        # while still rejecting a significant closure defect in that direction.
        box = {"normals": [[1,0,0,0], [0,0,1,0], [-1,0,0,0], [0,0,-1,0],
                           [0,1,0,0], [0,-1,0,0], [0,0,0,1], [0,0,0,-1]],
               "heights": [1]*8}
        for epsilon in (1e-16, 1e-4):
            b = [.25-epsilon,.25,.25,.25,epsilon]
            q = sum(b[i]*b[j]*omega(box["normals"][i],box["normals"][j])
                    for i in range(5) for j in range(i))
            candidate = {"capacity": -1/(2*q), "minimum": q, "maximum": -q,
                         "order": list(range(5)), "beta": b,
                         "maximizing_order": list(range(4,-1,-1)), "maximizing_beta": b[::-1]}
            if epsilon < 1e-8:
                validate_witness(box, candidate)
            else:
                with self.assertRaises(ValueError):
                    validate_witness(box, candidate)

    def test_feasible_is_not_optimal(self):
        candidate = square_witness((.125, .375, .125, .375))
        validate_witness(SQUARE, candidate)
        self.assertGreater(candidate["capacity"], 4)
        with self.assertRaises(ValueError):
            near(candidate["capacity"], 4)
        # Ordering is solved globally for these frozen feasible weights, yet
        # the capacity is still wrong. Fixed-beta DP is not a general solver.
        b = list(map(F, candidate["beta"]))
        w = [[b[i]*b[j]*omega(x,y) for j,y in enumerate(SQUARE["normals"])]
             for i,x in enumerate(SQUARE["normals"])]
        self.assertEqual(1/(2*ordering_dp(w)), F(16,3))

    def test_dp_against_full_permutations(self):
        rng = Random(20261110)
        for size in range(2, 7):
            for _ in range(3):
                w = [[F(0) for _ in range(size)] for _ in range(size)]
                for i in range(size):
                    for j in range(i):
                        w[i][j] = F(rng.randint(-7,7), rng.randint(1,9))
                        w[j][i] = -w[i][j]
                self.assertEqual(ordering_dp(w), max(score(w,p) for p in permutations(range(size))))

    def test_cyclic_reduction_requires_closure(self):
        data = simplex()
        a, s = data["normals"], len(data["normals"])
        w = [[F(omega(x,y), s*s) for y in a] for x in a]
        self.assertTrue(all(sum(row) == 0 for row in w))
        full = max(score(w,p) for p in permutations(range(s)))
        self.assertEqual(full, ordering_dp([r[1:] for r in w[1:]]))
        p = tuple(range(s))
        self.assertEqual(score(w,p), score(w,p[2:]+p[:2]))
        self.assertEqual(score(w,p), -score(w,p[::-1]))
        bad = [[F(0),F(1),F(0)], [F(-1),F(0),F(0)], [F(0)]*3]
        self.assertNotEqual(ordering_dp(bad), max(score(bad,(0,)+p) for p in permutations((1,2))))

    def test_square_needs_four_directions(self):
        self.assertEqual(expected(SQUARE), 4)
        # Every proper feasible support has Q=0; the d+1 support bound is false.
        for omitted in range(4):
            small = {k: [x for i,x in enumerate(v) if i != omitted] for k,v in SQUARE.items()}
            with self.assertRaises(ValueError):
                theorem_capacity(**small)

    def test_exact_singular_plateau_and_duplicate_row(self):
        duplicate = copy.deepcopy(SQUARE)
        duplicate["normals"].append([2,0])
        duplicate["heights"].append(2)
        self.assertEqual(expected(duplicate), 4)
        # Splitting a weighted normal between duplicate inequalities creates a
        # flat feasible direction; smaller supports must retain its optimum.
        for t in (F(0), F(1,8), F(1,4)):
            beta = [t, F(1,4), F(1,4), F(1,4), (F(1,4)-t)/2]
            order = [0,4,1,2,3]
            q = sum(beta[i]*beta[j]*omega(duplicate["normals"][i], duplicate["normals"][j])
                    for r,i in enumerate(order) for j in order[:r])
            self.assertEqual(q, F(-1,8))

    def test_planar_area_independent_of_stationary_oracle(self):
        for vertices in ([[-2,-1],[3,-1],[-1,3]],
                         [[-2,-1],[2,-2],[3,1],[-1,3]],
                         [[-2,-1],[1,-2],[3,0],[1,3],[-2,2]]):
            data, area = polygon(vertices)
            self.assertEqual(expected(data), area)

    def test_exact_free_weight_fixtures(self):
        # A six-facet R4 body has a genuine free direction, unlike a simplex.
        for height in (1,2):
            data = truncated(height)
            a = list(map(list, zip(*data["normals"]))) + [data["heights"]]
            solution = affine_solution(a, [0,0,0,0,1], 6)
            self.assertEqual(len(solution[1]), 1)
            self.assertEqual(expected(data), F(21,4) if height == 1 else F(6))


class ItineraryTests(unittest.TestCase):
    def test_physical_lp_and_returned_segments(self):
        import numpy as np
        from check_itinerary_lp import itinerary_lp
        for data in (SQUARE, simplex(4, 31), truncated(1)):
            a, h = np.asarray(data["normals"], float), np.asarray(data["heights"], float)
            m, d = a.shape
            best = math.inf
            for tail in permutations(range(1,m)):
                order = (0,)+tail
                result = itinerary_lp(a,h,order)
                if not result.success:
                    self.assertEqual(result.status, 2, result.message)
                    continue
                x = result.x[:d].copy()
                start = x.copy()
                durations = result.x[d:]
                self.assertGreaterEqual(float(min(durations)), -1e-9)
                for i,t in zip(order,durations):
                    self.assertLessEqual(float(max(a@x-h)), 1e-8)
                    near(a[i]@x, h[i], atol=1e-8)
                    x += t*np.r_[a[i,d//2:], -a[i,:d//2]]
                self.assertLessEqual(float(np.linalg.norm(x-start)), 1e-8)
                near(result.fun, sum(durations*h[list(order)])/2)
                self.assertGreater(result.fun, 0)
                best = min(best, result.fun)
            near(best, expected(data))

    def test_unchecked_itinerary_can_accept_constant_loop(self):
        from check_itinerary_lp import itinerary_lp
        result = itinerary_lp(SQUARE["normals"], SQUARE["heights"], (0,1))
        self.assertTrue(result.success)
        near(result.fun, 0, atol=1e-10)
        # Positive normalized closure is impossible on these adjacent normals.
        self.assertIsNone(affine_solution([[1,0],[0,1],[1,1]], [0,0,1], 2))


class ProductionTests(unittest.TestCase):
    def run_solver(self, data, reject=False):
        p = subprocess.run([str(BINARY)], input=json.dumps(data), text=True,
                           capture_output=True, timeout=120,
                           env=dict(os.environ, EHZ_PYTHON=sys.executable))
        if reject:
            self.assertNotEqual(p.returncode, 0, p.stdout)
            self.assertIn("error:", p.stderr)
            return None
        self.assertEqual(p.returncode, 0, p.stderr)
        result = json.loads(p.stdout)["result"]
        validate_witness(data, result)
        return result

    def check_capacity(self, data, target):
        result = self.run_solver(data)
        near(result["capacity"], target)
        # Witness and optimum comparisons are separate assertions.
        near(result["maximum"], 1/(2*F(target)))
        return result

    def test_exact_quality_and_metamorphic_cases(self):
        for data in (simplex(4,11), simplex(4,32), truncated(1), truncated(2)):
            target = expected(data)
            self.check_capacity(data, target)
            for kind in ("relabel", "row-scale", "translate", "shear", "antisymplectic"):
                with self.subTest(data=data, transform=kind):
                    self.check_capacity(transform(data,kind), target)

    def test_singular_redundancy_and_full_support(self):
        r = self.check_capacity(SQUARE, 4)
        self.assertEqual(len(r["order"]), 4)
        self.assertEqual(r["supports"], 11)
        duplicate = copy.deepcopy(SQUARE)
        duplicate["normals"].append([2,0])
        duplicate["heights"].append(2)
        self.check_capacity(duplicate, 4)
        redundant = copy.deepcopy(SQUARE)
        redundant["normals"].append([1,1])
        redundant["heights"].append(3)
        self.check_capacity(redundant, 4)

    def test_product_formula_and_multiple_optima(self):
        # Symmetric products: independent inradius/pairing identity.
        for u,v,h,g in (([[1,0],[0,1]], [[1,0],[0,1]], [1,1], [1,1]),
                        ([[1,1],[1,-1]], [[1,1],[1,-1]], [1,1], [1,1]),
                        ([[1,0],[0,1]], [[1,2],[2,-1]], [2,3], [3,2])):
            reps = [r+[0,0] for r in u] + [[0,0]+r for r in v]
            data = {"normals": reps+[[-x for x in r] for r in reps], "heights": (h+g)*2}
            pairing = max(F(abs(sum(x*y for x,y in zip(a,b))),ha*gb)
                          for a,ha in zip(u,h) for b,gb in zip(v,g))
            self.check_capacity(data, 4/pairing)

    def test_scale_and_conditioning_are_different(self):
        for scale in (2.**-30, 1., 2.**30):
            data = copy.deepcopy(SQUARE)
            data["heights"] = [scale]*4
            self.check_capacity(data, 4*scale*scale)
        for stretch in (1., 8., 64.):
            # diag(k,1/k) is symplectic: increasing conditioning, fixed capacity.
            data = {"normals": [[1/stretch,0],[0,stretch],[-1/stretch,0],[0,-stretch]],
                    "heights": [1]*4}
            self.check_capacity(data, 4)

    def test_rejects_domain_errors(self):
        for data in ({"normals": [[1,0],[-1,0],[0,1]], "heights": [1]*3},
                     {"normals": [[1,0],[-1,0],[0,0]], "heights": [1]*3},
                     {"normals": SQUARE["normals"], "heights": [1,1,0,1]},
                     {"normals": [[1,0,0]]*4, "heights": [1]*4},
                     {"normals": [[1,0],[-1,0],[0,1],[0,-1,0]], "heights": [1]*4}):
            with self.subTest(data=data):
                self.run_solver(data, reject=True)


class ExactProductionTests(unittest.TestCase):
    """Exercise the same optimization gates without rounding expected fractions."""

    def run_exact(self, data, *limits):
        # Fraction(float) means that float's exact binary value. Encoding it as
        # a rational string prevents the production decoder seeing a float.
        encoded = {"schema": "ehz-exact-h-v1",
                   "normals": [[str(F(x)) for x in row] for row in data["normals"]],
                   "heights": [str(F(x)) for x in data["heights"]]}
        process = subprocess.run([str(BINARY), "--mode", "exact", *limits],
                                 input=json.dumps(encoded), text=True,
                                 capture_output=True, timeout=120)
        self.assertEqual(process.returncode, 0, process.stderr)
        result = json.loads(process.stdout)["result"]
        # Independent domain and witness checks use original coordinates; full
        # global comparison below uses the separately derived Fraction oracle.
        from verify_exact import verify
        verified = verify({"result": result}, original=encoded, witness_only=True)
        self.assertTrue(verified["verified"])
        self.assertFalse(verified["global_optimality_verified"])
        return result

    def check_exact(self, data, target):
        result = self.run_exact(data)
        self.assertEqual(result["status"], "exact-complete")
        self.assertTrue(result["completion"]["complete"])
        self.assertEqual(F(result["capacity"]), F(target))
        self.assertEqual(F(result["witness"]["objective"]), 1 / (2 * F(target)))
        self.assertEqual(F(result["witness"]["capacity_upper_bound"]), F(target))
        return result

    def test_rational_oracle_and_metamorphic_cases(self):
        for data in (simplex(4, 11), simplex(4, 32), truncated(1), truncated(2)):
            target = expected(data)
            self.check_exact(data, target)
            for kind in ("relabel", "row-scale", "translate", "shear", "antisymplectic"):
                with self.subTest(data=data, transform=kind):
                    self.check_exact(transform(data, kind), target)

    def test_planar_area_and_boundary_supports(self):
        for vertices in ([[-2, -1], [3, -1], [-1, 3]],
                         [[-2, -1], [2, -2], [3, 1], [-1, 3]],
                         [[-2, -1], [1, -2], [3, 0], [1, 3], [-2, 2]]):
            data, area = polygon(vertices)
            self.check_exact(data, area)
        square = self.check_exact(SQUARE, F(4))
        self.assertEqual(len(square["witness"]["order"]), 4)
        self.assertEqual(square["stats"]["supports"], 11)
        for row, height in (([2, 0], 2), ([1, 1], 2), ([1, 1], 3)):
            data = copy.deepcopy(SQUARE)
            data["normals"].append(row)
            data["heights"].append(height)
            result = self.check_exact(data, F(4))
            self.assertEqual(result["domain"]["retained"], [0, 1, 2, 3])

    def test_exact_bit_length_and_perturbation_sensitivity(self):
        for k in (20, 60, 160):
            for sign in (-1, 0, 1):
                with self.subTest(bits=k, sign=sign):
                    epsilon = F(sign, 2**k)
                    data = copy.deepcopy(SQUARE)
                    data["heights"][0] = 1 + epsilon
                    # The planar rectangle has width 2+epsilon, height 2.
                    self.check_exact(data, 4 + 2 * epsilon)
        data = copy.deepcopy(SQUARE)
        scales = [F(2**180, 3**70), F(5**90, 7**60), F(11**80, 13**70), F(17**50, 19**40)]
        data["normals"] = [[scale * x for x in row]
                           for scale, row in zip(scales, data["normals"])]
        data["heights"] = scales
        self.check_exact(data, F(4))

    def test_exact_cyclic_and_reversal_conventions(self):
        result = self.check_exact(truncated(1), F(21, 4))
        data, witness = truncated(1), result["witness"]
        order, beta = witness["order"], list(map(F, witness["weights"]))
        def objective(route, weights):
            return sum((weights[i] * weights[j] * omega(data["normals"][route[i]],
                                                        data["normals"][route[j]])
                        for i in range(len(route)) for j in range(i)), F(0))
        value = F(witness["objective"])
        for shift in range(len(order)):
            self.assertEqual(objective(order[shift:] + order[:shift],
                                       beta[shift:] + beta[:shift]), value)
        self.assertEqual(objective(order[::-1], beta[::-1]), -value)

    def test_r4_symmetric_product_pairing_formula(self):
        # Independent inradius/pairing identity for centrally symmetric
        # Lagrangian products, including unequal heights and tilted factors.
        u, v, h, g = [[1, 0], [0, 1]], [[1, 2], [2, -1]], [2, 3], [3, 2]
        reps = [row + [0, 0] for row in u] + [[0, 0] + row for row in v]
        data = {"normals": reps + [[-x for x in row] for row in reps],
                "heights": (h + g) * 2}
        pairing = max(F(abs(sum(x * y for x, y in zip(a, b))), ha * gb)
                      for a, ha in zip(u, h) for b, gb in zip(v, g))
        self.check_exact(data, F(4) / pairing)

    def test_limits_preserve_exact_witness_without_complete_claim(self):
        empty = self.run_exact(SQUARE, "--max-steps", "0")
        self.assertEqual(empty["status"], "incomplete")
        self.assertIsNone(empty["capacity"])
        self.assertFalse(empty["completion"]["complete"])
        self.assertIsNone(empty["witness"])
        data = truncated(1)
        complete = self.check_exact(data, F(21, 4))
        # Stop before the final order; this free-weight fixture already has
        # positive witnesses. Work counts choose a stop, never prove coverage.
        steps = str(complete["stats"]["steps"] - 1)
        stopped = self.run_exact(data, "--max-steps", steps)
        self.assertFalse(stopped["completion"]["complete"])
        self.assertIsNone(stopped["capacity"])
        self.assertIsNotNone(stopped["witness"])
        self.assertEqual(stopped["status"], "exact-witness/incomplete")
        interval = self.run_exact(data, "--max-steps", steps, "--certify-interval")
        self.assertEqual(interval["status"], "certified-interval")
        self.assertLessEqual(F(interval["interval"]["lower"]), F(21, 4))
        self.assertGreaterEqual(F(interval["interval"]["upper"]), F(21, 4))
        self.assertFalse(interval["completion"]["complete"])


def main():
    global BINARY
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--math-only", action="store_true", help="exact/checker/LP tests; explicitly excludes Rust")
    mode.add_argument("--binary", type=Path, help="required existing Rust CLI; runs all tests")
    parser.add_argument("--exact-only", action="store_true",
                        help="with --binary, run only the exact production comparisons")
    args = parser.parse_args()
    if args.exact_only and not args.binary:
        parser.error("--exact-only requires --binary")
    classes = [MathematicalTests, ItineraryTests]
    if args.binary:
        BINARY = args.binary.resolve()
        if not BINARY.is_file():
            parser.error(f"binary does not exist: {BINARY}")
        classes = [ExactProductionTests] if args.exact_only else classes + [ProductionTests, ExactProductionTests]
    suite = unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromTestCase(c) for c in classes)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())
