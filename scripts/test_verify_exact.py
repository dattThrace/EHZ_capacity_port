#!/usr/bin/env python3
"""Independent verifier lifecycle and adversarial certificate tests (stdlib)."""
import copy
import json
import os
from pathlib import Path
import subprocess
from fractions import Fraction as F
import unittest

from verify_exact import VerificationError, canonical, rational, verify


def square_certificate():
    return {"result": {
        "schema": "ehz-exact-result-v1", "status": "exact-complete",
        "input": {"schema": "ehz-exact-h-v1",
                  "normals": [["1", "0"], ["-1", "0"], ["0", "1"], ["0", "-1"]],
                  "heights": ["1"] * 4},
        "domain": {
            "normals": [["1", "0"], ["-1", "0"], ["0", "1"], ["0", "-1"]],
            "heights": ["1"] * 4,
            "rows": [{"original": i, "scale": "1", "representative": i} for i in range(4)],
            "retained": [0, 1, 2, 3], "removals": [], "bounded_weights": ["1/4"] * 4},
        "witness": {"order": [0, 3, 1, 2], "weights": ["1/4"] * 4,
                    "objective": "1/8", "capacity_upper_bound": "4"},
        "capacity": "4", "completion": {"complete": True, "reason": "exhaustive"},
        "stats": {"orders": 1},  # Deliberately not used to establish coverage.
    }}


class VerifierTests(unittest.TestCase):
    def test_complete_requires_independent_enumeration(self):
        report = verify(square_certificate())
        self.assertEqual(report["status"], "exact-complete")
        self.assertTrue(report["global_optimality_verified"])
        self.assertGreater(report["independent_orders"], 1)
        witness = verify(square_certificate(), witness_only=True)
        self.assertEqual(witness["status"], "exact-witness/incomplete")
        self.assertFalse(witness["global_optimality_verified"])

    def test_canonical_exact_arithmetic(self):
        for value in ("0", "123456789012345678901234567890", "-7/13", "0.125", ".5", "1e-40", "1/-2"):
            number = rational(value)
            self.assertEqual(rational(canonical(number), True), number)
        self.assertEqual(rational(["10", "20"]), F(1, 2))
        for malformed in (0.25, 1, True, "1/0", " 1", "1_000", ["1", "0"], [1, 2], "nan"):
            with self.subTest(value=malformed), self.assertRaises(VerificationError):
                rational(malformed)
        with self.assertRaises(VerificationError):
            rational("2/4", True)

    def test_arbitrary_size_certificate_integers(self):
        # Exceeds Python 3.11's default 4,300 digit decimal conversion ceiling.
        huge = "1" + "0" * 5000 + "3"
        integer = rational(huge, True)
        self.assertEqual(canonical(integer), huge)
        fraction = rational(huge + "/7", True)
        self.assertEqual(canonical(fraction), huge + "/7")
        self.assertEqual(rational(canonical(fraction), True), fraction)
        self.assertEqual(rational([huge, "7"]), fraction)

    def test_source_parser_limits_and_unknown_fields_match_solver(self):
        for value in ("1e1000001", "1e-1000001", "1.1e-1000000"):
            with self.subTest(value=value), self.assertRaisesRegex(VerificationError, "exponent"):
                rational(value)
        # Certificate canonical fractions have no decimal expansion limit.
        with self.assertRaisesRegex(VerificationError, "noncanonical"):
            rational("1e1000001", True)
        cert = square_certificate()
        original = copy.deepcopy(cert["result"]["input"])
        original["ignored"] = "not part of exact input v1"
        with self.assertRaisesRegex(VerificationError, "input fields"):
            verify(cert, original=original)
        cert["result"]["input"]["ignored"] = "not part of exact input v1"
        with self.assertRaisesRegex(VerificationError, "input fields"):
            verify(cert)

    def test_mutated_witness_and_mapping_rejected(self):
        mutations = [
            lambda r: r["witness"]["weights"].__setitem__(0, "1/3"),
            lambda r: r["witness"].__setitem__("objective", "1/7"),
            lambda r: r["witness"].__setitem__("capacity_upper_bound", "5"),
            lambda r: r["witness"].__setitem__("order", [0, 1, 2, 3]),
            lambda r: r["domain"]["rows"][0].__setitem__("scale", "2"),
            lambda r: r["domain"]["rows"][0].__setitem__("representative", 1),
            lambda r: r["domain"]["rows"][0].__setitem__("original", 1),
            lambda r: r["domain"]["bounded_weights"].__setitem__(0, "0"),
            lambda r: r["completion"].__setitem__("complete", False),
            lambda r: r.__setitem__("capacity", "5"),
        ]
        for mutate in mutations:
            cert = square_certificate()
            mutate(cert["result"])
            with self.subTest(mutation=mutate), self.assertRaises(VerificationError):
                verify(cert)

    def test_feasible_suboptimal_witness_does_not_establish_global_optimality(self):
        cert = square_certificate()
        result = cert["result"]
        result["witness"].update(weights=["1/8", "3/8", "1/8", "3/8"],
                                 objective="3/32", capacity_upper_bound="16/3")
        result["capacity"] = "16/3"
        with self.assertRaisesRegex(VerificationError, "independent enumeration"):
            verify(cert)
        self.assertFalse(verify(cert, witness_only=True)["global_optimality_verified"])
        result.update(status="exact-witness/incomplete", capacity=None,
                      completion={"complete": False, "reason": "interrupted"})
        self.assertEqual(verify(cert)["capacity_upper_bound"], "16/3")

    def test_input_binding_uses_exact_values(self):
        cert = square_certificate()
        original = copy.deepcopy(cert["result"]["input"])
        original["heights"][0] = ["2", "2"]
        self.assertTrue(verify(cert, original=original)["verified"])
        original["heights"][0] = "1.000000000000000000000000000001"
        with self.assertRaisesRegex(VerificationError, "input identity"):
            verify(cert, original=original)

    def test_duplicate_and_certified_redundancy(self):
        cert = square_certificate()
        result = cert["result"]
        result["input"]["normals"] += [["2", "0"], ["1", "1"]]
        result["input"]["heights"] += ["2", "3"]
        domain = result["domain"]
        domain["rows"] += [{"original": 4, "scale": "1/2", "representative": 0},
                           {"original": 5, "scale": "1", "representative": 5}]
        domain["bounded_weights"] = ["1/10", "2/5", "1/10", "1/5", "1/10", "1/10"]
        domain["removals"] = [{"row": 5, "using": [0, 2], "multipliers": ["1", "1"]}]
        self.assertTrue(verify(cert)["global_optimality_verified"])
        domain["removals"][0]["multipliers"][0] = "2"
        with self.assertRaisesRegex(VerificationError, "redundancy normal"):
            verify(cert)

    def test_unjustified_retention_is_rejected_geometrically(self):
        cert = square_certificate()
        result = cert["result"]
        result["input"]["normals"] += [["1", "1"]]
        result["input"]["heights"] += ["3"]
        domain = result["domain"]
        domain["normals"] += [["1", "1"]]
        domain["heights"] += ["3"]
        domain["retained"] += [4]
        domain["rows"] += [{"original": 4, "scale": "1", "representative": 4}]
        domain["bounded_weights"] = ["1/7", "2/7", "1/7", "2/7", "1/7"]
        with self.assertRaisesRegex(VerificationError, "facet"):
            verify(cert)

    def test_no_incumbent_no_capacity(self):
        cert = square_certificate()
        result = cert["result"]
        result.update(status="incomplete", witness=None, capacity=None, domain=None,
                      completion={"complete": False, "reason": "domain budget"})
        report = verify(cert)
        self.assertFalse(report["domain_verified"])
        self.assertFalse(report["global_optimality_verified"])
        result["capacity"] = "4"
        with self.assertRaises(VerificationError):
            verify(cert)


class IntervalTests(unittest.TestCase):
    def certificate(self):
        cert = square_certificate()
        cert["result"].update(status="certified-interval", capacity=None,
                              completion={"complete": False, "reason": "interrupted"},
                              interval={"lower": "1", "upper": "4",
                                        "objective_upper_bound": "1/2",
                                        "method": "pairing-envelope-v1"})
        return cert

    def test_independent_interval_bound(self):
        report = verify(self.certificate())
        self.assertEqual(report["status"], "certified-interval")
        self.assertEqual(report["capacity_interval"], ["1", "4"])
        self.assertFalse(report["global_optimality_verified"])
        self.assertEqual(verify(self.certificate(), witness_only=True)["capacity_interval"], ["1", "4"])

    def test_mutated_bounds_and_method_rejected(self):
        for field, value in [("lower", "2"), ("upper", "3"),
                             ("objective_upper_bound", "1/4"), ("method", "trusted-counter")]:
            cert = self.certificate()
            cert["result"]["interval"][field] = value
            with self.subTest(field=field), self.assertRaises(VerificationError):
                verify(cert)
        for field, value in [("status", "exact-witness/incomplete"), ("interval", None)]:
            cert = self.certificate()
            cert["result"][field] = value
            with self.subTest(field=field), self.assertRaises(VerificationError):
                verify(cert)


ROOT = Path(__file__).resolve().parents[1]
BINARY = Path(os.environ.get("EHZ_EXACT_BINARY", ROOT / "target/release/ehz-capacity"))


@unittest.skipUnless(BINARY.is_file(), "build the release solver for certificate integration checks")
class SolverCertificateTests(unittest.TestCase):
    def solve(self, name="square", flags=()):
        original = json.loads((ROOT / "examples/exact" / (name + ".json")).read_text())
        process = subprocess.run([str(BINARY), *flags], input=json.dumps(original), text=True,
                                 capture_output=True, check=True)
        return json.loads(process.stdout), original

    def test_actual_domain_certificates(self):
        for name in ("square", "decimal_square", "duplicate_square", "redundant_square", "large_row_scales"):
            with self.subTest(name=name):
                cert, original = self.solve(name)
                self.assertTrue(verify(cert, original=original)["global_optimality_verified"])

    def test_actual_arbitrary_size_certificate(self):
        original = square_certificate()["result"]["input"]
        scale = int("1" + "0" * 5000 + "3")
        original["normals"] = [[str(int(x) * scale) for x in row] for row in original["normals"]]
        original["heights"] = [str(scale)] * 4
        process = subprocess.run([str(BINARY)], input=json.dumps(original), text=True,
                                 capture_output=True, check=True)
        report = verify(json.loads(process.stdout), original=original)
        self.assertEqual(report["capacity"], "4")
        self.assertTrue(report["global_optimality_verified"])

    def test_actual_interruption(self):
        cert, original = self.solve(flags=("--max-steps", "0"))
        self.assertEqual(verify(cert, original=original)["status"], "incomplete")
        complete, _ = self.solve()
        steps = complete["result"]["stats"]["steps"]
        cert, _ = self.solve(flags=("--max-steps", str(steps - 1)))
        self.assertEqual(verify(cert)["status"], "exact-witness/incomplete")
        self.assertFalse(verify(cert)["global_optimality_verified"])

    def test_actual_certified_interval(self):
        complete, _ = self.solve()
        steps = complete["result"]["stats"]["steps"]
        cert, original = self.solve(flags=("--max-steps", str(steps - 1), "--certify-interval"))
        report = verify(cert, original=original)
        self.assertEqual(report["status"], "certified-interval")
        self.assertFalse(report["global_optimality_verified"])


if __name__ == "__main__":
    unittest.main()
