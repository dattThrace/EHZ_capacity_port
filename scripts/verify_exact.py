#!/usr/bin/env python3
"""Verify rational solver certificates without using the Rust optimizer.

A positive ordered witness certifies only an upper bound on capacity. The
exact-complete claim is checked by an independent, full-permutation search in
paper_oracle.py. Search counts and the producer's completion flag are never
used as a proof of coverage. Intended for small H-polytopes only.
"""
import argparse
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
import json
from math import gcd, lcm
from pathlib import Path
import re
import sys

# Exact certificates have arbitrary-size integer strings. Python 3.11's
# process-wide decimal conversion ceiling is inappropriate for this verifier.
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

from paper_oracle import affine_solution, maximum_on_orders, omega


class VerificationError(ValueError):
    """The supplied data do not establish the claimed result."""


def require(condition, message):
    if not condition:
        raise VerificationError(message)


def canonical(value):
    return str(value.numerator) if value.denominator == 1 else str(value)


def rational(value, canonical_only=False):
    """Parse the exact JSON scalar contract, never passing through a float."""
    if isinstance(value, list):
        require(not canonical_only and len(value) == 2, "invalid fraction pair")
        require(all(isinstance(x, str) and re.fullmatch(r"[+-]?\d+", x, re.ASCII)
                    for x in value), "fraction pairs require integer strings")
        require(int(value[1]) != 0, "zero denominator")
        return F(int(value[0]), int(value[1]))
    require(isinstance(value, str), "rational must be a string or integer-string pair")
    grammar = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    require(re.fullmatch(grammar, value, re.ASCII) is not None or
            re.fullmatch(r"[+-]?\d+/[+-]?\d+", value, re.ASCII) is not None,
            "malformed exact rational")
    if canonical_only:
        require(re.fullmatch(r"(?:0|-?[1-9]\d*)(?:/[1-9]\d*)?", value, re.ASCII) is not None,
                "noncanonical certificate fraction")
    elif "/" not in value:
        mantissa, separator, exponent = value.lower().partition("e")
        exponent = int(exponent) if separator else 0
        places = len(mantissa.partition(".")[2])
        require(abs(exponent - places) <= 1_000_000, "decimal exponent too large")
    try:
        result = F(*map(int, value.split("/"))) if "/" in value else F(value)
    except (ValueError, ZeroDivisionError) as error:
        raise VerificationError("malformed exact rational") from error
    require(not canonical_only or value == canonical(result), "noncanonical certificate fraction")
    return result


def vector(value, canonical_only=True):
    require(isinstance(value, list), "expected a rational vector")
    return [rational(x, canonical_only) for x in value]


def matrix(value, canonical_only=True):
    require(isinstance(value, list), "expected a rational matrix")
    return [vector(row, canonical_only) for row in value]


def exact_input(value, canonical_only=False):
    require(isinstance(value, dict) and value.get("schema") == "ehz-exact-h-v1",
            "unsupported exact input schema")
    require(set(value) == {"schema", "normals", "heights"}, "unexpected exact input fields")
    normals = matrix(value.get("normals"), canonical_only)
    heights = vector(value.get("heights"), canonical_only)
    d = len(normals[0]) if normals else 0
    require(d > 0 and d % 2 == 0, "dimension must be positive and even")
    require(len(normals) == len(heights), "normal/height count mismatch")
    require(all(len(row) == d and any(row) for row in normals), "invalid normal rows")
    require(all(h > 0 for h in heights), "heights must be strictly positive")
    return normals, heights


def rank(rows, columns):
    solution = affine_solution(rows, [F(0)] * len(rows), columns)
    return columns - len(solution[1])


def row_scale(normal, height):
    values = normal + [height]
    denominator = lcm(*(v.denominator for v in values))
    divisor = reduce(gcd, (abs(v.numerator * (denominator // v.denominator)) for v in values))
    return F(denominator, divisor)


def indices(value, total, name):
    require(isinstance(value, list) and all(type(i) is int and 0 <= i < total for i in value),
            f"invalid {name} indices")
    require(len(set(value)) == len(value), f"duplicate {name} indices")
    return value


def verify_domain(domain, normals, heights):
    """Verify row identities and Farkas deletions, then establish genuine facets."""
    require(isinstance(domain, dict), "missing domain certificate")
    m, d = len(normals), len(normals[0])
    weights = vector(domain.get("bounded_weights"))
    require(len(weights) == m and all(w > 0 for w in weights) and sum(weights) == 1,
            "invalid positive-spanning weights")
    require(all(sum(weights[i] * normals[i][j] for i in range(m)) == 0 for j in range(d)),
            "positive-spanning closure failed")
    require(rank(normals, d) == d, "normals lack full column rank")

    rows = domain.get("rows")
    require(isinstance(rows, list) and len(rows) == m, "incomplete row mapping")
    normalized, normalized_heights, first, active = [], [], {}, []
    for i, entry in enumerate(rows):
        require(isinstance(entry, dict) and type(entry.get("original")) is int and
                entry["original"] == i, "invalid original row identity")
        scale = rational(entry.get("scale"), True)
        require(scale > 0 and scale == row_scale(normals[i], heights[i]),
                "incorrect positive primitive row scale")
        row, height = [scale * x for x in normals[i]], scale * heights[i]
        key = tuple(row + [height])
        representative = first.setdefault(key, i)
        require(type(entry.get("representative")) is int and
                entry["representative"] == representative, "incorrect duplicate row representative")
        normalized.append(row)
        normalized_heights.append(height)
        if representative == i:
            active.append(i)

    removals = domain.get("removals")
    require(isinstance(removals, list), "missing redundancy certificates")
    for removal in removals:
        require(isinstance(removal, dict), "invalid redundancy certificate")
        row = removal.get("row")
        require(type(row) is int and row in active, "invalid removed row identity")
        using = indices(removal.get("using"), m, "redundancy")
        require(using and row not in using and all(i in active for i in using),
                "redundancy certificate references unavailable rows")
        multipliers = vector(removal.get("multipliers"))
        require(len(multipliers) == len(using) and all(v >= 0 for v in multipliers),
                "invalid redundancy multipliers")
        require(all(sum(v * normals[i][j] for i, v in zip(using, multipliers)) == normals[row][j]
                    for j in range(d)), "redundancy normal identity failed")
        require(sum(v * heights[i] for i, v in zip(using, multipliers)) <= heights[row],
                "redundancy height inequality failed")
        active.remove(row)
    retained = indices(domain.get("retained"), m, "retained")
    require(retained == active, "retained rows do not match certified deletions")
    a, h = matrix(domain.get("normals")), vector(domain.get("heights"))
    require(a == [normalized[i] for i in retained] and
            h == [normalized_heights[i] for i in retained], "retained coordinate mapping failed")

    # The original domain is bounded and all deletions preserve it. Verify that
    # every remaining inequality actually supports a facet: feasible vertices
    # on its hyperplane must affinely span dimension d-1. This differs from the
    # producer's exhaustive nonnegative dual-support redundancy check.
    facets = [[] for _ in retained]
    for support in combinations(range(len(a)), d):
        solution = affine_solution([a[i] for i in support], [h[i] for i in support], d)
        if solution is None or solution[1]:
            continue
        point = solution[0]
        slacks = [height - sum(x * y for x, y in zip(row, point)) for row, height in zip(a, h)]
        if min(slacks) < 0:
            continue
        for i, slack in enumerate(slacks):
            if slack == 0:
                facets[i].append(point)
    for i, points in enumerate(facets):
        require(points and rank([[x - y for x, y in zip(point, points[0])]
                                 for point in points[1:]], d) == d - 1,
                f"retained row {retained[i]} has no independently established facet")
    return retained, a, h


def verify_witness(witness, normals, heights, retained):
    require(isinstance(witness, dict), "missing ordered witness")
    order = indices(witness.get("order"), len(normals), "witness order")
    require(order and all(i in retained for i in order), "witness must use retained original rows")
    weights = vector(witness.get("weights"))
    require(len(weights) == len(order) and all(v >= 0 for v in weights), "invalid witness weights")
    require(sum(v * heights[i] for i, v in zip(order, weights)) == 1,
            "witness normalization failed")
    require(all(sum(v * normals[i][j] for i, v in zip(order, weights)) == 0
                for j in range(len(normals[0]))), "witness closure failed")
    objective = sum((weights[i] * weights[j] * omega(normals[order[i]], normals[order[j]])
                     for i in range(len(order)) for j in range(i)), F(0))
    require(objective > 0 and rational(witness.get("objective"), True) == objective,
            "witness objective failed")
    bound = F(1) / (2 * objective)
    require(rational(witness.get("capacity_upper_bound"), True) == bound,
            "witness reciprocal failed")
    return objective, bound


def verify(certificate, original=None, witness_only=False):
    require(isinstance(certificate, dict), "expected a certificate object")
    result = certificate.get("result", certificate)
    require(isinstance(result, dict) and result.get("schema") == "ehz-exact-result-v1",
            "unsupported exact result schema")
    normals, heights = exact_input(result.get("input"), True)
    if original is not None:
        require(exact_input(original) == (normals, heights), "certificate input identity differs")
    status = result.get("status")
    require(status in ("exact-complete", "exact-witness/incomplete", "incomplete", "certified-interval"),
            "result does not claim a verifiable exact outcome")
    completion = result.get("completion")
    require(isinstance(completion, dict) and type(completion.get("complete")) is bool,
            "missing completion declaration")
    require(completion["complete"] == (status == "exact-complete"), "contradictory completion status")
    require(isinstance(completion.get("reason"), str), "missing completion reason")
    interval = result.get("interval")
    require((interval is not None) == (status == "certified-interval"),
            "interval/status mismatch")
    witness = result.get("witness")
    require((witness is not None) == (status != "incomplete"), "witness/status mismatch")
    require(result.get("capacity") is None or status == "exact-complete",
            "incomplete search cannot claim an exact capacity")
    if result.get("domain") is None:
        require(status == "incomplete" and witness is None, "missing verified geometric domain")
        return {"verified": True, "status": "incomplete", "global_optimality_verified": False,
                "domain_verified": False, "reason": "No bound or exact capacity established."}
    retained, a, h = verify_domain(result["domain"], normals, heights)
    if witness is None:
        return {"verified": True, "status": "incomplete", "global_optimality_verified": False,
                "domain_verified": True, "reason": "No positive witness established."}
    objective, bound = verify_witness(witness, normals, heights, retained)
    report = {"verified": True, "status": "exact-witness/incomplete",
              "global_optimality_verified": False, "domain_verified": True,
              "capacity_upper_bound": canonical(bound)}
    if status == "certified-interval":
        require(isinstance(interval, dict) and interval.get("method") == "pairing-envelope-v1",
                "unsupported interval bound method")
        normalized = [[x / height for x in row] for row, height in zip(normals, heights)]
        envelope = max(abs(omega(a, b)) for a in normalized for b in normalized)
        require(envelope > 0, "nonpositive pairing envelope")
        upper_objective = envelope / 2
        lower_capacity = 1 / envelope
        require(rational(interval.get("objective_upper_bound"), True) == upper_objective and
                rational(interval.get("lower"), True) == lower_capacity and
                rational(interval.get("upper"), True) == bound,
                "incorrect certified interval bounds")
        require(objective <= upper_objective and lower_capacity <= bound,
                "invalid certified interval ordering")
        report.update(status="certified-interval",
                      capacity_interval=[canonical(lower_capacity), canonical(bound)])
    if status == "exact-complete":
        require(rational(result.get("capacity"), True) == bound, "exact capacity differs from witness")
        if not witness_only:
            maximum, orders = maximum_on_orders(a, h)
            require(maximum == objective, "independent enumeration contradicts exact-complete claim")
            report.update(status="exact-complete", global_optimality_verified=True,
                          capacity=canonical(bound), independent_orders=orders)
    return report


def load(path):
    def no_duplicates(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, "duplicate JSON object key")
            result[key] = value
        return result
    return json.loads(Path(path).read_text(), object_pairs_hook=no_duplicates)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", help="JSON result emitted by ehz-capacity --mode exact")
    parser.add_argument("--input", help="bind verification to this original exact input JSON")
    parser.add_argument("--witness-only", action="store_true",
                        help="check feasibility and the upper bound without certifying global optimality")
    args = parser.parse_args()
    try:
        report = verify(load(args.certificate), load(args.input) if args.input else None,
                        args.witness_only)
    except (ValueError, TypeError, KeyError, IndexError, OSError) as error:
        print(json.dumps({"verified": False, "error": str(error)}))
        return 1
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
