#!/usr/bin/env python3
"""Exact acceptance gates and reproducible, unpruned rational-solver benchmarks.

The expected values stay fractions.Fraction throughout. The paper oracle uses
free-coordinate stationary equations and full permutations, independently of
Rust's KKT elimination and cyclic anchor. Requires Python's standard library.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
from random import Random
import subprocess
import tempfile
import time

from paper_oracle import omega, symmetric_capacity, theorem_capacity
from verify_exact import verify

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "ehz-exact-h-v1"


def rational(value):
    if isinstance(value, list):
        return F(int(value[0]), int(value[1]))
    return F(value)


def encode(normals, heights):
    return {"schema": SCHEMA, "normals": [[str(F(x)) for x in row] for row in normals],
            "heights": [str(F(h)) for h in heights]}


def decode(data):
    return [[rational(x) for x in row] for row in data["normals"]], list(map(rational, data["heights"]))


def polygon(vertices):
    """CCW outward facets and shoelace area, independently of Qhull."""
    vertices = [[F(x) for x in row] for row in vertices]
    normals, heights, area = [], [], F(0)
    for a, b in zip(vertices, vertices[1:] + vertices[:1]):
        normal = [b[1] - a[1], a[0] - b[0]]
        normals.append(normal)
        heights.append(sum(x*y for x, y in zip(normal, a)))
        area += (a[0]*b[1] - a[1]*b[0])/2
    assert area > 0 and min(heights) > 0
    # All edges must support the convex polygon; no accidental nonconvex family.
    assert all(sum(x*y for x, y in zip(n, v)) <= h
               for n, h in zip(normals, heights) for v in vertices)
    return encode(normals, heights), area


def cases():
    square = encode([[1, 0], [-1, 0], [0, 1], [0, -1]], [1]*4)
    yield "square", square, F(4), "planar area; full four-row support; tied cyclic optimizers"
    decimal = {"schema": SCHEMA, "normals": [["1.0", "0"], ["-1.00", "0"],
                 ["0", ["2", "2"]], ["0", "-1"]], "heights": ["0.1"]*4}
    yield "decimal_square", decimal, F(1, 25), "exact decimal parsing and rational pairs"
    data, area = polygon([[-2, -1], [2, -2], [3, 1], [-1, 3]])
    yield "asymmetric_quadrilateral", data, area, "independent shoelace area and free weights"
    rng = Random(20261109)
    for i in range(4):
        points = [[F(x) + F(rng.randint(-2, 2), 16) for x in row]
                  for row in [[-2, -1], [1, -2], [3, 0], [1, 3], [-2, 2]]]
        data, area = polygon(points)
        yield f"seeded_rational_pentagon_{i}", data, area, "seeded nonsymmetric free weights; exact planar area"
    normals = [[-int(i == j) for j in range(4)] for i in range(4)] + [[1]*4]
    yield "simplex_r4", encode(normals, [1]*5), F(25, 4), "Theorem 1.1 exact simplex"
    for h in [F(1), F(3, 2)]:
        data = encode(normals + [[1, 1, 0, 0]], [1]*5 + [h])
        expected = F(21, 4) if h == 1 else F(91, 16)
        yield f"truncated_simplex_r4_{str(h).replace('/', '_')}", data, expected, "nonsymmetric R4 free-weight stationary branch"
    # Rational symplectic q'=q+Bp with symmetric B, including large denominators.
    a, b, c = F(2, 7), F(1, 11), F(3, 13)
    transformed = [[q1, q2, p1-a*q1-b*q2, p2-b*q1-c*q2] for q1, q2, p1, p2 in normals]
    yield "symplectic_shear_r4", encode(transformed, [1]*5), F(25, 4), "exact rational symplectic invariance"
    duplicate_n = [[1, 0], [-1, 0], [0, 1], [0, -1], [2, 0], [0, -3]]
    yield "duplicate_square", encode(duplicate_n, [1, 1, 1, 1, 2, 3]), F(4), "positive duplicate-row plateaus"
    redundant_n = duplicate_n + [[1, 1], [1, 0], [-1, -1]]
    yield "redundant_square", encode(redundant_n, [1, 1, 1, 1, 2, 3, 3, 2, 2]), F(4), "certified redundant inequalities, including a vertex-tight row"
    scales = [F(2**160+1, 3**30), F(3**80, 2**70), F(7, 2**180), F(2**120, 11)]
    n, h = decode(square)
    yield "large_row_scales", encode([[v*s for v in row] for row, s in zip(n, scales)],
                                      [v*s for v, s in zip(h, scales)]), F(4), "equivalent positive row scales; independent input bit length"
    for k in [0, 20, 60, 120]:
        for sign in [-1, 0, 1]:
            eps = sign*F(1, 2**(k+1))
            # Width q=2+eps and p=2 stay positive for the entire family.
            data = encode(n, [1+eps, 1, 1, 1])
            yield f"perturb_square_k{k}_{sign:+d}", data, F(4)+2*eps, "rational 0,+/-2^-k family; distinctions beyond binary64"
    # In R2, this is the symmetric product [-2,2] x [-3,3].
    reps, hs = [[1, 0], [0, 1]], [F(2), F(3)]
    product_exact, _ = symmetric_capacity(reps, hs)
    assert product_exact == 24
    yield "symmetric_product", encode(reps + [[-x for x in row] for row in reps], hs*2), F(24), "independent Corollary 1.2 and interval-product formula"


def call(binary, data, flags=()):
    """Measure one process, including its per-process peak RSS via wait4."""
    work = ROOT / "work"
    work.mkdir(exist_ok=True)
    with tempfile.TemporaryFile(dir=work) as stdout, tempfile.TemporaryFile(dir=work) as stderr:
        start = time.perf_counter()
        child = subprocess.Popen([str(binary), *flags], stdin=subprocess.PIPE, stdout=stdout, stderr=stderr)
        child.stdin.write(json.dumps(data).encode())
        child.stdin.close()
        if hasattr(os, "wait4"):
            _, status, usage = os.wait4(child.pid, 0)
            child.returncode = os.waitstatus_to_exitcode(status)
            rss = usage.ru_maxrss  # KiB on Linux, the benchmark platform.
            cpu = usage.ru_utime + usage.ru_stime
        else:
            child.wait()
            rss, cpu = None, None
        elapsed = time.perf_counter() - start
        stdout.seek(0)
        stderr.seek(0)
        out, err = stdout.read().decode(), stderr.read().decode()
    try:
        result = json.loads(out)["result"]
    except (ValueError, KeyError) as exc:
        raise AssertionError(f"solver did not produce a result: {out!r}; stderr={err!r}") from exc
    return child.returncode, result, {"wall_seconds": elapsed, "cpu_seconds": cpu, "peak_rss_kib": rss}


def witness(data, result):
    n, h = decode(data)
    w = result["witness"]
    order, beta = w["order"], list(map(rational, w["weights"]))
    assert len(order) == len(set(order)) == len(beta) and min(beta) >= 0
    assert all(0 <= i < len(n) for i in order)
    assert sum(h[i]*b for i, b in zip(order, beta)) == 1
    assert all(sum(n[i][j]*b for i, b in zip(order, beta)) == 0 for j in range(len(n[0])))
    def score(route, weights):
        return sum((weights[i]*weights[j]*omega(n[route[i]], n[route[j]])
                    for i in range(len(route)) for j in range(i)), F(0))
    objective = score(order, beta)
    assert objective > 0 and objective == rational(w["objective"])
    assert rational(w["capacity_upper_bound"]) == 1/(2*objective)
    # Closure justifies cyclic anchoring; reversal changes the sign exactly.
    assert score(order[1:]+order[:1], beta[1:]+beta[:1]) == objective
    assert score(list(reversed(order)), list(reversed(beta))) == -objective


def input_bits(data):
    n, h = decode(data)
    return max(max(abs(x.numerator).bit_length(), x.denominator.bit_length()) for row in n+[h] for x in row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "target/release/ehz-capacity")
    parser.add_argument("--benchmark", action="store_true", help="write validation benchmark JSON and Markdown")
    args = parser.parse_args()
    binary = args.binary.resolve()
    records, lifecycle = [], []
    fixture_dir = ROOT / "examples/exact"
    fixture_dir.mkdir(parents=True, exist_ok=True)
    for name, data, expected, source in cases():
        # Keep a compact fixed fixture subset; parameter families remain reproducible here.
        if name in {"square", "decimal_square", "asymmetric_quadrilateral", "seeded_rational_pentagon_0",
                    "simplex_r4", "truncated_simplex_r4_3_2", "symplectic_shear_r4",
                    "duplicate_square", "redundant_square", "large_row_scales", "perturb_square_k120_+1", "symmetric_product"}:
            (fixture_dir / f"{name}.json").write_text(json.dumps(data, indent=2)+"\n")
        normals, heights = decode(data)
        # Duplicates/redundancies have already been checked against the same square.
        # Full original-row permutation oracles would add cost without new coverage.
        oracle_orders = None
        if len(normals) <= 6:
            oracle, oracle_orders = theorem_capacity(normals, heights)
            assert oracle == expected, (name, oracle, expected)
        code, result, measures = call(binary, data)
        assert code == 0, (name, result)
        assert result["status"] == "exact-complete" and result["completion"]["complete"], (name, result)
        assert rational(result["capacity"]) == expected, (name, result["capacity"], expected)
        witness(data, result)
        verified = verify({"result": result}, original=data)
        assert verified["global_optimality_verified"], (name, verified)
        if name == "square":
            assert len(result["witness"]["order"]) == 4
        record = {"case": name, "dimension": len(normals[0]), "input_facets": len(normals),
                  "retained_facets": len(result["domain"]["retained"]), "input_peak_bits": input_bits(data),
                  "expected": str(expected), "actual": result["capacity"], "basis": source,
                  "oracle_full_permutations": oracle_orders, "status": result["status"],
                  "independent_verifier_orders": verified["independent_orders"],
                  "stats": result["stats"], **measures}
        records.append(record)
        print(f"PASS {name}: {expected}", flush=True)
    square = next(data for name, data, *_ in cases() if name == "square")
    for flags in [("--max-steps", "0"), ("--timeout-ms", "0")]:
        code, result, _ = call(binary, square, flags)
        assert result["status"] == "incomplete" and not result["completion"]["complete"]
        assert result.get("capacity") is None
        assert not verify({"result": result}, original=square)["global_optimality_verified"]
        lifecycle.append({"flags": list(flags), "status": result["status"], "stats": result["stats"]})
        print(f"PASS interruption {flags[0]}", flush=True)
    # Find a budget preserving a verified incumbent before enumeration completes.
    complete_steps = records[0]["stats"].get("steps")
    assert isinstance(complete_steps, int) and complete_steps > 0
    found = False
    if complete_steps:
        for limit in range(1, complete_steps):
            _, result, _ = call(binary, square, ("--max-steps", str(limit)))
            if result.get("witness"):
                assert result["status"] == "exact-witness/incomplete" and not result["completion"]["complete"]
                assert result.get("capacity") is None
                witness(square, result)
                assert not verify({"result": result}, original=square)["global_optimality_verified"]
                lifecycle.append({"flags": ["--max-steps", str(limit)], "status": result["status"],
                                  "stats": result["stats"], "capacity_upper_bound": result["witness"]["capacity_upper_bound"]})
                found = True
                break
        assert found, "square should expose a verified partial incumbent with an intermediate budget"
        print("PASS interruption preserves exact witness", flush=True)
        _, result, _ = call(binary, square, ("--max-steps", str(limit), "--certify-interval"))
        assert result["status"] == "certified-interval" and not result["completion"]["complete"]
        assert result.get("capacity") is None
        report = verify({"result": result}, original=square)
        lower, upper = map(F, report["capacity_interval"])
        assert lower <= F(4) <= upper and not report["global_optimality_verified"]
        lifecycle.append({"flags": ["--max-steps", str(limit), "--certify-interval"],
                          "status": result["status"], "stats": result["stats"], "interval": result["interval"]})
        print("PASS independently certified partial interval", flush=True)
    invalid_cases = [
        ("zero denominator", {**square, "heights": [["1", "0"]]*4}),
        ("malformed fraction", {**square, "heights": ["1//2"]*4}),
        ("machine JSON numbers", {**square, "heights": [1]*4}),
        ("zero height", {**square, "heights": ["0"]*4}),
        ("zero normal", encode([[1, 0], [-1, 0], [0, 0]], [1]*3)),
        ("unbounded full rank", encode([[1, 0], [-1, 0], [0, 1]], [1]*3)),
        ("unbounded rank deficient", encode([[1, 0], [-1, 0], [2, 0]], [1]*3)),
        ("odd dimension", encode([[1], [-1]], [1]*2)),
    ]
    for name, data in invalid_cases:
        code, result, _ = call(binary, data)
        assert code != 0 and result["status"] == "invalid-input", (name, code, result)
        print(f"PASS reject {name}", flush=True)
    if args.benchmark:
        report = {"schema": "ehz-rational-benchmark-v1", "binary": str(binary.relative_to(ROOT)),
                  "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
                  "measured_at_utc": datetime.now(timezone.utc).isoformat(),
                  "toolchain_spec": (ROOT / "rust-toolchain.toml").read_text(),
                  "python": os.sys.version.split()[0], "platform": os.uname().sysname + " " + os.uname().machine,
                  "memory_method": "wait4 ru_maxrss in KiB (Linux), including process-launch footprint; solver VmHWM separately in stats.peak_rss_kib",
                  "exact_checks": len(records), "records": records, "lifecycle": lifecycle}
        (ROOT / "validation/rational-results.json").write_text(json.dumps(report, indent=2)+"\n")
        lines = ["# Rational solver baseline benchmark", "",
                 f"Measured at `{report['measured_at_utc']}` on `{report['platform']}` with Python `{report['python']}`.",
                 "Binary path, SHA256, and the pinned toolchain specification are recorded in",
                 "[`rational-results.json`](rational-results.json).", "",
                 "Reproduce with `cargo build --release --locked` followed by",
                 "`python3 scripts/validate_exact.py --benchmark`. Every expected value is compared",
                 "as an exact fraction. The independent paper oracle uses nullspace coordinates and",
                 "full permutations; Rust uses KKT systems and anchored cyclic orders. Planar",
                 "shoelace areas and the symmetric corollary supply additional independent checks.", "",
                 "These are correctness baseline measurements, with no pruning or speedup target.",
                 "Runtime and memory include one fresh Rust CLI process per case. Launch RSS comes",
                 "from Linux `wait4`, independently of the solver, and includes the pre-exec",
                 "process-launch footprint. Solver RSS is its `/proc/self/status` VmHWM near",
                 "search completion; these two measurements have different scopes.",
                 "Intermediate bit lengths are",
                 "the solver's observed elimination/search arithmetic statistic; see the machine",
                 "report for complete counters, CPU times, fractions, and input bit lengths.",
                 "The arithmetic statistic excludes transient BigInt implementation internals",
                 "during rational normalization. All emitted certificates are also independently",
                 "verified, including input identity, row transformations, domain, and global",
                 "coverage re-enumeration. The binary uses the pinned Rust 1.85.1 release profile.", "",
                 "| Case | d | Input/retained rows | Input bits | Intermediate bits | Supports | Orders | Time (s) | Launch/solver RSS (KiB) | Status |",
                 "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |"]
        for row in records:
            s = row["stats"]
            lines.append(f"| {row['case']} | {row['dimension']} | {row['input_facets']}/{row['retained_facets']} | {row['input_peak_bits']} | {s.get('peak_rational_bits', 'see JSON')} | {s.get('supports', 'see JSON')} | {s.get('orders', 'see JSON')} | {row['wall_seconds']:.6f} | {row['peak_rss_kib']}/{s.get('peak_rss_kib')} | {row['status']} |")
        lines += ["", "The perturbation family holds facet count and dimension fixed while increasing",
                  "denominator bit length. All heights remain strictly positive, including the",
                  "negative perturbation cases. Values at 60 and 120 bits differ even when binary64",
                  "would identify the inputs. Large row scales preserve exactly the same square.", "",
                  "Budget exhaustion, zero-timeout behavior, preservation of a partial exact witness,",
                  "malformed fractions, zero denominators, JSON floats/numbers, and invalid geometric",
                  "domains are checked by the same command. Certificate mutation and independent",
                  "coverage re-enumeration are checked by `python3 scripts/test_verify_exact.py`.", ""]
        lines += ["Lifecycle checks also record incomplete domain checks, interrupted feasible",
                  "witnesses, and conservative certified intervals in the machine report. Partial",
                  "witnesses remain capacity upper bounds and never produce exact-complete status.", ""]
        (ROOT / "validation/RATIONAL_BENCHMARK.md").write_text("\n".join(lines))
    print(f"PASS {len(records)} exact capacities, lifecycle checks, and {len(invalid_cases)} invalid inputs", flush=True)


if __name__ == "__main__":
    main()
