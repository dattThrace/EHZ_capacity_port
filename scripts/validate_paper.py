#!/usr/bin/env python3
"""Validate Rust against the PDF mathematics, never scripts/reference.py."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from fractions import Fraction as F
from paper_oracle import omega, symmetric_capacity, theorem_capacity

ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "research/sources/1712.03494.pdf"
PDF_SHA256 = "37a977d8039d6fc7e772dba4dfb046f102d1693846a9f7f3bc8b318d68e5b5d5"
BINARY = ROOT / "target/release/ehz-capacity"


def call(data, check=True):
    p = subprocess.run([str(BINARY)], input=json.dumps(data), text=True, capture_output=True,
                       env=dict(os.environ, EHZ_PYTHON=sys.executable))
    if not check:
        assert p.returncode != 0, p.stdout
        return None
    assert p.returncode == 0, p.stderr
    return json.loads(p.stdout)


def close(actual, expected):
    assert abs(actual/float(expected)-1) < 2e-8, (actual, expected)


def witness(data, result):
    normals = data["normals"]
    heights = data["heights"]
    for prefix, expected_sign in (("", -1), ("maximizing_", 1)):
        order = result[prefix + "order"]
        beta = result[prefix + "beta"]
        assert len(order) == len(set(order)) == len(beta)
        assert min(beta) >= 0
        t = [F.from_float(v) * F(heights[i]) for i, v in zip(order, beta)]
        assert abs(float(sum(t))-1) < 1e-8
        for j in range(len(normals[0])):
            terms = [F.from_float(b)*F(normals[i][j]) for i, b in zip(order, beta)]
            assert abs(float(sum(terms))) <= 1e-8*max(1, float(sum(abs(v) for v in terms)))
        # Evaluate the printed witness directly in the paper's coordinates.
        q = sum((F.from_float(beta[i])*F.from_float(beta[j]) *
                 omega(normals[order[i]], normals[order[j]])
                 for i in range(len(order)) for j in range(i)), F(0))
        close(float(q), expected_sign*result["maximum"])
        close(result["capacity"], 1/(2*abs(q)))


def polygon_facets(vertices):
    """Exact CCW edges, independent of the production Qhull/cofactor adapter."""
    n, h = [], []
    for a, b in zip(vertices, vertices[1:]+vertices[:1]):
        normal = [b[1]-a[1], a[0]-b[0]]
        n.append(normal)
        h.append(sum(x*y for x, y in zip(normal, a)))
    return {"normals": n, "heights": h}


def simplex(d):
    # Translate the standard simplex so that all heights are positive integers:
    # x_i >= -1, sum x_i <= 1. No hull implementation involved.
    return {"normals": [[-int(i == j) for j in range(d)] for i in range(d)] + [[1]*d],
            "heights": [1]*(d+1)}


def main():
    assert hashlib.sha256(PDF.read_bytes()).hexdigest() == PDF_SHA256
    records = []

    def check(name, data, expected, source, exact_orders=0):
        output = call(data)
        result = output["result"]
        close(result["capacity"], expected)
        facets = data if "normals" in data else output["geometry"]
        witness(facets, result)
        records.append({"case": name, "expected": str(expected), "actual": result["capacity"],
                        "basis": source, "oracle_orders": exact_orders})
        print(f"PASS {name}: {result['capacity']:.12g}, exact={expected}", flush=True)

    # RREF/nullspace exact maximization of Formula (1), over all permutations.
    polygons = [
        ("square", [[-1,-1],[1,-1],[1,1],[-1,1]]),
        ("triangle", [[-2,-1],[3,-1],[-1,3]]),
        ("asymmetric quadrilateral", [[-2,-1],[2,-2],[3,1],[-1,3]]),
        ("pentagon", [[-2,-1],[1,-2],[3,0],[1,3],[-2,2]]),
    ]
    for name, vertices in polygons:
        data = polygon_facets(vertices)
        exact, orders = theorem_capacity(**data)
        area = abs(sum(F(a[0]*b[1]-a[1]*b[0], 2) for a, b in zip(vertices, vertices[1:]+vertices[:1])))
        assert exact == area
        check(name + " facets", data, exact, "Theorem 1.1 + planar action/area", orders)
        check(name + " vertices", {"vertices": vertices}, exact, "same independent polygon/area")
    for d in (4, 6):
        data = simplex(d)
        exact, orders = theorem_capacity(**data)
        assert exact == F((d+1)**2, d)
        check(f"simplex R{d}", data, exact, "Theorem 1.1 exact nullspace oracle", orders)
    # A nonsymmetric six-facet body has free weights. This checks the stationary
    # branch independently of the symmetric corollary and unique-weight simplices.
    for height in (1, 2):
        data = simplex(4)
        data["normals"].append([1,1,0,0])
        data["heights"].append(height)
        exact, orders = theorem_capacity(**data)
        check(f"truncated R4 simplex height {height}", data, exact,
              "Theorem 1.1 exact free-coordinate quadratic maximization", orders)
    # Symmetric bodies: use the distinct formula from Corollary 1.2 as oracle.
    symmetric = [
        ("four-dimensional box", [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], [1,2,3,4]),
        ("upstream README product", [[1,0,0,0],[0,1,0,0],[0,0,1,1],[0,0,1,-1]], [1]*4),
        ("skew symmetric polytope", [[1,2,0,1],[0,1,1,0],[1,0,2,1],[0,0,1,1]], [2,3,4,5]),
    ]
    for name, reps, heights in symmetric:
        exact, orders = symmetric_capacity(reps, heights)
        data = {"normals": reps + [[-v for v in r] for r in reps], "heights": heights*2}
        check(name, data, exact, "Corollary 1.2 exact signed-weight oracle", orders)
    square = polygon_facets(polygons[0][1])
    for scale in (1e-9, 1e-6, 1e6, 1e9):
        check(f"height scaling {scale}", {"normals": square["normals"], "heights": [v*scale for v in square["heights"]]},
              F(4)*F(str(scale))**2, "Remark 3.9 quadratic homogeneity")
        check(f"vertex scaling {scale}", {"vertices": [[v*scale for v in r] for r in polygons[0][1]]},
              F(4)*F(str(scale))**2, "Remark 3.9 quadratic homogeneity")
    for scales in ([1e-9]*4, [1e-4]*4, [1e9]*4, [1e-8,1e8,1e-4,1e4]):
        check(f"facet row scales {scales}", {"normals": [[v*s for v in r] for r, s in zip(square["normals"], scales)],
                                               "heights": [v*s for v, s in zip(square["heights"], scales)]},
              F(4), "Theorem 1.1 beta/row substitution")
    # q'=q+Bp with symmetric B, on a rational simplex in R4.
    data = simplex(4)
    transformed = {"normals": [[a,b,c-2*a-b,d-a-3*b] for a,b,c,d in data["normals"]], "heights": data["heights"]}
    exact, orders = theorem_capacity(**transformed)
    assert exact == F(25,4)
    check("R4 symplectic shear", transformed, exact, "Remark 3.10 + exact Theorem 1.1", orders)
    # Invalid theorem domain: unbounded, including full-rank normals with a finite
    # closed characteristic in only one symplectic plane.
    for name, data in [
        ("unbounded full rank", {"normals":[[1,0,0,0],[-1,0,0,0],[0,0,1,0],[0,0,-1,0],[0,1,0,0],[0,0,0,1]], "heights":[1]*6}),
        ("unbounded rank deficient", {"normals":[[1,0,0,0],[-1,0,0,0],[0,0,1,0],[0,0,-1,0],[1,0,1,0]], "heights":[1]*5}),
        ("zero normal", {"normals":[[1,0],[-1,0],[0,0]], "heights":[1]*3}),
    ]:
        call(data, check=False)
        records.append({"case": name, "expected": "reject: outside theorem domain"})
        print(f"PASS reject {name}", flush=True)
    report = {"source": "arXiv:1712.03494v3", "pdf_sha256": PDF_SHA256, "cases": records}
    destination = ROOT / "validation/paper-results.json"
    destination.parent.mkdir(exist_ok=True)
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS {len(records)} source-paper checks")


if __name__ == "__main__":
    main()
