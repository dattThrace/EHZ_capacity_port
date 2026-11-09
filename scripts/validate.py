#!/usr/bin/env python3
"""Cross-language regression checks; not a MATLAB-runtime comparison."""
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
from scipy.spatial import ConvexHull
from reference import capacity, geometry

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "target/release/ehz-capacity"


def rust(data):
    output = subprocess.run([str(BINARY)], input=json.dumps(data), text=True, capture_output=True, check=True,
                            env=dict(os.environ, EHZ_PYTHON=sys.executable))
    result = json.loads(output.stdout)
    normals = np.array(data.get("normals", result.get("geometry", {}).get("normals")))
    heights = np.array(data.get("heights", result.get("geometry", {}).get("heights")))
    r = result["result"]
    order, beta = r["order"], np.array(r["beta"])
    assert np.linalg.norm(normals[order].T @ beta) < 1e-7
    assert abs(heights[order] @ beta - 1) < 1e-7
    assert beta.min() > -1e-10
    n = normals.shape[1] // 2
    active = normals[order]
    pairing = active[:, :n] @ active[:, n:].T - active[:, n:] @ active[:, :n].T
    score = beta @ np.tril(pairing, -1) @ beta
    assert abs(score - r["minimum"]) < 1e-7
    assert abs(-0.5 / score - r["capacity"]) < 1e-7
    return r["capacity"]


def close(a, b):
    assert np.isclose(a, b, atol=1e-7, rtol=1e-7), (a, b)


def main():
    count = 0
    for path in sorted((ROOT / "examples").glob("*.json")):
        data = json.loads(path.read_text())
        n, h = geometry(data["vertices"])
        expected = capacity(n, h)
        actual = rust(data)
        close(actual, expected)
        if path.stem in ("square", "readme_product"):
            close(actual, 4)
        if len(data["vertices"][0]) == 2:
            close(actual, ConvexHull(data["vertices"]).volume)
        print(f"PASS {path.stem}: Rust={actual:.12g}, Python={expected:.12g}")
        count += 1
    square = json.loads((ROOT / "examples/square.json").read_text())["vertices"]
    for label, vertices, expected in [
        ("translation", [[x+7, y-3] for x, y in square], 4),
        ("scaling", [[3*x, 3*y] for x, y in square], 36),
        ("shear", [[x+2*y, y] for x, y in square], 4),
        ("redundant vertices", square + [[0, 0], square[0]], 4),
    ]:
        close(rust({"vertices": vertices}), expected)
        print(f"PASS {label}")
        count += 1
    rng = np.random.default_rng(20261109)
    for i in range(8):
        vertices = rng.normal(size=(6, 2))
        close(rust({"vertices": vertices.tolist()}), ConvexHull(vertices).volume)
        count += 1
    n, h = geometry(square)
    scales = np.array([1, 2, 3, 4])
    close(rust({"normals": (n*scales[:, None]).tolist(), "heights": (h*scales).tolist()}), 4)
    count += 1
    for data in ({"vertices": [[0, 0], [1, 0], [2, 0]]}, {"normals": [[1, 0]], "heights": [0]}):
        p = subprocess.run([str(BINARY)], input=json.dumps(data), text=True, capture_output=True)
        assert p.returncode != 0 and "error:" in p.stderr
        count += 1
    print(f"PASS {count} regression cases, including witness residuals")


if __name__ == "__main__":
    main()
