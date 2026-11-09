#!/usr/bin/env python3
"""Qhull simplex-index adapter. Numerical geometry remains in the Rust library."""
import json
import sys
import numpy as np
from scipy.spatial import ConvexHull, Delaunay


def triangulate(vertices):
    points = np.asarray(vertices, dtype=float)
    # Match MATLAB's explicit QJ option (Delaunay additionally requires Qt).
    simplices = (
        [list(range(len(points)))]
        if len(points) == points.shape[1] + 1
        else Delaunay(points, qhull_options="QJ Qt").simplices.tolist()
    )
    return {
        "delaunay": simplices,
        "hull": ConvexHull(points, qhull_options="QJ Qt").simplices.tolist(),
    }


if __name__ == "__main__":
    try:
        json.dump(triangulate(json.load(sys.stdin)), sys.stdout, allow_nan=False)
    except Exception as error:
        print(str(error), file=sys.stderr)
        sys.exit(1)
