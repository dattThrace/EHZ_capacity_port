#!/usr/bin/env python3
"""Independent small floating reference for EHZ_perms.m. No search reductions."""
import itertools
import json
import math
import sys
import numpy as np
from scipy.linalg import qr
from scipy.spatial import ConvexHull, Delaunay


def ktimes_t(k, t):
    return [list(q) + list(p) for p in t for q in k]


def barycenter(k):
    k = np.asarray(k, dtype=float)
    simplices = np.array([range(len(k))]) if len(k) == k.shape[1] + 1 else Delaunay(k, qhull_options="QJ Qt").simplices
    volumes = np.array([abs(np.linalg.det(k[s[1:]] - k[s[0]])) for s in simplices])
    return np.sum(volumes[:, None] * k[simplices].sum(axis=1), axis=0) / volumes.sum() / (k.shape[1] + 1)


def geometry(k):
    k = np.asarray(k, dtype=float)
    center = barycenter(k)
    k = k - center
    hull = ConvexHull(k, qhull_options="QJ Qt")
    used = sorted(set(hull.simplices.flat))
    mapping = {v: i for i, v in enumerate(used)}
    k = k[used]
    normals, faces = [], []
    d = k.shape[1]
    for f in sorted(sorted(mapping[int(i)] for i in f) for f in hull.simplices):
        a = k[f[1:]] - k[f[0]]
        if np.linalg.matrix_rank(a, tol=1e-10) < d - 1:
            continue
        normal = np.array([(-1)**(j+1) * np.linalg.det(np.delete(a, j, axis=1)) for j in range(d)])
        normal /= np.linalg.norm(normal)
        for v in k:
            side = (v-k[f[0]]) @ normal
            if abs(side) > 1e-10:
                if side > 0:
                    normal = -normal
                break
        if not any(np.all(np.abs(n-normal) <= 1e-8) for n in normals):
            normals.append(normal)
            faces.append(f)
    normals = np.array(normals)
    heights = np.array([k[f[0]] @ n for f, n in zip(faces, normals)])
    return normals, heights


def prepare(a, b):
    a = a.copy()
    a[abs(a) < 1e-10] = 0
    _, r, p = qr(a.T, pivoting=True)
    _, r2, _ = qr(np.column_stack([a, b]).T, pivoting=True)
    def rank(r):
        indices = np.flatnonzero(abs(np.diag(r)) > 1e-8)
        return int(indices[-1]+1) if len(indices) else 0
    rank_a = rank(r)
    if rank_a != rank(r2):
        return None
    return a[p[:rank_a]], b[p[:rank_a]]


def capacity(normals, heights):
    normals, heights = np.array(normals, dtype=float), np.array(heights, dtype=float)
    m, d = normals.shape
    n = d // 2
    omega = normals[:, :n] @ normals[:, n:].T - normals[:, n:] @ normals[:, :n].T
    minimum = math.inf
    for size in range(2, m+1):
        for support in itertools.combinations(range(m), size):
            selected = list(support)
            system = prepare(np.vstack([normals[selected].T, heights[selected]]), np.r_[np.zeros(d), 1.])
            if system is None:
                continue
            a, b = system
            unique = np.linalg.solve(a, b) if a.shape[0] == size else None
            if unique is not None and not np.all(unique > -1e-10):
                continue
            for tail in itertools.permutations(range(size-1, 0, -1)):
                route = [0, *tail]
                order = [support[i] for i in route]
                h = np.tril(omega[np.ix_(order, order)], -1)
                h = h + h.T
                if unique is None:
                    ar = a[:, route]
                    mat = np.block([[h, ar.T], [ar, np.zeros((len(b), len(b)))]])
                    if abs(np.linalg.det(mat)) < 1e-10:
                        continue
                    beta = np.linalg.solve(mat, np.r_[np.zeros(size), b])[:size]
                else:
                    beta = unique[route]
                if np.all(beta > -1e-10):
                    minimum = min(minimum, beta @ h @ beta / 2)
    if not math.isfinite(minimum) or minimum >= 0:
        raise ValueError("no negative feasible score")
    return -0.5 / minimum


if __name__ == "__main__":
    with open(sys.argv[1]) as source:
        data = json.load(source)
    normals, heights = geometry(data["vertices"]) if "vertices" in data else (data["normals"], data["heights"])
    print(json.dumps({"capacity": capacity(normals, heights)}))
