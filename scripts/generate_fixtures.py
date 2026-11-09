#!/usr/bin/env python3
"""Reproduce small inputs, including the original README example."""
import json
from pathlib import Path
from reference import ktimes_t


def generate(destination):
    destination.mkdir(parents=True, exist_ok=True)
    square = [[-1, -1], [-1, 1], [1, -1], [1, 1]]
    fixtures = {
        "square": square,
        "triangle": [[0, 0], [2, 0], [0, 3]],
        "simplex4": [[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
        "readme_product": ktimes_t(square, [[1, 0], [0, 1], [-1, 0], [0, -1]]),
    }
    for name, vertices in fixtures.items():
        (destination / f"{name}.json").write_text(json.dumps({"vertices": vertices}, indent=2) + "\n")


if __name__ == "__main__":
    generate(Path(__file__).resolve().parents[1] / "examples")
