#!/usr/bin/env python3
"""Build/run the Rust CLI, with the current Python interpreter for Qhull."""
import os
from pathlib import Path
import subprocess
import sys

if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    env = dict(os.environ, EHZ_PYTHON=sys.executable)
    raise SystemExit(subprocess.call(["cargo", "run", "--release", "--locked", "--manifest-path", str(root / "Cargo.toml"), "--", *sys.argv[1:]], env=env))
