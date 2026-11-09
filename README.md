# EHZ capacity port

Research and reproducible validation for a planned Rust implementation of the
Ekeland–Hofer–Zehnder capacity of convex polytopes.

This initial repository preserves the research, Python validation probes, and
the original MATLAB implementation. The Rust library and CLI have not been
implemented yet. The probes are small mathematical checks, not a general-purpose
capacity solver or a MATLAB compatibility test suite.

## Contents

| Location | Contents |
| --- | --- |
| [Research memo](research/RESEARCH.md) | Literature review, mathematical conditions, proposed algorithms, and implementation plan |
| [Python probes](research/probes/) | Exact ordering DP, exact symmetric-product witnesses, and numerical itinerary LP checks |
| [Recorded results](research/results/) | Original run output, dependency versions, seed, and script hashes |
| [Sources](research/sources/) | Downloaded papers, extracted text, and source URL/version/SHA-256 manifest |
| [Discovery records](research/discovery/) | Original search responses, arXiv queries, and repository metadata |
| [Original MATLAB source](upstream/EHZ-capacity/) | Unmodified tracked files at upstream commit `22dd473cc8b9b2397b44cd13acbceb763dd978e9` |
| [Original Git history](upstream/EHZ-capacity.bundle) | Restorable bundle containing the upstream `main` branch and its five commits |
| [Upstream provenance](upstream/README.md) | Source URL, pinned revision, attribution, and history restoration commands |

The research artifacts are preserved byte for byte from the original workspace.
The memo's status statements and absolute `/workspace/scratch/ehz-research` paths
describe that earlier research session. Use the repository-relative commands
below to run the probes in a new checkout. Generated Python bytecode was omitted.

## Run the validation probes

The recorded environment used Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, and
SciPy's bundled HiGHS 1.8.0. The ordering and symmetric-product probes use only the
Python standard library; the itinerary probe also needs NumPy and SciPy.

From the repository root, using Python 3.12:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r research/requirements.txt
python research/probes/check_ordering_dp.py
python research/probes/check_symmetric_product.py
python research/probes/check_itinerary_lp.py
```

On Windows, activate the environment with `.venv\Scripts\Activate.ps1` in
PowerShell. Run without `-O` or `PYTHONOPTIMIZE`: these probes use assertions.

Expected coverage:

- 22 exact comparisons of ordering DP with exhaustive permutations/cyclic orders.
- 23 exact symmetric-product formula and orbit-witness checks.
- Numerical itinerary LP checks for three simplices, the square, and a selected
  README-product orbit, plus rejection of its reversed itinerary.

The LP results use floating-point tolerances. The product itinerary check is not
an exhaustive search of the eight-facet product. The exact symmetric-product
formula relies on the published inradius theorem. See the memo for the scope and
mathematical justification of each probe.

## Download and preserve

GitHub's **Code → Download ZIP** downloads this whole snapshot, including the
papers and upstream bundle. Alternatively:

```sh
git clone --branch work https://github.com/dattThrace/EHZ_capacity_port.git
```

The files are ordinary Git contents, so obtaining them does not require access
to the original Codex workspace. There are no Git submodules or Git LFS objects.

## Attribution

The original MATLAB implementation is by Pazit Haim-Kislev and is preserved with
its original README and history. Papers and other source materials retain their
authors' notices and terms; citations and source URLs are in the research memo
and source manifest. No blanket license is assigned to these collected materials.
