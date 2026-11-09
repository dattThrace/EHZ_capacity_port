# Research archive

Start with [RESEARCH.md](RESEARCH.md) for the mathematical and implementation
recommendations. All 49 non-cache files from the original research directory are
preserved verbatim, with this layout:

- `probes/`: the three executable Python validation probes.
- `results/`: the three original stdout logs and `probe-environment.json`.
- `sources/`: eleven PDF papers, their extracted text, the arXiv HTML article,
  and `sources.json`. Manifest `file` paths are relative to this directory.
- `discovery/`: original search/API responses and repository metadata. Some
  responses were irrelevant or bot-challenge pages; they are retained as a record
  of retrieval attempts, not relied on as mathematical evidence.

In `results/probe-environment.json`, script filenames refer to files in
`probes/`. The recorded SHA-256 hashes still apply because the probes are
unchanged. The file's platform data describes the original run, not a requirement
for executing the probes elsewhere.

For current run instructions, see the [root README](../README.md). The old
workspace paths in the memo and script docstrings document their original
locations. `requirements.txt` pins the two numerical packages used by the LP
probe; the other two probes require only the Python standard library.

No additional code or conclusions were introduced by archiving these files.
MATLAB execution and production-port validation remain future work.
