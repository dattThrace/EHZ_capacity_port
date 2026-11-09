# Original MATLAB repository

- Author: Pazit Haim-Kislev
- Repository: <https://github.com/pazithaimkislev/EHZ-capacity>
- Branch: `main`
- Pinned commit: `22dd473cc8b9b2397b44cd13acbceb763dd978e9`
- Commit date: 2024-05-25
- Commit subject: `Update README.md`
- Commits in the preserved history: 5

`EHZ-capacity/` contains every tracked file extracted directly with `git archive`
at the pinned commit. File contents and original line endings are unchanged.
`EHZ-capacity.bundle` preserves the complete reachable history of the upstream
`main` branch, without copying local Git configuration or authentication data.

To obtain a standalone checkout with that history, run from this directory and
choose a destination that does not already exist:

```sh
git clone --branch main EHZ-capacity.bundle /path/to/matlab-reference
git -C /path/to/matlab-reference rev-parse HEAD
```

The printed revision should match the pinned commit above. A bundle clone's
`origin` refers to the bundle; to associate that separate checkout with the
upstream server afterward:

```sh
git -C /path/to/matlab-reference remote set-url origin https://github.com/pazithaimkislev/EHZ-capacity.git
```

The original repository has no license file at this revision. This snapshot
preserves its attribution and does not assign it a new license. Its mathematical
basis is Haim-Kislev, *On the Symplectic Size of Convex Polytopes*, GAFA 29 (2019),
440–463, <https://doi.org/10.1007/s00039-019-00486-4>.
