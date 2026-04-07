# Changelog

All notable changes to this project will be documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project adheres to [Semantic Versioning](https://semver.org/).

## [0.1.1] — 2026-04-07

This is an honest hotfix release. 0.1.0 shipped with three latent
data-corruption bugs and a README quickstart that no pip user could run.
0.1.1 fixes the bugs, adds regression tests, and ships the example
fixtures inside the wheel so the quickstart actually works.

### Added
- `chilmesh.examples` module with `annulus()`, `donut()`, `block_o()`,
  `structured()` factories that return ready-to-use `CHILmesh` instances
  loaded from bundled `.fort.14` fixtures.
- `chilmesh/data/` package shipping all four fixtures via
  `[tool.setuptools.package-data]`, accessible from any install.
- `chilmesh/utils/__init__.py` so `chilmesh.utils` is a real package
  and `plot_utils.py` always lands in the wheel (B6).
- Regression tests for every B1-B6 bug, parametrized over the four
  fixtures where applicable; layer-disjoint-cover invariants replacing
  the old value-pinned annulus assertions.
- `pyproject.toml` keywords, classifiers, `[project.urls]`, and a
  `[project.optional-dependencies] dev` group.
- GitHub Actions CI workflow (`.github/workflows/ci.yml`) running on
  Python 3.10/3.11/3.12 × ubuntu/macos.
- This `CHANGELOG.md`.

### Fixed
- **B1**: `CHILmesh.write_to_fort14` recursed into itself with the wrong
  arity. Now delegates to the module-level `write_fort14`. Roundtrip
  identity test catches any future regression.
- **B3**: `_elem_type` no longer treats vertex id `0` as a sentinel for
  "missing 4th vertex". Padded triangles are detected by row-pattern.
- **B4**: `_ensure_ccw_orientation` no longer applies the quad
  permutation `[0,3,2,1]` to padded-triangle rows in mixed-element
  meshes. Each CW row is classified before being flipped.
- **B5**: `interior_angles` quad branch guards against zero-length edges
  with the same `+1e-12` epsilon already used by the triangle branch and
  takes the real part of the result, mirroring MATLAB's `real(acosd(...))`.
- `read_from_fort14` accepts node indices written in float form
  (e.g. `1.000000`) so legacy fixtures load.

### Changed
- Bumped `requires-python` to `>=3.10` (drops EOL 3.8/3.9).
- `requirements.txt` reduced from a 60-line MADMESHR ML environment
  freeze to the three actual runtime dependencies (`numpy`, `scipy`,
  `matplotlib`).
- README quickstart rewritten around `chilmesh.examples.annulus()`;
  the broken `/kaggle/working/...` path is gone, the install command
  now includes `-r`, the stale "pip installation" TODO is removed,
  an alpha-status banner has been added, and the MADMESHR "Future Work"
  badge has been reframed as a "See also / Downstream projects"
  section that explicitly notes CHILmesh is agnostic to MADMESHR.

### Investigated, no fix needed
- **B2** (`signed_area` quad shoelace): the audit reported the formula
  used 4 cross terms instead of 8. On inspection the existing 4-row
  formula expands to the correct 8 terms, and the new random-quad
  cross-check test against an independent shoelace passes. Tests
  retained as a hardening guard.

### Known issues / deferred to 0.1.2
- The O(n²) Python loops in `_build_elem2edge` / `_build_edge2elem`
  make Block_O (~5,200 elements) slow to load. Vectorisation lands
  in 0.1.2.
- Coverage gating, Hypothesis property tests, mixed-element smoothing
  characterisation, Windows + Python 3.13 CI matrix.

[0.1.1]: https://github.com/domattioli/CHILmesh/releases/tag/v0.1.1
