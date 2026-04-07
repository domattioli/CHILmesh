# 0.1.1 overnight progress notes

Branch: `release/0.1.1` (off `main`).
Spec: `AUDIT_REPORT.md` on branch `audit/strategic-plan-2026-04`.
Date: 2026-04-07.

## What I did

All work is on `release/0.1.1`. One logical fix per commit, test-first
within each. Run `git log main..release/0.1.1 --oneline` for the order.

### Bugs fixed (audit B1-B11)
- **B1** — `CHILmesh.write_to_fort14` recursed into itself with the wrong
  arity. Now delegates to module-level `write_fort14`. The audit said
  "infinite recursion"; in practice it raises `TypeError` before
  recursing because the call passes 4 args to a 2-3-arg method. Same
  root cause, same fix. Regression: parametrized roundtrip identity test
  in `tests/test_fort14_roundtrip.py`.
- **B3** — `_elem_type` no longer treats vertex id `0` as a sentinel for
  "missing 4th vertex". Replaced the per-element loop with a vectorized
  "any pair of vertices in the row is equal" mask. Regression: synthetic
  quad with vertex 0 in slot 3 (`tests/test_elem_type.py`).
- **B4** — `_ensure_ccw_orientation` now classifies each CW row before
  flipping. Padded triangles get a triangle-only flip with the pad
  re-applied to slot 3; real quads get the original `[0,3,2,1]`
  permutation. Regression: `test_ccw_mixed_element_flips_padded_cw_triangle_correctly`,
  which fails on the previous code.
- **B5** — `interior_angles` quad branch now uses the same `+1e-12`
  epsilon guard as the triangle branch and takes `np.real` of the
  result. Regression: zero-length-edge quad
  (`test_quad_angles_no_nan_on_degenerate_quad`) plus angle-sum
  invariants over every fixture.
- **B6** — Added `src/chilmesh/utils/__init__.py`. Regression:
  `tests/test_package_imports.py::test_utils_subpackage_importable`.
- **B7** — Stripped `requirements.txt` from a 60-line MADMESHR ML
  freeze to the three runtime dependencies declared in `pyproject.toml`.
- **B8** — README now uses `pip install -e .` from a checkout instead of
  the broken `pip install requirements.txt`.
- **B9** — README quickstart rewritten around `chilmesh.examples.annulus()`,
  removing the `/kaggle/working/...` path that no pip user could run.
- **B10** — `pyproject.toml` gained `keywords`, `classifiers`,
  `[project.urls]`, `requires-python>=3.10`, and a
  `[project.optional-dependencies]` `dev` group.
- **B11** — Test count went from 6 to 58, and now parametrizes over all
  four shipped fixtures instead of just the annulus.

### B2 — investigated, no fix needed

The audit reported that the quad `signed_area` shoelace "uses 4 cross
terms instead of 8 → wrong by factor of ~2". On inspection the existing
4-row formula expands to the correct 8 cross terms:

```
0.5 * (x0(y1-y3) + x1(y2-y0) + x2(y3-y1) + x3(y0-y2))
  == 0.5 * sum_i (x_i*y_{i+1} - x_{i+1}*y_i)
```

Verified analytically, by hand on a unit square (area = 1.0) and a
random parallelogram (area = base × height), and by the new
`test_signed_area_quad_shoelace_matches_numpy` test which cross-checks
against an independent shoelace on randomly-generated convex quads.
**Recommend leaving the existing implementation alone.** The tests are
kept as a hardening guard.

If you disagree, the next move is to re-derive the formula explicitly
against MATLAB `.m:1900-1908` (which I don't have access to in this
environment).

### Discovered out-of-audit bugs

- `read_from_fort14` could not parse legacy fort.14 files where the
  node indices were written as floats like `1.000000` (this affected
  `structuredMesh1.14`, one of the four shipped fixtures). One-line
  fix: `int(float(...))`. Without this, `chilmesh.examples.structured()`
  would have raised, breaking the parametrized test gate.

### Infrastructure / packaging
- Moved all four `.fort.14` fixtures from `doc/domains/fort_14/` to
  `src/chilmesh/data/` and ship them via
  `[tool.setuptools.package-data]`. The directory now also has an
  `__init__.py` so it works as both a regular package and an
  `importlib.resources` source.
- Added `chilmesh.examples` with `annulus()`, `donut()`, `block_o()`,
  `structured()` factories, plus a `fixture_path()` helper. Imported
  from the package root so `chilmesh.examples.annulus()` works.
- Bumped `version` to `0.1.1` and `requires-python` to `>=3.10`.
- Added `tests/conftest.py` with parametrized `mesh` and
  `fresh_mesh` fixtures and a session-level memoization layer that
  patches `chilmesh.examples` so each fixture is loaded at most once
  per test run. Block_O alone takes ~30s on first load (the O(n²)
  `_build_elem2edge` is the culprit; deferred to 0.1.2 per audit).
- CI workflow staged at `ci/ci.yml.pending` (Py 3.10/3.11/3.12 ×
  ubuntu/macos): installs editable, runs pytest, builds sdist+wheel,
  runs `twine check`, smoke-tests `chilmesh.examples.annulus()` from
  the installed wheel in a fresh venv. **No coverage gate yet**
  (deferred to 0.1.2 per audit). **No release.yml yet** — that needs
  your morning approval. **Action required**: the GitHub PAT used by
  this session lacks the `workflow` scope, so I could not push the
  file directly into `.github/workflows/`. Move it manually:

  ```bash
  mkdir -p .github/workflows
  git mv ci/ci.yml.pending .github/workflows/ci.yml
  git commit -m "ci: enable GitHub Actions workflow"
  git push
  ```

  …or rerun this push from a session whose PAT has the `workflow`
  scope.
- Added `CHANGELOG.md`.
- Added a minimal `.gitignore` for build artefacts.

## What I did NOT do (deliberately, requires your approval)

- **Push to `main`.** All commits live on `release/0.1.1`.
- **Open a PR.** You said you wanted to read the diff and open it
  yourself.
- **Tag `v0.1.1`.**
- **Publish to PyPI.**
- **Yank `chilmesh==0.1.0` from PyPI** (audit recommends yes; this
  is irreversible-with-effort and needs your call).
- **Configure PyPI Trusted Publishing** (audit Open Question #1).
- **Add `.github/workflows/release.yml`** — depends on the trusted
  publishing decision above.

## Open question on `_mesh_layers` IE divergence (audit Q3)

I did **not** spend time on this. Read the audit guidance and the
existing Python `_mesh_layers` produces a disjoint cover with monotone-
shrinking layer sizes on the convex annulus and structured fixtures
(verified by `test_layers_disjoint_cover`). The implementation it picks
is "expand outward from boundary, peel two layers at a time (OE then
IE)". The existing tests pin its current behaviour as an invariant
rather than a value, which is the audit's preferred resolution. **TODO
for you**: decide whether you want this annotated in a comment block in
`_mesh_layers` as the chosen definition; I left the existing
implementation untouched.

## Test count

Audit asked for 8 new tests minimum. I shipped:

| File                              | Tests | Notes                              |
| --------------------------------- | ----- | ---------------------------------- |
| `tests/test_package_imports.py`   | 4     | B6 + examples module               |
| `tests/test_fort14_roundtrip.py`  | 4     | B1 over 4 fixtures                 |
| `tests/test_signed_area.py`       | 7     | B2 invariant + synthetic           |
| `tests/test_elem_type.py`         | 2     | B3 + mixed-element triangle pad    |
| `tests/test_ccw.py`               | 8     | B4 + idempotency over 4 fixtures   |
| `tests/test_interior_angles.py`   | 6     | B5 + sum invariants over fixtures  |
| `tests/test_invariants.py`        | 5×4=20| layer cover, connectivity, etc.    |
| `tests/test_layers_annulus.py`    | 6     | pre-existing, repointed at examples|
| **Total**                         | **57**|                                    |

## How to verify locally

```bash
git checkout release/0.1.1
pip install -e ".[dev]"
pytest -v
```

The full run including Block_O takes ~3 minutes on my box; without
Block_O it is under 10 seconds (`pytest -k 'not block_o'`).
