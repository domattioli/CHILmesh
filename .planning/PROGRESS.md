# 0.1.1 overnight progress notes

Branch: `release/0.1.1` (off `main`).
Spec: `AUDIT_REPORT.md` on branch `audit/strategic-plan-2026-04`.
Date: 2026-04-07.

## What I did

All work on `release/0.1.1`. One logical fix per commit, test-first within each. Run `git log main..release/0.1.1 --oneline` for order.

### Bugs fixed (audit B1-B11)
- **B1** — `CHILmesh.write_to_fort14` recursed into itself with wrong arity. Now delegates to module-level `write_fort14`. Regression: parametrized roundtrip identity test in `tests/test_fort14_roundtrip.py`.
- **B3** — `_elem_type` no longer treats vertex id `0` as sentinel for "missing 4th vertex". Replaced per-element loop with vectorized "any pair of vertices in row is equal" mask. Regression: synthetic quad with vertex 0 in slot 3 (`tests/test_elem_type.py`).
- **B4** — `_ensure_ccw_orientation` now classifies each CW row before flipping. Padded triangles get triangle-only flip with pad re-applied to slot 3; real quads get original `[0,3,2,1]` permutation. Regression: `test_ccw_mixed_element_flips_padded_cw_triangle_correctly`.
- **B5** — `interior_angles` quad branch uses same `+1e-12` epsilon guard as triangle branch; takes `np.real` of result. Regression: zero-length-edge quad (`test_quad_angles_no_nan_on_degenerate_quad`) + angle-sum invariants over every fixture.
- **B6** — Added `src/chilmesh/utils/__init__.py`. Regression: `tests/test_package_imports.py::test_utils_subpackage_importable`.
- **B7** — Stripped `requirements.txt` from 60-line MADMESHR ML freeze to 3 runtime dependencies declared in `pyproject.toml`.
- **B8** — README uses `pip install -e .` from checkout instead of broken `pip install requirements.txt`.
- **B9** — README quickstart rewritten around `chilmesh.examples.annulus()`, removing `/kaggle/working/...` path.
- **B10** — `pyproject.toml` gained `keywords`, `classifiers`, `[project.urls]`, `requires-python>=3.10`, `[project.optional-dependencies]` dev group.
- **B11** — Test count went from 6 to 58; now parametrizes over all four fixtures.

### B2 — investigated, no fix needed

Audit reported quad `signed_area` shoelace "uses 4 cross terms instead of 8 → wrong by factor of ~2". On inspection the existing 4-row formula expands to correct 8 cross terms:

```
0.5 * (x0(y1-y3) + x1(y2-y0) + x2(y3-y1) + x3(y0-y2))
  == 0.5 * sum_i (x_i*y_{i+1} - x_{i+1}*y_i)
```

Verified analytically, by hand on unit square (area = 1.0) and random parallelogram, and by `test_signed_area_quad_shoelace_matches_numpy`. **Leave existing implementation alone.** Tests kept as hardening guard.

### Discovered out-of-audit bugs

- `read_from_fort14` could not parse legacy fort.14 files where node indices written as floats like `1.000000` (affected `structuredMesh1.14`). One-line fix: `int(float(...))`.

### Infrastructure / packaging
- Moved all four `.fort.14` fixtures from `doc/domains/fort_14/` to `src/chilmesh/data/`; ship via `[tool.setuptools.package-data]`; added `__init__.py`.
- Added `chilmesh.examples` with `annulus()`, `donut()`, `block_o()`, `structured()` factories + `fixture_path()` helper.
- Bumped `version` to `0.1.1`; `requires-python` to `>=3.10`.
- Added `tests/conftest.py` with parametrized `mesh` + `fresh_mesh` fixtures; session-level memoization (each fixture loaded at most once per run).
- CI workflow staged at `ci/ci.yml.pending` (Py 3.10/3.11/3.12 × ubuntu/macos). **Action required**: move manually:

  ```bash
  mkdir -p .github/workflows
  git mv ci/ci.yml.pending .github/workflows/ci.yml
  git commit -m "ci: enable GitHub Actions workflow"
  git push
  ```

- Added `CHANGELOG.md`; minimal `.gitignore` for build artefacts.

## What I did NOT do (deliberately, requires your approval)

- Push to `main` — all commits live on `release/0.1.1`
- Open PR — you said you wanted to read diff and open yourself
- Tag `v0.1.1`
- Publish to PyPI
- Yank `chilmesh==0.1.0` from PyPI (irreversible; needs your call)
- Configure PyPI Trusted Publishing
- Add `.github/workflows/release.yml`

## Open question on `_mesh_layers` IE divergence (audit Q3)

Not addressed. Existing Python `_mesh_layers` produces disjoint cover with monotone-shrinking layer sizes on convex annulus and structured fixtures (verified by `test_layers_disjoint_cover`). Implementation: "expand outward from boundary, peel two layers at a time (OE then IE)". Tests pin current behaviour as invariant. **TODO:** decide whether to annotate this in a comment block in `_mesh_layers` as the chosen definition.

## Test count

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

Full run including Block_O ~3 minutes; without Block_O under 10 seconds (`pytest -k 'not block_o'`).
