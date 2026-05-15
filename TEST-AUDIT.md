# CHILmesh Test Suite Audit

**Date:** 2026-05-15
**Branch:** `daily-issue-fixing`
**Issue:** [#110](https://github.com/domattioli/CHILmesh/issues/110)
**Scope:** Read-only audit. No code changes.

---

## TL;DR

- 27 test files, 239 `test_*` functions, **457 collected items** (parametrization expands ~218 cases).
- Suite runs in **~24s** on the editable install (Python 3.11, `MPLBACKEND=Agg`). 448 pass, 9 skipped, 0 failed, 12 warnings.
- **Line coverage 77%, branch coverage acceptable but with one large dead zone** (`CHILmesh.py` lines 1377-1424 and 1451-1609 — `_ordered_vertex_ring` + `angle_based_smoother`, both completely untested).
- Two CI workflows; only `python-package.yml` gates PRs. It exercises the matrix (Ubuntu/macOS × py3.10/3.11/3.12), builds wheels, and smoke-tests an installed wheel — good gating.
- No `TESTING.md` / `CONTRIBUTING.md`. The README test badge points at a workflow file that does not exist (`test.yml`).
- Custom marks (`@pytest.mark.slow`) are **unregistered** → warnings on every run.

---

## 1. Inventory & Layout

### File counts and class structure

| File | Tests | Classes | LOC |
|---|---:|---:|---:|
| `tests/test_plot_utils.py` | 104 | 0 | 217 |
| `tests/test_bridge_adapters.py` | 44 | 5 | 410 |
| `tests/test_smoothing.py` | 43 | 9+ | 520 |
| `tests/test_admesh_metadata.py` | 43 | 1 | 144 |
| `tests/test_advancing_front.py` | 26 | — | 349 |
| `tests/test_edge_map.py` | 23 | — | 226 |
| `tests/test_invariants.py` | 20 | — | (5 fns × 4 fixt) |
| `tests/test_performance_edge_building.py` | 18 | — | 184 |
| `tests/test_admesh_warmstart.py` | 15 | 5 | 377 |
| `tests/test_fast_init.py` | 14 | — | — |
| `tests/test_metadata_validation.py` | 12 | — | 194 |
| `tests/test_from_admesh_domain.py` | 11 | — | — |
| `tests/test_2dm_reader.py` | 10 | — | 215 |
| `tests/test_skeletonization_matlab_parity_external.py` | 9 | — | 154 |
| `tests/test_ccw.py` | 9 | — | — |
| `tests/test_skeletonization_invariant.py` | 8 | — | — |
| `tests/test_degeneracy.py` | 8 | — | 200 |
| `tests/test_signed_area.py` | 7 | — | — |
| `tests/test_layers_annulus.py` | 6 | — | 160 |
| `tests/test_interior_angles.py` | 6 | — | — |
| `tests/test_skeletonization_matlab_parity.py` | 5 | — | — |
| `tests/test_fort14_quad_roundtrip.py` | 5 | — | — |
| `tests/test_package_imports.py` | 4 | — | — |
| `tests/test_fort14_roundtrip.py` | 4 | — | — |
| `tests/test_elem_type.py` | 2 | — | — |
| `tests/test_readme_quickstart.py` | 1 | — | 123 |
| **Total** | **239 / 457 items** | — | 4 450 |

### Discoverability

- `pytest --collect-only` produces clean output — no collection errors.
- Naming is consistent: `tests/test_<topic>.py` and `test_<name>` functions.
- No `tests/unit/`, `tests/integration/`, `tests/e2e/` segmentation; everything is flat. Categorisation is implicit in filenames (e.g. `*_invariant`, `*_roundtrip`, `*_parity*`).

### Fixture infrastructure

- `tests/conftest.py:15` defines `FIXTURE_NAMES = ["annulus", "donut", "block_o", "structured", "quad_2x2"]`.
- A module-level cache memoises mesh loaders (`tests/conftest.py:20-39`) — the block_o load (~14s) hits once per session, not per test.
- Three top-level fixtures exported: `mesh`, `fresh_mesh`, `fixture_name`. Only `tests/test_invariants.py:18` uses `scope="module"`; everywhere else relies on the cache.

---

## 2. Coverage

Branch coverage was collected via `pytest --cov=src/chilmesh --cov-branch`. All numbers reflect the suite as-is (i.e. with the 9 skipped external-parity tests).

| Module | Stmts | Miss | Branch | BrPart | Cover |
|---|---:|---:|---:|---:|---:|
| `src/chilmesh/CHILmesh.py` | 893 | 225 | 396 | 41 | **73%** |
| `src/chilmesh/__init__.py` | 9 | 0 | 0 | 0 | 100% |
| `src/chilmesh/_vendor_admesh_truss.py` | 127 | 18 | 44 | 8 | 81% |
| `src/chilmesh/admesh_warmstart.py` | 163 | 49 | 48 | 3 | **67%** |
| `src/chilmesh/bridge.py` | 84 | 0 | 26 | 0 | 100% |
| `src/chilmesh/examples.py` | 23 | 0 | 2 | 1 | 96% |
| `src/chilmesh/mesh_topology.py` | 36 | 0 | 8 | 0 | 100% |
| `src/chilmesh/utils/plot_utils.py` | 199 | 19 | 78 | 9 | 88% |
| **TOTAL** | **1 536** | **311** | **602** | **62** | **77%** |

### Lowest-coverage hotspots (line ranges from `--cov-report=term-missing`)

| Rank | File | Range | Symbol(s) | Notes |
|---|---|---|---|---|
| 1 | `src/chilmesh/CHILmesh.py` | 1377-1424 | `_ordered_vertex_ring` | Helper for angle-based smoother. **Zero tests hit this.** |
| 2 | `src/chilmesh/CHILmesh.py` | 1451-1609 | `angle_based_smoother` (Zhou & Shimada) | **Zero direct tests.** The recent FEM-smoother work (issue #105) added 7 tests in `TestMixedElementFEMSmoother`, but only exercises FEM path, not the standalone angle-based smoother. |
| 3 | `src/chilmesh/admesh_warmstart.py` | 77-136 | `_vbnd_from_mesh` / boundary helpers | Internal helpers behind public `optimize_with_admesh_truss*`. Tested only via the public API, leaving the boundary-build paths thinly covered. |
| 4 | `src/chilmesh/CHILmesh.py` | 698-699, 713 | `boundary_node_indices` early-returns | One branch (no boundary edges) untested. |
| 5 | `src/chilmesh/CHILmesh.py` | 2130-2139 | `pinch_points` (tail of method) | Width-threshold edge cases. |
| 6 | `src/chilmesh/CHILmesh.py` | 908-914 | `_mesh_layers` fallback | Likely degenerate input path. |
| 7 | `src/chilmesh/CHILmesh.py` | 652-660, 679-689 | `_build_edge2elem` / `_build_vert2edge` fallbacks | Exception/edge branches. |
| 8 | `src/chilmesh/utils/plot_utils.py` | 320-326, 340-350 | Plot helpers | Two tail branches. |
| 9 | `src/chilmesh/_vendor_admesh_truss.py` | 280-292 | Vendored numerics tail | Likely an unreachable-on-our-fixtures path. |
| 10 | `src/chilmesh/admesh_warmstart.py` | 388-389, 513-514, 600 | Misc guard clauses | Negative-input branches. |

### Public-API coverage map

Public symbols (per `src/chilmesh/__init__.py:14-22` and `CHILmesh` class methods) cross-referenced against test references:

| Symbol | Hits in `tests/` | Coverage status |
|---|---:|---|
| `CHILmesh.read_from_fort14` | 92 IO refs (combined `read_from_fort14` / `write_to_fort14` / `write_fort14`) | Strong (roundtrip × 4 fixtures) |
| `CHILmesh.write_to_fort14` | (above) | Strong |
| `CHILmesh.read_from_2dm` | 10 in `test_2dm_reader.py` | OK |
| `CHILmesh.from_admesh_domain` | 11 in `test_from_admesh_domain.py` | OK |
| `CHILmesh.smooth_mesh` | 43 in `test_smoothing.py` | Strong for `method='fem'`; **no test for `method='angle'` / direct path** |
| `CHILmesh.angle_based_smoother` | 0 direct | **Untested** (see Finding F1) |
| `CHILmesh.direct_smoother` | 0 direct | **Untested standalone** (only via FEM dispatcher fallback paths) |
| `CHILmesh.elem_quality` | covered via smoothing + degeneracy tests | OK |
| `CHILmesh.interior_angles` | 6 in `test_interior_angles.py` | OK |
| `CHILmesh.get_layer` / `elements_in_layer` | covered via layer / invariant tests | OK |
| `CHILmesh.get_vertex_edges` / `get_vertex_elements` | covered via invariants | OK |
| `CHILmesh.advancing_front_*` | 37 refs in `test_advancing_front.py` | Strong |
| `CHILmesh.pinch_points` | covered | OK (tail uncovered, F5 above) |
| `CHILmesh.copy` | indirect via `fresh_mesh` fixture | Weak — no dedicated invariants test |
| `CHILmesh.change_points` | a handful | OK |
| `CHILmesh.admesh_metadata` | 43 in `test_admesh_metadata.py` | Strong |
| `EdgeMap` | 23 in `test_edge_map.py` | Strong |
| `optimize_with_admesh_truss` / `_arrays` | 15 in `test_admesh_warmstart.py` | OK |
| `write_fort14` (module-level) | covered | OK |
| `bridge.*` adapters | 44 in `test_bridge_adapters.py` | Strong |

---

## 3. Quality Smells

- **Type-only / existence assertions** — at least 45 sites use `isinstance(...)` or `is not None` as the only assertion. Hotspot: `tests/test_plot_utils.py` (≈18 occurrences from line 30 onward) — many tests assert only that `plot()` returns `(Figure, Axes)` without checking *what* was drawn. Examples: `tests/test_plot_utils.py:30`, `tests/test_plot_utils.py:37`, `tests/test_edge_map.py:72-73`, `tests/test_admesh_warmstart.py:197-211`, `tests/test_degeneracy.py:114`. (F2)
- **Assertion messages** — ~50% of asserts (218 / 439) carry an explanatory message. The remainder fail with bare `assert <expr>`, harder to triage from CI logs. Global pattern across `tests/`. (F6)
- **No tautological asserts** (no `assert True` / `assert 1 == 1`) — clean.
- **Over-mocking** — no `unittest.mock` usage detected in `tests/`; everything goes through the real `CHILmesh` class. Good.

---

## 4. Speed & Flakiness

### Top 15 slowest tests (`pytest --durations=15`, fresh editable install)

| Time | Test |
|---:|---|
| 8.69s | `tests/test_readme_quickstart.py::test_readme_annulus_regenerates_and_validates` |
| 0.77s | `tests/test_admesh_warmstart.py::TestQualityAndPerformance::test_determinism` |
| 0.75s | `tests/test_smoothing.py::TestSmootherIntegration::test_smooth_mesh_fem_triangles[block_o]` |
| 0.73s | `tests/test_smoothing.py::TestMixedSmoother::test_fem_smoother_triangle_mesh_backward_compat[block_o]` |
| 0.68s | `tests/test_smoothing.py::TestTriangleSmoother::test_fem_smoother_triangle_preserves_boundary[block_o]` |
| 0.67s | `tests/test_smoothing.py::TestTriangleSmoother::test_fem_smoother_triangle_preserves_z[block_o]` |
| 0.67s | `tests/test_smoothing.py::TestTriangleSmoother::test_fem_smoother_triangle_interior_changes[block_o]` |
| 0.52s | `tests/test_plot_utils.py::test_plot_label_all[structured]` |
| 0.46s | `tests/test_admesh_warmstart.py::TestBoundaryPreservation::test_vbnd_multiple_domains` |
| 0.46s | `tests/test_plot_utils.py::test_plot_label_all[annulus]` |
| 0.43s | `tests/test_admesh_warmstart.py::TestBoundaryPreservation::test_vbnd_annulus_array_form` |
| 0.41s | `tests/test_admesh_warmstart.py::TestBoundaryPreservation::test_vbnd_annulus_chilmesh_form` |
| 0.41s | `tests/test_plot_utils.py::test_plot_label_edge_only[structured]` |
| 0.41s | `tests/test_admesh_warmstart.py::TestQualityAndPerformance::test_performance_annulus` |
| 0.40s | `tests/test_admesh_warmstart.py::TestNonDegradation::test_enforce_non_degradation_true` |

Single outlier (`test_readme_annulus_regenerates_and_validates`) drives ~8s/24s of wall time. Everything else is sub-second after the fixture cache warms.

### Flakiness candidates

| Location | Risk | Notes |
|---|---|---|
| `tests/test_admesh_warmstart.py:303` | **High** | `np.random.uniform(-0.95, 0.95, ...)` — unseeded. Test passes today but a `Delaunay` triangulation built from non-deterministic interior points can flip element IDs / quality outcomes across runs. **F3.** |
| `tests/test_admesh_warmstart.py:367-371` | Low | `time.time()` used for performance gate. Tolerance not visible in this audit — verify the threshold is generous. **F8.** |
| `tests/test_skeletonization_matlab_parity_external.py:117` | Low | Runs only when `CHILMESH_RUN_EXTERNAL_PARITY=1` + `admesh-domains` installed; otherwise skipped. Acceptable opt-in. |
| Mesh-IO tests | Low | `tests/test_2dm_reader.py` uses `tempfile.NamedTemporaryFile(delete=False)` with manual `unlink()` in `finally` (`tests/test_2dm_reader.py:25, 36, 52, 63, 80, 91, 104, 115, 123, 134, 142, 153, 162, 173, 183, 194, 205, 216`). Works but `tmp_path` fixture is more idiomatic and crash-safe. **F7.** |
| Floating-point tolerances | Low | Roundtrip tests (`tests/test_fort14_roundtrip.py`, `tests/test_fort14_quad_roundtrip.py`) compare geometry directly — confirm explicit `atol`/`rtol` (not audited line-by-line here). |

---

## 5. Redundancy & Drift

- **Duplicate `FIXTURES` constants** — `["annulus", "donut", "block_o", "structured"]` is defined inline in `tests/test_fort14_roundtrip.py:17`, `tests/test_interior_angles.py:16`, `tests/test_ccw.py:16` (and at least `test_signed_area.py`, `test_layers_annulus.py`). All four omit `quad_2x2`, which `tests/conftest.py:15` includes. Drift relative to conftest. **F4.**
- **Skip markers** — 12 hits total; all justified:
  - `tests/test_skeletonization_invariant.py:44` — trivially-holds skip for small fixtures (fine).
  - `tests/test_skeletonization_matlab_parity_external.py:110-150` — external-parity gates (opt-in via env var, fine).
  - `tests/test_advancing_front.py:109-110, 180` — runtime conditional skips (mesh-shape dependent, fine).
  - `tests/test_admesh_warmstart.py:35` — `admesh_unavailable` skip when the external `admesh` package is missing (fine).
- **TODO / FIXME** — only a single FIXME-style note: `tests/test_skeletonization_matlab_parity_external.py:145` flagging that `None` entries in `MATLAB_REFERENCE_LAYER_COUNTS` (8 of them, `:59-71`) need ground-truth captured. Tracked-but-pending coverage debt.
- **Test recency** — earliest test commit `2026-04-27`; nothing >19 days stale at audit time. No stagnant tests.

---

## 6. External-Dependency Markers

- **`@pytest.mark.slow` is unregistered.** Two sites: `tests/test_performance_edge_building.py:122`, `tests/test_admesh_warmstart.py:359`. Produces `PytestUnknownMarkWarning` on every run and silently ignores any `-m "not slow"` intent. **F9.**
- **No filesystem / network markers needed.** All tests run against in-tree fixtures (`src/chilmesh/data/*.fort.14`, `*.14`). Total fixture payload: ~415 KB.
- **External-parity test** properly env-gated via `CHILMESH_RUN_EXTERNAL_PARITY` (`tests/test_skeletonization_matlab_parity_external.py`).
- **ADMESH-Truss vendor** runs entirely in-tree (`src/chilmesh/_vendor_admesh_truss.py`); the `admesh_warmstart` tests gracefully skip when the external `admesh` package is absent (`tests/test_admesh_warmstart.py:35`).

---

## 7. CI Gating

| Workflow | Trigger | Gates PR? | Matrix | Notes |
|---|---|---|---|---|
| `.github/workflows/python-package.yml` | push (`main`, `release/**`), PR to `main` | **Yes** | `ubuntu-latest` × `macos-latest` × py `3.10/3.11/3.12` | Runs full `pytest -v`, builds sdist + wheel, runs `twine check`, smoke-tests the installed wheel by importing `chilmesh` and instantiating `examples.annulus()`. Strong gate. |
| `.github/workflows/publish-pypi.yml` | release `published`, manual dispatch | No | py 3.11 only | Build + `twine upload` only. Does **not** rerun tests on release — relies on the prior CI run on the release commit. Acceptable but worth flagging. |

- **No Windows runner** despite `pyproject` claims Python 3.10+ portability. Constitution §5.2 lists "Python 3.10/3.11/3.12 × Ubuntu/macOS" — but the constitution also says **§4.1 "all platforms"** for the release checklist. Decide: add Windows or update constitution. **F10.**
- **Coverage not gated.** No `--cov-fail-under=N` on `pytest` invocation; coverage drops would not break CI.
- **Nightly job:** none. Acceptable while suite stays sub-30s, but flag for future when `slow`-marked tests grow.

---

## 8. Test-Data Hygiene

| Fixture | File | Size | Notes |
|---|---|---:|---|
| annulus | `src/chilmesh/data/annulus_200pts.fort.14` | 24 KB | Small, fast. |
| donut | `src/chilmesh/data/donut_domain.fort.14` | 12 KB | Small, fast. |
| structured | `src/chilmesh/data/structuredMesh1.14` | 40 KB | Quad mesh. |
| block_o | `src/chilmesh/data/Block_O.14` | 336 KB | Large; cached via `_MESH_CACHE` (`tests/conftest.py:20-39`). |
| quad_2x2 | `src/chilmesh/data/quad_2x2.fort.14` | 0.4 KB | Toy fixture; **not parametrised by `FIXTURES` constants in several test files** (see F4). |

- All fixtures committed to repo; <1 MB total. Healthy.
- No golden-output files (no expected-rendering blobs). Numerical invariants are computed inline.

---

## 9. Framework Hygiene

- `tests/conftest.py` is small (62 LOC) and single-purpose — no sprawl.
- Module-level mesh cache is documented (`tests/conftest.py:17-22`) and explicit about the mutation contract ("Tests that mutate a mesh in place must call `.copy()` first"). Good.
- Only one explicit `scope="module"` fixture (`tests/test_invariants.py:18`); everywhere else relies on the cache to amortise.
- **No `pytest.ini` / `[tool.pytest.ini_options]` markers registered.** `pyproject.toml` `[tool.pytest.ini_options]` only sets `testpaths`. Missing `markers = [...]` is the proximate cause of F9.

---

## 10. Docs & Onboarding

- **No `TESTING.md`.** No `CONTRIBUTING.md` either. New contributors must infer the test entrypoint from `.claude/CLAUDE.md:138-150` or the CI workflow.
- **README badge is broken:** `README.md:18` links `actions/workflows/test.yml` but the workflow file is `python-package.yml`. The badge image will 404. **F11.**
- Local-vs-CI parity is otherwise tight: `MPLBACKEND=Agg` is set in CI (`python-package.yml`) and required locally for plot tests (`tests/test_plot_utils.py:7`).
- `pytest` runs cold for a new contributor with `pip install -e ".[dev]" && python -m pytest -v` — verified during this audit.

---

## Findings (by severity)

Each finding is keyed `F<n>` for cross-referencing in the backlog below.

### High

- **F1 — `angle_based_smoother` and `_ordered_vertex_ring` are entirely uncovered.** `src/chilmesh/CHILmesh.py:1370-1424` and `:1426-1609`. ~190 lines of production code, including the Zhou & Shimada angle-based smoother referenced explicitly by the constitution (`.planning/constitution.md:340-346`, "DOMsmooth Hybrid Fallback"). Recent FEM-smoother work landed tests for the FEM path; the angle path is the documented fallback strategy for mixed meshes and has no regression guard.
- **F2 — `tests/test_plot_utils.py` is shape-only.** 104 collected items, but the dominant assertion pattern is `isinstance(fig, plt.Figure)` (≈18 sites starting `tests/test_plot_utils.py:30`). The tests verify that `plot()` doesn't crash; they do not verify that anything correct is drawn (no artist counts, no axis-limit checks, no colormap checks). Half of the file's coverage value is essentially smoke testing.
- ~~**F3 — Unseeded `np.random` in `tests/test_admesh_warmstart.py:303`.**~~ **WITHDRAWN 2026-05-15:** re-inspection shows `np.random.seed(42)` is set on `tests/test_admesh_warmstart.py:289` inside the same method (`TestExtensibility.test_array_form_source_agnostic`), 14 lines before the `np.random.uniform` call on `:303`. The seed is in scope and effective. Audit was overzealous. The underlying universal-convention proposal (autouse fixture to fail on unseeded `np.random.*`) was nonetheless filed upstream on `domattioli/DomI#63` as a hardening lint.

### Medium

- **F4 — Drift between `tests/conftest.py:15` and per-file `FIXTURES` constants.** Conftest lists `["annulus", "donut", "block_o", "structured", "quad_2x2"]`; `tests/test_fort14_roundtrip.py:17`, `tests/test_interior_angles.py:16`, `tests/test_ccw.py:16` (and at least `test_signed_area.py`, `test_layers_annulus.py`) omit `quad_2x2`. Either harmonise (preferred — define once in conftest and re-export) or document the omission. The quad-only fixture is new (`v0.2.x`) and should be exercised wherever invariants are claimed to hold for "all fixtures".
- **F5 — No standalone test for `CHILmesh.copy()`.** The fixture cache relies on `.copy()` for mutation safety (`tests/conftest.py:21-22`), but there's no invariants test confirming a deep, independent copy. If `.copy()` ever degrades to a shallow copy, every mutation test silently begins corrupting the cache.
- **F6 — Assertion messages absent on ~half of assertions.** 218 of 439 assertions carry an explanatory message; the rest fail with bare `assert <expr>`. For numerical-equality assertions in particular (geometry roundtrip, quality invariants), `assert a == b, f"got {a}, expected {b}"` saves CI-log triage time.
- **F7 — `tempfile.NamedTemporaryFile(delete=False)` + manual `unlink()` in `tests/test_2dm_reader.py`.** 9 occurrences (`tests/test_2dm_reader.py:25, 52, 80, 104, 123, 142, 162, 183, 205`), all paired with `try/finally`. Works, but `tmp_path` is the pytest-idiomatic choice and is crash-safe if a test errors before the `finally`.
- **F8 — Wall-clock performance gate in `tests/test_admesh_warmstart.py:367-371`.** `time.time()` deltas in CI are noisy across runners. Either widen the gate to clearly cover macOS variance, or move the perf gate behind `@pytest.mark.slow` and run it on a single canonical runner.

### Low

- **F9 — `@pytest.mark.slow` is unregistered.** Two sites (`tests/test_performance_edge_building.py:122`, `tests/test_admesh_warmstart.py:359`) emit `PytestUnknownMarkWarning` on every run. Register in `[tool.pytest.ini_options]` `markers = ["slow: marks slow tests..."]`.
- **F10 — No Windows runner.** Matrix is Ubuntu × macOS only (`.github/workflows/python-package.yml:14-16`). Either add Windows or amend `.planning/constitution.md` to scope "all platforms" to Linux + macOS.
- **F11 — README test badge points to non-existent workflow.** `README.md:18` references `actions/workflows/test.yml`; the file is `python-package.yml`. Badge will render as "unknown".
- **F12 — 8 untested `MATLAB_REFERENCE_LAYER_COUNTS` entries.** `tests/test_skeletonization_matlab_parity_external.py:59-71` carries `None` placeholders for lake-erie, delaware-bay, chesapeake-bay, great-lakes, lake-michigan, baranja-hill (×2). Tracking comment is good (`:145-152`) but the gap should land in an issue, not just inline.
- **F13 — No `branches` / lines coverage gate in CI.** Dropping below 77% would not break a build today.
- **F14 — `_skeletonize` boundary path (`src/chilmesh/CHILmesh.py:908-914`) and tail of `pinch_points` (`:2130-2139`) are uncovered.** Small but algorithmically interesting branches.

---

## Backlog Status (2026-05-15)

Resolved in this session (autonomous batch on `daily-issue-fixing`):

| ID | Status | Commit | Notes |
|---|---|---|---|
| F3 | Withdrawn | — | False alarm; seed already set 14 lines upstream of the cited call. |
| F4 | ✅ Resolved | `63bf016` | `TRI_FIXTURE_NAMES` exported from `conftest.py`; 4 test files now import. |
| F9 | ✅ Resolved | `20ac025` | `slow` marker registered in `[tool.pytest.ini_options]`. Warning count 12 → 10. |
| F11 | ✅ Resolved | `2427103` | README Tests badge points at `python-package.yml`. |

Outstanding backlog (per-finding follow-up tickets pending):

## Prioritized Backlog (Top 10)

Each item is intended to become a follow-up issue. The finding ID points back to the section above.

1. **[F1] Add regression tests for `angle_based_smoother` and `_ordered_vertex_ring`.** Cover at minimum: convergence on annulus, boundary pinning, no-element-collapse on a mixed-element fixture, idempotency on a converged mesh. Target ≥80% line coverage of `CHILmesh.py:1370-1609`.
2. **[F2] Promote `tests/test_plot_utils.py` from smoke tests to assertion tests.** For each `plot()` variant, assert: artist count > 0, axis limits enclose mesh extent, colormap matches request when `elem_color=` is a scalar field. Keep parametrization over `SMALL_FIXTURES`.
3. **[F3] Seed every `np.random` call in tests.** Audit `tests/`, replace unseeded calls (currently `tests/test_admesh_warmstart.py:303`) with `np.random.default_rng(<fixed seed>)`. Add a `conftest.py` autouse fixture that fails the test if it observes unseeded numpy use (optional second deliverable).
4. **[F4] Centralise the fixtures list.** Export `FIXTURE_NAMES` from `tests/conftest.py` (or a `tests/_fixtures.py`) and replace per-file `FIXTURES = [...]` constants. Decide whether `quad_2x2` belongs in invariants/roundtrip/CCW/interior-angles tests; either include or document the exclusion in the constant's docstring.
5. **[F1/F5] Add `test_copy.py` with deep-copy invariants.** Assert: mutating `copy().points` does not alter the original; `copy()` preserves `n_elems`, `n_verts`, `n_layers`, `adjacencies` dict identity (not aliased).
6. **[F9] Register pytest markers.** Add `[tool.pytest.ini_options] markers = ["slow: ..."]` to `pyproject.toml`. Confirm `pytest -m "not slow"` filters as expected.
7. **[F11] Fix README test badge.** Update `README.md:18` to reference `python-package.yml` (or add a `test.yml` redirect workflow). Verify the rendered badge.
8. **[F6] Add assertion messages to numerical comparisons.** Pass over `tests/test_fort14_roundtrip.py`, `tests/test_fort14_quad_roundtrip.py`, `tests/test_ccw.py`, `tests/test_signed_area.py`, `tests/test_invariants.py`. Each numeric assertion should carry a `f"..."` message with the fixture name and the diff.
9. **[F7] Migrate `tests/test_2dm_reader.py` to `tmp_path`.** Replace the 9 `tempfile.NamedTemporaryFile(delete=False)` blocks with `tmp_path / "<name>.2dm"`. Removes the try/finally boilerplate and is crash-safe.
10. **[F12] Capture remaining MATLAB reference layer counts.** File a tracking issue per missing entry in `tests/test_skeletonization_matlab_parity_external.py:59-71` (or one umbrella issue) with a checklist. Out-of-scope for the audit itself but blocks closing the "MATLAB parity" promise in `.planning/constitution.md:347-352`.

---

## "Do Nothing" List

These were considered and intentionally not added to the backlog.

- **Add Windows CI runner now.** macOS already exercises a non-Linux path. Adding Windows triples runner-minutes for a library whose downstream consumers (MADMESHR/ADMESH/ADMESH-Domains) run on Linux + macOS. Revisit if a Windows user reports a real bug.
- **Coverage-fail-under gate in CI.** Current 77% line / 89% branch is reasonable; adding a hard gate before F1 lands would force either a coverage drop or a sham fix. Land F1 first, then set the gate at 80%.
- **Split `tests/` into `unit/`/`integration/`/`e2e/`.** The flat layout is small (27 files), discoverable, and naming is consistent. Reorg is busywork without a clear payoff. Revisit once the suite passes ~50 files.
- **Convert function-style tests to `unittest.TestCase`.** No mocks, no setup/teardown complexity; class-based tests in `test_smoothing.py` / `test_bridge_adapters.py` are organisational only. The mixed style is fine.
- **Add property-based testing (Hypothesis).** Numerical-mesh invariants would benefit *in principle*, but the cost of stable shrinkers for meshes is significant. Defer until a concrete bug shows that example-based parametrization missed a case.
- **Tighten the existing fort.14 roundtrip tolerance.** The current tolerance survived the v0.1.1 → v0.2.0 perf rewrite without floating-point drift. Tightening it pre-emptively trades stability for no observed bug.

---

## Success Criteria (issue #110)

- [x] `TEST-AUDIT.md` committed on `daily-issue-fixing`.
- [x] Each finding has a `file:line` reference (or `global` tag) — see F1-F14 above.
- [x] Backlog items each link back to a finding.
- [x] No code changes in this issue's scope.

---

## Audit Methodology Notes

- Suite executed once with `MPLBACKEND=Agg python -m pytest --durations=15 -q` and once with `--cov=src/chilmesh --cov-branch --cov-report=term-missing` against the editable install (`pip install -e ".[dev]"`) under Python 3.11.15.
- All findings are derived from static reads of `tests/`, `pyproject.toml`, `.github/workflows/`, and runtime output of pytest. No test was modified.
- The audit deliberately did **not** read `.planning/MODERNIZATION_LESSONS_LEARNED.md` or other planning docs beyond `.claude/CLAUDE.md` and `.planning/constitution.md` — to keep the audit grounded in the suite as it stands, not in what it was promised to be.
