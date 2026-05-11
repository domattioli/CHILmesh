# CHILmesh Strategic Audit & Phased Plan

**Date:** 2026-04-07
**Method:** 8 independent auditor agents (Fidelity, Release Engineer, QA, Downstream Integration, Big Picture, Optimizer, Realist, Salesman) ran two rounds of independent analysis followed by mediated convergence.
**Output:** This document.

---

## Executive summary

CHILmesh 0.1.0 on PyPI ships **three latent bugs that corrupt user data on first use**: recursive `write_to_fort14` (infinite loop), wrong quad shoelace in `signed_area` (corrupts CCW orientation), and `_elem_type` misclassification of quads with 4th vertex == index 0. README has broken install command and kaggle example path.

Audit converged on **three-release plan** with hard test gates:

- **0.1.1 — "Honest hotfix"** (~12-15h): fix bugs, ship test fixtures via `importlib.resources`, add 8 tests, set up CI, fix README. Yank 0.1.0.
- **0.1.2 — "Test gate + perf"** (~20-30h): expand to ~25 tests including Hypothesis property layer, vectorize adjacency hotspots, add coverage gating.
- **0.2.0 — "General-purpose editing API"** (~30-50h): add `boundary_loops()`, `compact()`, `is_geometrically_valid()`, `vert2vert()`, `elem_centroid()`, `edge_length()`, `locate_point()` — all multi-consumer, library-neutral.

**Strategic identity (Big Picture):** CHILmesh = research-grade Python port of QuADMESH layer abstraction, maintained as substrate for downstream research, with first-class ADCIRC `.fort.14` I/O. Not `gmsh` replacement. Not 3D library. Not general mesh generator. Delaunay-from-scratch and native GMSH I/O cut; thin `meshio` adapter replaces both.

**Test policy:** No release tagged without all tests green and bug-regression suite covering every prior shipped bug.

---

## Convergent findings — bugs in 0.1.0

Independently surfaced by multiple auditors. No debate.

| # | File:line | Bug | Auditors | Severity |
|---|---|---|---|---|
| B1 | `src/chilmesh/CHILmesh.py:646` | `write_to_fort14` calls itself instead of module-level `write_fort14` at line 904 → infinite recursion on first call | Fidelity, QA, Realist | **Release blocker** |
| B2 | `src/chilmesh/CHILmesh.py:171` | `signed_area` quad shoelace uses 4 cross terms instead of 8 → wrong by factor of ~2 → corrupts CCW orientation for any quad mesh | Fidelity | **Release blocker** |
| B3 | `src/chilmesh/CHILmesh.py:208` | `_elem_type` misclassifies any quad whose 4th vertex is index 0 (sentinel collision with valid vertex id) | Fidelity | **Release blocker** |
| B4 | `src/chilmesh/CHILmesh.py:119` | `_ensure_ccw_orientation` flips mixed-element meshes incorrectly (applies quad permutation `[0,3,2,1]` to padded triangles) | Fidelity | High |
| B5 | `src/chilmesh/CHILmesh.py:716` | `interior_angles` quad branch has no `real()` / NaN guard; degenerate concave quads emit NaN that bypasses `elem_quality`'s `<= 359.99` zero-out (NaN comparisons are False) | Fidelity | High |
| B6 | `src/chilmesh/utils/` | Missing `__init__.py` → wheel may drop `plot_utils.py` depending on setuptools version | Release Engineer | High |
| B7 | `requirements.txt` | Contaminated with `torch`, `stable_baselines3`, `gymnasium`, `pandas`, `networkx`, `shapely`, `rtree` (~2 GB of MADMESHR ML deps) — none are runtime deps per `pyproject.toml` | Release, Realist | High |
| B8 | `README.md:67` | `pip install requirements.txt` (missing `-r`) | Release, Realist | Medium |
| B9 | `README.md:97` | Example path `/kaggle/working/CHILmesh/...` broken for any pip user | Release, Realist, Salesman | Medium |
| B10 | `pyproject.toml` | Missing `keywords`, `classifiers`, `[project.urls]` → PyPI page barren and unsearchable | Release, Salesman | Medium |
| B11 | `tests/` | Only 1 test file, only 1 of 4 available `.fort.14` fixtures used; 0 tests for write, quality, smoothing, mixed-element, I/O round-trip | QA, Realist | High |

**Note on B1:** Recursive `write_to_fort14` never successfully executed by any test. QA's `test_fort14_roundtrip_identity` would have caught it on day one — single highest-ROI test in entire plan.

---

## Strategic identity (Big Picture)

Big Picture auditor surfaced unstated tension: CHILmesh performing three identities simultaneously — (a) faithful Python port of 2017 thesis, (b) general-purpose 2D mesh library competing with `meshio`/`pygmsh`/`gmsh`, (c) substrate for author's downstream research line (MADMESHR, hydrodynamic ML). These imply incompatible priorities.

**Decision (adopted across audit):** identity = **(a) preserved + (c) actively maintained**. CHILmesh = research-grade Python implementation of QuADMESH layer abstraction, packaged cleanly enough to serve downstream research and ADCIRC users. Will *not* compete with `gmsh` on Delaunay generation or with `meshio` on format breadth.

**Implications:**

1. **Cut from roadmap:** Delaunay-from-scratch (`scipy.spatial.Delaunay` exists), native GMSH I/O (delegate to `meshio`), full MATLAB plotting parity, units conversion, `delaunayTriangulation` object I/O.
2. **Stay in scope:** layer-based decomposition (unique selling point), tri/quad/mixed first-class support, ADCIRC `.fort.14` I/O, FEM + angle-based smoothing, mesh quality/skewness metrics, clean editing API.
3. **README reframing:** "Future Work → MADMESHR" badge inverts dependency arrow. Reframe as "See also / Downstream projects." Add non-coupling statement.
4. **No MADMESHR-specific concepts** in CHILmesh APIs. Library-neutral naming. Anything serving only RL agents: out.

**3-year north star:** "MADMESHR (and 1-2 other downstream research projects) cite and depend on `chilmesh>=1.0` from PyPI; layer-based mesh abstraction from thesis available to anyone who reads QuADMESH paper; API has not broken in 18 months."

---

## Joint Phase 1 — 0.1.1 "Honest hotfix"

**Estimated work:** 12-15 hours. Realistically one weekend.

**Goal:** Ship release where `pip install chilmesh` followed by README quickstart works on clean machine, with regression tests for every shipped bug, and CI gates future merges.

### Code changes (correctness)

1. **Fix B1** — `write_to_fort14` (`CHILmesh.py:646`): delegate to module-level `write_fort14` at line 904. Five-line change.
2. **Fix B2** — `signed_area` quad shoelace (`CHILmesh.py:171`): replace with full 8-term cross-product expansion matching MATLAB `.m:1900-1908`. Vectorize — no Python loop, no `elem_id in tri_elems` membership.
3. **Fix B3** — `_elem_type` (`CHILmesh.py:208`): detect quads via row-pattern mask (`v[:,0]==v[:,3]` or similar sentinel) rather than `vertices[3] == 0`.
4. **Fix B4** — `_ensure_ccw_orientation` (`CHILmesh.py:119`): branch on `_elem_type` per element rather than `shape[1]`. Vectorize flip with fancy indexing.
5. **Fix B5** — `interior_angles` quad branch (`CHILmesh.py:716`): clip `cos` to `[-1, 1]` before `arccos`; use `np.real` on `arccos` output; mirror MATLAB `.m:1602` `real(acosd(...))`.
6. **Fix B6** — Add `src/chilmesh/utils/__init__.py`.

### Code changes (test fixture infrastructure)

7. **Move fort.14 fixtures** from `doc/domains/fort_14/` to `src/chilmesh/data/`; add `[tool.setuptools.package-data]`. Access via `importlib.resources.files("chilmesh.data") / "annulus_200pts.fort.14"`.
8. **Add `chilmesh.examples` module** exposing `examples.annulus()`, `examples.donut()`, `examples.block_o()`, `examples.structured()`. README quickstart becomes copy-pasteable.

### Tests (release gate — must pass before tag)

1. `test_fort14_roundtrip_identity[fixture]` — load each of 4 fixtures, write to tmp_path, reload, assert connectivity bit-identical and points within `np.allclose`. **Catches B1.** Parametrized over all 4 fixtures.
2. `test_signed_area_positive_ccw[fixture]` — for each fixture, `mesh.signed_area().min() > 0` after `_ensure_ccw_orientation`. **Catches B2 and B4** for quad-bearing fixtures.
3. `test_elem_type_classification` — synthetic quad mesh with elements including vertex 0 in 4th slot; assert all classified as quads. **Catches B3.**
4. `test_interior_angles_sum[fixture]` — for every triangle element, angles sum to π ± 1e-9; for every quad, sum to 2π ± 1e-9. **Catches B5.**
5. `test_layers_invariants[fixture]` — disjoint cover of all elements, layer-to-layer adjacency, monotone non-increasing layer size on convex domains. **Replaces** value-pinned assertions in existing `test_layers_annulus.py`.
6. `test_connectivity_wellformed[fixture]` — no duplicate vertex per element, all indices in range, no degenerate elements (zero area).
7. `test_ensure_ccw_orientation_idempotent[fixture]` — applying twice equals applying once.
8. `test_package_imports` — `import chilmesh; chilmesh.CHILmesh; chilmesh.examples.annulus(); from chilmesh.utils import plot_utils` all succeed. **Catches B6.**

Keep valid existing `test_layers_annulus.py` tests; rewrite value-pinned ones as invariants.

### Packaging

9. **Strip `requirements.txt`** to runtime deps only (`numpy>=1.23`, `scipy>=1.10`, `matplotlib>=3.6`) OR delete it and rely on `pyproject.toml`. Add `[project.optional-dependencies]` with `dev = [pytest, pytest-cov, build, twine, hypothesis]`, `examples = [jupyter]`.
10. **Add to `pyproject.toml`**:
    - `keywords = ["mesh", "adcirc", "fort.14", "hydrodynamics", "ocean", "smoothing", "quadrilateral", "finite-element", "mesh-quality", "computational-geometry"]`
    - `classifiers` including `Topic :: Scientific/Engineering :: Hydrology`, Oceanography, Mathematics, Python versions, `License :: OSI Approved :: MIT License`, `Operating System :: OS Independent`.
    - `[project.urls]`: Homepage, Documentation, Source, Issues, Changelog.
    - Bump `requires-python = ">=3.10"`.
11. **Bump version** to `0.1.1`.

### CI

12. Add `.github/workflows/ci.yml`:
    - Matrix: Python 3.10, 3.11, 3.12 × ubuntu-latest, macos-latest.
    - Steps: `pip install -e .[dev]`, `pytest -x --cov=chilmesh --cov-report=term-missing`, `python -m build`, `twine check dist/*`, install built wheel into fresh venv and `python -c "import chilmesh; chilmesh.examples.annulus()"`.
    - **Required check** on PRs to `main`.
13. Add `.github/workflows/release.yml` triggered on tag `v*`: build sdist+wheel, verify tag matches `pyproject.toml` version, run full test suite, publish via PyPI Trusted Publishing.

### Documentation

14. **Fix README example** — replace `/kaggle/working/...` with `chilmesh.examples.annulus()`. Verify end-to-end on clean venv.
15. **Fix README install snippet** — `pip install -r requirements.txt` (add `-r`).
16. **Delete stale "pip installation" To-Do** (done).
17. **Add `CHANGELOG.md`** documenting 0.1.0 → 0.1.1 with all bug fixes, `requires-python` bump, new `chilmesh.examples` API.
18. **Add "Status: alpha, API unstable pre-1.0" banner** to README.
19. **Reframe "Future Work → MADMESHR"** as "See also / Downstream projects" with explicit "CHILmesh remains agnostic to MADMESHR" sentence.

### Release engineering

20. **Yank 0.1.0 from PyPI** simultaneously with 0.1.1 publish. Yank reason: "Contains data-corrupting bugs in `signed_area` and `write_to_fort14`; upgrade to 0.1.1." (PEP 592 — pinned installs still resolve, fresh resolution skips.)

**0.1.1 sign-off criteria:**
- All 11 bugs B1-B11 fixed.
- All 8 tests green on CI matrix.
- Wheel installs cleanly on fresh venv; `chilmesh.examples.annulus()` returns working mesh.
- README example runs end-to-end without modification.
- 0.1.0 yanked from PyPI.

---

## Joint Phase 2 — 0.1.2 "Test gate + perf"

**Estimated work:** 20-30 hours over 2-4 weeks after 0.1.1 ships.

**Goal:** Make test suite real release gate with coverage enforcement; pay down O(N²) Python loops blocking scaling past ~10⁴ elements.

### Tests (expand from 8 to ~25)

21. **Geometry math characterization** (5-6 tests): equilateral triangle skew == 1.0, unit square quad area == 1.0, known-good quality on hand-built mesh, `interior_angles` column ordering documented and tested.
22. **Smoothing characterization** (3-4 tests): `angle_based_smoother` preserves topology and boundary points, `direct_smoother` does not increase max skew, idempotent on already-smooth structured mesh, no element inversions post-smoothing. Land **after** Phase 1 correctness fixes.
23. **I/O edge cases** (3-4 tests): rejects non-tri-non-quad explicitly, handles CRLF/LF line endings, refuses malformed fort.14 with clean error, write→read on smoothed mesh.
24. **Mixed-element tests** (2-3 tests): quality on mixed mesh, smoothing on mixed mesh, layer decomposition on mixed mesh.
25. **Hypothesis property layer** (3 tests): random CCW triangulations satisfy signed_area>0, Euler characteristic, edge↔elem symmetry. Strategy: 3-50 random points in unit square, scipy Delaunay → CHILmesh round-trip.
26. **Single-element edge cases** (2 tests): 1-tri mesh, 2-tri "diamond" mesh — adjacencies build, layer count == 1, all edges marked boundary.

### Performance (Optimizer Phase A)

27. **Vectorized adjacency build** — replace `_identify_edges` + `_build_elem2edge` + `_build_edge2elem` with single `np.unique`-driven pass. ~10⁴× speedup at 10⁵ elements.
28. **Vectorized `_elem_type`** — pure boolean ops, cache `tri_mask` on instance, invalidate on mutation. ~100× speedup.
29. **Vectorized `signed_area`** — vectorized shoelace using cached `tri_mask`.
30. **Fix quadratic membership in `interior_angles` and `elem_quality`** — replace `[elem_id in tri_elems for ...]` with cached mask. ~500× at 10⁵ elements.
31. **CSR adjacencies** — convert `vert2edge`, `vert2elem`, `elem2elem` from list-of-lists to CSR int32 arrays (`indptr`, `indices`). ~10× memory drop, zero-allocation neighbor queries.
32. **Cached fast paths** — `signed_area(elem_id)`, `interior_angles(elem_id)`, `elem_quality(elem_id)` should not allocate full arrays for single-element query.

### CI

33. **Add coverage gate** — `pytest --cov=chilmesh --cov-fail-under=60`. Ratchet to 80% by 0.2.0.
34. **Add Python 3.13 and Windows** to CI matrix.

### Documentation

35. **Add `examples/` directory** with three runnable scripts: `01_load_fort14.py`, `02_quality_and_smooth.py`, `03_layer_analysis.py`. CI runs them as smoke tests with `Agg` backend.
36. **Add API docstring pass** on public methods only.

---

## Joint Phase 3 — 0.2.0 "General-purpose editing API"

**Estimated work:** 30-50 hours over 1-2 months after 0.1.2 ships.

**Goal:** Add editing primitives adaptive remeshers, interactive editors, FE solvers, and downstream research projects need from 2D mesh library. **Library-neutral naming throughout. Anything serving only one consumer: cut.**

### New API (multi-consumer justified)

| API | Mirrors | Consumers served |
|---|---|---|
| `boundary_loops() -> List[np.ndarray]` | `boundary_edges()` | FE solver (BCs), editor, refiner, plotting |
| `compact() -> dict[str, np.ndarray]` returning `{vert, elem, edge}` permutations | new but explicit | refiner, editor, quad generator |
| `is_geometrically_valid(elem_ids=None) -> np.ndarray[bool]` (CCW + non-self-intersecting + positive area) | `signed_area` | refiner, editor, quad generator, RL agent |
| `vert2vert(vert_ids=None)` (one-ring vertex neighbors) | `vert2elem` | smoother, refiner, FE assembly |
| `edge_length(edge_ids=None)`, `elem_centroid(elem_ids=None)` | `edge2vert` | every consumer |
| `locate_point(xy, hint_elem=None) -> int` (containing element via walk) | MATLAB-flavored | refiner, editor, FE interpolation |
| `add_element(vertex_ids)`, `remove_element(elem_id)`, `flip_edge(edge_id)`, `split_edge(edge_id, t=0.5)` | new | refiner, editor, quad generator, RL agent |
| `elem_quality_one(elem_id)` scalar fast path | `elem_quality` | smoother, refiner, RL agent |

**Explicitly NOT added** (rejected as MADMESHR-coupling):
- `pending_loops` / `pending_elements` tracking
- `would_overlap(candidate_quad)` (RL-flavored — replaced by general `is_geometrically_valid` over `mesh.copy()`)
- `snapshot()` / `restore()` (`copy()` already exists; Pythonic pattern: `trial = mesh.copy(); trial.mutate(...); if ok: mesh = trial`)
- `commit()` / `rollback()` (heavy abstraction; `copy()` covers it)
- `valid_mask` tombstoning + stable handles (`compact()` + refresh handles instead)

### Tests for editing API

37. **Fuzz harness** — perform N random mutations, then re-run all Phase 1 invariants (`signed_area>0`, Euler, edge↔elem symmetry). Single test, high coverage.
38. **Per-primitive correctness tests** — `add_element` then `remove_element` is identity, `flip_edge` is involutive, `split_edge` increases vertex/element counts by expected amounts.
39. **Compact round-trip** — `m.copy().compact()` returns permutations such that re-indexing original yields same mesh.

### Tests as gate (mandatory)

40. **Coverage floor: 80%** of `topology + geometry + io` modules.
41. **No release tag** without all tests green AND coverage at floor.

---

## Cut list — explicitly deferred or dropped

| Item | Decision | Reason |
|---|---|---|
| Delaunay-from-scratch / zero-input `CHILmesh()` | **Cut** | `scipy.spatial.Delaunay` and `triangle` exist; no differentiation. Multi-month sink. |
| Native GMSH I/O | **Cut, replace with `meshio` adapter** | One ~50-line `meshio` adapter handles GMSH and 30 other formats. |
| `delaunayTriangulation` object I/O | **Cut** | Vanishingly small audience. |
| Full MATLAB plotting parity | **Cut** | Cosmetic. MATLAB plots not a spec. Existing matplotlib coverage sufficient. |
| Units conversion | **Cut** | Scope creep. fort.14 is unit-agnostic by convention. |
| Property-based fuzzing past Phase 2 | **Defer to 0.3** | Hypothesis is maintenance surface; revisit when test count > 50. |
| 3 OS × 5 Python CI matrix | **Defer to 0.1.2 (partial)** | ubuntu+macos × 3.10/3.11/3.12 in 0.1.1, add Windows + 3.13 in 0.1.2. |
| JOSS submission | **Defer past 1.0** | Get working+tested first. |
| Read the Docs / Sphinx site | **Defer to 0.3** | API docstring pass in 0.1.2 enough for now. |
| Mutation API | **0.2.0** (not earlier) | New surface on broken primitives compounds risk. |
| Renaming "CHILmesh" | **Keep name** | Lab provenance is citation story. Decouple from CHIL Lab in masthead instead. |
| MADMESHR-specific concepts | **Forever rejected** | CHILmesh is general-purpose library. MADMESHR builds on it, not reverse. |

---

## Open questions for maintainer

Require human decision before 0.1.1 tag.

1. **Is PyPI Trusted Publishing setup blocked?** If yes, 0.1.1 ships via manual `twine upload` and TP lands in 0.1.2.
2. **Do you want 0.1.0 yanked?** Strong audit recommendation: yes, with yank reason linking changelog. Yanking is public signal.
3. **What is actual `_mesh_layers` IE intent?** Fidelity found Python's IE diverges from MATLAB. Audit resolved to test invariants only. You wrote original — your call which is "right."
4. **Salesman's Colab notebook** — "5 minutes with CHILmesh" notebook (~4h work, lives in `notebooks/`). Slot in 0.1.2 or 0.2.0?
5. **Comparison table vs `meshio` / `pygmsh` / `oceanmesh`** — Three rows enough. Slot in 0.1.1 (README) or 0.1.2?

---

## Test plan summary

Tests are release gate per user directive "lots of tests written so that we can ensure they're working before publishing new versions."

| Phase | Test count (cumulative) | Coverage gate | New scope |
|---|---|---|---|
| 0.1.1 | **8 + existing** | None (track only) | Bug regression + invariants + parametrization over 4 fixtures + package imports |
| 0.1.2 | **~25** | 60% line, no decrease | Geometry math, smoothing characterization, I/O edge cases, mixed-element, Hypothesis properties, single-element edge cases |
| 0.2.0 | **~40** | 80% line, no decrease | Editing API fuzz harness, per-primitive correctness, compact round-trip |

**Hard rule:** every shipped bug gets regression test. Every new public API ships with at least one invariant test. CI refuses to publish if tests fail or coverage drops.

**Test fixtures:** all four `.fort.14` files in `src/chilmesh/data/`, accessible via `importlib.resources`. Tests parametrize over full set; hand-rolled fixtures only for synthetic edge cases (single triangle, two-tri diamond, vertex-0 quad).

---

## Auditor positions appendix

### Where auditors agreed

- B1-B11 real bugs, must fix before 0.1.1 tag.
- 0.1.0 should be yanked.
- Fort.14 fixtures must ship in wheel.
- CI must be required before any new release.
- `chilmesh.examples.annulus()` should be 0.1.1 deliverable.
- MADMESHR-specific concepts must not enter CHILmesh.

### Where auditors disagreed (and how resolved)

| Tension | Positions | Resolution |
|---|---|---|
| **0.1.1 scope** | Realist: 6.5h, recursion fix + paperwork, 1 test. Everyone else: 12-15h with 8 tests + correctness fixes + CI. | Bigger scope adopted (user directive: "lots of tests"). Realist's hour budgeting honored by deferring perf, plotting, mutation API. |
| **`_mesh_layers` IE ground truth** | Fidelity: port MATLAB literally. QA + Big Picture + Realist: invariants only, MATLAB not spec. | Invariant approach wins (4 vs 1). |
| **CI matrix size** | Release: 3 OS × 5 Python. QA: 2 OS × 4 Python. Realist: 1 OS × 1 Python. | Compromise: 2 OS × 3 Python in 0.1.1, expand in 0.1.2. |
| **Coverage gate timing** | Release: 80% by 0.1.2. QA: 60% in 0.1.1, 80% by 0.2.0. Realist: no gate. | QA's ratcheting: track-only in 0.1.1, 60% in 0.1.2, 80% in 0.2.0. |
| **Adjacency vectorization timing** | Optimizer + Fidelity: bundle with 0.1.1. Downstream + Realist: separate to 0.1.2. | 0.1.2 (separated). Hotfix surface stays minimal. |
| **Mutation API timing** | Downstream: 0.2.0 after correctness foundation. Fidelity: orthogonal extension. | 0.2.0 unanimous. |
| **Plotting parity** | Fidelity round 1: Phase 5. Everyone in round 2: cut/defer. | Cut from roadmap. |
| **MATLAB→Python full port** | Fidelity round 1: explicit goal. Big Picture + Realist: not goal. | Cut. Port what's needed; behavioral equivalence on documented use cases only. |

### What each auditor signed off on for joint Phase 1

All eight auditors converged on same 0.1.1 deliverable list, with nuances:

- **Fidelity:** signs off provided B1-B5 all fixed and existing layer test rewritten as invariants.
- **Release Engineer:** signs off provided 0.1.0 yanked, packaging metadata added, CI mandatory going forward. Acknowledged missing recursive write bug in round 1.
- **QA:** signs off provided all 8 tests exist and pass in CI before tagging. Tests must ship in same PR as fixes, test-first within PR.
- **Downstream (reframed):** signs off provided `boundary_loops()` ships in 0.1.1 and rest of editing API waits for 0.2.0. Drops all MADMESHR-coupled proposals from round 1.
- **Big Picture:** signs off provided README reframed around strategic identity, MADMESHR "Future Work" badge rewritten as "See also," and Delaunay-from-scratch / GMSH-from-scratch ambitions explicitly cut.
- **Optimizer:** signs off on 0.1.1 as scoped, but flags O(N²) `_build_elem2edge` / `_build_edge2elem` as single largest perf liability — must land in 0.1.2 or any consumer working on meshes >10⁴ elements hits wall.
- **Realist:** signs off reluctantly. Pushed back on test count but conceded to user directive. Insisted README example actually runs on clean venv as sign-off criterion.
- **Salesman:** signs off provided `chilmesh.examples.annulus()` ships in 0.1.1. The broken example path in 0.1.0 is worst credibility leak in repo and must die in 0.1.1.

---

## Sequencing recommendation

One weekend: do 0.1.1. Plan is sized for it.

If weekend overruns: minimum-viable-honest cut = **B1 + B2 + B3 + B6 + #7 + #8 + #14 + #15 + #19 + tests #1-#3 + #8 + minimal CI (#12 ubuntu only) + yank 0.1.0**. That is Realist's "make it honest" floor at ~6 hours. Ship that, then come back for rest of 0.1.1 the following week.

Two things cannot ship without: **(a)** recursive `write_to_fort14` bug fixed with regression test, and **(b)** README example actually running on clean venv. Everything else negotiable.

---

*This report is output of an 8-agent audit run on 2026-04-07. Treat as planning document, not binding commitment. Adapt as new information surfaces. Update or supersede when 0.1.1 ships.*
