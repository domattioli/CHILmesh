# CHILmesh Strategic Audit & Phased Plan

**Date:** 2026-04-07
**Method:** 8 independent auditor agents (Fidelity, Release Engineer, QA, Downstream Integration, Big Picture, Optimizer, Realist, Salesman) ran two rounds of independent analysis followed by mediated convergence.
**Output:** This document.

---

## Executive summary

CHILmesh 0.1.0 is on PyPI but ships **three latent bugs that corrupt user data on first use**: a recursive `write_to_fort14` (infinite loop on first call), a mathematically wrong quad shoelace in `signed_area` (corrupts CCW orientation for any quad mesh), and a `_elem_type` misclassification of any quad whose 4th vertex is index 0. The README also contains a broken install command and an example path pointing at `/kaggle/working/...` that no pip user can run.

The audit converged on a **three-release plan** with hard test gates:

- **0.1.1 — "Honest hotfix"** (~12-15 hours of work): fix the bugs, ship test fixtures via `importlib.resources`, add 8 tests covering the bugs and core invariants, set up minimal CI, fix README. Yank 0.1.0.
- **0.1.2 — "Test gate + perf"** (~20-30 hours): expand to ~25 tests including a Hypothesis property layer, vectorize the `np.unique`-class adjacency hotspots, add coverage gating to CI.
- **0.2.0 — "General-purpose editing API"** (~30-50 hours): add `boundary_loops()`, `compact()`, `is_geometrically_valid()`, `vert2vert()`, `elem_centroid()`, `edge_length()`, `locate_point()` — all multi-consumer, library-neutral.

**Strategic identity (Big Picture):** CHILmesh is a research-grade Python implementation of the QuADMESH layer abstraction, maintained as a substrate for downstream research, with first-class ADCIRC `.fort.14` I/O. It is *not* a `gmsh` replacement, *not* a 3D library, and *not* a general mesh generator. Delaunay-from-scratch and native GMSH I/O are cut from the roadmap; a thin `meshio` adapter replaces both.

**Test policy (per user directive):** No release tagged without all tests green and the bug-regression suite covering every prior shipped bug.

---

## Convergent findings — bugs in 0.1.0

These were independently surfaced by multiple auditors and have no debate.

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
| B9 | `README.md:97` | Example path `/kaggle/working/CHILmesh/...` is broken for any pip user | Release, Realist, Salesman | Medium |
| B10 | `pyproject.toml` | Missing `keywords`, `classifiers`, `[project.urls]` → PyPI page is barren and unsearchable | Release, Salesman | Medium |
| B11 | `tests/` | Only 1 test file, only 1 of 4 available `.fort.14` fixtures used; 0 tests for write, quality, smoothing, mixed-element, I/O round-trip | QA, Realist | High |

**Note on B1:** The recursive `write_to_fort14` has never been successfully executed by any test. QA's proposed `test_fort14_roundtrip_identity` would have caught it on day one — and is the single highest-ROI test in the entire plan.

---

## Strategic identity (Big Picture)

The Big Picture auditor surfaced an unstated tension: CHILmesh is currently performing three identities simultaneously and committing to none — (a) faithful Python port of a 2017 thesis, (b) general-purpose 2D mesh library competing with `meshio`/`pygmsh`/`gmsh`, (c) substrate for the author's downstream research line (MADMESHR, hydrodynamic ML). These imply incompatible priorities.

**Decision (adopted across the audit):** identity is **(a) preserved + (c) actively maintained**. CHILmesh is a research-grade Python implementation of the QuADMESH layer abstraction, packaged cleanly enough to serve downstream research and ADCIRC users. It will *not* compete with `gmsh` on Delaunay generation or with `meshio` on format breadth.

**Implications:**

1. **Cut from roadmap:** Delaunay-from-scratch (`scipy.spatial.Delaunay` exists), native GMSH I/O (delegate to `meshio`), full plotting parity with MATLAB (the MATLAB plots are not a spec), units conversion (scope creep), `delaunayTriangulation` object I/O (vanishingly small audience).
2. **Stay in scope:** layer-based decomposition (the unique selling point), tri/quad/mixed first-class support, ADCIRC `.fort.14` I/O, FEM + angle-based smoothing, mesh quality/skewness metrics, a clean editing API (general-purpose, not RL-flavored).
3. **README reframing:** the "Future Work → MADMESHR" badge inverts the dependency arrow (MADMESHR depends on CHILmesh, not the reverse). Reframe as "See also / Downstream projects." Add a non-coupling statement so future contributors don't try to merge concerns.
4. **No MADMESHR-specific concepts** in CHILmesh APIs. Library-neutral naming. Anything that ONLY serves an RL agent is out.

**3-year north star:** "MADMESHR (and 1-2 other downstream research projects) cite and depend on `chilmesh>=1.0` from PyPI; the layer-based mesh abstraction from the thesis is available to anyone who reads the QuADMESH paper; the API has not broken in 18 months."

---

## Joint Phase 1 — 0.1.1 "Honest hotfix"

**Estimated work:** 12-15 hours of focused effort. Realistically one weekend.

**Goal:** Ship a release where `pip install chilmesh` followed by the README quickstart actually works on a clean machine, with regression tests for every shipped bug, and CI gates future merges.

### Code changes (correctness)

1. **Fix B1** — `write_to_fort14` (`CHILmesh.py:646`): delegate to module-level `write_fort14` at line 904. Five-line change.
2. **Fix B2** — `signed_area` quad shoelace (`CHILmesh.py:171`): replace with full 8-term cross-product expansion matching MATLAB `.m:1900-1908`. Vectorize while you're there (Optimizer A4) — no Python loop, no `elem_id in tri_elems` membership.
3. **Fix B3** — `_elem_type` (`CHILmesh.py:208`): detect quads via row-pattern mask (`v[:,0]==v[:,3]` or similar sentinel) rather than `vertices[3] == 0`.
4. **Fix B4** — `_ensure_ccw_orientation` (`CHILmesh.py:119`): branch on `_elem_type` per element rather than `shape[1]`. Vectorize the flip with fancy indexing (Optimizer B6).
5. **Fix B5** — `interior_angles` quad branch (`CHILmesh.py:716`): clip `cos` to `[-1, 1]` before `arccos`; use `np.real` on `arccos` output if needed; mirror MATLAB's `.m:1602` `real(acosd(...))`.
6. **Fix B6** — Add `src/chilmesh/utils/__init__.py` (re-exporting `plot_utils` symbols if appropriate, or empty).

### Code changes (test fixture infrastructure)

7. **Move fort.14 fixtures** from `doc/domains/fort_14/` to `src/chilmesh/data/` and add `[tool.setuptools.package-data]` so they ship in the wheel. Access via `importlib.resources.files("chilmesh.data") / "annulus_200pts.fort.14"`.
8. **Add `chilmesh.examples` module** exposing `examples.annulus()`, `examples.donut()`, `examples.block_o()`, `examples.structured()` returning `CHILmesh` instances. This is a one-day fix with massive UX payoff (Salesman's #1 priority): the README quickstart becomes copy-pasteable.

### Tests (must exist, must pass before tag)

These are the **release gate**. CI refuses to publish if any of these fail.

1. `test_fort14_roundtrip_identity[fixture]` — load each of 4 fixtures, write to tmp_path, reload, assert connectivity bit-identical and points within `np.allclose`. **Catches B1.** Parametrized over all 4 fixtures.
2. `test_signed_area_positive_ccw[fixture]` — for each fixture, `mesh.signed_area().min() > 0` after `_ensure_ccw_orientation`. **Catches B2 and B4** for quad-bearing fixtures.
3. `test_elem_type_classification` — synthetic quad mesh whose elements include vertex 0 in the 4th slot; assert all classified as quads. **Catches B3.**
4. `test_interior_angles_sum[fixture]` — for every triangle element, angles sum to π ± 1e-9; for every quad element, angles sum to 2π ± 1e-9. **Catches B5.**
5. `test_layers_invariants[fixture]` — disjoint cover of all elements, layer-to-layer adjacency, monotone non-increasing layer size on convex domains. **Replaces** the existing `test_layers_annulus.py` value-pinned assertions; does not commit either MATLAB or Python's exact IE definition as ground truth. (QA's resolution, ratified by Big Picture and Realist; Fidelity's "port MATLAB literally" position was overruled.)
6. `test_connectivity_wellformed[fixture]` — no duplicate vertex per element, all indices in range, no degenerate elements (zero area).
7. `test_ensure_ccw_orientation_idempotent[fixture]` — applying twice equals applying once.
8. `test_package_imports` — `import chilmesh; chilmesh.CHILmesh; chilmesh.examples.annulus(); from chilmesh.utils import plot_utils` all succeed. **Catches B6.**

Plus: keep the existing `test_layers_annulus.py` tests that are still valid; rewrite the value-pinned ones as invariants.

### Packaging

9. **Strip `requirements.txt`** to runtime deps only (`numpy>=1.23`, `scipy>=1.10`, `matplotlib>=3.6`) OR delete it entirely and rely on `pyproject.toml`. Add `[project.optional-dependencies]` with `dev = [pytest, pytest-cov, build, twine, hypothesis]`, `examples = [jupyter]`.
10. **Add to `pyproject.toml`**:
    - `keywords = ["mesh", "adcirc", "fort.14", "hydrodynamics", "ocean", "smoothing", "quadrilateral", "finite-element", "mesh-quality", "computational-geometry"]`
    - `classifiers` including `Topic :: Scientific/Engineering :: Hydrology`, `Topic :: Scientific/Engineering :: Oceanography`, `Topic :: Scientific/Engineering :: Mathematics`, Python version classifiers, `License :: OSI Approved :: MIT License`, `Operating System :: OS Independent`.
    - `[project.urls]`: Homepage, Documentation, Source, Issues, Changelog.
    - Bump `requires-python = ">=3.10"` (drop EOL 3.8 and 3.9; honest about what's tested).
11. **Bump version** to `0.1.1`.

### CI

12. Add `.github/workflows/ci.yml`:
    - Matrix: Python 3.10, 3.11, 3.12 × ubuntu-latest, macos-latest. (Drop Windows and 3.13 to 0.1.2 — Realist correctly notes Windows fort.14 line endings will eat hours.)
    - Steps: `pip install -e .[dev]`, `pytest -x --cov=chilmesh --cov-report=term-missing`, `python -m build`, `twine check dist/*`, install the built wheel into a fresh venv and `python -c "import chilmesh; chilmesh.examples.annulus()"`.
    - **Required check** on PRs to `main` going forward. No coverage gate yet (added in 0.1.2).
13. Add `.github/workflows/release.yml` triggered on tag `v*`: build sdist+wheel, verify tag matches `pyproject.toml` version, run full test suite, publish via PyPI Trusted Publishing (configure once on PyPI; manual `twine upload` is acceptable fallback if OIDC setup blocks the hotfix).

### Documentation

14. **Fix README example** — replace `/kaggle/working/...` with `chilmesh.examples.annulus()` (now possible via #8). Verify end-to-end on a clean venv.
15. **Fix README install snippet** — `pip install -r requirements.txt` (add the `-r`).
16. **Delete the stale "pip installation" To-Do** (it's done).
17. **Add `CHANGELOG.md`** documenting 0.1.0 → 0.1.1 with all bug fixes, the `requires-python` bump as breaking, and the new `chilmesh.examples` API.
18. **Add a "Status: alpha, API unstable pre-1.0" banner** to the README.
19. **Reframe "Future Work → MADMESHR"** as "See also / Downstream projects" with an explicit "CHILmesh remains agnostic to MADMESHR" sentence.

### Release engineering

20. **Yank 0.1.0 from PyPI** (`pip yank chilmesh==0.1.0`) simultaneously with 0.1.1 publish. Yank reason: "Contains data-corrupting bugs in `signed_area` and `write_to_fort14`; upgrade to 0.1.1." Yank does not delete — pinned installs still resolve, fresh resolution skips. PEP 592.

**0.1.1 sign-off criteria:**

- All 11 bugs B1-B11 fixed.
- All 8 tests above green on CI matrix.
- Wheel installs cleanly on a fresh venv; `chilmesh.examples.annulus()` returns a working mesh.
- README example runs end-to-end without modification.
- 0.1.0 yanked from PyPI.

---

## Joint Phase 2 — 0.1.2 "Test gate + perf"

**Estimated work:** 20-30 hours over 2-4 weeks after 0.1.1 ships.

**Goal:** Make the test suite a real release gate with coverage enforcement, and pay down the O(N²) Python loops that block scaling past ~10⁴ elements.

### Tests (expand from 8 to ~25)

21. **Geometry math characterization** (5-6 tests): equilateral triangle skew == 1.0, unit square quad area == 1.0, known-good quality on a hand-built mesh, `interior_angles` column ordering documented and tested.
22. **Smoothing characterization** (3-4 tests): `angle_based_smoother` preserves topology and boundary points, `direct_smoother` does not increase max skew, idempotent on already-smooth structured mesh, no element inversions post-smoothing. Land **after** the Phase 1 correctness fixes — characterizing broken behavior is bad.
23. **I/O edge cases** (3-4 tests): rejects non-tri-non-quad explicitly, handles CRLF/LF line endings, refuses malformed fort.14 with a clean error, write→read on a smoothed mesh.
24. **Mixed-element tests** (2-3 tests): quality on mixed mesh, smoothing on mixed mesh, layer decomposition on mixed mesh.
25. **Hypothesis property layer** (3 tests): random CCW triangulations satisfy signed_area>0, Euler characteristic, edge↔elem symmetry. Strategy: 3-50 random points in unit square, scipy Delaunay → CHILmesh round-trip.
26. **Single-element edge cases** (2 tests): 1-tri mesh, 2-tri "diamond" mesh — adjacencies build, layer count == 1, all edges marked boundary.

### Performance (Optimizer Phase A)

27. **Vectorized adjacency build** — replace `_identify_edges` + `_build_elem2edge` + `_build_edge2elem` with a single `np.unique`-driven pass. ~10⁴× speedup at 10⁵ elements. Optimizer A1+A2.
28. **Vectorized `_elem_type`** — pure boolean ops, cache `tri_mask` on the instance, invalidate on mutation. ~100× speedup. Optimizer A3.
29. **Vectorized `signed_area`** — vectorized shoelace using cached `tri_mask`. Optimizer A4.
30. **Fix quadratic membership in `interior_angles` and `elem_quality`** — replace `[elem_id in tri_elems for ...]` with the cached mask. ~500× at 10⁵ elements. Optimizer A5.
31. **CSR adjacencies** — convert `vert2edge`, `vert2elem`, `elem2elem` from list-of-lists to CSR int32 arrays (`indptr`, `indices`). ~10× memory drop, zero-allocation neighbor queries. Optimizer B1+B2.
32. **Cached fast paths** — `signed_area(elem_id)`, `interior_angles(elem_id)`, `elem_quality(elem_id)` should not allocate full arrays for a single-element query. Optimizer B4.

### CI

33. **Add coverage gate** — `pytest --cov=chilmesh --cov-fail-under=60`. Ratchet to 80% by 0.2.0.
34. **Add Python 3.13 and Windows** to the CI matrix now that the cross-platform fort.14 line-ending edge cases are tested in #23.

### Documentation

35. **Add `examples/` directory** with three runnable scripts: `01_load_fort14.py`, `02_quality_and_smooth.py`, `03_layer_analysis.py`. CI runs them as smoke tests with the `Agg` matplotlib backend.
36. **Add API docstring pass** on the public methods only.

---

## Joint Phase 3 — 0.2.0 "General-purpose editing API"

**Estimated work:** 30-50 hours over 1-2 months after 0.1.2 ships.

**Goal:** Add the editing primitives that adaptive remeshers, interactive editors, FE solvers, and downstream research projects (including but not limited to MADMESHR) need from a 2D mesh library. **Library-neutral naming throughout. Anything that only serves one consumer is cut.**

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

**Explicitly NOT added** (rejected from the Downstream auditor's round-1 plan as MADMESHR-coupling):

- `pending_loops` / `pending_elements` tracking
- `would_overlap(candidate_quad)` (RL-flavored — replaced by general `is_geometrically_valid` over `mesh.copy()`)
- `snapshot()` / `restore()` (`copy()` already exists at line 899; the Pythonic trial-move pattern is `trial = mesh.copy(); trial.mutate(...); if ok: mesh = trial`)
- `commit()` / `rollback()` (heavy abstraction; `copy()` covers it)
- `valid_mask` tombstoning + stable handles (philosophy change touching every accessor; consumers can call `compact()` and refresh handles instead)

### Tests for editing API

37. **Fuzz harness** — perform N random mutations from the editing primitive set, then re-run all Phase 1 invariants (`signed_area>0`, Euler, edge↔elem symmetry, layer cover). Single test, high coverage. QA's recommendation.
38. **Per-primitive correctness tests** — `add_element` then `remove_element` is identity, `flip_edge` is involutive, `split_edge` increases vertex/element counts by expected amounts.
39. **Compact round-trip** — `m.copy().compact()` returns permutations such that re-indexing the original yields the same mesh.

### Tests as gate (now mandatory)

40. **Coverage floor: 80%** of `topology + geometry + io` modules. Smoothing and plotting can drag the project denominator; gate on subdirectories.
41. **No release tag** without all tests green AND coverage at floor.

---

## Cut list — explicitly deferred or dropped

These are out of scope for 0.1.x and 0.2.x. Most are deferred indefinitely.

| Item | Decision | Reason |
|---|---|---|
| Delaunay-from-scratch / zero-input `CHILmesh()` | **Cut** | `scipy.spatial.Delaunay` and `triangle` exist; no differentiation. Multi-month sink. |
| Native GMSH I/O | **Cut, replace with `meshio` adapter** | Don't reimplement format parsers; one ~50-line `meshio` adapter handles GMSH and 30 other formats. |
| `delaunayTriangulation` object I/O | **Cut** | Vanishingly small audience (MATLAB refugees). |
| Full MATLAB plotting parity (`axisCHILmesh`, `plotEdge`, `plotElem`, `plotLabel`, `plotPoint`, `plotQuality`) | **Cut** | Cosmetic. The MATLAB plots are not a spec. Existing matplotlib coverage is sufficient for headless use. |
| Units conversion (`convert2LL/SI/USC`, `reverseZdir`) | **Cut** | Scope creep. fort.14 is unit-agnostic by convention. |
| Property-based fuzzing past Phase 2 scope | **Defer to 0.3** | Hypothesis is a maintenance surface; revisit when test count > 50. |
| 3 OS × 5 Python CI matrix | **Defer to 0.1.2 (partial)** | Realist correctly notes cross-platform burns hours. ubuntu+macos × 3.10/3.11/3.12 in 0.1.1, add Windows + 3.13 in 0.1.2. |
| JOSS submission | **Defer past 1.0** | JOSS bar is "software that works and is tested." Get there first. |
| Read the Docs / Sphinx site | **Defer to 0.3** | API docstring pass in 0.1.2 is enough for now. |
| Mutation API | **0.2.0** (not earlier) | New surface area on top of broken primitives compounds risk. |
| Renaming "CHILmesh" | **Keep name** | Lab provenance is the citation story; cost-benefit is real. Decouple from CHIL Lab in masthead instead. |
| MADMESHR-specific concepts | **Forever rejected** | CHILmesh is a general-purpose library. MADMESHR builds on it, not the reverse. |

---

## Open questions for the maintainer

These need a human decision before 0.1.1 is tagged. The audit cannot resolve them without you.

1. **Is the PyPI Trusted Publishing setup blocked?** If yes, 0.1.1 ships via manual `twine upload` and TP lands in 0.1.2. If no, configure TP first (15 min on PyPI web UI) and use it for 0.1.1.
2. **Do you want 0.1.0 yanked?** Strong audit recommendation: yes, with a yank reason linking the changelog. But yanking is a public signal and you may prefer to leave 0.1.0 visible with a "deprecated" note in the README.
3. **What is the actual `_mesh_layers` IE intent?** Fidelity found Python's IE definition diverges from MATLAB. The audit resolved to test invariants only (neither pinned), but the *implementation* still has to pick one. Recommendation: pick whichever produces a disjoint cover with monotone-shrinking layer sizes on convex domains, document the choice in a comment, and move on. You wrote the original — your call which is "right."
4. **Salesman's Colab notebook** — Salesman calls a "5 minutes with CHILmesh" Colab notebook the single highest-leverage adoption artifact. It's ~4 hours of work and lives in `notebooks/`. Slot in 0.1.2 or 0.2.0?
5. **Comparison table vs `meshio` / `pygmsh` / `oceanmesh`** — Salesman's #4. Three rows is enough. Slot in 0.1.1 (README) or 0.1.2?

---

## Test plan summary

Per the user directive "lots of tests written so that we can ensure they're working before publishing new versions," tests are the release gate. The convergence:

| Phase | Test count (cumulative) | Coverage gate | New scope |
|---|---|---|---|
| 0.1.1 | **8 + existing** | None (track only) | Bug regression + invariants + parametrization over 4 fixtures + package imports |
| 0.1.2 | **~25** | 60% line, no decrease | Geometry math, smoothing characterization, I/O edge cases, mixed-element, Hypothesis properties, single-element edge cases |
| 0.2.0 | **~40** | 80% line, no decrease | Editing API fuzz harness, per-primitive correctness, compact round-trip |

**Hard rule:** every shipped bug gets a regression test. Every new public API ships with at least one invariant test. CI refuses to publish if tests fail or coverage drops.

**Test fixtures:** all four `.fort.14` files relocated to `src/chilmesh/data/`, accessible via `importlib.resources`. Tests parametrize over the full set; no hand-rolled fixtures except for synthetic edge cases (single triangle, two-tri diamond, vertex-0 quad).

---

## Auditor positions appendix

### Where the auditors agreed

- B1-B11 are real bugs, must be fixed before 0.1.1 tag.
- 0.1.0 should be yanked.
- Fort.14 fixtures must ship in the wheel.
- CI must be required before any new release.
- `chilmesh.examples.annulus()` should be a 0.1.1 deliverable (UX + tests share the same packaging fix).
- MADMESHR-specific concepts must not enter CHILmesh.

### Where the auditors disagreed (and how it was resolved)

| Tension | Positions | Resolution |
|---|---|---|
| **0.1.1 scope** | Realist: 6.5h, recursion fix + paperwork, 1 test. Everyone else: 12-15h with 8 tests + correctness fixes + CI. | Bigger scope adopted (user directive: "lots of tests before publishing"). Realist's hour budgeting honored by deferring perf, plotting, mutation API. |
| **`_mesh_layers` IE ground truth** | Fidelity: port MATLAB literally. QA + Big Picture + Realist: invariants only, MATLAB is not the spec. | Invariant approach wins (4 vs 1; aligned with strategic identity decision). |
| **CI matrix size** | Release: 3 OS × 5 Python. QA: 2 OS × 4 Python. Realist: 1 OS × 1 Python. | Compromise: 2 OS × 3 Python in 0.1.1, expand in 0.1.2. |
| **Coverage gate timing** | Release: 80% by 0.1.2. QA: 60% in 0.1.1, 80% by 0.2.0. Realist: no coverage gate. | QA's ratcheting approach: track-only in 0.1.1, 60% in 0.1.2, 80% in 0.2.0. |
| **Adjacency vectorization timing** | Optimizer + Fidelity: bundle with 0.1.1. Downstream + Realist: separate to 0.1.2 to keep hotfix surface small. | 0.1.2 (separated). Hotfix surface stays minimal; perf gets its own release with its own tests. |
| **Mutation API timing** | Downstream: 0.2.0 after correctness foundation. Fidelity: orthogonal extension, no parity conflict. | 0.2.0 unanimous. |
| **Plotting parity** | Fidelity round 1: Phase 5. Everyone in round 2: cut/defer. | Cut from roadmap. |
| **MATLAB→Python full port** | Fidelity round 1: explicit goal. Big Picture + Realist: not the goal. | Cut. Port what's needed; behavioral equivalence on documented use cases only. |

### What each auditor "signed off" on for joint Phase 1

All eight auditors converged on essentially the same 0.1.1 deliverable list, with the following nuances:

- **Fidelity:** signs off provided B1-B5 are all fixed and the existing layer test is rewritten as invariants.
- **Release Engineer:** signs off provided 0.1.0 is yanked, packaging metadata is added, and CI is mandatory going forward. Acknowledged missing the recursive write bug in round 1 — lesson: a packaging-only hotfix would have been malpractice.
- **QA:** signs off provided all 8 tests above exist and pass in CI before tagging. Tests must ship in the same PR as the fixes, test-first within the PR.
- **Downstream (reframed):** signs off provided `boundary_loops()` ships in 0.1.1 (small scope, broadly useful) and the rest of the editing API waits for 0.2.0. Drops all MADMESHR-coupled proposals from round 1.
- **Big Picture:** signs off provided the README is reframed around the strategic identity (research-grade port + downstream substrate, not general mesh generator), the MADMESHR "Future Work" badge is rewritten as "See also," and the Delaunay-from-scratch / GMSH-from-scratch ambitions are explicitly cut.
- **Optimizer:** signs off on 0.1.1 as scoped, but flags that the O(N²) `_build_elem2edge` / `_build_edge2elem` pattern is the single largest perf liability and must land in 0.1.2 or any consumer working on meshes >10⁴ elements will hit a wall.
- **Realist:** signs off reluctantly. Pushed back on the test count but conceded to the user directive. Insisted on the README example actually running on a clean venv as a sign-off criterion, which is now item #14 above.
- **Salesman:** signs off provided `chilmesh.examples.annulus()` ships in 0.1.1 (the single change that converts the README from "broken demo" to "copy-pasteable"). The full README rewrite, badges, and Colab notebook can wait for 0.1.2/0.2.0, but the broken example path in 0.1.0 is the worst credibility leak in the repo and must die in 0.1.1.

---

## Sequencing recommendation

If you have one weekend: do 0.1.1. The plan is sized for it.

If the weekend overruns: the minimum-viable-honest cut is **B1 + B2 + B3 + B6 + #7 + #8 + #14 + #15 + #19 + tests #1-#3 + #8 + minimal CI (#12 ubuntu only) + yank 0.1.0**. That is the Realist's "make it honest" floor and it is ~6 hours. Ship that, then come back for the rest of 0.1.1 the following week.

The two things you cannot ship without are: **(a)** the recursive `write_to_fort14` bug fixed with a regression test, and **(b)** the README example actually running on a clean venv. Everything else is negotiable.

---

*This report is the output of an 8-agent audit run on 2026-04-07. It should be treated as a planning document, not a binding commitment. Adapt as new information surfaces. Update or supersede when 0.1.1 ships.*
