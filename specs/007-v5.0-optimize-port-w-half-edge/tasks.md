---
description: "Tasks for half-edge data structure investigation & language optimization path"
---

# Tasks: Half-Edge Data Structure Investigation & Language Optimization Path

**Input**: Design documents from `specs/007-v5.0-optimize-port-w-half-edge/`
**Prerequisites**: spec.md ✓, clarify.md ✓, plan.md ✓ (research.md / data-model.md deferred — handled inline in T1/T2)

**Tests**: Required. Spec FR-004 mandates all 439 existing tests pass on both backends + new equivalence/basic tests.

**Organization**: Tasks grouped by user story (US1, US2, US3 from spec.md). US1 is MVP — ships independently.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (US1, US2, US3, or SETUP for shared infra)

## Path Conventions

Single-project library layout: `src/chilmesh/`, `tests/`, `scripts/`, `docs/`, `specs/007-v5.0-optimize-port-w-half-edge/`.

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Verify baseline + pin external mesh before any backend work begins.

- [ ] **T001** [SETUP] Verify `pytest -v` passes on `main` baseline (all 439 tests green). Record runtime + machine info. Output stored as `output/baseline_pytest.log` for the PR record.
- [ ] **T002** [SETUP] Pin WNAT_Hagen via SHA-256 per clarify C-3. Compute hash of `/tmp/admesh-domains/registry_data/meshes/WNAT_Hagen.14`; write to `output/benchmark.json` field `mesh_sha256`. Write helper script `scripts/verify_mesh_pin.py` that exits non-zero on mismatch.
- [ ] **T003** [P] [SETUP] Add pytest fixture to `tests/conftest.py` that monkeypatches `CHILMESH_TOPOLOGY_BACKEND` and resets on teardown per clarify F-2. Fixture name: `topology_backend_env` (parametrized over `[None, "halfedge"]`).

**Checkpoint**: Baseline locked, mesh pinned, env-var leakage prevented. Foundation ready.

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: DCEL module + factory wiring. MUST complete before any user story can be benchmarked.

- [ ] **T004** [SETUP] Create `src/chilmesh/mesh_topology_halfedge.py`. Implements:
  - `HalfEdgeTopology` class with internal `half_edges: ndarray[n_halfedges, 4]` (cols = origin / twin_idx / next_idx / face_idx per clarify C-1).
  - `build_halfedge_from_connectivity(elem2vert: ndarray, n_verts: int) -> HalfEdgeTopology` factory.
  - Runs `_ensure_ccw_orientation` (existing function) BEFORE assigning `next_idx` per clarify F-1.
  - Reuses `_elem_type` mask from B3 to skip padded-triangle 4th slot per clarify F-4.
  - Boundary half-edges set `twin_idx == -1` (matches `Edge2Elem == -1` sentinel) per clarify C-1.
  - Conversion methods: `to_edge2vert() -> ndarray`, `to_edge2elem() -> ndarray`, `to_elem2edge() -> ndarray`, `to_vert2edge() -> dict[int, set[int]]`, `to_vert2elem() -> dict[int, set[int]]`, `to_edgemap_list() -> list[tuple[int, int]]` (canonical sorted form per clarify C-4).
  - Module docstring documents the schema + sentinel convention.
  - Estimated ~300–350 LOC.

- [ ] **T005** [SETUP] Modify `src/chilmesh/CHILmesh.py` `_build_adjacencies()` to dispatch on `topology_backend`:
  - Read kwarg first; fall back to `os.environ.get("CHILMESH_TOPOLOGY_BACKEND")`; final fallback `"edgemap"`.
  - `ValueError` on unknown value per spec FR-008.
  - For `"edgemap"`: existing code path unchanged.
  - For `"halfedge"`: call `build_halfedge_from_connectivity()`, then call its conversion methods to fill the same adjacency dict the rest of the codebase consumes.
  - Thread the kwarg through `__init__`, `read_from_fort14`, `read_from_2dm`, `from_admesh_domain` (default `topology_backend=None`).
  - `examples/*` factories accept `**kwargs` passthrough.
  - Estimated ~50–80 LOC modified across 1–2 files.

**Checkpoint**: Backend is dispatchable; the rest of the codebase consumes a uniform adjacency dict regardless of backend.

---

## Phase 3: User Story 1 — Maintainer evaluates half-edge for adoption (Priority: P1) 🎯 MVP

**Goal**: Half-edge backend produces bit-identical adjacency outputs on all 4 fixtures and passes all 439 existing tests with zero modifications.

**Independent Test**: `CHILMESH_TOPOLOGY_BACKEND=halfedge pytest -v` returns 0 with 439 passed. `python -c "from chilmesh.examples import annulus; m1 = annulus(); m2 = annulus(topology_backend='halfedge'); assert m1.adjacencies['Edge2Vert'].tolist() == m2.adjacencies['Edge2Vert'].tolist()"` succeeds (or the equivalence test in T007 covers this more rigorously).

### Tests for User Story 1 (write first, ensure they fail before T004/T005 implementation)

- [ ] **T006** [P] [US1] Write `tests/test_halfedge_basic.py`. ~150 LOC, ~15 tests:
  - HalfEdgeTopology construction succeeds on each of the 4 fixtures (annulus, donut, block_o, structured).
  - Every half-edge has a valid `next_idx` (not `-1` for interior; cyclic walk around a face returns to origin).
  - Every interior half-edge has a valid `twin_idx`; every boundary half-edge has `twin_idx == -1`.
  - Twin-of-twin is identity: `half_edges[half_edges[i, 1], 1] == i` for all interior `i`.
  - Padded triangle does NOT generate a 4th half-edge (regression test for F-4; uses block_o which has mixed elements).
  - Degenerate quad fixture (synthetic, copied from `test_interior_angles.py::test_quad_angles_no_nan_on_degenerate_quad`) does not crash construction (regression for F-1).
  - Unknown `CHILMESH_TOPOLOGY_BACKEND` value raises `ValueError` (FR-008).
  - Uses the `topology_backend_env` fixture from T003 for clean env-var handling.

- [ ] **T007** [P] [US1] Write `tests/test_halfedge_equivalence.py`. ~80 LOC, 4 fixture-parametrized tests:
  - For each fixture: build EdgeMap mesh, build HalfEdge mesh, run the canonical-form comparator (per clarify C-4) on every adjacency output. All comparators return True.
  - The comparator itself is tested with synthetic permuted inputs (1 self-check test) so we don't ship a false-green comparator.
  - Uses the `topology_backend_env` fixture from T003.

### Implementation for User Story 1

- [ ] **T008** [US1] Implement T004 (`mesh_topology_halfedge.py`) and T005 (factory wiring in `CHILmesh.py`). Run T006 + T007 — must pass.

- [ ] **T009** [US1] Run full existing pytest suite with EdgeMap backend (default). Confirm 439 tests pass; record runtime as `output/edgemap_pytest.log`. **Gate**: zero failures or US1 does not ship.

- [ ] **T010** [US1] Run full existing pytest suite with `CHILMESH_TOPOLOGY_BACKEND=halfedge`. Confirm 439 tests pass; record runtime as `output/halfedge_pytest.log`. **Gate**: zero failures or US1 does not ship.

- [ ] **T011** [US1] Diff `output/edgemap_pytest.log` vs `output/halfedge_pytest.log`. Test outcomes must match line-by-line on the `PASSED`/`FAILED` markers (runtimes differ — that's expected). Commit both logs to the branch.

**Checkpoint**: At this point, US1 ships. Half-edge backend is correctness-equivalent. Decision to adopt depends on benchmark data (US3); US1 alone is shippable as "experimental backend available, parity verified".

---

## Phase 4: User Story 2 — Downstream consumer reads updated benchmark report (Priority: P2)

**Goal**: `docs/BENCHMARK.md` is regenerated from live data; downstream readers see truthful current numbers and clearly-labeled experimental half-edge results.

**Independent Test**: Open `docs/BENCHMARK.md`. Every number under "Current Performance (v0.4.1)" matches the latest line in `output/benchmark.json` to within rounding. Half-edge results live in a section labeled "Experimental — not part of the public API".

### Implementation for User Story 2

- [ ] **T012** [US2] Extend `scripts/benchmark_wnat_hagen.py` (or wrap it) so that `output/benchmark.json` gains the additive fields per clarify Q-2:
  - `mesh_sha256` (from T002)
  - `tracemalloc_peak_bytes` per operation (per clarify F-3; uses `tracemalloc.start()` → `take_snapshot()`)
  - `backend_variants` (list of {backend, fast_init_s, full_init_s, quality_s, ...})
  - Existing fields unchanged so v0.4.1 consumers still parse it.
  - Estimated ~80–120 LOC.

- [ ] **T013** [P] [US2] Update `docs/BENCHMARK.md`:
  - "Current Performance (v0.4.1)" section: numbers come from the latest `output/benchmark.json` (no hand-edited estimates).
  - "Half-Edge Variant (Experimental)" subsection: clearly labeled "EXPERIMENTAL — not part of the public API. Pin against the EdgeMap backend until v5.1 decision."
  - Reference `output/benchmark.json` as the source of truth.
  - FR-005 gate: live regen, not estimates.

- [ ] **T014** [P] [US2] Update README.md "Performance Characteristics" section with a 2-paragraph summary linking to `docs/BENCHMARK.md`. Mention v5.0 experimental backend exists and how to opt in (env var + kwarg).

**Checkpoint**: At this point, US2 ships independently of US3. A downstream consumer can now safely pin against EdgeMap and read the experimental half-edge numbers for informational purposes.

---

## Phase 5: User Story 3 — Maintainer decides on language porting (Priority: P3)

**Goal**: 3-variant benchmark table makes the adopt / archive / investigate-port decision mechanical (no judgment call).

**Independent Test**: `python scripts/benchmark_halfedge_variants.py` runs in <60 s on the reference container, emits a markdown table to stdout with three rows (current, v1, v2) and four columns (fast init, full init, quality, query latency) with percent deltas. The "Decision Trigger" column reads "adopt", "archive", or "investigate port" per spec SC-005.

### Implementation for User Story 3

- [ ] **T015** [US3] Write `scripts/benchmark_halfedge_variants.py`. ~200 LOC:
  - Imports CHILmesh, runs all three variants on WNAT_Hagen.
  - Variant 1 = current EdgeMap (baseline).
  - Variant 2 = half-edge v1 (Python `for`-loop walk).
  - Variant 3 = half-edge v2 with NumPy vectorized walk per clarify C-2 (`np.take(half_edges[:,2], face_starts)`).
  - Median-of-3 protocol per clarify F-3 to bound NFR-003 variance.
  - Outputs: stdout markdown table + `output/benchmark.json` update (additive fields per Q-2).
  - WNAT_Hagen SHA-256 verified at start; exits non-zero on mismatch (calls `scripts/verify_mesh_pin.py` from T002).

- [ ] **T016** [US3] Implement half-edge v2 vectorized walk in `mesh_topology_halfedge.py` (extends T004's module). Add internal `_walk_face_vectorized()` method. Document in module docstring that v2 is selected via an internal flag, NOT exposed in the public kwarg (the kwarg only switches edgemap ↔ halfedge; v1 vs v2 is internal benchmark plumbing).

- [ ] **T017** [US3] Run T015. Capture the markdown table. Insert it into `docs/BENCHMARK.md` under a new "Half-Edge Variants (3-Way Comparison)" section. Add a footnote describing the median-of-3 + tracemalloc methodology.

- [ ] **T018** [US3] Compute the decision-trigger value per spec SC-005:
  - Full init half-edge v2 / EdgeMap > 2.0 → "investigate port"
  - 0.9 ≤ ratio ≤ 1.1 (within 10% per NFR-001 + Q-4) → "adopt"
  - ratio > 1.1 → "archive (experimental, kept in tree per Q-4)"
  - Print the trigger word on the last line of the markdown table.

**Checkpoint**: At this point, US3 ships. The PR for this branch carries enough data for v5.1's DECISION record to be a 1-page document — no further benchmarking required.

---

## Phase 6: Cross-Cutting Polish & Audit

**Purpose**: Final tightening before /speckit-analyze hands off to PR review.

- [ ] **T019** [P] [POLISH] Audit `mesh_topology_halfedge.py` for inline comments that violate constitution principle VIII (docstrings normative, code self-documenting). Remove WHAT-comments; keep WHY-comments only.
- [ ] **T020** [P] [POLISH] Run `mypy src/chilmesh/mesh_topology_halfedge.py` (if mypy available in env). Zero errors required.
- [ ] **T021** [P] [POLISH] Confirm no test file was modified to accommodate the half-edge backend (FR-004). `git diff main..HEAD -- tests/` should show only `test_halfedge_basic.py`, `test_halfedge_equivalence.py`, and `conftest.py` (the latter only for the env-var fixture from T003).
- [ ] **T022** [POLISH] Run `/speckit-analyze` (Phase 3 of plan.md). Produces `specs/007-v5.0-optimize-port-w-half-edge/analyze.md` with PASS/PASS-with-followups/BLOCK verdict. Open PR only if verdict is PASS or PASS-with-followups.

**Checkpoint**: PR-ready.

---

## Dependency Summary

```
T001 ──┐
T002 ──┼─→ T003 ──→ T004 ──→ T005 ──→ T006, T007 ──→ T008 ──→ T009, T010 ──→ T011 ──→ ...
                                                                                  ↓
                                                                    T012 ──→ T013, T014
                                                                                  ↓
                                                                    T015, T016 ──→ T017 ──→ T018
                                                                                                ↓
                                                                                  T019, T020, T021
                                                                                                ↓
                                                                                              T022
```

Parallelism opportunities:
- T013 + T014 (different files)
- T015 + T016 (different files; T015 imports T016 — actually T016 must precede T015 — adjust)
- T019 + T020 + T021 (different files)

Effort estimates (T-shirt):
- T001–T003: XS (15 min each)
- T004: L (4–6 hours)
- T005: M (2–3 hours)
- T006: M (2–3 hours)
- T007: S (1–2 hours)
- T008: S (1 hour integration + debug)
- T009–T011: XS (run + commit)
- T012: M (2–3 hours)
- T013: S (1 hour)
- T014: XS (15 min)
- T015: M (2–3 hours)
- T016: S (1–2 hours; depends on T004 quality)
- T017–T018: XS each
- T019–T022: XS each

**Total estimated effort**: ~20–25 hours focused work; ~3–5 wall-days at typical pace.
