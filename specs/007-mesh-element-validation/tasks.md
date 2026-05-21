---
description: "Task list for spec 007 — mesh element validity test suite"
---

# Tasks: Mesh Element Validity Test Suite

**Input**: `/specs/007-mesh-element-validation/{spec,plan,research,data-model,quickstart,contracts}`
**Prerequisites**: spec.md, plan.md (both complete)

**Tests**: Tests ARE the deliverable in this spec. Every implementation task ships its test FIRST and verifies failure before writing implementation code (constitution principle III).

**Organization**: Tasks grouped by user story. US1 + US2 + US3 are all P1; US4 is P2 polish.

## Format: `[ID] [P?] [Story] Description`

- `[P]` = can run in parallel (different files, no dependencies)
- `[US#]` = user story tag (US1 = element-type, US2 = geometric, US3 = degenerate, US4 = pytest UX)

---

## Phase 1: Setup

**Purpose**: Scaffold `tests/_validity/` helper package and pin the contract.

- [ ] **T001** Create `tests/_validity/` directory with empty `__init__.py` that re-exports `validate_mesh_elements`, `MeshValidityReport`, `Violation`, `InformationalNote`.
- [ ] **T002** [P] Create `tests/_validity/README.md` (≤ 30 lines) pointing at `specs/007-mesh-element-validation/spec.md` and noting "test-suite-only, not public API (see issue #142)".
- [ ] **T003** [P] Copy frozen contract from `specs/007-mesh-element-validation/contracts/validator.py` (signatures only) into `tests/_validity/types.py` (data classes) and a stub `tests/_validity/validator.py` that raises `NotImplementedError`.

**Checkpoint**: `from tests._validity import validate_mesh_elements, MeshValidityReport` succeeds. Stub raises on call.

---

## Phase 2: Foundational (BLOCKS all user stories)

**Purpose**: Core predicates + broadphase that every per-story rule depends on.

- [ ] **T010** [P] [foundation] Write `tests/_validity/test_predicates.py` (UNIT-level, not in main suite) covering: `_orient` sign on cw/ccw/collinear triples; `segment_proper_cross` on shared-endpoint exemption, collinear-no-overlap, proper cross, t-junction (touching at interior point — counts as crossing); `point_strictly_in_polygon` on square + bowtie + degenerate triangle. **Tests MUST fail at this point** (no implementation yet).
- [ ] **T011** [P] [foundation] Write `tests/_validity/test_broadphase.py`: build a 4-element grid, assert all 6 unique pairs are emitted exactly once; assert pair count stays under `n * 8` on a 100-element random fixture.
- [ ] **T012** [foundation] Implement `tests/_validity/predicates.py` (`bbox_diag`, `classify_element`, `_orient`, `segment_proper_cross`, `point_strictly_in_polygon`, `is_self_intersecting_quad`). Pure numpy. Make T010 pass.
- [ ] **T013** [foundation] Implement `tests/_validity/broadphase.py` (`UniformGridIndex` build + `candidate_pairs`). Make T011 pass.
- [ ] **T014** [foundation] Implement `tests/_validity/fixtures.py`: `corrupt_to(mesh, row_id, new_row)` helper + each synthetic fixture from data-model.md (`bowtie_quad_mesh`, `interior_triangle_mesh`, `pentagon_mesh` via test-double, `overlapping_quads_mesh`, `edge_crossing_mesh`). Verify each fixture loads without exception.

**Checkpoint**: Predicates, broadphase, and fixtures all exist and have unit-level test coverage. No `validate_mesh_elements` yet.

---

## Phase 3: User Story 1 — Element-type composition (Priority: P1) MVP

**Goal**: Detect interior triangles, unsupported arities, layer-membership violations.

**Independent test**: Run validator on `annulus` (all triangles in layer 0 → pass), `interior_triangle_mesh` (interior tri planted → `INTERIOR_TRIANGLE_FORBIDDEN`), `pentagon_mesh` (5-vertex element → `UNSUPPORTED_ELEMENT_ARITY`).

### Tests first

- [ ] **T020** [P] [US1] Add `tests/test_mesh_element_validity.py` with two parametrized tests: `test_builtins_element_type[fixture]` over the four built-in fixtures (expect `ok=True` for the type-composition subset) and `test_synthetic_negatives_element_type[case]` over `interior_triangle_mesh`, `pentagon_mesh` (expect specific violation category). Mark `pytest.mark.skip(reason="pending validator")`.

### Implementation

- [ ] **T021** [US1] Implement element-type rules in `validator.py`:
  - FR-003: arity classification, emit `UNSUPPORTED_ELEMENT_ARITY` for ≥5-vertex rows.
  - FR-004: emit `DEGENERATE_QUAD_DUPLICATE_VERTEX` as informational note (NOT violation).
  - FR-005: zero-distinct-vertex triangle → handled by ID-distinctness check.
  - FR-006: each TRI must have ≥1 vertex in `mesh.boundary_vertices()`; else `INTERIOR_TRIANGLE_FORBIDDEN`.
  - FR-007: auto-trigger `_skeletonize()` if `mesh.layers` not yet populated; each TRI must belong to `layers["OE"][0] ∪ layers["IE"][0]`; else `INTERIOR_LAYER_TRIANGLE_FORBIDDEN`. Emit `LAYERS_AUTO_TRIGGERED` informational note when validator triggered the computation.
- [ ] **T022** [US1] Unskip T020. Verify all four built-in fixtures pass and both synthetic negative fixtures fail with the correct category.

**Checkpoint**: US1 acceptance scenarios 1-4 from spec.md pass.

---

## Phase 4: User Story 2 — Geometric validity (Priority: P1)

**Goal**: Detect bowtie quads, interior-overlap pairs, edge-crossing pairs. Planarity check.

**Independent test**: Run validator on `bowtie_quad_mesh`, `overlapping_quads_mesh`, `edge_crossing_mesh` and verify each produces its expected violation category.

### Tests first

- [ ] **T030** [P] [US2] Extend `tests/test_mesh_element_validity.py` with `test_synthetic_negatives_geometric[case]` over the three geometric negative fixtures. Mark skip until T031-T033 land.
- [ ] **T031** [P] [US2] Add `test_planarity` parametrized over: (a) built-in fixtures (`mesh.points.shape[1] == 2` → pass), (b) a synthetic 3D-coord fixture with varying `z` → expect `NON_PLANAR_MESH`.

### Implementation

- [ ] **T032** [US2] Implement planarity (FR-008) and bowtie (FR-009) rules in `validator.py`. Vectorize the bowtie check via `is_self_intersecting_quad` over all QUAD rows.
- [ ] **T033** [US2] Implement element-element non-overlap (FR-011) and broadphase (FR-012). Algorithm:
  1. Build edge-sharing map per research.md.
  2. If `n_elems > 5000`, build `UniformGridIndex`; else use direct `O(n²)` pair iteration.
  3. For each candidate pair `(i, j)` that does NOT share an edge:
     - Test `EDGE_CROSSING`: for each edge of `i` × each edge of `j`, call `segment_proper_cross` (with shared-vertex exemption).
     - Test `INTERIOR_OVERLAP`: any vertex of `i` strictly inside `j` polygon, or vice versa.
- [ ] **T034** [US2] Unskip T030 and T031. Verify all three geometric synthetic negatives fail with the correct category and all four built-in fixtures pass.

**Checkpoint**: US2 acceptance scenarios 1-5 from spec.md pass.

---

## Phase 5: User Story 3 — Degenerate-quad acceptance (Priority: P1)

**Goal**: Confirm degenerate ≠ violation. Cross-check against `tests/test_degeneracy.py` scenarios.

**Independent test**: Build the same fixtures that `tests/test_degeneracy.py` uses (padded-triangle, duplicate-vertex quad, mixed-element mesh) and assert `report.ok == True` and `DEGENERATE_*` notes are present.

### Tests first

- [ ] **T040** [US3] Add `test_degenerate_accepted` to `tests/test_mesh_element_validity.py` covering:
  - Padded-triangle quad: pass, expect `DEGENERATE_QUAD_DUPLICATE_VERTEX` note (it's the `d==a||b||c` form).
  - Zero-area collinear quad: pass, expect `DEGENERATE_ZERO_AREA` note.
  - Bowtie with near-zero area: FAIL with `SELF_INTERSECTING_QUAD` (verifies area is NOT used as a self-intersection proxy).

### Implementation

- [ ] **T041** [US3] Add zero-area detection in `validator.py`: compute signed area per element; if `|area| < tol_effective * bbox_diag`, emit `DEGENERATE_ZERO_AREA` note.
- [ ] **T042** [US3] Re-verify all existing `tests/test_degeneracy.py` tests still pass after the new module is importable (no spillover side effects).

**Checkpoint**: US3 acceptance scenarios 1-4 from spec.md pass.

---

## Phase 6: User Story 4 — Pytest UX & runtime budget (Priority: P2)

**Goal**: Pytest output is actionable; runtime fits the budget.

### Tests first

- [ ] **T050** [US4] Add `test_failure_message_format` (planted bowtie in `annulus`): assert assertion message contains element ID, violation category `SELF_INTERSECTING_QUAD`, two edge IDs, four vertex coords (per SC-006).
- [ ] **T051** [US4] Add `test_runtime_budget` parametrized over the four built-in fixtures. Assert `report.runtime_s` < `{annulus: 5, donut: 5, structured: 5, block_o: 45}` (SC-004).

### Implementation

- [ ] **T052** [US4] Implement FR-016 / FR-017 in `validator.py`: build assertion message that aggregates all violation categories (up to 10 per category) when invoked from the pytest harness. Wire into `test_mesh_element_validity.py` via `assert report.ok, _format_failures(report)`.
- [ ] **T053** [US4] Profile `block_o` if T051 fails. If broadphase is the bottleneck, switch to `scipy.spatial.cKDTree`-based broadphase per plan.md risk register.

**Checkpoint**: US4 acceptance scenarios 1-3 from spec.md pass.

---

## Phase 7: Polish

- [ ] **T060** [P] Add `LAYERS_AUTO_TRIGGERED` to the violation-category table in spec.md (was missing from the table; covered in FR-007).
- [ ] **T061** [P] `CHANGELOG.md` entry under "Internal" (NOT "Added"): "Internal test-suite helper `tests/_validity/` (see #142 for promotion discussion)".
- [ ] **T062** [P] Run `pytest -v` and confirm all 288 pre-existing tests still pass (SC-005).
- [ ] **T063** Final commit + push.

---

## Dependencies & Execution Order

- T001-T003 (Setup) → no dependencies.
- T010-T014 (Foundational) → blocks Phase 3+. T010-T011 run in parallel; T012 needs T010, T013 needs T011, T014 is independent.
- T020-T022 (US1) → blocks demo of MVP. T020 first, T021 implements, T022 enables.
- T030-T034 (US2) → independent of US1 logic but shares the same `validator.py`. Land US1 first to keep the diff reviewable.
- T040-T042 (US3) → independent of US1/US2; cross-tests existing degeneracy tests.
- T050-T053 (US4) → depends on US1+US2+US3 being landed (needs a full `report`).
- T060-T063 (Polish) → last.

### Parallel opportunities

- T002, T003 in parallel.
- T010, T011 in parallel.
- T020, T030, T040 (test scaffolding for each story) in parallel, but DO NOT unskip until the corresponding implementation lands.

---

## Mapping: FR / SC → Tasks

| FR/SC | Task(s) |
|-------|---------|
| FR-001 entry point | T003, T021 |
| FR-002 report shape | T003 (types.py) |
| FR-003 arity | T021 |
| FR-004 duplicate-vertex note | T021 |
| FR-005 tri distinctness | T021 |
| FR-006 boundary-vertex tri | T021 |
| FR-007 layer-membership + auto-skeletonize | T021 |
| FR-008 planarity | T031, T032 |
| FR-009 bowtie | T032 |
| FR-010 tri non-self-intersect (no-op) | T032 |
| FR-011 overlap + edge-cross | T033 |
| FR-012 broadphase | T013, T033 |
| FR-013 tolerance | T012 |
| FR-014 degenerate notes | T041 |
| FR-015 synthetic fixtures | T014 |
| FR-016 message cap | T052 |
| FR-017 pytest shape | T052 |
| SC-001 built-ins pass | T022, T034, T042 |
| SC-002 synthetic negatives fail correctly | T022, T034, T040 |
| SC-003 file location | T001, T003 |
| SC-004 runtime budget | T051 |
| SC-005 existing tests pass | T062 |
| SC-006 failure-message format | T050 |
| SC-007 ≥100× broadphase reduction | research.md; verified empirically in T053 |
| SC-008 issue opened | DONE (issue #142) |

---

## Notes

- All work lands on branch `claude/mesh-quad-triangle-spec-IIvKw` per user override in this session.
- Each task ≤ 1 logical commit. Commit messages: `feat(007):` for new code, `test(007):` for test additions, `chore(007):` for polish.
- No new package deps. numpy + pytest only.
- If a task expands beyond its checkpoint scope, STOP and update tasks.md before continuing.
