# Tasks: Fix Skeletonization Layer Separation Invariant

**Input**: Design documents from `/specs/006-skeletonize-matlab-port/`
**Prerequisites**: plan.md, spec.md (both complete)

## Format: `[ID] [P?] [Story] Description with file path`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: US1=correctness, US2=backward compat, US3=visualization
- Tests are explicitly required (FR-006 SC-006); TDD order applies

---

## Phase 1: Setup

- [ ] T001 Verify spec is complete: read `specs/006-skeletonize-matlab-port/spec.md` and `specs/006-skeletonize-matlab-port/plan.md` to confirm all 3 clarifications are recorded and FR-001 through FR-013 are present.
- [ ] T002 Capture baseline metrics: run `python -m pytest tests/test_layers_annulus.py -v` and record current pass/fail status; run `python -c "from chilmesh import examples; m=examples.annulus(); print(m.n_layers)"` and similar for donut/structured/block_o; save these to a temporary scratch file for reference.

## Phase 2: Foundational (blocking prerequisites)

No foundational refactoring needed — the existing `_skeletonize()` method is a single function and the dict structure is preserved.

## Phase 3: User Story 1 — Correctness (P1)

**Goal**: The layer separation invariant holds across all fixtures.

**Independent Test**: A regression test that enumerates layer pairs (k, m) with |k-m| ≥ 2 and asserts disjoint vertex sets, parametrized over all 4 fixtures. The test must pass.

### TDD Phase: Write Tests First

- [ ] T003 [P] [US1] Create `tests/test_layer_separation.py` with a parametrized test `test_layer_separation_invariant` over fixtures (annulus, donut, structured, block_o). For each mesh, enumerate all (k, m) layer pairs with |k-m| ≥ 2 and assert that `set(vertices in OE[k] ∪ IE[k]) ∩ set(vertices in OE[m] ∪ IE[m])` is empty. Mark the test `@pytest.mark.xfail(reason="Issue #74: pre-fix")` initially; will be removed after T005.

### Implementation Phase

- [ ] T004 [US1] Read the current `_skeletonize()` method in `src/chilmesh/CHILmesh.py` (lines ~767-872) for context.
- [ ] T005 [US1] Replace `_skeletonize()` body in `src/chilmesh/CHILmesh.py` with a faithful Python port of the MATLAB `meshLayers` function. Specifically:
   1. Initialize `Edge2VertIDs = self.adjacencies["Edge2Vert"].copy()` and `Edge2ElemIDs = self.adjacencies["Edge2Elem"].copy()`.
   2. Reset `self.layers = {"OE": [], "IE": [], "OV": [], "IV": [], "bEdgeIDs": []}`.
   3. Loop while `(Edge2ElemIDs >= 0).any()`. (Python uses -1 sentinel, not 0.)
   4. **Per-iteration steps** (matching MATLAB exactly, with line-comment references to the MATLAB source):
      - If `iL == 0`: `iLbEdgeIDs = boundary_edges()`; else `iLbEdgeIDs = np.where(np.sum(Edge2ElemIDs >= 0, axis=1) == 1)[0]`
      - Compute `OV = np.unique(filtered_vertices(Edge2VertIDs[iLbEdgeIDs]))` (filter out -1 padding per Q3)
      - Append OV, bEdgeIDs to layers
      - Compute `OE = np.unique(active_elements(Edge2ElemIDs[iLbEdgeIDs]))` (active = `>= 0`)
      - Append OE; flag `Edge2ElemIDs[np.isin(Edge2ElemIDs, OE)] = -1`
      - Compute "edges associated with OV": `iLbEdgeIDs2 = np.where(np.any(np.isin(Edge2VertIDs, OV), axis=1))[0]`
      - Compute `IE = np.unique(active_elements(Edge2ElemIDs[iLbEdgeIDs2]))`
      - Append IE; flag `Edge2VertIDs[np.isin(Edge2VertIDs, OV)] = -1`; flag `Edge2ElemIDs[np.isin(Edge2ElemIDs, IE)] = -1`
      - Compute `IV = setdiff(unique_filtered(connectivity_list[OE ∪ IE]), OV)` (filter -1 padding)
      - Append IV; increment `iL`
   5. After loop: `self.n_layers = iL`
- [ ] T006 [US1] Remove the `@pytest.mark.xfail` from `tests/test_layer_separation.py::test_layer_separation_invariant` and verify it now passes for all 4 fixtures.

## Phase 4: User Story 2 — Backward Compatibility (P1)

**Goal**: All existing tests pass; layer counts are pinned to MATLAB-correct values.

**Independent Test**: Full pytest suite passes, including `test_layers_annulus.py` and `test_invariants.py`.

- [ ] T007 [US2] Run `python -m pytest tests/ --tb=line` and capture all failures. Categorize each failure as: (a) genuine regression to fix, (b) test asserting a buggy layer count to update per Q2 Option A, or (c) pre-existing failure unrelated to skeletonization (e.g., `_detect_element_type`).
- [ ] T008 [US2] For each test in category (b) — typically `tests/test_layers_annulus.py::test_structured_grid_layers` and possibly `test_annulus_layers` — update the assertion to the new MATLAB-correct value. Add a comment explaining the change references issue #74.
- [ ] T009 [P] [US2] Create `tests/test_matlab_layer_counts.py` with a parametrized test `test_layer_count_matches_matlab_reference` that pins the new per-fixture layer counts. Initial values are captured from running `mesh.n_layers` after T005. Mark the test as the source of truth for layer counts going forward.
- [ ] T010 [US2] Run `python -m pytest tests/ -v` and verify all previously passing tests still pass (or are documented as updated per category-b).

## Phase 5: User Story 3 — Visualization (P2)

**Goal**: Public README image shows medially-correct layer rings.

**Independent Test**: Regenerate `tests/output/annulus_quickstart.png` and visually inspect column 2 of all 4 rows.

- [ ] T011 [US3] Run `python generate_4row_admesh.py` to regenerate `tests/output/annulus_quickstart.png` with the new skeletonization.
- [ ] T012 [US3] Visually inspect the regenerated image (display via Read tool) and verify column 2 (Layers) of all 4 rows shows clean concentric ring transitions with no Layer-N triangles touching Layer-0 triangles for N ≥ 2.

## Phase 6: Polish & Cross-Cutting

- [ ] T013 Update CHANGELOG.md (or create one if absent) with a brief note: "Fixed #74: skeletonization layer separation invariant (replaced with MATLAB-equivalent algorithm). Layer counts may change for existing meshes; this is expected and corrects a long-standing bug."
- [ ] T014 Update issue #74 with a comment linking to the spec, plan, and final commit; close the issue.
- [ ] T015 Commit all changes with a single atomic commit message (or two: one for the algorithm fix + tests, one for the visualization update).

## Dependencies

- T001, T002 → T003 (need spec confirmed before writing tests)
- T003 → T004 → T005 → T006 (TDD: test first, then implementation, then test verification)
- T005 → T007 (need new behavior to know what tests fail)
- T007 → T008 (need failure categorization to update tests)
- T005, T008 → T009 (need final layer counts to pin them)
- T010 (after all test updates) → T011 → T012 (visualization after tests pass)
- T010 + T012 → T013 → T014 → T015 (final cleanup and documentation)

## Parallel Execution Examples

- T003 and T009 can be created in parallel (different test files, no shared state). However, T009 depends on T005's output for the actual values, so create T003 first, run through to T005, then create T009.
- T013 (CHANGELOG update) can run in parallel with T011 (image regeneration).

## Implementation Strategy

**MVP scope = User Story 1 only** (T001-T006). Once the layer separation invariant is fixed and the regression test passes, the core correctness goal is achieved. US2 (backward compat / test updates) and US3 (visualization) follow naturally as cleanup.

If US1 alone takes longer than expected (e.g., performance regression on block_o), defer US3 to a follow-up commit.
