# Tasks: ADMESH Warm-Start Truss Optimization

**Feature**: `005-admesh-warm-start-truss`
**Plan**: [plan.md](plan.md)
**Date**: 2026-05-02
**Phase Status**: Phase 0 (Research) ✓ | Phase 1 (Design) ✓ | Phase 2 (Tasks) → Implementation

---

## Overview

Implementation tasks for the ADMESH warm-start truss adapter. Divided into 6 phases:
- **Phase 1**: Setup & prerequisites
- **Phase 2**: Core adapter implementation (array form)
- **Phase 3**: Vendored truss loop integration
- **Phase 4**: High-level CHILmesh wrapper
- **Phase 5**: Comprehensive test suite
- **Phase 6**: Demo script & visualization

**User Stories**:
- US1: Boundary-preserving optimization (MVP — array form + CHILmesh wrapper)
- US2: Four-row demo visualization (new layout per Q2=d)
- US3: Generic input support (array form + domain-agnostic SDF)

**Dependencies**: None blocking; all work is parallelizable after US1 core implementation.

**Estimated Effort**: ~30-40 implementation tasks; ~20 test/validation tasks.

**Success Criteria** (from plan.md):
- SC-001: Annulus median quality ≥ 0.60 after warm-start
- SC-002: Bit-exact boundary preservation (`np.array_equal` verification)
- SC-003: Adapter completes < 30s on 580-element annulus
- SC-007/SC-008: Domain-agnostic & source-agnostic (donut + array form)

---

## Phase 1: Setup & Prerequisites

- [ ] T001 [P] Create skeleton modules: `src/chilmesh/admesh_warmstart.py` and `src/chilmesh/_vendor_admesh_truss.py` with placeholder docstrings
- [ ] T002 [P] Verify ADMESH import chain: Test that `from admesh.distmesh import distmesh2d` works when ADMESH is pip-installed or available at pinned SHA `05bc68f`
- [ ] T003 Create module docstring for `admesh_warmstart.py` documenting: FR-001a, FR-001b, FR-002, FR-005-007, FR-016 with cross-references to contracts/api-contract.md
- [ ] T004 Add header comment to `_vendor_admesh_truss.py`: "Vendored from admesh.distmesh.distmesh2d at commit 05bc68f. Byte-identical except for warm-start preamble. Remove when ADMESH-B (public warm-start entry point) lands."

---

## Phase 2: Core Adapter Implementation (Array Form)

### Boundary Identification & Validation (FR-002, FR-005-007)

- [ ] T005 [P] Implement `_infer_boundary_from_triangles(triangles: ndarray) -> ndarray` helper: Edges appearing in exactly one triangle → boundary vertex indices; ensure deterministic ordering via `np.unique` (FR-007)
- [ ] T006 [P] Implement validation block in `optimize_with_admesh_truss_arrays()`:
  - V_TRI: Raise `NotImplementedError` if `triangles.shape[1] != 3`
  - V_AREA: Raise `ValueError` if any triangle has signed area ≤ 0; list offending indices (FR-005)
  - V_BND_SDF: Raise `ValueError` if any boundary point has `|sdf(p)| ≥ sdf_tolerance`; include max |sdf| and offending index (FR-006)
  - V_BBOX: Raise `ValueError` if any input point falls outside bbox; include offending indices
  - See contracts/api-contract.md for exact error message formats
- [ ] T007 Compute h0 if None: Mean of input edge lengths (from points and triangles)
- [ ] T008 Infer bbox if None: Min/max of all points with 5% margin on all sides; handle degenerate cases

### Core Algorithm (FR-003, FR-004, FR-008)

- [ ] T009 [P] Implement main function signature `optimize_with_admesh_truss_arrays(points, triangles, sdf, size_fn=None, boundary_indices=None, h0=None, bbox=None, dptol=1e-3, ttol=0.1, Fscale=1.2, deltat=0.2, geps_factor=1e-3, niter=500, seed=0, sdf_tolerance=1e-6, enforce_non_degradation=True) -> (points_out, triangles_out)`
- [ ] T010 Prepare normalized ValidatedInput: Extract `points_xy = points[:, :2]`, identify boundary and interior indices, compute `boundary_xy` and `interior_xy` subsets
- [ ] T011 Call vendored truss loop: `distmesh2d_warmstart(boundary_xy, interior_xy, sdf, size_fn, h0, bbox, dptol, ttol, Fscale, deltat, geps_factor, niter, seed)` → `(points_opt, triangles_opt)`
- [ ] T012 Implement non-degradation guard (FR-011):
  - Compute `median_q_in = np.median(elem_quality(triangles, points[:, :2]))`
  - Compute `median_q_out = np.median(elem_quality(triangles_opt, points_opt))`
  - If `enforce_non_degradation=True` and `median_q_out < median_q_in`: emit `RuntimeWarning` and return input unchanged
  - If `enforce_non_degradation=False`: always return optimized output
- [ ] T013 Return boundary-preserving output: `(points_out, triangles_out)` where `points_out[:B] == input_points[boundary_indices]` bit-exactly (FR-008)

### Determinism & RNG (FR-009)

- [ ] T014 Wire `seed` parameter through to vendored truss loop's RNG initialization; document that determinism holds for identical inputs + identical seed

---

## Phase 3: Vendored Truss Loop Integration

- [ ] T015 Copy inner loop from `admesh.distmesh.distmesh2d` (lines 140-220 of upstream at commit 05bc68f) into `_vendor_admesh_truss.distmesh2d_warmstart()` function
- [ ] T016 Modify preamble (only change): Replace `_initial_distribution() + _rejection_method()` with:
  ```python
  p = np.vstack([pfix_boundary, interior_initial])
  nfix = len(pfix_boundary)
  ```
- [ ] T017 Verify vendored loop structure matches upstream exactly: Same force computation, same Delaunay call, same re-triangulation cadence (ttol check), same convergence check (dptol on interior movement)
- [ ] T018 Implement `distmesh2d_warmstart(pfix_boundary, interior_initial, fd, fh, h0, bbox, dptol, ttol, Fscale, deltat, geps_factor, niter, seed)` wrapper:
  - Combine pfix and interior into working point array `p`
  - Run truss iterations with boundary pinned
  - Return `(p_final, t_final)` where triangles are from final Delaunay

---

## Phase 4: High-Level CHILmesh Wrapper (FR-001b)

- [ ] T019 Implement `optimize_with_admesh_truss(mesh: CHILmesh, sdf: Callable, size_fn: Optional[Callable]=None, **kwargs) -> CHILmesh`:
  - Extract `points` from `mesh.points`
  - Extract `triangles` from `mesh.connectivity_list`
  - Identify boundary via `mesh.boundary_edges() → Edge2Vert → np.unique` (FR-002)
  - Call `optimize_with_admesh_truss_arrays(points, triangles, sdf, size_fn, boundary_indices, **kwargs)`
  - Wrap output in fresh `CHILmesh(connectivity=triangles_out, points=column_stack([points_out, zeros]))`
  - Return new CHILmesh instance
- [ ] T020 Add docstring to high-level form linking to array form and documenting the thin-wrapper relationship
- [ ] T021 Update `src/chilmesh/__init__.py`: Import and re-export both `optimize_with_admesh_truss` and `optimize_with_admesh_truss_arrays`

---

## Phase 5: Comprehensive Test Suite

### Unit Tests: Boundary Preservation & Non-Degradation (FR-008, FR-011)

- [ ] T022 [P] Create `tests/test_admesh_warmstart.py` with test fixtures: annulus (simple), donut (domain-agnostic), structured (grid)
- [ ] T023 [P] Test V_BND (bit-exact boundary preservation): Run adapter on annulus, verify `np.array_equal(output[:B], input[boundary_indices])`
- [ ] T024 [P] Test V_BND on multiple domains: Repeat V_BND test on donut and structured
- [ ] T025 [P] Test non-degradation guard: Run adapter with known-good input; verify output quality ≥ input quality (or return input if not)
- [ ] T026 [P] Test non-degradation override: Run with `enforce_non_degradation=False`; verify optimizer output returned even if worse

### Unit Tests: Error Handling (FR-005-007, FR-015)

- [ ] T027 [P] V_TRI: Pass quad mesh (4-vertex elements) → expect `NotImplementedError` with exact message from contracts/api-contract.md
- [ ] T028 [P] V_AREA: Pass triangles with zero or negative area → expect `ValueError` listing indices
- [ ] T029 [P] V_BND_SDF: Pass boundary point off SDF zero set → expect `ValueError` with max |sdf| and index
- [ ] T030 [P] V_BBOX: Pass points outside inferred bbox → expect `ValueError` listing indices
- [ ] T031 [P] ADMESH import error: Mock missing ADMESH import → expect `ImportError` mentioning pinned SHA

### Unit Tests: Input Source Agnosticity (SC-007)

- [ ] T032 [P] Array form test (no CHILmesh): Use scipy.spatial.Delaunay to create raw triangulation, call array form with explicit boundary_indices
- [ ] T033 [P] CHILmesh form test: Load annulus fixture, call CHILmesh form, verify output is fresh CHILmesh instance with adjacencies recomputed

### Unit Tests: Domain Agnosticity (SC-008, Q3=d)

- [ ] T034 [P] Test on annulus with constant size_fn (same as Row 2 demo)
- [ ] T035 [P] Test on donut with graded size_fn: Finer near boundaries, coarser in middle (spec 004 pattern)
- [ ] T036 [P] Custom domain test: Unit square with user-supplied SDF, verify boundary preservation on square corners

### Integration Tests: Determinism & Performance (FR-009, SC-003)

- [ ] T037 [P] Determinism test: Run same adapter call twice with `seed=0`; verify `np.array_equal(out1[0], out2[0])` and `np.array_equal(out1[1], out2[1])`
- [ ] T038 [P] Performance baseline: Time adapter on annulus (580 elements); verify completes in < 30s (SC-003)

### Quality Metrics Tests (SC-001, SC-005)

- [ ] T039 [P] Annulus quality improvement: Verify median quality ≥ 0.60 after warm-start (SC-001, Q5=b)
- [ ] T040 [P] Quality does not regress: For all test domains, verify `median_q_out >= median_q_in` when `enforce_non_degradation=True`

---

## Phase 6: Demo Script & Visualization

### Restructure generate_4row_admesh.py (Q2=d, FR-012-013)

- [ ] T041 Modify `generate_4row_admesh.py` to new 4-row layout:
  - Row 1: `chilmesh.examples.annulus()` (raw Delaunay)
  - Row 2: `optimize_with_admesh_truss(row1, ANNULUS_SDF, size_fn=None, seed=0)` (warm-start of row 1)
  - Row 3: `row2.smooth_mesh(method='fem', acknowledge_change=True)` (FEM smoother applied to row 2)
  - Row 4: `admesh.quad_prep.smooth_for_quadrangulation(row2.points, row2.triangles, ANNULUS_SDF, ...)` (right-isoceles applied to row 2)
  - Document that rows 3 and 4 are siblings off row 2, not sequential
- [ ] T042 Remove old Row 3 (fresh ADMESH from bbox) — this pattern is no longer needed
- [ ] T043 Update subplot titles per contracts/visualization-output.md:
  - Row 1: "Raw Delaunay"
  - Row 2: "+ ADMESH Truss (warm-start)"
  - Row 3: "Row 2 + FEM Smoother"
  - Row 4: "Row 2 + Right-Isoceles Smoother"

### Demo Assertions (FR-013, V_BND-V_TRUSS_INVOKED)

- [ ] T044 Replace spec-004 V1-V7 assertions with new spec-005 assertions:
  - **V_BND**: `np.array_equal(row2.points[boundary_indices], row1.points[boundary_indices])` — warm-start preserved boundary
  - **V_BND_PROP**: Row 3 and Row 4 boundary points equal Row 2's within each smoother's tolerance
  - **V_QI**: `median_q(row2) > median_q(row1)` — warm-start improved quality
  - **V_CONN**: All four rows have positive triangle areas, no duplicate vertex indices
  - **V_CHAIN**: Row 3 input was row 2 (verify `id(row2)` or hash equality); Row 4 input was row 2 (NOT a fresh ADMESH)
  - **V_TRUSS_INVOKED**: `distmesh2d_warmstart` was called exactly once for Row 2 (track via module-level flag)
- [ ] T045 If any assertion fails, raise `RuntimeError` and DO NOT write PNG (fail-loud, per spec 004 pattern)
- [ ] T046 On all assertions passing, render 12 subplots (4 rows × 3 columns) and save PNG to `tests/output/annulus_quickstart.png`

### Colormaps (from spec 004, preserved in spec 005)

- [ ] T047 Column 1 (mesh): Black wireframe, no fill, no colorbar (triplot only)
- [ ] T048 Column 2 (layers): Parula 64-stop ListedColormap (`matplotlib.cm.get_cmap('parula')`), discrete colorbar with n_layers ticks labeled "Layer"
- [ ] T049 Column 3 (quality): cool_r (reversed cool: red=0, blue=1.0), continuous colorbar [0,1] labeled "Quality"

### Regenerate PNG Artifact (FR-012)

- [ ] T050 Run restructured `generate_4row_admesh.py` to regenerate `tests/output/annulus_quickstart.png`
- [ ] T051 Verify PNG dimensions ~1500×1800 px (figsize 15×18 @ 100 dpi), file size < 5 MB, format 8-bit RGBA non-interlaced
- [ ] T052 Commit regenerated PNG to branch

---

## Phase 7: Documentation & Polish

### Verify Quickstart Examples (FR-018)

- [ ] T053 [P] Example 1 (annulus CHILmesh form): Copy example from quickstart.md, run interactively, verify boundary preservation assertion passes
- [ ] T054 [P] Example 2 (donut with graded size_fn): Copy example from quickstart.md, run interactively, verify output quality and boundary
- [ ] T055 [P] Example 3 (raw arrays from non-CHILmesh source): Copy example from quickstart.md, run interactively, verify 4 corners of square preserved

### Update README & Cross-References

- [ ] T056 Update `README.md` figure caption to reference new 4-row pipeline: "Row 1: Raw Delaunay, Row 2: ADMESH warm-start, Row 3: FEM smoother, Row 4: Right-isoceles smoother"
- [ ] T057 Add inline references in README or a "How It Works" section linking to quickstart.md and api-contract.md
- [ ] T058 Verify README image link (`tests/output/annulus_quickstart.png`) is correct and image renders

### Cross-Repo Issue Filing (per spec 004 pattern)

- [ ] T059 File ADMESH issue ADMESH-A: "`admesh.routine` imports undefined `MeshOutput` on main HEAD" — document broken import, link to daily-issue-fixing branch fix, note that CHILmesh currently bypasses routine
- [ ] T060 File ADMESH issue ADMESH-B: "Add public warm-start entry point (`distmesh2d_warmstart` or `initial_points` parameter)" — note that CHILmesh currently vendors the inner loop pending this
- [ ] T061 File ADMESH issue ADMESH-C: "Document `pfix` bit-exact preservation guarantee in distmesh2d docstring" — cite source inspection evidence

### Final Verification

- [ ] T062 Run full test suite: `pytest tests/test_admesh_warmstart.py -v` — all tests pass
- [ ] T063 Run demo script: `python generate_4row_admesh.py` — PNG generated without assertion failures
- [ ] T064 Type hints: Verify `optimize_with_admesh_truss` and `optimize_with_admesh_truss_arrays` have complete type annotations
- [ ] T065 Import check: Verify `from chilmesh import optimize_with_admesh_truss, optimize_with_admesh_truss_arrays` works from clean Python session

---

## Parallel Opportunities

**Phase 2 parallelization** (after T004 skeleton created):
- T005-T008: Boundary identification & validation
- T009-T013: Core algorithm (all depend on T005-T008)

**Phase 5 parallelization** (after T022 fixtures created):
- T023-T026: Boundary preservation tests
- T027-T031: Error handling tests
- T032-T035: Domain/source agnosticity tests
- T037-T040: Integration tests

**Phase 6 parallelization** (after T043 layout created):
- T044-T049: Assertions & colormaps (independent of PNG rendering)
- T053-T055: Quickstart examples (independent of demo restructuring)

---

## Dependencies

- T001: no deps
- T002: depends on ADMESH being available (environmental)
- T003-T004: depend on T001
- T005-T008: depend on T001
- T009-T013: depend on T005-T008
- T015-T018: depend on T001 (skeleton exists)
- T019-T021: depend on T009-T018 (core functions exist)
- T022: no deps on implementation (test fixtures already available)
- T023-T040: depend on T019-T021 (functions ready to test)
- T041-T052: depend on T019-T021 (functions ready for demo)
- T053-T065: depend on T050 (all work complete)

**Sequential blockers**:
1. T001 (create skeleton)
2. T005-T021 (implement all public functions)
3. T022-T052 (tests & demo)
4. T053-T065 (final verification)

**Can run in parallel after each sequential blocker**:
- After T001: All phases can start skeleton review
- After T005-T008: Error handling tests (T027-T031) can run independently
- After T009-T021: All unit tests (T023-T040) can run in parallel
- After T041-T049: All demo verification (T053-T055) can run in parallel

---

## Success Criteria (Acceptance)

All tasks complete when:

1. **V_BND**: `np.array_equal(output_boundary, input_boundary)` passes for all test domains
2. **V_QI**: Annulus median quality ≥ 0.60 after warm-start
3. **SC-003**: Adapter completes < 30s on 580-element annulus
4. **V_TRUSS_INVOKED**: Demo script verifies truss was actually invoked for Row 2
5. **All tests pass**: `pytest tests/test_admesh_warmstart.py -v` with 40+ test cases
6. **PNG regenerated**: New 4-row visualization in `tests/output/annulus_quickstart.png`
7. **Examples work**: All three quickstart.md examples run without error
8. **Imports work**: `from chilmesh import optimize_with_admesh_truss, optimize_with_admesh_truss_arrays`

---

**Status**: Ready for implementation.  
**Date**: 2026-05-02  
**Phase**: 2 → Implementation

