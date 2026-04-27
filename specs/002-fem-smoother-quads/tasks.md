# Implementation Tasks: FEM Smoother for Quad & Mixed-Element Meshes

**Feature**: Extend FEM Smoother for Quad & Mixed-Element Meshes (Issue #4)  
**Specification**: [spec.md](spec.md)  
**Plan**: [plan.md](plan.md)  
**Branch**: `claude/youthful-goldberg-ueQ9R`  
**Created**: 2026-04-27

---

## Overview

This document breaks Issue #4 into executable tasks organized by user story. Each story (triangle, quad, mixed) is independently testable. Implement in order: Phase 2 (foundational), then Phase 3-5 (user stories), then Phase 6 (polish).

**MVP Scope**: Complete Phase 2 + Phase 3 (backward compatibility + tests) = minimal viable product  
**Full Scope**: Phases 2-6 (all three mesh types + validation + performance)

---

## Dependency Graph

```
Phase 1: Setup (none needed - library project)
  ↓
Phase 2: Foundational (blocking all user stories)
  ├→ T001: Understand Zhou & Shimada formulation
  ├→ T002: Design quad stiffness matrices
  └→ T003: Refactor element type detection logic
  ↓
Phase 3: [US1] Triangle Backward Compat (independent)
  ├→ T004: [P] Write regression tests (triangle)
  └→ T005: Refactor direct_smoother to element-type dispatcher
  ↓
Phase 4: [US2] Quad Mesh Support (independent, depends on Phase 2)
  ├→ T006: [P] Write tests for quad meshes (structured fixture)
  ├→ T007: Implement quad stiffness matrix calculation
  └→ T008: Validate quad smoother output
  ↓
Phase 5: [US3] Mixed-Element Support (independent, depends on Phase 2)
  ├→ T009: [P] Write tests for mixed-element meshes (synthetic)
  ├→ T010: Implement element-type-specific stiffness assembly
  └→ T011: Validate mixed-element smoother output
  ↓
Phase 6: Polish & Validation (final)
  ├→ T012: [P] Run full test suite on all 4 fixtures
  ├→ T013: Performance benchmarking (target: <2s quad, <3s mixed)
  └→ T014: Document quad stiffness derivation in code comments
```

---

## Parallel Execution

**Phase 3-5 can run in parallel** (after Phase 2 completes):
- Story 1, Story 2, Story 3 are independent
- Suggest: Implement Phase 3 first (simplest, validates foundation), then Phase 4 & 5 in parallel

---

## Phase 1: Setup

No setup tasks required. This is a library modification (no new project structure).

---

## Phase 2: Foundational (Blocking Prerequisites)

These tasks establish shared infrastructure for all user stories. Must complete before any story-specific work.

- [ ] T001 Research Zhou & Shimada triangle formulation and quad extension via analogy in `src/chilmesh/CHILmesh.py` (lines 1319-1373)

- [ ] T002 Design quad stiffness matrices (D_quad, T_quad analogs) as mathematical derivation comment in `src/chilmesh/CHILmesh.py`

- [ ] T003 Implement `_detect_element_type(connectivity_list)` helper function in `src/chilmesh/CHILmesh.py` to identify tri (3 cols) vs quad (4 cols) elements

---

## Phase 3: User Story 1 - Triangle Backward Compat [US1]

**Goal**: Ensure refactored `direct_smoother()` preserves triangle-only behavior  
**Independent Test**: Run `mesh.smooth_mesh('fem')` on annulus, donut, block_o, structured (tri portion) → compare outputs to v0.2.0 baseline  
**Acceptance**: All existing tests pass without modification, output numerical stability maintained

- [ ] T004 [P] [US1] Write regression tests for triangle smoother in `tests/test_smoothing.py` (fixtures: annulus, donut, block_o for triangles; include boundary condition verification)

- [ ] T005 [US1] Refactor `direct_smoother()` to use element-type dispatcher pattern in `src/chilmesh/CHILmesh.py`:
  - Extract triangle-specific stiffness assembly into `_tri_stiffness_assembly()`
  - Create dispatcher that calls `_tri_stiffness_assembly()` for all-triangle meshes
  - Preserve existing numerical behavior exactly (no algorithm changes)

- [ ] T006 [US1] Run regression tests on T004 to verify backward compatibility

---

## Phase 4: User Story 2 - Quad Mesh Support [US2]

**Goal**: Enable FEM smoothing on pure quad meshes  
**Independent Test**: Run `mesh.smooth_mesh('fem')` on structured fixture → verify quad elements smoothed, shapes valid, <2s execution  
**Acceptance**: Quad meshes can be smoothed without manual triangle conversion; quad element validity (no inverted elements); performance <2s

- [ ] T007 [P] [US2] Write tests for quad-only mesh smoother in `tests/test_smoothing.py` (structured fixture; include edge cases: degenerate quads, boundary conditions)

- [ ] T008 [US2] Implement `_quad_stiffness_assembly()` function in `src/chilmesh/CHILmesh.py`:
  - Derive D_quad and T_quad matrices following Zhou & Shimada analogy (documented math in code comment)
  - Assemble global stiffness matrix for quad elements
  - Apply boundary conditions (kinf for boundary nodes)

- [ ] T009 [US2] Update `direct_smoother()` element-type dispatcher to call `_quad_stiffness_assembly()` for quad meshes in `src/chilmesh/CHILmesh.py`

- [ ] T010 [US2] Add quad element validity check after smoothing in `src/chilmesh/CHILmesh.py` (verify no inverted elements, positive area)

- [ ] T011 [US2] Run quad tests on T007 to verify performance targets (<2s) and element validity (95%+)

---

## Phase 5: User Story 3 - Mixed-Element Mesh Support [US3]

**Goal**: Enable FEM smoothing on meshes with both triangles and quads  
**Independent Test**: Create synthetic mixed mesh (combine triangles + quads) → Run `mesh.smooth_mesh('fem')` → verify each element type uses correct stiffness, <3s execution  
**Acceptance**: Mixed meshes can be smoothed in single call; element-type-specific formulations applied; performance <3s; element validity maintained

- [ ] T012 [P] [US3] Write tests for mixed-element mesh smoother in `tests/test_smoothing.py` (synthetic mesh combining annulus triangles + structured quads; include element transition boundaries)

- [ ] T013 [US3] Implement `_mixed_stiffness_assembly()` function in `src/chilmesh/CHILmesh.py`:
  - Detect element types per-element (check connectivity_list row size)
  - Assemble global stiffness matrix with element-type-specific blocks
  - Maintain consistent global numbering for nodes

- [ ] T014 [US3] Update `direct_smoother()` dispatcher to handle mixed meshes (call `_mixed_stiffness_assembly()` when both tri and quad elements present) in `src/chilmesh/CHILmesh.py`

- [ ] T015 [US3] Add mixed-element validity check after smoothing in `src/chilmesh/CHILmesh.py` (verify no inverted triangles or quads, positive areas)

- [ ] T016 [US3] Run mixed-element tests on T012 to verify performance targets (<3s) and element validity (95%+)

---

## Phase 6: Polish & Cross-Cutting Concerns

Final validation, performance tuning, and documentation.

- [ ] T017 [P] Run full test suite (`pytest -v`) on all 4 fixtures in `tests/test_smoothing.py` and legacy test files to verify no regressions

- [ ] T018 [P] Benchmark performance:
  - Measure execution time for quad-only (structured fixture)
  - Measure execution time for mixed-element (synthetic mesh)
  - Document results in `BENCHMARK.md` with target confirmation (<2s quad, <3s mixed)

- [ ] T019 Document quad stiffness matrix derivation as inline comment in `src/chilmesh/CHILmesh.py` near `_quad_stiffness_assembly()` (math, rationale, references to Zhou & Shimada)

- [ ] T020 Update API documentation if needed (`API.md`) to clarify that `smooth_mesh('fem')` now supports all mesh types

---

## Summary

**Total Tasks**: 20 (T001-T020)  
**Tasks per User Story**:
- Setup: 0
- Foundational: 3
- US1 (Triangle): 3
- US2 (Quad): 5
- US3 (Mixed): 5
- Polish: 4

**MVP Scope** (Phases 2-3): 6 tasks (T001-T006)  
**Full Scope** (Phases 2-6): 20 tasks (T001-T020)

**Parallel Opportunities**:
- T004, T007, T012 (test writing) can start immediately after Phase 2
- T017, T018 (validation) can run in parallel in Phase 6
- Phases 3, 4, 5 can run in parallel (stories are independent)

**Independent Test Criteria**:
- **US1**: Triangle regression tests pass, all existing tests pass
- **US2**: Quad tests pass, <2s performance, 95%+ element validity
- **US3**: Mixed-element tests pass, <3s performance, 95%+ element validity

---

**Next Step**: Execute Phase 2 (T001-T003) to establish foundation, then run Phases 3-5 (MVP first: Phase 3, then add Phases 4-5)
