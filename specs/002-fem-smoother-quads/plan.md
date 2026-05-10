# Implementation Plan: FEM Smoother for Quad & Mixed-Element Meshes

**Branch**: `claude/youthful-goldberg-ueQ9R` | **Date**: 2026-04-27 | **Spec**: [specs/002-fem-smoother-quads/spec.md](spec.md)
**Issue**: GitHub #4

## Summary

Extend `direct_smoother()` method in CHILmesh to support quadrilateral and mixed-element meshes (currently only handles triangles). Implementation uses Zhou & Shimada analogy approach: extend existing triangle stiffness matrices (D and T) to quad geometry, maintaining energy-minimization principle and direct solver pattern. Target: all three mesh types (tri, quad, mixed) smoothed via same `smooth_mesh('fem')` interface with <3s execution on large meshes.

## Technical Context

**Language/Version**: Python 3.10+ (CHILmesh v0.2.0 baseline)
**Primary Dependencies**: numpy, scipy (sparse linear algebra), pytest
**Testing**: pytest with 4 built-in fixtures (annulus, donut, block_o, structured)
**Performance Goals**: <2s for quad-only init, <3s for mixed-element init (SC-002, SC-003)
**Constraints**: Backward compatibility (all existing tests must pass), element validity (95%+, SC-004), numeric stability
**Scale/Scope**: Support meshes up to 100k+ elements

## Constitution Check

✅ **Backward Compatibility**: Public API must be stable until v1.0
  - Current signature: `smooth_mesh(method='fem', acknowledge_change=False, *kwargs)`
  - No breaking changes to this signature
  - Will extend `direct_smoother()` internally; public API unchanged

✅ **Test-Driven Development**: All changes require parametrized tests on all 4 fixtures
  - Tests written before implementation

✅ **Skeletonization Preservation**: Core medial axis algorithm untouched
  - Only modifying `direct_smoother()` method, not mesh topology or layers

✅ **Type Hints Required**: New methods and extended signatures include type hints

**Constitution Status**: ALL GATES PASS

## Project Structure

### Documentation (this feature)

```text
specs/002-fem-smoother-quads/
├── plan.md                      # This file (Phase 1 output)
├── spec.md                       # Feature specification (input)
├── research.md                  # Phase 0 output
├── data-model.md                # Phase 1 output (if needed)
├── quickstart.md                # Phase 1 output (if needed)
├── contracts/                   # Phase 1 output (if applicable)
└── checklists/
    └── requirements.md          # Spec quality checklist
```

### Source Code

```text
src/chilmesh/
├── CHILmesh.py                  # Main mesh class (direct_smoother method here)
└── utils/
    └── plot_utils.py

tests/
├── test_smoothing.py            # NEW: FEM smoother tests
├── conftest.py                  # 4 fixtures: annulus, donut, block_o, structured
└── test_invariants.py           # Existing regression tests
```

**Structure Decision**: Modification localized to `src/chilmesh/CHILmesh.py` (`direct_smoother` method) + new test coverage in `tests/test_smoothing.py`. No structural changes to package layout.

## Phase 0: Research & Clarification

**Gate**: Constitution Check passed

**Research Tasks**:
1. **Quad stiffness formulation**: Clarified — Use Zhou & Shimada analogy (not bilinear isoparametric)
2. **Element connectivity validation**: How to detect and handle degenerate elements (quads with zero area)
3. **Mixed-element FEM assembly**: Verify proper element-type-specific stiffness in global system

**Deliverables**:
- `research.md` with findings on quad FEM formulation, degenerate handling, assembly patterns

## Phase 1: Design & Implementation Planning

**Prerequisites**: research.md complete

**Design Tasks**:
1. **Define quad stiffness matrices**: Create D_quad and T_quad analogs to triangle (D, T) following Zhou & Shimada principle
2. **Element dispatcher**: Modify `direct_smoother()` to:
   - Detect element type from connectivity_list column count (3=tri, 4=quad)
   - Apply correct stiffness based on type
   - Handle mixed-element assembly
3. **Validation logic**: Check element validity post-smoothing (no inversion, positive area)

**Data Model**:
- Element: Vertices indexed into point cloud
- Stiffness matrix: Sparse matrix relating node DOF (2D: x, y per node)
- Boundary detection: Edge-based (from `boundary_edges()`)
- Mixed mesh: Both tri and quad elements in single connectivity list

**Testing Strategy**:
- Regression: Ensure tri-only behavior unchanged
- New coverage: quad-only (structured fixture) and mixed (synthetic test case)
- Performance: <2s quad, <3s mixed on CHILmesh test meshes
