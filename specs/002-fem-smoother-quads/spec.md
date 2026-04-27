# Feature Specification: Extend FEM Smoother for Quad & Mixed-Element Meshes

**Feature Branch**: `claude/youthful-goldberg-ueQ9R`  
**Created**: 2026-04-27  
**Status**: Draft  
**Issue**: #4  
**Input**: FEM smoother python implementation does not yet work for quad or mixed-element meshes

## User Scenarios & Testing

### User Story 1 - Direct FEM Smoother Works on Triangle Meshes (Priority: P1)

A mesh researcher uses the direct FEM smoother on existing triangle-only meshes and expects it to continue working without breaking changes.

**Why this priority**: Current functionality must not regress. This is the baseline requirement.

**Independent Test**: Can be tested by running `mesh.smooth_mesh('fem')` on any triangle mesh fixture (annulus, donut, block_o, structured) and verifying smoothed points satisfy FEM energy minimization.

**Acceptance Scenarios**:

1. **Given** a triangle mesh is loaded, **When** `smooth_mesh('fem')` is called, **Then** boundary nodes remain fixed and interior nodes are smoothed according to FEM formulation

2. **Given** a triangle mesh with mixed boundary/interior, **When** smoothing completes, **Then** all boundary edges remain on mesh boundary

---

### User Story 2 - Direct FEM Smoother Works on Quad-Only Meshes (Priority: P1)

A mesh researcher has a pure quad mesh (all elements are quadrilaterals) and needs FEM smoothing without manually converting to triangles.

**Why this priority**: Pure quad meshes are common in structured grid generation and CFD. This is a critical capability gap.

**Independent Test**: Can be tested with the `structured` fixture (pure quad mesh) by calling `mesh.smooth_mesh('fem')` and verifying smoothed points satisfy FEM energy minimization for quad elements.

**Acceptance Scenarios**:

1. **Given** a quadrilateral mesh is loaded, **When** `smooth_mesh('fem')` is called, **Then** interior quad nodes are repositioned via FEM stiffness matrix

2. **Given** a quad mesh, **When** smoothing completes, **Then** quad element shapes improve or remain valid (no inversion)

---

### User Story 3 - Direct FEM Smoother Works on Mixed-Element Meshes (Priority: P1)

A mesh researcher has a mixed mesh (both triangles and quads) and needs FEM smoothing that respects element heterogeneity.

**Why this priority**: Mixed-element meshes are essential for hybrid meshing strategies. CHILmesh supports mixed elements throughout its API; smoother must follow.

**Independent Test**: Can be created by combining triangle and quad elements in a single mesh, then calling `smooth_mesh('fem')` and verifying smoothed points satisfy element-type-specific FEM formulations.

**Acceptance Scenarios**:

1. **Given** a mesh with both triangles and quads, **When** `smooth_mesh('fem')` is called, **Then** each element type uses its respective stiffness formulation

2. **Given** a mixed mesh, **When** smoothing completes, **Then** elements remain valid and interior nodes are repositioned consistently

---

### Edge Cases

- What happens when a mesh contains only boundary nodes (no interior nodes to smooth)?
- How does the smoother handle degenerate quads (e.g., zero area or very thin)?
- How does the smoother handle isolated triangles in a predominantly quad mesh?
- What happens at element transitions (tri-quad boundaries)?

## Requirements

### Functional Requirements

- **FR-001**: System MUST handle triangle-only meshes with existing behavior (no breaking changes)

- **FR-002**: System MUST handle quadrilateral-only meshes by applying quad-specific stiffness matrices

- **FR-003**: System MUST handle mixed-element meshes (triangles + quads) by applying element-type-specific stiffness matrices

- **FR-004**: System MUST identify element types from `connectivity_list` (triangles: 3 columns, quads: 4 columns)

- **FR-005**: System MUST construct FEM stiffness matrices appropriate to each element type:
  - Triangles: Use existing Zhou & Shimada formulation (D matrix and T rotation matrix)
  - Quads: Apply appropriate quad stiffness formulation (e.g., bilinear isoparametric or Balendran quad extension)

- **FR-006**: System MUST apply boundary conditions correctly for all element types (fixed boundary nodes with kinf stiffness)

- **FR-007**: System MUST preserve mesh topology (element-node connectivity) after smoothing

- **FR-008**: System MUST return smoothed point coordinates with preserved z-coordinates

### Key Entities

- **Element**: An individual mesh element (triangle or quad) defined by vertex indices
- **Stiffness Matrix**: FEM-derived matrix relating node displacements to energy minimization
- **Boundary Nodes**: Nodes on mesh boundary edges (marked with large stiffness kinf)
- **Interior Nodes**: Nodes not on boundary, subject to optimization via FEM solve

## Success Criteria

### Measurable Outcomes

- **SC-001**: All existing tests pass without modification (backward compatibility maintained)

- **SC-002**: Direct FEM smoother executes on quad-only mesh (structured fixture) in <2 seconds

- **SC-003**: Direct FEM smoother executes on mixed-element mesh in <3 seconds

- **SC-004**: Smoothed mesh maintains element validity (no inverted/degenerate elements) for 95%+ of test cases

- **SC-005**: Interior node displacements are measurable and consistent with FEM energy minimization principle

## Assumptions

- The existing triangle stiffness formulation (Zhou & Shimada via Balendran 2006) remains correct and unchanged

- Quad stiffness matrices can be derived via analogy to triangle formulation or standard FEM bilinear isoparametric element theory

- Boundary condition enforcement (kinf constraint) applies consistently across all element types

- Mixed-element meshes have contiguous element regions (no interleaving that would require special handling)

- The `connectivity_list` array always has consistent padding: triangles use 3 columns, quads use 4 columns (no padding with -1 or 0)

- Backward compatibility is required: existing code calling `smooth_mesh('fem')` on triangle meshes must work identically

- `scipy.sparse.linalg.spsolve` remains the solver; no performance optimization via direct linear algebra is required in this phase
