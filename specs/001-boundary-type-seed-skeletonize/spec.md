# Feature Specification: Boundary-Type Seeding for Skeletonization

**Feature Branch**: `001-boundary-type-seed-skeletonize`
**Created**: 2026-06-02
**Status**: Draft
**Input**: Issue #129 — layers/skeletonization should consider boundary types

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Filter skeletonization to land boundaries only (Priority: P1)

A mesh user has a coastal domain (e.g. WNAT) with both open-ocean boundaries and land (mainland/island) boundaries. They want skeletonization to peel only from land boundaries so that layer 0 starts at the coastline, not at the open-ocean edge. This lets downstream tools (QuADMesh) place quads in the open-ocean region.

**Why this priority**: Core use case from issue #129. Without this, WNAT-style workflows must post-process layers manually.

**Independent Test**: Construct a mesh with two distinct boundary segments (one `kind=open`, one `kind=flow`); pass `seed_boundary_kinds=['flow']`; verify layer 0 OV only contains nodes from the flow segment.

**Acceptance Scenarios**:

1. **Given** a mesh with `boundary_segments=[{kind:'open',...}, {kind:'flow',...}]`, **When** `_skeletonize(seed_boundary_kinds=['flow'])` is called, **Then** layer 0 seeds only from nodes on the `flow` segment.
2. **Given** the same mesh, **When** `seed_boundary_kinds=None` (default), **Then** behavior is identical to current (all boundary edges seed layer 0) — backward compatible.
3. **Given** a mesh with no `boundary_segments`, **When** any `seed_boundary_kinds` is provided, **Then** all boundary edges are used (graceful fallback) and a warning is emitted.

---

### User Story 2 — Filter by IBTYPE value (Priority: P2)

A user wants to seed only from IBTYPE=0 (mainland) boundaries and exclude island (IBTYPE=1–3) and weir boundaries. They pass `seed_ibtypes=[0]` and get layers that propagate from mainland only.

**Why this priority**: More granular than kind-based filtering; needed for complex ADCIRC meshes.

**Independent Test**: Construct mesh with segments of IBTYPE=0 and IBTYPE=1; pass `seed_ibtypes=[0]`; verify layer 0 nodes match only IBTYPE=0 segment.

**Acceptance Scenarios**:

1. **Given** a mesh with segments of mixed IBTYPE, **When** `seed_ibtypes=[0]` is passed, **Then** only IBTYPE=0 edges seed layer 0.
2. **Given** `seed_ibtypes` and `seed_boundary_kinds` both provided, **Then** the filter is the intersection (both must match).

---

### User Story 3 — Public API passthrough (Priority: P3)

The filtering parameters must be reachable through the public `CHILmesh.__init__()` and via the `full_init` / `fast_init` paths so users do not need to call private methods.

**Why this priority**: Without a public entry point the feature is unusable in practice.

**Independent Test**: Instantiate `CHILmesh(connectivity, points, seed_boundary_kinds=['flow'])`; verify `layers['OV'][0]` contains only flow-segment nodes.

**Acceptance Scenarios**:

1. **Given** `CHILmesh(..., seed_boundary_kinds=['flow'])`, **Then** skeletonization is seeded correctly from construction.
2. **Given** `compute_layers=False`, **Then** `seed_boundary_kinds` is accepted but silently ignored (no error).

---

### Edge Cases

- What if `seed_boundary_kinds` names a kind that no segment has? → Layer 0 seed set is empty → raise `ValueError` with clear message listing available kinds.
- What if filtered seed nodes do not cover any boundary edges? → Same `ValueError` path.
- What if `boundary_segments` is an empty list and filter is specified? → Graceful fallback: use all boundary edges + emit warning.
- What if mesh has no fort.14 boundaries at all (node+element-only mesh)? → Filter is ignored; all boundary edges used; no error.
- Mixed filter: `seed_boundary_kinds=['flow']` + `seed_ibtypes=[0, 1]` → intersection semantics: flow AND (ibtype 0 or 1).

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: `_skeletonize()` MUST accept `seed_boundary_kinds: list[str] | None = None` and `seed_ibtypes: list[int] | None = None` parameters.
- **FR-002**: When both parameters are `None`, behavior MUST be identical to the current implementation (all boundary edges seed layer 0).
- **FR-003**: When `seed_boundary_kinds` is provided, layer 0 MUST seed only from boundary edges whose nodes belong to segments with matching `kind`.
- **FR-004**: When `seed_ibtypes` is provided, layer 0 MUST seed only from boundary edges whose nodes belong to segments with matching `ibtype`.
- **FR-005**: When both are provided, MUST apply intersection (kind AND ibtype must both match).
- **FR-006**: MUST raise `ValueError` when the filter resolves to zero seed nodes AND `boundary_segments` is non-empty.
- **FR-007**: MUST emit a `warnings.warn` and fall back to all boundary edges when `boundary_segments` is empty and a filter is specified.
- **FR-008**: Parameters MUST be accessible via the public `CHILmesh.__init__()` constructor and forwarded through `_initialize_mesh()`.
- **FR-009**: MATLAB parity tests MUST still pass with default (`None`) parameters.

### Key Entities

- **BoundarySegment**: `{kind: "open"|"flow", ibtype: int|None, nodes: ndarray[int]}` — 0-based node indices.
- **SeedNodeSet**: the union of `nodes` arrays from all matching boundary segments; used to filter boundary edges for layer 0.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: All existing skeletonization tests pass unchanged with default parameters.
- **SC-002**: A WNAT-style mesh with mixed boundary kinds produces disjoint layer-0 node sets when filtered vs. unfiltered.
- **SC-003**: New tests covering kind-filter, ibtype-filter, combined filter, fallback, and error paths all pass.
- **SC-004**: No measurable change in skeletonization runtime for the default (unfiltered) path.

## Assumptions

- `boundary_segments` is populated only for meshes loaded from fort.14; node+element-only meshes have an empty list — filter gracefully no-ops.
- Edge semantics: a boundary edge "belongs" to a segment if BOTH its vertices are in that segment's `nodes` array; this matches the topological intent (edges are defined by consecutive node pairs on a boundary segment).
- `seed_boundary_kinds` values are case-sensitive and match the stored `kind` values (`"open"`, `"flow"`) exactly.
- The intersection of `seed_boundary_kinds` and `seed_ibtypes` is computed at the segment level before deriving the seed node set.
- MATLAB oracle tests (`test_skeletonization_matlab_parity_external.py`) use fort.14 meshes with `boundary_segments`; default parameters ensure parity is preserved.
