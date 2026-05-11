# Feature Specification: Fix Skeletonization Layer Separation Invariant

**Feature Branch**: `006-skeletonize-matlab-port`
**Created**: 2026-05-03
**Status**: Draft
**Input**: User description: "Port the original MATLAB QuADMesh+ meshLayers skeletonization algorithm to Python to fix issue #74 (layer separation invariant violations)."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Researchers see medially-correct layer assignments (Priority: P1)

Researcher loads any 2D mesh, visualizes skeletonization layers. Expects concentric rings at increasing depth from boundary — no node from boundary-layer element shares space with elements 2+ layers away.

**Why this priority**: Skeletonization foundational for downstream consumers (MADMESHR, ADMESH-Domains, hydrodynamic models). Incorrect layer assignment cascades into every analysis.

**Independent Test**: Load `chilmesh.examples.annulus()`, call `mesh.layers`, inspect via `mesh.plot_layer()`. Verify gradual color transitions (Layer 0 → 1 → 2 → ...) with no Layer-N triangle touching Layer-0 for any N ≥ 2.

**Acceptance Scenarios**:

1. **Given** an annulus mesh fixture, **When** I compute `mesh.layers`, **Then** for every pair of layers (k, k+2 or beyond), no element in layer k shares any vertex with any element in layer k+2.
2. **Given** any of the four built-in fixtures (annulus, donut, structured, block_o), **When** I render `mesh.plot_layer()`, **Then** the visualization shows visually-coherent concentric ring structure consistent with the medial-axis property.
3. **Given** the same mesh under both the new Python implementation and the original MATLAB `meshLayers` function, **When** I compare layer counts and OE/IE/OV/IV cardinalities, **Then** they match exactly (modulo 0-vs-1 indexing).

---

### User Story 2 - Downstream consumers continue to function (Priority: P1)

Developer using CHILmesh as dependency accesses `mesh.layers["OE"]`, `["IE"]`, `["OV"]`, `["IV"]`, `["bEdgeIDs"]`, and `mesh.n_layers`. Expects data structure shape and key names identical to v0.2.0.

**Why this priority**: Backward compatibility mandated by CLAUDE.md ("Public API stable until v1.0").

**Independent Test**: Run existing CHILmesh test suite (288 tests). All non-skeletonization tests pass without modification.

**Acceptance Scenarios**:

1. **Given** any pre-existing code that reads `mesh.layers["OE"]`, **When** the new implementation runs, **Then** the returned structure is a list of numpy arrays of element IDs (same type/shape as before).
2. **Given** any pre-existing code that reads `mesh.n_layers`, **When** the new implementation runs, **Then** an integer is returned (numerically may differ from before, but type contract is preserved).
3. **Given** the existing 288-test suite, **When** the fix is applied, **Then** all tests except those that hardcoded the buggy layer count assumption continue to pass.

---

### User Story 3 - Visualization regression check (Priority: P2)

README's `tests/output/annulus_quickstart.png` shows warm-start pipeline with skeletonization layers in column 2. After fix, must display medially-coherent rings.

**Why this priority**: Public-facing artifact on GitHub README. Visual incorrectness undermines trust.

**Independent Test**: Regenerate `tests/output/annulus_quickstart.png` via `python generate_4row_admesh.py` and visually inspect column 2 of each row.

**Acceptance Scenarios**:

1. **Given** the regenerated visualization, **When** inspecting column 2 of any row, **Then** no high-numbered (e.g., yellow) layer triangles are seen directly adjacent to Layer-0 (purple) boundary triangles.

### Edge Cases

- Single-element mesh: layer 0 OE; nLayers = 1.
- Mesh with internal holes (annulus, donut): both outer and inner boundaries seed BFS simultaneously.
- Mixed-element mesh: padding entries (-1) excluded from vertex enumeration.
- Multiple disconnected components: each skeletonized independently; indices global.
- Open ocean boundaries: out of scope.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The `_skeletonize()` method MUST assign every mesh element to exactly one layer (no duplicates, no omissions).
- **FR-002**: The layer separation invariant MUST hold: for every pair of layers (k, m) with |k - m| ≥ 2, the set of vertices in elements of layer k MUST be disjoint from the set of vertices in elements of layer m.
- **FR-003**: Layer counts MUST match the original MATLAB `meshLayers` function for the same input mesh (within ±0 layers; exact match expected).
- **FR-004**: The `self.layers` dict MUST contain exactly five keys: `OE`, `IE`, `OV`, `IV`, `bEdgeIDs`. Each value MUST be a list of length `n_layers`.
- **FR-005**: `self.n_layers` MUST equal the number of completed iterations of the algorithm.
- **FR-006**: For layer 0 (Python) / iL=1 (MATLAB): `bEdgeIDs[0] = boundary_edges()`, `OV[0] = unique vertices of boundary edges`, `OE[0] = elements adjacent to boundary edges`.
- **FR-007**: For layer k > 0: `bEdgeIDs[k]` MUST be the set of edges whose adjacent-element bookkeeping shows exactly one remaining (unflagged) element neighbor at the start of iteration k.
- **FR-008**: `IE[k]` MUST be computed as: elements still active in the bookkeeping that are adjacent to ANY edge whose `Edge2Vert` references at least one vertex in `OV[k]`. (Critical: this must include elements that share only a vertex — not necessarily an edge — with elements in `OE[k]`.)
- **FR-009**: `IV[k]` MUST be the set of vertices belonging to elements in `OE[k] ∪ IE[k]` minus the vertices in `OV[k]`.
- **FR-010**: After processing layer k, the algorithm MUST flag elements in `OE[k]` and `IE[k]` as "consumed" and flag vertices in `OV[k]` as "consumed" so they are excluded from future iterations.
- **FR-011**: The algorithm MUST terminate when no active elements remain in the bookkeeping.
- **FR-012**: The implementation MUST handle the Python ↔ MATLAB indexing translation correctly: 0-indexed elements/vertices/edges in Python; -1 sentinel for "no neighbor" in `Edge2Elem` (MATLAB uses 0/empty).
- **FR-013**: Padding entries in `connectivity_list` (typically `-1` for triangles padded to 4 columns) MUST be excluded from vertex enumeration in OV/IV computation.

### Key Entities

- **Layer**: Single iteration output — OE (outer elements), IE (inner elements), OV (outer vertices), IV (inner vertices), bEdgeIDs (boundary edges defining outer frontier).
- **Edge2VertIDs (working copy)**: `Edge2Vert` copy mutated during skeletonization (OV[k] vertices zeroed as consumed).
- **Edge2ElemIDs (working copy)**: `Edge2Elem` copy mutated during skeletonization (OE[k] ∪ IE[k] elements zeroed as consumed).

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Layer separation invariant holds with **0 violations** across all four built-in fixtures (annulus, donut, structured, block_o), measured by a regression test that enumerates all (k, m) layer pairs with |k-m| ≥ 2 and asserts disjoint vertex sets.
- **SC-002**: All 288 currently-passing tests continue to pass (or, if a test asserted a buggy layer count, it is updated to match the MATLAB-correct count and the update is documented in the commit message with a justification).
- **SC-003**: The annulus visualization (`tests/output/annulus_quickstart.png`, column 2 of all 4 rows) shows visually-coherent concentric layer rings.
- **SC-004**: For meshes where MATLAB reference output is known, the new Python implementation produces a layer count within ±0 of the MATLAB `meshLayers` function output. **Known reference values** (provided by the project maintainer from external MATLAB runs):
  - Italy domain: **15 layers**
  - Lake Erie domain: **17 layers**
  - Delaware Bay domain: **17 layers**
  - WNAT (Western North Atlantic): **~39 layers** (source variant unconfirmed; one of WNAT_Hagen / WNAT_Onur / WNAT_Test / WNAT_NC_inundation_v6c)
  - Wetting-and-drying test mesh: **15 layers**
  - The annulus, donut, and structured fixtures should also match MATLAB output, though the exact reference layer counts are not yet captured.
- **SC-005**: Implementation completes skeletonization in less than 5× the runtime of the previous implementation on the block_o fixture (~5,200 elements). Correctness > performance.
- **SC-006**: A new regression test `tests/test_skeletonization_invariant.py` is added that programmatically verifies SC-001 across all fixtures. Additionally, `tests/test_skeletonization_matlab_parity_external.py` documents and (when opted in) verifies MATLAB layer-count parity for catalog meshes (Italy=15, Lake Erie=17, Delaware Bay=17, Wetting-Drying=15, WNAT≈39).

## Clarifications

**Q1** (MATLAB equivalence without MATLAB): Use MATLAB source as literal translation reference with line-by-line comments. Verify primarily via layer separation invariant (SC-001). Cross-check against externally-provided ground-truth counts (Italy: 15, Lake Erie: 17) when available.

**Q2** (test failures from layer count changes): Option A — update broken tests to MATLAB-correct values; document in commit message; add `test_matlab_layer_counts` pinning new expected values per fixture.

**Q3** (mixed-element padding): Option A — filter `-1` padding via `v >= 0` check in OV/IV computation. MATLAB logic otherwise unchanged.

## Assumptions

- MATLAB `meshLayers` function is canonical. Python port matches bit-for-bit (modulo indexing).
- Italy and Lake Erie meshes not bundled; reference counts (15 and 17) from external MATLAB runs.
- "Open ocean boundary" branching out of scope; all boundary edges treated uniformly.
- `_skeletonize()` called once per mesh load (not hot path). Mutating copies of `Edge2Vert`/`Edge2Elem` is acceptable overhead.
- Four built-in fixtures representative of real-world inputs.
