# Feature Specification: Fix Skeletonization Layer Separation Invariant

**Feature Branch**: `006-skeletonize-matlab-port`
**Created**: 2026-05-03
**Status**: Draft
**Input**: User description: "Port the original MATLAB QuADMesh+ meshLayers skeletonization algorithm to Python to fix issue #74 (layer separation invariant violations)."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Researchers see medially-correct layer assignments (Priority: P1)

A computational hydrodynamics researcher loads any 2D triangular or mixed-element mesh into CHILmesh and visualizes its skeletonization layers. They expect each "layer" to represent a concentric ring at increasing depth from the boundary, with the property that no node from a boundary-layer element shares space with elements that are 2+ layers away.

**Why this priority**: Skeletonization is foundational for downstream consumers (MADMESHR, ADMESH-Domains, hydrodynamic models). Incorrect layer assignment cascades into every analysis: medial axis extraction, layer-based smoothing, advancing-front meshing, boundary-layer refinement.

**Independent Test**: Load `chilmesh.examples.annulus()`, call `mesh.layers`, inspect via `mesh.plot_layer()`. Verify that layer color transitions are gradual (Layer 0 → 1 → 2 → ...) with no Layer-N triangle directly touching a Layer-0 triangle for any N ≥ 2.

**Acceptance Scenarios**:

1. **Given** an annulus mesh fixture, **When** I compute `mesh.layers`, **Then** for every pair of layers (k, k+2 or beyond), no element in layer k shares any vertex with any element in layer k+2.
2. **Given** any of the four built-in fixtures (annulus, donut, structured, block_o), **When** I render `mesh.plot_layer()`, **Then** the visualization shows visually-coherent concentric ring structure consistent with the medial-axis property.
3. **Given** the same mesh under both the new Python implementation and the original MATLAB `meshLayers` function, **When** I compare layer counts and OE/IE/OV/IV cardinalities, **Then** they match exactly (modulo 0-vs-1 indexing).

---

### User Story 2 - Downstream consumers continue to function (Priority: P1)

A developer using CHILmesh as a dependency (e.g., MADMESHR, ADMESH-Domains) accesses `mesh.layers["OE"]`, `mesh.layers["IE"]`, `mesh.layers["OV"]`, `mesh.layers["IV"]`, `mesh.layers["bEdgeIDs"]`, and `mesh.n_layers`. They expect the data structure shape and key names to remain identical to v0.2.0.

**Why this priority**: Backward compatibility is mandated by CLAUDE.md ("Public API stable until v1.0").

**Independent Test**: Run the existing CHILmesh test suite (288 tests). All non-skeletonization tests must continue to pass without modification.

**Acceptance Scenarios**:

1. **Given** any pre-existing code that reads `mesh.layers["OE"]`, **When** the new implementation runs, **Then** the returned structure is a list of numpy arrays of element IDs (same type/shape as before).
2. **Given** any pre-existing code that reads `mesh.n_layers`, **When** the new implementation runs, **Then** an integer is returned (numerically may differ from before, but type contract is preserved).
3. **Given** the existing 288-test suite, **When** the fix is applied, **Then** all tests except those that hardcoded the buggy layer count assumption continue to pass.

---

### User Story 3 - Visualization regression check (Priority: P2)

The README's `tests/output/annulus_quickstart.png` shows the warm-start truss optimization pipeline with skeletonization layers in column 2. After the fix, these visualizations must display medially-coherent layer rings.

**Why this priority**: This is a public-facing artifact (rendered on GitHub README). Visual incorrectness undermines user trust.

**Independent Test**: Regenerate `tests/output/annulus_quickstart.png` via `python generate_4row_admesh.py` and visually inspect column 2 of each row.

**Acceptance Scenarios**:

1. **Given** the regenerated visualization, **When** inspecting column 2 of any row, **Then** no high-numbered (e.g., yellow) layer triangles are seen directly adjacent to Layer-0 (purple) boundary triangles.

### Edge Cases

- **Single-element mesh** (n_elems == 1): The single element is layer 0 OE; nLayers = 1.
- **Mesh with internal holes** (annulus, donut): Both outer and inner boundaries seed BFS; layers grow inward from both rings simultaneously.
- **Mixed-element mesh** (triangles + quads with -1 padding): Padding entries (-1) excluded from vertex enumeration.
- **Multiple disconnected components**: Each is skeletonized independently; layer indices are global.
- **Open ocean boundaries** (an ADCIRC concept noted as a future edit in the MATLAB source): Out of scope.

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

### Key Entities *(include if feature involves data)*

- **Layer**: A single iteration's output, comprising OE (outer elements), IE (inner elements), OV (outer vertices), IV (inner vertices), bEdgeIDs (boundary edges that defined this layer's outer frontier).
- **Edge2VertIDs (working copy)**: A copy of `Edge2Vert` mutated during skeletonization (vertex IDs in OV[k] are zeroed out as layers are consumed).
- **Edge2ElemIDs (working copy)**: A copy of `Edge2Elem` mutated during skeletonization (element IDs in OE[k] ∪ IE[k] are zeroed out as layers are consumed).

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Layer separation invariant holds with **0 violations** across all four built-in fixtures (annulus, donut, structured, block_o), measured by a regression test that enumerates all (k, m) layer pairs with |k-m| ≥ 2 and asserts disjoint vertex sets.
- **SC-002**: All 288 currently-passing tests continue to pass (or, if a test asserted a buggy layer count, it is updated to match the MATLAB-correct count and the update is documented in the commit message with a justification).
- **SC-003**: The annulus visualization (`tests/output/annulus_quickstart.png`, column 2 of all 4 rows) shows visually-coherent concentric layer rings.
- **SC-004**: For the annulus mesh, the new Python implementation produces a layer count within ±0 of the MATLAB `meshLayers` function output.
- **SC-005**: Implementation completes skeletonization in less than 5× the runtime of the previous implementation on the block_o fixture (~5,200 elements). Correctness > performance.
- **SC-006**: A new regression test `tests/test_layer_separation_invariant.py` is added that programmatically verifies SC-001 across all fixtures.

## Assumptions

- The MATLAB reference at `domattioli/QuADMesh-MATLAB/blob/master/00_CHILMesh_Class/@CHILmesh/CHILmesh.m`, function `meshLayers`, is the canonical correct algorithm. The Python port should match its behavior bit-for-bit (modulo indexing).
- The "open ocean boundary" branching mentioned in the MATLAB code (commented out) is out of scope; all boundary edges are treated uniformly.
- The `_skeletonize()` method is called once per mesh load (not in a hot path). Memory overhead of mutating copies of `Edge2Vert`/`Edge2Elem` is acceptable.
- Mesh fixtures (annulus, donut, structured, block_o) are representative of all real-world inputs that need to work correctly.
