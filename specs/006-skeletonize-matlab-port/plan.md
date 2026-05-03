# Implementation Plan: Fix Skeletonization Layer Separation Invariant

**Branch**: `main` (per project policy) | **Date**: 2026-05-03 | **Spec**: [spec.md](./spec.md)
**Input**: Feature specification from `/specs/006-skeletonize-matlab-port/spec.md`

## Summary

Replace the current `_skeletonize()` method in `src/chilmesh/CHILmesh.py` with a faithful Python port of the MATLAB `meshLayers` function from `domattioli/QuADMesh-MATLAB/blob/master/00_CHILMesh_Class/@CHILmesh/CHILmesh.m`. The current Python implementation uses element-element adjacency to compute Inner Elements (IE), which omits elements that share only a vertex (not an edge) with Outer Elements (OE). The MATLAB algorithm computes IE via edges that touch OV vertices — a strictly broader set. This bug fix restores the medial-axis layer separation invariant: a vertex in any element of layer k cannot appear in elements of layer k+2 or beyond.

The fix preserves the existing `self.layers` dict structure (keys: OE, IE, OV, IV, bEdgeIDs) and `self.n_layers` integer for full backward compatibility. Layer counts will change because they were buggy before — these are documented and tests are updated to MATLAB-correct values.

## Technical Context

**Language/Version**: Python 3.10+
**Primary Dependencies**: numpy (already required); no new dependencies
**Storage**: N/A (in-memory mesh state on the `CHILmesh` instance)
**Testing**: pytest (288 existing tests; will add 2 new tests for layer separation invariant and MATLAB-equivalent layer counts)
**Target Platform**: Cross-platform Python library (Linux, macOS, Windows; CI runs on Linux)
**Project Type**: Single-package Python library
**Performance Goals**: Skeletonization runtime < 5× the previous (buggy) implementation on the block_o fixture (~5,200 elements). Correctness is the primary goal.
**Constraints**:
- Public API stable until v1.0 — no changes to `mesh.layers` dict shape or `mesh.n_layers` type contract
- Mixed-element meshes (triangles + quads with `-1` padding) must continue to work
- All non-skeletonization tests must continue to pass
**Scale/Scope**: Largest fixture is block_o with 5,200 elements; production use cases (Italy, Lake Erie) reach 10,000-100,000+ elements

## Constitution Check

CHILmesh has no formal `.specify/memory/constitution.md`. Applicable principles from `.claude/CLAUDE.md`:

- **Backward Compatibility**: Public API stable until v1.0 — ✅ this fix preserves the `mesh.layers` dict shape and `mesh.n_layers` integer contract
- **Test-Driven Development**: ✅ Layer separation invariant test (SC-006) and MATLAB-equivalent layer count test will be added
- **Code Standards**: Type hints, minimal comments (only WHY) — ✅ the algorithm port will include MATLAB-source line references in comments because they document WHY each step exists
- **Branch Policy**: ✅ Working on `main` per active session direction (CLAUDE.md is contradictory between `daily-issue-fixing` and `planning-optimize_modernize`; user has been working on `main` throughout this session)

**Constitutional Status**: PASS — no violations.

## Project Structure

### Documentation (this feature)

```text
specs/006-skeletonize-matlab-port/
├── plan.md              # This file
├── spec.md              # Feature specification
├── research.md          # Phase 0 output (this plan inlines research; no separate file needed)
├── data-model.md        # Phase 1 output (this plan inlines data model; no separate file needed)
├── quickstart.md        # Phase 1 output (this plan inlines quickstart; no separate file needed)
├── contracts/           # N/A — no external interfaces being added or changed
├── checklists/
│   └── requirements.md  # Spec quality checklist (already complete)
└── tasks.md             # Phase 2 output (created by /speckit-tasks)
```

### Source Code (repository root)

```text
src/chilmesh/
├── CHILmesh.py          # MODIFY: replace _skeletonize() with MATLAB port (~100 lines)
├── utils/
│   └── plot_utils.py    # No changes (plot_layer reads self.layers, structure unchanged)
└── __init__.py          # No changes

tests/
├── test_layers_annulus.py        # MODIFY: update test_structured_grid_layers expectation if needed
├── test_layer_separation.py      # NEW: regression test for SC-001 invariant across all fixtures
├── test_matlab_layer_counts.py   # NEW: regression test pinning per-fixture layer counts
└── conftest.py                   # No changes
```

## Phase 0: Research

### MATLAB Reference Algorithm

Source: `domattioli/QuADMesh-MATLAB/blob/master/00_CHILMesh_Class/@CHILmesh/CHILmesh.m`, function `meshLayers` (~lines 1082-1245).

**Algorithm pseudocode** (verbatim translation, see spec FR-001 through FR-013):

```
INPUT: CM (mesh class with adjacencies)
OUTPUT: CM.Layers (struct with OE, IE, OV, IV, bEdgeIDs lists), CM.nLayers (int)

Edge2VertIDs := COPY(CM.edge2Vert)        # Working copy, will be mutated
Edge2ElemIDs := COPY(CM.edge2Elem)        # Working copy, will be mutated

CM.Layers := { OE: [], IE: [], OV: [], IV: [], bEdgeIDs: [] }
iL := 1                                    # Current layer (1-indexed in MATLAB)

WHILE any(Edge2ElemIDs > 0):              # Until all elements consumed
    # Step 1: identify outer boundary edges of layer iL
    IF iL == 1:
        iLbEdgeIDs := CM.boundaryEdges            # Mesh boundary
        OV[iL] := unique(Edge2VertIDs[iLbEdgeIDs, :])
    ELSE:
        # Edges with exactly one remaining active element neighbor
        iLbEdgeIDs := find(sum(Edge2ElemIDs > 0, axis=1) == 1)
        OV[iL] := unique(Edge2VertIDs[iLbEdgeIDs, :])

    bEdgeIDs[iL] := iLbEdgeIDs

    # Step 2: OE = elements adjacent to those boundary edges
    OE[iL] := unique(Edge2ElemIDs[iLbEdgeIDs, :] AND > 0)

    # Step 3: flag OE elements as consumed in Edge2ElemIDs (set to 0)
    Edge2ElemIDs[ismember(Edge2ElemIDs, OE[iL])] := 0

    # Step 4: get all edges that touch ANY OV[iL] vertex
    iLbEdgeIDs2 := find(sum(ismember(Edge2VertIDs, OV[iL]), axis=1) > 0)

    # Step 5: IE = elements still active that are adjacent to those edges
    IE[iL] := unique(Edge2ElemIDs[iLbEdgeIDs2, :] AND > 0)

    # Step 6: flag OV vertices and IE elements as consumed
    Edge2VertIDs[ismember(Edge2VertIDs, OV[iL])] := 0
    Edge2ElemIDs[ismember(Edge2ElemIDs, IE[iL])] := 0

    # Step 7: IV = vertices of (OE ∪ IE) connectivity NOT in OV[iL]
    IV[iL] := setdiff(unique(connectivity_list[OE[iL] ∪ IE[iL], :]), OV[iL])

    iL += 1

CM.nLayers := iL - 1
```

### Python Translation Notes

| MATLAB | Python |
|--------|--------|
| `iL = 1` (1-indexed) | `iL = 0` (0-indexed); `n_layers = iL` after loop |
| `Edge2ElemIDs > 0` (0 = invalid) | `Edge2ElemIDs >= 0` (-1 = invalid sentinel) |
| `Edge2ElemIDs(...) = 0` (consumed) | `Edge2ElemIDs[...] = -1` (consumed) |
| `Edge2VertIDs(...) = 0` (vertex 0 = invalid in MATLAB) | `Edge2VertIDs[...] = -1` (vertex 0 is valid in Python; use -1 sentinel) |
| `unique(...)` | `np.unique(...)` |
| `ismember(A, B)` | `np.isin(A, B)` |
| `setdiff(A, B)` | `np.setdiff1d(A, B)` |
| `sum(X > 0, 2)` (row sums) | `np.sum(X > 0, axis=1)` |

### Critical Bug In Current Python Implementation

Current code (lines 824-829):
```python
inner_elems = []
for elem in outer_elems:
    for neighbor in elem2elem[elem]:        # ← Uses ELEMENT adjacency
        if neighbor in remaining_elements:
            inner_elems.append(neighbor)
```

This computes IE as element-adjacency neighbors of OE. But MATLAB computes IE as elements adjacent to **edges that touch OV vertices**, which is broader: it includes elements that share only a vertex (not an edge) with OE. In a triangular mesh's vertex fan, this distinction matters for elements that touch a vertex at a non-shared edge.

### Mixed-Element Padding (per Q3 Decision)

Where `connectivity_list` may contain `-1` padding for triangles in mixed meshes:
- In `unique(connectivity_list[..., :])` calls, filter out `-1` entries: `vertices = vertices[vertices >= 0]`
- The `Edge2Vert` adjacency does not contain `-1` (edges are always between two real vertices), so no filtering needed there

### Performance Considerations

- The MATLAB algorithm is O(n_elems × n_edges) per iteration in the worst case, with O(n_layers) iterations.
- For the block_o mesh (~5,200 elements), this should complete in well under 1 second.
- The bookkeeping uses two `np.copy()` calls (one each for Edge2Vert and Edge2Elem) — memory is O(n_edges), negligible.

## Phase 1: Design

### Data Model (in-memory only, no persistent storage)

```python
self.layers = {
    "OE": List[np.ndarray[int]],       # OE[k] = element IDs in layer k's outer ring
    "IE": List[np.ndarray[int]],       # IE[k] = element IDs in layer k's inner ring
    "OV": List[np.ndarray[int]],       # OV[k] = vertex IDs on layer k's outer boundary
    "IV": List[np.ndarray[int]],       # IV[k] = vertex IDs on layer k's inner boundary
    "bEdgeIDs": List[np.ndarray[int]], # bEdgeIDs[k] = edge IDs that defined layer k's outer frontier
}
self.n_layers: int                     # Number of completed iterations
```

Layer separation invariant (formally):

```
∀ layers k, m where |k - m| ≥ 2:
    vertices(OE[k] ∪ IE[k]) ∩ vertices(OE[m] ∪ IE[m]) = ∅
```

### Quickstart (manual verification)

```python
import chilmesh
mesh = chilmesh.examples.annulus()
print(f"n_layers: {mesh.n_layers}")
for k in range(mesh.n_layers):
    print(f"Layer {k}: OE={len(mesh.layers['OE'][k])}, IE={len(mesh.layers['IE'][k])}, "
          f"OV={len(mesh.layers['OV'][k])}, IV={len(mesh.layers['IV'][k])}")

# Visual confirmation
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
mesh.plot_layer(ax=ax)
plt.show()
```

Expected: layered concentric rings, no Layer-N color directly touching Layer-0 for N ≥ 2.

### Contracts

No new external interfaces. Internal `_skeletonize()` method signature unchanged: `_skeletonize(self) -> None`.

## Phase 2: Tasks

To be generated by `/speckit-tasks` (see `tasks.md`).

## Risks and Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|-----------|
| Layer count change breaks downstream consumer (MADMESHR, ADMESH-Domains) that hardcoded a specific n_layers | Medium | High | Survey downstream consumers; add deprecation note in CHANGELOG; layer counts are now MATLAB-correct so consumers should reconcile to that |
| Performance regression on block_o (the algorithm now uses np.isin which can be O(n*m)) | Low | Low | Benchmark before/after; use np.unique + sorted lookup if needed |
| Mixed-element mesh (quads) breaks because the algorithm assumes triangles | Medium | Medium | Per Q3 decision: filter -1 padding; test on quad fixtures (structured, quad_2x2) explicitly |
| MATLAB algorithm has an undocumented edge case (e.g., disconnected mesh) that the Python port doesn't handle | Low | Medium | Test on all 4 fixtures; document any deviation in spec assumptions |

## Phase 3: Implementation Order

1. **Add layer separation invariant test first (TDD)**: Create `tests/test_layer_separation.py` with parametrized test across all 4 fixtures. Initially XFAIL (since current code violates it).
2. **Capture current (buggy) layer counts**: Record current `mesh.n_layers` for each fixture so we know what's changing.
3. **Port MATLAB algorithm**: Replace `_skeletonize()` body with the verbatim port; include MATLAB code as comments for traceability.
4. **Verify layer separation invariant passes**: The XFAIL test should now PASS — remove XFAIL marker.
5. **Capture new layer counts**: Run tests, observe new `mesh.n_layers` values, add `tests/test_matlab_layer_counts.py` pinning these values.
6. **Update or delete buggy assertions in `test_layers_annulus.py`**: Per Q2 decision (Option A), update tests with documented justification.
7. **Regenerate visualization**: Run `python generate_4row_admesh.py`; verify column 2 shows clean concentric rings.
8. **Run full test suite**: All 288+ tests should pass.

## Success Validation

Before merge:
- ✅ Layer separation invariant test passes for all 4 fixtures (SC-001, SC-006)
- ✅ All previously passing tests still pass; only tests with documented buggy assumptions are updated (SC-002)
- ✅ Visual inspection of `tests/output/annulus_quickstart.png` shows clean concentric rings (SC-003)
- ✅ Performance on block_o is within 5× of previous runtime (SC-005)
- ⚠️ MATLAB layer count match (SC-004): cannot verify directly (no MATLAB locally); rely on layer separation invariant + Italy/Lake Erie ground truth when those meshes are available
