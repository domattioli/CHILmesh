# Phase 5: Mesh Mutation API

**Goal:** Add topology-modifying operations (split, merge, swap, insert, remove) to enable adaptive meshing  
**Scope:** Design spec for API, algorithms, invariants  
**Impact:** Unblocks MADMESHR adaptive refinement (~50% of adaptive use cases)  
**Effort:** Design: 2 hours. Implementation: 8-12 hours (Phase 5 execution)

---

## Problem

CHILmesh is read-only post-init. Adaptive workflows require full rebuild:
- Split triangle (refinement) → create new CHILmesh from scratch (~3s per edit on 100k mesh)
- 1000 adaptive iterations × 3s = 50 minutes overhead just rebuilding

**Standard mesh libraries** (Triangle, meshpy, gmsh, pyvista) support on-the-fly topology edits.

---

## Core Operations

### 1. split_triangle(elem_id: int, point: Optional[ndarray]) → ndarray
**Semantics:** Triangle → 4 smaller triangles (1-to-4 subdivision)

**Input:**
- `elem_id`: Triangle to split
- `point`: Optional point inside; if None, use barycenter

**Output:** Array of 4 new element IDs

**Use case:** Refinement near features (shorelines, steep gradients)

```python
m = CHILmesh(...)
new_elem_ids = m.split_triangle(elem_id=100, point=np.array([0.5, 0.5]))
# new_elem_ids = [401, 402, 403, 404]
assert m.n_elems == original_count + 3  # 3 new elements added
```

### 2. swap_edge(edge_id: int) → tuple[int, int]
**Semantics:** Flip edge diagonal in quad (two adjacent tris)

**Input:**
- `edge_id`: Interior edge shared by exactly 2 triangles

**Output:** Tuple of 2 new element IDs (the two post-swap triangles)

**Precondition:** Edge must be interior (not boundary)

**Use case:** Quality improvement (Lawson flip, anti-clockwise restoration)

```python
m = CHILmesh(...)
new_elem_ids = m.swap_edge(edge_id=150)
assert m.n_elems == original_count  # Element count unchanged
```

### 3. merge_elements(elem_a: int, elem_b: int) → int
**Semantics:** Merge two adjacent triangles → 1 quad or 2 tris (depends on configuration)

**Input:**
- `elem_a`, `elem_b`: Adjacent elements

**Output:** Single element ID (merged result) or raises if non-adjacent

**Use case:** Coarsening in low-gradient regions (deep water)

**Precondition:** Elements must share exactly one edge

```python
m = CHILmesh(...)
merged_id = m.merge_elements(elem_a=10, elem_b=11)
assert m.n_elems < original_count
```

### 4. insert_vertex(point: ndarray) → int
**Semantics:** Add vertex to mesh; automatically triangulates containing element or edge

**Input:**
- `point`: 2D coordinate (x, y)

**Output:** New vertex ID

**Precondition:** Point must be inside or on boundary of mesh

**Algorithms:** Bowyer-Watson for interior, edge-split for boundary

**Use case:** Delaunay refinement, incremental mesh improvement

```python
m = CHILmesh(...)
new_vert_id = m.insert_vertex(point=np.array([0.3, 0.7]))
assert m.n_verts == original_verts + 1
```

### 5. remove_vertex(vert_id: int) → None
**Semantics:** Delete vertex; adjacent elements re-triangulated

**Input:**
- `vert_id`: Vertex to remove

**Precondition:** Vertex must not be on boundary (or special handling needed)

**Use case:** Degenerate vertex elimination, coarsening

```python
m = CHILmesh(...)
m.remove_vertex(vert_id=42)
assert m.n_verts < original_verts
```

### 6. collapse_edge(edge_id: int) → int
**Semantics:** Merge two vertices via edge collapse; one survives, elements updated

**Input:**
- `edge_id`: Edge to collapse

**Output:** Surviving vertex ID

**Use case:** Aggressive coarsening, degenerate element removal

**Precondition:** Collapse must not violate element orientation invariants

```python
m = CHILmesh(...)
surviving_vert_id = m.collapse_edge(edge_id=123)
assert m.n_verts < original_verts
```

---

## Bulk Operations (Transactional)

### smooth_topology(quality_threshold: float)
**Semantics:** Automatically swap low-quality edges until all meet threshold

**Input:**
- `quality_threshold`: Min acceptable quality (e.g., 0.3)

**Output:** Count of swaps performed

**Algorithm:** Iterate all interior edges, swap if quality below threshold

**Use case:** Post-refinement quality recovery

```python
m = CHILmesh(...)
swaps = m.smooth_topology(quality_threshold=0.3)
print(f"Performed {swaps} edge swaps")
assert all(q >= 0.3 for q in m.elem_quality())
```

---

## Invariants (Must Preserve)

✅ **No duplicate edges** (canonical (min, max) form)  
✅ **Adjacency consistency** (`edge2elem`, `vert2edge`, `vert2elem` synchronized)  
✅ **Boundary markers** (boundary edges remain -1; interior edges have 2 neighbors)  
✅ **Element orientation** (all positive area, counterclockwise)  
✅ **Mixed-element padding** (triangles padded to `[v0, v1, v2, v0]`)  
✅ **Layer invalidation** (skeletonization invalidated on topology change; must recompute)  
✅ **Coordinate continuity** (no changes to existing vertex coords)

---

## Implementation Strategy

### Data Structure Changes

**None required.** Existing `Edge2Elem`, `Elem2Edge`, adjacency dicts sufficient.

**Transaction Pattern:**
```python
def split_triangle(self, elem_id: int, point: Optional[ndarray]) -> ndarray:
    # Snapshot adjacency state (for rollback)
    snapshot = self._snapshot_adjacencies()
    
    try:
        # Perform edits
        new_elem_ids = self._split_triangle_impl(elem_id, point)
        # Validate invariants
        self._validate_adjacencies()
        # Invalidate layers (skeletonization stale)
        self.layers = {}
        return new_elem_ids
    except Exception as e:
        # Rollback on failure
        self._restore_adjacencies(snapshot)
        raise
```

### Performance

**Target:** ≤ 100 μs per single operation (split, swap, merge)

- Hash lookups: O(1) via `EdgeMap`
- Adjacency updates: O(degree of affected vertex) ≈ O(6-8)
- Total: O(1) amortized

---

## API Design

### Current (Read-Only)
```python
from chilmesh import CHILmesh
m = CHILmesh(connectivity, points)
# Can only read + smooth (coords only)
m.smooth_mesh(iterations=10)
```

### Post-Phase 5 (Mutable)
```python
from chilmesh import CHILmesh, MutableMesh

# Explicit opt-in to mutability
m = MutableMesh(connectivity, points)

# Topology edits
new_elems = m.split_triangle(elem_id=10)
m.swap_edge(edge_id=42)
m.merge_elements(elem_a=5, elem_b=6)

# Save edited mesh
m.write_to_fort14("refined.14")
```

**Design rationale:** Mutability is explicit. Immutable CHILmesh remains default for stability.

---

## Blockers & Concerns

### 1. Mixed-Element Mutation
**Concern:** How to split a quad? (Quad → 4 quads? Or 4 tris?)

**Design decision (TBD in Phase 5):** 
- Quad → 4 quads (paving)
- OR Quad → 4 tris (simpler)
- Recommend: paving (preserve mesh regularity)

### 2. Boundary Vertex Constraints
**Concern:** Vertices on boundary have rigid position constraints (coastline)

**Design decision:** Remove boundary vertices only via collapse (safe). insert_vertex on boundary re-meshes safely.

### 3. Layer Recomputation Cost
**Concern:** Every edit invalidates skeletonization; recompute is ~3s

**Design decision:** Accept cost for now. Phase 5.2 (incremental skeletonization) addresses this.

### 4. API Stability
**Concern:** Mutation API = major version change. Breaking API risk.

**Design decision:** MutableMesh = new class (no breaking changes to CHILmesh). Users opt-in.

---

## Testing Strategy

### Unit Tests
```python
def test_split_triangle_creates_four():
    """split_triangle → 4 sub-triangles"""
    m = MutableMesh(...)
    original_elems = m.n_elems
    new_ids = m.split_triangle(elem_id=0)
    assert len(new_ids) == 4
    assert m.n_elems == original_elems + 3

def test_split_preserves_invariants():
    """After split, all invariants hold"""
    m = MutableMesh(...)
    m.split_triangle(elem_id=5)
    assert_no_duplicate_edges(m)
    assert_all_elements_ccw(m)
    assert_adjacency_consistent(m)
```

### Regression Tests
- Compare quality metrics before/after (should improve or stay same)
- Verify 1000 random splits + swaps on block_o don't crash
- Determinism: same seed → same result

---

## Success Criteria

- [ ] API surface fully defined (this spec)
- [ ] 6 single operations implemented + tested
- [ ] 1 bulk operation (smooth_topology) implemented
- [ ] All invariants preserved (verified via assertions)
- [ ] Performance target met: ≤100 μs/op on 100k mesh
- [ ] MutableMesh class non-breaking to CHILmesh
- [ ] Documentation: usage examples, invariants, performance trade-offs

---

## Phase 5 Timeline (Rough)

- **Phase 5.0 (Mutation API):** This spec + implementation (8-12 hours)
- **Phase 5.1 (Mutation Tests):** Regression suite (4-6 hours)
- **Phase 5.2 (Incremental Skeletonization):** Local layer updates (6-8 hours)
- **Phase 5.3 (Spatial Indexing):** Point-location queries (4-6 hours)

**v2.0 release candidate:** After Phases 5.0-5.1 complete (~1-2 weeks)

---

## Related

- Issue #94: Mesh mutation API (this spec)
- Issue #93: Incremental skeletonization (depends on #94)
- Issue #92: Spatial indexing (complements #94)
- MADMESHR: Primary consumer (adaptive refinement)
- Constitution Principle II: Immutability by default (must be opt-in)
