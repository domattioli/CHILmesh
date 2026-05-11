# CHILmesh Data Structure Modernization Plan
**Phase: Specification & Planning (Spec-Kit Process)**
**Date: 2026-04-26**
**Branch: claude/zen-fermi-NGYbR**
**Status: Planning only - no implementation yet**

---

## 1. SPECIFICATION & CLARIFICATION

### 1.1 Current State Assessment

#### Existing Data Structures
```
CHILmesh Class Core:
├── points: ndarray[Nx3]              # Vertex coordinates (x, y, z)
├── connectivity_list: ndarray[Nx3|4] # Element vertex indices (tri or quad)
├── adjacencies: Dict[str, ...]       # Multiple adjacency representations
│   ├── Elem2Vert: connectivity_list (reference)
│   ├── Edge2Vert: ndarray[Mx2]       # M = number of unique edges
│   ├── Elem2Edge: ndarray[Nx3|4]     # Per-element edge IDs
│   ├── Vert2Edge: List[List[int]]    # Per-vertex incident edge IDs
│   ├── Vert2Elem: List[List[int]]    # Per-vertex incident element IDs
│   └── Edge2Elem: ndarray[Mx2]       # Up to 2 elements per edge (-1 if boundary)
└── layers: Dict[str, List]           # Skeletonization output
    ├── OE: Outer elements per layer
    ├── IE: Inner elements per layer
    ├── OV: Outer vertices per layer
    ├── IV: Inner vertices per layer
    └── bEdgeIDs: Boundary edges per layer
```

#### Critical Algorithms
1. **_skeletonize()**: Medial axis extraction via boundary peeling — iterative layer identification: boundary → outer elements → inner elements; used for mesh decomposition, visualization, mesh adaptation; depends on Edge2Elem, element-to-element connectivity via edges
2. **_build_adjacencies()**: Constructs all adjacency representations
   - `_identify_edges()`: O(n_elems × n_verts_per_elem) set operations
   - `_build_elem2edge()`: **O(n_elems × n_verts × n_edges)** ← BOTTLENECK
   - `_build_vert2edge()`: O(n_edges)
   - `_build_vert2elem()`: O(n_elems × n_verts_per_elem)
   - `_build_edge2elem()`: O(n_elems × n_verts × n_edges) ← BOTTLENECK
3. **Mesh Operations**: `smooth_mesh()`, `elem_quality()`, `interior_angles()`; mixed-element handling adds complexity (padded triangles as 4-col arrays)

#### Known Performance Issues
| Issue | Location | Complexity | Impact |
|-------|----------|-----------|--------|
| Edge search via linear scan | `_build_elem2edge`, `_build_edge2elem` | O(n) per lookup | O(n²) overall |
| Edge identification via set() | `_identify_edges` | O(1) set ops, but iteration | Acceptable |
| No spatial indexing | All operations | N/A | No point location, nearest-neighbor |
| Mixed-element sentinel handling | `_elem_type`, `_ensure_ccw_orientation` | O(n) iteration | Intricate logic, risk of bugs |
| List-of-lists adjacencies | `Vert2Edge`, `Vert2Elem` | O(1) access, O(k) traversal | Inconsistent with numpy arrays |

### 1.2 Downstream Dependencies

- **MADMESHR**: Likely needs mesh traversal, layer info, quality metrics; may need rapid neighbor queries, element/edge insertion/deletion
- **ADMESH**: Historical mesh refinement/coarsening/quality; may need efficient modification, topological queries
- **ADMESH-Domains**: Likely manages multiple domains/boundaries; may need subgraph queries, domain tagging

### 1.3 Functional Requirements

**Must Preserve:** skeletonization, Fort.14 I/O, public API, mixed-element support, plotting

**Must Support:** rapid traversal O(1)/O(log n); efficient node/edge insertion/deletion; bridge interface; fast plot conversion

**Must Improve:** edge building O(n²) → O(n log n) or O(n); code clarity

---

## 2. STATE-OF-THE-ART RESEARCH

### 2.1 Candidate Data Structures

#### A. Half-Edge (DCEL)
O(1) ring/face traversal; standard in CGAL. Higher memory overhead; complex insertion/deletion; not numpy-compatible. **Verdict: ⚠️ Promising but heavy**

#### B. Winged-Edge
Similar to half-edge, no significant advantage. **Verdict: ❌ Less preferable**

#### C. Hybrid: Explicit Hash Map Edge Lookup
Keep numpy arrays; replace linear edge search with dict `(v1,v2)→edge_id`. Minimal code changes; O(1) lookup; numpy-compatible. Does NOT fix list-of-lists. **Verdict: ✅ Quick win**

#### D. Explicit Sparse CSR/CSC Matrices
O(1) lookup; scipy integration; memory efficient. Overkill for small meshes; CSR inefficient for insertion. **Verdict: ⚠️ Good for large meshes**

#### E. NetworkX Graph
Rich graph algorithms; O(1) neighbor queries. Heavy dependency; Python objects slow; not geometry-aware; may not preserve element order. **Verdict: ⚠️ Overkill**

#### F. Refined Numpy + Dict Hybrid
Keep numpy for dense data; sorted arrays + binary search for edges; dict for Vert2Edge/Vert2Elem; CSR for large sparse ops. Scales small→large; clear separation; incremental migration. **Verdict: ✅ Best fit**

### 2.2 Skeletonization-Specific Considerations

```
repeat:
  1. Find boundary edges (edge2elem[:, 1] == -1)
  2. Get boundary vertices
  3. Find outer elements (elements touching boundary edges)
  4. Record OE, OV, bEdgeIDs
  5. Find inner elements (neighbors of outer elements)
  6. Record IE, IV
  7. Compute new boundary (edges between processed/unprocessed)
until no remaining elements
```

Key: doesn't depend on complex topological ops; needs efficient boundary identification + neighbor finding. Any structure supporting O(1) edge→elements sufficient.

### 2.3 Best Practice from Literature

CGAL: DCEL/half-edge. Triangle: edge-based lists + pointers. Gmsh: hybrid (connectivity arrays + lookup maps). VTK: cell/point arrays + incidence structures. DEAL.II: explicit adjacency matrices + sparse storage.

**Key insight:** No universal solution — depends on size (small: dense, large: sparse), operations (read-heavy: optimize lookup, write-heavy: optimize insertion), integration.

---

## 3. STRATEGIC PLAN

### 3.1 Recommended Approach: Layered Modernization

**Phase 1 (Low risk):** Hash map edge lookup O(n)→O(1); refactor `_build_elem2edge` and `_build_edge2elem`; ~2× speedup expected

**Phase 2 (Medium):** Migrate Vert2Edge/Vert2Elem to explicit dicts; update traversal patterns; add type hints + validation

**Phase 3 (Downstream coordination required):** Define CAI (CHILmesh Access Interface); adapter methods for MADMESHR/ADMESH/ADMESH-Domains; integration tests

**Phase 4 (Later if needed):** Spatial indexing (KD-tree); CSR matrices; caching; parallel ops

### 3.2 Design Principles

1. **Backward Compatibility First** — public API unchanged until v1.0; internal refactoring hidden; tests pass without modification
2. **Incremental Migration** — one adjacency at a time; each phase adds methods, keeps old deprecated
3. **Explicit Over Implicit** — replace magic list-of-lists with documented dicts; validation at build time; clear error messages
4. **Cache Locality & CPU Efficiency** — prefer numpy; structured arrays where appropriate; avoid unnecessary copies
5. **Documentation-First Design** — document invariants before code; every adjacency has docstring

### 3.3 Implementation Strategy

#### Data Structure Specification

```python
# NEW: Explicit edge mapping (replaces linear search)
edge_map: Dict[Tuple[int, int], int]
  # Maps (v1, v2) → edge_id (v1 < v2 by convention)
  # Built once in __init__, maintained on insertion/deletion

# UPGRADE: Keep numpy arrays for main adjacencies
Elem2Vert: ndarray[n_elems, 3|4]    # Current connectivity_list
Edge2Vert: ndarray[n_edges, 2]      # Sorted edge vertices
Elem2Edge: ndarray[n_elems, 3|4]    # Element→edge indices

# UPGRADE: Convert to explicit dicts for sparse adjacencies
Vert2Edge: Dict[int, Set[int]]      # vertex → set of edge indices
Vert2Elem: Dict[int, Set[int]]      # vertex → set of element indices

# KEEP: Most reliable parts
Edge2Elem: ndarray[n_edges, 2]      # (elem1, elem2) or elem2=-1 if boundary
layers: Dict[str, List[ndarray]]    # Skeletonization layers (unchanged)
```

#### API Evolution Plan

**Current API (stays valid):**
```python
mesh.boundary_edges()           # returns ndarray
mesh.elem2edge(elem_ids)        # returns ndarray
mesh.vert2elem(vert_id)         # returns list (WILL CHANGE RETURN TYPE)
```

**New API (added alongside old):**
```python
# Explicit accessors returning consistent types
mesh.get_vertex_edges(vert_id: int) → Set[int]      # New: explicit return
mesh.get_vertex_elements(vert_id: int) → Set[int]   # New: explicit return
mesh.find_edge(v1: int, v2: int) → int | None       # New: map lookup
mesh.add_element(vertices: ndarray) → int            # New: mutation
mesh.remove_element(elem_id: int) → bool            # New: mutation
```

---

## 4. DETAILED TASK BREAKDOWN

### 4.1 Phase 1: Hash Map Edge Lookup

1. [ ] **TASK-P1-01**: Create `EdgeMap` class — `add_edge`, `find_edge`, `remove_edge`; tests
2. [ ] **TASK-P1-02**: Refactor `_identify_edges()` to use EdgeMap — return edge list + edge_map
3. [ ] **TASK-P1-03**: Optimize `_build_elem2edge()` O(n²)→O(n log n) — use `edge_map.find_edge()`; benchmark before/after
4. [ ] **TASK-P1-04**: Optimize `_build_edge2elem()` — same approach; benchmark
5. [ ] **TASK-P1-05**: Store EdgeMap in adjacencies dict — `adjacencies['EdgeMap'] = edge_map`; backward compat
6. [ ] **TASK-P1-06**: Performance regression tests — all fixtures faster than baseline; ≥1.5× speedup on large fixtures

### 4.2 Phase 2: Adjacency Modernization

1. [ ] **TASK-P2-01**: Migrate `Vert2Edge` → `Dict[int, Set[int]]`
2. [ ] **TASK-P2-02**: Migrate `Vert2Elem` → `Dict[int, Set[int]]`
3. [ ] **TASK-P2-03**: Add explicit accessors — `get_vertex_edges()`, `get_vertex_elements()`
4. [ ] **TASK-P2-04**: Update internal traversal patterns
5. [ ] **TASK-P2-05**: Type hints + validation for adjacency methods
6. [ ] **TASK-P2-06**: Adjacency structure docs — format, invariants, access patterns

### 4.3 Phase 3: Bridge Infrastructure

1. [ ] **TASK-P3-01**: Define CAI — `docs/CHILmesh_Access_Interface.md`; stable API; backward compat guarantee
2. [ ] **TASK-P3-02**: Create bridge module — `src/chilmesh/bridge.py`; `MeshAdapterForMADMESHR`, etc.
3. [ ] **TASK-P3-03**: Integration tests — simulate MADMESHR/ADMESH/ADMESH-Domains; clear contracts
4. [ ] **TASK-P3-04**: Downstream migration docs — `docs/DOWNSTREAM_MIGRATION_GUIDE.md`; code examples

---

## 5. RISK ASSESSMENT & MITIGATION

| Risk | Impact | Probability | Mitigation |
|------|--------|-------------|-----------|
| Skeletonization breaks | Critical | Low | Retain original _skeletonize, add comprehensive tests |
| Fort.14 I/O breaks | Critical | Low | All existing tests must pass unchanged |
| Performance regression | High | Medium | Benchmark every step, revert if slower |
| API breakage | High | Medium | Deprecation warnings, parallel APIs until v1.0 |
| Memory bloat from extra maps | Medium | Low | Profile memory; optimize if >10% overhead |
| Downstream incompatibility | High | High | Early coordination with MADMESHR/ADMESH authors |

---

## 6. SUCCESS CRITERIA

### Phase 1
- [ ] `_build_elem2edge` and `_build_edge2elem` O(n log n) or O(n)
- [ ] All existing tests pass without modification
- [ ] ≥1.5× speedup on block_o fixture
- [ ] No new dependencies

### Phase 2
- [ ] Vert2Edge and Vert2Elem as explicit dicts
- [ ] Traversal patterns updated; type hints complete
- [ ] Docs updated with invariants

### Phase 3
- [ ] CAI published; bridge adapters for all three projects
- [ ] Integration tests passing; no hidden dependencies

### Overall
- [ ] 100% test pass rate; no breaking API changes
- [ ] ≥1.5× improvement on large meshes
- [ ] Code ready for v0.2.0

---

## 7. ESTIMATED EFFORT

| Phase | Tasks | Est. Effort | Risk Level |
|-------|-------|------------|-----------|
| 1: Edge Mapping | 6 | 8-12 hours | Low |
| 2: Adjacency Modernization | 6 | 12-16 hours | Medium |
| 3: Bridge Infrastructure | 4 | 20-24 hours | High |
| **Total** | **16** | **40-52 hours** | **Medium** |

---

## 8. NEXT STEPS

1. Stakeholder buy-in — share with MADMESHR/ADMESH/ADMESH-Domains authors; identify gaps; adjust timeline
2. Create GitHub Issues — one per task (16 total); link to this doc; assign milestones
3. Branch structure — `planning-optimize_modernize` for planning; `feature/phase1-edge-mapping` for Phase 1 work
4. Begin Phase 1 — start TASK-P1-01 (EdgeMap class); benchmark as you go

---

## 9. RELATED DOCUMENTS

- `README.md`: CHILmesh overview
- `PROGRESS.md`: Historical notes on 0.1.1 fixes
- `CHANGELOG.md`: Release history
- `src/chilmesh/CHILmesh.py`: Current implementation
- `tests/`: Existing test suite

---

## 10. GLOSSARY

| Term | Definition |
|------|-----------|
| DCEL | Doubly-Connected Edge List (half-edge data structure) |
| CAI | CHILmesh Access Interface (stable API for downstream projects) |
| Skeletonization | Medial axis extraction via boundary peeling |
| Adjacency | Topological relationship (e.g., Elem2Vert, Vert2Edge) |
| O(n) | Linear time complexity in n |
| O(n²) | Quadratic time complexity in n |
| O(n log n) | Linearithmic time complexity (merge sort, binary search) |

---

**Document Status: DRAFT - PLANNING PHASE**
**Next Review: After Phase 1 complete**
