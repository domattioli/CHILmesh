# CHILmesh Data Structure Modernization Plan
**Using Spec-Kit Methodology: Specify → Clarify → Plan → Tasks**

---

## SPECIFICATION: Current State & Constraints

### Current Architecture

CHILmesh: 2D mixed-element mesh library supporting triangular, quadrilateral, and hybrid geometries. Core responsibilities:
- Represent unstructured meshes (points, connectivity)
- Compute mesh properties (quality, angles, area)
- Perform mesh smoothing (FEM-based, angle-based)
- Implement mesh **skeletonization** (medial axis extraction via boundary peeling)
- I/O: `.fort.14` format (ADCIRC hydrodynamic models)
- Visualization via matplotlib

**Current Data Structures:**
1. `connectivity_list` (N×3 or N×4 ndarray): Element→Vertex adjacency
2. `points` (N×3 ndarray): Vertex coordinates (x, y, z)
3. `adjacencies` dict:
   - `Elem2Vert`: Connectivity list reference
   - `Edge2Vert` (N×2 ndarray): Edge→Vertex map
   - `Elem2Edge` (N×M ndarray): Element→Edges
   - `Vert2Edge` (list of lists): Vertex→Edges
   - `Vert2Elem` (list of lists): Vertex→Elements
   - `Edge2Elem` (N×2 ndarray): Edge→[elem1, elem2] with sentinel -1 for boundary
4. `layers` dict: OE, IE, OV, IV, bEdgeIDs (lists per layer)

**Known Performance Issues:**
- `_build_elem2edge()`: O(n²) edge discovery (Block_O takes ~30s)
- Skeletonization: Iterates over remaining elements with set membership checks
- Dynamic operations: Adding/removing edges/nodes requires full rebuild
- Mixed-element handling: Padded-triangle sentinel (vertex3 == vertex0) confusing in 0-indexed Python

**Hard Constraints (Must Preserve):**
1. **Skeletonization output semantics**: Layer decomposition must remain valid (disjoint cover, monotone-shrinking per audit Q3)
2. **API compatibility**: Public method signatures unchanged
3. **.fort.14 roundtrip**: Exact I/O behavior preserved (audit B1)
4. **Mixed-element support**: Triangles + quads in same mesh (audit B4)
5. **Coordinate system**: (x, y, z) with z=0 for 2D meshes

**Soft Constraints:**
- Backwards compatibility for deprecated methods (_mesh_layers)
- No external dependencies beyond scipy, numpy, matplotlib

---

## CLARIFICATION: Downstream & Upstream Needs

### Downstream (Consumers)
1. **MADMESHR**: RL-based advancing-front mesh generator
   - Needs: Efficient element insertion, edge insertion/deletion, domain splitting
   - Current blocker: No API for incremental element/edge addition

2. **ADMESH-Domains**: Mesh registry (public coastal/riverine domain catalog)
   - Needs: Fast mesh loading, metadata queries, batch operations
   - Challenge: May need 100K–1M node meshes

3. **Future: ADMESH** (GitHub 404 — assumed mesh adaptation library)
   - Likely needs: Edge-swapping, node repositioning, refinement/coarsening

### Bridge Opportunities
- CHILmesh as "glue layer": Central mesh representation for MADMESHR→ADMESH→ADMESH-Domains pipeline
- Skeletonization as feature: Pinch-point detection, topology-aware placement

---

## MODERNIZATION GOALS

### Performance
1. **O(n log n) edge discovery**: Replace O(n²) brute-force
2. **Constant-time adjacency lookups**: Graph structure instead of list-of-lists
3. **Efficient skeletonization**: Reduce redundant set operations

### Functionality
1. **Dynamic mesh alterations**: Add nodes, edges, elements; remove nodes
2. **Extended graph traversal**: BFS, DFS, shortest paths, connected-component analysis
3. **Pinch-point detection**: For MADMESHR domain-splitting
4. **Incremental layer updates**: Rebuild only affected layers after changes

### Maintainability
1. **Unified graph representation**: Single source of truth for topology
2. **Type safety**: Dataclass or NamedTuple for edges, elements, vertices
3. **Comprehensive testing**: 100% coverage on graph operations

---

## CANDIDATE DATA STRUCTURES

### Option A: NetworkX Graph
**Verdict:** Good for prototyping. Not ideal for 1M-node production meshes.

### Option B: Custom Compact Graph (Recommended)
```
class MeshGraph:
    vertices: np.ndarray (N×3)        # Coordinates
    edges: List[Tuple[int, int]]      # Unique, unordered edges
    edge2elem: np.ndarray (E×2)       # Edge→[elem_left, elem_right]
    elem2edge: List[List[int]]        # Element→[edge_ids]
    vert2edge: List[List[int]]        # Vertex→[edge_ids]
    vert2elem: List[List[int]]        # Vertex→[element_ids]
    elem2vert: np.ndarray (M×4)       # Element→Vertex (padded to 4)
    elem_type: np.ndarray             # [0/1] for [tri/quad]
```
**Best For:** Production systems needing O(1) lookups and vectorized ops.

### Option C: CSR Sparse Matrix
**Best For:** Physics-based smoothing (FEM), spectral analysis.

### Option D: Half-Edge / Doubly-Linked Representation
**Best For:** Advanced mesh editing tools, topological queries.

---

## RECOMMENDED HYBRID APPROACH

**Phase 1 (Immediate):** Upgrade adjacencies to Option B (Compact Graph)
- Keep `.fort.14` I/O intact
- Replace Elem2Edge, Vert2Elem list-of-lists with efficient dicts/arrays
- Add O(n log n) edge discovery

**Phase 2 (Integration):** Dynamic alteration API
- `add_element(vertices, elem_type)`, `remove_element(elem_id)`
- `add_vertex(x, y, z)`, `remove_vertex(vert_id, strategy='merge')`
- Transactional: Batch operations with single consistency check

**Phase 3 (Algorithms):** Extended graph traversal
- BFS/DFS, connected-component analysis
- Pinch-point detection (MADMESHR bottleneck identification)

**Phase 4 (Downstream Integration):** Bridge to MADMESHR/ADMESH
- Advancing-front API, domain-splitting API, batch import/export

---

## SUCCESS CRITERIA

### Functional
- [ ] All existing tests pass (regression)
- [ ] `_build_elem2edge()` runtime < 1s for Block_O (vs. ~30s current)
- [ ] New tests for dynamic operations
- [ ] Skeletonization produces identical layer decomposition
- [ ] `.fort.14` roundtrip byte-identical

### Performance
- [ ] Edge discovery: O(n log n) instead of O(n²)
- [ ] Adjacency lookup: O(1) instead of O(k)
- [ ] Memory: ≤10% overhead vs. current

### Architectural
- [ ] Single graph data structure (no dual Elem2Edge + Vert2Elem lists)
- [ ] Type hints throughout

### Integration
- [ ] MADMESHR can use CHILmesh for advancing-front placement
- [ ] ADMESH-Domains can bulk-load meshes without slowdowns

---

## RESEARCH ARTIFACTS TO CREATE

1. **Graph Comparison Benchmark** (`research/graph_benchmarks.py`)
2. **Skeletonization Algorithm Analysis** (`research/skeletonization_analysis.md`)
3. **Dynamic Mesh Operation Design** (`research/dynamic_ops_design.md`)
4. **Pinch-Point Detection Algorithm** (`research/pinch_point_detection.md`)
5. **API Design Document** (`research/api_design.md`)

---

## TIMELINE & MILESTONES

| Phase | Effort | Timeline | Blocker |
|-------|--------|----------|---------|
| Research | 2–3 days | Week 1 | None; parallelizable |
| Phase 1 | 3–5 days | Week 2 | Research completion |
| Phase 2 | 2–3 days | Week 3 | Phase 1 completion |
| Phase 3 | 3–5 days | Week 3–4 | Phase 2 completion |
| Phase 4 | 2–3 days | Week 4–5 | Coordination with MADMESHR/ADMESH owners |
| **Total** | ~13–19 days | ~5 weeks | — |

**Gate:** Each phase has full test coverage (unit + integration) before proceeding.

---

## Appendix: Current Code Paths

### Initialization Flow
```
CHILmesh.__init__()
  → _initialize_mesh()
    → _ensure_ccw_orientation()
    → _build_adjacencies()
      → _identify_edges()        [O(n²) brute-force]
      → _build_elem2edge()       [O(n²) brute-force]
      → _build_vert2edge()       [O(n) per vertex]
      → _build_vert2elem()       [O(n) per vertex]
      → _build_edge2elem()       [O(n) per edge]
    → _skeletonize()
      → _mesh_layers()           [deprecated wrapper]
      → [Iterative boundary peeling]
```

### Performance Hotspots
```python
# _identify_edges(): O(n²) edge enumeration
for elem in elements:
    for pair in combinations(vertices, 2):
        edges.add(pair)

# _build_elem2edge(): O(n²) edge matching
for elem in elements:
    for edge in edges:
        if elem has edge:
            elem2edge[elem].append(edge_id)
```

---

**Document Status:** DRAFT - Ready for Specification Review & Clarification
**Last Updated:** 2026-04-26
**Author:** Research Agent (Claude)
