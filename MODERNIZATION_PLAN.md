# CHILmesh Data Structure Modernization Plan
**Using Spec-Kit Methodology: Specify → Clarify → Plan → Tasks**

---

## SPECIFICATION: Current State & Constraints

### Current Architecture

**CHILmesh** is a 2D mixed-element mesh library supporting triangular, quadrilateral, and hybrid geometries. Core responsibilities:
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
4. `layers` dict (skeletonization result):
   - `OE`: List of arrays (outer elements per layer)
   - `IE`: List of arrays (inner elements per layer)
   - `OV`: List of arrays (outer vertices per layer)
   - `IV`: List of arrays (inner vertices per layer)
   - `bEdgeIDs`: List of arrays (boundary edge IDs per layer)

**Known Performance Issues:**
- `_build_elem2edge()`: O(n²) edge discovery (audit note: Block_O takes ~30s)
- Skeletonization: Iterates over remaining elements with set membership checks
- Layer structure: Multiple separate lists → cache misses, consistency burden
- Dynamic operations: Adding/removing edges/nodes requires full rebuild
- Mixed-element handling: Padded-triangle sentinel (vertex3 == vertex0) is confusing in 0-indexed Python

**Hard Constraints (Must Preserve):**
1. **Skeletonization output semantics**: Layer decomposition must remain valid (disjoint cover, monotone-shrinking layer sizes per audit Q3)
2. **API compatibility**: Public methods like `signed_area()`, `interior_angles()`, `elem_quality()`, etc. must maintain signatures
3. **.fort.14 roundtrip**: Must preserve exact I/O behavior (audit B1)
4. **Mixed-element support**: Triangles + quads in same mesh (audit B4)
5. **Coordinate system**: (x, y, z) with z=0 for 2D meshes

**Soft Constraints (Nice to Have):**
- Backwards compatibility for deprecated methods (_mesh_layers)
- Plotting via matplotlib (plot utilities)
- No external dependencies beyond scipy, numpy, matplotlib

---

## CLARIFICATION: Downstream & Upstream Needs

### Upstream (Dependencies)
- **scipy.sparse**: LIL matrix for FEM smoother
- **scipy.spatial.Delaunay**: Random mesh generation
- **numpy**: All array operations
- **matplotlib**: Visualization

### Downstream (Consumers)
1. **MADMESHR**: Reinforcement learning–based advancing-front mesh generator
   - Needs: Efficient element insertion, edge insertion/deletion, domain splitting
   - Operations: Continuous placement of quads/triangles at advancing front
   - Challenge: Dynamic mesh growth during element placement
   - Current blocker: No API for incremental element/edge addition

2. **ADMESH-Domains**: Mesh registry (public coastal/riverine domain catalog)
   - Needs: Fast mesh loading, metadata queries, batch operations
   - Operations: Bulk load, filter by quality/size, export subsets
   - Challenge: May need to handle 100K–1M node meshes

3. **Future: ADMESH** (GitHub 404 — assumed mesh adaptation library)
   - Likely needs: Edge-swapping, node repositioning, refinement/coarsening
   - Challenge: Requires efficient point location and edge traversal

### Bridge Opportunities
- **CHILmesh as "glue layer"**: Central mesh representation for MADMESHR→ADMESH→ADMESH-Domains pipeline
- **Skeletonization as feature**: Could be valuable for MADMESH domain analysis (pinch-point detection, topology-aware placement)

---

## MODERNIZATION GOALS

### Performance
1. **O(n log n) edge discovery**: Replace O(n²) brute-force with spatial index or sorted iteration
2. **Constant-time adjacency lookups**: Use graph structure instead of list-of-lists
3. **Efficient skeletonization**: Reduce redundant set operations, use priority queue for layer frontier

### Functionality
1. **Dynamic mesh alterations**: Add nodes, edges, elements; remove nodes (preserving topology)
2. **Extended graph traversal**: BFS, DFS, shortest paths, connected-component analysis
3. **Pinch-point detection**: For MADMESHR domain-splitting (identify narrow bottlenecks)
4. **Incremental layer updates**: Rebuild only affected layers after node insertion/deletion

### Maintainability
1. **Unified graph representation**: Single source of truth for topology (avoid Elem2Edge, Vert2Elem dual lists)
2. **Type safety**: Use dataclass or NamedTuple for edges, elements, vertices
3. **Clear separation of concerns**: Topology layer vs. geometry layer vs. algorithm layer
4. **Comprehensive testing**: 100% coverage on graph operations before shipping

---

## CANDIDATE DATA STRUCTURES & RESEARCH SUMMARY

### Option A: NetworkX Graph (Reference Implementation)
**Pros:**
- Industry-standard, battle-tested, excellent docs
- BFS/DFS/shortest-paths built-in
- Node/edge attribute storage
- Community support

**Cons:**
- Python-only (no NumPy vectorization)
- Memory overhead for general graphs
- Slower than optimized sparse matrices for large meshes

**Verdict:** Good for prototyping, validation, algorithm development. Not ideal for 1M-node production meshes.

---

### Option B: Custom Compact Graph (Recommended)
**Structure:**
```
class MeshGraph:
    vertices: np.ndarray (N×3)        # Coordinates
    edges: List[Tuple[int, int]]      # Unique, unordered edges (edge list)
    edge2elem: np.ndarray (E×2)       # Edge→[elem_left, elem_right]
    elem2edge: List[List[int]]        # Element→[edge_ids]
    vert2edge: List[List[int]]        # Vertex→[edge_ids]
    vert2elem: List[List[int]]        # Vertex→[element_ids]
    elem2vert: np.ndarray (M×4)       # Element→Vertex (padded to 4)
    elem_type: np.ndarray             # [0/1] for [tri/quad]
```

**Pros:**
- Explicit, minimal, NumPy-friendly
- Fast node/edge iteration
- Direct mapping for spatial algorithms
- Cache-friendly for linear traversals

**Cons:**
- Must maintain invariants manually (consistency)
- No built-in traversal algorithms (must code BFS/DFS)
- Requires careful synchronization when adding/removing elements

**Best For:** Production systems needing O(1) lookups and vectorized ops.

---

### Option C: CSR Sparse Matrix (Adjacency Matrix)
**Structure:**
```
adj_matrix: scipy.sparse.csr_matrix (N×N)  # Vertex×Vertex adjacency
elem_data: Separate storage for element properties
```

**Pros:**
- Excellent for spectral methods, Laplacian matrices
- Linear algebra library support (sparse.linalg)
- Built-in matrix operations (products, solves)

**Cons:**
- Poor for mixed-element geometries (quads vs. triangles)
- Harder to track edge multiplicity (if any)
- Sparse matrix construction is O(n log n) (sorting overhead)

**Best For:** Physics-based smoothing (FEM), spectral analysis.

---

### Option D: Half-Edge / Doubly-Linked Representation (Academic Gold Standard)
**Structure:** Each edge has two "half-edges" (directed), linked in circular chains around vertices.

**Pros:**
- Fast local queries (neighbors of vertex in CCW order)
- Elegant for manifold topology operations
- Standard in computational geometry

**Cons:**
- Complex to implement correctly (off-by-one errors common)
- Requires special handling for non-manifold boundaries (mixed triangles/quads)
- Overkill for static meshes; mainly useful for dynamic topology operations

**Best For:** Advanced mesh editing tools, topological queries.

---

## RECOMMENDED HYBRID APPROACH

**Phase 1 (Immediate):** Upgrade adjacencies to **Option B (Compact Graph)**
- Keep existing `.fort.14` I/O intact
- Incrementally replace Elem2Edge, Vert2Elem list-of-lists with efficient dicts/arrays
- Add O(n log n) edge discovery
- Benchmarks: Compare against current implementation

**Phase 2 (Integration):** API for dynamic alterations
- `add_element(vertices, elem_type)`: Insert element, update edges
- `remove_element(elem_id)`: Delete element, clean up orphaned edges
- `add_vertex(x, y, z)`: Insert node, update vertex lists
- `remove_vertex(vert_id, strategy='merge')`: Delete node (merge strategy TBD)
- Transactional: Batch operations with single consistency check

**Phase 3 (Algorithms):** Extended graph traversal
- BFS/DFS from any starting vertex/element
- Connected-component analysis
- Pinch-point detection (bottleneck identification for MADMESHR)
- Incremental skeletonization (update layers after changes, not full rebuild)

**Phase 4 (Downstream Integration):** Bridge to MADMESHR/ADMESH
- Advancing-front API (expose frontier edges, support element placement)
- Domain-splitting API (detect pinch points, split at narrowest edge)
- Batch import/export for ADMESH-Domains registry

---

## SUCCESS CRITERIA

### Functional
- [ ] All existing tests pass (regression)
- [ ] `_build_elem2edge()` runtime < 1s for Block_O (vs. ~30s current)
- [ ] New tests for dynamic operations (add/remove element, vertex)
- [ ] Skeletonization produces identical layer decomposition
- [ ] `.fort.14` roundtrip byte-identical

### Performance
- [ ] Edge discovery: O(n log n) instead of O(n²)
- [ ] Adjacency lookup: O(1) instead of O(k) where k = degree
- [ ] Memory: ≤10% overhead vs. current implementation

### Architectural
- [ ] Single graph data structure (no dual Elem2Edge + Vert2Elem lists)
- [ ] Clear public/private separation (no leaking implementation details)
- [ ] Comprehensive docstrings (algorithm descriptions, complexity, examples)
- [ ] Type hints throughout (no `Any` without justification)

### Integration
- [ ] MADMESHR can use CHILmesh for advancing-front placement
- [ ] ADMESH-Domains can bulk-load meshes without slowdowns
- [ ] ADMESH can use pinch-point detection API (when available)

---

## RESEARCH ARTIFACTS TO CREATE

Before implementation:

1. **Graph Comparison Benchmark** (`research/graph_benchmarks.py`)
   - Implement all four options (NetworkX, Compact, CSR, Half-Edge)
   - Measure: construction time, adjacency lookup, traversal speed
   - Test on annulus, structured, block_o, and synthetic 1M-node mesh
   - Generate comparison table

2. **Skeletonization Algorithm Analysis** (`research/skeletonization_analysis.md`)
   - Current implementation: describe algorithm (BFS layer peeling?)
   - Complexity analysis: steps, set operations, redundancy
   - Opportunities: priority queue, memoization, incremental updates
   - Invariants: disjoint cover, monotone-shrinking sizes (per audit Q3)

3. **Dynamic Mesh Operation Design** (`research/dynamic_ops_design.md`)
   - Use cases from MADMESHR (advancing-front placement)
   - Consistency rules: When can we add/remove without rebuilding?
   - Transactional model: Batch operations, rollback semantics

4. **Pinch-Point Detection Algorithm** (`research/pinch_point_detection.md`)
   - Definition: Narrow bottlenecks in mesh (e.g., Voronoi diagram medial axis)
   - Algorithm options: Medial axis pruning, bottleneck sampling, distance field
   - Integration with skeletonization (already compute layers—can we reuse?)

5. **API Design Document** (`research/api_design.md`)
   - Proposed new methods: `add_element()`, `remove_vertex()`, etc.
   - Signatures, docstrings, examples
   - Error handling: What if topology breaks?
   - Backwards compatibility: Deprecation path for old methods

---

## ISSUE TRACKING & WORK BREAKDOWN

Will create GitHub issues for each phase:

- [ ] **Research Phase** (Research artifacts, no code changes)
  - Issue: Graph Structure Research & Benchmarking
  - Issue: Skeletonization Algorithm Analysis
  - Issue: Dynamic Mesh Operations Design
  - Issue: Pinch-Point Detection Design

- [ ] **Phase 1: Graph Modernization** (Upgrade adjacencies)
  - Issue: Replace Elem2Edge/Vert2Elem lists with efficient dicts
  - Issue: Implement O(n log n) edge discovery
  - Issue: Add efficiency tests (benchmark suite)

- [ ] **Phase 2: Dynamic Operations** (Add/remove elements)
  - Issue: Implement `add_element()` with consistency checks
  - Issue: Implement `remove_element()` with edge cleanup
  - Issue: Transactional batch API
  - Issue: Regression tests for dynamic operations

- [ ] **Phase 3: Graph Algorithms** (BFS, DFS, pinch-point detection)
  - Issue: Implement BFS/DFS traversal
  - Issue: Implement connected-component analysis
  - Issue: Implement pinch-point detection
  - Issue: Incremental skeletonization updates

- [ ] **Phase 4: Integration** (Bridge to MADMESHR/ADMESH)
  - Issue: MADMESHR advancing-front API
  - Issue: Domain-splitting API for pinch points
  - Issue: ADMESH-Domains bulk-load optimization

---

## TIMELINE & MILESTONES

**All work on branch:** `planning-optimize_modernize` (no immediate shipping)

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

## GOVERNANCE & DOCUMENTATION UPDATES

**Will revise:**
1. **CLAUDE.md** (project charter): Add "Modernize data structures" as standing task
2. **constitution.md** (architectural decisions): Document chosen graph structure, API stability guarantees
3. **PROJECT_PLAN.md** (roadmap): Add phases, timelines, success criteria
4. **CHANGELOG.md**: Version 0.2.0 entry (data structure refactor, new APIs)
5. **API.md** (to be created): Public interface, breaking changes, migration guide

**New files:**
- `GRAPH_DESIGN.md`: Chosen graph structure, invariants, maintenance guide
- `research/` directory with benchmark results, algorithm analyses, design docs

---

## Questions for Stakeholders

1. **Backwards compatibility urgency?** Can we ship 0.2.0 with breaking changes (e.g., moving `layers` dict structure) if migration guide is clear?
2. **Downstream priority?** Which is more urgent: MADMESHR advancing-front support or ADMESH-Domains bulk-load optimization?
3. **ADMESH status:** Is ADMESH (GitHub 404) still under development? Affects Phase 4 planning.
4. **Pinch-point definition:** For MADMESHR domain splitting, what metric defines a "pinch point"? (Voronoi width? Edge length? Bottleneck sampling?)
5. **Testing infrastructure:** Can we enable coverage gates for 0.2.0+? (Audit deferred coverage to 0.1.2)

---

## Appendix: Current Code Paths (For Reference)

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

### Performance Hotspot: _identify_edges()
```python
# Current: O(n²) edge enumeration
for elem in elements:
    for pair in combinations(vertices, 2):
        edges.add(pair)
```

### Performance Hotspot: _build_elem2edge()
```python
# Current: O(n²) edge matching
for elem in elements:
    for edge in edges:
        if elem has edge:
            elem2edge[elem].append(edge_id)
```

---

**Document Status:** DRAFT - Ready for Specification Review & Clarification
**Last Updated:** 2026-04-26
**Author:** Research Agent (Claude)
