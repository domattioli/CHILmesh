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
1. **_skeletonize()**: Medial axis extraction via boundary peeling
   - Iterative layer identification: boundary → outer elements → inner elements
   - Used for mesh decomposition, visualization, and potentially mesh adaptation
   - Depends on: Edge2Elem, element-to-element connectivity via edges

2. **_build_adjacencies()**: Constructs all adjacency representations
   - `_identify_edges()`: O(n_elems × n_verts_per_elem) set operations
   - `_build_elem2edge()`: **O(n_elems × n_verts × n_edges)** ← BOTTLENECK
   - `_build_vert2edge()`: O(n_edges)
   - `_build_vert2elem()`: O(n_elems × n_verts_per_elem)
   - `_build_edge2elem()`: O(n_elems × n_verts × n_edges) ← BOTTLENECK

3. **Mesh Operations**
   - `smooth_mesh()`: FEM or angle-based smoothing
   - `elem_quality()`: Angular skewness calculation
   - `interior_angles()`: Per-vertex angle computation
   - Mixed-element handling adds complexity (padded triangles as 4-col arrays)

#### Known Performance Issues
| Issue | Location | Complexity | Impact |
|-------|----------|-----------|--------|
| Edge search via linear scan | `_build_elem2edge`, `_build_edge2elem` | O(n) per lookup | O(n²) overall |
| Edge identification via set() | `_identify_edges` | O(1) set ops, but iteration | Acceptable |
| No spatial indexing | All operations | N/A | No point location, nearest-neighbor |
| Mixed-element sentinel handling | `_elem_type`, `_ensure_ccw_orientation` | O(n) iteration | Intricate logic, risk of bugs |
| List-of-lists adjacencies | `Vert2Edge`, `Vert2Elem` | O(1) access, O(k) traversal | Inconsistent with numpy arrays |

### 1.2 Downstream Dependencies

#### MADMESHR (Mesh Adaptation Research)
- **Unknown specifics** (external repo)
- Likely uses: Mesh traversal, layer information, quality metrics
- May need: Rapid neighbor queries, element/edge insertion/deletion

#### ADMESH (Mesh Adaptation)
- **Unknown specifics** (external repo)
- Historical: Mesh refinement, coarsening, quality improvement
- May need: Efficient modification operations, topological queries

#### ADMESH-Domains (Domain Handling)
- **Unknown specifics** (external repo)
- Likely manages: Multiple domains, domain boundaries
- May need: Subgraph queries, domain tagging, cross-domain connectivity

### 1.3 Functional Requirements

**Must Preserve:**
1. Skeletonization (boundary peeling layer identification)
2. Fort.14 I/O (ADCIRC format compatibility)
3. Public API surface (user-facing methods)
4. Mixed-element support (triangles + quads)
5. Plotting/visualization

**Must Support:**
1. Rapid traversal: adjacency queries in O(1) or O(log n)
2. Efficient alteration: node/edge insertion/deletion
3. Bridge interface: clear API for downstream projects
4. Plotting: fast conversion to plottable format

**Must Improve:**
1. Edge building from O(n²) to O(n log n) or O(n)
2. Spatial queries (if needed by downstream projects)
3. Code clarity and maintainability

---

## 2. STATE-OF-THE-ART RESEARCH

### 2.1 Candidate Data Structures

#### A. Half-Edge (Doubly-Connected Edge List - DCEL)
**What it is:** Each edge stored bidirectionally with pointers to next/prev edges in face cycle.

**Pros:**
- O(1) ring traversal around vertex
- O(1) face boundary traversal
- Elegant handling of 2D topology
- Standard in computational geometry (CGAL, JChem)

**Cons:**
- Higher memory overhead per edge (4-6 pointers vs 2 integers)
- Complex insertion/deletion (must maintain face cycles)
- May be overkill for read-heavy skeletonization
- Not directly compatible with numpy arrays

**Verdict for CHILmesh:** ⚠️ Promising but operationally heavy

#### B. Winged-Edge Structure
**What it is:** Similar to half-edge but edges store references to adjacent faces and boundary edges.

**Pros:**
- Slightly simpler than half-edge in some use cases
- Still O(1) local topology queries

**Cons:**
- Similar memory overhead
- Similar complexity for insertion/deletion
- No significant advantage over half-edge for our use case

**Verdict for CHILmesh:** ❌ Less preferable than half-edge

#### C. Hybrid: Explicit Hash Map Edge Lookup
**What it is:** Keep numpy arrays but replace linear edge search with dict/hash table mapping (v1, v2) → edge_id.

**Pros:**
- Minimal code changes needed
- O(1) edge lookup (amortized)
- Stays compatible with numpy arrays
- Easy to implement incrementally
- Preserves current adjacency structure

**Cons:**
- Does NOT fix list-of-lists adjacencies (Vert2Edge, Vert2Elem)
- Doesn't improve memory layout for cache efficiency

**Verdict for CHILmesh:** ✅ Quick win, low risk

#### D. Explicit Sparse CSR/CSC Matrices
**What it is:** Store adjacencies as scipy.sparse matrices (already used in smoothing code).

**Pros:**
- O(1) lookup via .getrow() / .getcol()
- Memory efficient for sparse graphs
- Integrates with scipy ecosystem
- Fast matrix operations (transposition, multiplication)
- Good for graph algorithms (spectral clustering, etc.)

**Cons:**
- Overkill for small meshes (<100k elements)
- Different access pattern from current code
- CSR format inefficient for insertion

**Verdict for CHILmesh:** ⚠️ Good for large meshes, adds dependency on existing code

#### E. NetworkX Graph
**What it is:** Use NetworkX DiGraph or MultiGraph to represent mesh topology.

**Pros:**
- Rich graph algorithms library (BFS, DFS, shortest path, connected components)
- O(1) neighbor queries
- Natural representation of topology
- Excellent for skeletonization (can use BFS/DFS)

**Cons:**
- Heavy dependency (~2500 lines, many features we don't use)
- Python objects slower than numpy for tight loops
- Adds layer of indirection
- May not preserve element order for fort.14 I/O
- Not designed for geometric meshes (stores only topology)

**Verdict for CHILmesh:** ⚠️ Overkill, but worth considering if bridge needs complex graph algorithms

#### F. Refined Numpy + Dict Hybrid
**What it is:** Extend current approach:
- Keep numpy arrays for dense data (points, connectivity, most adjacencies)
- Use sorted arrays + binary search for edge lookups
- Use dict for Vert2Edge and Vert2Elem
- Use CSR matrices for large sparse operations

**Pros:**
- Minimal dependency additions
- Leverages existing ecosystem (numpy, scipy)
- Scales well from small to large meshes
- Clear separation of concerns (dense vs sparse)
- Incremental migration path

**Cons:**
- More complex implementation than single approach
- Requires careful design to avoid mixing paradigms

**Verdict for CHILmesh:** ✅ Best match for requirements

### 2.2 Skeletonization-Specific Considerations

The current `_skeletonize()` algorithm:
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

**Key observations:**
- Doesn't depend on complex topological operations
- Needs efficient: boundary identification, neighbor finding
- Layer data structure (lists of arrays per layer) is separate concern

**Data structure implications:**
- Any structure supporting O(1) edge→elements is sufficient
- BFS/layer-by-layer processing doesn't require pointer-based rings
- Current approach actually reasonable for this specific algorithm

### 2.3 Best Practice from Literature

**Computational Geometry Landmarks:**
- CGAL (Computational Geometry Algorithms Library): Uses DCEL/half-edge
- Triangle (J.R. Shewchuk): Uses edge-based lists + pointer structures
- Gmsh: Uses hybrid approach (connectivity arrays + lookup maps)
- VTK: Uses cell/point arrays + incidence structures

**Mesh Optimization Tools:**
- DEAL.II (FEM): Explicit adjacency matrices + sparse storage
- PETSc: Distributed mesh with CSR-based topology
- Exodus II: Optimized I/O of implicit connectivity

**Key insight:** No single universal solution—depends on:
1. Size (small: dense arrays, large: sparse structures)
2. Operations (read-heavy: optimize lookup, write-heavy: optimize insertion)
3. Integration (match ecosystem conventions)

---

## 3. STRATEGIC PLAN

### 3.1 Recommended Approach: Layered Modernization

**Phase 1: Foundation (Immediate, low risk)**
- Implement hash map edge lookup (O(n) → O(1) amortized)
- Refactor `_build_elem2edge` and `_build_edge2elem`
- Target: Eliminate O(n²) bottleneck, ~2x speedup expected

**Phase 2: Adjacency Modernization (Medium complexity)**
- Migrate Vert2Edge, Vert2Elem to explicit dicts
- Update all traversal patterns (for-in-loops → dict lookups)
- Add type hints and validation

**Phase 3: Bridge Infrastructure (Requires coordination with downstream)**
- Define clear "CHILmesh Access Interface" (CAI)
- Document invariants and traversal patterns
- Create adapter methods for MADMESHR/ADMESH/ADMESH-Domains
- Add integration tests

**Phase 4: Advanced Optimizations (Later, if needed)**
- Spatial indexing (KD-tree for point queries)
- CSR matrices for large sparse operations
- Caching strategies
- Parallel operations (numpy vectorization improvements)

### 3.2 Design Principles

1. **Backward Compatibility First**
   - Public API unchanged until v1.0
   - Internal refactoring hidden behind same methods
   - Tests must pass without modification

2. **Incremental Migration**
   - One adjacency at a time
   - Each phase adds new methods, keeps old ones deprecated
   - Clear migration path documented

3. **Explicit Over Implicit**
   - Replace magic list-of-lists with documented dict structures
   - Add validation at adjacency building time
   - Clear error messages on invariant violations

4. **Cache Locality & CPU Efficiency**
   - Prefer numpy arrays over Python objects
   - Use structured arrays where appropriate
   - Avoid unnecessary copying

5. **Documentation-First Design**
   - Document invariants before code
   - Every adjacency structure has a docstring
   - Graph representation documented in comments

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

**Tasks:**
1. [ ] **TASK-P1-01**: Create `EdgeMap` class with insertion/lookup/removal
   - Type: Implementation
   - Scope: `src/chilmesh/mesh_topology.py` (new file)
   - Deliverable: Class with `add_edge`, `find_edge`, `remove_edge` methods
   - Test: Unit tests for edge mapping

2. [ ] **TASK-P1-02**: Refactor `_identify_edges()` to use EdgeMap
   - Type: Refactoring
   - Scope: `CHILmesh._identify_edges()`
   - Change: Return both edge list and edge_map
   - Test: Same semantics, verify output unchanged

3. [ ] **TASK-P1-03**: Optimize `_build_elem2edge()` from O(n²) to O(n log n)
   - Type: Performance
   - Scope: `CHILmesh._build_elem2edge()`
   - Change: Use edge_map.find_edge() instead of linear search
   - Benchmark: Time before/after on all fixtures

4. [ ] **TASK-P1-04**: Optimize `_build_edge2elem()` similarly
   - Type: Performance
   - Scope: `CHILmesh._build_edge2elem()`
   - Benchmark: Time before/after

5. [ ] **TASK-P1-05**: Store EdgeMap in adjacencies dict
   - Type: Integration
   - Scope: `CHILmesh._build_adjacencies()`
   - Store: `adjacencies['EdgeMap'] = edge_map`
   - Maintain backward compat: old methods still work

6. [ ] **TASK-P1-06**: Performance regression tests
   - Type: Testing
   - Scope: `tests/test_performance_edge_building.py`
   - Verify: All fixtures build faster than baseline
   - Accept: Only if ≥1.5x speedup on large fixtures

### 4.2 Phase 2: Adjacency Modernization

**Tasks:**
1. [ ] **TASK-P2-01**: Migrate `Vert2Edge` to explicit dict
   - Type: Refactoring
   - Scope: `CHILmesh._build_vert2edge()`
   - Change: Return `Dict[int, Set[int]]` instead of `List[List[int]]`
   - Maintain: `adjacencies['Vert2Edge']` format

2. [ ] **TASK-P2-02**: Migrate `Vert2Elem` to explicit dict
   - Type: Refactoring
   - Scope: `CHILmesh._build_vert2elem()`
   - Change: Return `Dict[int, Set[int]]`

3. [ ] **TASK-P2-03**: Add explicit accessor methods
   - Type: API
   - Scope: `CHILmesh` class
   - Add: `get_vertex_edges(vert_id)`, `get_vertex_elements(vert_id)`
   - Return: Consistent `Set[int]` type

4. [ ] **TASK-P2-04**: Update internal traversal patterns
   - Type: Refactoring
   - Scope: All methods iterating over Vert2Edge, Vert2Elem
   - Search: `for edge_id in self.adjacencies['Vert2Edge'][...]`
   - Replace: Use new dict-based iteration

5. [ ] **TASK-P2-05**: Type hints and validation
   - Type: Code quality
   - Scope: `CHILmesh` class
   - Add: Full type annotations for adjacency methods
   - Add: Validation (keys exist, values are sets, etc.)

6. [ ] **TASK-P2-06**: Adjacency structure documentation
   - Type: Documentation
   - Scope: `CHILmesh` class docstring
   - Document: Each adjacency format, invariants, access patterns
   - Example: Illustrate correct usage

### 4.3 Phase 3: Bridge Infrastructure

**Tasks:**
1. [ ] **TASK-P3-01**: Define CHILmesh Access Interface (CAI)
   - Type: Design
   - Scope: `docs/CHILmesh_Access_Interface.md` (new)
   - Content: Stable API for downstream projects
   - Guarantee: Backward compatibility guarantee

2. [ ] **TASK-P3-02**: Create bridge module
   - Type: Implementation
   - Scope: `src/chilmesh/bridge.py` (new)
   - Implement: Adapter classes for MADMESHR/ADMESH/ADMESH-Domains
   - Example: `MeshAdapterForMADMESHR`, etc.

3. [ ] **TASK-P3-03**: Integration tests with mock downstream projects
   - Type: Testing
   - Scope: `tests/test_bridge_integration.py`
   - Simulate: MADMESHR, ADMESH, ADMESH-Domains usage patterns
   - Verify: Clear contracts, no hidden dependencies

4. [ ] **TASK-P3-04**: Document migration for downstream projects
   - Type: Documentation
   - Scope: `docs/DOWNSTREAM_MIGRATION_GUIDE.md` (new)
   - Content: How MADMESHR/ADMESH/ADMESH-Domains should integrate
   - Examples: Code snippets for common operations

---

## 5. RISK ASSESSMENT & MITIGATION

| Risk | Impact | Probability | Mitigation |
|------|--------|-------------|-----------|
| Skeletonization breaks | Critical | Low | Retain original _skeletonize, add comprehensive tests |
| Fort.14 I/O breaks | Critical | Low | All existing tests must pass unchanged |
| Performance regression | High | Medium | Benchmark every step, revert if slower |
| API breakage | High | Medium | Deprecation warnings, parallel APIs until v1.0 |
| Memory bloat from extra maps | Medium | Low | Profile memory usage, optimize if >10% overhead |
| Downstream project incompatibility | High | High | Early coordination with MADMESHR/ADMESH authors |

---

## 6. SUCCESS CRITERIA

### Phase 1 (Edge Mapping)
- [ ] `_build_elem2edge` and `_build_edge2elem` are O(n log n) or O(n)
- [ ] All existing tests pass without modification
- [ ] Benchmarks show ≥1.5x speedup on block_o fixture (largest)
- [ ] No new dependencies added

### Phase 2 (Adjacency Modernization)
- [ ] `Vert2Edge` and `Vert2Elem` are explicit dicts
- [ ] All traversal patterns updated
- [ ] Type hints complete for adjacency methods
- [ ] Documentation updated with invariants

### Phase 3 (Bridge Infrastructure)
- [ ] CAI document published
- [ ] Bridge module has adapters for all three downstream projects
- [ ] Integration tests cover common usage patterns
- [ ] No external dependencies on implementation details

### Overall
- [ ] 100% test pass rate on all fixtures
- [ ] No breaking changes to public API
- [ ] Performance improved by ≥1.5x on large meshes
- [ ] Code ready for v0.2.0 release

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

1. **Get Stakeholder Buy-In**
   - Share this plan with MADMESHR/ADMESH/ADMESH-Domains authors
   - Identify any additional requirements
   - Adjust timeline if needed

2. **Create GitHub Issues**
   - One issue per task (16 issues total)
   - Link to this planning document
   - Assign to sprints/milestones

3. **Set Up Branch Structure**
   - Keep `claude/zen-fermi-NGYbR` for planning
   - Create `feature/phase1-edge-mapping` for Phase 1 work
   - Create separate branches for phases 2 and 3

4. **Begin Phase 1**
   - Start with TASK-P1-01 (EdgeMap class)
   - Get code review before proceeding
   - Benchmark as you go

---

## 9. RELATED DOCUMENTS

- `README.md`: Overview of CHILmesh functionality
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
