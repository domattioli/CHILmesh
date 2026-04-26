# Phase 3: Graph Algorithms & Optimization Specification

**Branch:** `planning-optimize_modernize`  
**Date:** 2026-04-26  
**Phase:** 3 – Optimization & Algorithm Suite  
**Methodology:** Spec-Kit  
**Prerequisite:** Phase 2 complete (dynamic ops API stable)

---

## Executive Summary

Phase 3 optimizes CHILmesh's graph operations and adds advanced traversal algorithms. The focus is achieving <1s initialization for Block_O and supporting algorithms for MADMESHR (BFS, connected components, pinch-point detection).

**Deliverables:**
- O(n log n) edge discovery (replace O(n²) _identify_edges)
- BFS/DFS graph traversal methods
- Connected-component analysis
- Pinch-point detection (via skeletonization reuse)
- Optional incremental skeletonization (Phase 3B, if time allows)

**Success:** Block_O initialization <1s, all algorithm tests passing.

---

## Specification

### 3.1 Performance Goals

| Metric | Target | Current | Improvement |
|--------|--------|---------|-------------|
| Block_O init time | <1s | ~30s | >30× speedup |
| Edge discovery | O(n log n) | O(n²) | √n improvement |
| Adjacency lookup | O(1) | O(k) | Better worst-case |
| Memory overhead | <10% | ~20% | Leaner impl |

### 3.2 In-Scope Optimizations

#### **3.2.1 O(n log n) Edge Discovery**

**Current (O(n²)):**
```python
def _identify_edges(self):
    edges = set()
    for elem_id in range(n_elems):
        vertices = connectivity_list[elem_id]
        for i, j in combinations(vertices, 2):
            edges.add((min(i,j), max(i,j)))  # O(n²) across all elements
```

**Optimized (O(n log n)):**
```python
def _identify_edges_optimized(self):
    # Collect all (vertex_pair, elem_id) tuples, sort, deduplicate
    candidates = []
    for elem_id, vertices in enumerate(connectivity_list):
        for i, j in combinations(vertices, 2):
            candidates.append((min(i,j), max(i,j), elem_id))
    candidates.sort()  # O(n log n)
    edges = [tuple(c[:2]) for c in candidates]  # O(n) dedup
    # Build edge2elem adjacency in single pass: O(n)
```

**Success Metric:** Sorting + dedup overhead < 100ms for 100K vertices.

#### **3.2.2 BFS/DFS Graph Traversal**

**New Methods:**
```python
def bfs(self, start_vertex: int, callback=None) -> List[int]:
    """Breadth-first search from vertex. Optionally apply callback to each."""
    visited = set()
    queue = deque([start_vertex])
    order = []
    while queue:
        v = queue.popleft()
        if v in visited: continue
        visited.add(v)
        order.append(v)
        if callback: callback(v)
        for neighbor in self.vert2edge[v]:  # O(1) lookup
            queue.append(neighbor)
    return order

def dfs(self, start_vertex: int, callback=None) -> List[int]:
    """Depth-first search. Useful for topological analysis."""
    # Similar to BFS but use stack instead of deque
```

**Test Cases:**
- [ ] BFS from vertex 0, verify all reachable vertices visited
- [ ] DFS returns same vertex set as BFS (for connected mesh)
- [ ] BFS order respects distance (layers match skeletonization if applicable)

#### **3.2.3 Connected-Component Analysis**

**New Method:**
```python
def connected_components(self) -> Dict[int, List[int]]:
    """Partition vertices into connected components."""
    components = {}
    comp_id = 0
    unvisited = set(range(self.n_verts))
    while unvisited:
        seed = unvisited.pop()
        component = set(self.bfs(seed))
        components[comp_id] = sorted(component)
        unvisited -= component
        comp_id += 1
    return components
```

**Use Case:** Detect if mesh is singly-connected (for MADMESHR domain splitting).

#### **3.2.4 Pinch-Point Detection (Via Skeletonization)**

**Hypothesis:** Skeletonization (boundary peeling) already identifies narrow regions. Can we reuse layer data?

**Algorithm (Sketch):**
```python
def pinch_points(self, width_threshold: float) -> List[Tuple[int, int]]:
    """Identify edges that form bottlenecks (narrow layer transitions)."""
    if not self.layers or len(self.layers['OE']) < 2:
        return []
    
    pinch_edges = []
    for layer_id in range(1, len(self.layers['OE'])):
        outer_edges = self.adjacencies['Edge2Elem'][self.layers['bEdgeIDs'][layer_id]]
        for edge_id, (e1, e2) in outer_edges:
            # Check if edge transitions between two narrow layers
            if self.edge_width(edge_id) < width_threshold:
                pinch_edges.append(edge_id)
    return pinch_edges
```

**Success Metric:** Pinch-point detection completes in <100ms for Block_O.

### 3.3 Optional Phase 3B: Incremental Skeletonization

If time allows after primary goals, optimize layer rebuilds:

- **Current:** Full rebuild every time topology changes
- **Optimized:** Track "frontier" edges (boundary of current layer), update incrementally

**Risk:** Incremental updates are complex; full rebuilds are simpler and likely fast enough.

**Decision:** Implement only if Phase 3A complete + benchmarks show benefit.

### 3.4 Test Coverage

**Performance Tests (Benchmarks):**
- [ ] block_o init time <1s (vs. 30s current)
- [ ] edge discovery on 100K mesh <500ms
- [ ] BFS/DFS on large mesh <100ms

**Correctness Tests:**
- [ ] BFS produces same order as manual traversal (on small mesh)
- [ ] Connected components partition all vertices (disjoint union)
- [ ] Pinch-point detection matches manual inspection (on test fixtures)
- [ ] All Phase 2 regression tests still passing

### 3.5 Acceptance Criteria

- [ ] `_identify_edges()` refactored to O(n log n)
- [ ] `bfs()`, `dfs()` methods implemented and tested
- [ ] `connected_components()` method implemented and tested
- [ ] `pinch_points()` method implemented (reuses skeletonization)
- [ ] Block_O initialization <1s
- [ ] All benchmark tests pass
- [ ] All Phase 2 regression tests still passing
- [ ] Docstrings + examples for all new methods
- [ ] No performance regression for existing methods

### 3.6 Out of Scope (Phase 3)

- ❌ Parallel graph algorithms
- ❌ Spectral methods (Laplacian eigenvalues, etc.)
- ❌ Force-directed layout or mesh visualization improvements
- ❌ Incremental skeletonization (Phase 3B, TBD)

---

## Timeline

**3–5 days** (O(n log n) edge discovery is critical path; algorithms are secondary)

---

## Handoff to Phase 4

Once Phase 3 is complete:
1. All graph operations are O(1)–O(n log n)
2. Algorithm suite available for MADMESHR integration
3. Performance goals met (Block_O <1s)
4. Phase 4 focuses on MADMESHR/ADMESH integration only

