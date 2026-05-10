# Phase 3: Graph Algorithms & Optimization Specification

**Branch:** `planning-optimize_modernize`  
**Date:** 2026-04-26  
**Phase:** 3 – Optimization & Algorithm Suite  
**Prerequisite:** Phase 2 complete (dynamic ops API stable)

---

## Executive Summary

Optimize CHILmesh graph operations and add traversal algorithms. Target: <1s initialization for Block_O; BFS, connected components, pinch-point detection for MADMESHR.

**Deliverables:** O(n log n) edge discovery; BFS/DFS; connected-component analysis; pinch-point detection (via skeletonization reuse); optional incremental skeletonization (Phase 3B).

**Success:** Block_O initialization <1s; all algorithm tests passing.

---

## Specification

### Performance Goals

| Metric | Target | Current | Improvement |
|--------|--------|---------|-------------|
| Block_O init time | <1s | ~30s | >30× |
| Edge discovery | O(n log n) | O(n²) | √n |
| Adjacency lookup | O(1) | O(k) | constant |
| Memory overhead | <10% | ~20% | leaner |

### In-Scope Optimizations

#### O(n log n) Edge Discovery

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
    candidates = []
    for elem_id, vertices in enumerate(connectivity_list):
        for i, j in combinations(vertices, 2):
            candidates.append((min(i,j), max(i,j), elem_id))
    candidates.sort()  # O(n log n)
    edges = [tuple(c[:2]) for c in candidates]  # O(n) dedup
    # Build edge2elem adjacency in single pass: O(n)
```

#### BFS/DFS Graph Traversal

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
    # Similar to BFS but use stack
```

#### Connected-Component Analysis

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

#### Pinch-Point Detection (Via Skeletonization)

```python
def pinch_points(self, width_threshold: float) -> List[Tuple[int, int]]:
    """Identify edges forming bottlenecks (narrow layer transitions)."""
    if not self.layers or len(self.layers['OE']) < 2:
        return []
    pinch_edges = []
    for layer_id in range(1, len(self.layers['OE'])):
        outer_edges = self.adjacencies['Edge2Elem'][self.layers['bEdgeIDs'][layer_id]]
        for edge_id, (e1, e2) in outer_edges:
            if self.edge_width(edge_id) < width_threshold:
                pinch_edges.append(edge_id)
    return pinch_edges
```

**Success Metric:** Pinch-point detection completes in <100ms for Block_O.

### Optional Phase 3B: Incremental Skeletonization

Full rebuild likely fast enough; incremental updates are complex. Implement only if Phase 3A complete + benchmarks show benefit.

### Test Coverage

**Performance:**
- [ ] block_o init time <1s; edge discovery on 100K mesh <500ms; BFS/DFS on large mesh <100ms

**Correctness:**
- [ ] BFS produces same vertex set as DFS (connected mesh)
- [ ] Connected components partition all vertices (disjoint union)
- [ ] All Phase 2 regression tests still passing

### Acceptance Criteria

- [ ] `_identify_edges()` refactored to O(n log n)
- [ ] `bfs()`, `dfs()`, `connected_components()`, `pinch_points()` implemented + tested
- [ ] Block_O initialization <1s; all benchmark tests pass
- [ ] All Phase 2 regression tests still passing

### Out of Scope

- ❌ Parallel graph algorithms
- ❌ Spectral methods (Laplacian eigenvalues, etc.)
- ❌ Incremental skeletonization (Phase 3B, TBD)

---

**Timeline:** 3–5 days (O(n log n) edge discovery is critical path)

**Handoff to Phase 4:** All graph ops O(1)–O(n log n); algorithm suite available for MADMESHR integration; Phase 4 focuses on MADMESHR/ADMESH integration only.

