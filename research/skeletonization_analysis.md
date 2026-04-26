# Skeletonization Algorithm Analysis

## Overview

Mesh skeletonization (a.k.a. "medial axis extraction via boundary peeling") is a key CHILmesh feature that decomposes the mesh into hierarchical layers by iteratively removing concentric "shells" from the boundary inward. This document analyzes the current algorithm, identifies inefficiencies, and proposes optimizations.

---

## Current Algorithm (Python Implementation)

**Location:** `CHILmesh.py:_skeletonize()` lines 511–616

### Algorithm Pseudocode

```
skeletonize():
  layers = {OE: [], IE: [], OV: [], IV: [], bEdgeIDs: []}
  remaining_elements = set(all element IDs)
  boundary_edges = set(edges with only 1 adjacent element)

  while remaining_elements non-empty and boundary_edges non-empty:
    # Outer layer (elements adjacent to boundary)
    outer_verts = vertices(boundary_edges)
    outer_elems = {elem : any edge of elem in boundary_edges and elem in remaining}
    
    # Inner layer (neighbors of outer elements)
    inner_elems = {elem : neighbor of outer_elem and elem in remaining}

    # Store layer
    layers.OE.append(outer_elems)
    layers.IE.append(inner_elems)
    layers.OV.append(outer_verts)
    layers.IV.append(vertices(inner_elems) - outer_verts)
    layers.bEdgeIDs.append(boundary_edges)

    # Remove from remaining
    remaining -= outer_elems
    remaining -= inner_elems

    # Recompute boundary edges for next iteration
    boundary_edges = {edge : exactly one endpoint in remaining}

  return layers
```

### Key Data Structures Used

| Structure | Type | Complexity | Notes |
|-----------|------|-----------|-------|
| `remaining_elements` | `set(int)` | O(1) membership | OK |
| `boundary_edges` | `np.ndarray(int)` or `list` | O(n) iteration | OK |
| `edge2elem` | `np.ndarray(E×2)` | O(1) lookup | OK (precomputed) |
| `edge2vert` | `np.ndarray(E×2)` | O(1) lookup | OK (precomputed) |
| `elem2elem` | `list[list]` | Built in loop | **Problem: O(n²) to build** |

---

## Complexity Analysis

### Current Implementation Complexity

Let:
- N = number of vertices
- M = number of elements
- E = number of edges
- L = number of layers (typically O(√M) for 2D mesh)

**Initialization (per layer):**
1. Outer vertices: O(|boundary_edges| × 2) = O(B) where B is boundary size
2. Outer elements: O(B) (each edge has ≤2 adjacent elements)
3. Inner elements: O(|outer_elems| × degree) = O(|outer_elems| × 8) ≈ O(|outer_elems|)
4. Inner vertices: O(|outer_elems| + |inner_elems|) × 4 vertices/element = O(|OE| + |IE|)

**Boundary recomputation:**
```python
for edge_idx, (e1, e2) in enumerate(edge2elem):
  if (e1 in remaining) != (e2 in remaining):
    boundary_edges.append(edge_idx)
```
- O(E) iterations per layer
- O(1) set membership check (set is fast)
- **Total per layer: O(E)**

**Sum over all layers:**
- Outer/inner element discovery: O(M) total (each element processed once)
- Boundary recomputation: O(L × E)
- Set membership checks: O(M) total (one per element)

**Total:** O(L × E) = O(√M × E)

For typical 2D mesh: **O(√M × M) = O(M^1.5)** (since E ≈ 3M for triangular mesh)

### Inefficiencies

1. **O(E) boundary recomputation per layer:** Iterates all edges even though only "frontier" edges change. Could maintain frontier incrementally.
2. **Set membership in loop:** While O(1) per check, the constant is non-trivial for large sets.
3. **Repeated vertex gathering:** `inner_vertices = all_vertices - outer_vertices` can be computed during element traversal.
4. **`elem2elem` construction on-the-fly:** Built inside `_skeletonize()` but already partially exists in `adjacencies['Edge2Elem']`.

---

## Audit Notes (Q3)

**Q3 (from PROGRESS.md):** Investigate IE (inner elements) divergence in skeletonization.

**Finding:** Existing Python implementation produces a disjoint cover with monotone-shrinking layer sizes on convex annulus and structured fixtures. The algorithm:
- Expands outward from the boundary
- Peels two layers at a time (OE then IE)
- Produces valid skeletonization as measured by `test_layers_disjoint_cover`

**Resolution:** Current behavior is preserved as an invariant in tests. No algorithmic change needed (per audit guidance).

---

## Optimization Opportunities

### Opportunity 1: Incremental Frontier Tracking (O(E) → O(frontier))

Instead of iterating all edges, maintain a frontier queue:

```python
frontier = PriorityQueue()  # Edges with exactly one endpoint in remaining
for edge in boundary_edges:
  frontier.put(edge, priority=...)

while frontier not empty:
  edge = frontier.pop()
  if exactly_one_endpoint_in_remaining(edge):
    yield edge
  else:
    frontier.discard(edge)  # Both endpoints processed
```

**Benefit:** Avoids O(E) iteration; processes only changed edges.
**Complexity:** O(frontier size) per layer, typically O(|boundary|) ≈ O(√M).
**Cost:** Requires tracking edge endpoints carefully; must handle removals from frontier.

### Opportunity 2: Layer-Aware Skeletonization (Incremental Updates)

If topology changes (add/remove element), instead of full re-skeletonization:
1. Identify affected layers (layers containing added/removed elements)
2. Recompute only affected layers and downstream

**Benefit:** Enables incremental mesh updates (for MADMESHR advancing-front).
**Complexity:** O(affected layers × E) instead of O(L × E).
**Cost:** Requires careful invalidation of layer relationships.

### Opportunity 3: Medial Axis Pruning (For Pinch-Point Detection)

Current skeletonization produces full medial axis. For domain splitting (MADMESHR use case), we only need bottlenecks:

```python
def pinch_points(mesh):
  """Find narrow regions (bottlenecks) in mesh."""
  skeleton = skeletonize(mesh)
  for layer in skeleton.layers:
    frontier_width = measure_layer_width(layer)  # Voronoi vertex spacing?
    if frontier_width < threshold:
      yield layer.boundary_edges  # These are potential split points
```

**Benefit:** Provides MADMESHR domain-splitting API.
**Complexity:** O(L × perimeter) to measure widths.
**Cost:** Metric for "width" is domain-specific (geometry vs. graph metric).

---

## Proposed Implementation (Phase 3)

### Stage 1: Profiling (No code change)
1. Measure `_skeletonize()` runtime on annulus, structured, block_o
2. Identify actual bottleneck (frontier recomputation? set membership? vertex gathering?)
3. Verify that **O(E) per layer** is the limiting factor (likely is)

### Stage 2: Refactor for Clarity (O(M^1.5) → O(M^1.5), same complexity, clearer code)
1. Extract `elem2elem` construction into separate method (cache result)
2. Document invariants: layer disjointness, monotone shrinking
3. Add type hints: `remaining_elements: Set[int]`, etc.

### Stage 3: Optimize Frontier (O(E/L) expected speedup)
1. Implement frontier queue tracking
2. Benchmark against current implementation
3. Gate: Must preserve output (layer structure identical)

### Stage 4: Incremental Layers (Phase 2 dependency: requires add/remove element API)
1. Design incremental update semantics
2. Implement `add_element() → update_layers()`
3. Regression: All tests must pass; skeletonization output unchanged for static meshes

---

## Invariants to Preserve

### Hard Invariants (Must Hold After Modernization)

1. **Disjoint cover:** Every element appears in exactly one layer (OE[i] or IE[i])
   ```python
   all_elems = set()
   for i in range(n_layers):
     all_elems.update(layers['OE'][i])
     all_elems.update(layers['IE'][i])
   assert len(all_elems) == n_elems
   ```

2. **Monotone-shrinking layer sizes (on convex meshes):**
   ```python
   for i in range(n_layers - 1):
     assert len(layers['OE'][i]) >= len(layers['OE'][i+1])
   ```
   (Note: Not guaranteed for non-convex meshes; add comment)

3. **Layer boundary validity:** Each layer's boundary edges are exactly those with one element in the layer and one outside
   ```python
   for i, edge in enumerate(layers['bEdgeIDs'][i]):
     e1, e2 = edge2elem[edge]
     assert (e1 in layers['OE'][i] or e1 in layers['IE'][i]) XOR (e2 in layers['OE'][i] or e2 in layers['IE'][i])
   ```

4. **Vertex containment:** All vertices of elements in layer i must be in OV[i] ∪ IV[i]
   ```python
   for elem in layers['OE'][i] + layers['IE'][i]:
     for vert in elem2vert[elem]:
       assert vert in layers['OV'][i] or vert in layers['IV'][i]
   ```

### Soft Invariants (Currently Hold, Nice to Preserve)

5. **Outer vertices = boundary vertices:** OV[i] = vertices of boundary_edges[i]
   (Could be optimized away if carefully recomputed)

6. **Layer peeling order:** Layers numbered from boundary inward (layer 0 touches original mesh boundary)
   (Natural consequence of BFS; document as desired property)

---

## Testing Strategy

### Current Tests
- `test_layers_disjoint_cover`: Validates invariant #1 (disjoint cover)
- `test_layers_annulus.py`: Visual/numerical validation on specific meshes

### New Tests Needed

```python
def test_skeletonization_invariants():
  """Test all hard invariants on parametrized fixtures."""
  for mesh in [annulus, structured, donut, block_o]:
    # Invariant 1: Disjoint cover
    # Invariant 2: Monotone shrinking (if convex)
    # Invariant 3: Boundary validity
    # Invariant 4: Vertex containment

def test_skeletonization_incremental_update():
  """Test that add_element/remove_element preserve layers correctly."""
  # Add element → recompute affected layers
  # Remove element → recompute affected layers
  # Result matches full re-skeletonization

def test_pinch_point_detection():
  """Test pinch-point API for MADMESHR use case."""
  # Load mesh, run skeletonization
  # Identify bottlenecks
  # Verify detected points are geometrically narrow
```

---

## References

- [Medial Axis Wikipedia](https://en.wikipedia.org/wiki/Medial_axis)
- [Image Skeletonization (Zhang-Suen)](https://en.wikipedia.org/wiki/Skeletonization)
- Mattioli Thesis: "Mesh Layers" concept (deprecated name) ≈ medial axis extraction via BFS
- Audit Note Q3 (PROGRESS.md): Disjoint cover and monotone-shrinking sizes preserved

---

## Next Steps

1. Profile current implementation (identify actual hotspot)
2. Review proposed frontier optimization with team
3. Design incremental update semantics (post Phase 2)
4. Implement stages 2–4 with full test coverage

**Status:** DRAFT - Ready for Phase 3 planning
