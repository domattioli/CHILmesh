# Dynamic Mesh Operations Design

## Purpose

Enable CHILmesh to support **incremental mesh modification** (adding/removing elements, vertices, edges) while maintaining consistency of adjacency structures and skeletonization. This is essential for:

1. **MADMESHR integration:** Advancing-front mesh generation places elements one at a time
2. **Mesh adaptation:** Future ADMESH use cases (refinement, coarsening, edge swapping)
3. **Interactive mesh editing:** Domain splitting, boundary modification

---

## Use Cases

### Use Case 1: Advancing-Front Element Placement (MADMESHR)

```python
# Pseudocode: RL agent places elements iteratively
mesh = CHILmesh(...)  # Start with empty or seed mesh

while frontier.has_elements():
  state = extract_state(frontier)  # Boundary configuration
  action = agent.select_action(state)  # Learn: quad vs triangle, position, angle
  
  new_element = place_element(action)
  mesh.add_element(
    vertices=[v1, v2, v3, (v4 if quad else v3)],  # Padded if triangle
    elem_type='quad' if quad else 'tri'
  )
  # Mesh must update:
  # - connectivity_list
  # - adjacencies (Elem2Edge, Vert2Elem, Edge2Elem, etc.)
  # - skeletonization (affected layers only, not full rebuild)
  
  frontier = update_frontier(mesh)  # Find next placement site
```

**Requirements:**
- `add_element(vertices, elem_type)` → integer elem_id
- Incremental layer updates (avoid full `_skeletonize()` per element)
- Frontier management (expose boundary edges, advance pointers)
- Domain splitting (detect pinch points, split domain)

---

### Use Case 2: Adaptive Mesh Refinement (ADMESH)

```python
# Refine element based on error estimate
mesh = load_mesh(...)
error = solve_pde_and_estimate_error(mesh)

for elem_id in high_error_elements:
  # Option A: Split quad into 4 quads
  # Option B: Split triangle into 4 triangles
  # Option C: Insert vertex at centroid + remesh surrounding
  
  mesh.split_element(elem_id, strategy='split_4')
  # Updates adjacencies + skeletonization
```

**Requirements:**
- `split_element(elem_id, strategy)`
- `insert_vertex(x, y, z) → vert_id`
- `remove_vertex(vert_id, strategy='collapse_edges')`
- `swap_edge(edge_id)` (flip diagonal in quad pair)

---

### Use Case 3: Domain Splitting (MADMESHR)

```python
# Detect bottlenecks in domain, split for multi-region meshing
mesh = CHILmesh(...)
pinch_points = mesh.pinch_points(threshold=0.1)

for pinch_edge in pinch_points:
  v1, v2 = mesh.edge2vert(pinch_edge)
  # Remove pinch edge, separating mesh into two domains
  mesh.remove_edge(pinch_edge)
  # Returns list of connected components
  domains = mesh.connected_components()
```

**Requirements:**
- `pinch_points(threshold)` → [edge_ids]
- `remove_edge(edge_id)` → [component_ids]
- `connected_components()` → [[vert_ids], ...]

---

## API Design (Proposed)

### Element Operations

#### `add_element(vertices, elem_type='auto')`

```python
def add_element(
    self,
    vertices: np.ndarray | List[int],  # Shape (3,) or (4,) 
    elem_type: str = 'auto'  # 'tri', 'quad', 'auto'
) -> int:
  """
  Add a new element to the mesh.
  
  Parameters:
    vertices: Vertex indices. Shape (3,) for triangle or (4,) for quad.
              If (3,), padded to (4,) with vertices[2].
    elem_type: 'tri' (force triangle), 'quad' (force quad), or 'auto' (infer from length).
  
  Returns:
    New element ID.
  
  Raises:
    ValueError: If vertices not in [0, n_verts) or invalid elem_type.
    TopologyError: If elem already exists (duplicate edge).
  
  Notes:
    - Updates connectivity_list, edge structures, adjacencies
    - Incremental layer update (only affected layers recomputed)
    - Returned elem_id can be used for queries (elem2edge, elem2vert, etc.)
  
  Example:
    >>> new_elem_id = mesh.add_element([1, 2, 3, 4], elem_type='quad')
    >>> mesh.elem2vert(new_elem_id)
    array([1, 2, 3, 4])
  """
```

**Consistency Checks:**
- All vertex IDs must exist: `all(0 <= v < n_verts for v in vertices)`
- If triangle (len=3): pad to 4 automatically
- Elements must not be degenerate (allow 0-area for now, flag in quality check)
- New edges: Check if already exist (consistency)

**Adjacency Updates:**
```python
# Pseudocode for internal updates
def _add_element_internal(vertices, elem_type):
  elem_id = n_elems
  n_elems += 1
  
  # 1. Update connectivity
  connectivity_list[elem_id] = vertices
  
  # 2. Update/create edges
  for i, (v1, v2) in enumerate(edges_of_elem(vertices)):
    edge_key = tuple(sorted([v1, v2]))
    if edge_key not in edge_id_map:
      edge_id = create_new_edge(edge_key)
    else:
      edge_id = edge_id_map[edge_key]
    
    elem2edge[elem_id].append(edge_id)
    vert2edge[v1].add(edge_id)
    vert2edge[v2].add(edge_id)
  
  # 3. Update elem2elem (for skeletonization)
  for edge_id in elem2edge[elem_id]:
    e1, e2 = edge2elem[edge_id]
    if e1 != elem_id and e1 >= 0:
      elem2elem[elem_id].append(e1)
      elem2elem[e1].append(elem_id)
    if e2 != elem_id and e2 >= 0:
      elem2elem[elem_id].append(e2)
      elem2elem[e2].append(elem_id)
  
  # 4. Update vert2elem
  for v in vertices:
    vert2elem[v].add(elem_id)
  
  # 5. Incremental layer update
  update_layers_for_new_element(elem_id)
```

---

#### `remove_element(elem_id)`

```python
def remove_element(self, elem_id: int) -> List[int]:
  """
  Remove element from mesh, cleaning up orphaned edges.
  
  Parameters:
    elem_id: Element to remove.
  
  Returns:
    List of edge IDs that became isolated (no adjacent elements).
  
  Raises:
    IndexError: If elem_id not in [0, n_elems).
    
  Notes:
    - Orphaned edges (no adjacent element) are NOT automatically removed
      (caller can optionally remove via remove_edge).
    - This preserves vertex ordering, but gaps form in element list.
      Use defragment() if dense numbering needed (not recommended; use sparse IDs).
    - Incremental layer update only.
  
  Example:
    >>> orphaned_edges = mesh.remove_element(5)
    >>> for edge_id in orphaned_edges:
    ...   mesh.remove_edge(edge_id)  # Optional cleanup
  """
```

**Consistency Checks:**
- Element must exist: `0 <= elem_id < n_elems`
- Update adjacencies: remove from edge2elem, elem2edge, vert2elem, elem2elem

---

### Vertex Operations

#### `add_vertex(x, y, z=0.0)`

```python
def add_vertex(self, x: float, y: float, z: float = 0.0) -> int:
  """
  Add a new vertex to the mesh (initially unconnected).
  
  Parameters:
    x, y, z: Coordinates.
  
  Returns:
    New vertex ID.
  
  Notes:
    - Vertex is created but has no adjacent elements (isolated node).
    - Caller must use add_element() to connect it.
  
  Example:
    >>> v_new = mesh.add_vertex(1.5, 2.3, 0.0)
    >>> # Later: mesh.add_element([v_new, v1, v2])
  """
```

---

#### `remove_vertex(vert_id, strategy='error')`

```python
def remove_vertex(
    self,
    vert_id: int,
    strategy: str = 'error'  # 'error', 'collapse_edges', 'remove_orphaned'
) -> List[int]:
  """
  Remove vertex from mesh.
  
  Parameters:
    vert_id: Vertex to remove.
    strategy:
      - 'error': Raise error if vertex has adjacent elements (default, safest).
      - 'collapse_edges': Merge vertex into first neighbor (edge collapse in CAD).
      - 'remove_orphaned': Remove only if isolated (no adjacent elements/edges).
  
  Returns:
    List of element IDs that were deleted (due to edge collapse or orphaning).
  
  Raises:
    ValueError: If vertex has adjacent elements and strategy='error'.
    IndexError: If vert_id not in [0, n_verts).
  
  Notes:
    - Vertex IDs are NOT renumbered (sparse indexing). Gaps form in vertex list.
    - For collapse_edges: Merge all edges containing vert_id into neighbor.
      This deletes elements sharing edges with vert_id (topologically unsafe!).
      Use with caution; may violate manifold assumptions.
  """
```

---

### Edge Operations

#### `remove_edge(edge_id, strategy='error')`

```python
def remove_edge(
    self,
    edge_id: int,
    strategy: str = 'error'  # 'error', 'delete_elements', 'split_domain'
) -> List[int]:
  """
  Remove an edge from the mesh.
  
  Parameters:
    edge_id: Edge to remove.
    strategy:
      - 'error': Raise if edge is interior (two adjacent elements) (default).
      - 'delete_elements': Remove all elements incident to this edge.
      - 'split_domain': Mark edge as removed; may result in disconnected regions.
  
  Returns:
    List of element IDs deleted (empty if boundary edge or strategy='split_domain').
  
  Raises:
    ValueError: If edge is interior and strategy='error'.
    IndexError: If edge_id not in [0, n_edges).
  
  Example (domain splitting):
    >>> pinch_edges = mesh.pinch_points(threshold=0.1)
    >>> for edge in pinch_edges:
    ...   mesh.remove_edge(edge, strategy='split_domain')
    >>> components = mesh.connected_components()
  """
```

---

#### `swap_edge(edge_id)`

```python
def swap_edge(self, edge_id: int) -> bool:
  """
  Swap diagonal of a quad pair (flip shared edge).
  
  For a quad pair sharing edge (v1, v2), this flips the edge to (v3, v4)
  where [v1, v2, v3, v4] form a convex quadrilateral.
  
  Parameters:
    edge_id: Edge to swap.
  
  Returns:
    True if swap successful, False if edge is boundary or non-quad-pair.
  
  Raises:
    IndexError: If edge_id not in [0, n_edges).
  
  Notes:
    - Used for mesh quality improvement (angle-based edge swapping).
    - Only valid for quads; raises error if adjacent elements are triangles.
    - Preserves vertex counts, reorders connectivity.
  
  Example:
    >>> for edge_id in mesh.edges:
    ...   if mesh.swap_improves_quality(edge_id):
    ...     mesh.swap_edge(edge_id)
  """
```

---

### Query Operations (New)

#### `connected_components()`

```python
def connected_components(self) -> List[np.ndarray]:
  """
  Identify connected components of mesh.
  
  Returns:
    List of vertex ID arrays, one per component.
  
  Notes:
    - Useful for domain splitting (remove_edge can create components).
    - Components are sorted by size (largest first).
  
  Example:
    >>> mesh.remove_edge(pinch_edge, strategy='split_domain')
    >>> domains = mesh.connected_components()
    >>> assert len(domains) == 2  # Mesh split into 2 regions
    >>> mesh1 = mesh.subgraph(domains[0])
    >>> mesh2 = mesh.subgraph(domains[1])
  """
```

---

#### `pinch_points(threshold)`

```python
def pinch_points(self, threshold: float = 0.1) -> np.ndarray:
  """
  Identify bottleneck edges (narrow regions suitable for domain splitting).
  
  Parameters:
    threshold: Relative width threshold (0–1). Lower = narrower regions.
               Metric: frontier_width / average_mesh_width.
  
  Returns:
    Array of edge IDs representing pinch points.
  
  Notes:
    - Uses skeletonization result (layer frontier widths).
    - Intended for MADMESHR domain splitting (divide narrow regions).
    - Metric is heuristic; may need tuning for specific applications.
  
  Example:
    >>> pinch_edges = mesh.pinch_points(threshold=0.15)
    >>> # Suitable edges for domain splitting
  """
```

---

## Transactional Batch API (Future)

```python
class MeshTransaction:
  """Batch operations with rollback semantics."""
  
  def add_element(...) -> int: ...
  def remove_element(...) -> None: ...
  def add_vertex(...) -> int: ...
  # ... etc
  
  def commit(self) -> None:
    """Apply all changes, update adjacencies once."""
    
  def rollback(self) -> None:
    """Discard all changes."""

# Usage:
with mesh.transaction() as txn:
  for (vertices, elem_type) in elements_to_add:
    txn.add_element(vertices, elem_type)
  # If exception: rollback automatically
  # If success: commit with single adjacency rebuild
```

**Benefit:** Avoid O(n) adjacency updates per operation; batch to single O(n) rebuild.

---

## Incremental Skeletonization

### Approach 1: Full Rebuild (Simple, O(M^1.5))

```python
def add_element(self, vertices, elem_type):
  # ... add to connectivity_list, edges, adjacencies
  self._skeletonize()  # Full rebuild
```

**Pros:** Simple, guaranteed consistency.
**Cons:** Slow for MADMESHR (1000s of elements).

---

### Approach 2: Layer-Aware Invalidation (Complex, O(affected layers))

```python
def add_element(self, vertices, elem_type):
  elem_id = self._add_element_internal(vertices, elem_type)
  
  # Find affected layers (layers whose frontier might change)
  affected = self._find_affected_layers(elem_id)
  
  # Recompute only affected layers
  for layer_idx in affected:
    self._recompute_layer(layer_idx)
```

**Pros:** Faster for incremental updates (amortized O(frontier × iterations)).
**Cons:** Complex logic; easy to introduce bugs.

---

### Approach 3: Lazy Invalidation (Hybrid)

```python
def add_element(self, vertices, elem_type):
  elem_id = self._add_element_internal(vertices, elem_type)
  self._invalidate_layers = True  # Mark as dirty
  
def get_layer(self, layer_idx):
  if self._invalidate_layers:
    self._skeletonize()  # Rebuild only when queried
    self._invalidate_layers = False
  return self.layers[...]
```

**Pros:** Simple, lazy evaluation, good for MADMESHR (only check frontier when needed).
**Cons:** Non-deterministic timing (first access triggers rebuild).

---

## Recommended Approach

**Phase 2A (Add/Remove API):** Approach 1 (Full Rebuild)
- Simple to implement and test
- Adequate for initial MADMESHR integration (100–1000 elements)

**Phase 2B (Performance Optimization):** Approach 3 (Lazy Invalidation)
- Add `_invalidate_layers` flag
- Introduce `get_layer()` method (rebuilds if dirty)
- Minimal disruption to existing API
- Benchmark: Compare speed vs. full rebuild

**Phase 3 (Advanced):** Approach 2 (Layer-Aware) if profiling shows bottleneck
- Implement `_find_affected_layers()`, `_recompute_layer()`
- Requires careful testing (layer invariants)

---

## Testing Strategy

### Unit Tests

```python
def test_add_element_triangular():
  mesh = CHILmesh(...)
  v_new = mesh.add_element([1, 2, 3], elem_type='tri')
  assert v_new in mesh.connectivity_list
  assert mesh.elem_quality()[v_new] is valid
  
def test_add_element_quad():
  mesh = CHILmesh(...)
  v_new = mesh.add_element([1, 2, 3, 4], elem_type='quad')
  assert mesh.elem_quality()[v_new] is valid

def test_add_element_updates_adjacencies():
  # Verify Elem2Edge, Vert2Elem, Edge2Elem all updated correctly
  
def test_remove_element_consistency():
  # Verify all adjacencies cleaned up
  
def test_add_remove_roundtrip():
  # Add element → remove element → original state
  
def test_dynamic_ops_preserve_skeletonization_invariants():
  # Disjoint cover, monotone shrinking
```

### Integration Tests

```python
def test_advancing_front_simulation():
  # Simulate MADMESHR-like element placement
  mesh = CHILmesh(seed=boundary_loop)
  for _ in range(100):
    frontier_edges = mesh.boundary_edges()
    new_elem = select_element_from_frontier(frontier_edges)
    mesh.add_element(new_elem)
  # Verify mesh is valid (no overlaps, valid connectivity)
  
def test_domain_splitting():
  # Load mesh, find pinch point, split
  mesh = CHILmesh.read_from_fort14(...)
  pinch = mesh.pinch_points()[0]
  mesh.remove_edge(pinch, strategy='split_domain')
  domains = mesh.connected_components()
  assert len(domains) == 2
```

---

## Error Handling

### Exceptions

```python
class TopologyError(Exception):
  """Raised when operation violates mesh topology invariants."""
  
class DegenerateElementError(TopologyError):
  """Element has zero area or invalid geometry."""
  
class InconsistencyError(TopologyError):
  """Adjacency structures are inconsistent (internal bug)."""
```

### Validation Hooks

```python
def _validate_adjacencies(self):
  """Check invariants: edge2elem is symmetric, no duplicates, etc."""
  # Called after each operation in debug mode
  # Silent in production (control via environment variable)
```

---

## Performance Targets

| Operation | Target | Notes |
|-----------|--------|-------|
| `add_element()` | < 10 ms | Amortized over 1000 adds |
| `remove_element()` | < 5 ms | |
| Layer invalidation | O(1) | Lazy rebuild |
| Full layer rebuild | < 100 ms | Block_O in <1s total |
| Pinch-point detection | < 50 ms | Reuses skeletonization |

---

## Dependencies & Compatibility

- No external library changes required (numpy, scipy still sufficient)
- Backwards compatible: Old API (static mesh) unchanged
- New API opt-in (users call `add_element()` explicitly)

---

## References

- [Edge Collapse in Mesh Simplification](https://en.wikipedia.org/wiki/Mesh_simplification)
- [Advancing Front Meshing](https://en.wikipedia.org/wiki/Advancing_front_technique)
- MADMESHR documentation (when available)

---

**Status:** DRAFT - Ready for Phase 2 planning
**Next Review:** After Phase 1 graph modernization complete
