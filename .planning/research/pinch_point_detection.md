# Pinch-Point Detection Algorithm

## Problem Statement

For MADMESHR advancing-front mesh generation, we need to identify **bottlenecks** (narrow regions) in a domain where the mesh width becomes critically small. These "pinch points" are ideal candidates for **domain splitting**: divide the domain into independent sub-regions, mesh each separately, and merge results.

**Example:** A domain shaped like a dumbbell (two large lobes connected by a thin neck). MADMESHR can:
1. Detect the pinch point (narrow neck)
2. Split the domain at that point
3. Mesh left and right lobes independently
4. Merge the two sub-meshes

---

## Use Case: MADMESHR Integration

```python
# Pseudocode: Domain-aware advancing-front
def adaptive_mesh_generation(domain_polygon, target_quality=0.9):
  """Generate mesh with adaptive domain splitting."""
  
  # Stage 1: Detect bottlenecks
  mesh = CHILmesh(domain_polygon)  # Coarse seed mesh
  pinch_edges = mesh.pinch_points(threshold=0.15)  # 15% relative width
  
  if pinch_edges:
    # Stage 2: Split domain
    for edge in pinch_edges:
      mesh.remove_edge(edge, strategy='split_domain')
    
    domains = mesh.connected_components()
    
    # Stage 3: Mesh each domain independently
    sub_meshes = []
    for domain_verts in domains:
      domain_mesh = adaptive_mesh_generation(domain_verts, target_quality)
      sub_meshes.append(domain_mesh)
    
    # Stage 4: Merge sub-meshes (align at shared boundaries)
    result = merge_meshes(sub_meshes)
    return result
  else:
    # No pinch points: mesh entire domain using advancing front
    agent = RL_Agent.load('advancing_front_policy.ckpt')
    while frontier.has_elements():
      action = agent.act(frontier_state())
      place_element(action)
    return mesh
```

---

## Candidate Algorithms

### Algorithm A: Skeletonization-Based (Recommended)

**Idea:** Reuse CHILmesh's skeletonization (medial axis extraction) to identify narrow layer frontiers.

```python
def pinch_points_via_skeleton(mesh, threshold=0.1):
  """
  Identify pinch points using medial axis layers.
  
  Parameters:
    mesh: CHILmesh instance (already skeletonized)
    threshold: Relative width threshold (0–1)
              pinch if frontier_width / avg_mesh_width < threshold
  
  Returns:
    Array of edge IDs representing bottlenecks.
  """
  
  skeleton = mesh.layers
  avg_width = estimate_mesh_width(mesh)  # See below
  
  pinch_edges = []
  for layer_idx in range(mesh.n_layers):
    frontier_edges = skeleton['bEdgeIDs'][layer_idx]
    frontier_width = measure_frontier_width(frontier_edges)
    
    if frontier_width < threshold * avg_width:
      pinch_edges.extend(frontier_edges)
  
  return np.array(pinch_edges)
```

**Pros:**
- Reuses existing `_skeletonize()` computation
- No new algorithms needed (framework already present)
- Elegant: pinch points are naturally narrow frontiers
- O(layers) complexity, which is O(√M) for 2D meshes

**Cons:**
- Frontier width metric is geometric (needs Voronoi analysis or proxy)
- Threshold tuning needed (domain-specific)

**Complexity:** O(√M × frontier size) = O(M) expected

---

### Algorithm B: Voronoi-Based Medial Axis

**Idea:** Compute Voronoi diagram, extract medial axis skeleton, identify narrow regions.

```python
from scipy.spatial import Voronoi

def pinch_points_via_voronoi(mesh, threshold=0.1):
  """Use Voronoi diagram to compute exact medial axis."""
  
  vor = Voronoi(mesh.points[:, :2])  # 2D points
  
  # Voronoi vertices → medial axis
  # Voronoi edges → topology
  # Voronoi ridge half-length → distance to boundary
  
  # Identify narrow Voronoi regions (small ridge lengths)
  pinch_verts = []
  for ridge_idx, ridge in enumerate(vor.ridge_vertices):
    if ridge_half_length(ridge) < threshold * avg_width:
      pinch_verts.append(ridge_idx)
  
  # Map back to mesh edges (intersections of Voronoi edges with boundary)
  pinch_edges = map_voronoi_to_mesh_edges(pinch_verts, mesh)
  return pinch_edges
```

**Pros:**
- Mathematically rigorous (exact medial axis)
- Width metric is well-defined (Voronoi ridge distance)

**Cons:**
- Voronoi computation O(n log n) (expensive for large meshes)
- Complex mapping back to mesh edges
- Overkill if approximate detection sufficient

**Complexity:** O(n log n + n × frontier size)

---

### Algorithm C: Distance Field–Based

**Idea:** Compute distance field from boundary, identify critical points (Hessian analysis).

```python
from scipy.ndimage import distance_transform_edt

def pinch_points_via_distance_field(mesh, threshold=0.1):
  """
  Rasterize mesh onto grid, compute distance field, identify bottlenecks.
  
  This is a grid-based approximation of the medial axis.
  """
  
  # Rasterize mesh to binary image (inside/outside)
  grid = rasterize_mesh(mesh, resolution=...)
  
  # Distance transform: each pixel = distance to boundary
  dist_field = distance_transform_edt(grid)
  
  # Identify ridges (local maxima in distance field)
  # These are candidates for pinch points
  ridges = identify_ridges(dist_field)
  
  # Find narrow ridges (distance < threshold × max_distance)
  narrow_ridges = [r for r in ridges if dist_field[r] < threshold * dist_field.max()]
  
  # Map back to mesh edges
  pinch_edges = map_grid_to_mesh_edges(narrow_ridges, mesh)
  return pinch_edges
```

**Pros:**
- Simple to implement (scipy provides distance_transform_edt)
- Good for quick approximation
- Scales well (linear in grid size, not mesh size)

**Cons:**
- Grid resolution affects accuracy (tradeoff: memory vs. precision)
- Requires rasterization (adds memory, may lose geometry details)
- Not suitable for non-convex domains (distance field can be misleading)

**Complexity:** O(grid^2) = O(domain area / resolution^2)

---

## Recommended: Algorithm A (Skeletonization-Based)

### Frontier Width Metric

**Option 1: Voronoi-Inspired Proxy**

Approximate the "width" of a frontier as the **minimum Euclidean distance** from any frontier edge's centroid to the mesh boundary:

```python
def frontier_width_euclidean(frontier_edges, mesh):
  """Width = distance from frontier to nearest boundary point."""
  
  boundary_edges = mesh.boundary_edges()
  boundary_points = []
  for edge in boundary_edges:
    v1, v2 = mesh.edge2vert(edge)
    p1, p2 = mesh.points[v1], mesh.points[v2]
    boundary_points.extend([p1, p2])
  
  width = float('inf')
  for edge in frontier_edges:
    v1, v2 = mesh.edge2vert(edge)
    centroid = (mesh.points[v1] + mesh.points[v2]) / 2
    
    # Distance from centroid to nearest boundary point
    dist = min(np.linalg.norm(centroid - p) for p in boundary_points)
    width = min(width, dist)
  
  return width
```

**Complexity:** O(|frontier| × |boundary|) = O(frontier × perimeter)

**Accuracy:** Good proxy for narrow regions; easy to compute.

---

**Option 2: Element Count Ratio**

Simpler heuristic: Frontier is narrow if the number of elements in the layer drops sharply:

```python
def frontier_width_element_ratio(layer_idx, mesh):
  """Width proxy = ratio of outer to inner element counts."""
  
  n_outer = len(mesh.layers['OE'][layer_idx])
  n_inner = len(mesh.layers['IE'][layer_idx])
  
  if n_inner == 0:
    return float('inf')  # Last layer, boundary
  
  ratio = n_outer / n_inner
  
  # Narrow frontier: many outer elements, few inner
  # (advancing front is bottlenecked)
  # Ratio > 2 suggests pinch point
  return ratio
```

**Complexity:** O(1) per layer = O(n_layers) = O(√M)

**Accuracy:** Heuristic; works well for convex domains, less reliable for complex shapes.

---

### Implementation (Phase 3)

```python
def pinch_points(
    self,
    threshold: float = 0.15,
    metric: str = 'euclidean'  # 'euclidean', 'element_ratio'
) -> np.ndarray:
  """
  Identify bottleneck edges for domain splitting.
  
  Parameters:
    threshold: Relative width threshold (0–1).
               If metric='euclidean': edges with width < threshold × avg_width.
               If metric='element_ratio': layers with ratio > 1/threshold.
    metric: Width computation method.
  
  Returns:
    Array of edge IDs representing pinch points.
  
  Notes:
    - Requires skeletonization (automatic if not done).
    - Intended for MADMESHR domain splitting.
    - Threshold is domain-specific; recommend experimenting with 0.1–0.3.
  
  Example:
    >>> mesh = CHILmesh.read_from_fort14('domain.14')
    >>> pinch_edges = mesh.pinch_points(threshold=0.2, metric='euclidean')
    >>> print(f"Found {len(pinch_edges)} pinch points")
  """
  
  if self.n_layers == 0:
    self._skeletonize()  # Compute layers if needed
  
  if metric == 'euclidean':
    return self._pinch_points_euclidean(threshold)
  elif metric == 'element_ratio':
    return self._pinch_points_element_ratio(threshold)
  else:
    raise ValueError(f"Unknown metric: {metric}")
```

---

## Testing Strategy

### Unit Tests

```python
def test_pinch_points_empty_mesh():
  mesh = CHILmesh(...)
  pinch = mesh.pinch_points()
  assert len(pinch) == 0

def test_pinch_points_dumbbell_domain():
  mesh = create_dumbbell_mesh()
  pinch = mesh.pinch_points(threshold=0.2)
  assert len(pinch) > 0

def test_pinch_points_annulus():
  mesh = CHILmesh.examples.annulus()
  pinch = mesh.pinch_points(threshold=0.2)
  assert len(pinch) == 0

def test_pinch_points_threshold_sensitivity():
  mesh = create_test_mesh()
  pinch_1 = mesh.pinch_points(threshold=0.1)
  pinch_2 = mesh.pinch_points(threshold=0.3)
  assert len(pinch_1) <= len(pinch_2)
```

---

## Performance Targets

| Operation | Target | Notes |
|-----------|--------|-------|
| `pinch_points()` (Euclidean) | < 100 ms | O(layers × frontier) ≈ O(M) |
| `pinch_points()` (Element ratio) | < 10 ms | O(layers) ≈ O(√M) |

---

**Status:** DRAFT - Ready for Phase 3 planning
**Next Review:** After Phase 1–2 completion
