# Phase 5.3: Spatial Indexing for Point Queries

**Goal:** Fast point-location + nearest-neighbor queries (O(log n) instead of O(n))  
**Impact:** Observation ingestion, mesh-to-mesh interpolation, particle tracking  
**Performance:** 100k point lookups from hours → seconds  
**Effort:** Design: 1 hour. Implementation: 4-6 hours

---

## Problem

Current: Brute-force scan every element for point-in-element.  
Example: WNAT mesh (98k elements) + 100k point lookups = minutes

---

## API

### find_element(point: ndarray) → int
Return containing element ID.

```python
m = CHILmesh(...)
elem_id = m.find_element(point=np.array([0.5, 0.5]))
```

### find_elements_in_radius(point: ndarray, radius: float) → ndarray
All elements within radius.

```python
nearby = m.find_elements_in_radius(point=np.array([0.5, 0.5]), radius=0.1)
```

### nearest_vertices(point: ndarray, k: int) → ndarray
k nearest mesh vertices.

```python
neighbors = m.nearest_vertices(point=np.array([0.5, 0.5]), k=3)
```

### nearest_elements(point: ndarray, k: int) → ndarray
k nearest by centroid.

```python
nearby_elems = m.nearest_elements(point=np.array([0.5, 0.5]), k=5)
```

---

## Implementation

**Backend:** `scipy.spatial.cKDTree` (fast C implementation)

```python
class CHILmesh:
    def __init__(self, ...):
        # ... existing code ...
        # Build spatial index at init (lazy option available)
        self._vertex_tree = cKDTree(self.points)
        self._centroid_tree = cKDTree(self.centroids)
```

**Point-in-element:** KD-tree narrows candidates → barycentric check

---

## Success Criteria

- [ ] `find_element()` implemented
- [ ] `find_elements_in_radius()` implemented
- [ ] `nearest_vertices()` and `nearest_elements()` implemented
- [ ] Performance: 100k lookups < 5 seconds (vs. hours brute-force)
- [ ] Build cost: <0.5s on WNAT (98k elements)
- [ ] Edge cases: points outside mesh handled (return -1 or nearest)

---

## Use Cases

- NOAA gauge ingestion → find mesh element for observation
- Mesh-to-mesh interpolation → find nearest source vertices
- Particle tracking → locate drifter at each timestep
- Local refinement → all elements within feature region

---

## Related

- Issue #92: Spatial indexing (this spec)
- Issue #94, #93: Mesh mutation API, incremental skeletonization
- MADMESHR: Point observations for data assimilation
- Constitution IV: Correctness (edge cases must be robust)

---

