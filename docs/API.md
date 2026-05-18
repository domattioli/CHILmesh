# CHILmesh Public API Reference

**Version:** 0.2.0  
**Stability:** Stable (guaranteed through v1.0)  
**Access Interface:** CHILmesh Access Interface (CAI)

---

## Overview

CHILmesh provides an API for 2D mesh processing, smoothing, and analysis on triangular, quadrilateral, and mixed-element meshes (with limited mutation primitives such as `insert_vertex` and the advancing-front addition path). All public methods below are guaranteed stable through v1.0.

**Import:**
```python
from chilmesh import CHILmesh, write_fort14
from chilmesh.bridge import (
    MeshAdapterForMADMESHR,
    MeshAdapterForADMESH,
    MeshAdapterForADMESHDomains,
)
```

---

## Core Class: CHILmesh

### File I/O Methods

#### `read_from_fort14(path: Path, compute_layers: bool = True) -> CHILmesh`
**Stability:** Stable | **Complexity:** O(n log n)

Load mesh from ADCIRC `.fort.14` file.

**Parameters:** `path` — fort.14 file path; `compute_layers` — if False, skip skeletonization (fast bulk loading)

```python
mesh = CHILmesh.read_from_fort14(Path("mesh.fort.14"))
mesh_fast = CHILmesh.read_from_fort14(Path("mesh.fort.14"), compute_layers=False)
```

---

#### `read_from_2dm(path: Path, compute_layers: bool = True) -> CHILmesh`
**Stability:** Stable | **Complexity:** O(n log n)

Load mesh from SMS `.2dm` file.

---

#### `write_to_fort14(filename: str, grid_name: str = "CHILmesh Grid") -> bool`
**Stability:** Stable | **Complexity:** O(n)

Write mesh to ADCIRC `.fort.14` format. Returns True on success.

```python
success = mesh.write_to_fort14("output.fort.14", grid_name="My Mesh")
```

---

#### `from_admesh_domain(record, compute_layers: bool = True) -> CHILmesh`
**Stability:** Stable | **Complexity:** O(n log n)

Load mesh from ADMESH-Domains catalog record (duck-typed).

---

#### `copy() -> CHILmesh`
**Stability:** Stable | **Complexity:** O(n)

Deep copy of mesh.

---

### Adjacency & Topology Methods

#### `boundary_edges() -> np.ndarray`
**Stability:** Stable | **Complexity:** O(n_edges)

Returns array of boundary edge IDs.

```python
boundary = mesh.boundary_edges()
for edge_id in boundary:
    v1, v2 = mesh.edge2vert(edge_id)[0]
```

---

#### `get_vertex_edges(vert_id: int) -> Set[int]`
**Stability:** Stable | **Complexity:** O(k) where k = vertex degree

Returns set of edge IDs incident to vertex. Raises `ValueError` if vert_id out of range.

---

#### `get_vertex_elements(vert_id: int) -> Set[int]`
**Stability:** Stable | **Complexity:** O(k)

Returns set of element IDs incident to vertex. Raises `ValueError` if vert_id out of range.

---

#### `edge2vert(edge_ids: Optional[int|List[int]]) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(m)

Returns array of shape (m, 2) with vertex indices for specified edges (None = all).

```python
v1, v2 = mesh.edge2vert(5)[0]
all_verts = mesh.edge2vert()
```

---

#### `elem2edge(elem_ids: Optional[int|List[int]]) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(m)

Returns array of shape (m, 3|4) with edge IDs for specified elements (None = all).

---

#### `edge2elem(edge_ids: Optional[int|List[int]]) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(m)

Returns array of shape (m, 2) — col 0: first adjacent element; col 1: second element (-1 for boundary).

---

### Analysis & Quality Methods

#### `elem_quality(elem_ids: Optional[List[int]] = None, quality_type: str = 'skew') -> Tuple[np.ndarray, np.ndarray, Dict]`
**Stability:** Stable | **Complexity:** O(m)

Returns `(quality_array, angles_array, stats_dict)` where stats_dict has keys 'mean', 'min', 'max', 'std'.

```python
quality, angles, stats = mesh.elem_quality()
poor_count = np.sum(quality < 0.3)
```

---

#### `interior_angles(elem_ids: Optional[List[int]] = None) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(m)

Returns array of angles in degrees for specified elements.

---

#### `get_layer(layer_idx: int) -> Dict[str, np.ndarray]`
**Stability:** Stable | **Complexity:** O(1)

Returns dict with keys 'OE' (outer edges), 'IE' (inner edges), 'OV' (outer verts), 'IV' (inner verts) for layer at `layer_idx`.

```python
layer = mesh.get_layer(0)
outer_verts = layer['OV']
```

---

#### `admesh_metadata() -> Dict[str, Any]`
**Stability:** Stable | **Complexity:** O(1)

Returns dict: `node_count`, `element_count`, `element_type` ('tri'/'quad'/'mixed'), `bounding_box` ((x_min, y_min), (x_max, y_max)), `element_types`.

---

### Mesh Modification Methods (v0.2.0+)

#### `add_advancing_front_element(vertices: List[int], elem_type: str = "tri") -> int`
**Stability:** Stable | **Complexity:** O(k)

Add element during advancing-front generation. Returns new element ID.

**Parameters:** `vertices` — [v0, v1, v2] for tri or [v0, v1, v2, v3] for quad; `elem_type` — 'tri' or 'quad'

```python
new_elem_id = mesh.add_advancing_front_element([v1, v2, v3], 'tri')
```

---

#### `remove_boundary_loop(edge_ids: List[int]) -> None`
**Stability:** Stable | **Complexity:** O(m + n)

Remove boundary loop elements (residual closure).

```python
boundary = mesh.advancing_front_boundary_edges()[:4]
mesh.remove_boundary_loop(boundary)
```

---

#### `advancing_front_boundary_edges() -> List[int]`
**Stability:** Stable | **Complexity:** O(n_edges)

Returns sorted list of boundary edge IDs for advancing-front placement.

```python
boundary = mesh.advancing_front_boundary_edges()
for edge_id in boundary[:10]:
    v1, v2 = mesh.edge2vert(edge_id)[0]
```

---

#### `pinch_points(width_threshold: float = 0.5) -> List[int]`
**Stability:** Stable | **Complexity:** O(n_verts * k)

Returns sorted list of bottleneck vertex IDs for domain splitting.

**Parameters:** `width_threshold` — ratio threshold [0, 1]

```python
pinches = mesh.pinch_points(width_threshold=0.3)
```

---

### Smoothing & Optimization Methods

#### `smooth_mesh(method: str, acknowledge_change: bool = False, **kwargs) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(iterations * n)

Smooth mesh using 'direct' or 'angle_based' method. `acknowledge_change` must be True.

---

#### `angle_based_smoother(angle_limit: float = 30.0) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(n)

Enforce minimum interior angle. Returns new point coordinates.

---

#### `direct_smoother(kinf: float = 1e12) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(n)

Direct smoothing method. Returns new point coordinates.

---

### Properties & Data Access

| Property | Type | Description |
|----------|------|-------------|
| `n_verts` | int (read-only) | Number of vertices |
| `n_elems` | int (read-only) | Number of elements |
| `n_edges` | int (read-only) | Number of edges |
| `points` | ndarray (n_verts, 3) (read-only) | Vertex coordinates (x, y, z) |
| `connectivity_list` | ndarray (n_elems, 3\|4) (read-only) | Element connectivity |
| `layers` | Dict (read-only) | Skeletonization layers |
| `grid_name` | str (read/write) | Grid name for fort.14 output |
| `adjacencies` | Dict (read-only) | Low-level structures (internal use) |

---

## Bridge Adapters (v0.2.0+)

### MeshAdapterForMADMESHR

| Method | Description |
|--------|-------------|
| `get_element_neighbors(elem_id: int) -> Set[int]` | All elements adjacent to specified element |
| `get_element_quality_neighborhood(elem_id: int) -> Dict` | Quality metrics for element + neighbors |
| `get_refinement_region(elem_ids: List[int], include_neighbors: bool = True) -> Set[int]` | Refinement region from seed elements |

### MeshAdapterForADMESH

| Method | Description |
|--------|-------------|
| `get_mesh_quality_report(quality_type: str = "skew") -> Dict` | Comprehensive mesh quality report |
| `get_element_angles_summary(elem_ids: Optional[List[int]] = None) -> Dict` | Angle-based quality summary |

### MeshAdapterForADMESHDomains

| Method | Description |
|--------|-------------|
| `get_domain_boundaries() -> Dict[int, Set[int]]` | Mesh boundaries as domain boundaries |
| `get_mesh_connectivity_info() -> Dict` | High-level connectivity statistics |

---

## Module Functions

#### `write_fort14(filename: Path, points: np.ndarray, elements: np.ndarray, grid_name: str) -> bool`
Write mesh data directly to fort.14 (without CHILmesh instance). O(n).

---

## Deprecations (v0.1.1 → v0.2.0)

- ❌ `adjacencies['Vert2Edge']` as list-of-lists → ✅ `get_vertex_edges(v)`
- ❌ `adjacencies['Vert2Elem']` as list-of-lists → ✅ `get_vertex_elements(v)`

---

## Stability Guarantees

✅ All methods in this document stable through v1.0  
✅ Signatures will not change  
✅ Return types will not change  
✅ Behavior will not change  

**Breaking changes policy:** major version bump + 2-week deprecation notice + migration guidance

---

## Error Handling

| Exception | Cause | Resolution |
|-----------|-------|-----------|
| `ValueError` | Index out of range | Verify index in valid range |
| `ValueError` | Invalid parameter | Check enum/type constraints |
| `RuntimeError` | Mesh state violation | Ensure mesh valid after modification |
| `IOError` | File not found | Verify file path exists |

---

## Performance Notes

**O(1):** `get_vertex_edges(v)`, `get_vertex_elements(v)`, individual `elem2edge()`/`edge2vert()` calls

**O(n):** `boundary_edges()`, `elem_quality()`, `pinch_points()`

**O(n log n):** `read_from_fort14()`, `read_from_2dm()`

✅ Phase 1-3 optimizations eliminated O(n²) behaviors from v0.1.1

---

## Examples

### Load and Analyze
```python
from chilmesh import CHILmesh
from pathlib import Path

mesh = CHILmesh.read_from_fort14(Path("mesh.fort.14"))
print(f"Vertices: {mesh.n_verts}, Elements: {mesh.n_elems}")

quality, angles, stats = mesh.elem_quality()
print(f"Mean quality: {stats['mean']:.3f}")
```

### Use Bridge Adapter
```python
from chilmesh.bridge import MeshAdapterForADMESH

adapter = MeshAdapterForADMESH(mesh)
report = adapter.get_mesh_quality_report()
print(f"Poor elements: {report['poor_count']}")
```

### Advancing-Front Generation
```python
boundary = mesh.advancing_front_boundary_edges()
while len(boundary) > 4:
    v1, v2 = mesh.edge2vert(boundary[0])[0]
    v3 = mesh.n_verts  # New vertex
    new_elem = mesh.add_advancing_front_element([v1, v2, v3], 'tri')
    boundary = mesh.advancing_front_boundary_edges()
```

---

**Last Updated:** 2026-04-27  
**API Version:** 1.0
