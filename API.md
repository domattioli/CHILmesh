# CHILmesh Public API Reference

**Version:** 0.2.0  
**Stability:** Stable (guaranteed through v1.0)  
**Access Interface:** CHILmesh Access Interface (CAI)

---

## Overview

CHILmesh provides a comprehensive API for 2D mesh manipulation, analysis, and generation. This document describes all public methods guaranteed stable through version 1.0.

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

Load mesh from ADCIRC `.fort.14` file format.

**Parameters:**
- `path`: Path to fort.14 file
- `compute_layers`: If False, skip skeletonization (fast bulk loading)

**Returns:** CHILmesh instance

**Example:**
```python
mesh = CHILmesh.read_from_fort14(Path("mesh.fort.14"))
mesh_fast = CHILmesh.read_from_fort14(Path("mesh.fort.14"), compute_layers=False)
```

---

#### `read_from_2dm(path: Path, compute_layers: bool = True) -> CHILmesh`
**Stability:** Stable | **Complexity:** O(n log n)

Load mesh from SMS `.2dm` file format.

**Parameters:**
- `path`: Path to .2dm file
- `compute_layers`: If False, skip skeletonization

**Returns:** CHILmesh instance

---

#### `write_to_fort14(filename: str, grid_name: str = "CHILmesh Grid") -> bool`
**Stability:** Stable | **Complexity:** O(n)

Write mesh to ADCIRC `.fort.14` format.

**Parameters:**
- `filename`: Output file path
- `grid_name`: Header string for the mesh

**Returns:** True on success, False on failure

**Example:**
```python
success = mesh.write_to_fort14("output.fort.14", grid_name="My Mesh")
```

---

#### `from_admesh_domain(record, compute_layers: bool = True) -> CHILmesh`
**Stability:** Stable | **Complexity:** O(n log n)

Load mesh from ADMESH-Domains catalog record (duck-typed).

**Parameters:**
- `record`: Catalog record with path or content
- `compute_layers`: If False, skip skeletonization

**Returns:** CHILmesh instance

---

#### `copy() -> CHILmesh`
**Stability:** Stable | **Complexity:** O(n)

Create a deep copy of the mesh.

**Returns:** Independent CHILmesh copy

---

### Adjacency & Topology Methods

#### `boundary_edges() -> np.ndarray`
**Stability:** Stable | **Complexity:** O(n_edges)

Identify boundary edges of the mesh.

**Returns:** Array of boundary edge IDs

**Example:**
```python
boundary = mesh.boundary_edges()
for edge_id in boundary:
    v1, v2 = mesh.edge2vert(edge_id)[0]
```

---

#### `get_vertex_edges(vert_id: int) -> Set[int]`
**Stability:** Stable | **Complexity:** O(k) where k = vertex degree

Get all edges incident to a vertex.

**Parameters:**
- `vert_id`: Vertex index [0, n_verts)

**Returns:** Set of edge IDs

**Raises:** ValueError if vert_id out of range

**Example:**
```python
edges = mesh.get_vertex_edges(5)
for edge_id in edges:
    v1, v2 = mesh.edge2vert(edge_id)[0]
```

---

#### `get_vertex_elements(vert_id: int) -> Set[int]`
**Stability:** Stable | **Complexity:** O(k) where k = vertex degree

Get all elements incident to a vertex.

**Parameters:**
- `vert_id`: Vertex index [0, n_verts)

**Returns:** Set of element IDs

**Raises:** ValueError if vert_id out of range

---

#### `edge2vert(edge_ids: Optional[int|List[int]]) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(m) where m = number of edges

Get vertex endpoints for specified edges.

**Parameters:**
- `edge_ids`: Edge index or list of indices (None for all)

**Returns:** Array of shape (m, 2) with vertex indices

**Example:**
```python
v1, v2 = mesh.edge2vert(5)[0]
all_verts = mesh.edge2vert()  # All edges
```

---

#### `elem2edge(elem_ids: Optional[int|List[int]]) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(m)

Get edge IDs for specified elements.

**Parameters:**
- `elem_ids`: Element index or list (None for all)

**Returns:** Array of shape (m, 3|4) with edge IDs

**Example:**
```python
edges = mesh.elem2edge(0)[0]  # Edges of element 0
```

---

#### `edge2elem(edge_ids: Optional[int|List[int]]) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(m)

Get element IDs adjacent to specified edges.

**Parameters:**
- `edge_ids`: Edge index or list (None for all)

**Returns:** Array of shape (m, 2) with element IDs
- Column 0: First adjacent element
- Column 1: Second adjacent element (or -1 for boundary)

**Example:**
```python
e1, e2 = mesh.edge2elem(5)[0]
```

---

### Analysis & Quality Methods

#### `elem_quality(elem_ids: Optional[List[int]] = None, quality_type: str = 'skew') -> Tuple[np.ndarray, np.ndarray, Dict]`
**Stability:** Stable | **Complexity:** O(m) where m = num elements

Compute element quality metrics.

**Parameters:**
- `elem_ids`: Element indices (None for all)
- `quality_type`: Quality metric ('skew' recommended)

**Returns:** Tuple of (quality_array, angles_array, stats_dict)
- `quality_array`: Quality values [0, 1] for each element
- `angles_array`: Interior angles in degrees
- `stats_dict`: {'mean', 'min', 'max', 'std'} statistics

**Example:**
```python
quality, angles, stats = mesh.elem_quality()
poor_count = np.sum(quality < 0.3)
```

---

#### `interior_angles(elem_ids: Optional[List[int]] = None) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(m)

Get interior angles for specified elements.

**Parameters:**
- `elem_ids`: Element indices (None for all)

**Returns:** Array of angles in degrees

---

#### `get_layer(layer_idx: int) -> Dict[str, np.ndarray]`
**Stability:** Stable | **Complexity:** O(1)

Retrieve a specific skeletonization layer.

**Parameters:**
- `layer_idx`: Layer index [0, n_layers)

**Returns:** Dict with keys 'OE' (outer edges), 'IE' (inner edges), 'OV' (outer verts), 'IV' (inner verts)

**Example:**
```python
layer = mesh.get_layer(0)  # First layer
outer_verts = layer['OV']
```

---

#### `admesh_metadata() -> Dict[str, Any]`
**Stability:** Stable | **Complexity:** O(1)

Get mesh metadata for ADMESH-Domains catalog.

**Returns:** Dict with:
- `node_count`: Number of vertices
- `element_count`: Number of elements
- `element_type`: 'tri', 'quad', or 'mixed'
- `bounding_box`: ((x_min, y_min), (x_max, y_max))
- `element_types`: List of element sizes present

**Example:**
```python
meta = mesh.admesh_metadata()
print(f"Type: {meta['element_type']}, Nodes: {meta['node_count']}")
```

---

### Mesh Modification Methods (v0.2.0+)

#### `add_advancing_front_element(vertices: List[int], elem_type: str = "tri") -> int`
**Stability:** Stable | **Complexity:** O(k) where k = element degree

Add element to mesh during advancing-front generation.

**Parameters:**
- `vertices`: Vertex indices
  - Triangle: [v0, v1, v2]
  - Quad: [v0, v1, v2, v3]
- `elem_type`: 'tri' or 'quad'

**Returns:** ID of newly added element

**Raises:** ValueError if vertices invalid

**Example:**
```python
new_elem_id = mesh.add_advancing_front_element([v1, v2, v3], 'tri')
```

---

#### `remove_boundary_loop(edge_ids: List[int]) -> None`
**Stability:** Stable | **Complexity:** O(m + n)

Remove boundary loop elements (residual closure).

**Parameters:**
- `edge_ids`: Boundary edge IDs to process

**Raises:** ValueError if edge_id out of range

**Example:**
```python
boundary = mesh.advancing_front_boundary_edges()[:4]
mesh.remove_boundary_loop(boundary)
```

---

#### `advancing_front_boundary_edges() -> List[int]`
**Stability:** Stable | **Complexity:** O(n_edges)

Get boundary edge IDs for advancing-front placement.

**Returns:** Sorted list of boundary edge IDs

**Example:**
```python
boundary = mesh.advancing_front_boundary_edges()
for edge_id in boundary[:10]:
    v1, v2 = mesh.edge2vert(edge_id)[0]
    # Place next element along this edge
```

---

#### `pinch_points(width_threshold: float = 0.5) -> List[int]`
**Stability:** Stable | **Complexity:** O(n_verts * k) where k = avg degree

Identify bottleneck vertices for domain splitting.

**Parameters:**
- `width_threshold`: Ratio threshold for pinch detection [0, 1]

**Returns:** Sorted list of pinch point vertex IDs

**Example:**
```python
pinches = mesh.pinch_points(width_threshold=0.3)
for v in pinches:
    # Consider splitting domain at this vertex
```

---

### Smoothing & Optimization Methods

#### `smooth_mesh(method: str, acknowledge_change: bool = False, **kwargs) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(iterations * n)

Smooth mesh using specified method.

**Parameters:**
- `method`: Smoothing type ('direct' or 'angle_based')
- `acknowledge_change`: Safety flag (must be True)
- `**kwargs`: Method-specific parameters

**Returns:** New point coordinates

---

#### `angle_based_smoother(angle_limit: float = 30.0) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(n)

Smooth mesh by enforcing minimum interior angle.

**Parameters:**
- `angle_limit`: Minimum interior angle in degrees

**Returns:** New point coordinates

---

#### `direct_smoother(kinf: float = 1e12) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(n)

Smooth mesh using direct method.

**Parameters:**
- `kinf`: Penalty parameter

**Returns:** New point coordinates

---

### Properties & Data Access

#### `n_verts: int` (read-only)
Number of vertices in the mesh.

#### `n_elems: int` (read-only)
Number of elements in the mesh.

#### `n_edges: int` (read-only)
Number of edges in the mesh.

#### `points: np.ndarray` (read-only)
Vertex coordinates, shape (n_verts, 3)
- Column 0: x coordinate
- Column 1: y coordinate
- Column 2: z coordinate (often 0 for 2D meshes)

#### `connectivity_list: np.ndarray` (read-only)
Element connectivity, shape (n_elems, 3|4)
- Triangles: [v0, v1, v2, v0] (4-column) or [v0, v1, v2] (3-column)
- Quads: [v0, v1, v2, v3]

#### `layers: Dict` (read-only)
Skeletonization layers indexed by layer number.

#### `grid_name: str` (read/write)
Grid name for fort.14 output.

#### `adjacencies: Dict` (read-only)
Low-level adjacency structures (internal use only)
- `Edge2Vert`: Edge endpoints
- `Elem2Edge`: Element edges
- `Edge2Elem`: Adjacent elements per edge
- `Vert2Edge`: Vertex incident edges (dict)
- `Vert2Elem`: Vertex incident elements (dict)

---

## Bridge Adapters (v0.2.0+)

### MeshAdapterForMADMESHR

Provides convenience methods for mesh adaptation research.

#### `get_element_neighbors(elem_id: int) -> Set[int]`
Get all elements adjacent to specified element.

#### `get_element_quality_neighborhood(elem_id: int) -> Dict`
Get quality metrics for element and neighbors.

#### `get_refinement_region(elem_ids: List[int], include_neighbors: bool = True) -> Set[int]`
Get refinement region from seed elements.

---

### MeshAdapterForADMESH

Provides quality assessment methods for adaptive mesh refinement.

#### `get_mesh_quality_report(quality_type: str = "skew") -> Dict`
Get comprehensive mesh quality report.

#### `get_element_angles_summary(elem_ids: Optional[List[int]] = None) -> Dict`
Get angle-based quality summary.

---

### MeshAdapterForADMESHDomains

Provides domain-level queries for decomposition.

#### `get_domain_boundaries() -> Dict[int, Set[int]]`
Get mesh boundaries as domain boundaries.

#### `get_mesh_connectivity_info() -> Dict`
Get high-level connectivity statistics.

---

## Module Functions

#### `write_fort14(filename: Path, points: np.ndarray, elements: np.ndarray, grid_name: str) -> bool`
**Stability:** Stable | **Complexity:** O(n)

Write mesh data directly to fort.14 format (without CHILmesh instance).

---

## Deprecations

### v0.1.1 → v0.2.0

The following changes are backward compatible but discouraged:

- ❌ Direct access to `adjacencies['Vert2Edge']` as list-of-lists
- ✅ Use `get_vertex_edges(v)` instead

- ❌ Direct access to `adjacencies['Vert2Elem']` as list-of-lists
- ✅ Use `get_vertex_elements(v)` instead

---

## Stability Guarantees

### Through v1.0

✅ All methods in this document guaranteed stable  
✅ Method signatures will not change  
✅ Return types will not change  
✅ Behavior will not change  

### Breaking Changes Policy

- Major version bump (0.2 → 1.0)
- 2-week deprecation notice before removal
- Clear migration guidance provided

---

## Error Handling

### Common Exceptions

| Exception | Cause | Resolution |
|-----------|-------|-----------|
| `ValueError` | Index out of range | Verify index is in valid range |
| `ValueError` | Invalid parameter | Check enum/type constraints |
| `RuntimeError` | Mesh state violation | Ensure mesh is valid after modification |
| `IOError` | File not found | Verify file path exists |

---

## Performance Notes

### O(1) Operations
- `get_vertex_edges(v)` - Dict lookup
- `get_vertex_elements(v)` - Dict lookup
- Individual `elem2edge()` / `edge2vert()` calls

### O(n) Operations
- `boundary_edges()` - Scan all edges
- `elem_quality()` - Compute for all elements
- `pinch_points()` - Analyze all vertices

### O(n log n) Operations
- `read_from_fort14()` - Sort edges during discovery
- `read_from_2dm()` - Build adjacencies

### O(n²) Avoided
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
