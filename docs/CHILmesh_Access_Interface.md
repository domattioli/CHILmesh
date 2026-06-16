# CHILmesh Access Interface (CAI) Specification

**Version:** 1.0  
**Effective:** CHILmesh 0.2.0+  
**Stability:** Guaranteed stable through v1.0  

---

## Executive Summary

**CAI** defines stable, documented public API that downstream projects (MADMESHR, ADMESH, Valence) can depend on. All methods listed here maintain signatures, return types, and semantics through v1.0. Breaking changes require major version bump + minimum 2-week advance notice.

---

## 1. Core Mesh Properties

```python
# All read-only public attributes

mesh.points: np.ndarray                  # Shape: [n_verts, 3], dtype: float64
mesh.connectivity_list: np.ndarray      # Shape: [n_elems, 3|4], dtype: int32
# Triangles: 3 unique vertices + 1 repeated; Quads: 4 vertices; Mixed: 4-column layout

mesh.n_verts: int                       # Number of vertices
mesh.n_elems: int                       # Number of elements
mesh.n_edges: int                       # Number of edges

mesh.type: str                          # 'Triangular', 'Quadrilateral', or 'Mixed-Element'
mesh.grid_name: str                     # Optional mesh name
mesh.boundary_condition: dict | None    # Optional boundary conditions (fort.14 format)
```

**Invariants:** `points.shape[0] == n_verts`; `connectivity_list.shape[0] == n_elems`; all vertex indices ∈ [0, n_verts); `type` matches actual composition

---

## 2. Topology Query Methods

```python
def boundary_edges(self) -> np.ndarray:
    """Identify boundary edges. Returns ndarray of edge IDs. O(n_edges)."""

def edge2vert(self, edge_ids: Optional[Union[int, List[int], np.ndarray]]) -> np.ndarray:
    """
    Get vertex endpoints for edges. Returns ndarray[2] (scalar) or ndarray[m, 2] (array).
    Vertices in canonical form [v_min, v_max]. O(1) per edge.
    
    >>> v1, v2 = mesh.edge2vert(5)
    >>> all_edges = mesh.edge2vert(None)  # All edge endpoints
    """

def elem2edge(self, elem_ids: Optional[Union[int, List[int], np.ndarray]]) -> np.ndarray:
    """
    Get edge IDs for elements. Returns ndarray[3|4] (scalar) or ndarray[m, 3|4] (array).
    Triangles padded to 4 columns (repeated vertex in col 3). O(1) per element.
    
    >>> edges = mesh.elem2edge(0)  # 3 or 4 edge IDs
    """

def get_vertex_edges(self, vert_id: int) -> Set[int]:
    """
    Get all edges incident to vertex. Returns Set[int]. O(1) + O(k).
    Raises ValueError if vert_id out of range.
    
    >>> edges = mesh.get_vertex_edges(0)
    """

def get_vertex_elements(self, vert_id: int) -> Set[int]:
    """
    Get all elements incident to vertex. Returns Set[int]. O(1) + O(k).
    Raises ValueError if vert_id out of range.
    """

def edge2elem(self, edge_ids: Optional[Union[int, List[int], np.ndarray]]) -> np.ndarray:
    """
    Get elements adjacent to edges. Returns ndarray[2] (scalar) or ndarray[m, 2] (array).
    Boundary edges: [elem, -1]. Interior: [elem1, elem2]. O(1) per edge.
    
    >>> e1, e2 = mesh.edge2elem(10)
    >>> neighbors = set([e1, e2]) - {-1}
    """
```

---

## 3. Quality & Geometry Methods

```python
def signed_area(self, elem_ids=None) -> np.ndarray:
    """
    Signed area of elements. Positive=CCW, negative=CW. O(k).
    
    >>> areas = mesh.signed_area()
    >>> cw_elements = np.where(areas < 0)[0]
    """

def interior_angles(self, elem_ids=None) -> np.ndarray:
    """
    Interior angles at each vertex. Returns radians. O(k).
    
    >>> angles = mesh.interior_angles(0)  # ndarray[3|4]
    """

def elem_quality(self, elem_ids=None, quality_type: str = 'skew') -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Element quality metrics. Returns (quality[0,1], angles, stats). O(k).
    stats keys: 'mean', 'min', 'max', 'std'.
    
    >>> quality, angles, stats = mesh.elem_quality()
    >>> print(f"Mean quality: {stats['mean']:.3f}")
    """
```

---

## 4. File I/O Methods

```python
@staticmethod
def read_from_fort14(full_file_name: Path, compute_layers: bool = True) -> "CHILmesh":
    """Load from ADCIRC .fort.14. Raises FileNotFoundError, ValueError."""

@staticmethod
def read_from_2dm(full_file_name: Path, compute_layers: bool = True) -> "CHILmesh":
    """Load from SMS .2dm. Raises FileNotFoundError, ValueError."""

@staticmethod
def from_admesh_domain(record, compute_layers: bool = True) -> "CHILmesh":
    """
    Create from Valence catalog record (duck-typed, zero imports).
    Record must have 'fort14_path' attribute.
    """

def write_to_fort14(self, filename: str, grid_name: Optional[str] = None) -> bool:
    """Save to ADCIRC .fort.14. Returns True on success."""
```

---

## 5. Utility Methods

```python
def copy(self) -> "CHILmesh":
    """Deep copy of mesh. O(n_verts + n_elems)."""

def admesh_metadata(self) -> Dict[str, Any]:
    """
    Valence compatible metadata.
    Returns: node_count, element_count, element_type ('tri'/'quad'/'mixed'), bounding_box.
    """
```

---

## 6. Stability Guarantees

**Guaranteed stable (until v2.0):**
- All method signatures listed above
- Return types (numpy arrays, Sets, Dicts maintain structure)
- Behavior and semantics
- Public attribute types

**NOT guaranteed (implementation details):**
- Internal adjacency data structures
- Exact layer ordering in skeletonization
- Memory layout and caching strategies
- Specific error message text
- Performance within order of magnitude

**Breaking change policy:**
- Before v1.0: breaking changes only in minor versions (0.3.0, 0.4.0)
- From v1.0+: breaking changes only in major versions (2.0.0)
- Notice: minimum 2 weeks via GitHub issue tagged `breaking-change`
- Transition: old API deprecated for at least 2 releases before removal

---

## 7. Common Usage Patterns

### Reading and Inspecting
```python
from chilmesh import CHILmesh
from pathlib import Path

mesh = CHILmesh.read_from_fort14(Path("domain.fort.14"))
print(f"Type: {mesh.type}")
print(f"Vertices: {mesh.n_verts}, Elements: {mesh.n_elems}")
coords = mesh.points[0:10]
```

### Finding Element Neighbors
```python
def get_element_neighbors(mesh, elem_id):
    neighbors = set()
    edges = mesh.elem2edge(elem_id)
    for edge_id in edges:
        e1, e2 = mesh.edge2elem(edge_id)
        if e1 == elem_id and e2 >= 0: neighbors.add(e2)
        elif e2 == elem_id and e1 >= 0: neighbors.add(e1)
    return neighbors
```

### Finding Vertex Neighbors
```python
def get_vertex_neighbors(mesh, vert_id):
    neighbors = set()
    for edge_id in mesh.get_vertex_edges(vert_id):
        v1, v2 = mesh.edge2vert(edge_id)
        neighbors.add(v2 if v1 == vert_id else v1)
    return neighbors
```

### Quality Assessment
```python
quality, angles, stats = mesh.elem_quality()
poor_elems = np.where(quality < 0.3)[0]
print(f"Mean quality: {stats['mean']:.3f}, Poor elements: {len(poor_elems)}")
```

---

## 8. Error Handling

| Exception | When Raised |
|-----------|------------|
| `ValueError` | Invalid input bounds |
| `FileNotFoundError` | File I/O fails |
| `TypeError` | Wrong argument type |
| `AssertionError` | Internal validation fails (rare) |

---

## 9. Backward Compatibility (v0.1.x → 0.2.0+)

All existing code continues to work. New features are additive (new methods added, existing methods unchanged). **No code changes required**, but migration to new methods recommended for clarity.

---

## Summary Table

| Feature | Signature | Guarantee |
|---------|-----------|-----------|
| `mesh.points` | ndarray[n,3] | Stable |
| `mesh.connectivity_list` | ndarray[m,3\|4] | Stable |
| `mesh.n_verts, n_elems, n_edges` | int | Stable |
| `mesh.type` | str | Stable |
| `mesh.boundary_edges()` | ndarray | Stable |
| `mesh.edge2vert()` | ndarray | Stable |
| `mesh.elem2edge()` | ndarray | Stable |
| `mesh.edge2elem()` | ndarray | Stable |
| `mesh.get_vertex_edges()` | Set[int] | Stable |
| `mesh.get_vertex_elements()` | Set[int] | Stable |
| `mesh.elem_quality()` | Tuple | Stable |
| `mesh.signed_area()` | ndarray | Stable |
| `mesh.interior_angles()` | ndarray | Stable |
| `CHILmesh.read_from_fort14()` | CHILmesh | Stable |
| `mesh.write_to_fort14()` | bool | Stable |
| `mesh.copy()` | CHILmesh | Stable |
| `mesh.admesh_metadata()` | Dict | Stable |

---

**Last Updated:** 2026-04-27  
**CAI Version:** 1.0
