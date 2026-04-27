# CHILmesh Access Interface (CAI) Specification

**Version:** 1.0  
**Effective:** CHILmesh 0.2.0+  
**Stability:** Guaranteed stable through v1.0  

---

## Executive Summary

The **CHILmesh Access Interface (CAI)** defines the stable, documented public API that downstream projects (MADMESHR, ADMESH, ADMESH-Domains) can depend on. All methods listed here will maintain their signatures, return types, and semantics through v1.0. Breaking changes require a major version bump and minimum 2-week advance notice.

---

## 1. Core Mesh Properties

### Basic Metadata

```python
# All read-only public attributes

mesh.points: np.ndarray                  # Shape: [n_verts, 3], dtype: float64
# Vertex coordinates (x, y, z)
# Example: mesh.points[0] → [x, y, z]

mesh.connectivity_list: np.ndarray      # Shape: [n_elems, 3|4], dtype: int32
# Element vertex indices
# Triangles: row has 3 unique vertices + 1 repeated
# Quads: row has 4 vertices
# Mixed: Both formats with consistent 4-column layout
# Example: mesh.connectivity_list[0] → [v0, v1, v2, v3]

mesh.n_verts: int                       # Number of vertices
mesh.n_elems: int                       # Number of elements
mesh.n_edges: int                       # Number of edges

mesh.type: str                          # Mesh type classification
# Values: 'Triangular', 'Quadrilateral', 'Mixed-Element'
# Determined from connectivity_list analysis

mesh.grid_name: str                     # Optional mesh name (from file or set manually)
mesh.boundary_condition: dict | None    # Optional boundary conditions (fort.14 format)
```

**Invariants:**
- `points.shape[0] == n_verts`
- `connectivity_list.shape[0] == n_elems`
- All vertex indices in connectivity_list ∈ [0, n_verts)
- `type` matches actual element composition
- `n_edges` computed from edge discovery (not user-settable)

---

## 2. Topology Query Methods

### Edge Queries

```python
def boundary_edges(self) -> np.ndarray:
    """
    Identify boundary edges of the mesh.
    
    Returns:
        ndarray[m, 1] where m = number of boundary edges
        Contains edge IDs that form the mesh boundary
    
    Complexity: O(n_edges)
    
    Example:
        >>> boundary = mesh.boundary_edges()
        >>> for edge_id in boundary:
        ...     v1, v2 = mesh.edge2vert(edge_id)
    """
```

```python
def edge2vert(self, edge_ids: Optional[Union[int, List[int], np.ndarray]]) -> np.ndarray:
    """
    Get vertex endpoints for edges.
    
    Parameters:
        edge_ids: Single edge ID, list of IDs, or None (all edges)
        
    Returns:
        If edge_ids is scalar: ndarray[2]
        If edge_ids is array: ndarray[len(edge_ids), 2]
        Vertices always in canonical form [v_min, v_max]
    
    Complexity: O(1) per edge queried
    
    Example:
        >>> v1, v2 = mesh.edge2vert(5)
        >>> all_edges = mesh.edge2vert(None)  # All edge endpoints
    """
```

### Element Queries

```python
def elem2edge(self, elem_ids: Optional[Union[int, List[int], np.ndarray]]) -> np.ndarray:
    """
    Get edge IDs for elements.
    
    Parameters:
        elem_ids: Single element ID, list of IDs, or None (all elements)
        
    Returns:
        If elem_ids is scalar: ndarray[3|4] for triangle|quad
        If elem_ids is array: ndarray[len(elem_ids), 3|4]
        Triangles padded to 4 columns (repeated vertex in col 3)
    
    Complexity: O(1) per element queried
    
    Example:
        >>> edges = mesh.elem2edge(0)  # 3 or 4 edge IDs
        >>> all_edges = mesh.elem2edge(None)  # All element edges
    """
```

### Vertex Adjacency Queries

```python
def get_vertex_edges(self, vert_id: int) -> Set[int]:
    """
    Get all edges incident to a vertex.
    
    Parameters:
        vert_id: Vertex index ∈ [0, n_verts)
        
    Returns:
        Set[int]: Edge IDs touching this vertex
        Empty set if vertex is isolated
    
    Raises:
        ValueError: If vert_id out of range
        
    Complexity: O(1) lookup + O(k) copy, where k = incident edges
    
    Example:
        >>> edges = mesh.get_vertex_edges(0)
        >>> for edge_id in edges:
        ...     v1, v2 = mesh.edge2vert(edge_id)
    """
```

```python
def get_vertex_elements(self, vert_id: int) -> Set[int]:
    """
    Get all elements incident to a vertex.
    
    Parameters:
        vert_id: Vertex index ∈ [0, n_verts)
        
    Returns:
        Set[int]: Element IDs containing this vertex
        Empty set if vertex is isolated
    
    Raises:
        ValueError: If vert_id out of range
        
    Complexity: O(1) lookup + O(k) copy, where k = incident elements
    
    Example:
        >>> elems = mesh.get_vertex_elements(0)
        >>> for elem_id in elems:
        ...     verts = mesh.connectivity_list[elem_id]
    """
```

### Edge-Element Queries

```python
def edge2elem(self, edge_ids: Optional[Union[int, List[int], np.ndarray]]) -> np.ndarray:
    """
    Get elements adjacent to edges.
    
    Parameters:
        edge_ids: Single edge ID, list of IDs, or None (all edges)
        
    Returns:
        If edge_ids is scalar: ndarray[2] = [elem1, elem2]
        If edge_ids is array: ndarray[len(edge_ids), 2]
        Boundary edges: [elem, -1]
        Interior edges: [elem1, elem2]
    
    Complexity: O(1) per edge queried
    
    Example:
        >>> e1, e2 = mesh.edge2elem(10)
        >>> neighbors = set([e1, e2]) - {-1}  # Skip -1 boundary marker
    """
```

---

## 3. Quality & Geometry Methods

```python
def signed_area(self, elem_ids: Optional[Union[int, List[int], np.ndarray]] = None) -> np.ndarray:
    """
    Compute signed area of elements.
    
    Positive for CCW-oriented elements, negative for CW.
    
    Parameters:
        elem_ids: Element IDs to query, or None (all elements)
        
    Returns:
        ndarray[len(elem_ids)] of float64
        
    Complexity: O(k) where k = number of elements queried
    
    Example:
        >>> areas = mesh.signed_area()
        >>> cw_elements = np.where(areas < 0)[0]
    """
```

```python
def interior_angles(self, elem_ids: Optional[Union[int, List[int], np.ndarray]] = None) -> np.ndarray:
    """
    Compute interior angles at each vertex in elements.
    
    Parameters:
        elem_ids: Element IDs to query, or None (all elements)
        
    Returns:
        If elem_ids is scalar: ndarray[3|4] angles in radians
        If elem_ids is array: ndarray[len(elem_ids), 3|4]
        
    Complexity: O(k) where k = number of elements queried
    
    Example:
        >>> angles = mesh.interior_angles(0)
        >>> min_angle = np.min(angles)
        >>> print(f"Min angle: {np.degrees(min_angle):.1f}°")
    """
```

```python
def elem_quality(self, elem_ids: Optional[Union[int, List[int], np.ndarray]] = None,
                 quality_type: str = 'skew') -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Compute element quality metrics.
    
    Parameters:
        elem_ids: Element IDs to query, or None (all elements)
        quality_type: 'skew' (default), other types TBD
        
    Returns:
        Tuple of:
        - quality: ndarray[len(elem_ids)] ∈ [0, 1], 1=perfect
        - angles: ndarray[len(elem_ids), 3|4] in radians
        - stats: Dict with keys 'mean', 'min', 'max', 'std'
        
    Complexity: O(k) where k = number of elements queried
    
    Example:
        >>> quality, angles, stats = mesh.elem_quality()
        >>> print(f"Mean quality: {stats['mean']:.3f}")
        >>> poor = np.where(quality < 0.3)[0]
    """
```

---

## 4. File I/O Methods

### Reading

```python
@staticmethod
def read_from_fort14(full_file_name: Path, compute_layers: bool = True) -> "CHILmesh":
    """
    Load mesh from ADCIRC .fort.14 file.
    
    Parameters:
        full_file_name: Path to .fort.14 file
        compute_layers: If True, compute skeletonization layers (default: True)
        
    Returns:
        CHILmesh instance
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format invalid
        
    Example:
        >>> mesh = CHILmesh.read_from_fort14(Path("domain.fort.14"))
    """
```

```python
@staticmethod
def read_from_2dm(full_file_name: Path, compute_layers: bool = True) -> "CHILmesh":
    """
    Load mesh from SMS .2dm file.
    
    Parameters:
        full_file_name: Path to .2dm file
        compute_layers: If True, compute skeletonization layers
        
    Returns:
        CHILmesh instance
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format invalid
        
    Example:
        >>> mesh = CHILmesh.read_from_2dm(Path("mesh.2dm"))
    """
```

```python
@staticmethod
def from_admesh_domain(record, compute_layers: bool = True) -> "CHILmesh":
    """
    Create CHILmesh from ADMESH-Domains catalog record.
    
    Duck-typed interface (zero dependencies on ADMESH-Domains).
    Record must have 'fort14_path' attribute.
    
    Parameters:
        record: ADMESH-Domains record object
        compute_layers: If True, compute skeletonization layers
        
    Returns:
        CHILmesh instance
        
    Example:
        >>> mesh = CHILmesh.from_admesh_domain(domain_record)
    """
```

### Writing

```python
def write_to_fort14(self, filename: str, grid_name: Optional[str] = None) -> bool:
    """
    Save mesh to ADCIRC .fort.14 file.
    
    Parameters:
        filename: Output file path
        grid_name: Optional grid name (uses self.grid_name if not provided)
        
    Returns:
        True if successful, False if write failed
        
    Example:
        >>> success = mesh.write_to_fort14("domain_modified.14")
    """
```

---

## 5. Utility Methods

```python
def copy(self) -> "CHILmesh":
    """
    Create a deep copy of the mesh.
    
    Returns:
        New CHILmesh instance with all data copied
        
    Complexity: O(n_verts + n_elems)
    
    Example:
        >>> mesh2 = mesh.copy()
        >>> mesh2.points[0] = [0, 0, 0]  # Won't affect original
    """
```

```python
def admesh_metadata(self) -> Dict[str, Any]:
    """
    Get ADMESH-Domains compatible metadata.
    
    Returns:
        Dict with keys:
        - 'node_count': int
        - 'element_count': int
        - 'element_type': str ('tri', 'quad', 'mixed')
        - 'bounding_box': Dict with 'min' and 'max'
        
    Example:
        >>> metadata = mesh.admesh_metadata()
        >>> print(f"Nodes: {metadata['node_count']}")
    """
```

---

## 6. Stability Guarantees

### Guaranteed Stable (until v2.0)

✅ **Will NOT change:**
- All method signatures listed above
- Return types (numpy arrays, Sets, Dicts maintain structure)
- Behavior and semantics of all methods
- Public attribute types (points, connectivity_list, n_verts, etc.)

✅ **Only deprecation path:**
- If a method must change, it will be marked `@deprecated` first
- Minimum 2-release deprecation period before removal
- Clear migration path provided

### NOT Guaranteed (Implementation Details)

❌ **May change in patch/minor releases:**
- Internal adjacency data structures
- Exact layer ordering in skeletonization
- Memory layout and caching strategies
- Skipping validation in specific cases
- Specific error message text (only error type guaranteed)
- Performance characteristics (within order of magnitude)
- Numerical precision beyond what tests verify

**Do NOT rely on implementation details in your code.**

### Breaking Change Policy

1. **Before v1.0**: Breaking changes only in minor versions (0.3.0, 0.4.0, etc.)
2. **From v1.0+**: Breaking changes only in major versions (2.0.0, 3.0.0, etc.)
3. **Notice**: Minimum 2 weeks via GitHub issue with tag `breaking-change`
4. **Migration**: Detailed migration guide provided
5. **Transition period**: Old API deprecated for at least 2 releases before removal

---

## 7. Common Usage Patterns

### Reading and Inspecting

```python
from chilmesh import CHILmesh
from pathlib import Path

# Load mesh
mesh = CHILmesh.read_from_fort14(Path("domain.fort.14"))

# Inspect properties
print(f"Type: {mesh.type}")
print(f"Vertices: {mesh.n_verts}, Elements: {mesh.n_elems}, Edges: {mesh.n_edges}")

# Access data
coords = mesh.points[0:10]          # First 10 vertices
elem_verts = mesh.connectivity_list[0]  # First element vertices
```

### Finding Element Neighbors

```python
def get_element_neighbors(mesh, elem_id):
    """Find all elements adjacent to a given element."""
    neighbors = set()
    edges = mesh.elem2edge(elem_id)
    
    for edge_id in edges:
        e1, e2 = mesh.edge2elem(edge_id)
        if e1 == elem_id and e2 >= 0:
            neighbors.add(e2)
        elif e2 == elem_id and e1 >= 0:
            neighbors.add(e1)
    
    return neighbors
```

### Finding Vertex Neighbors

```python
def get_vertex_neighbors(mesh, vert_id):
    """Find all vertices adjacent to a given vertex."""
    neighbors = set()
    edges = mesh.get_vertex_edges(vert_id)
    
    for edge_id in edges:
        v1, v2 = mesh.edge2vert(edge_id)
        other = v2 if v1 == vert_id else v1
        neighbors.add(other)
    
    return neighbors
```

### Quality Assessment

```python
import numpy as np

# Get quality metrics
quality, angles, stats = mesh.elem_quality()

# Find poor quality elements
threshold = 0.3
poor_elems = np.where(quality < threshold)[0]

print(f"Mean quality: {stats['mean']:.3f}")
print(f"Poor elements: {len(poor_elems)}")

# Check angles
min_angles = np.min(angles, axis=1)  # Minimum angle per element
min_angle = np.min(min_angles)
print(f"Minimum angle: {np.degrees(min_angle):.1f}°")
```

---

## 8. Error Handling

All methods use standard Python exceptions:

| Exception | When Raised | Example |
|-----------|------------|---------|
| `ValueError` | Invalid input bounds | `mesh.get_vertex_edges(mesh.n_verts)` |
| `FileNotFoundError` | File I/O fails | `CHILmesh.read_from_fort14("missing.14")` |
| `TypeError` | Wrong argument type | `mesh.signed_area("invalid")` |
| `AssertionError` | Internal validation fails | Data corruption (rare) |

**Always handle exceptions in downstream code.**

---

## 9. Version Information

```python
import chilmesh
print(chilmesh.__version__)      # e.g., "0.2.0"
```

Check version to decide on compatibility:
- If using CHILmesh < 0.2.0: Some methods may not exist
- If using CHILmesh >= 0.2.0: All CAI methods guaranteed available

---

## 10. Backward Compatibility Notes

### Updating from v0.1.x to 0.2.0+

All existing code continues to work. New features are additive:
- New methods added (get_vertex_edges, get_vertex_elements)
- Existing methods unchanged
- New public API available

**No code changes required**, but migration to new methods recommended for clarity.

---

## 11. Support and Feedback

- **Documentation**: See `docs/ADJACENCY_STRUCTURES.md` for detailed topology info
- **Examples**: See `examples/` for typical usage patterns
- **Issues**: Report bugs or request clarification at GitHub issues
- **Integration**: See `docs/DOWNSTREAM_MIGRATION_GUIDE.md` for MADMESHR/ADMESH/ADMESH-Domains

---

## Summary Table

| Feature | Signature | Guarantee | Typical Use |
|---------|-----------|-----------|------------|
| `mesh.points` | ndarray[n,3] | Stable | Access vertex coordinates |
| `mesh.connectivity_list` | ndarray[m,3\|4] | Stable | Access element definitions |
| `mesh.n_verts, n_elems, n_edges` | int | Stable | Query mesh size |
| `mesh.type` | str | Stable | Check element type |
| `mesh.boundary_edges()` | ndarray | Stable | Find boundary |
| `mesh.edge2vert()` | ndarray | Stable | Edge endpoints |
| `mesh.elem2edge()` | ndarray | Stable | Element edges |
| `mesh.edge2elem()` | ndarray | Stable | Edge neighbors |
| `mesh.get_vertex_edges()` | Set[int] | Stable | Vertex incident edges |
| `mesh.get_vertex_elements()` | Set[int] | Stable | Vertex incident elements |
| `mesh.elem_quality()` | Tuple | Stable | Quality metrics |
| `mesh.signed_area()` | ndarray | Stable | Element areas |
| `mesh.interior_angles()` | ndarray | Stable | Vertex angles |
| `CHILmesh.read_from_fort14()` | CHILmesh | Stable | Load ADCIRC format |
| `mesh.write_to_fort14()` | bool | Stable | Save ADCIRC format |
| `mesh.copy()` | CHILmesh | Stable | Clone mesh |
| `mesh.admesh_metadata()` | Dict | Stable | ADMESH-Domains metadata |

---

**Last Updated:** 2026-04-27  
**CAI Version:** 1.0
