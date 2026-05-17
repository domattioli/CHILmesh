# CHILmesh Adjacency Structures Guide

## Overview

CHILmesh maintains multiple topology representations for efficient queries. This doc describes each structure, semantics, access patterns, and invariants.

## Mesh Representation

### Core Data

**points**: `ndarray[n_verts, 3]` — vertex coordinates (x, y, z); access: `mesh.points[vert_id]`

**connectivity_list**: `ndarray[n_elems, 3|4]` — element vertex indices; triangles: 3 unique + 1 repeated; quads: 4 vertices; mixed: consistent 4-column format; access: `mesh.connectivity_list[elem_id]`

### Adjacency Dictionary

All structures stored in `mesh.adjacencies`:

#### Dense Adjacencies (numpy arrays)

**Edge2Vert**: `ndarray[n_edges, 2]` — maps edge ID → `[v_min, v_max]` (canonical, always sorted); access: `v1, v2 = mesh.edge2vert(edge_id)`

**Elem2Edge**: `ndarray[n_elems, 3|4]` — maps element ID → edge IDs; triangles: 3 edges (col 3 unused); quads: 4 edges; all edge IDs valid [0, n_edges)

**Edge2Elem**: `ndarray[n_edges, 2]` — maps edge ID → `[elem1, elem2]`; elem2 = -1 for boundary; access: `elem1, elem2 = mesh.edge2elem(edge_id)`

#### Sparse Adjacencies (dicts + sets)

**Vert2Edge**: `Dict[int, Set[int]]` — vertex ID → incident edge IDs; invariant: if edge (v1,v2) exists, both vert2edge[v1] and vert2edge[v2] contain it; typical degree: 4-6 per vertex; access: `mesh.get_vertex_edges(vert_id)`

**Vert2Elem**: `Dict[int, Set[int]]` — vertex ID → incident element IDs; invariant: if elem contains v, then v ∈ vert2elem[v]; typical degree: 4-8 per vertex; access: `mesh.get_vertex_elements(vert_id)`

**EdgeMap**: Internal hash table for O(1) edge lookup by vertex pair — not public API; maps (v1,v2) → edge_id with automatic canonical form.

## Access Patterns

### Query vertex incident edges
```python
edges = mesh.get_vertex_edges(vert_id)  # Returns Set[int]
for edge_id in edges:
    v1, v2 = mesh.edge2vert(edge_id)
```

### Query vertex incident elements
```python
elems = mesh.get_vertex_elements(vert_id)  # Returns Set[int]
for elem_id in elems:
    verts = mesh.connectivity_list[elem_id]
```

### Find element neighbors
```python
edges = mesh.elem2edge(elem_id)
neighbors = set()
for edge_id in edges:
    elem1, elem2 = mesh.edge2elem(edge_id)
    if elem1 == elem_id:
        if elem2 >= 0: neighbors.add(elem2)
    else:
        neighbors.add(elem1)
```

### Find vertex neighbors
```python
edges = mesh.get_vertex_edges(vert_id)
neighbors = set()
for edge_id in edges:
    v1, v2 = mesh.edge2vert(edge_id)
    neighbors.add(v2 if v1 == vert_id else v1)
```

## Invariants

### Connectivity Validity
- All vertex indices ∈ [0, n_verts); edge indices ∈ [0, n_edges); element indices ∈ [0, n_elems)

### Edge Consistency
- Every edge in exactly one element's edge list
- Edge2Vert[e] = sorted([v1, v2]) (canonical form)
- No duplicate edges

### Adjacency Completeness
- All edges referenced by some element
- Vert2Edge: if v is in edge e, then e ∈ vert2edge[v]
- Vert2Elem: if v is in elem e, then e ∈ vert2elem[v]

### Orientation Consistency
- All elements CCW (enforced by `_ensure_ccw_orientation()`)
- Verified by positive signed area

### Skeletonization Partition (when layers computed)
- Layers form disjoint cover of all elements
- Layer sizes decrease monotonically
- No element missing or duplicated across layers

## Structure Comparison: Old vs New

| Aspect | Old (Phase 1) | New (Phase 2) |
|--------|---------------|---------------|
| Type | `List[List[int]]` | `Dict[int, Set[int]]` |
| Access | `vert2edge[v]` | `mesh.get_vertex_edges(v)` |
| Validation | Manual | Automatic via getter methods |
| Lookup | O(k) linear | O(1) dict lookup |
| Semantics | Implicit (may be missing) | Explicit (all vertices present) |

## Performance Characteristics

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| `edge2vert(e)` | O(1) | Direct array indexing |
| `elem2edge(e)` | O(1) | Direct array indexing |
| `edge2elem(e)` | O(1) | Direct array indexing |
| `get_vertex_edges(v)` | O(1) + O(k) copy | k = incident edges |
| `get_vertex_elements(v)` | O(1) + O(k) copy | k = incident elements |
| `_build_vert2edge()` | O(n_edges) | Linear |
| `_build_vert2elem()` | O(n_elems * avg_elem_size) | Linear |
| `_validate_adjacencies()` | O(n_verts + n_edges) | Full check |

## Debugging

### Validating Adjacencies
```python
try:
    mesh._validate_adjacencies()
    print("All adjacencies valid")
except AssertionError as e:
    print(f"Adjacency violation: {e}")
```

### Common Errors

**`AssertionError: Vert2Edge has X entries, expected Y`** — missing/extra vertices; fix: rebuild via `_build_adjacencies()`

**`AssertionError: Vert2Edge[V] is not a set`** — data corrupted (list instead of set); fix: check code modifying `adjacencies` dict directly

**`AssertionError: Invalid edge ID X in Vert2Edge[V]`** — edge ID out of range; fix: verify edge building

**`AssertionError: Vertex V not endpoint of edge E`** — wrong edge in vert2edge; fix: check canonical form enforcement

## References

Phase 2 Implementation:
- P2-01: Migrate Vert2Edge to explicit dict (#17)
- P2-02: Migrate Vert2Elem to explicit dict (#18)
- P2-03: Add getter methods for vertex adjacencies (#19)
- P2-04: Update traversal patterns (#20)
- P2-05: Add type hints and validation (#21)
- P2-06: Document adjacency structures (#22)
