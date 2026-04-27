# CHILmesh Adjacency Structures Guide

## Overview

CHILmesh maintains multiple representations of mesh topology to support efficient queries and algorithms. This document describes each adjacency structure, their semantics, access patterns, and invariants.

## Mesh Representation

### Core Data

**points**: `ndarray[n_verts, 3]`
- Vertex coordinates (x, y, z)
- Access: `mesh.points[vert_id]` → [x, y, z]

**connectivity_list**: `ndarray[n_elems, 3|4]`
- Element vertex indices
- Triangles: row has 3 unique vertices + 1 repeated
- Quads: row has 4 vertices
- Mixed: Both triangles and quads with consistent 4-column format
- Access: `mesh.connectivity_list[elem_id]` → [v0, v1, v2, v3]

### Adjacency Dictionary

All adjacency structures are stored in `mesh.adjacencies` dict:

#### Dense Adjacencies (numpy arrays)

**Edge2Vert**: `ndarray[n_edges, 2]`
- Maps edge ID to its vertex endpoints
- Format: `[v_min, v_max]` (canonical form for order-agnostic lookup)
- Invariant: Always sorted (v_min < v_max)
- Access: `v1, v2 = mesh.edge2vert(edge_id)`

**Elem2Edge**: `ndarray[n_elems, 3|4]`
- Maps element ID to its edge IDs
- Triangles: 3 edges (columns 0-2, column 3 unused)
- Quads: 4 edges (all columns used)
- Invariant: All edge IDs valid [0, n_edges)
- Access: `edge_ids = mesh.elem2edge(elem_id)`

**Edge2Elem**: `ndarray[n_edges, 2]`
- Maps edge ID to (up to) 2 adjacent elements
- Format: `[elem1, elem2]` where elem2 = -1 for boundary edges
- Invariant: elem indices valid [0, n_elems) or -1
- Access: `elem1, elem2 = mesh.edge2elem(edge_id)`

#### Sparse Adjacencies (dicts + sets)

**Vert2Edge**: `Dict[int, Set[int]]`
- Maps vertex ID to incident edge IDs
- Invariant: If edge (v1, v2) exists, both vert2edge[v1] and vert2edge[v2] contain that edge ID
- Typical degree: 4-6 edges per vertex (mesh-dependent)
- Access: `edges = mesh.get_vertex_edges(vert_id)` → Set[int]

**Vert2Elem**: `Dict[int, Set[int]]`
- Maps vertex ID to incident element IDs
- Invariant: If elem contains vertex v, then v ∈ vert2elem[v]
- Typical degree: 4-8 elements per vertex (mesh-dependent)
- Access: `elems = mesh.get_vertex_elements(vert_id)` → Set[int]

**EdgeMap**: Internal `EdgeMap` class
- Hash table for O(1) edge lookup by vertex pair
- Used internally during adjacency building
- Not part of public API
- Maps edge (v1, v2) → edge_id with automatic canonical form enforcement

## Access Patterns

### Query vertex incident edges

```python
# Get all edges touching a vertex
edges = mesh.get_vertex_edges(vert_id)  # Returns Set[int]

# Iterate and access edge endpoints
for edge_id in edges:
    v1, v2 = mesh.edge2vert(edge_id)
    # Process edge...
```

### Query vertex incident elements

```python
# Get all elements containing a vertex
elems = mesh.get_vertex_elements(vert_id)  # Returns Set[int]

# Iterate and access element vertices
for elem_id in elems:
    verts = mesh.connectivity_list[elem_id]
    # Process element...
```

### Find element neighbors

```python
# Get all neighbors of an element
elem_id = 0
edges = mesh.elem2edge(elem_id)
neighbors = set()

for edge_id in edges:
    elem1, elem2 = mesh.edge2elem(edge_id)
    if elem1 == elem_id:
        if elem2 >= 0:
            neighbors.add(elem2)
    else:
        neighbors.add(elem1)
```

### Find edge neighbors of a vertex

```python
# Get all vertices adjacent to a vertex
edges = mesh.get_vertex_edges(vert_id)
neighbors = set()

for edge_id in edges:
    v1, v2 = mesh.edge2vert(edge_id)
    other_vertex = v2 if v1 == vert_id else v1
    neighbors.add(other_vertex)
```

## Invariants

All adjacency structures must satisfy these invariants to maintain mesh validity:

### Connectivity Validity
- All vertex indices in `[0, n_verts)`
- All element indices in `[0, n_elems)`
- All edge indices in `[0, n_edges)`

### Edge Consistency
- Every edge appears in exactly one element's edge list
- Edge2Vert[e] = sorted([v1, v2]) (canonical form)
- No duplicate edges in mesh

### Adjacency Completeness
- All edges referenced by some element
- All vertices referenced by some element (may be isolated in theory)
- Vert2Edge: if v is in edge e, then e ∈ vert2edge[v]
- Vert2Elem: if v is in elem e, then e ∈ vert2elem[v]

### Orientation Consistency
- All elements counter-clockwise (CCW)
- Enforced by `_ensure_ccw_orientation()`
- Verified by positive signed area

### Skeletonization Partition (when layers computed)
- Layers form disjoint cover of all elements
- Layer sizes decrease monotonically (number of elements decreases)
- No element missing or duplicated across layers

## Structure Comparison: Old vs New

### Vert2Edge / Vert2Elem

| Aspect | Old (Phase 1) | New (Phase 2) |
|--------|---------------|---------------|
| Type | `List[List[int]]` | `Dict[int, Set[int]]` |
| Access | `vert2edge[v]` | `mesh.get_vertex_edges(v)` |
| Iteration | `for e in vert2edge[v]` | `for e in mesh.get_vertex_edges(v)` |
| Validation | Manual (user responsible) | Automatic via getter methods |
| Memory | O(n_verts * avg_degree) | Same |
| Lookup | O(k) linear search | O(1) dict lookup |
| Semantics | Implicit (vertices might be missing) | Explicit (all vertices present) |

## Performance Characteristics

### Access Times

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| `edge2vert(e)` | O(1) | Direct array indexing |
| `elem2edge(e)` | O(1) | Direct array indexing |
| `edge2elem(e)` | O(1) | Direct array indexing |
| `get_vertex_edges(v)` | O(1) + O(k) copy | k = incident edges |
| `get_vertex_elements(v)` | O(1) + O(k) copy | k = incident elements |

### Building Time

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| `_build_vert2edge()` | O(n_edges) | Linear in number of edges |
| `_build_vert2elem()` | O(n_elems * avg_elem_size) | Linear in connectivity |
| `_validate_adjacencies()` | O(n_verts + n_edges) | Full structure check |

## Debugging

### Validating Adjacencies

Call `_validate_adjacencies()` to check all invariants:

```python
try:
    mesh._validate_adjacencies()
    print("✓ All adjacencies valid")
except AssertionError as e:
    print(f"✗ Adjacency violation: {e}")
```

### Common Errors

**AssertionError: Vert2Edge has X entries, expected Y**
- Cause: Missing or extra vertices in adjacency dict
- Fix: Rebuild adjacencies via `_build_adjacencies()`

**AssertionError: Vert2Edge[V] is not a set**
- Cause: Data structure corrupted (contains list instead of set)
- Fix: Check code modifying `adjacencies` dict directly

**AssertionError: Invalid edge ID X in Vert2Edge[V]**
- Cause: Edge ID out of range
- Fix: Verify edge building didn't skip/corrupt edge numbering

**AssertionError: Vertex V not endpoint of edge E**
- Cause: Vert2Edge references wrong edge
- Fix: Check edge ordering (canonical form enforcement)

## Related Documentation

- `CHILmesh.get_vertex_edges()`: Public API for Vert2Edge queries
- `CHILmesh.get_vertex_elements()`: Public API for Vert2Elem queries
- `CHILmesh.edge2vert()`: Public API for Edge2Vert queries
- `CHILmesh.edge2elem()`: Public API for Edge2Elem queries
- `CHILmesh.elem2edge()`: Public API for Elem2Edge queries

## References

Phase 2 Implementation:
- P2-01: Migrate Vert2Edge to explicit dict structure (issue #17)
- P2-02: Migrate Vert2Elem to explicit dict structure (issue #18)
- P2-03: Add getter methods for vertex adjacencies (issue #19)
- P2-04: Update traversal patterns (issue #20)
- P2-05: Add type hints and validation (issue #21)
- P2-06: Document adjacency structures (issue #22)
