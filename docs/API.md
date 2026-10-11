# CHILmesh Public API Reference

CHILmesh processes, smooths, and analyzes 2D triangular, quadrilateral, and mixed-element meshes. Every public method in this file is stable through v1.0. This file lists the methods first, then the capabilities, the mesh registry, and the command-line tools.

**Version:** 0.2.0  
**Stability:** Stable (guaranteed through v1.0)  
**Access Interface:** CHILmesh Access Interface (CAI)  
**Last Updated:** 2026-04-27  
**API Version:** 1.0

## Contents

1. [Import](#import)
2. [Core class: CHILmesh](#core-class-chilmesh)
3. [Bridge adapters](#bridge-adapters-v020)
4. [Module functions](#module-functions)
5. [Stability, deprecations, and errors](#stability-guarantees)
6. [Performance notes](#performance-notes)
7. [Examples](#examples)
8. [Capabilities overview](#capabilities-overview)
9. [Fetching meshes](#fetching-meshes-with-chilmeshregistry)
10. [Command-line tools and examples](#examples-and-cli-moved-from-the-readme)
11. [Roadmap](#roadmap-items-moved-from-the-readme)

---

## Import

The API covers mutation only through limited primitives, such as `insert_vertex` and the advancing-front addition path.

```python
from chilmesh import CHILmesh, write_fort14
from chilmesh.bridge import (
    MeshAdapterForADMESH,
    MeshAdapterForADMESHDomains,
)
```

---

## Core class: CHILmesh

### File I/O methods

#### `read_from_fort14(path: Path, compute_layers: bool = True) -> CHILmesh`
**Stability:** Stable | **Complexity:** O(n log n)

Loads a mesh from an ADCIRC `.fort.14` file.

**Parameters:**

- `path`: the fort.14 file path.
- `compute_layers`: if False, the method skips skeletonization (fast bulk loading).

```python
mesh = CHILmesh.read_from_fort14(Path("mesh.fort.14"))
mesh_fast = CHILmesh.read_from_fort14(Path("mesh.fort.14"), compute_layers=False)
```

---

#### `read_from_2dm(path: Path, compute_layers: bool = True) -> CHILmesh`
**Stability:** Stable | **Complexity:** O(n log n)

Loads a mesh from an SMS `.2dm` file.

---

#### `write_to_fort14(filename: str, grid_name: str = "CHILmesh Grid") -> bool`
**Stability:** Stable | **Complexity:** O(n)

Writes the mesh in ADCIRC `.fort.14` format. Returns True on success.

```python
success = mesh.write_to_fort14("output.fort.14", grid_name="My Mesh")
```

---

#### `from_admesh_domain(record, compute_layers: bool = True) -> CHILmesh`
**Stability:** Stable | **Complexity:** O(n log n)

Loads a mesh from a Valence catalog record (duck-typed).

---

#### `copy() -> CHILmesh`
**Stability:** Stable | **Complexity:** O(n)

Returns a deep copy of the mesh.

---

### Adjacency and topology methods

#### `boundary_edges() -> np.ndarray`
**Stability:** Stable | **Complexity:** O(n_edges)

Returns an array of boundary edge IDs.

```python
boundary = mesh.boundary_edges()
for edge_id in boundary:
    v1, v2 = mesh.edge2vert(edge_id)[0]
```

---

#### `get_vertex_edges(vert_id: int) -> Set[int]`
**Stability:** Stable | **Complexity:** O(k) where k = vertex degree

Returns the set of edge IDs incident to a vertex. Raises `ValueError` if vert_id is out of range.

---

#### `get_vertex_elements(vert_id: int) -> Set[int]`
**Stability:** Stable | **Complexity:** O(k)

Returns the set of element IDs incident to a vertex. Raises `ValueError` if vert_id is out of range.

---

#### `edge2vert(edge_ids: Optional[int|List[int]]) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(m)

Returns an array of shape (m, 2) with the vertex indices of the specified edges (None = all).

```python
v1, v2 = mesh.edge2vert(5)[0]
all_verts = mesh.edge2vert()
```

---

#### `elem2edge(elem_ids: Optional[int|List[int]]) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(m)

Returns an array of shape (m, 3|4) with the edge IDs of the specified elements (None = all).

---

#### `edge2elem(edge_ids: Optional[int|List[int]]) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(m)

Returns an array of shape (m, 2). Column 0 holds the first adjacent element. Column 1 holds the second element (-1 for a boundary edge).

---

### Analysis and quality methods

#### `elem_quality(elem_ids: Optional[List[int]] = None, quality_type: str = 'skew') -> Tuple[np.ndarray, np.ndarray, Dict]`
**Stability:** Stable | **Complexity:** O(m)

Returns `(quality_array, angles_array, stats_dict)`. The stats_dict has the keys 'mean', 'min', 'max', 'std'.

```python
quality, angles, stats = mesh.elem_quality()
poor_count = np.sum(quality < 0.3)
```

---

#### `interior_angles(elem_ids: Optional[List[int]] = None) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(m)

Returns an array of angles in degrees for the specified elements.

---

#### `get_layer(layer_idx: int) -> Dict[str, np.ndarray]`
**Stability:** Stable | **Complexity:** O(1)

Returns a dict for the layer at `layer_idx`, with the keys 'OE' (outer edges), 'IE' (inner edges), 'OV' (outer verts), and 'IV' (inner verts).

```python
layer = mesh.get_layer(0)
outer_verts = layer['OV']
```

---

#### `admesh_metadata() -> Dict[str, Any]`
**Stability:** Stable | **Complexity:** O(1)

Returns a dict with `node_count`, `element_count`, `element_type` ('tri'/'quad'/'mixed'), `bounding_box` ((x_min, y_min), (x_max, y_max)), and `element_types`.

---

### Mesh modification methods (v0.2.0+)

#### `add_advancing_front_element(vertices: List[int], elem_type: str = "tri") -> int`
**Stability:** Stable | **Complexity:** O(k)

Adds an element during advancing-front generation. Returns the new element ID.

**Parameters:**

- `vertices`: [v0, v1, v2] for a tri or [v0, v1, v2, v3] for a quad.
- `elem_type`: 'tri' or 'quad'.

```python
new_elem_id = mesh.add_advancing_front_element([v1, v2, v3], 'tri')
```

---

#### `remove_boundary_loop(edge_ids: List[int]) -> None`
**Stability:** Stable | **Complexity:** O(m + n)

Removes the elements of a boundary loop (residual closure).

```python
boundary = mesh.advancing_front_boundary_edges()[:4]
mesh.remove_boundary_loop(boundary)
```

---

#### `advancing_front_boundary_edges() -> List[int]`
**Stability:** Stable | **Complexity:** O(n_edges)

Returns a sorted list of boundary edge IDs for advancing-front placement.

```python
boundary = mesh.advancing_front_boundary_edges()
for edge_id in boundary[:10]:
    v1, v2 = mesh.edge2vert(edge_id)[0]
```

---

#### `pinch_points(width_threshold: float = 0.5) -> List[int]`
**Stability:** Stable | **Complexity:** O(n_verts * k)

Returns a sorted list of bottleneck vertex IDs for domain splitting.

**Parameters:**

- `width_threshold`: the ratio threshold, in [0, 1].

```python
pinches = mesh.pinch_points(width_threshold=0.3)
```

---

### Smoothing and optimization methods

#### `smooth_mesh(method: str, acknowledge_change: bool = False, **kwargs) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(iterations * n)

Smooths the mesh with the 'direct' or 'angle_based' method. `acknowledge_change` must be True.

---

#### `angle_based_smoother(angle_limit: float = 30.0) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(n)

Enforces a minimum interior angle. Returns the new point coordinates.

---

#### `direct_smoother(kinf: float = 1e12) -> np.ndarray`
**Stability:** Stable | **Complexity:** O(n)

Applies the direct smoothing method. Returns the new point coordinates.

---

### Properties and data access

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

## Bridge adapters (v0.2.0+)

### Legacy research adapter

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

## Module functions

#### `write_fort14(filename: Path, points: np.ndarray, elements: np.ndarray, grid_name: str) -> bool`

Writes mesh data directly to fort.14, without a CHILmesh instance. Complexity: O(n).

---

## Stability guarantees

Every method in this document is stable through v1.0:

- Signatures do not change.
- Return types do not change.
- Behavior does not change.

**Breaking changes policy:** a major version bump, a 2-week deprecation notice, and migration guidance.

### Deprecations (v0.1.1 to v0.2.0)

| Deprecated | Use instead |
|------------|-------------|
| `adjacencies['Vert2Edge']` as list-of-lists | `get_vertex_edges(v)` |
| `adjacencies['Vert2Elem']` as list-of-lists | `get_vertex_elements(v)` |

### Error handling

| Exception | Cause | Resolution |
|-----------|-------|-----------|
| `ValueError` | Index out of range | Verify index in valid range |
| `ValueError` | Invalid parameter | Check enum/type constraints |
| `RuntimeError` | Mesh state violation | Ensure mesh valid after modification |
| `IOError` | File not found | Verify file path exists |

---

## Performance notes

Phase 1-3 optimizations eliminated the O(n²) behaviors of v0.1.1. Current complexity classes:

- **O(1):** `get_vertex_edges(v)`, `get_vertex_elements(v)`, individual `elem2edge()` and `edge2vert()` calls.
- **O(n):** `boundary_edges()`, `elem_quality()`, `pinch_points()`.
- **O(n log n):** `read_from_fort14()`, `read_from_2dm()`.

---

## Examples

### Load and analyze

```python
from chilmesh import CHILmesh
from pathlib import Path

mesh = CHILmesh.read_from_fort14(Path("mesh.fort.14"))
print(f"Vertices: {mesh.n_verts}, Elements: {mesh.n_elems}")

quality, angles, stats = mesh.elem_quality()
print(f"Mean quality: {stats['mean']:.3f}")
```

### Use a bridge adapter

```python
from chilmesh.bridge import MeshAdapterForADMESH

adapter = MeshAdapterForADMESH(mesh)
report = adapter.get_mesh_quality_report()
print(f"Poor elements: {report['poor_count']}")
```

### Advancing-front generation

```python
boundary = mesh.advancing_front_boundary_edges()
while len(boundary) > 4:
    v1, v2 = mesh.edge2vert(boundary[0])[0]
    v3 = mesh.n_verts  # New vertex
    new_elem = mesh.add_advancing_front_element([v1, v2, v3], 'tri')
    boundary = mesh.advancing_front_boundary_edges()
```

---

## Capabilities overview

| Area | What ships | Entry points |
|---|---|---|
| File formats | ADCIRC fort.14 (lossless round trip, boundary records kept), fort.13 nodal attributes, fort.15 control file; SMS Aquaveo .2dm; Gmsh .msh 2.2 and 4.1; header-only `summary()` for fort.14, .2dm, fort.13, fort.15, .msh, .npy and .npz without loading the mesh | `Mesh.read_from_fort14`, `read_fort13`, `read_fort15`, `Mesh.read_from_2dm`, `Mesh.read_from_msh`, `chilmesh.summary` |
| Topology | Triangular, quadrilateral and mixed meshes on one API; seven adjacency tables; concentric layer peel with per-layer `OE`, `IE`, `OV`, `IV` and boundary-edge IDs | `Mesh.peel_layers`, `mesh.layers` |
| Quality | Skew quality and interior angles per element; per-edge Courant number and a CFL gate from depths and a time step | `Mesh.elem_quality`, `courant_number`, `cfl_gate` |
| Smoothing | Balendran direct FEM; Zhou-Shimada angle-based; ADMESH spring-based truss against a signed distance function; plain Laplacian. FEM, angle-based and truss accept an optional target size function. Boundary nodes fixed, topology unchanged, mixed meshes accepted | `Mesh.smooth_mesh(method=...)`, `Mesh.smooth` |
| Mutation | 13 topology operations: split triangle, split edge, split boundary edge, swap edge, merge elements, remove vertex, collapse edge, insert vertex, move boundary node, split triangles, topological smoothing, local re-peel, layer diff | `chilmesh.mutations` |
| Boundary editing | Advancing-front element addition, boundary-loop removal, pinch-point detection in narrow channels | `Mesh.add_advancing_front_element`, `Mesh.remove_boundary_loop`, `Mesh.pinch_points` |
| Spatial queries | Point-in-element, radius search, k-nearest vertices, all O(log n) after a one-time index | `Mesh.find_element`, `Mesh.find_elements_in_radius`, `Mesh.nearest_vertices` |
| Geometry | Haversine distance, convex hull, Hausdorff distance, antimeridian-aware bounding boxes | `chilmesh.geometry` |
| Interop | ADMESH domain adapter and truss warm start; renumbering-tolerant node matching between two meshes with nodal-field deltas | `Mesh.from_admesh_domain`, `chilmesh.admesh_warmstart`, `match_nodes` |
| Plotting | Mesh, per-element quality map, matched-colormap quality histogram, layer paths; experimental GPU rendering through the `[gpu]` extra | `Mesh.plot`, `Mesh.plot_quality`, `Mesh.plot_quality_histogram`, `chilplotting.render_image(backend="gpu")` |

---

## Fetching meshes with `chilmesh.registry`

Meshes come through Valence. The `chilmesh.registry` module ships in 1.5.5 (unreleased).

1. The module uses only the standard library.
2. It checks each download against the SHA-256 hash that the service reports.
3. The client talks to https://valence-api.iqh766c2.workers.dev by default. Set `VALENCE_API_URL` only to override it.
4. The service allows 500 downloads per token per day.

```python
from chilmesh import registry

names = registry.list_meshes()                       # filenames, e.g. "Block_O.14"
mesh = registry.fetch_mesh("Block_O.14")            # verified download, returns a CHILmesh
registry.download_mesh("Block_O.14", "mesh.14")      # or save the verified bytes to a file
```

Errors derive from `registry.RegistryError`:

- `RegistryAuthError`: the token is missing, invalid, or expired.
- `RegistryNotFoundError` (also a `ValueError`): the mesh name is unknown or invalid.
- `RegistryUnavailableError`: no API URL is set, the network failed, or the API is rate limited or down.
- `RegistryIntegrityError`: the downloaded bytes do not match the reported hash, or a response is malformed.

---

## Examples and CLI (moved from the README)

Run the example scripts:

```bash
python ../examples/01_quickstart.py        # load, stats, plot
python ../examples/02_fort14_roundtrip.py  # fort.14 read and write
python ../examples/03_smoothing.py         # angle-based smoother
python ../examples/04_spatial_queries.py   # find_element, radius search, k-nearest
```

Use the command-line tool:

```bash
chilmesh info mesh.fort.14                                   # vertex, element and layer counts
chilmesh summary mesh.fort.14                                # header-only metadata, no full parse
chilmesh convert mesh.2dm mesh.fort.14                       # format conversion
chilmesh smooth mesh.fort.14 -o out.fort.14 --method fem     # smooth and write
chilmesh plot mesh.fort.14 -o mesh.png --quality             # render
```

`python -m chilmesh` is equivalent, and every subcommand takes `--help`.

Four fixtures ship with the package: `chilmesh.examples.{annulus, donut, block_o, structured}()`. `chilmesh.CHILmesh` remains an alias of `Mesh`.

---

## Roadmap items moved from the README

- Next: a faster CPU renderer and a full-lifecycle benchmark on the global STOFS-2D-Global mesh ([#276](https://github.com/domattioli/CHILmesh/issues/276)).
- Research: a medial-axis `skeletonize()` ([#223](https://github.com/domattioli/CHILmesh/issues/223)); size-controlled mesh cartograms ([#219](https://github.com/domattioli/CHILmesh/issues/219)).
