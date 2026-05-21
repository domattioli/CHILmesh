# Architecture

**Analysis Date:** 2026-05-21

## System Overview

CHILmesh is a single-class 2D mesh library built around `CHILmesh` — a monolithic mesh container that owns all topology, geometry, and layer (skeletonization) state. The design prioritizes O(1) adjacency queries and fast initialization via hash-mapped edge lookup.

```text
┌────────────────────────────────────────────────────────────────┐
│                     Public API Layer                            │
│  (read_from_fort14, read_from_2dm, from_admesh_domain, etc.)   │
└──────────┬────────────────────────────────────────┬────────────┘
           │                                        │
           ▼                                        ▼
┌───────────────────────────┐        ┌─────────────────────────────┐
│   CHILmesh Core           │        │   I/O Adapters               │
│   (CHILmesh.py)           │        │   - fort.14 parser           │
│                           │        │   - .2dm parser              │
│  - Mesh container         │        │   - ADMESH-Domains bridge    │
│  - Adjacencies dict       │        │   - bridge.py adapters       │
│  - Layers (skeletonize)   │        │                              │
│  - Smoothing methods      │        │   Located: CHILmesh.py       │
│  - Quality/analysis       │        │            bridge.py         │
│  - Visualization (mixin)  │        │            admesh_warmstart  │
│                           │        │                              │
│  Located: src/chilmesh/   │        │   src/chilmesh/bridge.py     │
│           CHILmesh.py     │        │   src/chilmesh/admesh_*.py   │
└──────────┬────────────────┘        └──────────┬──────────────────┘
           │                                   │
           ▼                                   ▼
┌───────────────────────────────────────────────────────────────┐
│                    Data Layer                                  │
│                                                               │
│  Core Mesh Representation:                                   │
│  - points: ndarray[n_verts, 3]       (x, y, z coordinates)   │
│  - connectivity_list: ndarray[n_elems, 3|4]  (element verts) │
│                                                               │
│  Adjacency Structures (in adjacencies dict):                 │
│  - Edge2Vert: ndarray[n_edges, 2]    (edge endpoints)        │
│  - Elem2Edge: ndarray[n_elems, 3|4]  (element edges)         │
│  - Edge2Elem: ndarray[n_edges, 2]    (adjacent elements)     │
│  - Vert2Edge: Dict[int, Set[int]]    (vertex incident edges) │
│  - Vert2Elem: Dict[int, Set[int]]    (vertex incident elems) │
│  - EdgeMap: hash table                (O(1) edge lookup)      │
│                                                               │
│  Layer Output (after skeletonize):                           │
│  - layers['OE']: List[ndarray] (outer elements per layer)    │
│  - layers['IE']: List[ndarray] (inner elements per layer)    │
│  - layers['OV']: List[ndarray] (outer vertices per layer)    │
│  - layers['IV']: List[ndarray] (inner vertices per layer)    │
│  - layers['bEdgeIDs']: List[ndarray] (boundary edges)        │
│                                                               │
│  Located: CHILmesh.py attributes                            │
└───────────────────────────────────────────────────────────────┘
           │
           ▼
┌───────────────────────────────────────────────────────────────┐
│              Supporting Infrastructure                         │
│                                                               │
│  mesh_topology.py: EdgeMap (hash-based edge ID lookup)       │
│  layer_paths.py: Traversal helpers for layer graph extraction│
│  mutations.py: Element add/remove primitives (MutableMesh)   │
│  plot_utils.py: CHILmeshPlotMixin for visualization         │
│  cli.py: Command-line interface (argparse-based)            │
│  examples.py: Built-in fixtures (annulus, donut, block_o,   │
│               structured)                                     │
│                                                               │
│  Located: src/chilmesh/                                      │
└───────────────────────────────────────────────────────────────┘
```

## Component Responsibilities

| Component | Responsibility | File |
|-----------|----------------|------|
| **CHILmesh** | Main container; owns all state (points, connectivity, adjacencies, layers); exports all public methods | `src/chilmesh/CHILmesh.py` (2,292 lines) |
| **EdgeMap** | Hash table for O(1) edge ID lookup by vertex pair (v1, v2) → edge_id; used only during adjacency building | `src/chilmesh/mesh_topology.py` |
| **Bridge adapters** | Duck-typed integration points for MADMESHR, ADMESH, ADMESH-Domains; layer on top of public API | `src/chilmesh/bridge.py` |
| **I/O parsers** | fort.14 and .2dm format readers/writers; both in CHILmesh.py as static methods | `src/chilmesh/CHILmesh.py` (read_from_fort14, read_from_2dm, write_fort14) |
| **Smoothing** | FEM, angle-based, and ADMESH truss smoothers; mutate points in-place via change_points() | `src/chilmesh/CHILmesh.py` (__direct_smoother, __angle_based_smoother, admesh_warmstart.py) |
| **Skeletonization** | Medial axis extraction via concentric boundary peeling (layer-based); stores output in layers dict | `src/chilmesh/CHILmesh.py` (_skeletonize) |
| **Visualization** | Plotting (wireframe, quality, layers); added as mixin from utils | `src/chilmesh/utils/plot_utils.py` |
| **Mutations** | Experimental: add/remove elements; tracked in MutableMesh; limited, not all operations supported | `src/chilmesh/mutations.py` |
| **CLI** | Command-line interface (info, convert, smooth, plot); wrapper over public API | `src/chilmesh/cli.py` |
| **Examples** | Test fixtures (annulus, donut, block_o, structured); used in tests and CLI | `src/chilmesh/examples.py` |

## Pattern Overview

**Overall:** Single-class monolithic container with hash-optimized adjacency structures

**Key Characteristics:**
- **One mesh class** — `CHILmesh` owns all state; no separate topology/geometry classes
- **Eager adjacency building** — All adjacencies built once during initialization (unless compute_layers=False)
- **Hash-mapped edge lookup** — O(1) amortized via EdgeMap; eliminates O(n²) bottleneck
- **Vectorized operations** — All geometry (signed_area, angles, quality) via numpy ufuncs
- **In-place mutations** — Smoothing, element addition, coordinate changes all mutate state
- **Optional skeletonization** — compute_layers=False skips expensive medial axis extraction for fast bulk loading

## Layers

**Initialization Layer:**
- Purpose: Parse input files, validate mesh, build adjacencies
- Location: `CHILmesh.__init__()`, `_initialize_mesh()`, `read_from_fort14()`, `read_from_2dm()`, `from_admesh_domain()`
- Contains: I/O parsing, CCW orientation enforcement, adjacency building, optional skeletonization
- Depends on: numpy arrays, EdgeMap hash table
- Used by: All users (entry point for mesh loading)

**Topology Layer:**
- Purpose: Manage mesh connectivity; provide O(1) adjacency queries
- Location: `_build_adjacencies()`, `_identify_edges()`, `_build_elem2edge()`, `_build_vert2edge()`, `_build_vert2elem()`, `_build_edge2elem()`
- Contains: EdgeMap, Edge2Vert, Elem2Edge, Vert2Edge, Vert2Elem, Edge2Elem structures
- Depends on: numpy arrays, dict/set operations
- Used by: Skeletonization, smoothing, mutations, spatial queries

**Geometry Layer:**
- Purpose: Compute mesh properties (signed area, interior angles, element quality)
- Location: `signed_area()`, `interior_angles()`, `elem_quality()`
- Contains: Vectorized numpy operations; fully exploits SIMD via numpy ufuncs
- Depends on: numpy arrays, connectivity_list, points
- Used by: Quality analysis, visualization, smoothing convergence checks

**Skeletonization Layer:**
- Purpose: Decompose mesh into concentric layers via boundary peeling
- Location: `_skeletonize()`
- Contains: Iterative boundary extraction; tracks OV, OE, IE, IV per layer; stores layer output in layers dict
- Depends on: Edge2Vert, Edge2Elem, connectivity_list
- Used by: Layer-based queries, visualization, downstream analysis

**Spatial Index Layer:**
- Purpose: Enable fast point-in-element and k-nearest-vertex queries
- Location: `_build_spatial_indices()` (lazy), `find_element()`, `nearest_vertices()`, `find_elements_in_radius()`
- Contains: cKDTree on element centroids; barycentric coordinate checking
- Depends on: scipy.spatial.cKDTree
- Used by: Spatial queries (optional, v0.4.0+)

**Smoothing Layer:**
- Purpose: Improve mesh quality via node movement
- Location: `smooth_mesh()`, `_direct_smoother()` (FEM), `_angle_based_smoother()` (angle-based), `admesh_warmstart.py` (spring-based)
- Contains: Sparse matrix assembly (scipy.sparse), stiffness solvers, iterative refinement
- Depends on: scipy.sparse.lil_matrix, scipy.sparse.linalg.spsolve, points array
- Used by: Mesh quality improvement workflows

**Visualization Layer:**
- Purpose: Render mesh (wireframe, quality, layers)
- Location: `src/chilmesh/utils/plot_utils.py` (CHILmeshPlotMixin)
- Contains: matplotlib patches, colormap selection, axis setup
- Depends on: matplotlib, points, connectivity_list, optional: layers, quality arrays
- Used by: CLI plotting, Jupyter notebooks

## Data Flow

### Primary Request Path: Load Mesh

1. **Entry** — `CHILmesh.read_from_fort14('mesh.14')` (`CHILmesh.py:1126`)
   - Calls `_fort14_to_arrays()` (lines ~1160–1180) to parse file
   - Returns (points, connectivity, grid_name)

2. **Initialization** — `CHILmesh.__init__(connectivity, points, ...)` (`CHILmesh.py:153`)
   - Stores points, connectivity_list, grid_name
   - Calls `_initialize_mesh(compute_layers=True)`

3. **Mesh Setup** — `_initialize_mesh()` (`CHILmesh.py:205`)
   - Ensures z-coordinates on points (pad if 2D)
   - Calls `_ensure_ccw_orientation()` to fix element winding
   - Calls `_build_adjacencies()` if compute_layers=True
   - Calls `_build_spatial_indices()` (always, for future spatial queries)

4. **Adjacency Building** — `_build_adjacencies()` (`CHILmesh.py:391`)
   - Calls `_identify_edges()` → returns EdgeMap and edge list (`src/chilmesh/mesh_topology.py`)
   - Calls `_build_elem2edge()` using EdgeMap for O(1) lookups
   - Calls `_build_vert2edge()` from Edge2Vert
   - Calls `_build_vert2elem()` from connectivity
   - Calls `_build_edge2elem()` to pair edges with adjacent elements
   - Calls `_validate_adjacencies()` to check invariants
   - Stores all structures in `adjacencies` dict

5. **Skeletonization** — `_skeletonize()` (`CHILmesh.py:958`)
   - Iteratively peels boundary layers inward
   - Tracks OV (outer vertices), OE (outer elements), IE (inner elements), IV (inner vertices), bEdgeIDs per layer
   - Uses working copies of Edge2Vert and Edge2Elem marked with -1 as consumed
   - Stores layers in `layers` dict with keys: OE, IE, OV, IV, bEdgeIDs

6. **Output** — Mesh ready for queries, visualization, analysis
   - `.n_verts`, `.n_elems`, `.n_edges`, `.n_layers` populated
   - `.points`, `.connectivity_list` available
   - `.adjacencies` dict (Vert2Edge, Edge2Vert, Elem2Edge, Edge2Elem, Vert2Elem, EdgeMap)
   - `.layers` dict (OE, IE, OV, IV, bEdgeIDs)

### Secondary Flow: Element Quality Analysis

1. **Entry** — `mesh.elem_quality()` (`CHILmesh.py:~750`)
2. **Signed Area** — `self.signed_area()` (fully vectorized; separate paths for tri vs quad)
3. **Interior Angles** — `self.interior_angles()` (vectorized; law of cosines on per-element geometry)
4. **Quality Metrics** — Min/max/skew quality per element
5. **Output** — Tuple of (quality_array, angles_array, stats_dict)

**State Management:**
- None — quality analysis is read-only; does not mutate mesh state
- Uses: points, connectivity_list, Edge2Vert, elem type
- Produces: temporary numpy arrays (no saved state)

### Tertiary Flow: Mesh Smoothing

1. **Entry** — `mesh.smooth_mesh(method='fem', acknowledge_change=True)` (`CHILmesh.py:~800`)
2. **Method Selection** — 'fem' → `_direct_smoother()`, 'angle-based' → `_angle_based_smoother()`
3. **FEM Smoother** — `_direct_smoother()` (`CHILmesh.py:~850`)
   - Builds sparse stiffness matrix (scipy.sparse.lil_matrix) from nodal adjacencies
   - Solves Kx = 0 for interior nodes (boundary nodes fixed)
   - Moves points via `change_points(new_points, acknowledge_change=True)`
4. **Angle-Based Smoother** — `_angle_based_smoother()` (iterative)
   - For each interior node, maximize minimum incident angle
   - Centroid-based repositioning
   - Repeats up to n_iter times or until convergence
5. **Mutation** — `change_points()` updates `self.points` in-place
6. **Post-Smooth** — Adjacencies unchanged (topology preserved); may need to update spatial index if used

## Key Abstractions

**EdgeMap:**
- Purpose: O(1) edge ID lookup by vertex pair
- Example: `edge_map.find_edge(v1, v2)` → edge_id (or None)
- Pattern: Hash table with canonical form (min, max) to avoid duplicates
- Location: `src/chilmesh/mesh_topology.py` (lines 9–125)
- Used during: `_build_adjacencies()` only (ephemeral, not stored in public API)

**Adjacency Dictionary:**
- Purpose: Centralized access to topology
- Structure: `adjacencies = {"Edge2Vert": ..., "Elem2Edge": ..., ...}`
- Access: Public getters like `mesh.edge2vert(edge_id)`, `mesh.elem2edge(elem_id)`
- Invariants: All vertex/element/edge IDs valid; no duplicates; edges in canonical form

**Layers (Skeletonization Output):**
- Purpose: Represent mesh structure as concentric rings
- Structure: `layers = {"OE": [layer0_outer_elems, ...], "IE": [...], ...}`
- Invariants: Layers form disjoint cover of all elements; layer sizes decrease monotonically
- Access: `mesh.layers[key][layer_idx]` or `mesh.elements_in_layer(layer_idx)`

## Entry Points

**File I/O:**
- `CHILmesh.read_from_fort14(path, compute_layers=True)` — Load ADCIRC mesh (`src/chilmesh/CHILmesh.py:1126`)
- `CHILmesh.read_from_2dm(path, compute_layers=True)` — Load SMS mesh (`src/chilmesh/CHILmesh.py:1244`)
- `CHILmesh.from_admesh_domain(record, compute_layers=True)` — Load from catalog (`src/chilmesh/CHILmesh.py:1195`)

**Examples:**
- `chilmesh.examples.annulus()` — Small triangular mesh
- `chilmesh.examples.donut()` — Medium donut-shaped mesh
- `chilmesh.examples.block_o()` — Large O-shaped mesh (100k elements)
- `chilmesh.examples.structured()` — Structured quad mesh

**CLI:**
- `chilmesh info <mesh>` — Print mesh statistics
- `chilmesh convert <in> <out>` — Format conversion
- `chilmesh smooth <mesh> -o <out> --method <m>` — Smooth mesh
- `chilmesh plot <mesh> -o <image>` — Generate visualization

## Architectural Constraints

- **Threading:** Single-threaded event loop (Python GIL); no explicit threading; scipy calls may release GIL for C operations
- **Global state:** None; all state owned by CHILmesh instance; EdgeMap is ephemeral (discarded after adjacency build)
- **Circular imports:** None detected; import order: `__init__` → CHILmesh → utils.plot_utils → examples
- **Memory model:** Dense arrays (points, connectivity_list, Edge2Vert, Elem2Edge) + sparse dicts (Vert2Edge, Vert2Elem); no memory pooling or custom allocators

## Anti-Patterns

### Iteration Over Vert2Edge/Vert2Elem

**What happens:** Code naively copies full set on every access:
```python
edges = mesh.adjacencies['Vert2Edge'][v]  # Returns a set
```

**Why it's wrong:** Defensive copy overhead if called in tight loop; potential for external code to mutate returned set and break invariants.

**Do this instead:** 
```python
edges = mesh.get_vertex_edges(v)  # Returns Set[int]; public getter with documented contract
```
Location: `src/chilmesh/CHILmesh.py` (getter methods added Phase 2)

### Manual Edge Canonicalization

**What happens:** 
```python
edge = (v1, v2)  # Not sorted
if edge in edge2vert:  # Fails if edge stored as (v2, v1)
```

**Why it's wrong:** Edge ordering matters; canonical form (min, max) enforced in EdgeMap and Edge2Vert construction, but naive code breaks.

**Do this instead:**
```python
edge_id = edge_map.find_edge(v1, v2)  # EdgeMap handles canonicalization
```
Location: `src/chilmesh/mesh_topology.py` (EdgeMap.find_edge, lines 41–53)

### Assuming Fixed 3-Column Connectivity

**What happens:**
```python
for elem_id in range(n_elems):
    v0, v1, v2 = connectivity_list[elem_id]  # Fails on quads or mixed
```

**Why it's wrong:** Mixed meshes have 4-column connectivity with padded triangles; manual unpacking breaks.

**Do this instead:**
```python
tri_elems, quad_elems = mesh._elem_type()  # Classify elements first
# Or:
connectivity = mesh.connectivity_list[elem_id]  # Array indexing; works for 3 or 4 columns
```
Location: `src/chilmesh/CHILmesh.py` (lines 347–389, _elem_type method)

## Error Handling

**Strategy:** Fail fast with descriptive messages; validate at boundaries (I/O, adjacency checking)

**Patterns:**

1. **I/O Validation** (fort.14 parsing):
   - Check file existence, format version, element count
   - Raise FileNotFoundError, ValueError on malformed input

2. **Adjacency Validation** (post-build):
   - `_validate_adjacencies()` checks: all vertices in Vert2Edge/Vert2Elem; edge IDs valid; no missing entries
   - Raises AssertionError with diagnostic message

3. **Orientation Checks** (CCW enforcement):
   - Compute signed area; flip if negative
   - Handles mixed-element case (padded triangles vs real quads separately)

## Cross-Cutting Concerns

**Logging:** None (pure library; no logging framework)

**Validation:** 
- Adjacency invariants checked post-build via `_validate_adjacencies()`
- File format validation in read_from_fort14/2dm (line count, element count, vertex references)

**Authentication:** Not applicable (file-based, no network I/O)

---

*Architecture analysis: 2026-05-21*
