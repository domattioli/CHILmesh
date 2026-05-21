# Technology Stack

**Analysis Date:** 2026-05-21

## Languages

**Primary:**
- Python 3.8+ - Core mesh library, all public APIs, CLI
  - **Location:** `src/chilmesh/` (2,292 lines in main class)
  - **Notes:** Pure Python with numpy/scipy core, actively developed

**Secondary:**
- MATLAB - Original legacy implementation (unmaintained reference)
  - **Location:** `src/@CHILmesh/CHILmesh.m` 
  - **Status:** No longer maintained; Python implementation is production target
  - **Deprecated:** MATLAB version documented in CITATION.cff; MSc thesis reference only

## Runtime

**Environment:**
- CPython 3.8–3.12 (tested via pyproject.toml classifiers)
- Platform: OS-independent (Linux, macOS, Windows)

**Package Manager:**
- pip (primary distribution channel via PyPI)
- conda-forge (planned support per README)
- setuptools (build backend) — configured in `pyproject.toml` with package-dir mapping to `src/`
- Lockfile: None (requirements.txt is advisory, not locked)

## Frameworks

**Core:**
- **numpy** ≥1.23 - Dense mesh arrays, vectorized geometry ops
  - Edge2Vert, Elem2Edge, element connectivity lists (n_elems × 3|4 arrays)
  - Vectorized signed_area(), interior_angles(), elem_quality() implementations
  - All batch operations on 100k+ element meshes
  
- **scipy** ≥1.10 - Sparse matrix solvers, spatial indexing
  - `scipy.sparse` (lil_matrix) — Finite element method (FEM) smoother stiffness assembly
  - `scipy.spatial.Delaunay` — Mesh I/O fallback triangulation, test fixtures
  - `scipy.spatial.cKDTree` — O(log n) spatial queries (point-in-element via centroid tree, nearest neighbors)
  - `scipy.sparse.linalg.spsolve` — Direct FEM smoother (Balendran method)

- **matplotlib** ≥3.6 - Visualization (plotting API not part of public mesh operations)
  - `matplotlib.cm` — Quality colormaps
  - `matplotlib.pyplot` — Per-element quality plots, layer wireframe rendering
  - Located in `src/chilmesh/utils/plot_utils.py` (CHILmeshPlotMixin)

**Testing:**
- **pytest** ≥7 - Test runner (439 passing, 26 test files, ~17.2s total runtime)
  - Fixtures: annulus, donut, block_o (100k elements), structured quad mesh
  - Parametrization across all 4 fixtures for geometry/topology invariants
  - Configuration: `pyproject.toml` [pytest.ini_options] with `slow` marker for block_o exclusion
  - Optional: pytest-cov (coverage reporting), pytest-xdist (parallel runs)

**Build/Dev:**
- **setuptools** ≥61 - Build backend; package discovery via `src/` layout
- **build** - PEP 517 wheel/sdist builder (dev dependency)
- **twine** - PyPI upload tool (dev dependency)

## Key Dependencies

**Critical (mesh core):**
- **numpy** ≥1.23
  - Why: Adjacency structures (Edge2Vert ndarray[n_edges, 2], Elem2Edge ndarray[n_elems, 3|4]) require vectorized ops for 100k+ element performance
  - Usage: Dense element arrays, batch geometric transforms, layer masks (np.isin, np.concatenate, np.unique)

- **scipy** ≥1.10
  - Why: Spatial indexing (cKDTree) enables O(log n) point-in-element queries; sparse solvers required for FEM smoothing stiffness matrix
  - Usage: `cKDTree` for centroid/vertex queries, `lil_matrix` + `spsolve` for direct FEM solver

- **matplotlib** ≥3.6
  - Why: Visualization API (quality plots, layer overlays) part of public API but optional (imports gated)
  - Usage: colormap selection, patch rendering in CHILmeshPlotMixin

**Infrastructure (optional/integration):**
- **scipy** (Delaunay) — Test fixture generation, example meshes
- **admesh-domains** (duck-typed, not a hard dependency) — Catalog integration via `from_admesh_domain()` class method
  - No import required; uses attribute introspection (`getattr(record, 'filename', None)`)
  - Failure mode: FileNotFoundError with helpful message if file path not found

**Vendored:**
- **_vendor_admesh_truss.py** (298 lines) — Spring-based smoother from ADMESH project
  - Standalone; no external ADMESH import required
  - Integrated via `optimize_with_admesh_truss()` public API

## Configuration

**Environment:**
- No required .env file
- Optional: External mesh files (ADMESH-Domains registry) via `from_admesh_domain()`
  - Search paths: explicit path → /tmp/admesh-domains/ → home directory
  - See `benchmark_wnat_hagen.py` for example

**Build:**
- `pyproject.toml` [build-system]: requires setuptools ≥61, backend "setuptools.build_meta"
- `pyproject.toml` [project.scripts]: CLI entry point `chilmesh = "chilmesh.cli:main"`
- Package layout: `src/` directory (mapping in [tool.setuptools] package-dir)
- Data files: `src/chilmesh/data/*.fort.14` bundled fixtures (annulus, donut, block_o, structured)

## Platform Requirements

**Development:**
- Python ≥3.8, ≤3.12
- pip / conda / uv (package installer)
- Virtual environment recommended (`python -m venv`)

**Production:**
- Python ≥3.8
- numpy, scipy, matplotlib (installed via pip/conda)
- No compilation required (pure Python + numpy C extensions)
- Tested on: Linux (CI), macOS, Windows (self-reported via setup.py classifiers)

**Performance Baseline (v0.4.1):**
| Operation | Mesh | Time | Optimization |
|-----------|------|------|--------------|
| Fast init (no layers) | WNAT_Hagen (52.7k verts) | 0.44s | Adjacency building deferred |
| Full init (with layers) | WNAT_Hagen | 3.26s | Hash-mapped EdgeMap O(1) edge lookup |
| Quality analysis | WNAT_Hagen | 0.07s | Vectorized numpy |
| Total workflow | WNAT_Hagen | 3.33s | 937× faster than v0.1.1 (was ~13.4s) |

**Scaling Profile:**
- Edge lookup: O(1) amortized via EdgeMap hash table (`find_edge()` <0.2 μs/call, 5k samples)
- Vert2Edge adjacency: O(1) dict access (~0.17 μs/call)
- Elem2Edge bulk: O(1) per-element ndarray access (~0.5 μs/call)
- Spatial queries: O(log n) via cKDTree (`find_element()` ~50 μs/call)
- Skeletonization: O(n) graph traversal with numpy set operations

## Dependencies Graph

```
CHILmesh (public API)
├── numpy ≥1.23 (mesh arrays, vectorized ops)
├── scipy ≥1.10
│   ├── scipy.spatial.cKDTree (spatial indices)
│   ├── scipy.sparse (FEM stiffness)
│   └── scipy.spatial.Delaunay (test fixtures)
├── matplotlib ≥3.6 (visualization, optional at import)
├── admesh_warmstart.py (ADMESH Spring-Based Truss, vendored)
│   └── _vendor_admesh_truss.py (298 lines, no external deps)
└── mesh_topology.py
    └── EdgeMap (hash-based edge O(1) lookup, local class)

CLI (chilmesh.cli:main)
└── argparse (stdlib)

Test Suite
└── pytest ≥7, pytest-cov, pytest-xdist (optional)
```

## External Integrations at Runtime

**File I/O:**
- `.fort.14` (ADCIRC fort.14 format) — Text-based, parsed line-by-line
- `.2dm` (SMS format) — Text-based, node + element sections
- Both: 1-based vertex indexing → 0-based Python arrays (offset handling in readers)

**Downstream Projects:**
- **ADMESH** — uses `CHILmesh` as topology backend; accesses via public API
- **MADMESHR** — quad/mixed-element generation; depends on CHILmesh data structures
- **ADMESH-Domains** — mesh catalog; `from_admesh_domain(record)` integrates without hard import

**Version & Metadata:**
- `__version__` pulled from importlib.metadata (looks up installed package)
- Fallback: "0.0.0" if not installed (e.g., local development before setup)

---

*Stack analysis: 2026-05-21*
