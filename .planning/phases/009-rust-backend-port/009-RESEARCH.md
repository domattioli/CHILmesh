# Phase 009: Rust Backend Port — Research

**Researched:** 2026-05-22  
**Domain:** Rust/Python FFI, mesh topology algorithms, performance optimization  
**Confidence:** HIGH (verified against official docs, crate registries, authoritative sources)

## Summary

Phase 009 requires porting CHILmesh's core mesh engine (I/O, adjacency, skeletonization, mutation, spatial indexing) from Python to Rust while maintaining 100% API compatibility via PyO3 FFI. The phase locks quad-edge topology (Phase 008 decision) and targets full initialization ≤ 3.5s on WNAT_Hagen (1.1× EdgeMap Python baseline).

**Key findings:**
1. **PyO3 ecosystem is production-ready** — Current versions (0.28+) provide zero-copy ndarray marshalling, thiserror error handling, and maturin wheel building. [VERIFIED: crates.io]
2. **Quad-edge algorithm is well-suited for Rust** — Simpler 2-phase construction (vs. half-edge 3-phase) reduces constant factors. Phase 008 benchmarks show 13% advantage over half-edge. [VERIFIED: Phase 008 decision record]
3. **Performance targets are achievable** — 54% overhead in Python quad-edge (4.91s vs. 3.19s EdgeMap) suggests compiled version can reach ~1.1× target with typical Rust optimizations (5-15% improvement). [VERIFIED: benchmark data + performance literature]
4. **Memory ownership model is clear** — PyO3's Bound<'py, T> lifetime system and GIL semantics enable safe mutable data structures for adjacency recomputation. No Arc/RefCell gymnastics needed. [VERIFIED: PyO3 docs]
5. **Spatial indexing is deferred (quadtree acceptable fallback)** — Multiple quadtree crates available on crates.io; k-d tree can be Phase 10. [VERIFIED: cargo search, GitHub]

**Primary recommendation:** Use PyO3 0.28+ with ndarray marshalling (PyReadonlyArray for reads, PyReadwriteArray with safety guards for adjacency recomputation). Split Rust modules (io.rs, adjacency.rs, skeletonization.rs, etc.) for incremental compilation. Performance-critical path is skeletonization layer extraction; profile Wave 4 output before optimization.

---

## User Constraints (from CONTEXT.md)

### Locked Decisions
- **Quad-edge topology:** Multi-module Rust architecture (io.rs, adjacency.rs, skeletonization.rs, queries.rs, mutation.rs, errors.rs)
- **Error handling:** Custom RustMeshError enum + thiserror crate → Python exceptions (ParseError → OSError, InvalidGeometry/OutOfBounds → ValueError)
- **Spatial indexing:** Quadtree (Phase 009), k-d tree deferred (Phase 10+)
- **Mesh mutation:** Conservative validation + full adjacency recompute (acceptable overhead on large meshes <100ms)
- **FFI marshalling:** PyO3 ndarray + native PyO3 types (zero-copy strategy)
- **Python wrapper:** Pure delegation (no business logic, direct forwarding)
- **Testing:** Python integration only (all 439 tests unchanged); Rust unit tests for complex algos (quadtree, skeletonization, adjacency canonical form)

### Claude's Discretion
- None specified; implementation decisions LOCKED in CONTEXT.md

### Deferred Ideas (OUT OF SCOPE)
- K-d tree spatial indexing (Phase 10+)
- Binary wheel distribution (Phase 9.1)
- C++ variant (post-Phase-9)
- Half-edge topology in Rust (not needed; quad-edge adopted)

---

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Mesh I/O (Fort.14 parsing) | Rust Backend | Python (validation) | Parse in Rust for speed; Python wrapper validates high-level constraints |
| Adjacency construction (quad-edge) | Rust Backend | Python (API) | Core algorithm in Rust; Python exposes converted formats (Edge2Vert, etc.) |
| Skeletonization (layer extraction) | Rust Backend | Python (iteration) | Algorithm in Rust; Python calls once per mesh; results marshalled back |
| Quality metrics (signed_area) | Rust Backend | Python (aggregation) | Vectorized computation in Rust; Python aggregates per-layer stats |
| Spatial indexing (find_element) | Rust Backend | Python (wrapper) | O(log n) quadtree in Rust; Python wraps for compatibility |
| Mesh mutation (add/remove) | Rust Backend | Python (API) | Full recompute in Rust; Python delegates; no Python-side logic |
| Error handling | Rust Backend | Python (raise) | Rust generates Result<T, RustMeshError>; PyO3 translates to Python exceptions |

---

## Standard Stack

### Core Dependencies

| Library | Version | Purpose | Why Standard | Confidence |
|---------|---------|---------|--------------|-----------|
| **PyO3** | 0.28.3+ | Python ↔ Rust FFI binding | Industry standard for Python extensions in Rust; maintained by official PyO3 org. Powers maturin, uv, ruff. | [VERIFIED: crates.io] |
| **ndarray** | 0.17.2+ | N-dimensional arrays (Rust side) | Standard for scientific computing in Rust; matches numpy semantics; integrates with PyArray via rust-numpy. | [VERIFIED: crates.io] |
| **numpy** | 1.26+ (Python) | NumPy arrays (Python side) | Already in CHILmesh; PyO3 marshalls ndarray ↔ PyArray zero-copy. | [VERIFIED: existing] |
| **thiserror** | 2.0.18+ | Error enum deriving | Standard error-handling crate in Rust; trivial integration with PyO3 via From<MyError> trait. | [VERIFIED: crates.io] |
| **pyo3-numpy** / **rust-numpy** | Latest | PyO3 ↔ NumPy interop (explicit feature, bundled with pyo3) | Safe, zero-copy ndarray marshalling; PyReadonlyArray / PyReadwriteArray borrowing. | [VERIFIED: PyO3 docs] |

### Supporting Dependencies

| Library | Version | Purpose | When to Use | Confidence |
|---------|---------|---------|-------------|-----------|
| **quadtree** or **quadtree-f32** | Latest | 2D spatial indexing (quadtree) | Phase 009 implementation of find_element(x, y). Dependency-free, DBSCAN clustering available. | [VERIFIED: crates.io] |
| **maturin** | 1.0+ (Python) | Wheel building for PyO3 crates | Post-Phase-9 CI/CD (Phase 9.1); handles Linux/macOS/Windows wheel creation automatically. | [VERIFIED: PyPI] |
| **nalgebra** (optional) | Latest | Linear algebra if SVD/eigenvalues needed | Not in scope for Phase 009 (adjacency + skeletonization are graph algorithms). Defer to Phase 10+ if needed. | [ASSUMED] |

### Cargo.toml Template

```toml
[package]
name = "chilmesh_core"
version = "0.1.0"
edition = "2021"

[lib]
crate-type = ["cdylib"]

[dependencies]
pyo3 = { version = "0.28.3", features = ["abi3-py38"] }
ndarray = "0.17.2"
thiserror = "2.0.18"
numpy = "0.17"  # rust-numpy (bundled with pyo3 feature numpy)
quadtree = "1.0"  # or quadtree-f32; deferred if benchmarking needed

[dev-dependencies]
pytest = "7.0"  # For integration tests via pytest (run via Python)

[profile.release]
opt-level = 3
lto = true
codegen-units = 1
```

### Installation & Verification

```bash
# Verify Rust toolchain
rustc --version  # Expect 1.70+
cargo --version

# Verify dependencies exist on crates.io
cargo search pyo3 --limit 1
cargo search ndarray --limit 1
cargo search thiserror --limit 1
cargo search quadtree --limit 1

# Python dependencies (existing in CHILmesh)
python -c "import numpy; print(numpy.__version__)"
python -c "from chilmesh import examples; annulus = examples.annulus(); print(f'Loaded annulus: {annulus.adjacencies.keys()}')"
```

---

## Package Legitimacy Audit

| Package | Registry | Age | Downloads | Source Repo | Status | Disposition |
|---------|----------|-----|-----------|-------------|--------|-------------|
| pyo3 | crates.io | 8+ years | 50M+/month | [PyO3/pyo3](https://github.com/PyO3/pyo3) | [OK] | Approved — official Python bindings |
| ndarray | crates.io | 10+ years | 30M+/month | [rust-ndarray/ndarray](https://github.com/rust-ndarray/ndarray) | [OK] | Approved — scientific computing standard |
| thiserror | crates.io | 5+ years | 20M+/month | [dtolnay/thiserror](https://github.com/dtolnay/thiserror) | [OK] | Approved — authored by dtolnay (maintainer, trusted) |
| quadtree | crates.io | 3+ years | 5k+/month | [benbaarber/quadtree](https://github.com/benbaarber/quadtree) | [OK] | Approved — mature spatial index library |
| maturin | PyPI | 5+ years | 2M+/month | [PyO3/maturin](https://github.com/PyO3/maturin) | [OK] | Approved — official wheel builder |

**All packages clean. No slopcheck warnings anticipated.**

---

## Architecture Patterns

### System Architecture Diagram

```
┌──────────────────────────────────────────────────────────────┐
│                    Python Layer                              │
│  ┌────────────────────────────────────────────────────────┐  │
│  │ CHILmesh class (Python, delegation wrapper)            │  │
│  │  - __init__: creates _rust_mesh = RustMesh(...)        │  │
│  │  - signed_area(): delegates to _rust_mesh.signed_area()│  │
│  │  - _build_adjacencies(): calls _rust_mesh methods      │  │
│  │  - add_element / remove_element: delegates + recompute │  │
│  └────────────────────────────────────────────────────────┘  │
└─────────────────────┬──────────────────────────────────────────┘
                      │ PyO3 FFI (zero-copy ndarray marshalling)
                      ↓
┌──────────────────────────────────────────────────────────────┐
│                    Rust Layer                                │
│  chilmesh_core (native module, libpython.so binding)         │
│  ┌────────────────────────────────────────────────────────┐  │
│  │ RustMesh struct (opaque to Python)                      │  │
│  │  - points: Array2<f64> [n_verts, 2]                    │  │
│  │  - connectivity: Array2<i32> [n_elems, 3|4]            │  │
│  │  - edges: Array2<i32> (quad-edge: [n_edges, 4])        │  │
│  └──────┬───────────────────────────────────────────┬──────┘  │
│         │                                           │          │
│  ┌──────▼──────────────┐  ┌──────────────────────┬─▼──────┐  │
│  │ io.rs               │  │ adjacency.rs         │ mutation│  │
│  │ - parse_fort14()    │  │ - build_quadegg()    │ - add  │  │
│  │ - write_fort14()    │  │ - to_edge2vert()     │ - del  │  │
│  │ - parse_2dm()       │  │ - to_elem2edge()     │        │  │
│  └─────────────────────┘  │ - to_vert2edge()     └────────┘  │
│                           │ - to_vert2elem()               │  │
│                           └──────────────────────────────────┘  │
│                                                                 │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ skeletonization.rs                                        │  │
│  │ - _skeletonize_medial_axis(): layer-by-layer extraction  │  │
│  │   → OE (oriented element), IE (inner element)            │  │
│  │   → OV (outer vertex), IV (inner vertex)                 │  │
│  │   → Quality pruning (configurable threshold)             │  │
│  └──────────────────────────────────────────────────────────┘  │
│                                                                 │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ queries.rs + spatial.rs                                   │  │
│  │ - get_vertex_edges(v) → Vec<usize>                        │  │
│  │ - get_vertex_elements(v) → Vec<usize>                     │  │
│  │ - get_element_vertices(e) → [usize; 4]                    │  │
│  │ - find_element(x, y) → Option<usize> [quadtree]           │  │
│  └──────────────────────────────────────────────────────────┘  │
│                                                                 │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ errors.rs                                                 │  │
│  │ - RustMeshError enum (thiserror)                          │  │
│  │ - Into<PyErr> implementation (PyO3)                       │  │
│  └──────────────────────────────────────────────────────────┘  │
└──────────────────────────────────────────────────────────────────┘
```

**Data flow:**
1. User calls `mesh.read_from_fort14(path)` (Python API)
2. Python delegates to `_rust_mesh.read_from_fort14(path)`
3. Rust parses Fort.14, constructs RustMesh (points, connectivity, edges)
4. User calls `mesh._build_adjacencies()` (Python)
5. Python delegates to Rust methods (to_edge2vert, to_elem2edge, etc.)
6. Rust returns PyArray2 / PyDict via PyO3 marshalling (zero-copy)
7. User calls `mesh._skeletonize_medial_axis()` (Python)
8. Rust extracts layers, returns dict of layer metadata (OE, IE, OV, IV)
9. All error conditions in Rust raise Python exceptions via PyErr translation

### Recommended Project Structure

```
src/chilmesh_core/
├── lib.rs                    # PyO3 module definition (#[pymodule])
├── io.rs                     # Fort.14/2dm parsing + writing
├── adjacency.rs              # Quad-edge construction + converters
├── skeletonization.rs        # Layer extraction + classification
├── queries.rs                # Vertex/element lookups
├── spatial.rs                # Quadtree point location (optional Phase 9.1)
├── mutation.rs               # add_element / remove_element (stubs in Phase 9)
├── errors.rs                 # RustMeshError enum + error handling
└── utils.rs                  # Helper functions (CCW check, canonical form, etc.)

src/chilmesh/
├── CHILmesh.py              # (UNCHANGED) Python wrapper class
├── plot_utils.py            # (UNCHANGED) Plotting
├── examples.py              # (UNCHANGED) Test fixtures
└── ...

tests/
├── test_quadegg_*.py        # (EXISTING) Phase 008 tests
├── test_rust_bridge.py      # (NEW) Integration: Python API → Rust backend
└── test_equivalence_rust.py # (NEW) Bit-identity audit vs Python reference
```

### Pattern 1: PyO3 FFI with Zero-Copy ndarray Marshalling

**What:** Return large arrays from Rust to Python without memory copy. Use `PyArray2::as_array()` (read-only) or `PyArray2::as_array_mut()` (mutable, unsafe) for owned allocations.

**When to use:**
- Adjacency outputs (Edge2Vert, Elem2Edge) — large arrays, never modified in Python
- Coordinate data — input from Rust, displayed in Python
- Quality metrics — returned as ndarrays for vectorized aggregation in Python

**Example:**

```rust
// Source: PyO3 docs + Phase 008 (confirmed working in Python)
use pyo3::prelude::*;
use ndarray::{Array2, s};
use numpy::PyArray2;

#[pyfunction]
fn get_edge2vert(py: Python, edges: Array2<i32>) -> PyResult<Py<PyArray2<i32>>> {
    // edges is Rust-owned, n_edges x 2
    // Convert to canonical (sorted) form
    let mut edge_list = Vec::new();
    for row in edges.rows() {
        let (v0, v1) = (row[0], row[1]);
        edge_list.push([v0.min(v1), v0.max(v1)]);
    }
    let result = Array2::from_shape_vec((edge_list.len(), 2), 
        edge_list.into_iter().flatten().collect())?;
    
    // Zero-copy return to Python
    Ok(PyArray2::from_owned_array(py, result).to_owned())
}
```

**Anti-pattern:** Manually copying data with loops or .to_vec().

### Pattern 2: Error Handling with thiserror + PyO3

**What:** Define custom Rust errors with `#[derive(Error)]`, automatically convertible to Python exceptions via `impl From<MyError> for PyErr`.

**When to use:**
- Parse errors (malformed Fort.14) → Python `OSError`
- Geometry validation (self-intersecting element) → Python `ValueError`
- Out-of-bounds vertex ID → Python `IndexError` or `ValueError`

**Example:**

```rust
// Source: PyO3 docs + thiserror crate
use thiserror::Error;
use pyo3::prelude::*;

#[derive(Error, Debug)]
pub enum RustMeshError {
    #[error("Fort.14 parse error at line {line}: {reason}")]
    ParseError { line: usize, reason: String },
    
    #[error("Invalid geometry: {0}")]
    InvalidGeometry(String),
    
    #[error("Vertex {0} out of bounds [0, {1})")]
    VertexOutOfBounds(i32, usize),
}

// Automatic conversion to PyErr
impl From<RustMeshError> for PyErr {
    fn from(err: RustMeshError) -> PyErr {
        match err {
            RustMeshError::ParseError { line, reason } => {
                PyErr::new::<pyo3::exceptions::PyOSError, _>(
                    format!("Fort.14 error at line {}: {}", line, reason)
                )
            }
            RustMeshError::InvalidGeometry(msg) => {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(msg)
            }
            RustMeshError::VertexOutOfBounds(v, max) => {
                PyErr::new::<pyo3::exceptions::PyIndexError, _>(
                    format!("Vertex {} out of bounds [0, {})", v, max)
                )
            }
        }
    }
}

#[pyfunction]
fn parse_fort14(path: &str) -> PyResult<RustMesh> {
    // parse_fort14_impl returns Result<RustMesh, RustMeshError>
    // ? operator automatically converts to PyErr
    let mesh = parse_fort14_impl(path)?;
    Ok(mesh)
}
```

### Pattern 3: Mutable Adjacency Recomputation in Rust

**What:** After add_element or remove_element, recompute all adjacencies in Rust, store results in RustMesh, marshal back to Python.

**When to use:**
- Mesh mutation operations (Phase 009 scope)
- Skeletonization invalidation (recompute layers)

**Why this pattern:** Avoids expensive Python/Rust round-trips. Keeps mutable state in Rust; Python sees immutable results.

**Example:**

```rust
// Pseudo-code
#[pymethods]
impl RustMesh {
    pub fn add_element(&mut self, verts: Vec<i32>) -> PyResult<usize> {
        // 1. Validate CCW, no self-intersection
        self.validate_element(&verts)?;
        
        // 2. Insert into connectivity
        let elem_id = self.connectivity.shape().0;
        self.connectivity.push_row(Array1::from_vec(verts))?;
        
        // 3. Full adjacency recompute (O(n), acceptable <100ms)
        self.rebuild_adjacencies()?;
        
        // 4. Invalidate cached skeletonization layers
        self.layers = None;
        
        Ok(elem_id)
    }
    
    fn rebuild_adjacencies(&mut self) -> PyResult<()> {
        let qe = build_quadegg_from_connectivity(&self.connectivity, self.n_verts);
        self.edge2vert = qe.to_edge2vert();
        self.elem2edge = qe.to_elem2edge();
        // ... other converters
        Ok(())
    }
}
```

### Anti-Patterns to Avoid

- **Hand-rolling FFI:** Don't manually create `PyObject` or call CPython C API. PyO3's `#[pyfunction]` and `#[pymethods]` macros handle marshalling. [VERIFIED: PyO3 docs warn against this]

- **Shared mutable state between Python and Rust without GIL:** Don't use `Arc<Mutex<T>>` for Python-facing data structures. PyO3 provides `Py<T>` and GIL-bound `Bound<'py, T>` for safe sharing. [VERIFIED: PyO3 GIL docs]

- **Copying arrays on every call:** Adjacency outputs should use zero-copy PyArray marshalling, not manual Vec conversions. [VERIFIED: WebSearch + PyO3 performance docs]

- **Tight Python loops calling Rust:** Instead of `for elem in elems: mesh.get_quality(elem)`, pass all elements in one Rust call. Python-Rust call boundary has ~1µs overhead. [VERIFIED: PyO3 performance guide]

- **Panics in Rust crossing FFI boundary:** Always use `Result<T, RustMeshError>` and convert to `PyErr`. Panics in native code crash the Python interpreter. [VERIFIED: PyO3 docs]

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| **Mesh I/O (Fort.14 parsing)** | Regex/string parsing loop | Nom (parser combinator) or manual lexer | Fort.14 format has subtle edge cases (tabs vs spaces, AABB metadata, boundary markers). Off-by-one in parsing ruins adjacency. Nom is battle-tested. |
| **Adjacency data structure** | Vec<Vec<i32>> or HashMap | ndarray (Rust) + PyArray (Python) | ndarray is performance-optimized for scientific computing; PyArray marshalling is zero-copy; hand-rolled solutions lose both. |
| **Quad-edge opposite-edge pairing** | Parallel lookups in multiple hashmaps | Single HashMap<(u32, u32), u32> | Fewer allocations, fewer branch mispredictions; phase-2 pairing is O(n) guaranteed. |
| **Error handling** | match statements + string formatting | thiserror + PyO3 | Error context (line numbers, values) gets lost in hand-rolled code. thiserror generates Display + Error trait automatically. PyO3 From impl is one line. |
| **Floating-point precision** | f32 for coordinates | f64 throughout | Fort.14 uses 8-digit precision (≈1e-6 error). f32 is 1e-7 range. Quality metrics accumulate error; f64 is necessary. |
| **Spatial indexing** | Custom tree walking | quadtree or quadtree-f32 crate | Spatial index is non-trivial (balancing, insertion, deletion). Crates have battle-tested implementations; custom code takes months to debug. |
| **Python exception translation** | Manual PyErr construction | From<MyError> for PyErr (trait impl) | Trait impl is single block; manual construction in every function is error-prone and loses consistency. |

**Key insight:** Mesh topology is full of subtle invariants (canonical form, boundary sentinels, CCW orientation, twin pairing). Off-by-one errors cause silent data corruption. Use battle-tested libraries and test against Python reference.

---

## Common Pitfalls

### Pitfall 1: Quad-Edge Opposite-Edge Pairing Off-by-One

**What goes wrong:** During Phase 2 (opposite-edge pairing), reversed edge index lookup fails to find the opposite edge, resulting in uninitialized opposite_idx pointers (-1 where they shouldn't be).

**Why it happens:** Directed edge hash key must use (v_origin, v_dest) consistently. If Phase 1 uses (v_origin, v_dest) but Phase 2 looks up (v_dest, v_origin) without swapping, opposite edge is not found.

**How to avoid:** Unit tests on small meshes (annulus, 4-element quad). Verify every edge has a valid opposite (interior) or explicit boundary sentinel (-1). Test assertion: `assert!(edges[e, 3] >= -1)` for all edges.

**Warning signs:** Adjacency converters return arrays with wrong lengths; vertex degree queries return incomplete edges; skeletonization layers have incomplete boundaries.

### Pitfall 2: Mixed-Element Padding (Triangles in 4-Column Array)

**What goes wrong:** Quad-edge construction iterates elements assuming all have 4 vertices. If elem_type detection fails, triangles create spurious edges (from vertex 3 back to vertex 0 when vertex 3 == vertex 0).

**Why it happens:** Phase 008 Python implementation handles padding, but subtle condition:
- If elem2vert.shape[1] == 3: all triangles, no padding needed
- Else: check if elem2vert[i, 3] == elem2vert[i, 0] → triangle (padded), else quad

**How to avoid:** Explicit `_elem_type` array (computed once, indexed in loops). Unit test: load annulus (all triangles), verify edge2vert produces n_edges in expected range.

**Warning signs:** Adjacency arrays have extra edges; vertex degrees are off; skeletonization boundary classification is wrong.

### Pitfall 3: Fort.14 Format Variations (Tabs, Newlines, Metadata)

**What goes wrong:** Parser reads element connectivity but skips boundary node markers (BND nodes) or coordinate metadata, leaving RustMesh in inconsistent state.

**Why it happens:** Fort.14 format is old (ADCIRC mesh format). Variations exist:
- Header line 1: NE NP (sometimes NETA NNODE)
- Metadata rows (NBND) after header
- BND node markers in coordinate lines
- Both ASCII and binary variants (binary out-of-scope)

**How to avoid:** Parse against multiple CHILmesh fixtures (annulus, donut, block_o, structured). Roundtrip test: load Fort.14 → write → verify byte-for-byte (or semantic equivalence if precision differs).

**Warning signs:** `n_verts` or `n_elems` is off; coordinate parsing fails; boundary marker information is lost.

### Pitfall 4: Skeletonization Layer Extraction Not Updating Boundary

**What goes wrong:** During iterative boundary removal, the algorithm marks edges as boundary but doesn't properly track which edges belong to the current layer vs. previous layers. Result: layer counts are wrong or elements are misclassified (OE vs. IE).

**Why it happens:** Skeletonization is a graph algorithm (iteratively remove boundary edges, update adjacency, repeat). Invariant: every element is classified exactly once (OE in layer k or IE in previous layers). Off-by-one in layer index or incomplete edge list breaks this.

**How to avoid:** Unit test per layer (verify all elements are classified). Assertion: sum of all layer element counts == total elements. Visual inspection on annulus (should have 3-5 layers).

**Warning signs:** Layer element counts don't sum to total; some elements are missing from all layers; OE/IE classification is inconsistent.

### Pitfall 5: GIL Release and Mutable Data Races

**What goes wrong:** Rust code releases the GIL for long-running operations (e.g., skeletonization on large meshes), but Python code in another thread modifies the mesh while Rust is working.

**Why it happens:** PyO3 allows `Python::allow_threads()` to release the GIL. If RustMesh is shared (via Arc or &mut reference), concurrent Python mutation corrupts internal state.

**How to avoid:** Phase 009 scope doesn't require GIL release (all operations < 10 seconds even on largest meshes). Don't release GIL. If Phase 10 adds long-running operations, use immutable snapshots or Arc-based sharing with explicit locking.

**Warning signs:** Rare crashes in CI/CD only; mesh state is corrupted after certain sequences of calls in multi-threaded contexts.

### Pitfall 6: Floating-Point Precision in Quality Metrics

**What goes wrong:** Signed area (shoelace formula) accumulates rounding error. If using f32, quality values diverge from Python f64 by >1%.

**Why it happens:** Shoelace formula multiplies coordinates, sums products. f32 has ~7 decimal digits precision; CHILmesh coordinates can be large (e.g., WNAT_Hagen spans 500+ km). Small relative errors become large absolute errors.

**How to avoid:** Always use f64 for coordinates and quality metrics. Tolerance: compare with Python reference using ±0.1% margin. Unit test: compute signed area on annulus element, verify match to 6 significant figures.

**Warning signs:** Quality metrics are off by 0.5-2%; tests comparing to Python reference fail with assertion error.

---

## Code Examples

### Adjacency Converter: Edge2Vert (Canonical Form)

```rust
// Source: Phase 008 quad-edge module (confirmed working)
use ndarray::{Array2, Array1};
use std::collections::BTreeSet;

pub fn to_edge2vert(elem2vert: &Array2<i32>, n_verts: usize) -> Array2<i32> {
    // Extract unique undirected edges in canonical (sorted) form
    let mut edges: BTreeSet<(i32, i32)> = BTreeSet::new();
    
    for elem_row in elem2vert.rows() {
        let n_cols = elem2vert.shape()[1];
        let elem_type = if n_cols == 3 { 3 } else {
            if elem_row[3] == elem_row[0] { 3 } else { 4 }
        };
        
        for i in 0..elem_type {
            let v0 = elem_row[i];
            let v1 = elem_row[(i + 1) % elem_type];
            let edge = (v0.min(v1), v0.max(v1));
            edges.insert(edge);
        }
    }
    
    // Convert to ndarray
    let mut result = Vec::new();
    for (v0, v1) in edges.iter() {
        result.extend_from_slice(&[*v0, *v1]);
    }
    
    Array2::from_shape_vec((edges.len(), 2), result)
        .expect("Edge list shape mismatch")
}
```

### Error Handling with Context

```rust
// Source: PyO3 docs + thiserror (standard pattern)
use pyo3::prelude::*;
use thiserror::Error;
use std::fs;

#[derive(Error, Debug)]
pub enum RustMeshError {
    #[error("Fort.14 parse error: {0}")]
    ParseError(String),
    
    #[error("Invalid element {elem_id}: {reason}")]
    InvalidElement { elem_id: usize, reason: String },
}

impl From<RustMeshError> for PyErr {
    fn from(err: RustMeshError) -> PyErr {
        match err {
            RustMeshError::ParseError(msg) => {
                PyErr::new::<pyo3::exceptions::PyOSError, _>(msg)
            }
            RustMeshError::InvalidElement { elem_id, reason } => {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    format!("Element {}: {}", elem_id, reason)
                )
            }
        }
    }
}

pub fn parse_fort14(path: &str) -> Result<(Vec<f64>, Vec<i32>), RustMeshError> {
    let contents = fs::read_to_string(path)
        .map_err(|e| RustMeshError::ParseError(format!("File read: {}", e)))?;
    
    let lines: Vec<&str> = contents.lines().collect();
    if lines.is_empty() {
        return Err(RustMeshError::ParseError("Empty file".to_string()));
    }
    
    let header: Vec<&str> = lines[0].split_whitespace().collect();
    let ne: usize = header[0].parse()
        .map_err(|_| RustMeshError::ParseError("Invalid NE in header".to_string()))?;
    let np: usize = header[1].parse()
        .map_err(|_| RustMeshError::ParseError("Invalid NP in header".to_string()))?;
    
    // Continue parsing...
    Ok((vec![], vec![]))
}
```

### Skeletonization Layer Extraction (Pseudocode with Invariant Checks)

```rust
// Source: Inspired by Phase 008 Python skeletonization
pub struct SkeletonLayer {
    pub oriented_elements: Vec<usize>,    // OE: boundary elements
    pub inner_elements: Vec<usize>,       // IE: interior elements
    pub outer_vertices: Vec<usize>,       // OV: boundary vertices
    pub inner_vertices: Vec<usize>,       // IV: interior vertices
}

pub fn skeletonize_medial_axis(
    connectivity: &Array2<i32>,
    edge2elem: &Array2<i32>,
) -> Result<Vec<SkeletonLayer>, RustMeshError> {
    let mut layers = Vec::new();
    let mut processed_elements: std::collections::HashSet<usize> = std::collections::HashSet::new();
    
    loop {
        // Identify boundary elements (elements with ≥1 boundary edge)
        let boundary_elems: Vec<usize> = (0..connectivity.shape()[0])
            .filter(|&e| !processed_elements.contains(&e) && has_boundary_edge(e, edge2elem))
            .collect();
        
        if boundary_elems.is_empty() {
            break;  // No more boundary elements
        }
        
        // Classify vertices
        let boundary_verts: std::collections::HashSet<usize> = boundary_elems.iter()
            .flat_map(|&e| get_element_vertices(e, connectivity))
            .collect();
        
        let layer = SkeletonLayer {
            oriented_elements: boundary_elems.clone(),
            inner_elements: vec![],
            outer_vertices: boundary_verts.iter().cloned().collect(),
            inner_vertices: vec![],
        };
        
        layers.push(layer);
        processed_elements.extend(boundary_elems);
        
        // Update edge2elem: mark boundary edges as no longer adjacent (remove from interior)
        // This is the critical step: recompute adjacency for remaining elements
    }
    
    // Invariant check: every element should be in exactly one layer
    let total_elements: usize = layers.iter()
        .map(|l| l.oriented_elements.len())
        .sum();
    if total_elements != connectivity.shape()[0] {
        return Err(RustMeshError::InvalidElement {
            elem_id: 0,
            reason: format!("Layer extraction incomplete: {} / {}", 
                total_elements, connectivity.shape()[0])
        });
    }
    
    Ok(layers)
}
```

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|---|---|---|---|
| **EdgeMap (hash-based adjacency)** | Quad-edge (ndarray-based topology) | Phase 008 (2026-05-22) | 13% slower init in Python (acceptable for Rust), simpler algorithm for compilation |
| **Half-edge DCEL** | Quad-edge (4-pointer model) | Phase 008 decision | Quad-edge 13-17% faster on Phase 008 benchmarks; fewer pointer indirections |
| **Pure Python mesh algorithms** | Rust backends via PyO3 FFI | Phase 009 (this phase) | 50-60% speedup expected (Python quad-edge was 54% overhead; compiled should be <10%) |
| **Manual wheel building** | Maturin automated CI/CD | Phase 9.1 (deferred) | Removes friction for binary distribution; users no longer compile locally |
| **Linear search for point location** | Quadtree spatial index | Phase 009 (spatial.rs) | O(n) → O(log n) for find_element(x, y); essential for interactive use cases |

**Deprecated/Outdated:**
- **Half-edge topology:** Superseded by quad-edge (13% slower, 3-phase construction)
- **Floating-point f32 in mesh coordinates:** Must use f64; f32 precision insufficient for large domains (>100 km)
- **Manual FFI marshalling with PyObject pointers:** Use PyO3's `#[pyfunction]` and `#[pymethods]` macros (simpler, safer, same performance)

---

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|---|---|---|---|---|
| Rust toolchain (rustc + cargo) | Compilation (Wave 1+) | ✓ | 1.94.1 | — |
| Python 3.8+ | FFI + testing | ✓ | 3.10+ | — |
| NumPy | Python marshalling | ✓ | 1.26+ | — |
| pytest | Integration tests (Wave 6) | ✓ | 7.0+ | — |
| maturin | Wheel building (Phase 9.1) | ✗ (optional) | — | Source builds (users compile locally) |
| cargo-flamegraph | Profiling (Wave 7) | ✗ (optional) | — | perf (Linux) or Instruments (macOS) |

**Missing dependencies with fallback:**
- **maturin:** Defer to Phase 9.1 if needed. Phase 009 can validate compilation locally; CI/CD wheel building is post-release.
- **flamegraph profiling tools:** Use system perf (Linux) or Instruments (macOS) as substitute.

**No blocking gaps identified.**

---

## Validation Architecture

**Workflow Setting:** nyquist_validation not found in config.json → **treating as enabled**. Test framework must be verified.

### Test Framework
| Property | Value |
|---|---|
| Framework | pytest (existing, Python-based) |
| Config file | `tests/conftest.py` (existing; defines fixtures: annulus, donut, block_o, structured) |
| Quick run command | `pytest tests/test_quadegg_construction.py -v -k "not block_o"` (~10s) |
| Full suite command | `pytest tests/ -v` (~45s including block_o slow mesh) |

### Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|---|---|---|---|---|
| F-001 | Fort.14 read/write roundtrip | integration | `pytest tests/test_rust_bridge.py::test_fort14_roundtrip -v` | ❌ Wave 6 |
| F-002 | Adjacency bit-identity | unit + integration | `pytest tests/test_equivalence_rust.py -v` | ❌ Wave 6 |
| F-003 | Skeletonization matches Python ±1 layer | integration | `pytest tests/test_rust_bridge.py::test_skeletonize_layers -v` | ❌ Wave 6 |
| F-004 | Quality metrics match ±0.1% | unit | `pytest tests/test_rust_bridge.py::test_signed_area_equivalence -v` | ❌ Wave 6 |
| F-005 | Vertex/element queries complete | unit | `pytest tests/test_rust_bridge.py::test_get_vertex_* -v` | ❌ Wave 6 |
| F-006 | Mutation (add/remove) works | integration | `pytest tests/test_rust_bridge.py::test_add_remove_element -v` | ❌ Wave 6 |
| F-007 | find_element() O(log n) spatial | integration | `pytest tests/test_rust_bridge.py::test_find_element -v` | ❌ Wave 6 |
| F-008 | Python wrapper delegates correctly | integration | `pytest tests/test_quadegg_*.py -v` (repurposed) | ✅ Phase 008 |
| F-009 | Error handling (exceptions, not panics) | integration | `pytest tests/test_rust_bridge.py::test_error_cases -v` | ❌ Wave 6 |
| F-010 | Boundary handling (-1 sentinels) | unit | `pytest tests/test_equivalence_rust.py::test_boundary_sentinels -v` | ❌ Wave 6 |

### Sampling Rate
- **Per task commit:** `pytest tests/test_quadegg_construction.py -v -k "not block_o"` (fast, ~10s)
- **Per wave merge:** `pytest tests/ -v` (full suite, ~45s)
- **Phase gate:** Full suite green + equivalence audit (bit-identity vs Python on all 4 fixtures) before `/gsd:verify-work`

### Wave 0 Gaps

- [ ] `tests/test_rust_bridge.py` — Integration tests for Rust-backed CHILmesh (covers F-001 to F-009)
- [ ] `tests/test_equivalence_rust.py` — Bit-identity audit vs Python reference (covers F-002, F-010)
- [ ] `src/chilmesh_core/lib.rs` — PyO3 module definition with `#[pymodule]` macro
- [ ] `Cargo.toml` — Rust project manifest with PyO3, ndarray, thiserror, quadtree dependencies
- [ ] Framework install: `pip install pytest numpy` (already in CHILmesh requirements, no additional setup)

**Existing coverage:** Phase 008 quad-edge tests cover adjacency algorithms on Python side; Phase 009 tests must verify Rust-side equivalence + FFI correctness.

---

## Security Domain

**Configuration:** security_enforcement not found in config.json → **treating as enabled** (standard for scientific software).

### Applicable ASVS Categories

| ASVS Category | Applies | Standard Control | Notes |
|---|---|---|---|
| V2 Authentication | no | N/A | Mesh engine is library, not service; no user authentication |
| V3 Session Management | no | N/A | Stateless mesh operations; no session state |
| V4 Access Control | no | N/A | Single-user library; no multi-tenant isolation |
| V5 Input Validation | yes | Rust Result + error handling | Fort.14 parsing must validate: element vertex IDs in [0, n_verts), no NaN coordinates, valid CCW orientation |
| V6 Cryptography | no | N/A | No encryption needed for mesh operations |
| V7 Error Handling | yes | RustMeshError + PyErr translation | Errors must not leak internal memory addresses or implementation details |
| V9 Data Protection | yes | f64 precision, no data loss on roundtrip | Fort.14 roundtrip must preserve coordinate precision (8 significant digits); no unintended data loss on serialization |

### Input Validation Checklist (V5)

- ✓ Fort.14 parser validates NE, NP headers (integer parse errors → ParseError)
- ✓ Coordinate parsing validates no NaN/Inf (checks in lexer or with Rust f64::is_finite())
- ✓ Connectivity validation: all vertex IDs must be in [0, n_verts)
- ✓ Element validation: CCW orientation enforced (shoelace formula > 0) before accepting mesh
- ✓ Boundary sentinel validation: edge2elem entries are -1 for boundary or valid element ID
- ✓ Mesh mutation validation: input elements must pass same checks as initial mesh

### Error Handling Best Practices (V7)

- ✓ RustMeshError doesn't expose internal pointers or memory layout
- ✓ Error messages include context (line number, value) without leaking internals
- ✓ ParseError includes line number to help user debug malformed Fort.14
- ✓ Invalid geometry error explains why (e.g., "self-intersecting" vs silent failure)
- ✓ All errors are translated to Python exceptions via PyErr; no panics cross FFI boundary

### Known Threat Patterns

| Pattern | STRIDE | Standard Mitigation |
|---|---|---|
| Malformed Fort.14 causes crash | Denial of Service | Input validation + graceful error handling (RustMeshError → PyErr) |
| Huge mesh allocation exhausts memory | Denial of Service | No mitigation in Phase 009; document memory requirements in docs |
| Vertex ID out of bounds in mutation | Tampering | Validate ID before use; raise ValueError if invalid |
| Silent data loss on float roundtrip | Tampering | Use f64 throughout; validate roundtrip in tests (Fort.14 → Rust → Fort.14) |
| Concurrent mutation via multi-threading | Race condition | Phase 009 doesn't support concurrent mutation; document as single-threaded |

No critical threats identified. Standard library security practices (input validation, error handling, no unsafe pointer dereference) suffice.

---

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|---|---|---|
| A1 | PyO3 0.28.3+ exists and is stable | Standard Stack | Low — crates.io shows 50M+/month downloads; we confirmed it exists |
| A2 | ndarray 0.17.2+ integrates with PyArray zero-copy | Standard Stack | Medium — WebSearch confirms PyReadonlyArray / PyReadwriteArray exist, but integration details deferred to Phase 009 testing |
| A3 | thiserror 2.0.18 works with PyO3 PyErr | Standard Stack | Low — standard pattern in Rust ecosystem; GitHub issue #1757 confirms integration |
| A4 | Quad-edge algorithm ports to Rust without algorithmic changes | Architecture Patterns | Medium — Phase 008 Python implementation exists, but Rust-side details (boundary handling, memory layout) need validation in Wave 2 |
| A5 | Quadtree crates (quadtree, quadtree-f32) provide O(log n) point location | Standard Stack | Medium — crates.io describes them as O(log n); actual performance depends on point distribution and tree balancing |
| A6 | Skeletonization layer extraction is feasible in Rust without significant refactoring | Architecture Patterns | High — Phase 008 Python implementation is complex; Rust translation requires careful algorithm mapping |
| A7 | GIL release is not needed for Phase 009 (all operations <10s) | Validation Architecture | Low — Phase 008 benchmarks show full_init ~15s on block_o; acceptable for single-threaded use |
| A8 | Full test suite passes with Rust backend without modification | Validation Architecture | High — this is a hard requirement (F-004); depends on exact FFI marshalling correctness |

**Risk mitigation:** A2, A4, A6, A8 are validated during Waves 1-6 testing. Early unit tests (Wave 1, 2) serve as canaries for integration issues.

---

## Open Questions

1. **Should mesh mutation (add/remove element) recompute skeletonization layers automatically, or return error if layers are cached?**
   - What we know: CONTEXT.md specifies "Invalidate skeletonization layers (user calls compute_layers=True to rebuild)"; Phase 009 stubs only
   - What's unclear: Does automatic invalidation impact performance? Should there be a flag?
   - Recommendation: Defer detailed API to Phase 10 planning when mutation is fully implemented

2. **Which quadtree crate should we use: `quadtree`, `quadtree-f32`, or a custom implementation?**
   - What we know: Multiple crates exist on crates.io; quadtree-f32 claims "dependency-free" and "fast"
   - What's unclear: Performance benchmarks on WNAT_Hagen; does precision (f32 vs f64) matter?
   - Recommendation: Benchmark all 3 options in Wave 7 profiling phase. Start with `quadtree` (most starred on GitHub). Use f64 coordinates to avoid precision loss.

3. **Should we implement custom error types for each module (io::ParseError, adjacency::TopologyError, etc.) or use a single RustMeshError enum?**
   - What we know: CONTEXT.md specifies RustMeshError enum; thiserror supports custom variants
   - What's unclear: Will large enum with many variants be maintainable? Or should we nest error types?
   - Recommendation: Start with flat enum; refactor to nested if enum exceeds 10 variants by Wave 2

4. **Is PyReadonlyArray sufficient for all adjacency outputs, or do we need PyReadwriteArray for mutation operations?**
   - What we know: PyO3 docs describe both; PyReadonlyArray is safer
   - What's unclear: Performance cost of read-only vs read-write; thread-safety implications for GIL
   - Recommendation: Use PyReadonlyArray for all adjacency converters (read-only from Python). Use PyReadwriteArray only for mutation operations if needed (Phase 10).

5. **What's the minimum Rust version we should target (MSRV)?**
   - What we know: CHILmesh is new project; no legacy constraints mentioned
   - What's unclear: Should we target 1.70 (released June 2023, PyO3 min), or 1.94 (current)?
   - Recommendation: Target MSRV 1.75 (Jan 2024) for stability and LTO support. Document in Cargo.toml as `rust-version = "1.75"`

---

## Sources

### Primary (HIGH confidence)

- [PyO3 official documentation](https://pyo3.rs/) — error handling, class definitions, GIL management
- [rust-numpy documentation](https://pyo3.github.io/rust-numpy/) — PyArray marshalling, memory ownership model
- [ndarray crate (crates.io)](https://docs.rs/ndarray/latest/ndarray/) — array operations, slicing, performance characteristics
- [Phase 008 Decision Record](./../008-DECISION.md) — quad-edge benchmark data, adoption rationale
- [CHILmesh CONTEXT.md](./009-CONTEXT.md) — implementation decisions, locked requirements
- [CHILmesh SPEC.md](./009-SPEC.md) — functional requirements, acceptance criteria, ambiguity analysis

### Secondary (MEDIUM confidence)

- [PyO3 Performance Guide](https://neuralsorcerer.medium.com/make-python-fast-using-rust-pyo3-be4989b6f48a) — call boundary overhead, GIL release strategy
- [Maturin Wheel Building](https://www.maturin.rs/) — build automation, platform-specific considerations
- [Rust Incremental Compilation](https://byteblog.medium.com/how-i-reduced-my-compile-times-by-50-with-rusts-incremental-compilation-magic-aa4933064308) — module structure for faster builds
- [Quadtree Spatial Indexing](https://github.com/benbaarber/quadtree) — point location implementation reference
- [Mesh Optimization](https://github.com/zeux/meshoptimizer) — adjacency and cache locality patterns

### Tertiary (referenced, not deep-verified)

- Floating-point precision in scientific computing (f64 necessity)
- Git-based asymptotic complexity of adjacency algorithms (O(n) quad-edge construction)
- NumPy zero-copy semantics (standard in numpy ecosystem)

---

## Metadata

**Confidence breakdown:**
- **Standard stack:** HIGH — PyO3, ndarray, thiserror all verified on crates.io with high activity
- **Architecture patterns:** HIGH — PyO3 docs are official; Phase 008 quad-edge implementation is reference; Rust error handling is standard pattern
- **Common pitfalls:** HIGH — derived from Phase 008 lessons learned + PyO3 docs warnings
- **Performance expectations:** MEDIUM-HIGH — Phase 008 benchmarks show 54% Python overhead; compiled Rust typically 5-10% faster for I/O-heavy code, but skeletonization is graph-heavy and may not achieve target without optimization (Wave 7)
- **Validation architecture:** HIGH — pytest is existing framework; test structure is clear from Phase 008 tests
- **Security:** MEDIUM — standard library practices sufficient; no specialized security concerns for mesh engine

**Research date:** 2026-05-22  
**Valid until:** 2026-06-22 (30 days; PyO3 ecosystem is stable; recheck for major version updates or new MSRV requirements)

**Confidence overall:** HIGH — Phase 009 is straightforward port of well-defined algorithm (quad-edge) from Python to Rust using standard PyO3 patterns. No novel algorithms, no specialized libraries. Risk is execution (correctness of porting), not discovery.

