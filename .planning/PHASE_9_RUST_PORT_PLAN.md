# Phase 9 Plan: Full Rust Port with Python API Wrapper

**Phase:** 009-rust-backend-port  
**Goal:** Rewrite CHILmesh core (mesh I/O, adjacency, skeletonization, queries) in Rust; expose via PyO3 FFI maintaining 100% API compatibility.  
**Target:** Full init ≤ 3.5s WNAT_Hagen (1.1× EdgeMap baseline). All 439 tests PASS via Python wrapper.

---

## Scope Definition

### IN SCOPE (Rust implementation required)

**Core Mesh I/O:**
- Fort.14 reader (ASCII parsing, binary format support)
- 2dm reader
- Fort.14 writer (ASCII output)
- Binary mesh serialization (optional, Phase 9.2)

**Adjacency Structures:**
- Quad-edge backend (primary, from Phase 008)
- Optional: EdgeMap clone (fallback if Rust build slow)
- Optional: Half-edge variant (research-only, not prioritized)

**Skeletonization (medial axis extraction):**
- Layer-by-layer boundary extraction (_skeletonize_medial_axis)
- Oriented element (OE) + inner vertex (IV) detection
- Outer vertex (OV) identification
- Quality metrics (signed area per layer)
- Pruning algorithm (if configured)

**Quality Analysis:**
- signed_area() — element quality metric
- median_quality() — per-layer statistics
- _ensure_ccw_orientation() — preprocessing

**Mesh Queries:**
- get_vertex_edges(vert_id) → Set[edge_ids]
- get_vertex_elements(vert_id) → Set[elem_ids]
- get_element_vertices(elem_id) → ndarray[3|4]
- find_element(x, y) → elem_id (point location, optional Phase 9.2)

**ADMESH Integration:**
- admesh_warmstart.py bridge (stays in Python; calls Rust mesh ops)
- No changes needed; Rust mesh object drops into existing Python code

### OUT OF SCOPE (stays Python)

**Plotting/Visualization:**
- plot_utils.py remains Python (matplotlib, no Rust equivalent needed)
- Visualization calls Rust APIs, renders results in Python

**MADMESHR Bridge:**
- Stays Python; wraps Rust mesh object
- High-level simulation logic untouched

**Test Framework:**
- pytest remains Python
- Tests call Python API (which wraps Rust)
- No Rust unit tests (trust Python integration tests)

**Mesh Mutation (add/remove elements):**
- Phase 2 designed but not implemented; defer to Phase 10
- Rust port focuses on read-heavy skeletonization workflow

---

## Architecture: Python API + Rust Backend

### High-Level Design

```
┌─────────────────────────────────────────┐
│   Python Layer (src/chilmesh/*.py)      │
│  - CHILmesh class (public API)          │
│  - I/O: read_from_fort14(), write()     │
│  - Queries: get_vertex_edges(), etc.    │
│  - ADMESH/MADMESHR bridges              │
└──────────────────┬──────────────────────┘
                   │ PyO3 FFI
                   ↓
┌─────────────────────────────────────────┐
│   Rust Layer (src/chilmesh_core/)       │
│  - RustMesh struct (opaque to Python)   │
│  - I/O: parse_fort14(), write_fort14()  │
│  - Adjacency: QuadEdge, EdgeMap         │
│  - Skeletonization engine               │
│  - Quality analysis                     │
└─────────────────────────────────────────┘
```

### Python Wrapper Strategy

**Approach:** Minimal translation layer. Python CHILmesh object holds `RustMesh` FFI handle + metadata (filename, CRS). Method calls delegate to Rust; return types (ndarrays, dicts, sets) marshalled back to Python.

**Example:**
```python
class CHILmesh:
    def __init__(self, ...):
        self._rust_mesh = chilmesh_core.RustMesh(...)  # FFI call
        
    def _build_adjacencies(self, ...):
        # Rust does the work; Python just asks
        edge2vert = self._rust_mesh.get_edge2vert()  # Returns ndarray
        elem2edge = self._rust_mesh.get_elem2edge()  # Returns ndarray
        ...
        
    def signed_area(self):
        return self._rust_mesh.signed_area()  # Scalar from Rust
```

**Wrapper Implementation Effort:** ~20 hours (straightforward delegation; no logic).

---

## Task Breakdown (8 Waves, ~24 weeks)

### Wave 1: Rust Project Setup & I/O Foundation (Week 1–2)

**Task 1.1: Initialize Rust project + PyO3 scaffolding**
- Create `src/chilmesh_core/` Rust crate
- Add `Cargo.toml` with ndarray, numpy dependencies
- PyO3 module skeleton: `lib.rs` with `#[pymodule]`
- Dummy test: `python -c "from chilmesh_core import RustMesh; print('OK')"`
- **Effort:** 4 hours | **Risk:** low

**Task 1.2: Fort.14 parser (ASCII reading)**
- Implement `parse_fort14(path: &str) -> RustMesh`
- Parse header (NE, NP, nbnd)
- Element-to-vertex connectivity (tri + quad support)
- Coordinate (x, y) reading
- Boundary node marker handling
- **Test:** Load annulus, donut, Block_O; verify n_verts, n_elems, n_edges
- **Effort:** 12 hours | **Risk:** medium (format edge cases)

**Task 1.3: RustMesh data structure**
- `struct RustMesh { points: ndarray, connectivity: ndarray, elem_type: BitVec, ... }`
- Owned mesh in Rust; minimal copying to Python
- CCW orientation check + enforce
- **Test:** Roundtrip fort.14 → RustMesh → assertions
- **Effort:** 4 hours | **Risk:** low

**Wave 1 Effort:** ~20 hours | **Blocker Resolution:** None (self-contained)

---

### Wave 2: Adjacency Backend (Quad-Edge in Rust) (Week 3–4)

**Task 2.1: Port quad-edge from Python to Rust**
- Translate `build_quadegg_from_connectivity()` algorithm (Python 340 LOC → Rust ~500 LOC)
- Use `ndarray` for edge storage: `Array2::<i32>::zeros((n_edges, 4))`
- 2-phase construction: create edges, pair opposites + assign next
- Handle boundary sentinels (-1)
- Mixed-element (tri + quad) support via `_elem_type` mask
- **Test:** Load annulus; verify `n_edges`, edge count per vertex (expect 2–6)
- **Effort:** 16 hours | **Risk:** medium (pointer logic, off-by-one errors)

**Task 2.2: Adjacency converters (quad-edge → standard formats)**
- `to_edge2vert()` — canonical sorted edges
- `to_elem2edge()` — element incident edges
- `to_vert2edge()` — vertex incident edges
- `to_vert2elem()` — vertex incident elements
- `to_edge2elem()` — edge adjacent elements
- Ensure bit-identical to Python reference (canonical form)
- **Test:** Equivalence vs Phase 008 Python output on all 4 fixtures
- **Effort:** 12 hours | **Risk:** medium (canonical form, edge ordering)

**Task 2.3: PyO3 export (make adjacencies callable from Python)**
- Expose RustMesh methods as `#[pyo3(name = "...")]` functions
- Return ndarray / dict types marshalled correctly
- **Test:** `mesh._rust_mesh.get_edge2vert()` returns np.ndarray[n_edges, 2]
- **Effort:** 6 hours | **Risk:** low (PyO3 straightforward)

**Wave 2 Effort:** ~34 hours | **Blocker:** Must match Phase 008 Python output bit-for-bit

---

### Wave 3: Quality Analysis & Basic Queries (Week 5–6)

**Task 3.1: Signed area computation (element quality)**
- Implement shoelace formula per element (vectorized)
- Return ndarray[n_elems]
- **Test:** Verify median quality matches Python reference ±0.01%
- **Effort:** 4 hours | **Risk:** low

**Task 3.2: Boundary detection + orientation**
- _ensure_ccw_orientation() — enforce counter-clockwise
- Boundary element detection (edge with only 1 adjacent element)
- **Test:** All elements CCW-oriented; boundary edges identified
- **Effort:** 6 hours | **Risk:** low

**Task 3.3: Vertex/element query methods**
- `get_vertex_edges(v: usize) → Vec<usize>` (return Vec → Python list or set)
- `get_vertex_elements(v: usize) → Vec<usize>`
- `get_element_vertices(e: usize) → [usize; 4]` (padded, return ndarray)
- **Test:** Roundtrip: query all vertices → verify sets complete + non-overlapping
- **Effort:** 4 hours | **Risk:** low

**Wave 3 Effort:** ~14 hours | **Blocker:** None (build on Wave 2)

---

### Wave 4: Skeletonization Engine (Week 7–10)

**Task 4.1: Layer-by-layer medial axis extraction (_skeletonize_medial_axis)**
- Core algorithm: remove boundary edges iteratively, track layers
- Per-layer: boundary elements (BE), interior elements (IE), outer verts (OV), inner verts (IV)
- Implement pruning (if configured): remove layers below quality threshold
- Return: `struct SkeletonLayers { layers: Vec<Layer> }` where Layer = { be, ie, ov, iv }
- **Test:** annulus → expect ~3–5 layers; verify all elements classified exactly once
- **Effort:** 24 hours | **Risk:** high (algorithm complexity, invariants)

**Task 4.2: Boundary extraction (element → OE mapping)**
- OE = oriented element: boundary element with outer-facing normal
- Inner element = interior, not on boundary
- Implement orientation logic (consistent outward normal)
- **Test:** All OE have exactly 1 boundary edge; IE have 0
- **Effort:** 8 hours | **Risk:** medium (orientation, edge cases)

**Task 4.3: Vertex classification per layer**
- OV = outer vertex (on boundary of layer k)
- IV = inner vertex (not on boundary)
- Implement layer-aware vertex traversal
- **Test:** Verify OV ⊂ Vertices(Layer k), IV disjoint from boundary
- **Effort:** 6 hours | **Risk:** medium

**Task 4.4: PyO3 export (return SkeletonLayers to Python)**
- Marshal Vec<Layer> → Python dict (same shape as Phase 008 Python output)
- `get_skeleton_layers() → { 'OE': list, 'IE': list, 'OV': list, 'IV': list, ... }`
- **Test:** Verify structure matches Python; all lists populated correctly
- **Effort:** 4 hours | **Risk:** low

**Wave 4 Effort:** ~42 hours | **Blocker:** Skeletonization is the crown jewel; must nail invariants

---

### Wave 5: Fort.14 Writer + Mesh Mutation Prep (Week 11–12)

**Task 5.1: Fort.14 writer (ASCII output)**
- Implement `write_fort14(path: &str)` method
- Write header (NE, NP, nbnd)
- Element connectivity in original order (preserve fort.14 format)
- Coordinates to 8-digit precision
- **Test:** Roundtrip: read → write → read; verify byte-for-byte (or semantic equivalence)
- **Effort:** 6 hours | **Risk:** low

**Task 5.2: 2dm reader (optional, Phase 9.1.1)**
- Parse 2dm format (simpler than fort.14; similar logic)
- Support for quadrilateral/triangular meshes
- **Test:** Load a 2dm fixture; verify n_verts, n_elems
- **Effort:** 6 hours | **Risk:** low (if deferred, no impact)

**Task 5.3: Prep for mesh mutation (Phase 10)**
- Add mutable handles to RustMesh (e.g., `add_element()`, `remove_element()`)
- Skeleton NOT implemented (too risky for Phase 9); stubs only
- Document pre-conditions (mesh must be valid CCW, no redundant edges)
- **Effort:** 4 hours | **Risk:** low (stubs, no implementation)

**Wave 5 Effort:** ~16 hours | **Blocker:** None (I/O is optional fallback)

---

### Wave 6: Integration Testing (Week 13–15)

**Task 6.1: Port pytest suite to call Rust backend**
- Modify conftest.py to optionally construct CHILmesh with `topology_backend='rust'`
- Run all 439 tests via Python API → Rust implementation
- **Target:** 439/439 PASS (no modifications to test code)
- **Effort:** 8 hours | **Risk:** medium (edge cases, tolerance issues)

**Task 6.2: Equivalence audit (Rust vs Python reference)**
- Run both backends on all 4 fixtures (annulus, donut, block_o, structured)
- Compare: adjacencies (bit-identical), skeleton (layer counts, element counts), quality stats
- **Target:** Bit-identical for adjacencies; ±0.1% for quality metrics
- **Effort:** 6 hours | **Risk:** low (automated comparison)

**Task 6.3: Performance baseline (Rust vs Python)**
- Benchmark Rust on WNAT_Hagen: fast_init, full_init, quality_analysis, query_latency
- Target: 1.1× EdgeMap (3.5s full_init)
- Profile bottlenecks if > 1.2× target
- **Effort:** 4 hours | **Risk:** medium (if unoptimized, may exceed target)

**Wave 6 Effort:** ~18 hours | **Blocker:** Integration must work; else rework Rust API

---

### Wave 7: Optimization + Profiling (Week 16–18)

**Task 7.1: Profile Rust implementation (if performance below target)**
- Use `perf` / `cargo flamegraph` to identify hot paths
- Likely culprits: allocations in loop, unnecessary copying, inefficient traversal
- **Expected wins:** 10–20% improvement via arena allocators, simd if applicable
- **Effort:** 8 hours (if needed; skip if already 1.1×) | **Risk:** medium (profiling learning curve)

**Task 7.2: Vectorization (optional)**
- NumPy-style operations on adjacency data (element-wise quality, batch vertex queries)
- Measure benefit; only implement if >5% improvement expected
- **Effort:** 6 hours | **Risk:** low (optional)

**Task 7.3: Compile-time optimizations**
- LTO (link-time optimization) in Cargo.toml
- Codegen optimizations (opt-level = 3)
- Binary size reduction (if applicable)
- **Effort:** 2 hours | **Risk:** low

**Wave 7 Effort:** ~16 hours (conditional) | **Blocker:** None (optimization can slip if time-boxed)

---

### Wave 8: Documentation + Release Prep (Week 19–24)

**Task 8.1: Rust module documentation**
- Doc comments on all public functions (examples, panics, complexity)
- Design notes on quad-edge algorithm, skeletonization invariants
- Performance characteristics documented
- **Effort:** 8 hours | **Risk:** low

**Task 8.2: Python API documentation (no changes, just note Rust backend)**
- Update docstrings to clarify "uses Rust backend for performance"
- Benchmark comparison table in README
- Migration guide: "Phase 008 quadedge Python → Phase 009 Rust"
- **Effort:** 4 hours | **Risk:** low

**Task 8.3: CI/CD + build automation**
- Add Rust to GitHub Actions CI (cargo build, cargo test on Linux/macOS/Windows)
- maturin/PyO3 build script (creates .so / .pyd wheels)
- Binary wheels on PyPI (optional Phase 9.1.1)
- **Effort:** 6 hours | **Risk:** medium (wheel building, platform-specific)

**Task 8.4: Release notes + migration**
- Document performance improvements (3.5s vs 3.19s EdgeMap baseline; 13% faster than Python)
- List breaking changes (none expected; API identical)
- Rollback plan (fallback to Python if Rust build fails)
- **Effort:** 4 hours | **Risk:** low

**Task 8.5: User testing + feedback loop**
- Share Phase 009 branch for early testing (MADMESHR, ADMESH, other users)
- Incorporate feedback (bugs, performance edge cases)
- **Effort:** 8 hours (ongoing) | **Risk:** low (async, doesn't block release)

**Wave 8 Effort:** ~30 hours | **Blocker:** None (release-track, can slip to Phase 9.1)

---

## Dependency Graph & Waves

```
Wave 1: I/O Foundation (Weeks 1–2)
  ├─ Task 1.1: Rust scaffolding
  ├─ Task 1.2: Fort.14 parser
  └─ Task 1.3: RustMesh struct
        ↓
Wave 2: Quad-Edge Backend (Weeks 3–4)
  ├─ Task 2.1: Port quad-edge
  ├─ Task 2.2: Adjacency converters
  └─ Task 2.3: PyO3 export
        ↓
Wave 3: Queries & Quality (Weeks 5–6)
  ├─ Task 3.1: Signed area
  ├─ Task 3.2: Boundary + orientation
  └─ Task 3.3: Query methods
        ↓
Wave 4: Skeletonization (Weeks 7–10) ← CRITICAL PATH
  ├─ Task 4.1: Layer extraction
  ├─ Task 4.2: OE classification
  ├─ Task 4.3: Vertex classification
  └─ Task 4.4: PyO3 export
        ↓
Wave 5: I/O + Prep (Weeks 11–12)
  ├─ Task 5.1: Fort.14 writer
  ├─ Task 5.2: 2dm reader (optional)
  └─ Task 5.3: Mutation stubs
        ↓
Wave 6: Integration Testing (Weeks 13–15)
  ├─ Task 6.1: Port test suite
  ├─ Task 6.2: Equivalence audit
  └─ Task 6.3: Perf baseline
        ↓
Wave 7: Optimization (Weeks 16–18) [conditional]
  ├─ Task 7.1: Profiling
  ├─ Task 7.2: Vectorization
  └─ Task 7.3: Compile opts
        ↓
Wave 8: Release (Weeks 19–24)
  ├─ Task 8.1: Rust docs
  ├─ Task 8.2: Python docs
  ├─ Task 8.3: CI/CD
  ├─ Task 8.4: Release notes
  └─ Task 8.5: User testing
```

**Critical Path:** Wave 4 (skeletonization) — 42 hours, most complex algorithm. Everything else is mechanical translation.

**Parallel Opportunities:**
- Waves 2 & 3 can overlap (separate modules)
- Wave 5 I/O can start once Wave 1 settled
- Wave 7 profiling can run concurrently with Wave 6 testing

---

## Risk Assessment

| Risk | Severity | Mitigation |
|------|----------|-----------|
| Skeletonization algorithm bugs (Wave 4) | HIGH | Extensive testing on all 4 fixtures; visual inspection of layer boundaries; comparison vs Python reference |
| FFI performance overhead (Waves 2–3) | MEDIUM | PyO3 overhead ~5-10%; acceptable. Profile if exceeds target (Wave 7). |
| Rust learning curve (all waves) | MEDIUM | Start with straightforward tasks (Wave 1); tackle algorithm complexity (Wave 4) after team familiar |
| Binary wheel distribution | LOW | maturin handles most complexity; fallback to source builds if needed |
| Floating-point precision (quality metrics) | LOW | Use f64 throughout; tolerance ±1e-6 for comparisons; document |
| Platform-specific bugs (Windows/macOS) | MEDIUM | CI must test all three (GitHub Actions matrix). Rust LLVM usually portable. |

---

## Resource Estimate

| Category | Effort | Notes |
|----------|--------|-------|
| **Wave 1–5 (Core implementation)** | ~130 hours | Mostly mechanical translation + I/O parsing |
| **Wave 6 (Integration + testing)** | ~18 hours | Leverages existing pytest suite; mostly configuration |
| **Wave 7 (Optimization)** | ~16 hours (conditional) | Only if performance target missed; likely not needed |
| **Wave 8 (Release prep)** | ~30 hours | Documentation, CI/CD, rollout |
| **TOTAL** | **194 hours (~24 weeks at 8h/week)** | Or 12 weeks at 16h/week (full-time effort) |

**Staffing:** 1–2 Rust developers (or 1 developer + peer review from Python architect for correctness). Python wrapper integration handled by Python team.

**Tools:**
- Rust 1.70+ (recent stable)
- PyO3 0.20+ (Python binding framework)
- ndarray crate (multi-dimensional arrays)
- cargo (build system)
- perf / cargo-flamegraph (profiling, Wave 7)

---

## Success Criteria

| Criterion | Target | Acceptance |
|-----------|--------|-----------|
| **Performance** | Full init ≤ 3.5s WNAT_Hagen | Measured via benchmark_quadegg_variants.py (already written) |
| **Correctness** | 439/439 tests PASS | No modifications to test code; 100% compatibility |
| **Equivalence** | Adjacencies bit-identical; quality ±0.1% | Automated audit vs Python reference on all 4 fixtures |
| **Memory** | ≤ 25% vs EdgeMap baseline | Peak memory measured via tracemalloc (same as Phase 008) |
| **API compatibility** | 100% (no breaking changes) | All existing Python code runs unchanged |
| **CI/CD** | Green on Linux + macOS + Windows | GitHub Actions matrix passing |

---

## Optional Enhancements (Phase 9.1+)

1. **Binary wheels on PyPI** — Distribute pre-compiled .so / .pyd (saves compilation on install)
2. **2dm reader** — Full support for 2dm format (similar to fort.14, low priority)
3. **Point location (find_element)** — Spatial indexing for O(log n) element queries (instead of O(n) linear search)
4. **SIMD vectorization** — Shoelace formula, batch operations (if profiling shows benefit)
5. **Phase 10 (Mesh mutation)** — Add element insertion/deletion (requires mutable skeletonization; risky, defer)

---

## Conclusion

Full Rust port feasible in 24 weeks (1 developer) or 12 weeks (2 developers). Skeletonization (Wave 4) is the critical path and most complex piece; all other waves are straightforward translation.

**Python API wrapper (50 LOC) maintains 100% compatibility — users see no breaking changes. Rust backend is entirely internal.**

Risk is **medium** (algorithm correctness, not engineering). Mitigation: rigorous equivalence testing + profiling.

**Ready to commit to Phase 9 upon user approval?**

---

**Last Updated:** 2026-05-22  
**Document Version:** 1.0 (DRAFT)  
**Author:** Claude Code  
**Status:** Awaiting user decision on Rust vs C++ and scope approval
