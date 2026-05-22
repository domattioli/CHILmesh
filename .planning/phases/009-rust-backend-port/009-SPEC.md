# SPEC: Phase 009 — Rust Backend Port with Python API Wrapper

**Phase:** 009-rust-backend-port  
**Title:** Full Rust Backend Port (Quad-Edge Topology, Mesh Mutation, Spatial Indexing)  
**Created:** 2026-05-22  
**Ambiguity Score:** 0.172 (gate: ≤ 0.20) ✓  
**Status:** LOCKED

---

## Goal

Rewrite CHILmesh core (mesh I/O, adjacency, skeletonization, mutation, spatial indexing) in Rust; expose via PyO3 FFI maintaining 100% Python API compatibility. Deliver production-ready Rust backend with all functionality except plotting and ADMESH integration (separate repo).

**Performance Target:** Full init ≤ 3.5s WNAT_Hagen (1.1× EdgeMap Python baseline, 13% improvement over half-edge).

**Scope Expansion (Phase 009 vs draft):** Includes mesh mutation (add/remove elements) + spatial indexing (point location) — ~230 hours total.

## In Scope

1. **Mesh I/O (Fort.14 + 2dm)**
   - Fort.14 ASCII reader (parsing, validation)
   - Fort.14 ASCII writer (roundtrip-compatible)
   - 2dm reader (triangular + quadrilateral support)
   - Binary serialization (optional, Phase 9.1)

2. **Adjacency Backend (Quad-Edge)**
   - `build_quadegg_from_connectivity()` — 2-phase O(n) construction
   - All 6 converters: Edge2Vert, Edge2Elem, Elem2Edge, Vert2Edge, Vert2Elem
   - Boundary handling (-1 sentinels)
   - Mixed-element support (triangles + quads)
   - Bit-identical output to Python reference (canonical form)

3. **Skeletonization (Medial Axis Extraction)**
   - Layer-by-layer boundary removal algorithm
   - OE (oriented element) classification
   - IV (inner vertex) / OV (outer vertex) identification per layer
   - Quality pruning (remove layers below threshold)
   - Complexity: O(n) per layer × L layers

4. **Quality Analysis**
   - `signed_area()` per element (shoelace formula)
   - Median quality across all elements
   - Layer-wise quality statistics
   - CCW orientation enforcement + validation

5. **Mesh Queries**
   - `get_vertex_edges(v) → Vec<usize>` (all incident edges)
   - `get_vertex_elements(v) → Vec<usize>` (all incident elements)
   - `get_element_vertices(e) → [usize; 4]` (padded quad)

6. **Mesh Mutation (NEW — Phase 009 commitment)**
   - `add_element(verts: [usize; 3|4], elem_type) → elem_id`
   - `remove_element(elem_id)`
   - Adjacency invalidation + recomputation
   - Skeletonization recomputation (expensive but necessary)
   - Invariant: Mesh remains CCW-oriented, no self-intersections
   - Precondition: User must validate input mesh before mutation

7. **Spatial Indexing (NEW — Phase 009 commitment)**
   - `find_element(x: f64, y: f64) → Option<elem_id>` — point location
   - Quadtree or k-d tree implementation (TBD during planning)
   - O(log n) lookup complexity
   - Supports all 4 fixtures (annulus, donut, block_o, structured)
   - Handles boundary cases (point outside domain → None)

8. **Python Wrapper (FFI + API)**
   - PyO3 module `chilmesh_core` exposing RustMesh
   - Python CHILmesh class holds opaque Rust handle
   - Method delegation (e.g., `mesh.signed_area()` → `_rust_mesh.signed_area()`)
   - ndarray/dict marshalling (adjacencies returned as native Python types)
   - Error translation (Rust Result → Python exceptions with context)
   - 100% backward compatibility — all 439 tests PASS unmodified

## Out of Scope

- ✗ Plotting/visualization (stays Python, `src/chilmesh/plot_utils.py` unchanged)
- ✗ ADMESH integration (ADMESH is separate repo; does not touch CHILmesh core)
- ✗ MADMESHR bridge (stays Python; wraps Rust mesh)
- ✗ Binary wheels (source-only; users build locally. Phase 9.1+ adds pre-built wheels)
- ✗ Language porting to C++ (Rust chosen; C++ optional post-launch if needed)
- ✗ Half-edge topology variant (Rust port uses quad-edge only)
- ✗ Public API beyond `topology_backend` kwarg (internal feature, unchanged)

## Constraints

1. **Backward compatibility:** All 439 existing tests PASS without modification. Python API identical pre/post-Phase-009.

2. **Performance ceiling:** Full init ≤ 3.5s WNAT_Hagen (1.1× Python EdgeMap baseline). If profiling shows > 1.2× overage, implement optimization (Wave 7).

3. **Memory overhead:** ≤ 25% vs EdgeMap baseline (same as half-edge Python).

4. **Benchmark variance:** < 5% across 3 trials per operation.

5. **Mixed-element handling:** Triangles remain padded to 4-column arrays with `_elem_type` mask. No spurious edges.

6. **Error semantics:** Rust panics become Python exceptions with context. Invalid input (e.g., self-intersecting element) raises `ValueError` or `RuntimeError` (TBD).

7. **Fort.14 format:** Roundtrip reading → writing must produce ASCII output compatible with existing tools. Binary precision: f64 throughout (no f32).

## Requirements

### Functional Requirements

| Req | Current State | Target State | Acceptance Criterion |
|-----|---|---|---|
| **F-001: Quad-edge I/O** | Python quad-edge module exists (Phase 008) | Fort.14 reader/writer in Rust + Python wrapper | `CHILmesh.read_from_fort14(path)` succeeds; `mesh.write_fort14(path)` roundtrips |
| **F-002: Adjacency converters** | Python converters exist (Phase 008) | Rust converters (6 types) + FFI export | Bit-identical output vs Python reference; 0 failures on canonical-form audit |
| **F-003: Skeletonization** | Implemented in Python | Rewritten in Rust; layer-by-layer extraction | All 439 tests PASS; median quality matches Python ±0.1% |
| **F-004: Quality metrics** | Python `signed_area()` | Rust `signed_area()` + layer stats | Answers match Python ±1e-6 per element |
| **F-005: Vertex queries** | Python dict-based (Vert2Edge, Vert2Elem) | Rust Vec-based + FFI to Python | `get_vertex_edges(v)` returns list; O(1) for typical vertices |
| **F-006: Mesh mutation** | Not implemented (Phase 2 designed only) | Add/remove elements in Rust + recompute adjacency | `add_element([v0, v1, v2]) → elem_id`; mesh remains valid |
| **F-007: Spatial indexing** | None (linear search fallback) | Quadtree/k-d tree in Rust | `find_element(x, y)` O(log n); returns elem_id or None |
| **F-008: Python wrapper** | None (will create) | PyO3 module + CHILmesh class delegates | All Python imports work; `type(mesh._rust_mesh).__module__ == "chilmesh_core"` |
| **F-009: Error handling** | None (will add) | Rust Result → Python exceptions with context | Invalid mesh raises `ValueError`; I/O errors raise `OSError` |
| **F-010: Boundary handling** | Python (-1 sentinels) | Rust (-1 sentinels) | Edge2Elem boundary entries = -1; unmarshal to Python dict |

### Non-Functional Requirements

| Req | Metric | Acceptance |
|-----|--------|-----------|
| **N-001: Performance** | Full init WNAT_Hagen | ≤ 3.5s (1.1× Python EdgeMap) |
| **N-002: Memory** | Peak memory vs EdgeMap | ≤ 25% increase |
| **N-003: Variance** | Std dev across 3 trials | < 5% per operation |
| **N-004: Compilation** | Cargo build time | < 60s (debug), < 120s (release) |
| **N-005: Binary size** | Compiled .so/.pyd | < 50 MB (unstripped) |

## Boundaries

### In Scope (Yes, Phase 009 delivers)
✓ Full Rust rewrite of core mesh (I/O, adjacency, skeletonization, mutation, indexing)  
✓ Python API wrapper (100% backward compatible)  
✓ All 439 tests PASS unmodified  
✓ Bit-identical adjacency vs EdgeMap reference  
✓ Mesh mutation + spatial indexing  
✓ Source-only distribution (users compile locally)

### Out of Scope (No, deferred or adjacent)
✗ Plotting/visualization  
✗ ADMESH/MADMESHR integration  
✗ Binary wheels (PyPI pre-built; Phase 9.1+)  
✗ C++ variant  
✗ Half-edge topology  

**Boundary reasoning:** Phase 009 is "core topology + operations in Rust." Visualization, downstream integration, and distribution mechanics are separate concerns. Mutation + indexing added to Phase 009 per explicit request.

## Acceptance Criteria

### Pass/Fail Checkboxes

- [ ] `src/chilmesh_core/` Rust crate exists with `Cargo.toml`, `lib.rs`, PyO3 scaffolding
- [ ] Fort.14 reader implemented; loads annulus, donut, block_o, structured without error
- [ ] Fort.14 writer implemented; roundtrip read → write → read produces byte-compatible ASCII
- [ ] Quad-edge adjacency converters ported from Python; output bit-identical on all 4 fixtures
- [ ] PyO3 module `chilmesh_core` builds and imports successfully: `python -c "from chilmesh_core import RustMesh"`
- [ ] Python wrapper class (CHILmesh) delegates to Rust; `topology_backend` parameter accepted but internally unused (Rust is now default)
- [ ] Skeletonization layer extraction working in Rust; layer counts match Python ±1 layer
- [ ] Quality metrics (`signed_area()`) match Python reference ±0.1% per element
- [ ] Mutation operations (`add_element`, `remove_element`) work; adjacency recomputed correctly
- [ ] Spatial indexing (`find_element`) O(log n); returns correct elem_id for point in domain
- [ ] All 439 pytest tests PASS with `CHILMESH_TOPOLOGY_BACKEND=rust` (env var, Python still recognizes it)
- [ ] Equivalence audit: adjacency outputs (Edge2Vert, Elem2Edge, etc.) bit-identical to Python reference on all 4 fixtures
- [ ] Performance: Full init on WNAT_Hagen ≤ 3.5s (measured via benchmark_quadegg_variants.py extended for Rust)
- [ ] Memory: Peak memory ≤ 25% vs EdgeMap baseline
- [ ] Variance: < 5% across 3 trials per operation
- [ ] CI/CD: GitHub Actions matrix (Linux + macOS + Windows) all green
- [ ] Error handling: Invalid input raises Python exception with context; no panics visible to user

## Success Metrics

1. **Code complete** — All 10 functional requirements (F-001 to F-010) verified
2. **Tests passing** — 439/439 existing tests PASS; new integration tests for mutation + indexing
3. **Benchmarked** — 4 operations measured: fast_init, full_init, quality_analysis, query_latency. Target ≤ 3.5s full_init
4. **Decided** — Performance analysis in DECISION.md: "Rust quad-edge meets target, ready for production" OR "optimization needed (Wave 7)" OR "Fallback to Python recommended"

## Ambiguity Report

**Final Score: 0.172** (gate: ≤ 0.20) ✓

| Dimension | Score | Minimum | Status |
|---|---|---|---|
| Goal Clarity | 0.85 | 0.75 | ✓ Exceeded |
| Boundary Clarity | 0.88 | 0.70 | ✓ Exceeded |
| Constraint Clarity | 0.80 | 0.65 | ✓ Exceeded |
| Acceptance Criteria | 0.75 | 0.70 | ✓ Met |

**Interview findings:**
- Scope fully specified (quad-edge + mutation + indexing, no half-edge fallback, source-only distribution)
- Error handling strategy locked (Result → Python exceptions with context)
- Performance ceiling locked (≤ 3.5s WNAT_Hagen)
- Backward compatibility guaranteed (all 439 tests, no API changes)

**Unresolved (deferred to planning):** Spatial indexing algorithm (quadtree vs k-d tree) — will be decided during Phase 009 planning based on complexity/performance tradeoffs.

---

## Next Step

Run `/gsd:discuss-phase 009` to lock implementation decisions:
- Rust architecture (module structure, error types, FFI bindings)
- Spatial indexing algorithm choice (quadtree vs k-d tree)
- Mutation invariant enforcement (validation, error recovery)
- Integration strategy (Python wrapper delegation pattern)

This SPEC.md is now locked. discuss-phase will treat all requirements above as frozen.

---

**Last Updated:** 2026-05-22  
**Document Version:** 1.0  
**Status:** LOCKED (ready for discuss-phase)
