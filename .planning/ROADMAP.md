# CHILmesh Roadmap — Phases 008–009

**Project:** CHILmesh 0.2.0+ (Modernization)  
**Timeline:** 2026-05-22 onward  
**Maintained by:** gsd-spec-phase, gsd-discuss-phase, gsd-plan-phase

---

## Phase 008: Quad-Edge Topology Optimization ✅ COMPLETE

**Status:** COMPLETE (2026-05-22)  
**Goal:** Benchmark quad-edge adjacency backend vs EdgeMap and half-edge; decide on Phase 9 porting target.  
**Scope:** Quad-edge implementation (340 LOC), 104 unit tests, 4-backend benchmark on WNAT_Hagen.

**Decision:** ADOPT QUAD-EDGE for Phase 9+ — 13% faster than half-edge, simpler 2-phase construction.

**Deliverables:**
- ✅ `src/chilmesh/mesh_topology_quadegg.py` — Quad-edge backend
- ✅ `tests/test_quadegg_construction.py` — 72 construction tests
- ✅ `tests/test_quadegg_equivalence.py` — 32 equivalence tests
- ✅ `scripts/benchmark_quadegg_variants.py` — 4-backend benchmark
- ✅ `.planning/008-DECISION.md` — decision record with benchmark results
- ✅ Branch: `008-optimize-port-w-quad-edge`

---

## Phase 009: Rust Backend Port with Python API Wrapper (IN PROGRESS)

**Status:** DISCUSSING IMPLEMENTATION DECISIONS  
**Goal:** Rewrite CHILmesh core (mesh I/O, adjacency, skeletonization, mutation, spatial indexing) in Rust; expose via PyO3 FFI maintaining 100% Python API compatibility. Deliver production-ready Rust backend with all functionality except plotting and ADMESH integration (separate repo).

**Performance Target:** Full init ≤ 3.5s WNAT_Hagen (1.1× EdgeMap Python baseline).

**Scope Expansion (Phase 009 vs draft):** Includes mesh mutation (add/remove elements) + spatial indexing (point location) — ~230 hours total.

### In Scope

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

7. **Spatial Indexing (NEW — Phase 009 commitment)**
   - `find_element(x: f64, y: f64) → Option<elem_id>` — point location
   - Quadtree or k-d tree implementation (TBD during planning)
   - O(log n) lookup complexity

8. **Python Wrapper (FFI + API)**
   - PyO3 module `chilmesh_core` exposing RustMesh
   - Python CHILmesh class holds opaque Rust handle
   - Method delegation with 100% backward compatibility
   - All 439 tests PASS unmodified

### Out of Scope

- ✗ Plotting/visualization (stays Python)
- ✗ ADMESH integration (ADMESH is separate repo)
- ✗ Binary wheels (source-only; Phase 9.1+ adds pre-built wheels)
- ✗ Language porting to C++ (Rust chosen)
- ✗ Half-edge topology variant (Rust port uses quad-edge only)

**Canonical Refs:**
- `.planning/PHASE_9_RUST_PORT_PLAN.md` — comprehensive 8-wave breakdown
- `specs/009-rust-backend-port/009-SPEC.md` — locked requirements (ambiguity ≤ 0.20)

**Branch:** `claude/optimize-new-lang-bc7la` (created from 008-optimize-port-w-quad-edge)

---

## Phase 010+: Deferred

- Spatial indexing optimization (if Phase 9 underperforms)
- Binary wheel distribution
- C++ variant (optional)
- Half-edge topology variant (if Rust quad-edge insufficient)
- Additional mutation operations (mesh smoothing, subdivision)

---

**Last Updated:** 2026-05-22  
**Version:** 1.0
