# Context: Phase 009 — Rust Backend Port (Implementation Decisions)

**Phase:** 009-rust-backend-port  
**Status:** LOCKED — Ready for planning  
**Date:** 2026-05-22

---

## Domain

Rewrite CHILmesh mesh engine (I/O, adjacency, skeletonization, mutation, indexing) in Rust. Expose via PyO3 FFI. Keep Python API 100% backward compatible.

Target: Full init ≤3.5s WNAT_Hagen. All 439 tests PASS unchanged.

---

## Locked Requirements (from SPEC.md)

**SPEC.md is authoritative.** Downstream agents: read `specs/009-rust-backend-port/009-SPEC.md` before planning. 10 functional reqs (F-001 to F-010), 5 non-functional (N-001 to N-005), 18 acceptance criteria.

Key gates:
- Full init ≤ 3.5s WNAT_Hagen (1.1× EdgeMap Python baseline)
- Performance variance < 5% across 3 trials
- All 439 tests PASS unmodified
- Adjacency bit-identical to Python reference (canonical form)
- Mesh mutation + spatial indexing in Phase 009 (not deferred)

---

## Implementation Decisions

### 1. Rust Module Architecture: Multi-module

```
src/chilmesh_core/
├── lib.rs                 # PyO3 FFI + module imports
├── io.rs                  # Fort.14/2dm parser + writer
├── adjacency.rs           # Quad-edge construction + converters
├── skeletonization.rs     # Layer extraction + quality metrics
├── queries.rs             # Vertex/element lookups
├── mutation.rs            # add_element + remove_element
└── errors.rs              # Error type definitions
```

**Why:** Mirrors Python structure. Clear separation of concerns. Allows parallel dev on Wave 2–5. Simpler imports.

**Impact:** Slightly longer compilation (~5–10s per Wave), offset by clarity + concurrent work.

---

### 2. Error Handling: Custom enum + thiserror crate

```rust
use thiserror::Error;

#[derive(Error, Debug)]
pub enum RustMeshError {
    #[error("Parse error: {0}")]
    ParseError(String),
    
    #[error("Invalid geometry: {0}")]
    InvalidGeometry(String),
    
    #[error("Vertex out of bounds: {0}")]
    OutOfBounds(String),
}
```

PyO3 bridges:
- `ParseError` → Python `OSError`
- `InvalidGeometry`, `OutOfBounds` → Python `ValueError`
- Error message includes context (line, input description)

**Why:** Type-safe, composable, explicit semantics. Matches Python API design.

---

### 3. Spatial Indexing: Quadtree (Phase 009), k-d tree deferred (Phase 10+)

Quadtree construction: recursive axis-aligned subdivision. O(log n) average case.

**Why:** Simpler impl, fewer edge cases, works on annulus/donut/block_o. If profiling shows bottleneck → defer k-d tree to Phase 10.

---

### 4. Mesh Mutation: Conservative validation + full recompute

`add_element([v0, v1, v2])`  
- Validate: no self-intersection, vertices in bounds, CCW orientation
- Insert into mesh
- Recompute adjacency (O(n), acceptable)
- Invalidate skeletonization layers (user calls `compute_layers=True` to rebuild)

`remove_element(elem_id)`  
- Check: no orphaned vertices
- Mark edges as boundary (-1 sentinel in Edge2Elem)
- Recompute adjacency
- Invalidate layers

**Why:** Simplifies invariant checking. Recomputation still <100ms on large meshes.

---

### 5. FFI & Marshalling: PyO3 ndarray + native PyO3 types

Return ndarrays: `PyArray2::as_slice_mut()` (zero-copy)  
Return dicts: `PyDict`  
Return sets: `PyFrozenSet` (immutable) or `PySet` (if mutable)

**Overhead:** <5% vs Python baseline (measured in Phase 008 quad-edge).

**Why:** Standard PyO3 pattern, well-tested, minimal custom code.

---

### 6. Python Wrapper: Pure delegation

```python
class CHILmesh:
    def __init__(self, ...):
        self._rust_mesh = chilmesh_core.RustMesh(...)
    
    def signed_area(self):
        return self._rust_mesh.signed_area()
```

No business logic in Python. Direct method forwarding.

**Why:** Zero overhead, clear ownership, easy audit.

---

### 7. Testing: Python integration only, Rust unit tests for complex algos

All 439 pytest tests unchanged (call Python API).  
Rust unit tests for:
- Quadtree insertion + point location
- Skeletonization layer extraction (OE/IV/OV classification)
- Adjacency canonical form (vs Python reference)

No separate Rust binaries. Tests linked into libpython.

**Why:** Single test suite. Python passing = no regressions.

---

## Code Context

**Key reference implementations (Python → Rust):**

- `src/chilmesh/mesh_topology_quadegg.py` (340 LOC) — Quad-edge canonical form (use as Rust reference)
- `src/chilmesh/CHILmesh.py._build_adjacencies()` — Adjacency API (replicate in Rust)
- `src/chilmesh/CHILmesh.py._skeletonize_medial_axis()` — Skeletonization engine (critical algorithm)
- `src/chilmesh/examples.py` — Test fixtures (annulus, donut, block_o, structured)
- `scripts/benchmark_quadegg_variants.py` — Benchmark framework for performance validation

**Reusable patterns:**
- Canonical-form adjacency comparator (Phase 008, in test suite)
- Fixture parametrization over 4 meshes (pytest pattern)
- Fort.14 round-trip testing (existing test suite)

---

## Canonical Refs

**MUST READ before planning:**
- `.planning/ROADMAP.md` — Phase 008/009 scope + status
- `specs/009-rust-backend-port/009-SPEC.md` — Locked requirements (ambiguity ≤ 0.20)
- `.planning/PHASE_9_RUST_PORT_PLAN.md` — 8-wave breakdown, effort estimates, risks
- `src/chilmesh/mesh_topology_quadegg.py` — Quad-edge reference impl (340 LOC)
- `scripts/benchmark_quadegg_variants.py` — Performance measurement framework

**Optional context:**
- `.planning/008-DECISION.md` — Why quad-edge adopted (13% faster than half-edge)
- `output/benchmark.json` — WNAT_Hagen baseline times (EdgeMap 3.19s, QuadEdge 4.91s)

---

## Deferred Ideas

- **K-d tree for spatial indexing** — Phase 10+ if quadtree underperforms
- **Binary wheel distribution** — Phase 9.1 (users build locally for now)
- **C++ variant** — Post-Phase-9 if Rust team requests
- **Half-edge topology in Rust** — Research-only, deferred unless quad-edge insufficient

---

## Success Gate for Planning

Planner reads:
1. SPEC.md (requirements locked)
2. CONTEXT.md (implementation decisions locked)
3. Phase 9 plan + 8-wave breakdown

Planner outputs PLAN.md with:
- Detailed task breakdown per wave
- Wave 1–2 critical path (I/O + adjacency)
- Dependency analysis
- Risk mitigation (compilation time, PyO3 FFI pitfalls)
- Effort allocation per wave

---

**Last Updated:** 2026-05-22  
**Document Version:** 1.0 (LOCKED)
