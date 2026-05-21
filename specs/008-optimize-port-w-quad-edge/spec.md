# Feature Specification: Quad-Edge Data Structure Investigation & Benchmarking

**Feature Branch**: `008-optimize-port-w-quad-edge`
**Created**: 2026-05-21
**Status**: Draft (awaiting specify → clarify)
**Input**: User description: "Investigate quad-edge data structure as optimization alternative; benchmark against current EdgeMap and half-edge variants. Evaluate whether 4-connected edge topology yields measurable performance gains on WNAT_Hagen reference mesh."

**Related Issues**: #68 (Quad-edge investigation), #137 (Half-edge alternative), comparison study

---

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Library maintainer compares quad-edge to half-edge baseline (Priority: P1)

CHILmesh maintainer wants to compare the quad-edge (4-connected) data structure against the current EdgeMap baseline and the previously-benchmarked half-edge (DCEL) implementation to determine which topology backend yields the best performance on the WNAT_Hagen reference mesh.

**Why this priority**: Half-edge benchmark (Phase 007) showed 80% slowdown; quad-edge is a candidate for superior performance due to 4-connected symmetry. P1 blocks the adoption decision until all three variants are compared head-to-head on the same benchmark rig.

**Independent Test**: Run `python scripts/benchmark_quad_edge_variants.py` (modeled on half-edge benchmark) and confirm output table with four columns: EdgeMap (baseline), Half-Edge-v1, Half-Edge-v2 (vectorized), Quad-Edge (native). Table MUST show all four operations (fast_init, full_init, quality_analysis, query_latency) with median-of-3 trials and percent deltas vs. EdgeMap baseline.

**Acceptance Scenarios**:

1. **Given** the codebase on `008-optimize-port-w-quad-edge` branch, **When** `pytest -v` runs with `CHILMESH_TOPOLOGY_BACKEND=quadegg`, **Then** all existing tests pass with zero modifications.
2. **Given** a benchmark script runs quad-edge alongside the three existing variants (EdgeMap, HE-v1, HE-v2) on WNAT_Hagen (52.7k vertices), **When** results are collected, **Then** a markdown table is produced showing operation-by-operation deltas; quad-edge result MUST be present and either "meets or exceeds half-edge" or "slower than half-edge".
3. **Given** the quad-edge backend is selected, **When** mesh topology is built on any of the four fixtures (annulus, donut, block_o, structured), **Then** adjacency outputs (`Edge2Vert`, `Elem2Edge`, `Vert2Edge`, `Vert2Elem`, `Edge2Elem`) are bit-identical to the EdgeMap baseline.

---

### User Story 2 — Comparison drives adoption decision (Priority: P2)

After benchmarking quad-edge, maintainer has clear data to decide whether to:
- Adopt quad-edge over half-edge
- Stick with half-edge
- Archive both and optimize EdgeMap further
- Pursue language porting (Rust/C++) only after choosing the best topology

**Why this priority**: P1 produces the numbers; P2 is using them. Lower priority because it gates the next phase, not this one.

**Independent Test**: Read `output/benchmark.json` and verify it contains all four backend variants with complete metrics. Decision trigger MUST be unambiguous: "quad-edge fastest" → pursue adoption; "quad-edge slower than HE but faster than EdgeMap" → pursue HE; "quad-edge and HE both slower than EdgeMap" → archive topology work, focus elsewhere.

**Acceptance Scenarios**:

1. **Given** the benchmark data in `output/benchmark.json`, **When** a maintainer reads it, **Then** the quad-edge row is complete and comparable to all three other backends on the same hardware.
2. **Given** analysis of the results, **When** written into `.planning/008-DECISION.md`, **Then** the decision is data-driven (cite exact numbers from JSON) and unambiguous.

---

### User Story 3 — Quad-edge implementation guides future ports (Priority: P3)

If quad-edge outperforms half-edge, the pure-Python quad-edge implementation serves as a reference implementation for a potential Rust/C++ port in Phase 9+. Implementation clarity (well-documented construction algorithm, clear invariants) guides porting effort.

**Why this priority**: Deferred to post-analysis. Only becomes relevant if quad-edge wins the comparison.

**Independent Test**: Read `src/chilmesh/mesh_topology_quadegg.py` and verify construction algorithm is documented with complexity analysis (expected O(n) like half-edge). Code MUST have docstrings explaining the four-field tuple structure per half-edge precedent.

**Acceptance Scenarios**:

1. **Given** quad-edge implementation, **When** reviewed for clarity, **Then** construction algorithm is documented and the four-field-per-edge tuple is explained.
2. **Given** a future porting effort, **When** the pure-Python reference is read, **Then** the logic is clear enough to port to Rust/C++ without guessing.

---

### Edge Cases

- **Mixed-element padded triangles**: Like half-edge, quad-edge MUST handle padded triangles correctly without storing spurious edges for the padding slot. Reuse existing `_elem_type` mask.
- **Boundary edges**: Quad-edge natively supports undefined neighbors; converter MUST preserve `-1` sentinel for `Edge2Elem` output compatibility.
- **Disconnected meshes**: Quad-edge construction MUST handle multiple connected components independently.
- **Degenerate elements**: Run `_ensure_ccw_orientation()` before constructing quad-edge pointers (same precondition as half-edge).
- **Environment variable mis-set**: `CHILMESH_TOPOLOGY_BACKEND=quadegg` (typo or misspelling) MUST raise `ValueError` at mesh construction, not silently fall back.
- **Memory measurement variance**: Use `tracemalloc` snapshots; record median-of-3 per NFR-003.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST provide a quad-edge topology backend selectable via either a constructor argument (`topology_backend='quadegg'`) or the environment variable `CHILMESH_TOPOLOGY_BACKEND`. Quad-edge stored as `ndarray[n_edges, 4]` with columns `[edge_id, next_cw, next_ccw, opposite_edge_id]` (or canonical form per Wikipedia quad-edge definition); all integer-index pointers, no per-edge Python objects.
- **FR-002**: System MUST keep public CHILmesh API surface unchanged (only new optional `topology_backend` kwarg). No downstream import paths break.
- **FR-003**: System MUST produce bit-identical adjacency outputs (`Edge2Vert`, `Elem2Edge`, `Vert2Edge`, `Vert2Elem`, `Edge2Elem`) between EdgeMap, half-edge, and quad-edge backends on every built-in fixture. Canonical-form comparator (sorted edges, per-row ID sets).
- **FR-004**: System MUST pass all 439 existing tests on the quad-edge backend without test-code modifications.
- **FR-005**: Benchmark script MUST emit markdown table with four backends: EdgeMap, HE-v1, HE-v2, Quad-Edge. Four operations per backend. Median-of-3 trials per operation.
- **FR-006**: System MUST preserve mixed-element handling without spurious edges.
- **FR-007**: System MUST raise `ValueError` on unknown `CHILMESH_TOPOLOGY_BACKEND` value (including typos like "quadegg" misspelling).

### Non-Functional Requirements

- **NFR-001**: Quad-edge backend MUST initialize WNAT_Hagen in ≤ 3.6 s (match half-edge baseline, <10% degradation from v0.4.0 baseline of 3.26s).
- **NFR-002**: Quad-edge backend MUST NOT increase peak memory by more than 25% vs. EdgeMap baseline.
- **NFR-003**: Benchmark variance across three runs MUST be < 5% per operation.

### Key Entities

- **QuadEdge**: One directed edge as part of a 4-tuple (dual representation of edge topology). Knows: next clockwise edge, next counter-clockwise edge, opposite edge (twin), edge ID.
- **QuadEdgeTopology**: Container for all quad-edges + lookup indices. Equivalent in role to EdgeMap and HalfEdgeTopology; emits same adjacency outputs.
- **TopologyBackend (enum-like)**: Three values — `"edgemap"` (default), `"halfedge"` (Phase 007), `"quadegg"` (this feature).

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: 100% of existing tests (439 tests) pass with `CHILMESH_TOPOLOGY_BACKEND=quadegg` set, no test code modified.
- **SC-002**: WNAT_Hagen full init on quad-edge backend completes in ≤ 3.6 s (matches NFR-001).
- **SC-003**: Benchmark table shows quad-edge alongside EdgeMap, HE-v1, HE-v2; all four backends have complete metrics.
- **SC-004**: Adjacency-equivalence test (quad-edge vs. EdgeMap, fixture-by-fixture) shows zero mismatched entries.
- **SC-005**: `.planning/008-DECISION.md` cites benchmark numbers and recommends one of: "adopt quad-edge", "adopt half-edge", "stick with EdgeMap", or "archive topology work".

## Assumptions

- **A-001**: Reference hardware and Python version match Phase 007 (half-edge) benchmarks for fair comparison. Recorded in `output/benchmark.json`.
- **A-002**: All 439 tests pass on main branch baseline before this work begins.
- **A-003**: WNAT_Hagen mesh available via ADMESH-Domains; SHA-256 verified (same pin file as Phase 007).
- **A-004**: Pure-Python quad-edge is the prototype target. Rust/C++ porting is Phase 9+ deferred work.
- **A-005**: No new public API beyond `topology_backend='quadegg'` kwarg — backward compatible.
- **A-006**: Quad-edge construction algorithm achieves O(n) complexity like half-edge (not O(n²) like naive implementations).

---

## Out of Scope (Deferred Work)

- Rust/C++ quad-edge implementation (deferred to Phase 9+).
- Compiled optimizations (Numba, Cython, C extensions).
- Spatial indexing or mesh mutation.
- Public API additions beyond `topology_backend` kwarg.
- Skeletonization algorithm refactor (must remain bit-identical).
