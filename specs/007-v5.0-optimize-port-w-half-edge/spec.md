# Feature Specification: Half-Edge Data Structure Investigation & Language Optimization Path

**Feature Branch**: `v5.0-optimize-port-w-half-edge`
**Created**: 2026-05-21
**Status**: Draft (post-specify, awaiting clarify)
**Input**: User description: "Investigate half-edge / doubly-connected edge list as optimization path; evaluate language porting (Rust/C++) as follow-on. Recompile benchmarks and update reader."

**Related Issues**: #137 (Half-edge investigation), #68 (Quad-edge alternative)

---

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Library maintainer evaluates half-edge for adoption (Priority: P1)

CHILmesh maintainer wants to know whether replacing the current EdgeMap + dict-based adjacency model with a half-edge (DCEL) structure will yield measurable performance improvement on the WNAT_Hagen reference mesh, without regressing correctness on any of the four built-in fixtures.

**Why this priority**: Phase 4 already achieved 937× speedup; further gains require structural change. Half-edge is the canonical alternative cited in CGAL/Triangle/Gmsh literature. Without P1 the maintainer cannot decide whether to invest further effort in DCEL or move to Rust/C++ porting (P3).

**Independent Test**: Run `pytest -v` on a half-edge backend and confirm 100% pass rate; run `python scripts/benchmark_halfedge_variants.py` and confirm output table with three columns (current, half-edge v1, optimized variant) and percent deltas for fast_init, full_init, quality_analysis, query latency. P1 ships as soon as those two artifacts exist.

**Acceptance Scenarios**:

1. **Given** the codebase on `v5.0-optimize-port-w-half-edge` branch, **When** `pytest -v` runs with `CHILMESH_TOPOLOGY_BACKEND=halfedge`, **Then** all existing tests pass with zero modifications.
2. **Given** a benchmark script runs the three variants on WNAT_Hagen (52.7k vertices), **When** results are collected, **Then** a markdown table is produced showing operation-by-operation deltas vs. the v0.4.1 baseline.
3. **Given** the half-edge backend is selected, **When** mesh topology is built on any of the four fixtures (annulus, donut, block_o, structured), **Then** adjacency outputs (`Edge2Vert`, `Elem2Edge`, `Vert2Edge`, `Vert2Elem`, `Edge2Elem`) are bit-identical to the EdgeMap baseline.

---

### User Story 2 — Downstream library consumer reads updated benchmark report (Priority: P2)

Downstream consumers (MADMESHR, ADMESH, ADMESH-Domains) need to see the current performance characteristics of CHILmesh on the reference mesh before deciding whether to upgrade their pinned version.

**Why this priority**: Documentation lag was a P2 surfaced in the Phase 4 lessons-learned (`docs/BENCHMARK.md` was stale). Consumers cannot pin to a CHILmesh version without trustworthy numbers. Lower priority than P1 because consumers don't see the half-edge prototype directly — they see the reader.

**Independent Test**: Open `docs/BENCHMARK.md` on the branch and confirm the v0.4.1 figures match what `python scripts/benchmark_wnat_hagen.py` reports today, and that any half-edge results are clearly tagged "prototype / not adopted yet" until the DECISION record is created.

**Acceptance Scenarios**:

1. **Given** `docs/BENCHMARK.md` on this branch, **When** a downstream consumer reads it, **Then** they see a clearly labeled "Current Performance (v0.4.1)" section with timings that match the live benchmark script output.
2. **Given** a half-edge prototype has been benchmarked, **When** a consumer reads the report, **Then** the half-edge column is labeled "experimental — not part of the public API" so they do not pin against it prematurely.

---

### User Story 3 — Library maintainer decides on language porting (Priority: P3)

If half-edge in pure Python shows >2× speedup, the maintainer wants justification to invest in a Rust/C++ port (via PyO3, pybind11, or cffi). If half-edge shows ≤baseline performance, the maintainer wants the half-edge branch archived with a recorded "not beneficial" decision so future contributors don't redo the same analysis.

**Why this priority**: This is the entire reason the feature is called "optimize-port-w-half-edge" — DCEL is a stepping stone, not the destination. P3 because it's a downstream gate: it only fires after P1's benchmark numbers are in.

**Independent Test**: Read `.planning/v5.0-DECISION.md` (created in a follow-on phase). It must contain either "Recommend Rust/C++ port in Phase 6 because [data]" or "Archive half-edge; pure-Python DCEL provides no benefit because [data]" — with the actual measured percent deltas inline.

**Acceptance Scenarios**:

1. **Given** the benchmark table from P1, **When** the maintainer reads it, **Then** the speedup column unambiguously triggers either "pursue port" or "archive" — no judgment call needed beyond reading the number.
2. **Given** "pursue port" is triggered, **When** Phase 6 opens, **Then** the DECISION record specifies which language (Rust preferred for memory safety; C++ as fallback) and which binding library (PyO3 for Rust, pybind11 for C++).

---

### Edge Cases

- **Mixed-element padded triangles**: The current adjacency model uses a padded-quad convention (vertex 0 in slot 3 marks a triangle in a 4-column array). Half-edge must preserve this without storing fake edges for the padding slot, otherwise `_elem_type` regresses (B3 from v0.1.1).
- **Boundary edges**: `Edge2Elem` uses `-1` as the sentinel for the "outside" element on a boundary edge. Half-edge representations natively use `None` / `null` for the twin half-edge of a boundary; the converter MUST emit `-1` to preserve downstream contracts.
- **Disconnected meshes**: ADMESH-Domains catalog includes meshes with multiple connected components. Half-edge construction MUST handle each component independently — no global assumption that every vertex has degree ≥ 1.
- **Degenerate elements**: Zero-area triangles / zero-length edges (e.g., the quad angle test from B5) must not crash half-edge construction — they must produce the same warnings/values the current code does.
- **Environment variable mis-set**: If `CHILMESH_TOPOLOGY_BACKEND` is set to an unknown value (e.g., typo "halfedage"), the system MUST raise a clear `ValueError` at mesh construction — not silently fall back to EdgeMap and confuse the benchmark.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST provide an internal half-edge (DCEL) topology backend selectable via either a constructor argument (`topology_backend='halfedge'`) or the environment variable `CHILMESH_TOPOLOGY_BACKEND`.
- **FR-002**: System MUST keep the public CHILmesh API surface (constructor signature minus the new optional kwarg, public methods, public attributes) unchanged. No downstream import paths break.
- **FR-003**: System MUST produce bit-identical adjacency outputs (`Edge2Vert`, `Elem2Edge`, `Vert2Edge`, `Vert2Elem`, `Edge2Elem`, `EdgeMap.to_list()`) between the EdgeMap backend and the half-edge backend on every built-in fixture.
- **FR-004**: System MUST pass all 439 existing tests on both backends without test-code modifications.
- **FR-005**: System MUST recompile `docs/BENCHMARK.md` figures from a live benchmark run (no estimated numbers); the document MUST link the JSON output (`output/benchmark.json`) as the source of truth.
- **FR-006**: System MUST emit a markdown comparison table (current vs. half-edge v1 vs. optimized variant) for at least four operations: fast init, full init, quality analysis, query latency.
- **FR-007**: System MUST preserve mixed-element (triangle + padded-quad) handling without storing dummy half-edges for the padding slot.
- **FR-008**: System MUST raise `ValueError` at mesh construction if `CHILMESH_TOPOLOGY_BACKEND` holds an unknown value.

### Non-Functional Requirements

- **NFR-001**: Half-edge backend MUST initialize WNAT_Hagen in ≤ 3.6 s on the reference hardware (≤ 10% degradation from the v0.4.0 baseline of 3.26 s).
- **NFR-002**: Half-edge backend MUST NOT increase peak memory by more than 25% vs. the EdgeMap backend, measured on WNAT_Hagen.
- **NFR-003**: Benchmark variance across three runs MUST be < 5% per operation (rules out one-off noise driving the adoption decision).

### Key Entities

- **HalfEdge**: One directed half of an undirected edge. Knows: origin vertex, twin half-edge, next half-edge around its face, incident face. In CHILmesh's 2D mixed-element setting, a "face" is an element (triangle or quad).
- **HalfEdgeTopology**: Container for all HalfEdges + lookup indices. Equivalent in role to EdgeMap; emits the same dict/ndarray adjacency outputs.
- **TopologyBackend (enum-like)**: Two values — `"edgemap"` (current default) and `"halfedge"` (this feature). Used by the factory inside `_build_adjacencies()`.
- **BenchmarkVariant**: Tuple of (backend_name, mesh_path, operation_label, mean_seconds, std_seconds). Three variants per operation drive the table in `docs/BENCHMARK.md`.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: 100% of existing tests (439 tests) pass with `CHILMESH_TOPOLOGY_BACKEND=halfedge` set, no test code modified.
- **SC-002**: WNAT_Hagen full init (with layers) on half-edge backend completes in ≤ 3.6 s on reference hardware (matches NFR-001).
- **SC-003**: `docs/BENCHMARK.md` is regenerated from live benchmark output; the JSON file referenced from the doc parses cleanly and contains all four operation columns.
- **SC-004**: Adjacency-equivalence test (half-edge vs. EdgeMap, fixture-by-fixture) shows zero mismatched entries.
- **SC-005**: The decision-trigger column in the benchmark table unambiguously points to one of three outcomes: "adopt" / "archive" / "investigate port" — no further analysis required to choose.

## Assumptions

- **A-001**: Reference hardware for benchmark numbers is whichever machine ran the v0.4.1 baseline (current container counts; we record `platform.platform()` + `platform.python_version()` in the JSON output).
- **A-002**: All 439 existing tests are currently passing on the `main` branch baseline before this work begins (verified pre-flight).
- **A-003**: The WNAT_Hagen mesh remains available via ADMESH-Domains; if the upstream registry moves the file, the benchmark script's fallback paths catch it.
- **A-004**: Pure-Python half-edge is the prototype target. If P3's "investigate port" trigger fires, that's a separate phase — this spec does not commit to delivering Rust/C++ code.
- **A-005**: No new public API enters the package. The `topology_backend` constructor argument is the only addition, and it's a kwarg with a default — backward compatible.
- **A-006**: Phase 4 lessons-learned conventions hold: one logical change per commit, regression tests for every fix, and `chilmesh.examples` is the canonical entry point for fixtures.

---

## Out of Scope (Deferred Work)

The following are explicitly NOT delivered by this spec:

- Full Rust/C++ implementation (deferred to a Phase 6 if P3 trigger fires).
- New language binding infrastructure (PyO3, pybind11, cffi) — also Phase 6.
- Spatial indexing improvements (#94) — independent Phase 5+ work.
- Incremental / dynamic mesh mutation (#93) — independent Phase 5+ work.
- Public API additions beyond `topology_backend` kwarg.
- Skeletonization algorithm refactor (behavior must remain bit-identical).
