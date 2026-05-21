# SPEC: Phase 008 — Quad-Edge Data Structure Investigation

**Phase:** 008-optimize-port-w-quad-edge  
**Title:** Quad-Edge Data Structure Investigation & Benchmarking  
**Created:** 2026-05-21  
**Ambiguity Score:** 0.170 (gate: ≤ 0.20) ✓  

---

## Goal

Implement a quad-edge (4-connected) topology backend for CHILmesh, starting from the Wikipedia quad-edge definition and adapting as needed for 2D mixed-element (triangle + quad) meshes. Benchmark quad-edge alongside existing backends (EdgeMap baseline, Half-Edge v1/v2) on the WNAT_Hagen reference mesh. Produce a comparison table and a data-driven decision record (DECISION.md) recommending which backend to pursue for language porting (Phase 9+).

## In Scope

1. **Quad-edge definition & adaptation** — Start with Wikipedia's quad-edge (4-tuple per edge; next-clockwise, next-counter-clockwise, opposite-edge pointers). Adapt for CHILmesh constraints (mixed-element triangles/quads, existing `_elem_type` mask, boundary sentinels). Document any departures from wiki definition with rationale.

2. **Quad-edge implementation** — Pure-Python ndarray-based storage (like half-edge); no compiled extensions. Required adjacency conversions (to_edge2vert, to_elem2edge, to_vert2edge, to_vert2elem, to_edge2elem) MUST produce bit-identical outputs vs. EdgeMap baseline.

3. **Backend integration** — Add `topology_backend='quadegg'` kwarg to CHILmesh constructor. Integrate into `_build_adjacencies()` dispatch (same pattern as half-edge). Preserve public API — no breaking changes.

4. **Equivalence validation** — All 439 existing tests MUST pass with `CHILMESH_TOPOLOGY_BACKEND=quadegg`; no test modifications. Bit-identical adjacency audit (canonical-form comparator) on all four fixtures.

5. **Benchmarking** — Run quad-edge alongside EdgeMap, Half-Edge-v1, Half-Edge-v2 on WNAT_Hagen (52.7k vertices). Median-of-3 trials per operation. Four operations: fast_init, full_init, quality_analysis, query_latency. Update `output/benchmark.json` and emit markdown comparison table.

6. **Decision record** — Write `008-DECISION.md` citing exact benchmark numbers. Recommendation: adopt quad-edge, adopt half-edge, stick with EdgeMap, or archive topology work. Decision is data-driven — no pre-set speed threshold.

## Out of Scope

- Rust/C++ quad-edge port (Phase 9+ deferred work)
- Compiled optimizations (Numba, Cython, C extensions)
- Spatial indexing, mesh mutation, public API beyond `topology_backend` kwarg
- Skeletonization algorithm refactor (must remain bit-identical)
- Changes to existing test files or test methodology

## Requirements

### Functional Requirements

| Req | Current State | Target State | Acceptance Criterion |
|-----|---|---|---|
| **F-001: Quad-edge definition** | No quad-edge in CHILmesh | Quad-edge documented, adapted for mixed-element meshes, stored as ndarray | `src/chilmesh/mesh_topology_quadegg.py` exists; docstring explains the 4-tuple structure and how it differs from wiki (if at all) |
| **F-002: Backend selection** | Only "edgemap" and "halfedge" supported | Support `topology_backend='quadegg'` kwarg + env var | `CHILmesh(topology_backend='quadegg')` succeeds; `CHILMESH_TOPOLOGY_BACKEND=quadegg pytest -v` passes |
| **F-003: Adjacency equivalence** | Half-edge shows equivalence gaps in Elem2Edge ordering | Quad-edge produces bit-identical adjacency outputs (sorted edges, canonical form) | `tests/test_quad_edge_equivalence.py` (8 test cases, one per adjacency type) all pass |
| **F-004: Test pass rate** | 439 tests pass on EdgeMap + HalfEdge backends | 439 tests pass on quad-edge backend without modification | `CHILMESH_TOPOLOGY_BACKEND=quadegg pytest -v` → all green |
| **F-005: Benchmark inclusion** | Benchmark shows EdgeMap, HE-v1, HE-v2 only | Benchmark shows all four backends: EdgeMap, HE-v1, HE-v2, Quad-Edge | `output/benchmark.json` contains a "Quad-Edge" entry with four operations; markdown table has four columns |
| **F-006: Decision record** | No DECISION.md | `008-DECISION.md` exists, cites benchmark data, recommends one clear path | File exists and states "Adopt quad-edge" OR "Adopt half-edge" OR "Stick with EdgeMap" OR "Archive topology work"; each statement backed by inline benchmark numbers |

### Non-Functional Requirements

| Req | Metric | Acceptance |
|-----|---|---|
| **N-001: Performance parity** | Quad-edge WNAT_Hagen full init time | ≤ 3.6s (same as half-edge NFR-001; <10% regression from v0.4.0 baseline of 3.26s) |
| **N-002: Memory overhead** | Quad-edge peak memory vs. EdgeMap | ≤ 25% increase (same as half-edge) |
| **N-003: Benchmark variance** | 3-trial std dev per operation | < 5% variance (rules out one-off noise) |

## Boundaries

### In Scope (Yes, this phase delivers)
✓ Quad-edge implementation, definition documentation, equivalence testing, benchmarking, decision record

### Out of Scope (No, deferred or adjacent)
✗ Rust/C++ porting (Phase 9+)
✗ Compiled optimizations  
✗ Public API beyond `topology_backend` kwarg  
✗ Skeletonization refactor  
✗ Spatial indexing  
✗ Test methodology changes  

**Boundary reasoning:** Quad-edge is a structural alternative; porting, compilation, and skeletonization are separate concerns. This phase is purely "Python prototype + evaluation". Public API stability is a production concern; `topology_backend` is optional and internal.

## Constraints

1. **Wiki definition as starting point** — Quad-edge MUST begin with Wikipedia's 4-tuple model; deviations from wiki are documented and justified.
2. **Mixed-element support** — Quad-edge MUST handle triangles (padded quads with `_elem_type` mask) and quads without storing spurious edges.
3. **Boundary sentinels** — Quad-edge boundary representation (undefined twin) MUST convert to `-1` in `Edge2Elem` output for downstream compatibility.
4. **Bit-identical adjacency** — Canonical-form comparator applies (edges sorted by (min_vert, max_vert); per-row ID sets); naive array equality does NOT apply.
5. **No breaking changes** — Public methods, attribute names, fort.14 I/O MUST remain unchanged.

## Acceptance Criteria (Pass/Fail Checkboxes)

- [ ] `src/chilmesh/mesh_topology_quadegg.py` implemented with clear construction algorithm docstring
- [ ] `topology_backend='quadegg'` kwarg works in CHILmesh constructor and via `CHILMESH_TOPOLOGY_BACKEND` env var
- [ ] `CHILMESH_TOPOLOGY_BACKEND=quadegg pytest -v` → 439/439 tests pass (no test modifications)
- [ ] Adjacency equivalence tests (quad-edge vs. EdgeMap) all pass on four fixtures
- [ ] `output/benchmark.json` contains "Quad-Edge" backend with four operations complete
- [ ] Benchmark markdown table shows four backends × four operations with percent deltas
- [ ] `008-DECISION.md` exists and recommends clear next step with benchmark data cited
- [ ] Quad-edge WNAT_Hagen full init ≤ 3.6s (N-001)
- [ ] Memory overhead ≤ 25% vs. EdgeMap (N-002)
- [ ] Benchmark variance < 5% (N-003)

## Ambiguity Report

**Final Ambiguity Score: 0.170** (gate: ≤ 0.20) ✓

| Dimension | Score | Minimum | Status |
|---|---|---|---|
| Goal Clarity | 0.85 | 0.75 | ✓ Exceeded |
| Boundary Clarity | 0.80 | 0.70 | ✓ Exceeded |
| Constraint Clarity | 0.78 | 0.65 | ✓ Exceeded |
| Acceptance Criteria | 0.85 | 0.70 | ✓ Exceeded |

**Interview findings:**
- Quad-edge definition is IN SCOPE (to be determined during phase, starting from wiki)
- Benchmark is NEUTRAL (no pre-set speed threshold; data drives decision)
- Adaptation to CHILmesh is PERMITTED (wiki + customization allowed)

---

## Success Metrics (How We Know It's Done)

1. **Code complete** — Quad-edge backend implemented, passes all tests, integrates into dispatch
2. **Benchmarked** — Four-backend comparison table produced; all operations measured
3. **Decided** — DECISION.md recommends clear adoption path (adopt QE, adopt HE, stick with EM, archive)
4. **Documented** — Quad-edge structure documented; wiki vs. adapted differences explained

---

**Next Step:** Run `/gsd:discuss-phase 008` to lock implementation approach (data structure layout, construction algorithm, conversion strategy).
