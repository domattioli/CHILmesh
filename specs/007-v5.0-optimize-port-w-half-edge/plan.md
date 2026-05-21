# Implementation Plan: Half-Edge Data Structure Investigation & Language Optimization Path

**Branch**: `v5.0-optimize-port-w-half-edge` | **Date**: 2026-05-21 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `specs/007-v5.0-optimize-port-w-half-edge/spec.md` + clarifications from `clarify.md`

## Summary

Deliver an internal DCEL (doubly-connected edge list) topology backend for CHILmesh as a feasibility prototype. Run a 3-variant benchmark on WNAT_Hagen (current EdgeMap, half-edge v1, half-edge v2 with NumPy vectorized face walk). Recompile `docs/BENCHMARK.md` from live data. Defer the adopt/archive/port decision to a v5.1 follow-on phase.

The half-edge backend is selected via either `topology_backend='halfedge'` kwarg or `CHILMESH_TOPOLOGY_BACKEND` env var. Public API surface stays unchanged. All 439 existing tests run unmodified on both backends.

## Technical Context

**Language/Version**: Python 3.10+ (pure Python; no compiled extensions in this phase per clarify Q-1 and spec out-of-scope)
**Primary Dependencies**: numpy ≥1.23 (DCEL storage as `ndarray[n_halfedges, 4]` per clarify C-1), scipy ≥1.10 (existing), matplotlib ≥3.6 (existing, visualization only)
**Storage**: N/A — mesh data only; benchmark JSON written to `output/benchmark.json`
**Testing**: pytest ≥7 with parametrized fixtures (annulus, donut, block_o, structured); session-level memoization stays
**Target Platform**: Linux/macOS/Windows (CI matrix unchanged)
**Project Type**: Library (Python package, src/ layout)
**Performance Goals**:
- Half-edge full init on WNAT_Hagen ≤ 3.6 s (NFR-001; 10% degradation budget over v0.4.0's 3.26 s)
- Memory overhead ≤ 25% vs. EdgeMap measured via tracemalloc (NFR-002, F-3)
- Variance across 3 runs < 5% per operation (NFR-003)
**Constraints**:
- Public API surface frozen — only one new kwarg permitted (`topology_backend`)
- All 439 existing tests pass on both backends without test-code modification (FR-004)
- Bit-identical adjacency outputs under canonical-form comparator (FR-003, C-4)
- Mixed-element (padded triangle) handling preserved (FR-007, F-4)
**Scale/Scope**: ~700 LOC new code (300 DCEL module + 200 tests + 200 benchmark script + factory wiring); 4 fixtures × 2 backends = 8 test runs in CI

## Constitution Check

*Gate: must pass before Phase 0 research; re-check after Phase 1 design.*

Against `.specify/memory/constitution.md` v1.1 — CHILmesh project constitution:

| Principle | Check | Status |
|---|---|---|
| I. Library-First | Half-edge is an internal module; no new public exports beyond one kwarg. Importable as `from chilmesh.mesh_topology_halfedge import HalfEdgeTopology` (private use). | ✅ Pass |
| II. Mesh Immutability | Half-edge built once during `_build_adjacencies()`; no incremental updates in this phase (spec out-of-scope; clarify Q-3). | ✅ Pass |
| III. Test-First | All 5 tasks have regression tests; equivalence test in Task 2 written BEFORE Task 3 (factory integration). | ✅ Pass |
| IV. Geometric Correctness Over Performance | Half-edge construction runs `_ensure_ccw_orientation` first (F-1); bit-identical adjacency is the gate (FR-003), not raw speed. | ✅ Pass |
| V. Coordinate-System Agnostic | No new coordinate assumptions; half-edge operates on topology, not geometry. | ✅ Pass |
| VI. Format Pluralism (ADCIRC + 2DM) | Both readers route through `_build_adjacencies()`; both backends auto-supported. | ✅ Pass |
| VII. Public API Stability | One new optional kwarg (`topology_backend`); default value preserves existing behavior; no breaking changes. | ✅ Pass |
| VIII. Documentation = Contract | `mesh_topology_halfedge.py` ships with docstrings; `docs/BENCHMARK.md` updated with live numbers (FR-005). | ✅ Pass |
| IX. Mixed-Element Matrices + Correctness | Padded-triangle handling preserved via F-4 mitigation (reuse `_elem_type` mask). | ✅ Pass |
| X. MATLAB Parity Where It Matters | Skeletonization output unchanged (spec out-of-scope: "Skeletonization algorithm refactor"). | ✅ Pass |

**No violations. No entries needed in Complexity Tracking section.**

## Project Structure

### Documentation (this feature)

```text
specs/007-v5.0-optimize-port-w-half-edge/
├── spec.md              # Feature spec (3 user stories, 8 FRs, 3 NFRs, 5 SCs)
├── clarify.md           # 3-round Spec-Kit clarify (C-1..C-4, F-1..F-4, Q-1..Q-4)
├── plan.md              # This file
├── research.md          # Phase 0 output — DCEL design analysis (created next)
├── data-model.md        # Phase 1 output — HalfEdgeTopology schema, canonical comparator
├── quickstart.md        # Phase 1 output — how to flip backend, how to run benchmark
└── tasks.md             # Phase 2 output (/speckit-tasks command)
```

### Source Code (repository root)

```text
src/chilmesh/
├── CHILmesh.py                  # MODIFIED: _build_adjacencies() gains backend factory dispatch
├── mesh_topology.py             # UNCHANGED: existing EdgeMap (current default backend)
├── mesh_topology_halfedge.py    # NEW: HalfEdgeTopology, build_halfedge_from_connectivity()
├── __init__.py                  # MODIFIED: re-exports nothing new (kwarg only; private module)
└── examples/                    # UNCHANGED: factories accept topology_backend via **kwargs passthrough

tests/
├── conftest.py                  # MODIFIED: add topology_backend fixture with env-var monkeypatch (F-2)
├── test_halfedge_basic.py       # NEW: ~15 tests — construction on 4 fixtures, twin pairing, boundary
├── test_halfedge_equivalence.py # NEW: canonical-form comparator vs. EdgeMap on all 4 fixtures (C-4)
└── [existing 26 test files]     # UNCHANGED: run on both backends via env var

scripts/
├── benchmark_wnat_hagen.py      # UNCHANGED: continues to produce v0.4.1-shape JSON
├── benchmark_halfedge_variants.py # NEW: runs current + v1 + v2 variants; writes 3-column markdown
└── verify_mesh_pin.py           # NEW: SHA-256 check for WNAT_Hagen (C-3)

docs/
└── BENCHMARK.md                 # MODIFIED: live regen from output/benchmark.json (FR-005)

output/
├── benchmark.json               # MODIFIED: gains backend_variants[], mesh_sha256, tracemalloc_peak (Q-2)
└── benchmark.log                # MODIFIED: textual transcript of last benchmark run
```

**Structure Decision**: Single-project library layout (Option 1 from template). One new internal module (`mesh_topology_halfedge.py`), one factory edit in `CHILmesh.py`, two new test files, two new scripts. No new top-level directories; no submodules; no packaging changes.

## Phase 0 — Research (next deliverable: `research.md`)

Research deliverables for `research.md`:

1. **DCEL storage layout** — confirm `ndarray[n_halfedges, 4]` (origin/twin/next/face) is the cache-friendliest pure-Python option. Benchmark a per-object alternative (HalfEdge dataclass) on a 1k-edge synthetic mesh; document the 3–10× slowdown that justifies C-1's decision.
2. **Canonical-form comparator algorithm** — write the sort-key for `EdgeMap.to_list()` output and demonstrate it's stable across both backends. Reference C-4 explicitly.
3. **Tracemalloc methodology** — confirm `tracemalloc.start()` / `take_snapshot()` / filter to `chilmesh.*` modules gives reproducible numbers per F-3; document the median-of-3 protocol.
4. **NumPy vectorized face walk** — prototype `np.take(half_edges[:,2], face_starts)` against the Python `for`-loop on a 10k-element synthetic mesh; document expected speedup range for the v2 "optimized variant" per C-2.
5. **Prior-art survey (1 paragraph each)** — CGAL HalfedgeDS, Triangle's edge representation, Gmsh's hybrid model, OpenMesh's property-based DCEL. Document what we're NOT borrowing and why (out-of-scope: properties system, lazy iterators, geometric predicates).

## Phase 1 — Design (deliverables: `data-model.md`, `quickstart.md`)

`data-model.md` documents:
- HalfEdgeTopology schema (the 4-column ndarray)
- TopologyBackend enum (`"edgemap" | "halfedge"`)
- Canonical-form comparator (sort keys, frozenset construction)
- BenchmarkVariant tuple structure → JSON schema (Q-2 additive fields)
- Boundary half-edge sentinel convention (`twin_idx == -1`, maps to `Edge2Elem == -1`)

`quickstart.md` documents:
- One-liner backend switch (env var + kwarg)
- How to run the 3-variant benchmark
- How to interpret the markdown table (decision-trigger column rule)
- Where to find the WNAT_Hagen SHA-256 pin

## Phase 2 — Implementation Tasks (deliverable: `tasks.md` via `/speckit-tasks`)

5 tasks. Sequential dependency chain (no parallelism viable — each task gates the next).

1. **T1**: Design + implement `mesh_topology_halfedge.py` (DCEL with ndarray-of-pointers per C-1; padded-triangle handling per F-4; CCW orientation pre-pass per F-1). ~300–350 LOC.
2. **T2**: Write `test_halfedge_basic.py` + `test_halfedge_equivalence.py` covering construction on all 4 fixtures, twin pairing invariants, boundary edge sentinels, and the canonical-form comparator vs. EdgeMap. ~200 LOC, ~15 tests + 4 fixture-parametrized equivalence tests.
3. **T3**: Factory integration in `CHILmesh._build_adjacencies()` — accept `topology_backend` kwarg threaded from `__init__` / `read_from_fort14` / `read_from_2dm` / `from_admesh_domain`. Add env-var fallback. Wire into `examples/*` factories via **kwargs passthrough. ValueError on unknown backend value (FR-008). ~50–80 LOC modified.
4. **T4**: Run full `pytest -v` suite on both backends. Baseline (default EdgeMap), variant (`CHILMESH_TOPOLOGY_BACKEND=halfedge pytest -v`). Capture both logs in `output/` for the PR record. Zero failures gate.
5. **T5**: `benchmark_halfedge_variants.py` extends `benchmark_wnat_hagen.py` with backend cycling + 3-run median + tracemalloc memory probe. Writes additive JSON fields (Q-2) + new markdown comparison table inserted into `docs/BENCHMARK.md`. Verifies WNAT_Hagen SHA-256 pin (C-3) — fails loudly on mismatch.

Each task ends with a verification step that re-runs the prior tasks' tests to confirm no regression.

## Phase 3 — Analysis & Audit (post-tasks; deliverable: `analyze.md`)

After T5 completes, `/speckit-analyze` runs against the actually-shipped artifacts:
- Cross-check every FR in `spec.md` against the code that satisfies it
- Cross-check every NFR against the measured numbers in `output/benchmark.json`
- Cross-check every clarification (C-1..C-4, F-1..F-4, Q-1..Q-4) against its mitigation site in source
- Produce one of three verdicts: PASS / PASS-with-followups / BLOCK
- BLOCK triggers a /speckit-clarify re-run; PASS-with-followups files individual issues; PASS hands off to `v5.1-decision` follow-on

## Complexity Tracking

> Only fill if Constitution Check has violations that must be justified.

**No violations.** All 10 constitution principles pass without exceptions.

## Risk Register

| Risk | Likelihood | Impact | Mitigation | Owner |
|---|---|---|---|---|
| Half-edge slower than EdgeMap in pure Python | Medium | Low (negative result still valuable per Q-4) | Document as experimental; archive backend; the negative result IS the deliverable | T5 author |
| Tracemalloc noise exceeds 5% (NFR-003 fails) | Low | Medium (NFR not provable) | Median-of-3 protocol per F-3; document the methodology in `research.md` so the choice is defensible | T5 author |
| ADMESH-Domains mesh file changes mid-phase | Low | High (benchmark results not comparable) | SHA-256 pin per C-3; benchmark script fails loudly on mismatch; pin checked into JSON | T5 author |
| Env-var leakage between tests | Medium | High (false-green tests) | monkeypatch fixture per F-2; documented in test file header | T2 author |
| Padded triangle gets a spurious 4th half-edge | Medium | High (regresses B3 from v0.1.1) | F-4 mitigation: reuse `_elem_type` mask before iterating element vertices | T1 author |
| Equivalence comparator passes on row-permuted but semantically wrong data | Low | High (false-green correctness) | Canonical-form sort BEFORE comparison; test the comparator itself with synthetic permutations | T2 author |
