# CONTEXT: Phase 008 — Quad-Edge Implementation Decisions

**Phase:** 008-optimize-port-w-quad-edge  
**Date:** 2026-05-21  
**Status:** Ready for research and planning  

---

## Spec Lock

**SPEC.md locked 6 functional requirements:**
1. Quad-edge definition & adaptation (start wiki, adapt for mixed-element CHILmesh)
2. Backend integration (`topology_backend='quadegg'` kwarg + env var)
3. Adjacency equivalence (bit-identical to EdgeMap on all fixtures)
4. Test pass rate (439/439 without modifications)
5. Benchmark inclusion (four backends: EdgeMap, HE-v1, HE-v2, Quad-Edge)
6. Decision record (DECISION.md with clear recommendation)

SPEC.md also locked NFR-001 (≤ 3.6s WNAT_Hagen full init), NFR-002 (≤ 25% memory overhead), NFR-003 (< 5% benchmark variance).

---

## Implementation Decisions

### Data Structure

**Decision:** Wiki semantics + CHILmesh naming convention

- Follow Wikipedia quad-edge 4-tuple structure (origin, next-clockwise, next-counter-clockwise, opposite-edge)
- Adopt CHILmesh field names for code clarity (e.g., `origin`, `twin_idx`, `next_cw`, `next_ccw` or similar)
- Store as `ndarray[n_edges, 4]` with integer-index pointers throughout
- Directionality (undirected vs. directed edges) TBD after algorithm design — let construction algorithm determine whether edges are (1-per-pair undirected) or (2-per-pair directed like half-edge)

**Rationale:** Wiki structure is proven; CHILmesh naming keeps codebase consistent and code maintainable. Directionality deferred because it's a consequence of the construction strategy, not an independent choice.

### Construction Algorithm

**Decision:** Invent new algorithm optimized for quad-edge (not reusing half-edge 3-phase)

- Half-edge Phase 007 proved 3-phase O(n) works, but is it optimal for quad-edge?
- Quad-edge's 4-pointer structure may allow direct (next_cw, next_ccw) assignment in a single pass
- Research phase MUST investigate whether a specialized quad-edge algorithm beats the 3-phase adaptation
- Final implementation complexity MUST be O(n) per SPEC.md A-006 (half-edge reference)

**Rationale:** Quad-edge is fundamentally different from half-edge (undirected vs. semi-directed). Specialized algorithm may yield better cache locality or fewer pointer chases.

### Boundary Representation

**Decision:** Start wiki-style, adapt for hydro domains' multiple boundaries

- Wiki quad-edge naturally leaves boundary edges' neighbors undefined
- Initial implementation: use -1 sentinel like half-edge for simplicity (convert to downstream -1 sentinel)
- **Constraint:** CHILmesh is used in hydro simulations with multiple disconnected boundary components; be prepared to extend boundary encoding beyond single -1 sentinel if ADMESH-Domains mesh catalog requires it
- If multiple boundaries emerge as a blocker, DECISION.md will note this as a design constraint for Phase 9+ (language porting)

**Rationale:** Keeps implementation simple initially; acknowledges real-world use case without over-designing.

### Adjacency Conversion Strategy

**Decision:** Write specialized native quad-edge → adjacency converter (not reusing half-edge converters)

- Do NOT convert quad-edge → half-edge → adjacency (would double overhead)
- Implement direct quad-edge traversal to emit Edge2Vert, Elem2Edge, Vert2Edge, Vert2Elem, Edge2Elem
- Converter MUST understand quad-edge semantics natively (not treat it as a half-edge variant)
- Canonical-form comparator applies: edges sorted by (min_vert, max_vert); per-row ID sets compared

**Rationale:** Reusing half-edge converters defeats the purpose of benchmarking quad-edge separately. Native converter also isolates quad-edge performance bottlenecks.

### Benchmark Scope

**Decision:** Quad-Edge v1 only (no v2 vectorized variant, TBD for v1.1)

- Benchmark quad-edge v1 (base Python implementation) against EdgeMap, HE-v1, HE-v2
- If quad-edge v1 shows promise, v1.1 can add NumPy vectorized walk (like HE-v2)
- This phase focuses on algorithm + adjacency strategy; optimization deferred

**Rationale:** Don't over-optimize prematurely. If quad-edge is slower than both half-edge variants, vectorization is wasted effort. If it's promising, v1.1 can invest in optimization.

---

## Reusable Patterns & Code Context

**Half-Edge Reference Implementation (Phase 007):**
- `src/chilmesh/mesh_topology_halfedge.py` — 3-phase construction (create, pair twins, assign next)
- `src/chilmesh/CHILmesh.py` — backend dispatch pattern in `_build_adjacencies()`
- `scripts/benchmark_halfedge_variants.py` — benchmark scaffold (4 operations, median-of-3, WNAT_Hagen)
- `tests/test_halfedge_equivalence.py` — equivalence test pattern (canonical-form comparator)
- Result: HE 80% slower than EdgeMap; archived as experimental, off by default

**Test Fixtures:**
- `tests/conftest.py` — parametrized over [annulus, donut, block_o, structured]
- All 439 existing tests reuse these; expect to pass without modification on quad-edge backend

**Adjacency Output Contract:**
- `Edge2Vert`, `Elem2Edge`, `Vert2Edge`, `Vert2Elem`, `Edge2Elem` outputs MUST match EdgeMap bit-for-bit (canonical form)
- Precedent: half-edge failed this due to edge ID ordering; quad-edge MUST nail it on first try

---

## Deferred Ideas (Not This Phase)

- Quad-edge v2 vectorization → Phase 008-v1.1 (only if v1 shows promise)
- Language porting (Rust/C++ via PyO3/pybind11) → Phase 9+ (only if quad-edge wins benchmark)
- Spatial indexing or mesh mutation → independent phases
- Skeletonization refactor → independent phase (must remain bit-identical)

---

## Canonical Refs

Downstream agents MUST read these files:

- **specs/008-optimize-port-w-quad-edge/008-SPEC.md** — locked requirements (6 FR + 3 NFR + 9 AC)
- **specs/007-optimize-port-w-half-edge/007-SPEC.md** — half-edge reference (to understand what NOT to replicate)
- **specs/007-optimize-port-w-half-edge/007-ANALYZE.md** — why half-edge was 80% slower (inform quad-edge algorithm design)
- **src/chilmesh/mesh_topology_halfedge.py** — working proof of adjacency conversion + 3-phase construction
- **src/chilmesh/CHILmesh.py** — backend dispatch integration pattern
- **scripts/benchmark_halfedge_variants.py** — benchmark scaffold to extend for quad-edge
- **tests/conftest.py** — test fixture parametrization (use as-is for quad-edge)
- **tests/test_halfedge_equivalence.py** — canonical-form equivalence test pattern

---

## Next Step

`/gsd:plan-phase 008` — planner will use these decisions to create a detailed implementation plan (task breakdown, dependencies, verification strategy).

Plan will likely include:
1. **Research:** Investigate optimal quad-edge algorithm for 2D meshes (literature review, compare to 3-phase adaptation)
2. **Implementation:** Build quad-edge backend step-by-step (data structure, construction, conversion)
3. **Testing:** Equivalence + test pass-rate validation
4. **Benchmarking:** Run quad-edge v1 against three other backends on WNAT_Hagen
5. **Decision:** Analyze results, write DECISION.md (adopt QE / adopt HE / stick with EM / archive)
