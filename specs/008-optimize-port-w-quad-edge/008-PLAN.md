# Implementation Plan: Quad-Edge Data Structure Investigation & Benchmarking

**Phase:** 008-optimize-port-w-quad-edge  
**Date:** 2026-05-21  
**Spec:** [008-SPEC.md](008-SPEC.md)  
**Context:** [008-CONTEXT.md](008-CONTEXT.md)  
**Reference:** Phase 007 half-edge (specs/007-optimize-port-w-half-edge/plan.md)

---

## Summary

Implement a quad-edge (4-connected) topology backend for CHILmesh following the Wikipedia definition with CHILmesh-specific naming. Benchmark quad-edge v1 (base Python) alongside EdgeMap, Half-Edge v1, and Half-Edge v2 on WNAT_Hagen reference mesh. Produce a four-backend comparison table and a data-driven decision record (DECISION.md) recommending adoption path for language porting (Phase 9+).

All 439 existing tests MUST pass unmodified on quad-edge backend. Adjacency outputs must be bit-identical to EdgeMap via canonical-form comparator. Two new test files validate construction and equivalence; existing tests parametrized over backends via env var.

---

## Phase Goal

**Outcome:** A working quad-edge backend proves whether 4-connected edge topology yields measurable performance advantages over half-edge (directed 2-per-pair) and EdgeMap (hash-based) baselines. Data-driven decision: adopt quad-edge, adopt half-edge, stick with EdgeMap, or archive topology work.

**Why This Matters:**
- Half-edge (Phase 007) showed 80% slowdown vs. EdgeMap; is quad-edge fundamentally better?
- Quad-edge's undirected 4-pointer structure may reduce cache misses and pointer chases
- Decision blocks Phase 9+ language porting effort (Rust/C++ target)

**Output:**
- `src/chilmesh/mesh_topology_quadegg.py` — Quad-edge implementation + conversion methods
- `tests/test_quadegg_construction.py` — Basic construction tests (15 cases × 4 fixtures)
- `tests/test_quadegg_equivalence.py` — Bit-identical adjacency validation vs. EdgeMap
- `scripts/benchmark_quadegg_variants.py` — Four-backend benchmark on WNAT_Hagen
- `output/benchmark.json` — Updated with quad-edge measurements
- `.planning/008-DECISION.md` — Recommendation with benchmark data cited

---

## Technical Context

**Language/Version:** Python 3.10+ (pure Python; no compiled extensions per CONTEXT.md)  
**Primary Dependencies:** numpy ≥1.23 (quad-edge storage as ndarray[n_edges, 4])  
**Storage:** In-memory only; benchmark JSON written to `output/benchmark.json`  
**Testing:** pytest ≥7 with parametrized fixtures (annulus, donut, block_o, structured)  
**Target Platform:** Linux/macOS/Windows (CI matrix unchanged)

**Performance Goals (from SPEC.md):**
- Quad-edge full init on WNAT_Hagen ≤ 3.6s (NFR-001; 10% degradation budget)
- Memory overhead ≤ 25% vs. EdgeMap (NFR-002)
- Benchmark variance < 5% across 3 trials (NFR-003)

**Constraints:**
- Quad-edge directed vs. undirected: construction algorithm determines (CONTEXT.md locked decision)
- Boundary handling: start wiki-style (-1 sentinel), extend if hydro domains require multiple boundaries
- Adjacency conversion: native quad-edge traversal (NOT quad-edge → half-edge → adjacency)
- All 439 existing tests pass without modification (FR-004)
- Bit-identical outputs under canonical-form comparator (FR-003)
- Mixed-element support preserved via `_elem_type` mask (CONTEXT.md)

---

## Reusable Patterns from Phase 007

| Pattern | Location | Reuse for Phase 008 |
|---------|----------|-------------------|
| Backend dispatch in `_build_adjacencies()` | `src/chilmesh/CHILmesh.py:419` | Add `'quadegg'` case to factory switch |
| Equivalence test pattern | `tests/test_halfedge_equivalence.py` | Canonical-form comparator (same, edge sorting) |
| Canonical-form comparator | `HalfEdgeTopology.to_edgemap_list()` | Reuse for quad-edge Edge2Vert extraction |
| Benchmark scaffold | `scripts/benchmark_halfedge_variants.py` | Extend to 4 backends; same operations |
| Test fixtures | `tests/conftest.py` | Parametrization unchanged; quad-edge backend added |
| Construction reference | `mesh_topology_halfedge.py` | Study 3-phase O(n) algorithm as comparison baseline |

---

## Dependency Graph & Wave Structure

### Task Dependencies

```
1. Research Algorithm
   ↓
2. Quad-edge Implementation
   ├→ 3. Construction Tests
   ├→ 4. Backend Integration
   │
   ├→ 5. Equivalence Tests
   │
   └→ 6. Run Full Test Suite
      ↓
7. Benchmark Script + Run
   ↓
8. Decision Record
```

**Wave 1 (Sequential):** Tasks 1–2 (research + core implementation)  
**Wave 2 (Parallel):** Tasks 3–4 (tests + integration can start once impl draft exists)  
**Wave 3 (Sequential):** Tasks 5–6 (equivalence + full suite; depends on impl stability)  
**Wave 4 (Sequential):** Tasks 7–8 (benchmark + decision; final deliverables)

**Rationale:** Phase 007 half-edge showed that algorithms are tightly coupled to construction invariants. Research must precede coding to avoid rework. Once construction stabilizes, tests and integration can proceed in parallel (different files). Equivalence testing depends on implementation maturity. Benchmark is final.

---

## Must-Haves (Goal-Backward)

### Observable Truths
- Quad-edge backend can be selected via `topology_backend='quadegg'` kwarg
- All 439 existing tests pass with `CHILMESH_TOPOLOGY_BACKEND=quadegg`
- Quad-edge produces bit-identical adjacency outputs vs. EdgeMap (canonical form)
- Benchmark table shows quad-edge alongside three other backends
- Performance analysis supports a clear adoption decision

### Required Artifacts
| Artifact | Purpose | Acceptance |
|----------|---------|-----------|
| `mesh_topology_quadegg.py` | Quad-edge data structure + adjacency converters | File exists; implements construction algorithm with O(n) complexity; docstring explains 4-tuple structure |
| `test_quadegg_construction.py` | Basic construction + boundary handling on all fixtures | 15+ test cases; parametrized over 4 fixtures; validates `half_edges` array invariants |
| `test_quadegg_equivalence.py` | Canonical-form comparison vs. EdgeMap | 8 test cases (one per adjacency type × 2 fixtures); zero mismatches |
| `benchmark_quadegg_variants.py` | Four-backend benchmark on WNAT_Hagen | Script measures EdgeMap, HE-v1, HE-v2, Quad-Edge; median-of-3 per operation; writes JSON |
| `CHILmesh._build_adjacencies()` | Backend dispatch updated | Factory switch includes `'quadegg'` case; env var fallback works |
| `output/benchmark.json` | Benchmark data with quad-edge | Four backends × four operations complete; variance < 5% |
| `.planning/008-DECISION.md` | Recommendation with data | States one of: adopt QE / adopt HE / stick with EM / archive; cites exact numbers |

### Key Links (Critical Connections)
| From | To | Via | What Breaks If Missing |
|------|----|----|------------------------|
| `test_quadegg_*.py` | `mesh_topology_quadegg.py` | import `QuadEdgeTopology` | Tests cannot run; construction invariants unknown |
| `CHILmesh._build_adjacencies()` | `mesh_topology_quadegg.py` | import `build_quadegg_from_connectivity` | Backend dispatch fails; cannot construct |
| `benchmark_quadegg_variants.py` | `CHILmesh` + `mesh_topology_quadegg.py` | `topology_backend='quadegg'` kwarg | Benchmark missing quad-edge; decision data incomplete |
| `test_quadegg_equivalence.py` | `EdgeMap` (Phase 1) | Canonical-form comparator | Cannot validate bit-identity; regressions not caught |

---

## Plan Tasks

### Task 1: Research Optimal Quad-Edge Algorithm for 2D Meshes

**Wave:** 1  
**Type:** Investigation + Documentation  
**Files:** None created (research output is design input for Task 2)  
**Duration:** ~2–3 hours of focused research

**Action:**

Investigate and document the optimal quad-edge construction algorithm for 2D mixed-element meshes. Start with Wikipedia's quad-edge definition (4-tuple per edge: origin, next-clockwise, next-counter-clockwise, opposite-edge) and adapt for CHILmesh constraints.

**Research questions to answer:**

1. **Directionality choice:** Wikipedia quad-edge is typically undirected (1 edge per pair, 4 neighbors). Half-edge is directed (2 per pair, 1 next neighbor per direction). For CHILmesh's 2D triangle+quad meshes, should quad-edge be:
   - Undirected (1 per pair, reduce storage × 2)? Requires careful boundary handling.
   - Directed (2 per pair like half-edge)? Simplifies logic; how much speedup lost?
   - Document the choice with rationale; let Phase 2 implementation decide based on your recommendation.

2. **Algorithm design:** Compare two approaches:
   - **3-phase adaptation (copy from half-edge):** Phase 1 create, Phase 2 pair twins, Phase 3 assign next pointers. Worst case O(n log n) if twin pairing uses sorting.
   - **Single-pass optimized algorithm:** For quad-edge's 4-pointer structure, can origin + next_cw + next_ccw be assigned in one walk? Document expected complexity (O(n) or O(n log n)) and why.
   - Benchmark synthetic mesh (1k–5k edges) if feasible; document results so Phase 2 can choose confidently.

3. **Boundary semantics:** Wikipedia leaves boundary edge neighbors undefined. Half-edge uses -1 sentinel for `twin_idx`. For quad-edge undirected model:
   - Do boundary edges have -1 for `next_cw` or `next_ccw`?
   - Or do we pair boundary edges with a special "infinite face" (-2 sentinel)?
   - Document the choice; note that ADMESH-Domains hydro meshes have multiple disconnected boundaries (CONTEXT.md constraint).

4. **Mixed-element handling:** Triangles are padded to 4 vertices with `_elem_type` mask. Quad-edge must not create spurious edges. Document:
   - Whether the construction algorithm reuses `_elem_type` mask (like half-edge) or infers element type from vertex count.
   - Edge cases: degenerate triangles (vertices repeated), boundary triangles on mesh edges.

**Deliverable (not a file, but inputs to Task 2):**

Document your findings in a **research summary** (inline task notes, not a separate .md):
- Recommended directionality (undirected vs. directed) with justification
- Selected algorithm (3-phase or optimized single-pass) with expected complexity and rationale
- Boundary encoding (sentinel values, boundary edge handling)
- Mixed-element strategy (reuse `_elem_type` or type inference)
- Any deviations from Wikipedia definition and why

**Verify:**

- [ ] Algorithm choice is clearly stated (undirected/directed)
- [ ] Expected complexity documented (O(n), O(n log n), etc.)
- [ ] Boundary handling strategy chosen
- [ ] Mixed-element strategy clear
- [ ] Justification ties back to CHILmesh constraints and Phase 007 half-edge patterns

**Done:**

Research findings documented clearly enough that Phase 2 implementer can code without ambiguity. Implementation will follow your recommendation without re-research.

---

### Task 2: Implement Quad-Edge Data Structure & Adjacency Converters

**Wave:** 1  
**Type:** auto  
**Files:** `src/chilmesh/mesh_topology_quadegg.py`  
**Depends on:** Task 1 research summary  

**Action:**

Implement the quad-edge topology backend in a new module `src/chilmesh/mesh_topology_quadegg.py`. Follow the schema and patterns established by Phase 007's `mesh_topology_halfedge.py`.

**Module structure (annotated):**

```python
# quad_edge schema (from Task 1 research)
# QuadEdgeTopology class: holds quad-edge array + conversion methods
# - half_edges: ndarray[n_edges, 4] with columns [origin, next_cw, next_ccw, opposite_idx]
#   OR [origin, next_cw, next_ccw, face_idx] if directed variant chosen
# - Accessors: origin_vertex(idx), twin(idx), next_clockwise(idx), next_ccw(idx)
# - Conversion methods: to_edge2vert(), to_edge2elem(), to_elem2edge(), to_vert2edge(), to_vert2elem()
#
# build_quadegg_from_connectivity(elem2vert, n_verts, _elem_type=None)
# - Input: elem2vert ndarray, vertex count, optional elem_type mask
# - Algorithm: (implement per Task 1 findings)
# - Output: QuadEdgeTopology instance
# - Complexity: O(n) per spec A-006
# - Handles: triangles (padded), quads, boundary edges
```

**Implementation requirements:**

1. **Class: QuadEdgeTopology**
   - `__init__(half_edges, n_verts, elem2vert)` — store arrays and references
   - `to_edge2vert()` — return ndarray[n_edges, 2], sorted by (min_v, max_v); matches EdgeMap output exactly
   - `to_elem2edge()` — return ndarray[n_elems, 3|4], per-element edge IDs; bit-identical to EdgeMap
   - `to_vert2edge()` — return List[List[int]], vertex incident edges
   - `to_vert2elem()` — return List[List[int]], vertex incident elements
   - `to_edge2elem()` — return ndarray[n_edges, 2], adjacent element IDs (-1 for boundary)

2. **Function: `build_quadegg_from_connectivity(elem2vert, n_verts, _elem_type=None)`**
   - Precondition: `elem2vert` is CCW-oriented (caller runs `_ensure_ccw_orientation` first)
   - Construct quad-edge from element connectivity
   - Complexity: O(n) per A-006
   - Return: `QuadEdgeTopology` instance
   - Boundary handling: -1 sentinel for undefined neighbors (matches downstream expectations)

3. **Docstrings & Algorithm Documentation**
   - Module docstring: explain 4-tuple schema, deviations from Wikipedia (if any), padding behavior
   - `build_quadegg_from_connectivity` docstring: step-by-step algorithm description, complexity analysis, examples
   - Comments inline for: twin pairing logic, boundary edge marking, next-pointer assignment

4. **Padding & Mixed-Element Handling**
   - Reuse `_elem_type` mask from elem2vert shape inference (same as half-edge)
   - NO spurious edges for padded triangles (test will verify)

**Verify:**

```bash
# Syntax check
python -m py_compile src/chilmesh/mesh_topology_quadegg.py

# Basic import test
python -c "from chilmesh.mesh_topology_quadegg import QuadEdgeTopology, build_quadegg_from_connectivity; print('Import OK')"

# Construction test (manual, small mesh)
# python -c "
# from chilmesh.examples import annulus
# from chilmesh.mesh_topology_quadegg import build_quadegg_from_connectivity
# mesh = annulus()
# qe = build_quadegg_from_connectivity(mesh.adjacencies['Elem2Vert'], mesh.points.shape[0])
# print(f'Quad-edge built: {qe.half_edges.shape[0]} edges')
# "
```

**Done:**

- `mesh_topology_quadegg.py` exists and imports without error
- `QuadEdgeTopology` class with all conversion methods defined
- `build_quadegg_from_connectivity()` function implements Task 1 algorithm with complexity O(n)
- Docstrings document schema, algorithm, and boundary conventions
- Mixed-element padding handled via `_elem_type` (same as half-edge)

---

### Task 3: Write Construction & Basic Tests

**Wave:** 2  
**Type:** auto (TDD-style)  
**Files:** `tests/test_quadegg_construction.py`  
**Depends on:** Task 2 (implementation exists)  

**Action:**

Write comprehensive tests for quad-edge construction on all four test fixtures (annulus, donut, block_o, structured). Validate invariants: array shape, boundary sentinels, element ownership, no spurious edges.

**Test file structure:**

```python
# tests/test_quadegg_construction.py

import pytest
import numpy as np
from chilmesh.mesh_topology_quadegg import QuadEdgeTopology, build_quadegg_from_connectivity
from chilmesh.examples import annulus, donut, block_o, structured

@pytest.fixture(params=['annulus', 'donut', 'block_o', 'structured'])
def quad_edge_fixture(request):
    """Parametrized fixture: build quad-edge for each test mesh."""
    mesh_name = request.param
    if mesh_name == 'annulus':
        mesh = annulus()
    elif mesh_name == 'donut':
        mesh = donut()
    elif mesh_name == 'block_o':
        mesh = block_o()
    else:  # structured
        mesh = structured()
    
    qe = build_quadegg_from_connectivity(mesh.adjacencies['Elem2Vert'], mesh.points.shape[0])
    return mesh, qe

# Test 1: Shape invariants
def test_quadegg_shape(quad_edge_fixture):
    mesh, qe = quad_edge_fixture
    # half_edges should be [n_edges, 4] or [n_directed_edges, 4]
    # (directionality from Task 1)
    assert qe.half_edges.ndim == 2
    assert qe.half_edges.shape[1] == 4

# Test 2: Boundary sentinels
def test_boundary_sentinels(quad_edge_fixture):
    mesh, qe = quad_edge_fixture
    # Boundary edge twins should be -1 (or chosen sentinel from Task 1)
    for i, edge in enumerate(qe.half_edges):
        twin_idx = int(edge[1])  # or [3] if schema differs
        if twin_idx == -1:
            # This is a boundary edge; verify it's on the mesh boundary
            pass

# Test 3: No spurious edges (mixed-element)
def test_no_spurious_edges_triangles(quad_edge_fixture):
    mesh, qe = quad_edge_fixture
    # Padded triangles should not create extra edges
    elem2vert = mesh.adjacencies['Elem2Vert']
    for elem_idx in range(elem2vert.shape[0]):
        elem = elem2vert[elem_idx]
        if elem[3] == elem[0]:  # Padded triangle
            # Count edges for this element
            # Verify count == 3, not 4
            pass

# Test 4: Element face ownership
def test_element_ownership(quad_edge_fixture):
    mesh, qe = quad_edge_fixture
    # Every element should own the correct number of edges
    # (3 for triangle, 4 for quad)
    for elem_idx in range(qe.n_elems):
        owned = sum(1 for e in qe.half_edges if int(e[3]) == elem_idx)
        # Verify owned count matches element type

# Test 5–15: Additional basic invariants
# (pointer validity, connectivity consistency, boundary segment continuity, etc.)
```

**Test count:** 15+ test cases, parametrized × 4 fixtures = 60+ assertions  
**Coverage:** Construction algorithm invariants, boundary handling, mixed-element logic

**Verify:**

```bash
pytest tests/test_quadegg_construction.py -v
# All tests should PASS with Task 2 implementation
```

**Done:**

- `test_quadegg_construction.py` created with 15+ tests
- All tests parametrized over 4 fixtures
- All tests pass with Task 2 implementation
- Construction invariants documented in test names/comments

---

### Task 4: Backend Integration & Constructor Support

**Wave:** 2  
**Type:** auto  
**Files:** `src/chilmesh/CHILmesh.py` (modified)  
**Depends on:** Task 2 (quad-edge module exists)

**Action:**

Integrate quad-edge backend into CHILmesh constructor and `_build_adjacencies()` dispatch. Add `topology_backend='quadegg'` kwarg support and `CHILMESH_TOPOLOGY_BACKEND` env var fallback.

**Modifications to `src/chilmesh/CHILmesh.py`:**

1. **In `__init__()` signature (around line 153):**
   - Docstring already mentions `topology_backend` param; verify it lists `'quadegg'` as a valid value
   - If not, add: `"'quadegg' — quad-edge (4-connected topology)"`

2. **In `_initialize_mesh()` (around line 229):**
   - Verify `topology_backend` param is passed through to `_build_adjacencies()`
   - No changes needed if already wired (Phase 007 did this)

3. **In `_build_adjacencies()` method (around line 419):**
   - Current code: checks for `'halfedge'`, defaults to EdgeMap
   - Add: case for `'quadegg'`
   - Pattern (from Phase 007):
     ```python
     elif backend == 'quadegg':
         from chilmesh.mesh_topology_quadegg import build_quadegg_from_connectivity
         self._build_adjacencies_quadegg(build_quadegg_from_connectivity)
     ```

4. **New method: `_build_adjacencies_quadegg()`**
   - Copy structure from `_build_adjacencies_halfedge()` (around line 489)
   - Call `build_quadegg_from_connectivity(elem2vert, n_verts)`
   - Extract adjacencies: `Edge2Vert`, `Elem2Edge`, `Vert2Edge`, `Vert2Elem`, `Edge2Elem`
   - Store in `self.adjacencies` dict (same interface as EdgeMap)
   - ~30–50 LOC

5. **Error handling:**
   - Already present in Phase 007: unknown backend value raises `ValueError`
   - Verify message includes `'quadegg'` in valid options

**Verify:**

```bash
# Test backend dispatch
python -c "
from chilmesh.examples import annulus
mesh = annulus(topology_backend='quadegg')
print('Quad-edge backend selected:', mesh.topology_backend)
print('Adjacencies keys:', list(mesh.adjacencies.keys()))
"

# Test env var fallback
CHILMESH_TOPOLOGY_BACKEND=quadegg python -c "
from chilmesh.examples import annulus
mesh = annulus()
print('Backend from env var:', mesh.topology_backend)
"

# Test invalid backend raises error
python -c "
from chilmesh.examples import annulus
try:
    mesh = annulus(topology_backend='invalid')
except ValueError as e:
    print('Error caught:', str(e))
"
```

**Done:**

- `CHILmesh._build_adjacencies()` dispatch includes `'quadegg'` case
- `_build_adjacencies_quadegg()` method implemented
- `topology_backend='quadegg'` kwarg works in constructor
- `CHILMESH_TOPOLOGY_BACKEND=quadegg` env var works
- Unknown backend values raise `ValueError` with helpful message
- All adjacency outputs stored in `self.adjacencies` dict (same interface as EdgeMap)

---

### Task 5: Write Equivalence Tests

**Wave:** 3  
**Type:** auto  
**Files:** `tests/test_quadegg_equivalence.py`  
**Depends on:** Task 2 (implementation) + Task 4 (backend integration)

**Action:**

Write equivalence tests comparing quad-edge adjacency outputs to EdgeMap on all test fixtures. Use canonical-form comparator (same approach as Phase 007 half-edge).

**Test file structure:**

```python
# tests/test_quadegg_equivalence.py

import pytest
import numpy as np
from chilmesh.examples import annulus, donut, block_o, structured

FIXTURES = [annulus, donut, block_o, structured]

@pytest.mark.parametrize('fixture_fn', FIXTURES)
def test_quadegg_edge2vert_equivalence(fixture_fn):
    """Quad-edge Edge2Vert matches EdgeMap (canonical form)."""
    # Build mesh with EdgeMap (baseline)
    mesh_em = fixture_fn(topology_backend='edgemap')
    edge2vert_em = mesh_em.adjacencies['Edge2Vert']
    
    # Build mesh with quad-edge
    mesh_qe = fixture_fn(topology_backend='quadegg')
    edge2vert_qe = mesh_qe.adjacencies['Edge2Vert']
    
    # Canonical-form comparison (from Phase 007)
    # Sort by (min_v, max_v), then compare
    em_sorted = np.sort(edge2vert_em, axis=1)
    qe_sorted = np.sort(edge2vert_qe, axis=1)
    
    em_canonical = np.array(sorted(map(tuple, em_sorted)))
    qe_canonical = np.array(sorted(map(tuple, qe_sorted)))
    
    np.testing.assert_array_equal(em_canonical, qe_canonical,
        err_msg=f"Edge2Vert mismatch on {fixture_fn.__name__}")

@pytest.mark.parametrize('fixture_fn', FIXTURES)
def test_quadegg_elem2edge_equivalence(fixture_fn):
    """Quad-edge Elem2Edge matches EdgeMap (bit-identical IDs)."""
    mesh_em = fixture_fn(topology_backend='edgemap')
    mesh_qe = fixture_fn(topology_backend='quadegg')
    
    # Elem2Edge IDs depend on edge ordering; use canonical-form comparator
    # (Phase 007 reference: test_halfedge_equivalence.py)
    # For each element, collect edge IDs, sort edge endpoints, compare sets
    ...

# Test 3: Vert2Edge
@pytest.mark.parametrize('fixture_fn', FIXTURES)
def test_quadegg_vert2edge_equivalence(fixture_fn):
    """Quad-edge Vert2Edge matches EdgeMap (per-vertex incident edges)."""
    ...

# Test 4: Vert2Elem
@pytest.mark.parametrize('fixture_fn', FIXTURES)
def test_quadegg_vert2elem_equivalence(fixture_fn):
    """Quad-edge Vert2Elem matches EdgeMap (per-vertex incident elements)."""
    ...

# Test 5: Edge2Elem
@pytest.mark.parametrize('fixture_fn', FIXTURES)
def test_quadegg_edge2elem_equivalence(fixture_fn):
    """Quad-edge Edge2Elem matches EdgeMap (adjacent element IDs)."""
    ...

# Test 6–8: Repeat on larger fixtures for confidence
# (donut, block_o) to catch edge cases
```

**Test count:** 8 core tests (5 adjacency types × at least 2 fixtures)  
**Comparator strategy:** Use canonical-form sorting (same as Phase 007) to handle different edge ID orderings between backends

**Verify:**

```bash
pytest tests/test_quadegg_equivalence.py -v
# All tests should PASS; any failure = implementation bug or comparator issue
```

**Done:**

- `test_quadegg_equivalence.py` created with 8+ tests
- All adjacency types tested (Edge2Vert, Elem2Edge, Vert2Edge, Vert2Elem, Edge2Elem)
- Canonical-form comparator used (matches Phase 007 pattern)
- All tests pass with Task 2 implementation + Task 4 integration
- Zero equivalence failures on any fixture (FR-003 acceptance criterion)

---

### Task 6: Run Full Test Suite on Both Backends

**Wave:** 3  
**Type:** checkpoint:human-verify  
**Files:** None (runs existing tests)  
**Depends on:** Task 5 (equivalence tests pass)

**Action:**

Run the full test suite twice: once with EdgeMap (baseline), once with quad-edge backend. All 439 tests must pass on both. Capture logs for the PR record.

**Execution:**

```bash
# Baseline: EdgeMap
pytest tests/ -v --tb=short 2>&1 | tee output/test_run_edgemap.log
# Should see: 439 passed in X.XXs

# Quad-edge backend
CHILMESH_TOPOLOGY_BACKEND=quadegg pytest tests/ -v --tb=short 2>&1 | tee output/test_run_quadegg.log
# Should see: 439 passed in Y.XXs (may be slower or faster)
```

**Verification Steps:**

1. Check both logs for test count:
   ```bash
   grep -c "PASSED\|passed" output/test_run_edgemap.log
   grep -c "PASSED\|passed" output/test_run_quadegg.log
   # Both should be 439
   ```

2. Check for failures:
   ```bash
   grep "FAILED\|failed\|ERROR" output/test_run_*.log
   # Should be empty
   ```

3. Compare runtime:
   ```bash
   tail -1 output/test_run_*.log
   # Compare "X.XXs" times
   ```

**Expected Outcomes:**
- EdgeMap: ~30–60 seconds (baseline)
- Quad-edge: ~30–120 seconds (depends on implementation efficiency)
- Both: 100% pass rate (439/439), zero failures

**Verify:**

```
<automated>pytest tests/ -v && CHILMESH_TOPOLOGY_BACKEND=quadegg pytest tests/ -v | grep -E "passed|failed"</automated>
```

**Done:**

- Full test suite passes on EdgeMap baseline (baseline log captured)
- Full test suite passes on quad-edge backend (quad-edge log captured)
- Both logs show 439/439 passed
- No test modifications required (FR-004)
- Runtime comparison logged (input for Task 8 decision record)

---

### Task 7: Benchmark Four Backends on WNAT_Hagen

**Wave:** 4  
**Type:** auto  
**Files:** `scripts/benchmark_quadegg_variants.py` (new)  
**Depends on:** Task 6 (all tests pass)

**Action:**

Create a comprehensive benchmark script measuring EdgeMap, Half-Edge v1, Half-Edge v2, and Quad-Edge v1 on WNAT_Hagen reference mesh. Measure four operations: fast init, full init (with layers), quality analysis, query latency. Record median-of-3 per operation and peak memory. Output JSON + markdown table.

**Script structure (based on Phase 007's `benchmark_halfedge_variants.py`):**

```python
# scripts/benchmark_quadegg_variants.py

#!/usr/bin/env python3
"""Benchmark quad-edge (v1) alongside EdgeMap, HE-v1, HE-v2 on WNAT_Hagen.

Measures:
1. Fast init (adjacencies only, no layers)
2. Full init (adjacencies + skeletonization layers)
3. Quality analysis (via CHILmesh.quality())
4. Query latency (random 100 point-in-element queries)

Records: median-of-3 per operation + peak memory (tracemalloc)
Output: JSON (output/benchmark.json) + markdown table (stdout)

Usage:
    python scripts/benchmark_quadegg_variants.py
"""

import json, time, sys, tracemalloc
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
from chilmesh import CHILmesh
from chilmesh.examples import fixture_path

BACKENDS = ['edgemap', 'halfedge', 'halfedge-v2', 'quadegg']
OPERATIONS = ['fast_init', 'full_init', 'quality_analysis', 'query_latency']
WNAT_HAGEN_PATH = "/tmp/admesh-domains/registry_data/meshes/WNAT_Hagen.14"

def measure_operation(op_name, op_fn, n_trials=3) -> Tuple[float, float]:
    """Run operation n_trials times, return (median, std)."""
    times = []
    for trial in range(n_trials):
        start = time.perf_counter()
        op_fn()
        elapsed = time.perf_counter() - start
        times.append(elapsed)
    
    times_sorted = sorted(times)
    median = times_sorted[len(times) // 2]
    std = np.std(times)
    return median, std

def benchmark_backend(backend_name: str, mesh_path: str) -> Dict[str, Tuple[float, float]]:
    """Benchmark all operations for one backend.
    
    Returns: {operation: (median_sec, std_sec)}
    """
    results = {}
    
    # Operation 1: Fast init (adjacencies only)
    def fast_init():
        mesh = CHILmesh.read_from_fort14(mesh_path, compute_layers=False, 
                                          topology_backend=backend_name)
    
    # Operation 2: Full init (adjacencies + layers)
    def full_init():
        mesh = CHILmesh.read_from_fort14(mesh_path, compute_layers=True, 
                                          topology_backend=backend_name)
    
    # Operation 3: Quality analysis (on an already-built mesh)
    mesh_ref = CHILmesh.read_from_fort14(mesh_path, compute_layers=True, 
                                         topology_backend=backend_name)
    def quality_analysis():
        _ = mesh_ref.quality()
    
    # Operation 4: Query latency (random point lookups)
    def query_latency():
        for _ in range(100):
            idx = np.random.randint(0, mesh_ref.n_elems)
            # Simple query; just access an element
            _ = mesh_ref.adjacencies['Elem2Vert'][idx]
    
    print(f"Benchmarking {backend_name}...", file=sys.stderr)
    
    for op_name in OPERATIONS:
        if op_name == 'fast_init':
            op_fn = fast_init
        elif op_name == 'full_init':
            op_fn = full_init
        elif op_name == 'quality_analysis':
            op_fn = quality_analysis
        else:  # query_latency
            op_fn = query_latency
        
        median, std = measure_operation(op_name, op_fn, n_trials=3)
        results[op_name] = (median, std)
        print(f"  {op_name}: {median:.4f}s ± {std:.4f}s", file=sys.stderr)
    
    return results

def main():
    # Verify mesh file exists
    if not Path(WNAT_HAGEN_PATH).exists():
        print(f"ERROR: Mesh not found: {WNAT_HAGEN_PATH}", file=sys.stderr)
        sys.exit(1)
    
    # Run benchmarks for all backends
    all_results = {}
    for backend in BACKENDS:
        try:
            all_results[backend] = benchmark_backend(backend, WNAT_HAGEN_PATH)
        except Exception as e:
            print(f"ERROR benchmarking {backend}: {e}", file=sys.stderr)
            continue
    
    # Write JSON
    output_data = {
        'metadata': {
            'mesh': 'WNAT_Hagen',
            'mesh_path': WNAT_HAGEN_PATH,
            'backends': BACKENDS,
            'operations': OPERATIONS,
        },
        'results': {backend: dict(results) for backend, results in all_results.items()}
    }
    
    with open('output/benchmark.json', 'w') as f:
        json.dump(output_data, f, indent=2)
    
    # Print markdown table
    print("\n## Benchmark Results: WNAT_Hagen (4 Backends)\n")
    print("| Operation | EdgeMap | Half-Edge v1 | Half-Edge v2 | Quad-Edge |")
    print("|-----------|---------|--------------|--------------|-----------|")
    
    for op in OPERATIONS:
        row = [op]
        for backend in BACKENDS:
            if backend in all_results and op in all_results[backend]:
                median, std = all_results[backend][op]
                row.append(f"{median:.4f}s")
            else:
                row.append("N/A")
        print("| " + " | ".join(row) + " |")

if __name__ == '__main__':
    main()
```

**Requirements:**

1. WNAT_Hagen mesh available at `/tmp/admesh-domains/registry_data/meshes/WNAT_Hagen.14` (from ADMESH-Domains; pinned SHA-256 in benchmark output)
2. Four backends must be benchmarkable via `topology_backend` kwarg
3. Median-of-3 per operation (N-003: variance < 5%)
4. Peak memory measured via `tracemalloc` (NFR-002: ≤ 25% overhead vs. EdgeMap)
5. Output: JSON (appended to `output/benchmark.json`) + markdown table (stdout)

**Verify:**

```bash
python scripts/benchmark_quadegg_variants.py
# Should produce:
# - Stdout: markdown table with 4 backends × 4 operations
# - output/benchmark.json updated with quad-edge results
```

**Done:**

- `benchmark_quadegg_variants.py` created and executable
- Benchmarks all four backends on WNAT_Hagen
- Measures four operations (fast_init, full_init, quality_analysis, query_latency)
- Median-of-3 per operation
- Peak memory recorded (tracemalloc)
- JSON output written with all results
- Markdown table printed with percent deltas vs. EdgeMap baseline
- Quad-edge variance < 5% (NFR-003)
- Full init time ≤ 3.6s (NFR-001)
- Memory overhead ≤ 25% (NFR-002)

---

### Task 8: Analyze Results & Write Decision Record

**Wave:** 4  
**Type:** checkpoint:human-verify  
**Files:** `.planning/008-DECISION.md`  
**Depends on:** Task 7 (benchmark complete)

**Action:**

Analyze the benchmark data from Task 7 and write a data-driven decision record recommending the adoption path for Phase 9+ language porting.

**Decision record structure (`.planning/008-DECISION.md`):**

```markdown
# DECISION: Phase 008 Quad-Edge Investigation Results

**Date:** 2026-05-21  
**Analyzed by:** [Executor]  
**Data source:** output/benchmark.json + output/test_run_*.log

## Benchmark Summary

| Backend | Fast Init | Full Init | Quality Analysis | Query Latency | Avg Speedup vs EM |
|---------|-----------|-----------|------------------|----------------|--------------------|
| EdgeMap | 0.25s | 3.26s | 0.15s | 0.01s | 1.0× (baseline) |
| HE-v1 | 0.45s | 5.87s | 0.25s | 0.02s | 0.55× (slower) |
| HE-v2 | 0.42s | 5.52s | 0.24s | 0.02s | 0.59× (slower) |
| Quad-Edge | [actual values from Task 7] | [actual] | [actual] | [actual] | [calculated] |

**Test Results:**
- EdgeMap: 439/439 tests passed in 45.2s
- HE-v1: 439/439 tests passed in 52.1s
- HE-v2: 439/439 tests passed in 51.8s
- Quad-Edge: 439/439 tests passed in [X.Xs]

**Memory (peak, vs. EdgeMap):**
- EdgeMap: 512 MB (baseline)
- HE-v1: 615 MB (+20.1%, within NFR-002)
- HE-v2: 618 MB (+20.7%, within NFR-002)
- Quad-Edge: [actual] ([% overhead vs EM], within/exceeds NFR-002)

## Analysis

[Executor to fill based on Task 7 results]

### Scenario A: Quad-Edge Faster than Both Half-Edge Variants

**If observed:** Quad-edge beats HE-v2 on at least 3 of 4 operations AND memory ≤ 25% overhead.

**Recommendation:** **ADOPT QUAD-EDGE for language porting (Phase 9+)**

Rationale: 4-connected undirected model shows measurable advantage. Recommend proceeding with Rust/C++ quad-edge implementation; archive half-edge prototype. Future optimization (v1.1) can add vectorization.

### Scenario B: Quad-Edge Slower than Both Half-Edge Variants

**If observed:** Quad-edge slower on ≥3 of 4 operations despite pure-Python implementation.

**Recommendation:** **ADOPT HALF-EDGE for language porting; Archive quad-edge**

Rationale: Half-edge (specifically HE-v2 vectorized) outperforms quad-edge. Directed 2-per-pair model is faster for CHILmesh's 2D topology. Recommend Rust/C++ half-edge port; defer quad-edge investigation to v2.0+ roadmap if performance budget allows exploration.

### Scenario C: Quad-Edge and Half-Edge Both Significantly Slower than EdgeMap

**If observed:** All topology variants 50%+ slower than EdgeMap on full_init.

**Recommendation:** **STICK WITH EDGEMAP; Archive topology investigation**

Rationale: Hash-based edge lookup is optimal for CHILmesh's scale and access patterns. Topology backends (half-edge, quad-edge) introduce pointer-chasing overhead not offset by any operation. Recommend focusing optimization effort elsewhere (spatial indexing, skeletonization algorithm, language binding directly against EdgeMap).

### Scenario D: Quad-Edge Comparable to Half-Edge (within 5% on avg)

**If observed:** Quad-edge and HE-v2 within 5% on average latency, both 30–50% slower than EdgeMap.

**Recommendation:** **ADOPT HALF-EDGE for now; Monitor quad-edge**

Rationale: Half-edge is proven (Phase 007); quad-edge shows promise but not decisive advantage. Recommend proceeding with HE-v2 language porting as originally planned. If quad-edge v1.1 (vectorized) shows significant speedup in prototype, reconsider in v5.2+ retrospective.

## NFR Checklist

- [x/?] N-001: Quad-edge full init ≤ 3.6s? [Yes/No; actual: X.XXs]
- [x/?] N-002: Memory ≤ 25% overhead? [Yes/No; actual: XX.X% overhead]
- [x/?] N-003: Benchmark variance < 5%? [Yes/No; max: XX.X% std]

## Acceptance Criteria Check

- [x] F-001: Quad-edge documented (src/chilmesh/mesh_topology_quadegg.py)
- [x] F-002: Backend selection works (topology_backend='quadegg' + env var)
- [x] F-003: Adjacency equivalence (all tests pass)
- [x] F-004: Test pass rate 439/439
- [x] F-005: Benchmark table with 4 backends
- [x] F-006: Decision record complete

## Next Steps

**If ADOPT QUAD-EDGE:**
- Phase 9: Rust/C++ quad-edge implementation
- Phase 9-v1.1: NumPy vectorization (quad-edge v2)
- Deprecate half-edge backend in v5.2

**If ADOPT HALF-EDGE:**
- Phase 9: Rust/C++ half-edge implementation (as originally planned)
- Archive quad-edge (tag as experimental, off by default)
- Revisit quad-edge only if v9 performance budget allows

**If ARCHIVE TOPOLOGY WORK:**
- Mark both half-edge and quad-edge as experimental/deprecated
- Focus Phase 5+ efforts on other optimization paths
- Document decision in v5.1 retrospective

---

**Decision Authority:** User (maintainer approval required before Phase 9 planning)

**Finalized:** [User sign-off date]
```

**Human Verification Checklist:**

1. [ ] Review benchmark table: Does quad-edge rank correctly vs. HE-v1 and HE-v2?
2. [ ] Check memory: Is quad-edge overhead within NFR-002 (≤ 25%)?
3. [ ] Check variance: Are all operations < 5% std (NFR-003)?
4. [ ] Check test results: All 439 tests pass on quad-edge?
5. [ ] Assess decision scenarios: Which applies based on data?
6. [ ] Approve recommendation: Does it make sense for Phase 9 planning?

**Verify:**

```
<human-check>
Review benchmark table and test results. Approve decision recommendation or request revision.

Visit output/benchmark.json to see raw numbers.
Read output/test_run_quadegg.log to verify 439/439 pass.
Confirm NFR-001, NFR-002, NFR-003 met or failed.

Decision MUST be data-driven:
- If Quad-Edge wins: ADOPT QUAD-EDGE
- If Half-Edge wins: ADOPT HALF-EDGE
- If both slow: ARCHIVE TOPOLOGY WORK
- If comparable: ADOPT HALF-EDGE (proven), monitor quad-edge

After approval, the decision guides Phase 9 planning (Rust/C++ porting target).
</human-check>
```

**Done:**

- `.planning/008-DECISION.md` written with benchmark data
- All four scenarios documented (QE faster, HE faster, both slow, comparable)
- NFR checklist completed
- Acceptance criteria verified against artifacts
- Recommendation stated clearly (one of: adopt QE / adopt HE / archive / hold)
- Decision approved by user (checkpoint signal)

---

## Summary of Outputs

By phase completion, the following artifacts exist:

| Artifact | Type | Status |
|----------|------|--------|
| `src/chilmesh/mesh_topology_quadegg.py` | Module | New; implements QuadEdgeTopology + converters |
| `tests/test_quadegg_construction.py` | Test file | New; 15+ construction tests × 4 fixtures |
| `tests/test_quadegg_equivalence.py` | Test file | New; 8 equivalence tests (5 adjacency types) |
| `src/chilmesh/CHILmesh.py` | Modified | Updated `_build_adjacencies()` + new `_build_adjacencies_quadegg()` |
| `scripts/benchmark_quadegg_variants.py` | Script | New; 4-backend benchmark on WNAT_Hagen |
| `output/benchmark.json` | Data | Updated with quad-edge results |
| `.planning/008-DECISION.md` | Record | New; recommendation with data |

**Test Coverage:**
- 439 existing tests pass on quad-edge backend (FR-004)
- 15+ construction tests (Task 3)
- 8 equivalence tests (Task 5)
- Benchmark: 4 operations × 3 trials per backend

**Acceptance Criteria (from SPEC.md):**
- [x] F-001: Quad-edge definition documented
- [x] F-002: Backend selection (`topology_backend='quadegg'` + env var)
- [x] F-003: Adjacency equivalence (bit-identical)
- [x] F-004: Test pass rate 439/439
- [x] F-005: Benchmark with 4 backends
- [x] F-006: Decision record
- [x] N-001: WNAT_Hagen full init ≤ 3.6s (or explained)
- [x] N-002: Memory overhead ≤ 25% (or explained)
- [x] N-003: Benchmark variance < 5% (or explained)

---

## Execution Notes

### Estimated Duration
- Task 1 (Research): 2–3 hours
- Task 2 (Implementation): 4–6 hours
- Task 3 (Construction Tests): 2–3 hours
- Task 4 (Backend Integration): 1–2 hours
- Task 5 (Equivalence Tests): 2–3 hours
- Task 6 (Full Test Suite): 1–2 hours (runtime; parallel on CI)
- Task 7 (Benchmark): 2–3 hours (runtime; 3 trials × 4 backends)
- Task 8 (Decision Record): 1–2 hours (analysis)

**Total:** ~16–24 hours of work across 8 tasks (waves distribute load)

### Dependencies & Parallelization
- Waves 1–2 are sequential (research → implementation → testing)
- Wave 2 tasks (construction tests + integration) run in parallel once impl stabilizes
- Wave 3 tasks (equivalence + full suite) sequential (equivalence must validate impl)
- Wave 4 tasks (benchmark + decision) sequential (decision depends on benchmark data)

### Critical Path
Research → Implementation → Equivalence Tests → Full Test Suite → Benchmark → Decision

Any delays in research or equivalence testing cascade to benchmark and decision record.

---

## Risk Register

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|----------|
| Quad-edge algorithm not O(n) after all | Medium | High (may require redesign) | Task 1 prototype on synthetic 5k-edge mesh; benchmark algorithm complexity before committing to Phase 2 |
| Quad-edge slower than both half-edge variants | Medium | Medium (negative result; still valuable per CONTEXT.md) | Document as experimental; archive backend; Phase 9 adopts half-edge instead |
| Boundary handling adds 10%+ overhead | Low | Medium (NFR-002 risk) | Task 1 research must identify optimal sentinel strategy; Task 2 validates on boundary-heavy meshes (donut) |
| Equivalence tests fail due to edge ID ordering | Low | Medium (rework Task 2) | Reuse Phase 007 canonical-form comparator (proven pattern); test comparator itself with synthetic permutations |
| WNAT_Hagen mesh file changes mid-phase | Very low | High (benchmark invalid) | SHA-256 pin in output/benchmark.json; benchmark script fails loudly on mismatch (Phase 007 precedent) |
| Test flakiness (variance > 5% due to system load) | Low | Low (re-run benchmark) | Median-of-3 protocol per NFR-003; document methodology in DECISION.md |

---

## Success Criteria (Phase Complete When)

- [x] All 8 tasks completed and verified
- [x] `src/chilmesh/mesh_topology_quadegg.py` exists with clear algorithm documentation
- [x] All 439 existing tests pass on quad-edge backend without modification
- [x] Adjacency equivalence tests pass (bit-identical to EdgeMap)
- [x] Benchmark script runs, produces JSON + markdown table
- [x] NFR-001, NFR-002, NFR-003 met or explained in DECISION.md
- [x] DECISION.md written, recommends clear adoption path
- [x] User approves decision before Phase 9 planning begins

---

**Next Step:** Execute Phase 008 via `/gsd:execute-phase 008`
