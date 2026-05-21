# Spec-Kit: Clarify Phase
# "Half-Edge Data Structure Investigation & Language Optimization Path"

**Branch:** `v5.0-optimize-port-w-half-edge`
**Date:** 2026-05-21
**Prerequisite:** `specs/007-v5.0-optimize-port-w-half-edge/spec.md` (Specify phase)
**Related Issues:** #137 (Half-edge investigation), #68 (Quad-edge alternative)

---

## CLARIFY ROUND 1 — Boundary Keeper

### Finding C-1: "Half-edge" without DCEL boundary contract is under-specified

**Evidence:** spec.md FR-001 names the backend "half-edge / DCEL" but does not specify whether the implementation is a half-edge list (HEL — every directed edge stored independently) or a true doubly-connected edge list (DCEL — twin pointers between half-edges encoded as references, not lookups).

**Impact:** Memory footprint, query latency, and conversion-to-adjacency code paths all differ. NFR-002 (≤25% memory overhead) cannot be verified against an unspecified target.

**Clarification:** Use DCEL with explicit twin / next / prev pointers stored as integer indices into a single contiguous array of half-edges. NumPy-friendly, cache-friendly, avoids Python object overhead per half-edge. Specifically:
- `half_edges`: `ndarray[n_halfedges, 4]` columns = `[origin_vertex, twin_idx, next_idx, face_idx]`
- Boundary half-edges have `twin_idx == -1` (consistent with `Edge2Elem` sentinel for boundaries).

This is the only memory model that meets NFR-002 in pure Python; per-object models (`HalfEdge` dataclass per directed edge) blow the 25% budget on WNAT_Hagen (300k half-edges × ~200B per Python object = 60MB just for half-edges).

---

### Finding C-2: "Optimized variant" in the benchmark table is unspecified

**Evidence:** spec.md FR-006 mandates three columns — current, half-edge v1, "optimized variant" — but does not say what optimization is being measured.

**Impact:** The decision-trigger column (SC-005) becomes unanchored. We could optimize v1 → v2 along multiple axes (Numba JIT, array layout, cython compile) and end up with no consistent comparison.

**Clarification:** "Optimized variant" is defined as half-edge v1 with one specific optimization applied: replacing the Python `for`-loop that walks half-edges around a face with a NumPy vectorized scatter/gather. Other optimization candidates (Numba, Cython, Rust) are explicitly out of scope for this spec — they belong to P3's "investigate port" trigger.

Concretely: variant 2 = "v1 with the `next_idx` walk replaced by `np.take(half_edges[:,2], face_starts)`". Pure NumPy. No new dependencies.

---

### Finding C-3: "WNAT_Hagen" is an external dependency without pin

**Evidence:** spec.md A-003 assumes WNAT_Hagen "remains available via ADMESH-Domains". The actual mesh file lives in a separate repo with no version pin. If ADMESH-Domains rebuilds the registry, the file may differ across benchmark runs.

**Impact:** SC-002 (≤ 3.6 s on WNAT_Hagen) is not reproducible across time. NFR-003 (variance < 5%) only makes sense if we're benchmarking the same input.

**Clarification:** Record the SHA-256 of the WNAT_Hagen.14 file in the benchmark JSON output. The benchmark script must verify this hash and fail loudly if it changes. Pin recorded once during the first run on v5.0; checked in to `output/benchmark.json` and read on every subsequent run.

---

### Finding C-4: "Bit-identical" adjacency requires a precise equality definition

**Evidence:** spec.md FR-003 says half-edge backend MUST produce "bit-identical" adjacency outputs. But `EdgeMap` returns edges in insertion order; half-edge construction emits edges in face-traversal order. The arrays come out re-ordered.

**Impact:** Naive `np.array_equal` comparison fails on row-permuted output even when the topology is semantically identical. Downstream tests that index by edge ID break.

**Clarification:** FR-003's "bit-identical" applies to **values, not row order**. Specifically:
- `EdgeMap.to_list()` output sorted by `(min_vert, max_vert)` MUST match half-edge backend's equivalent output, row-by-row.
- `Edge2Vert`, `Edge2Elem`: compared as sets-of-frozensets-of-tuples (order-agnostic).
- `Elem2Edge`, `Vert2Edge`, `Vert2Elem`: each row's set of edge IDs must match after re-indexing edges to the canonical sorted ordering.

The equivalence test (Task 2 in plan.md) implements this canonical-form comparator. Naive `np.array_equal` is documented as NOT the contract.

---

## CLARIFY ROUND 2 — Failure Analyst

### Finding F-1: Half-edge construction on a degenerate quad with collinear vertices

**Evidence:** B5 from v0.1.1 (the `+1e-12` epsilon guard in `interior_angles`) confirms degenerate quads exist in real meshes. Half-edge construction relies on consistent orientation (next half-edge in CCW order). A collinear quad has no defined CCW order.

**Failure mode:** Half-edge `next_idx` becomes ambiguous; downstream `_skeletonize()` produces different layer counts.

**Prevention:** Half-edge construction MUST run `_ensure_ccw_orientation()` (existing B4 fix) BEFORE assigning `next_idx`. Add a regression test using the synthetic degenerate quad fixture committed in v0.1.1 (`tests/test_interior_angles.py::test_quad_angles_no_nan_on_degenerate_quad`).

---

### Finding F-2: `CHILMESH_TOPOLOGY_BACKEND` env var set in pytest fixture but not unset

**Evidence:** `tests/conftest.py` uses session-level memoization. If a test sets `CHILMESH_TOPOLOGY_BACKEND=halfedge` mid-suite, the change persists into all subsequent tests.

**Failure mode:** Tests appear to pass under one backend but actually exercise the other one. False-green benchmark + test correlation.

**Prevention:** Add a `monkeypatch.setenv` / `monkeypatch.delenv` pattern. All backend-switching MUST go through a pytest fixture that resets the env var on teardown. Document this in `tests/test_halfedge_basic.py` header.

---

### Finding F-3: Memory measurement under garbage collection variance

**Evidence:** NFR-002 specifies "≤25% memory overhead". Python's GC behavior makes peak RSS measurement noisy — a single test fixture's lifetime can swing peak by 10MB.

**Failure mode:** Benchmark reports half-edge as 30% over budget on one run, 5% over on the next.

**Prevention:** Use `tracemalloc` (snapshot before / after construction) instead of `resource.getrusage` (RSS). Record allocated bytes, not RSS. Run three trials; report median.

---

### Finding F-4: Mixed-element padded triangle → spurious half-edge

**Evidence:** Mixed-element model uses vertex 0 in slot 3 to mark a triangle in a 4-column array. Naive half-edge construction iterates over 4 vertices per element → creates a 4th half-edge for triangles that doesn't exist.

**Failure mode:** Triangle elements get 4 half-edges, two of which point at vertex 0 erroneously. Adjacency outputs become inconsistent (regresses B3 from v0.1.1).

**Prevention:** Half-edge construction MUST check element type per row before iterating. Skip the padded slot for triangles. Reuse the existing `_elem_type` vectorized mask from B3 — do not re-implement.

---

## CLARIFY ROUND 3 — Seed Closer

### Q-1: Backend selection — env var OR constructor kwarg OR both?

**Options:**
- A: Env var only (`CHILMESH_TOPOLOGY_BACKEND=halfedge pytest`)
- B: Kwarg only (`CHILmesh(..., topology_backend='halfedge')`)
- C: Both (env var as default, kwarg overrides) ← **selected**

**Decision:** Option C. Env var lets the whole test suite run on the variant without modifying any test code (per spec FR-004, NO test modifications). Kwarg lets the benchmark script flip backends per-variant explicitly. Both is canonical for "library default vs. caller override" patterns.

---

### Q-2: Benchmark JSON schema — backward-compatible or breaking?

**Options:**
- A: Add new fields, keep old fields (`v0.4.1` shape compatible) ← **selected**
- B: New schema version, breaking change

**Decision:** Option A. Existing tooling (CI archival, future regression alerts) already parses the v0.4.1 shape. Add new fields (`backend_variants[]`, `mesh_sha256`, `tracemalloc_peak`) without removing existing keys. JSON consumers ignore unknown fields.

---

### Q-3: When does the DECISION record get created — this phase or next?

**Options:**
- A: This phase — DECISION written immediately on benchmark completion
- B: Next phase (`v5.1` follow-on) — DECISION is a separate deliverable ← **selected**

**Decision:** Option B. Aligns with CONTEXT.md decision #4 ("benchmark-driven decision record, defer publish"). This phase delivers the *data*; the next phase produces the *recommendation*. Cleaner audit trail: each phase has exactly one deliverable type (this phase = code + benchmark; next = decision document + roadmap update).

---

### Q-4: What if benchmark shows half-edge is *slower* by >10%?

**Options:**
- A: Treat as failure; do not commit half-edge backend
- B: Commit half-edge backend as documented "experimental, off by default" + DECISION record archives it ← **selected**

**Decision:** Option B. Even a negative result is valuable — it closes #137 with data. The backend stays in the tree but is documented as "not recommended; archived for reference". A future contributor can pick it up if a Rust/C++ port becomes viable (different performance characteristics in compiled languages).

---

## UPDATED AMBIGUITY SCORES (Post-Clarify)

```
Goal Clarity:        0.93  ✓ (was 0.85 — C-1, C-2 sharpened "half-edge" and "optimized")
Boundary Clarity:    0.95  ✓ (was 0.85 — Q-3, Q-4 lock the deliverable boundary)
Constraint Clarity:  0.92  ✓ (was 0.80 — F-3 specified tracemalloc; NFR-002 now measurable)
Acceptance Criteria: 0.94  ✓ (was 0.85 — C-4 defined "bit-identical" precisely)

Ambiguity: 1.0 − (0.35×0.93 + 0.25×0.95 + 0.20×0.92 + 0.20×0.94) = 0.066  ← gate ≤ 0.20 ✓
```

**Gate passed. Spec + Clarify complete. Ready for /speckit-plan.**

---

## Spec Updates Required (Folded back into spec.md before plan.md)

- [ ] FR-001: Add "DCEL with integer-index pointers" detail per C-1
- [ ] FR-006: Specify the "optimized variant" as the NumPy vectorized walk per C-2
- [ ] A-003: Add SHA-256 pin requirement per C-3
- [ ] FR-003: Reference canonical-form comparator per C-4
- [ ] Edge cases: add F-1 (degenerate quad), F-2 (env var leak), F-3 (GC variance), F-4 (padded triangle)
- [ ] Out of Scope: explicitly list "Numba/Cython optimization of v1" per C-2

These updates should be applied during /speckit-plan as the planner reads spec.md + clarify.md together.
