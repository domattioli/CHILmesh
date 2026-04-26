# Phase 2: Dynamic Mesh Operations Implementation Specification

**Branch:** `planning-optimize_modernize`  
**Date:** 2026-04-26  
**Phase:** 2 – Dynamic Operations API (Implementation Preview)  
**Methodology:** Spec-Kit  
**Prerequisite:** Phase 1 complete (decision memo + recommendations)

---

## Executive Summary

Phase 2 translates the decision memo from Phase 1 into a working dynamic mesh API. The focus is implementing `add_element()`, `remove_element()`, and supporting operations with **zero impact to skeletonization semantics**.

**Deliverables:** 
- Updated `CHILmesh` class with 5+ new public methods
- Comprehensive test suite (regression + dynamic ops)
- Documentation (docstrings + API guide)

**Success:** All tests pass, skeletonization output byte-identical on all fixtures.

---

## Specification

### 2.1 Goals

**Primary:** Implement a clean, efficient API for incremental mesh modification that MADMESHR can invoke during advancing-front generation.

**Secondary:** Maintain skeletonization invariants perfectly (audit Q3: disjoint cover, monotone-shrinking).

### 2.2 In-Scope Operations

| Operation | Signature | Purpose | Complexity Target |
|-----------|-----------|---------|------------------|
| `add_element(vertices, elem_type)` | `int` (new elem ID) | Insert quad/triangle | O(k) where k = elem degree |
| `remove_element(elem_id)` | `None` | Delete element, clean orphaned edges | O(k) |
| `add_vertex(x, y, z)` | `int` (new vert ID) | Insert isolated vertex | O(1) |
| `remove_vertex(vert_id, strategy)` | `None` | Delete vertex (merge/split strategies) | O(k) |
| `split_edge(edge_id, position)` | `(int, int)` (2 new elems) | Subdivide edge into 2 elements | O(k) |
| `swap_edge(edge_id)` | `None` | Flip edge (quad refinement) | O(1) |

### 2.3 Consistency Rules (Per Operation)

**add_element(vertices, elem_type):**
- Pre: All vertices exist; no isolated vertices created
- Post: New element added to connectivity_list; edges registered in adjacency; element NOT added to layers (rebuilt on next _skeletonize() call)
- Invariant: Mesh topology remains valid (no self-intersections checked; assumed well-formed input)

**remove_element(elem_id):**
- Pre: Element exists; not in use by other operations
- Post: Element removed; orphaned edges deleted; layer info stale (rebuilds on next _skeletonize())
- Invariant: Remaining elements form valid mesh

**Transactional Model (Batch):**
```python
with mesh.batch_operations():
    mesh.add_element(...)
    mesh.add_element(...)
    mesh.remove_element(...)
# Consistency check + layer rebuild happens once at exit
```

### 2.4 Test Coverage

**Regression Suite (Existing):**
- [ ] All 1X existing tests pass (0 regressions)
- [ ] test_layers_annulus, test_elem_type, etc. unchanged

**New Dynamic Ops Tests:**
- [ ] test_add_element_triangle (add to annulus, verify adjacencies)
- [ ] test_add_element_quad (add to block_o, verify adjacencies)
- [ ] test_remove_element (remove from annulus, verify orphaned edge cleanup)
- [ ] test_add_vertex_isolated (add unused vertex, does nothing visible)
- [ ] test_batch_operations (multiple adds, single rebuild)
- [ ] test_skeletonization_unchanged (add/remove, then _skeletonize(), verify same layers)
- [ ] test_fort14_roundtrip_after_dynamic_ops (modify mesh, export, reimport, byte-identical)

### 2.5 Acceptance Criteria

- [ ] 5+ new public methods implemented with type hints
- [ ] All 6 regression tests passing
- [ ] All 7 dynamic ops tests passing
- [ ] Docstrings for all new methods (parameters, returns, examples)
- [ ] Zero performance regression (full init time unchanged for unmodified meshes)
- [ ] Skeletonization output identical on all 4 fixtures after dynamic ops

### 2.6 Out of Scope (Phase 2)

- ❌ Performance optimization of `add_element()` (that's Phase 3)
- ❌ Element quality validation (assume well-formed input)
- ❌ Parallel mesh modification
- ❌ Incremental skeletonization (full rebuild on demand; Phase 3 may optimize)

---

## Timeline

**3–5 days** (gated on Phase 1 completion)

---

## Handoff to Phase 3

Once Phase 2 is complete:
1. Dynamic ops API is stable (not changing in Phase 3)
2. Phase 3 focuses on optimization: O(n log n) edge discovery, incremental layers
3. All Phase 2 tests must continue to pass

