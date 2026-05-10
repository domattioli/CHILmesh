# Phase 2: Dynamic Mesh Operations Implementation Specification

**Branch:** `planning-optimize_modernize`  
**Date:** 2026-04-26  
**Phase:** 2 – Dynamic Operations API (Implementation Preview)  
**Prerequisite:** Phase 1 complete (decision memo + recommendations)

---

## Executive Summary

Translates Phase 1 decision memo into working dynamic mesh API. Implements `add_element()`, `remove_element()`, and supporting operations with **zero impact to skeletonization semantics**.

**Deliverables:** 5+ new public methods; regression + dynamic ops test suite; docstrings + API guide.

**Success:** All tests pass; skeletonization output byte-identical on all fixtures.

---

## Specification

### Goals

**Primary:** Implement clean, efficient API for incremental mesh modification MADMESHR can invoke during advancing-front generation.

**Secondary:** Maintain skeletonization invariants exactly (audit Q3: disjoint cover, monotone-shrinking).

### In-Scope Operations

| Operation | Signature | Purpose | Complexity |
|-----------|-----------|---------|------------|
| `add_element(vertices, elem_type)` | `int` (new elem ID) | Insert quad/triangle | O(k) where k = elem degree |
| `remove_element(elem_id)` | `None` | Delete element, clean orphaned edges | O(k) |
| `add_vertex(x, y, z)` | `int` (new vert ID) | Insert isolated vertex | O(1) |
| `remove_vertex(vert_id, strategy)` | `None` | Delete vertex (merge/split strategies) | O(k) |
| `split_edge(edge_id, position)` | `(int, int)` (2 new elems) | Subdivide edge | O(k) |
| `swap_edge(edge_id)` | `None` | Flip edge (quad refinement) | O(1) |

### Consistency Rules

**add_element(vertices, elem_type):**
- Pre: All vertices exist; no isolated vertices created
- Post: New element added to connectivity_list; edges registered; element NOT added to layers (rebuilt on next `_skeletonize()`)

**remove_element(elem_id):**
- Pre: Element exists
- Post: Element removed; orphaned edges deleted; layer info stale (rebuilds on next `_skeletonize()`)

**Transactional Model:**
```python
with mesh.batch_operations():
    mesh.add_element(...)
    mesh.add_element(...)
    mesh.remove_element(...)
# Consistency check + layer rebuild once at exit
```

### Test Coverage

**Regression Suite:**
- [ ] All existing tests pass (0 regressions)

**New Dynamic Ops Tests:**
- [ ] test_add_element_triangle (add to annulus, verify adjacencies)
- [ ] test_add_element_quad (add to block_o, verify adjacencies)
- [ ] test_remove_element (remove from annulus, verify orphaned edge cleanup)
- [ ] test_add_vertex_isolated
- [ ] test_batch_operations (multiple adds, single rebuild)
- [ ] test_skeletonization_unchanged (add/remove → _skeletonize() → same layers)
- [ ] test_fort14_roundtrip_after_dynamic_ops

### Acceptance Criteria

- [ ] 5+ new public methods with type hints
- [ ] All 13+ tests passing (regression + dynamic ops)
- [ ] Zero performance regression (full init unchanged for unmodified meshes)
- [ ] Skeletonization output identical on all 4 fixtures after dynamic ops

### Out of Scope (Phase 2)

- ❌ Performance optimization of `add_element()` (Phase 3)
- ❌ Element quality validation (assume well-formed input)
- ❌ Incremental skeletonization (Phase 3 may optimize)

---

**Timeline:** 3–5 days (gated on Phase 1 completion)

**Handoff to Phase 3:** Dynamic ops API stable; Phase 3 focuses on O(n log n) edge discovery and incremental layers.

