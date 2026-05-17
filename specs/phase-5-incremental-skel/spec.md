# Phase 5.2: Incremental Skeletonization for Dynamic Meshes

**Goal:** Update mesh layers in-place instead of full rebuild (3s → 50ms per edit)  
**Depends on:** Phase 5.0 (Mutation API)  
**Impact:** Interactive mesh editing, real-time layer feedback  
**Effort:** Design: 1 hour. Implementation: 6-8 hours

---

## Problem

Current `_skeletonize()` rebuilds all layers from scratch:
- Block_O (100k elements): ~3s per call
- Adaptive loop: 1000 edits × 3s = 50 minutes wasted

---

## API

### reskeletonize_local(elem_ids: ndarray, radius: int = 2)
Update layers only in neighborhood of edited elements.

**Performance target:** <50ms for radius-2 on 100k mesh

```python
m = MutableMesh(...)
new_elems = m.split_triangle(elem_id=5000)
m.reskeletonize_local(elem_ids=new_elems, radius=2)  # Fast local update
```

### skeletonize_diff(prev_state: dict) → dict
Return only layers that changed.

```python
prev_layers = deepcopy(m.layers)
m.reskeletonize_local(...)
deltas = m.skeletonize_diff(prev_state=prev_layers)
```

---

## Algorithm

**Local boundary detection:** Mark dirty elements → BFS outward until convergence

**Convergence heuristic:** Stop when 3 consecutive elements have unchanged layer assignment

**Fallback:** If layer count changes globally, full rebuild triggered (safe but slower)

---

## Success Criteria

- [ ] `reskeletonize_local()` implemented
- [ ] Performance: <50ms on 100k mesh
- [ ] Correctness: matches full rebuild
- [ ] Fallback triggers on global topology changes
- [ ] 1000 random edits without corruption
