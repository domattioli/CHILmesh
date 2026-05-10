# Quad-Edge Data Structure: CHILmesh Feasibility Study

**Date:** 2026-05-09  
**Author:** Claude Code (routine session)  
**Tracking:** GitHub issue #68  
**Status:** Research — no implementation committed

---

## 1. Background

The Quad-Edge data structure (Guibas & Stolfi, 1985) represents a 2-manifold subdivision by storing four directed half-edges per undirected edge.  The fundamental primitives are:

| Operation | Semantics |
|-----------|----------|
| `e.Org()` | Origin vertex of directed edge `e` |
| `e.Dest()` | Destination vertex |
| `e.Rot()` | Dual edge (rotates 90° into dual graph) |
| `e.Sym()` | Reverse of `e` (`e.Rot().Rot()`) |
| `e.Onext()` | Next edge around `e.Org()` CCW |
| `e.Lface()` | Left face of `e` |
| `splice(e, f)` | Atomic local topology modification |

Reference: Guibas, L. & Stolfi, J. (1985). "Primitives for the manipulation of general subdivisions and the computation of Voronoi diagrams." *ACM Trans. Graph.* 4(2):74-123.

---

## 2. Current CHILmesh Architecture

CHILmesh v0.2.0 uses flat numpy arrays and dicts:

```
Edge2Vert:  ndarray[n_edges, 2]       O(1) endpoint lookup (array index)
Elem2Edge:  ndarray[n_elems, 3|4]     O(1) element→edge
Edge2Elem:  ndarray[n_edges, 2]       O(1) edge→adjacent elements (-1 sentinel)
Vert2Edge:  Dict[int, Set[int]]       O(k) vertex star  
Vert2Elem:  Dict[int, Set[int]]       O(k) vertex star
EdgeMap:    Dict[tuple, int]          O(1) canonical (v0,v1)→edge_id
```

All Phase 1–4 optimizations landed in v0.2.0 and already eliminated O(n²) bottlenecks.

---

## 3. What Quad-Edge Would Give CHILmesh

### Advantages

| Feature | Benefit | Applicable CHILmesh use case |
|---------|---------|-----------------------------|
| Algebraic traversal (`Onext`, `Lnext`) | Elegant, error-free fan/ring traversal | Skeletonization layer peeling; advancing front |
| `splice` primitive | O(1) edge split / flip / insertion | Future edge-flip smoothing |
| Dual graph automatic | Voronoi vertices at no extra cost | Medial axis research |
| Topological invariants by construction | Manifold property maintained automatically | Fewer consistency checks needed |

### Disadvantages

| Factor | Impact |
|--------|--------|
| Memory: 4 directed half-edges per edge | ~4× edge storage vs current `Edge2Vert` array |
| Cache locality | Pointer-based QE is harder to vectorise than flat numpy arrays |
| Implementation complexity | `splice` is subtle; bugs are topologically destructive |
| API incompatibility | Would break `edge2vert(edge_ids)` array semantics |
| Python overhead | Pure Python QE traversal vs numpy vectorised ops is slower for batch ops |
| Mixed tri+quad meshes | QE assumes orientable 2-manifold; padded triangles in 4-col connectivity need special handling |

---

## 4. Operation-by-Operation Comparison

### Skeletonization (`_skeletonize`)

Current approach: boundary peeling using `Edge2Elem` sentinel scan + `Vert2Edge` dict lookup.

With QE: `Onext` ring traversal to collect boundary fans. Cleaner code, but skeletonization already runs in < 1s on all fixtures.  **No speed gain expected.**

### Advancing Front (`add_advancing_front_element`)

Current approach: appends to connectivity arrays, rebuilds adjacency on demand.

With QE: `splice(e, f)` would insert the new element in O(1) and automatically maintain topological consistency. **This is the strongest argument for QE** in CHILmesh.

### Quality / Smoothing (`elem_quality`, `smooth_mesh`)

Current approach: fully vectorised numpy on `points[connectivity_list]`.

With QE: traversal to gather vertex stars, then same numpy ops. **No benefit; likely slower** due to traversal overhead before the vectorised step.

### Edge Flipping (not yet implemented)

With QE: `splice` makes edge flipping a 5-line function.  Without QE, flipping requires surgical array edits across 5 arrays. **QE wins clearly** for any future edge-flip-based smoother.

---

## 5. Memory Impact (Block_O reference mesh: 52,774 vertices, 98,365 elements)

| Structure | Current (bytes) | With QE (bytes) | Delta |
|-----------|----------------|-----------------|-------|
| Edge records (~148k edges) | ~2.4 MB (Edge2Vert int32×2) | ~9.5 MB (4 half-edge ptrs × int32) | +7 MB |
| Adjacency dicts (Vert2Edge + Vert2Elem) | ~8 MB est. | eliminated | -8 MB |
| **Net** | | | ~0 MB — roughly neutral |

For the Block_O mesh, memory impact is neutral. For truly large meshes (1M+ elements), the QE pointer overhead may become significant.

---

## 6. Recommendation

### Short term (v0.2.x)

**Do not adopt Quad-Edge.** The performance wins are already realized (Phase 1–4). A full rewrite would:
- Break the stable public API
- Require re-certification of all 57+ tests
- Introduce subtle correctness risks in `splice`
- Not improve benchmark numbers for current use cases

### Medium term (v0.3.0)

**Targeted evaluation in the advancing-front module only.**  The advancing front (`add_advancing_front_element`, `remove_boundary_loop`) is the one subsystem where QE’s `splice` primitive would eliminate the most complexity. Consider:
1. Implement a minimal `HalfEdge` class (just `Org, Dest, Sym, Onext, splice`) as an internal structure used only within the advancing-front code path.
2. Keep the flat numpy arrays for all other operations (quality, smoothing, I/O).
3. Benchmark on Block_O advancing-front scenario before committing.

### Long term (v1.0+)

**Revisit if edge-flip smoothing is prioritized.** Edge flips are the canonical QE use case. If a Delaunay refinement or Lawson-flip smoother is ever added, adopting QE for that subsystem would be well-justified.

---

## 7. Open Questions for Maintainer

1. Is edge-flip smoothing planned for v0.3.0? If yes, QE becomes a stronger candidate.
2. What is the expected scale of advancing-front meshes? If front sizes routinely exceed 100k edges, `splice`-based insertion avoids `np.append` O(n) copies.
3. Are there plans to expose a dual graph (Voronoi) API? QE provides the dual automatically.

---

## References

- Guibas, L., Stolfi, J. (1985). Primitives for manipulation of general subdivisions. *ACM TOG* 4(2).
- Shewchuk, J. (1996). Triangle: Engineering a 2D quality mesh generator. *WACG*. (Uses QE internally)
- de Berg et al. (2008). *Computational Geometry: Algorithms and Applications*, Ch. 9 (half-edge)
- Wikipedia: https://en.wikipedia.org/wiki/Quad-edge
