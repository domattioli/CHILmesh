# Architecture — the graph engine

CHILmesh treats a mesh as a **graph**: vertices, edges, and elements are nodes,
and the adjacency tables are the edges between them. `_build_adjacencies`
(`src/chilmesh/CHILmesh.py`) assembles seven tables once; every query, smoother,
and the layerizer then reads them in constant or linear time.

## Adjacency tables

| Table | Shape | Maps |
|---|---|---|
| `Elem2Vert` | (n_elems, 3\|4) | element → vertices |
| `Edge2Vert` | (n_edges, 2) | edge → endpoints (canonical min,max) |
| `Elem2Edge` | (n_elems, 3\|4) | element → edges |
| `Edge2Elem` | (n_edges, 2) | edge → elements (−1 = boundary) |
| `Vert2Edge` | dict[int → set] | vertex → incident edges |
| `Vert2Elem` | dict[int → set] | vertex → incident elements |
| `EdgeMap` | hash | (v₀, v₁) → edge id |

Edge IDs use first-encounter (element-major, slot-minor) ordering to stay
consistent with the C++ half-edge backend, which reproduces these tables
bit-for-bit (`n_layers` parity holds across both — see
[`BENCHMARK.md`](BENCHMARK.md)).

## Complexities (n = element count)

| Operation | Cost | How |
|---|---|---|
| Edge lookup / dedup | **O(1)** | `EdgeMap` hash (was O(n²) pre-v0.2 — the 937× init speedup) |
| Adjacency build | **O(n log n)** | vectorized `np.unique` / argsort |
| Layerization (`_layerize`) | **O(n)** | concentric layer peel, ~1 s / 60k elems in Python; ~15× faster in C++ |
| Spatial query (`find_element`, `nearest_vertices`) | **O(log n)** | `cKDTree` |
| Vertex valence / 1-ring | **O(1) + O(degree)** | dict lookup |

## Layerization vs skeletonization

`_layerize` peels concentric element rings (the discrete, banded analogue of a
distance field). It is **not** medial-axis or skeleton extraction — those are
distinct constructs documented in [`CONCEPTS.md`](CONCEPTS.md). The ownership
boundary between CHILmesh (post-mesh analysis) and ADMESH (pre-mesh generation)
is set in [#221](https://github.com/domattioli/CHILmesh/issues/221).
