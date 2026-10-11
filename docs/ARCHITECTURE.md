# Architecture: the graph engine

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
bit-for-bit (`n_layers` parity holds across both; see
[`BENCHMARK.md`](BENCHMARK.md)).

## Complexities (n = element count)

| Operation | Cost | How |
|---|---|---|
| Edge lookup / dedup | **O(1)** | `EdgeMap` hash (was O(n²) pre-v0.2; the 937× init speedup) |
| Adjacency build | **O(n log n)** | vectorized `np.unique` / argsort |
| Layer peel (`_peel`) | **O(n)** | concentric layer peel, ~1 s / 60k elems in Python; ~15× faster in C++ |
| Spatial query (`find_element`, `nearest_vertices`) | **O(log n)** | `cKDTree` |
| Vertex valence / 1-ring | **O(1) + O(degree)** | dict lookup |

## Layerization vs skeletonization

`_peel` peels concentric element rings (the discrete, banded analogue of a
distance field). It is **not** medial-axis or skeleton extraction, those are
distinct constructs documented in [`CONCEPTS.md`](CONCEPTS.md). The ownership
boundary between CHILmesh (post-mesh analysis) and ADMESH (pre-mesh generation)
is set in [#221](https://github.com/domattioli/CHILmesh/issues/221).


---

## Backends

| Backend | Role | Status | How to get it |
|---|---|---|---|
| Python | Reference implementation; every other backend is validated against it | Default | `pip install chilmesh` |
| C++ (half-edge) | Accelerator, bit-identical output, 8.6× to 14.7× on full init | Recommended | `pip install "chilmesh[cpp]"` (prebuilt wheel), or `pip install ./src/chilmesh_cpp` (source build) |
| Rust (quad-edge) | Output-equivalent on all four fixtures (`rust-equivalence` CI); measured 2× to 5× slower than C++ on full init | Frozen, not developed further ([`docs/RUST_EVALUATION.md`](RUST_EVALUATION.md)) | `maturin build` in `src/chilmesh_core`, not recommended |
| MATLAB | Original 2017 implementation | Archived | [`src/@CHILmesh/CHILmesh.m`](../src/@CHILmesh/CHILmesh.m) |

With `CHILMESH_BACKEND` unset, CHILmesh selects the fastest available backend in the order C++, Rust, Python. Force one with `CHILMESH_BACKEND=python|cpp|rust` and inspect the choice with `chilmesh.backend_info()`:

```python
import chilmesh
chilmesh.backend_info()
# After `pip install "chilmesh[cpp]"` or a source build of the C++ extension:
# {'available': ['cpp', 'python'], 'selected': 'cpp',
#  'versions': {'cpp': '0.6.0', 'python': '1.5.0'}}
```

chilmesh 1.5.x requires `chilmesh-cpp>=0.6,<0.7`; both version numbers change together when the extension API changes. Platform coverage is listed in the Installation section of the README.

Parity is gated in CI. The `cpp-equivalence` and `rust-equivalence` jobs build each extension and run the equivalence suite on every pull request to `main` and every push to `main`. After each `chilmesh-cpp` release, a further CI job installs the published wheel from PyPI on Linux, macOS and Windows and runs the same suite.

---

## Why CHILmesh (moved from the README)

CHILmesh began in 2015-2017 at The Ohio State University as the data structure behind QuADMESH+, a layer-peel-driven indirect triangle-to-quadrilateral conversion heuristic for hydrodynamic models (Mattioli, MSc thesis, 2017; see [Citation](../README.md#8-citation)). The mesh generator needed three things no single tool provided: one adjacency model that treats triangles, quadrilaterals and mixed meshes identically; a concentric layer decomposition of the domain, from the boundary inward; and the ADCIRC and SMS file formats that coastal-ocean models read.

The Python library keeps that scope and adds the engineering the MATLAB code lacked. Seven adjacency tables are built once ([`ARCHITECTURE.md`](ARCHITECTURE.md)), giving O(1) edge lookup, O(n) layer peel and O(log n) spatial queries. Every algorithm has one reference implementation in Python. The C++ backend reproduces it bit-for-bit, verified by 76 parametrized equivalence cases in [`tests/test_backend_equivalence.py`](../tests/test_backend_equivalence.py). Terminology is fixed in [`CONCEPTS.md`](CONCEPTS.md): a layer is a concentric band of elements, distinct from a distance field, a medial axis or a skeleton ("vertex", "node" and "point" are interchangeable throughout this documentation).
