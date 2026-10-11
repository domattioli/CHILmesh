# Rust Backend Evaluation

Rust does not beat C++ on any measured operation, and no Python functionality should move to Rust for speed. Full init on Block_O takes 24 ms in Rust, 5 ms in C++ and 68 ms in Python. The Rust backend is frozen (2026-07-14): it stays output-equivalent and gets no further development.

**Date:** 2026-07-14
**Backend status:** ❄️ **FROZEN** (2026-07-14, operator-directed). The backend is kept and output-equivalent, and it is not developed further. See `src/chilmesh_core/STATUS.md`.
**Status:** Evaluation with evidence and a recommendation. The O(1) query fix was applied; the backend is otherwise frozen.
**Question:** Does a Rust implementation make sense for CHILmesh over the C++ or Python backends *in any regard*? Could Rust replace Python in any functionality to improve performance?

This document records the first *measured* head-to-head of all three backends. The README and `BENCHMARK.md` Rust perf cells were `tbd`, and Rust was excluded from the cross-language tables ([#163](https://github.com/domattioli/CHILmesh/issues/163)).

It complements the earlier design-level study [`.planning/research/quad_edge_feasibility.md`](../.planning/research/quad_edge_feasibility.md) (2026-05-09). That study reasoned about quad-edge *before* the Rust crate existed. This document measures the crate that shipped in PR #220.

---

## 1. What already exists

Rust is a wired-in third backend, not a greenfield idea.

- **Crate:** `src/chilmesh_core/` (PyO3 + `ndarray` + `numpy`). It exposes a `RustMesh` class with I/O (`fort.14`/`.2dm`), `build_adjacencies`, `skeletonize`, `compute_quality`, adjacency getters, vertex/element queries, and `add_element`/`remove_element`.
- **Python wrapper:** `src/chilmesh/backends/rust_backend.py` (`RUST_AVAILABLE`, `full_init`).
- **Selection:** `chilmesh.backend_info()` reports `rust` when the extension is importable and exposes `RustMesh`.
- **CI:** the `rust-equivalence` job (`.github/workflows/python-package.yml`) builds the crate with `maturin`, asserts the backend is active, and runs the cross-backend equivalence suite.

The question is whether the shipped Rust backend earns its place over C++ or Python.

---

## 2. Methodology

The harness built all three backends from source and measured them like-for-like. Every backend runs the same operation on the same in-memory `(connectivity, points)` arrays. No fort.14 parse happens inside the timed region. The harness is the committed `scripts/benchmark.py`.

**Terminology.** The layer decomposition is the **layer peel** (*onion peeling*). Per the #187 lexicon, "skeletonization" is deprecated for this operation and reserved for a future medial-axis op (#223). The backend extension method is still named `skeletonize()` pending a deferred cross-language backend rename. That symbol name persists in code, and prose says "peel".

- **Build:** Rust uses `maturin build --release` (crate profile: `opt-level=3`, `lto=true`, `codegen-units=1`). C++ uses `pip wheel ./src/chilmesh_cpp` (scikit-build-core + pybind11, `Release`).
- **Correctness gate:** `tests/test_backend_equivalence.py`, **76/76 pass**. The Rust build under test is output-correct: layer counts, layer-member sets, `bEdgeIDs`, edge ordering, and signed areas all match Python. The perf numbers below therefore come from a *correct* Rust backend.
- **Machine:** Intel Xeon @ 2.80 GHz (4 cores), 15 GiB RAM, Linux 6.18, Python 3.11.15, rustc/cargo 1.94.1, g++ 13.3. Absolute times depend on the machine. The **ratios** are the portable result.
- **Meshes:** the repo's bundled fixtures (`donut`, `annulus`, `Block_O`). The continental-scale ENPAC2003 reference mesh is not in the repo. Its README perf table keeps its `tbd` Rust cells, and this evaluation reports no numbers for a mesh it did not run.

---

## 3. Measured results

Each value is the median of 5–7 runs.

- `full-init` is adjacency build plus layer peel, the apples-to-apples cross-backend operation.
- `vert-edge` is one `get_vertex_edges` call, averaged over up to 2000 vertices.

### Full init (adjacency + layer peel)

| Mesh | Elements | C++ | Rust | Python | Rust vs C++ | Rust vs Python |
|---|---:|---:|---:|---:|---:|---:|
| donut | 276 | ~0 ms | 1 ms | 5 ms | slower | ~5× faster |
| annulus | 580 | 1 ms | 2 ms | 8 ms | ~2× slower | ~4× faster |
| Block_O | 5,214 | 5 ms | 24 ms | 68 ms | **~5× slower** | ~2.8× faster |

`n_layers` parity ✅ across all three on every mesh.

### Vertex-edge lookup (per call)

| Mesh | Elements | C++ | Rust (before fix) | Rust (after fix) | Python |
|---|---:|---:|---:|---:|---:|
| donut | 276 | 0.64 μs | 50 μs | **0.53 μs** | 0.41 μs |
| annulus | 580 | 0.81 μs | 103 μs | **0.40 μs** | 0.44 μs |
| Block_O | 5,214 | 0.64 μs | 954 μs | **0.32 μs** | 0.51 μs |

Before the fix (§4.2), the Rust query cost was **100–1800× worse** than C++ and Python. It also **grew with mesh size**, so the cause was a scaling defect and not a constant overhead.

After the fix, the query is O(1). It is on par with, or slightly faster than, C++ and Python. On Block_O it dropped from 954 μs to 0.32 μs (~3000×). Equivalence still holds (76/76 tests pass).

---

## 4. Findings

### 4.1 Full init: Rust loses to C++ and beats Python

Rust is ~3–5× faster than pure Python. It is ~2–5× *slower* than C++. The gap widens with mesh size, and Rust is ~5× behind C++ on Block_O. C++ remains the performance backend.

### 4.2 Queries were slow because of an algorithmic defect, now fixed

The root cause was in the shipped crate, not in the language.

1. `RustMesh.get_vertex_edges` (`src/chilmesh_core/lib.rs`) called `adjacency::to_edge2vert(&self.connectivity)` on **every call**.
2. That call rebuilt the entire canonical edge list (O(n_elems), with allocation).
3. The method then scanned **all** edges linearly (O(n_edges)).
4. Each lookup therefore cost O(n_elems + n_edges).
5. C++ and Python cache a `Vert2Edge` adjacency and answer in O(1). This is why the Rust query time tracked mesh size.

The fix has two parts:

1. `build_adjacencies` now precomputes the vertex→edge index once, from the same `to_edge2vert` source, so output stays bit-identical.
2. `get_vertex_edges` returns the cached row in O(1), and `set_connectivity` invalidates the cache.

The result: Block_O dropped from 954 μs to 0.32 μs (~3000×). That is on par with, or faster than, C++ and Python. All 76 equivalence tests still pass (see the "after fix" column in §3).

### 4.3 Correctness is not the issue

All 76 equivalence tests pass. Rust is slower than the shipped C++ backend on the hot path, and its query API had no cached adjacency.

### 4.4 The `tbd` and "excluded" doc state is resolved

The README Rust perf cells were `tbd`. `BENCHMARK.md` said "Rust is excluded, its skeletonization is incomplete (#163)". Both statements are stale. #163 is closed, the layer peel reaches parity, and the table in §3 supplies the missing perf evidence.

### 4.5 The backends differ in internal topology, not in output

- **Python** is the reference. It uses flat numpy arrays (`Edge2Vert`, `Elem2Edge`, `Edge2Elem`), dict adjacencies (`Vert2Edge`, `Vert2Elem`), and an `EdgeMap` hash.
- **C++** uses a **half-edge (DCEL)** structure (`src/chilmesh_cpp/src/halfedge.*`).
- **Rust** uses a **quad-edge** structure (Guibas–Stolfi 1985; four directed edges per undirected edge).

All three produce **bit-identical** layers, adjacency tables, and signed areas. `tests/test_backend_equivalence.py` guards that equivalence.

The quad-edge choice was an explicit experiment (`.planning/008-DECISION.md`). It tested whether a different topology structure would scale better than the flat and half-edge representations. The measurements show that it does **not** beat the C++ half-edge backend, so the data-structure bet did not pay off on performance.

That result and the doubled compiled-backend maintenance surface are **why the backend is frozen**. The backend is bit-exact, so correctness is not the reason. It is a slower second way to compute what C++ already computes faster.

---

## 5. Could Rust *replace Python* in any functionality to improve performance?

No. A functionality is a candidate only if it meets three conditions:

1. It is currently pure Python with **no** C++ acceleration.
2. It is actually perf-critical.
3. Rust would plausibly *win* at it.

| Functionality | Today | Perf-critical? | Would Rust replacing Python win? | Verdict |
|---|---|---|---|---|
| adjacency / layer-peel / signed-area / vertex queries | Python **+ C++** (+ Rust) | yes | C++ already replaces Python and beats Rust; Rust adds nothing | **No; C++ owns this** |
| FEM smoother (`method='fem'`, direct + iterative) | Python → `scipy.sparse` `spsolve` / MINRES | yes at scale | the solve is already in compiled SuiteSparse/LAPACK; Rust would reimplement a sparse solver; huge effort, unlikely to beat, high risk | **No** |
| angle-based smoother | Python, numpy-vectorized | moderate | numpy batch ops are already near-C; Rust gain marginal | **No** |
| ADMESH warm-start truss | Python, numpy-vectorized | at scale | force loop is vectorized; any win is marginal and better placed in the existing C++ backend | **Marginal; prefer C++** |
| fort.14/.13/.15 + gmsh I/O | Python | no (I/O-bound) | parse time is rarely the bottleneck; Rust I/O exists in the crate but shows no measured win | **No** |
| spatial indexing (point-location / NN; planned Phase 5) | not implemented | would be | greenfield, so no incumbent to beat; a *fair* Rust candidate on merit | **Only open niche, but see below** |
| topology mutations (split/swap/merge/collapse, #94) | Python | moderate | quad-edge `splice` gives O(1) local edits (the feasibility doc's strongest argument); the crate has `add_element`/`remove_element` stubs | **Compiled backend could help; but C++ half-edge is the vehicle** |

Everywhere Rust could help, one of two things is already true:

- **C++ already replaces Python there.** C++ is faster, produces bit-identical output, and is the maintained acceleration path. Nobody would *select* Rust as the active backend to get one accelerated function while paying a 5× penalty on the shared core.
- **The hot loop is already compiled** (numpy/scipy/BLAS/SuiteSparse). Rewriting those mature, well-tested kernels in Rust is high-effort and high-risk, and it is unlikely to beat them.

The only open niches are Phase-5 spatial indexing and topology-mutation `splice`. Both are *greenfield*, so Rust would compete on merit there. Even so, the sensible vehicle is the **existing C++ half-edge backend**. It is already ahead on the shared core, already bit-identical, and already in CI. A second compiled backend would be slower and less maintained.

The rule of thumb: **if a pure-Python path ever needs acceleration, extend C++, not Rust.**

### 5.1 Default backend and opt-in: can Rust replace Python as the default?

Users opt in to C++ today. The question is whether Rust could become the *default*, so users get speed without opting in. The answer turns on distribution, not language.

- **Both C++ and Rust are opt-in source builds today.** A plain `pip install chilmesh` from PyPI ships **pure Python only**, with no compiled extension ([#229](https://github.com/domattioli/CHILmesh/issues/229)). C++ needs a C++ toolchain and CMake (`pip install ./src/chilmesh_cpp`). Rust needs a Rust toolchain (`maturin build …`). Both need a compiler and a build step, so Rust is not a lighter opt-in than C++.
- **Auto-selection already prefers Rust over Python when Rust is built.** `backend_info()` (`src/chilmesh/__init__.py`) picks the fastest *available* backend in the order C++ → Rust → Python, and it honors `CHILMESH_BACKEND`. A user who has not opted into C++ already gets Rust *if the Rust extension is present*. That user gets pure Python only because the wheel ships nothing compiled.
- **A zero-opt-in compiled default needs prebuilt binary wheels** (manylinux/macOS/Windows via `cibuildwheel`/`maturin`). That work is already planned for **C++** ([#229](https://github.com/domattioli/CHILmesh/issues/229)). A prebuilt wheel should carry the **faster** backend, and C++ is ~5× faster than Rust on full init. Defaulting to Rust would give up that 5×.
- **Rust has one distribution advantage.** maturin builds an `abi3` wheel, so a *single* wheel per platform works across Python 3.8+. The pybind11 C++ extension needs a wheel per Python version. Fewer wheels is a real packaging convenience, but it does not outweigh a 5× slowdown on the hot path. It is at best a stopgap if the C++ wheel matrix stalls.

**Net:**

1. Keep **pure Python as the universal default and fallback**. It has zero dependencies, runs everywhere, and is the reference every backend is validated against.
2. Pursue **prebuilt C++ wheels** as the auto-selected accelerator (existing plan).
3. Do not make Rust the default backend.
4. Do not migrate any Python computation to Rust for speed.

---

## 6. Recommendation

**Do not invest further in Rust as a competing backend, and do not migrate any Python functionality to Rust for performance.**

Rationale:

- C++ already delivers the acceleration: ≈8.6–15× over Python on full init on the headline meshes, and ~5× over Rust here. Its output is bit-identical, and it is the documented recommended backend.
- The Rust backend is a second compiled-backend maintenance surface (toolchain, wheels, a dedicated CI job, equivalence tests). It is **slower than C++ on the hot path** and had a query-scaling defect (§4.2).
- Rust's theoretical advantages (memory safety, `cargo`, fearless concurrency) are real in general. They are **unrealized and unneeded** at CHILmesh's scale: a 5,000-element mesh full-inits in 24 ms. No workload here needs a compiled, memory-safe, parallel core that C++ does not already cover.

**Concretely:**

1. **The Rust backend is FROZEN** (operator-directed, 2026-07-14). It is kept and output-equivalent, but not developed further and not kept in lockstep as the Python API grows. The status banner is `src/chilmesh_core/STATUS.md`.
2. The crate, its `rust-equivalence` CI job, the equivalence tests, and the benchmark data all stay. No new feature work lands on the backend. Freezing, rather than deprecating, retains the reference and parity value.
3. **The acceleration path is prebuilt C++ binary wheels**, not a second backend. A plain `pip install` then gives C++ speed with no toolchain. The plan is [`docs/dev/PREBUILT_WHEELS_PLAN.md`](dev/PREBUILT_WHEELS_PLAN.md) (#229).
4. **The Rust query fix is applied.** `RustMesh` now caches `Vert2Edge`, so `get_vertex_edges` is O(1) and no longer rebuilds the edge list per call (§4.2). Rust queries now match C++ and Python.
5. The fix removes a real footgun. `backend_info()` auto-selects Rust over Python when Rust is built and C++ is not (`src/chilmesh/__init__.py`). Before the fix, a Rust-only build ran query-heavy workloads far slower than pure Python.
6. The fix does not change the top-line conclusion. Rust full-init is still ~5× behind C++, so C++ remains the acceleration path.

---

## 7. Reproduce

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]" maturin

# C++ backend
pip wheel ./src/chilmesh_cpp -w wheelhouse --no-deps
pip install --force-reinstall --no-deps wheelhouse/chilmesh_cpp-*.whl

# Rust backend
maturin build --release --manifest-path src/chilmesh_core/Cargo.toml --out wheelhouse
pip install --force-reinstall --no-deps wheelhouse/chilmesh_core-*.whl

# confirm all three, then measure
python -c "import chilmesh; print(chilmesh.backend_info())"   # -> ['cpp','rust','python']
python -m pytest tests/test_backend_equivalence.py -q         # 76 passed
CHILMESH_RUN_BENCH=1 python scripts/benchmark.py --mesh src/chilmesh/data/Block_O.14 --repeats 5
```

Run backend imports from a directory **outside** the repo, or install non-editable wheels. Then `chilmesh_cpp` and `chilmesh_core` resolve to the built extensions and not to the `src/` source dirs on `sys.path` (the namespace-stub shadowing noted in #163).
