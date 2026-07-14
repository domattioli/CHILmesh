# Rust Backend Evaluation

**Date:** 2026-07-14
**Author:** Claude Code (routine session)
**Status:** Evaluation — evidence + recommendation, no backend code changes
**Question:** Does a Rust implementation make sense for CHILmesh over the C++ or
Python backends *in any regard*? And specifically: **could Rust replace Python in
any functionality to improve performance?**

This document records the first *measured* head-to-head of all three backends
(the README/`BENCHMARK.md` Rust perf cells were `tbd` and Rust was excluded from
the cross-language tables — see [#163](https://github.com/domattioli/CHILmesh/issues/163)).
It complements the earlier design-level study
[`.planning/research/quad_edge_feasibility.md`](../.planning/research/quad_edge_feasibility.md)
(2026-05-09), which reasoned about quad-edge *before* the Rust crate existed;
this one measures the crate that shipped in PR #220.

---

## 1. What already exists

Rust is **not** greenfield here — it is a wired-in third backend:

- **Crate:** `src/chilmesh_core/` (PyO3 + `ndarray` + `numpy`), exposing a
  `RustMesh` class with I/O (`fort.14`/`.2dm`), `build_adjacencies`,
  `skeletonize`, `compute_quality`, adjacency getters, vertex/element queries,
  and `add_element`/`remove_element`.
- **Python wrapper:** `src/chilmesh/backends/rust_backend.py` (`RUST_AVAILABLE`,
  `full_init`).
- **Selection:** `chilmesh.backend_info()` reports `rust` when the extension is
  importable and exposes `RustMesh`.
- **CI:** the `rust-equivalence` job (`.github/workflows/python-package.yml`)
  builds the crate with `maturin`, asserts the backend is active, and runs the
  cross-backend equivalence suite.

So the question is not "should we start Rust" but "given the shipped Rust
backend, does it earn its place over C++ or Python?"

---

## 2. Methodology

All three backends were built from source **in this session** and measured
like-for-like — every backend runs the same operation on the same in-memory
`(connectivity, points)` arrays (no fort.14 parse inside the timed region),
using the committed harness `scripts/benchmark.py`.

- **Build:** `maturin build --release` (crate profile: `opt-level=3`, `lto=true`,
  `codegen-units=1`) for Rust; `pip wheel ./src/chilmesh_cpp` (scikit-build-core
  + pybind11, `Release`) for C++.
- **Correctness gate:** `tests/test_backend_equivalence.py` — **76/76 pass**, so
  the Rust build under test is output-correct (layer counts, layer-member sets,
  `bEdgeIDs`, edge ordering, and signed areas all match Python). Perf below is
  therefore measured on a *correct* Rust backend, not a broken one.
- **Machine:** Intel Xeon @ 2.80 GHz (4 cores), 15 GiB RAM, Linux 6.18,
  Python 3.11.15, rustc/cargo 1.94.1, g++ 13.3. Absolute times are
  machine-dependent; the **ratios** are the portable result.
- **Meshes:** the repo's bundled fixtures (`donut`, `annulus`, `Block_O`).
  The continental-scale ENPAC2003 reference mesh is **not** in the repo, so its
  README perf table keeps its `tbd` Rust cells — this evaluation does not
  fabricate numbers for a mesh it did not run.

---

## 3. Measured results

Medians of 5–7 runs. `full-init` = adjacency build + skeletonization (the
apples-to-apples cross-backend operation). `vert-edge` = one `get_vertex_edges`
call, averaged over up to 2000 vertices.

### Full init (adjacency + skeletonization)

| Mesh | Elements | C++ | Rust | Python | Rust vs C++ | Rust vs Python |
|---|---:|---:|---:|---:|---:|---:|
| donut | 276 | ~0 ms | 1 ms | 5 ms | slower | ~5× faster |
| annulus | 580 | 1 ms | 2 ms | 8 ms | ~2× slower | ~4× faster |
| Block_O | 5,214 | 5 ms | 24 ms | 68 ms | **~5× slower** | ~2.8× faster |

`n_layers` parity ✅ across all three on every mesh.

### Vertex-edge lookup (per call)

| Mesh | Elements | C++ | Rust | Python |
|---|---:|---:|---:|---:|
| donut | 276 | 0.64 μs | **50 μs** | 0.41 μs |
| annulus | 580 | 0.81 μs | **103 μs** | 0.44 μs |
| Block_O | 5,214 | 0.64 μs | **954 μs** | 0.51 μs |

The Rust query cost is **100–1800× worse** than C++/Python and **grows with
mesh size** — a scaling defect, not a constant overhead.

---

## 4. Findings

1. **Full init: Rust loses to C++, beats Python.** Rust is ~3–5× faster than
   pure Python but ~2–5× *slower* than C++, and the gap widens with mesh size
   (~5× behind C++ on Block_O). C++ remains the performance backend.

2. **Queries: Rust is catastrophically slow — worse than Python.** Root cause is
   an algorithmic defect in the shipped crate, not the language.
   `RustMesh.get_vertex_edges` (`src/chilmesh_core/lib.rs`) calls
   `adjacency::to_edge2vert(&self.connectivity)` on **every call**, rebuilding
   the entire canonical edge list (O(n_elems), with allocation), then
   `queries::get_vertex_edges` linearly scans **all** edges (O(n_edges)). So each
   lookup is O(n_elems + n_edges); C++ and Python cache a `Vert2Edge` adjacency
   and answer in O(1). This is why the Rust query time tracks mesh size.
   *(Reported as a finding; not fixed in this docs-only evaluation.)*

3. **Correctness is not the issue.** All 76 equivalence tests pass. The problem
   is purely that Rust is slower than the already-shipped C++ backend on the hot
   path, and its query API was never given a cached adjacency.

4. **Prior `tbd`/"excluded" doc state is now resolved.** README's Rust perf cells
   were `tbd`; `BENCHMARK.md` said "Rust is excluded — its skeletonization is
   incomplete (#163)". Both are stale: #163 is closed, skeletonization reaches
   parity, and the missing perf evidence is the table above.

---

## 5. Could Rust *replace Python* in any functionality to improve performance?

The pointed version of the question. A functionality is a candidate only if it is
(a) currently pure-Python with **no** C++ acceleration, (b) actually
perf-critical, and (c) something Rust would plausibly *win* at.

| Functionality | Today | Perf-critical? | Would Rust replacing Python win? | Verdict |
|---|---|---|---|---|
| adjacency / skeletonize / signed-area / vertex queries | Python **+ C++** (+ Rust) | yes | C++ already replaces Python and beats Rust; Rust adds nothing | **No — C++ owns this** |
| FEM smoother (`method='fem'`, direct + iterative) | Python → `scipy.sparse` `spsolve` / MINRES | yes at scale | the solve is already in compiled SuiteSparse/LAPACK; Rust would reimplement a sparse solver — huge effort, unlikely to beat, high risk | **No** |
| angle-based smoother | Python, numpy-vectorized | moderate | numpy batch ops are already near-C; Rust gain marginal | **No** |
| ADMESH warm-start truss | Python, numpy-vectorized | at scale | force loop is vectorized; any win is marginal and better placed in the existing C++ backend | **Marginal — prefer C++** |
| fort.14/.13/.15 + gmsh I/O | Python | no (I/O-bound) | parse time is rarely the bottleneck; Rust I/O exists in the crate but shows no measured win | **No** |
| spatial indexing (point-location / NN — planned Phase 5) | not implemented | would be | greenfield, so no incumbent to beat — a *fair* Rust candidate on merit | **Only open niche, but see below** |
| topology mutations (split/swap/merge/collapse, #94) | Python | moderate | quad-edge `splice` gives O(1) local edits (the feasibility doc's strongest argument); the crate has `add_element`/`remove_element` stubs | **Compiled backend could help — but C++ half-edge is the vehicle** |

**Answer: there is no functionality where replacing Python with Rust is the right
performance move.** Everywhere Rust could help, one of two things is already true:

- **C++ already replaces Python there, faster, with bit-identical output, and is
  the maintained acceleration path.** Adding the same work to the slower Rust
  backend helps no one — you would never *select* Rust as the active backend to
  get one accelerated function while paying a 5× penalty on the shared core.
- **The hot loop is already compiled** (numpy/scipy/BLAS/SuiteSparse). Rewriting
  those mature, well-tested kernels in Rust is high-effort, high-risk, and
  unlikely to beat them.

The only genuinely open niches — Phase-5 spatial indexing and topology-mutation
`splice` — are *greenfield*, where Rust would compete on merit. But even there
the sensible vehicle is the **existing C++ half-edge backend** (already ahead on
the shared core, already bit-identical, already in CI), not a second compiled
backend that is slower and less maintained. The rule of thumb: **if a
pure-Python path ever needs acceleration, extend C++, not Rust.**

---

## 6. Recommendation

**Do not invest further in Rust as a competing backend, and do not migrate any
Python functionality to Rust for performance.**

Rationale:

- C++ already delivers the acceleration (≈8.6–15× over Python on full init on the
  headline meshes; ~5× over Rust here), with bit-identical output, and is the
  documented recommended backend.
- The Rust backend is a second compiled-backend maintenance surface (toolchain,
  wheels, a dedicated CI job, equivalence tests) that is **slower than C++ on the
  hot path and has a query-scaling defect** (§4.2).
- Rust's theoretical advantages (memory safety, `cargo`, fearless concurrency)
  are real in general but **unrealized and unneeded** at CHILmesh's scale — a
  5,000-element mesh full-inits in 24 ms; there is no workload here that a
  compiled, memory-safe, parallel core would rescue that C++ has not already
  covered.

**Concretely:**

- Keep the Rust backend **clearly labeled experimental** (it already is). The
  README/`BENCHMARK.md` updates that accompany this evaluation replace the stale
  `tbd`/"incomplete" language with the measured result and a pointer here.
- If maintenance cost of the second compiled backend is a concern, a **freeze or
  deprecation** of `chilmesh_core` is a reasonable operator call — it removes a
  CI job and a wheel-build path for a backend with no performance niche over C++.
  (Not actioned here; flagged for operator sign-off.)
- If Rust is nonetheless kept active, the one worthwhile, well-scoped fix is to
  **cache `Vert2Edge` in `RustMesh`** so `get_vertex_edges` is O(1) instead of
  rebuilding the edge list per call (§4.2). That would bring queries in line with
  C++/Python — but it does not change the top-line conclusion that C++ is the
  right acceleration path.

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

> Run backend imports from a directory **outside** the repo (or install
> non-editable wheels) so `chilmesh_cpp`/`chilmesh_core` resolve to the built
> extensions rather than the `src/` source dirs on `sys.path` (the namespace-stub
> shadowing noted in #163).
