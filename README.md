<p align="center">
  <video src="videos/quadmesh_logo.mp4" autoplay loop muted playsinline width="680"></video>
</p>

<p align="center">
  <img src="output/readme_pipeline_annulus.gif" alt="CHILmesh pipeline — raw → smoothed → skeletonized" width="720">
</p>

<h1 align="center">CHILmesh</h1>

<p align="center">
  <strong>Fast 2D mesh processing, smoothing, and analysis for triangular, quadrilateral, and mixed-element meshes. Intended for hydrodynamic domains.</strong>
</p>

<p align="center">
  <strong><a href="https://scholar.google.com/citations?user=IBFSkOcAAAAJ&hl=en">Dominik Mattioli</a><sup>1†</sup>, <a href="https://scholar.google.com/citations?user=mYPzjIwAAAAJ&hl=en">Ethan Kubatko</a><sup>2</sup></strong><br>
  <sup>†</sup>Corresponding author | <sup>1</sup>Unaffiliated | <sup>2</sup>Ohio State University (CHIL)
</p>

<p align="center">
  <a href="https://ceg.osu.edu/computational-hydrodynamics-and-informatics-laboratory"><img src="https://img.shields.io/badge/CHIL%20Lab%20@%20OSU-a7b1b7?logo=academia&logoColor=ba0c2f&labelColor=ba0c2f" alt="CHIL Lab @ OSU"></a>
  <a href="https://pypi.org/project/chilmesh/"><img src="https://img.shields.io/github/v/release/domattioli/CHILmesh?include_prereleases&label=PyPI&logo=python&logoColor=white&color=blue&cacheSeconds=300" alt="PyPI"></a>
  <a href="https://github.com/domattioli/CHILmesh/actions/workflows/python-package.yml"><img src="https://img.shields.io/github/actions/workflow/status/domattioli/CHILmesh/python-package.yml?label=Tests&logo=github" alt="Tests"></a>
  <a href="https://doi.org/10.5281/zenodo.20263854"><img src=".github/badges/zenodo-doi.svg" alt="DOI"></a>
  <a href="https://www.mathworks.com/matlabcentral/fileexchange/135632-chilmesh"><img src=".github/badges/matlab-file-exchange.svg" alt="MATLAB File Exchange"></a>
  <a href="https://github.com/domattioli/CHILmesh/blob/main/LICENSE"><img src="https://img.shields.io/badge/License-PolyForm%20NC%20%2B%20No--AI-red.svg?style=flat-square" alt="License"></a>
</p>

> **MATLAB users:** This Python library is the actively-developed successor to the original MATLAB codebase. The original (no longer maintained) is at [`src/@CHILmesh/CHILmesh.m`](src/@CHILmesh/CHILmesh.m) and on [MathWorks](https://www.mathworks.com/matlabcentral/fileexchange/135632-chilmesh/files/src/@CHILmesh/CHILmesh.m).

---

## Why CHILmesh

**The stable backbone for hydrodynamic mesh tooling.** Sibling projects [ADMESH](https://github.com/domattioli/ADMESH), [ADMESH-Domains](https://github.com/domattioli/ADMESH-Domains), and [QuADMesh](https://github.com/domattioli/QuADMesh) build on top of it.

- **Pythonic API** — `from chilmesh import Mesh`; backwards-compatible `CHILmesh` alias preserved.
- **C++ acceleration, bit-identical output** — half-edge extension is **~24× faster than pure Python** on full init, verified bit-for-bit by [36 cross-backend equivalence tests](tests/test_backend_equivalence.py).
- **One interface for all topologies** — triangles, quadrilaterals, and mixed meshes share the same call surface.
- **Stable v1.x API** — sibling projects can pin `chilmesh>=1.0,<2`.

---

## Installation

```bash
pip install chilmesh                        # PyPI
uv pip install chilmesh                     # uv
conda install -c conda-forge chilmesh       # conda-forge (pending)
pip install -e .                            # from source
```

---

## Quick Start

```python
from chilmesh import Mesh

mesh = Mesh.read_from_fort14("ocean.14")
mesh.smooth_mesh(method="fem", acknowledge_change=True)
quality, angles, stats = mesh.elem_quality()
mesh.plot_quality()
```

The legacy `chilmesh.CHILmesh` import is preserved for backward compatibility. Built-in fixtures live at `chilmesh.examples.{annulus, donut, block_o, structured}()`. See [`examples/`](examples/) for runnable scripts.

---

## Features

- **Fast** — full init + quality analysis on a 98,365-element mesh in ~1.7 s (4.6× faster than v0.2.0)
- **Mixed-element** — triangles, quads, and mixed meshes share one API
- **Smoothing** — Balendran direct FEM, Zhou-Shimada angle-based, and ADMESH Spring-Based Truss
- **Analysis** — element quality, interior angles, layer-based skeletonization (medial axis)
- **I/O** — [ADCIRC](https://adcirc.org/) `.fort.14` and [SMS Aquaveo](https://www.aquaveo.com/sms) `.2dm` read/write
- **Spatial queries** — point-in-element, k-nearest vertices, radius search at O(log n)
- **Mesh alterations** — `insert_vertex`, coord moves, advancing-front element addition; full mutation suite tracked in [#94](https://github.com/domattioli/CHILmesh/issues/94)
- **ADMESH-Domains integration** — `from_admesh_domain()` adapter

<p align="center">
  <img src="output/wnat_hagen_showcase.png?v=3" alt="WNAT_Hagen quality plot and distribution">
  <br>
  <sub><em><strong>Figure 1.</strong> Scale demo on WNAT_Hagen (52,774 vertices · 98,365 elements). <code>plot_quality()</code> renders per-element skew quality; <code>plot_quality_histogram()</code> emits the matched-colormap distribution beneath. Reproduce: <code>python scripts/generate_wnat_showcase.py</code>.</em></sub>
</p>

### Performance

Reference workload: WNAT_Hagen (52,774 vertices · 98,365 elements). Median of 3 trials. **v1.0.0 backends are output-equivalent** — the C++ extension produces bit-identical skeletonization layers to Python, verified by [`tests/test_backend_equivalence.py`](tests/test_backend_equivalence.py).

| Metric | v0.1.0 MATLAB ‡ | v0.2.0 Python Port | v0.3.0 Python Optimized | v0.4.0 Rust † | v1.0.0 C++ |
|---|---:|---:|---:|---:|---:|
| Fast init (adj, no skeletonization) | 0.27 s | ~3.9 s | 1.31 s | 0.029 s | 0.036 s |
| Skeletonization only | 0.67 s | ~3.8 s | 0.32 s | 0.20 s | 0.033 s |
| Full init (adj + skeletonization) | 1.04 s | 7.7 s | 1.65 s | 0.23 s | 0.069 s |
| Quality analysis | 12 ms | 6.6 s | 6.4 ms | <1 ms | <1 ms |
| Vertex-edge lookup (per call) | ~2200 μs | ~700 μs | 0.34 μs | 0.02 μs | 0.04 μs |

**C++ is ~24× faster than Python on full init.** ‡ MATLAB v0.1.0 measured under GNU Octave 8.4 — treat as the original-algorithm baseline, not a MATLAB-vs-Octave claim. † Rust (v0.4.0) skeletonization is incomplete ([#163](https://github.com/domattioli/CHILmesh/issues/163)); its full-init figure reflects a partial peel. Full methodology and raw data: [`docs/BENCHMARK.md`](docs/BENCHMARK.md).

### Validation

Python, C++, and the original MATLAB/Octave implementation all produce identical `n_layers` (medial-axis skeletonization) across the ADMESH-Domains catalog, from 557 to 132k vertices. Identical connectivity + points are fed to both implementations; only the layering algorithm is compared.

| Mesh | Vertices | Elements | MATLAB | Python | C++ | Match |
|---|--:|--:|--:|--:|--:|:--:|
| Baranja Hill (ADMESH v2) | 557 | 1,011 | 10 | 10 | 10 | ✅ |
| Baranja Hill | 645 | 1,193 | 12 | 12 | 12 | ✅ |
| Wetting/Drying test | 2,716 | 4,978 | 15 | 15 | 15 | ✅ |
| Lake Erie (refined) | 5,095 | 9,688 | 20 | 20 | 20 | ✅ |
| Lake Erie (5k) | 13,266 | 24,910 | 17 | 17 | 17 | ✅ |
| Delaware Bay | 14,449 | 26,698 | 17 | 17 | 17 | ✅ |
| Delaware Bay (h 100–20000) | 14,449 | 26,697 | 17 | 17 | 17 | ✅ |
| Lake Michigan | 21,981 | 41,887 | 25 | 25 | 25 | ✅ |
| WNAT (Hagen) | 52,774 | 98,365 | 30 | 30 | 30 | ✅ |
| Chesapeake Bay | 83,388 | 160,734 | 55 | 55 | 55 | ✅ |
| Great Lakes | 132,162 | 250,905 | 46 | 46 | 46 | ✅ |

### Smoothing

Three algorithms — each preserves boundary nodes, leaves topology unchanged, and accepts mixed-element meshes.

| Algorithm | API call | Style | Best for |
|---|---|---|---|
| **Balendran direct FEM** | `smooth_mesh(method='fem')` | One-shot sparse solve | General-purpose default; stable on tri/quad/mixed |
| **Zhou-Shimada angle-based** | `smooth_mesh(method='angle-based')` | Iterative, angle-maximising | Difficult mixed meshes where FEM stalls |
| **ADMESH Spring-Based Truss** | `chilmesh.optimize_with_admesh_truss(mesh, sdf, ...)` | Spring/force relaxation against SDF | Quality gains with SDF-respecting boundary nodes |

**References.**
- Balendran (1999). *A direct smoothing method for surface meshes.* Proc. 8th IMR, pp. 189–193.
- Zhou & Shimada (2000). *An angle-based approach to two-dimensional mesh smoothing.* Proc. 9th IMR, pp. 373–384.
- Conroy et al. (2012). *ADMESH: An advanced, automatic unstructured mesh generator for shallow water models.* [doi:10.1007/s10236-012-0574-0](https://doi.org/10.1007/s10236-012-0574-0).

### Backends

`pip install chilmesh` gives you the pure-Python implementation — zero compiled dependencies, runs everywhere, and is the canonical reference every other backend is validated against. The C++ extension is the high-performance opt-in: same algorithms, bit-identical output, ~24× faster on full init.

| Language | Role | How to get it |
|---|---|---|
| **Python** | Reference implementation — the default | `pip install chilmesh` |
| **C++** | High-performance backend (half-edge) — bit-identical output | `pip install ./src/chilmesh_cpp` (build from source) |
| Rust | Experimental (quad-edge); skeletonization is incomplete — see [#163](https://github.com/domattioli/CHILmesh/issues/163) | source build, not recommended yet |
| MATLAB | Original 2017 implementation, archived & unmaintained | [`src/@CHILmesh/CHILmesh.m`](src/@CHILmesh/CHILmesh.m) |

```python
import chilmesh

chilmesh.backend_info()
# {'available': ['cpp', 'python'],
#  'selected': 'cpp',
#  'versions': {'cpp': '0.6.0.dev0', 'python': '1.1.0'}}
```

Force a specific backend with `CHILMESH_BACKEND` (`python` or `cpp`). When unset, the fastest available is picked. Pre-built binary wheels (`manylinux` / `macOS` / `Windows`) via `cibuildwheel` are planned — see [`docs/`](docs/) for build-from-source instructions.

### Examples

```bash
python examples/01_quickstart.py        # load, stats, plot
python examples/02_fort14_roundtrip.py  # fort.14 read/write
python examples/03_smoothing.py         # angle-based smoother
python examples/04_spatial_queries.py   # find_element, radius search, k-nearest
```

### CLI

```bash
chilmesh info mesh.fort.14                                      # stats
chilmesh convert mesh.2dm mesh.fort.14                         # format conversion
chilmesh smooth mesh.fort.14 -o out.fort.14 --method fem       # smooth in-place
chilmesh plot mesh.fort.14 -o mesh.png --quality               # render
```

Also available as `python -m chilmesh`. Each subcommand has `--help`.

---

## Ecosystem

CHILmesh is the core engine for the ADCIRC mesh ecosystem. Sibling projects build on it; CHILmesh depends on none of them.

| Repo | Role |
|---|---|
| [ADMESH](https://github.com/domattioli/ADMESH) | Unstructured triangle mesh generator; consumes CHILmesh for adjacency, smoothing, and quality analysis |
| [ADMESH-Domains](https://github.com/domattioli/ADMESH-Domains) | Curated ADCIRC mesh registry; `Mesh.from_admesh_domain()` reads from it directly |
| [QuADMesh](https://github.com/domattioli/QuADMesh) | Quad mesh generator (MATLAB → Python port, in progress); CHILmesh data structure descends from the original QuADMesh+ |
| [MADMESHing](https://github.com/domattioli/MADMESHing) | Benchmark harness comparing ADMESH triangulation vs quad generators; uses CHILmesh for quality analysis |

*[DomI](https://github.com/domattioli/DomI) provides dev-session skills and governance infrastructure for all repos.*

---

## Status & Roadmap

- **Shipped (v1.0.0)**: C++ half-edge backend (~24× faster on full init); bit-identical output verified; 36 cross-backend equivalence tests; fort.14 + .2dm I/O; mixed-element support.
- **In flight**: Pre-built binary wheels (cibuildwheel, manylinux/macOS/Windows) · Rust skeletonization completion ([#163](https://github.com/domattioli/CHILmesh/issues/163)) · Full mutation suite ([#94](https://github.com/domattioli/CHILmesh/issues/94))
- **Next**: conda-forge packaging · mkdocs API site · advancing-front element mutation

Open issues: [github.com/domattioli/CHILmesh/issues](https://github.com/domattioli/CHILmesh/issues)

---

## Documentation

- [`docs/API.md`](docs/API.md) — full API reference
- [`docs/BENCHMARK.md`](docs/BENCHMARK.md) — benchmark methodology and raw data
- [`TESTING.md`](TESTING.md) — test guide (pytest markers, local commands)
- [`examples/`](examples/) — runnable scripts (quickstart, fort.14 round-trip, smoothing, spatial queries)

---

## Citation

CHILmesh originated in MATLAB as the data structure backing a skeletonization-driven indirect tri-to-quad conversion heuristic (Mattioli, OSU MSc [thesis](https://github.com/user-attachments/files/19727573/QuADMESH__Thesis_Doc.pdf), 2017). This Python library is the actively-developed successor.

```bibtex
@software{mattioli_chilmesh,
  author    = {Mattioli, Dominik O. and Kubatko, Ethan J.},
  title     = {{CHILmesh}: a fast 2D mesh library for triangular,
               quadrilateral, and mixed-element grids},
  year      = {2026},
  publisher = {Zenodo},
  version   = {1.0.0},
  doi       = {10.5281/zenodo.20263854},
  url       = {https://github.com/domattioli/CHILmesh}
}
```

**MATLAB source (Mattioli, 2017).** [Read thesis (PDF)](https://github.com/user-attachments/files/19727573/QuADMESH__Thesis_Doc.pdf)

```bibtex
@mastersthesis{mattioli2017quadmesh,
  author = {Mattioli, Dominik O.},
  title  = {{QuADMESH+}: A Quadrangular ADvanced Mesh Generator
            for Hydrodynamic Models},
  school = {The Ohio State University},
  year   = {2017},
  url    = {http://rave.ohiolink.edu/etdc/view?acc_num=osu1500627779532088}
}
```

---

## Contributing

Issues and PRs welcome at [github.com/domattioli/CHILmesh](https://github.com/domattioli/CHILmesh). Run `pytest -v` before opening a PR — see [`TESTING.md`](TESTING.md).

---

## License

**Noncommercial / research use only.** Licensed under the PolyForm Noncommercial
License 1.0.0 **with an additional No-AI/ML-training restriction** — see
[LICENSE](LICENSE) and [AI-USAGE.md](AI-USAGE.md). No commercial use and no use
as AI/ML training data without a separate written license. Commercial or
AI-training licenses: domburner@duck.com
