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
  <a href="https://pypi.org/project/chilmesh/"><img src="https://img.shields.io/pypi/v/chilmesh?label=PyPI&logo=python&logoColor=white&cacheSeconds=300" alt="PyPI"></a>
  <a href="https://github.com/domattioli/CHILmesh/actions/workflows/python-package.yml"><img src="https://img.shields.io/github/actions/workflow/status/domattioli/CHILmesh/python-package.yml?label=Tests&logo=github" alt="Tests"></a>
  <a href="https://doi.org/10.5281/zenodo.20263854"><img src=".github/badges/zenodo-doi.svg" alt="DOI"></a>
  <a href="https://www.mathworks.com/matlabcentral/fileexchange/135632-chilmesh"><img src=".github/badges/matlab-file-exchange.svg" alt="MATLAB File Exchange"></a>
  <a href="https://github.com/domattioli/CHILmesh/blob/main/LICENSE"><img src="https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square" alt="License"></a>
</p>

> **MATLAB users:** This Python library is the actively-developed successor to the original MATLAB codebase. The original (no longer maintained) is at [`src/@CHILmesh/CHILmesh.m`](src/@CHILmesh/CHILmesh.m) and on [MathWorks](https://www.mathworks.com/matlabcentral/fileexchange/135632-chilmesh/files/src/@CHILmesh/CHILmesh.m).

---

## Table of Contents

- [Quick Start](#quick-start)
- [Gallery](#gallery)
- [Features](#features)
- [Installation](#installation)
- [Performance](#performance)
- [API Overview](#api-overview)
- [Mesh Smoothing](#mesh-smoothing)
- [Examples](#examples)
- [CLI](#cli)
- [Downstream Projects](#downstream-projects)
- [Contributing](#contributing)
- [Citation](#citation)

---

## Why CHILmesh

CHILmesh is the **stable backbone** for hydrodynamic mesh tooling — used by [ADMESH](https://github.com/domattioli/ADMESH), [MADMESHR](https://github.com/domattioli/MADMESHR), and [ADMESH-Domains](https://github.com/domattioli/ADMESH-Domains). v1.0.0 ships a Pythonic public API (`from chilmesh import Mesh`) backed by a half-edge C++ extension that is **66× faster than pure Python** on skeletonization while producing bit-identical layer output (verified by [36 cross-backend equivalence tests](tests/test_backend_equivalence.py)). Triangles, quadrilaterals, and mixed meshes share one interface; downstream projects can rely on a stable v1.x API.

---

## Quick Start

```bash
pip install chilmesh
```

```python
from chilmesh import Mesh

mesh = Mesh.read_from_fort14("ocean.14")
mesh.smooth_mesh(method="fem", acknowledge_change=True)
quality, angles, stats = mesh.elem_quality()
mesh.plot_quality()
```

The legacy `chilmesh.CHILmesh` import is preserved for backward compatibility. Built-in fixtures live at `chilmesh.examples.{annulus, donut, block_o, structured}()`. See [`examples/`](examples/) for runnable scripts.

---

## Gallery

<p align="center">
  <img src="output/wnat_hagen_showcase.png?v=3" alt="WNAT_Hagen quality plot and distribution">
  <br>
  <sub><em><strong>Figure 1.</strong> Scale demo on WNAT_Hagen (52,774 vertices · 98,365 elements). <code>plot_quality()</code> renders per-element skew quality; <code>plot_quality_histogram()</code> emits the matched-colormap distribution beneath. Reproduce: <code>python scripts/generate_wnat_showcase.py</code>.</em></sub>
</p>

<p align="center">
  <img src="output/mixed_mesh_showcase.png?v=2" alt="Mixed-element mesh: wireframe, layers, quality">
  <br>
  <sub><em><strong>Figure 2.</strong> Mixed-element pipeline — wireframe, skeletonization, and per-element quality on one tri+quad mesh. Reproduce: <code>python scripts/generate_mixed_truss_demo.py</code>.</em></sub>
</p>

<p align="center">
  <img src="output/annulus_quickstart.png?v=7" alt="Skeletonization + quality plotting across three smoothing states">
  <br>
  <sub><em><strong>Figure 3.</strong> <code>plot_layer()</code> and <code>plot_quality()</code> tracking skeletonization and quality across raw → truss → FEM smoothing. Reproduce: <code>python scripts/generate_3row_admesh.py</code>.</em></sub>
</p>

<p align="center">
  <img src="output/readme_pipeline_annulus.gif" alt="Animated 4-stage CHILmesh annulus pipeline">
  <br>
  <sub><em><strong>Figure 4.</strong> Figure 3 as a Manim animation. Higher-fidelity 1080p at <a href="output/readme_pipeline_annulus.mp4"><code>output/readme_pipeline_annulus.mp4</code></a>.</em></sub>
</p>

---

## Features

- **Fast** — full init + quality analysis on a 98,365-element mesh in ~3.3 s (4.3× faster than v0.2.0)
- **Mixed-element** — triangles, quads, and mixed meshes share one API
- **Smoothing** — Balendran direct FEM, Zhou-Shimada angle-based, and ADMESH Spring-Based Truss
- **Analysis** — element quality, interior angles, layer-based skeletonization (medial axis)
- **I/O** — [ADCIRC](https://adcirc.org/) `.fort.14` and [SMS Aquaveo](https://www.aquaveo.com/sms) `.2dm` read/write
- **Spatial queries** — point-in-element, k-nearest vertices, radius search at O(log n)
- **Mesh alterations** — `insert_vertex`, coord moves, advancing-front element addition; full mutation suite tracked in [#94](https://github.com/domattioli/CHILmesh/issues/94)
- **ADMESH-Domains integration** — `from_admesh_domain()` adapter

---

## Installation

```bash
pip install chilmesh                        # PyPI
uv pip install chilmesh                     # uv
conda install -c conda-forge chilmesh       # conda-forge (pending)
pip install -e .                            # from source
```

---

## Performance

Reference workload: WNAT_Hagen (52,774 vertices · 98,365 elements). Median of 3 trials. **v1.0.0 backends are output-equivalent** — the C++ extension produces bit-identical skeletonization layers to Python, verified by [`tests/test_backend_equivalence.py`](tests/test_backend_equivalence.py).

| Metric | v0.2.0 MATLAB ※ | v1.0.0 Python | v1.0.0 Rust † | v1.0.0 C++ |
|---|---:|---:|---:|---:|
| Fast init (adj, no skeletonization) | ~3.9 s | 1.01 s | 0.029 s | **0.036 s** |
| Skeletonization only | ~3.8 s | 2.20 s | 0.20 s | **0.033 s** |
| Full init (adj + skeletonization) | 7.7 s | 3.21 s | 0.23 s | **0.069 s** |
| Quality analysis | 6.6 s | 57 ms | <1 ms | **<1 ms** |
| Vertex-edge lookup (per call) | ~700 μs | **0.08 μs** | 0.02 μs | 0.04 μs |

**C++ is 46× faster than Python on full init and 66× faster on skeletonization** — at the same logical output. The Python implementation remains the canonical reference; C++ is opt-in via direct `chilmesh_cpp` import, Rust via `chilmesh_core`.

※ MATLAB v0.2.0 = direct Python port of original MATLAB implementation ([Mattioli, OSU MSc thesis, 2017](https://github.com/user-attachments/files/19727573/QuADMESH__Thesis_Doc.pdf)).  
† Rust fast init includes fort.14 file I/O; Python and C++ receive raw arrays.

Full methodology and raw data: [`docs/BENCHMARK.md`](docs/BENCHMARK.md).

---

## Backends

CHILmesh v1.0.0 ships pure-Python by default. The C++ half-edge extension (`chilmesh_cpp`) and the Rust quad-edge extension (`chilmesh_core`) are optional accelerators that produce identical mesh topology and layer output.

```python
import chilmesh

chilmesh.backend_info()
# {'available': ['cpp', 'rust', 'python'],
#  'selected': 'cpp',
#  'versions': {'cpp': '0.6.0.dev0', 'rust': '0.5.0.dev0', 'python': '1.0.0'}}
```

Force a specific backend with the `CHILMESH_BACKEND` environment variable (`python`, `cpp`, or `rust`). When in doubt, leave it unset — defaults pick the fastest available.

**Building the C++ extension** (from source):

```bash
cd src/chilmesh_cpp
pip install .                                       # uses scikit-build-core + pybind11
```

Pre-built binary wheels for `manylinux` / `macOS` / `Windows` arrive in v1.1.0 via `cibuildwheel`.

---

## API Overview

```python
from chilmesh import Mesh, examples

# Load
mesh = examples.annulus()
mesh = Mesh.read_from_fort14('mesh.14')
mesh = Mesh.read_from_2dm('mesh.2dm')

# Smooth, analyse, visualise
mesh.smooth_mesh(method='fem', acknowledge_change=True)
quality, angles, stats = mesh.elem_quality()
mesh.plot()             # wireframe
mesh.plot_quality()     # per-element quality
mesh.plot_layer()       # skeletonization layers

# Skeletonization output
layers = mesh.Layers    # {'OE', 'IE', 'OV', 'IV', 'bEdgeIDs'} per layer

# Spatial queries
elem_id = mesh.find_element([0.5, 0.0])
neighbors = mesh.nearest_vertices([0.5, 0.0], k=5)
in_radius = mesh.find_elements_in_radius([0.5, 0.0], radius=0.2)
```

Full reference: [`docs/API.md`](docs/API.md).

---

## Mesh Smoothing

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

---

## Examples

```bash
python examples/01_quickstart.py        # load, stats, plot
python examples/02_fort14_roundtrip.py  # fort.14 read/write
python examples/03_smoothing.py         # angle-based smoother
python examples/04_spatial_queries.py   # find_element, radius search, k-nearest
```

---

## CLI

```bash
chilmesh info mesh.fort.14                                      # stats
chilmesh convert mesh.2dm mesh.fort.14                         # format conversion
chilmesh smooth mesh.fort.14 -o out.fort.14 --method fem       # smooth in-place
chilmesh plot mesh.fort.14 -o mesh.png --quality               # render
```

Also available as `python -m chilmesh`. Each subcommand has `--help`.

---

## Downstream Projects

| Repo | Description |
|---|---|
| [**ADMESH**](https://github.com/domattioli/ADMESH) | Optimized 2D triangular mesh generation for hydrodynamic domains |
| [**MADMESHR**](https://github.com/domattioli/MADMESHR) | AI-based quad- and mixed-element generation for hydrodynamic domains |
| [**ADMESH-Domains**](https://github.com/domattioli/ADMESH-Domains) | Mesh catalog for hydrodynamic domains |

---

## Contributing

Issues and PRs welcome at [github.com/domattioli/CHILmesh](https://github.com/domattioli/CHILmesh). Run `pytest -v` before opening a PR — see [`TESTING.md`](TESTING.md).

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
  version   = {0.4.1},
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

## License

MIT — see [LICENSE](LICENSE).
