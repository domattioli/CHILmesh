<h1 align="center">
  CHILmesh
</h1>

<p align="center">
  <strong>Fast 2D mesh processing, smoothing, and analysis for triangular, quadrilateral, and mixed-element meshes.</strong>
</p>

<p align="center">
  <strong><a href="https://scholar.google.com/citations?user=IBFSkOcAAAAJ&hl=en">Dominik Mattioli</a><sup>1†</sup>, <a href="https://scholar.google.com/citations?user=mYPzjIwAAAAJ&hl=en">Ethan Kubatko</a><sup>2</sup></strong><br>
  <sup>†</sup>Corresponding author | <sup>1</sup>Penn State University | <sup>2</sup>Ohio State University (CHIL)
</p>

<p align="center">
  <a href="https://ceg.osu.edu/computational-hydrodynamics-and-informatics-laboratory"><img src="https://img.shields.io/badge/CHIL%20Lab%20@%20OSU-a7b1b7?logo=academia&logoColor=ba0c2f&labelColor=ba0c2f" alt="CHIL Lab @ OSU"></a>
  <a href="https://github.com/domattioli/ADMESH"><img src="https://img.shields.io/badge/OSU_CHIL-ADMESH-66bb33?logo=github&logoColor=ba0c2f&labelColor=ffffff" alt="ADMESH"></a>
  <a href="https://pypi.org/project/chilmesh/"><img src="https://img.shields.io/pypi/v/chilmesh?label=PyPI&logo=python&logoColor=white" alt="PyPI"></a>
  <a href="https://github.com/domattioli/CHILmesh/actions/workflows/python-package.yml"><img src="https://img.shields.io/github/actions/workflow/status/domattioli/CHILmesh/python-package.yml?label=Tests&logo=github" alt="Tests"></a>
  <a href="https://zenodo.org/badge/latestdoi/693749657"><img src="https://zenodo.org/badge/693749657.svg" alt="DOI"></a>
  <a href="https://github.com/domattioli/CHILmesh/blob/main/LICENSE"><img src="https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square" alt="License"></a>
</p>

> **Note for MATLAB users**: This Python implementation is the actively-developed successor to the original MATLAB codebase. It is **still in development** and the API may evolve. The original MATLAB code (no longer maintained) remains available for reference at [src/@CHILmesh/CHILmesh.m](src/@CHILmesh/CHILmesh.m).

---

## Table of Contents

- [Quick Start](#quick-start)
- [Gallery](#gallery)
- [Features](#features)
- [Installation](#installation)
- [Performance](#performance-v030)
- [API Overview](#api-overview)
- [Mesh Smoothing](#mesh-smoothing)
- [Examples](#examples)
- [CLI](#cli)
- [Contributing](#contributing)
- [Citation](#citation)

---

## Quick Start

```bash
pip install chilmesh
```

```python
import chilmesh
import matplotlib.pyplot as plt

mesh = chilmesh.examples.annulus()
mesh.smooth_mesh(method='fem', acknowledge_change=True)

quality, angles, stats = mesh.elem_quality()
mesh.plot_quality()
plt.show()
```

See [`examples/`](examples/) for more runnable scripts.

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
  <sub><em><strong>Figure 2.</strong> Mixed-element pipeline — wireframe, skeletonization, and per-element quality on one tri+quad mesh, all via the standard API. Reproduce: <code>python scripts/generate_mixed_truss_demo.py</code>.</em></sub>
</p>

<p align="center">
  <img src="output/annulus_quickstart.png?v=7" alt="Skeletonization + quality plotting across three smoothing states">
  <br>
  <sub><em><strong>Figure 3.</strong> Flagship plots <code>plot_layer()</code> and <code>plot_quality()</code> tracking how skeletonization and quality respond to smoothing (raw → truss → FEM). Reproduce: <code>python scripts/generate_3row_admesh.py</code>.</em></sub>
</p>

---

## Features

- **Fast** — full init + quality analysis on a 98,365-element mesh in **~3.3 s** (4.3× faster than v0.2.0).
  - Hash-mapped O(1) edge lookups, vectorised numpy core ops, kd-tree spatial queries at O(log n)
- **Mixed-element** — triangles, quads, and mixed meshes share one API
- **Smoothing** — three algorithms: Balendran direct FEM (one-shot solve), Zhou-Shimada angle-based (iterative), and the ADMESH Spring-Based Truss Smoother (force relaxation)
- **Mesh alterations** — `insert_vertex`, coord-only vertex moves, advancing-front element addition; topology-update primitives via the `MutableMesh` API (full mutation suite tracked in [#94](https://github.com/domattioli/CHILmesh/issues/94))
- **Analysis** — element quality, interior angles, layer-based skeletonization
- **I/O** — [ADCIRC](https://adcirc.org/) `.fort.14` and [SMS](https://www.aquaveo.com/sms) `.2dm` read/write. ([gmsh](https://gmsh.info/) Coming Soon)
- **Spatial queries** — point-in-element, k-nearest vertices, radius search.
- **Mesh traversal algorithms** (in-development)
- **ADMESH-Domains integration** — `from_admesh_domain()` adapter for catalog meshes

---

## Installation

From PyPI (pip):
```bash
pip install chilmesh
```

With [uv](https://docs.astral.sh/uv/) (faster, pip-compatible):
```bash
uv pip install chilmesh        # or:  uv add chilmesh
```

From conda-forge (once published):
```bash
conda install -c conda-forge chilmesh
# or: mamba install -c conda-forge chilmesh
```

From source:
```bash
git clone https://github.com/domattioli/CHILmesh && cd CHILmesh
pip install -e .
```

---

## Performance

CHILmesh is engineered for fast initialisation, query, and analysis on large unstructured 2D meshes. Hash-mapped edge adjacencies reduce topology build from `O(n²)` to amortised `O(n)`; core operations (`signed_area`, `interior_angles`, `elem_quality`) are fully vectorised over numpy arrays; a centroid kd-tree backs spatial queries (`find_element`, `nearest_vertices`) at `O(log n)` per call.

Reference workload: WNAT_Hagen (52,774 vertices · 98,365 elements).

| Stage | v0.2.0 | v0.4.0 |
|---|---:|---:|
| Fast init (no layers) | 3.9 s | **0.44 s** |
| Full init (with layers) | 7.7 s | **3.26 s** |
| Quality analysis | 6.6 s | **0.07 s** |
| **Total workflow** | **14.3 s** | **3.33 s** |
| `find_element` (per call) | n/a | **< 50 μs** |
| `Vert2Edge` lookup (per call) | 0.7 μs | **0.17 μs** |

Per-stage breakdown, methodology, and historical baselines in [`docs/BENCHMARK.md`](docs/BENCHMARK.md). Reproduce locally: `python scripts/benchmark_wnat_hagen.py --json results.json`.

---

## API Overview

```python
import chilmesh

# Load
mesh = chilmesh.examples.annulus()
mesh = chilmesh.CHILmesh.read_from_fort14('mesh.14')
mesh = chilmesh.CHILmesh.read_from_2dm('mesh.2dm')

# Smooth, analyse, visualise
mesh.smooth_mesh(method='fem', acknowledge_change=True)
quality, angles, stats = mesh.elem_quality()
mesh.plot()             # wireframe
mesh.plot_quality()     # per-element quality
mesh.plot_layer()       # skeletonization layers

# Skeletonization output
layers = mesh.layers    # {'OE', 'IE', 'OV', 'IV'} per layer

# Spatial queries (v0.3.0)
elem_id = mesh.find_element([0.5, 0.0])
neighbors = mesh.nearest_vertices([0.5, 0.0], k=5)
in_radius = mesh.find_elements_in_radius([0.5, 0.0], radius=0.2)
```

Full reference in [`docs/API.md`](docs/API.md). Optional ADMESH truss warm-start via `chilmesh.optimize_with_admesh_truss`.

---

## Mesh Smoothing

Three smoothing algorithms — pick by use case. Each preserves boundary nodes, leaves topology unchanged, and accepts mixed-element meshes.

| Algorithm | API | Style | When |
|---|---|---|---|
| **Balendran direct FEM** | `smooth_mesh(method='fem', ...)` → `direct_smoother(kinf=1e12)` | One-shot sparse solve | Best general-purpose default. Stable on tri / quad / mixed. |
| **Zhou-Shimada angle-based** | `smooth_mesh(method='angle-based', ...)` → `angle_based_smoother(n_iter, omega, tol)` | Iterative, angle-maximising | DOMsmooth hybrid fallback for difficult mixed meshes where FEM stalls. |
| **ADMESH Spring-Based Truss Smoother** | `chilmesh.optimize_with_admesh_truss(mesh, sdf, niter, Fscale)` | distmesh2d-style spring/force relaxation against a signed-distance field | When you want quality gains plus boundary nodes that respect a domain SDF (e.g., coastline). |

```python
mesh.smooth_mesh(method='fem', acknowledge_change=True)         # default
mesh.smooth_mesh(method='angle-based', acknowledge_change=True) # fallback
mesh = chilmesh.optimize_with_admesh_truss(mesh, sdf, niter=500, Fscale=0.5)
```

Stiffness assembly, convergence parameters, and algorithm details: [`docs/API.md`](docs/API.md).

**References.**
- FEM smoother: Balendran, B. (1999). *"A direct smoothing method for surface meshes."* Proc. 8th International Meshing Roundtable, pp. 189–193.
- Angle-based smoother: Zhou, M. & Shimada, K. (2000). *"An angle-based approach to two-dimensional mesh smoothing."* Proc. 9th IMR, pp. 373–384.
- ADMESH Spring-Based Truss Smoother: Conroy et al. (2012) *"ADMESH: An advanced, automatic unstructured mesh generator for shallow water models."* [doi:10.1007/s10236-012-0574-0](https://doi.org/10.1007/s10236-012-0574-0).

---

## Examples

Runnable scripts in [`examples/`](examples/) demonstrate common tasks against bundled fixtures — no external mesh files required:

- [`01_quickstart.py`](examples/01_quickstart.py) — load a mesh, print stats, save a plot
- [`02_fort14_roundtrip.py`](examples/02_fort14_roundtrip.py) — fort.14 read / write
- [`03_smoothing.py`](examples/03_smoothing.py) — angle-based smoother on perturbed interior
- [`04_spatial_queries.py`](examples/04_spatial_queries.py) — `find_element`, radius search, k-nearest vertices

```bash
python examples/01_quickstart.py
```

---

## CLI

`chilmesh` ships with a small shell entry point for inspection, conversion, smoothing, and plotting. No new dependencies — pure stdlib `argparse` over the existing public API.

```bash
# Mesh statistics (verts, elems, edges, layers, quality)
chilmesh info path/to/mesh.fort.14

# Format conversion (output format inferred from suffix)
chilmesh convert mesh.2dm mesh.fort.14

# In-place smoothing
chilmesh smooth mesh.fort.14 -o smoothed.fort.14 --method angle-based --iter 50

# Static figure (PNG / PDF / SVG by suffix; --layers or --quality for overlays)
chilmesh plot mesh.fort.14 -o mesh.png --quality
```

Each subcommand has its own `--help` with an example. Also available as `python -m chilmesh ...` when the script isn't on PATH.

---

## Downstream Projects

[**ADMESH**](https://github.com/domattioli/ADMESH) — Optimized 2D triangular mesh generation for hydrodynamic domains
[**MADMESHR**](https://github.com/domattioli/MADMESHR) — AI based quad- and mixed element generation for hydrodynamic domains.
[**ADMESH-Domains**](https://github.com/domattioli/ADMESH-Domains) — Mesh catalog for hydrodynamic domains.

---

## Contributing

Issues and pull requests welcome at [github.com/domattioli/CHILmesh](https://github.com/domattioli/CHILmesh). Run `pytest -v` before opening a PR — see [`TESTING.md`](TESTING.md) for the test-suite guide.

---

## Citation

CHILmesh originated in MATLAB as the mixed-element data structure backing a skeletonization-driven heuristic for indirect triangle-to-quad conversion that preserves the underlying size function (Mattioli, OSU MSc [thesis](https://github.com/user-attachments/files/19727573/QuADMESH__Thesis_Doc.pdf), 2017). This Python implementation is the actively-developed successor, with `.fort.14` I/O and a shared API for downstream projects (MADMESHR, ADMESH, ADMESH-Domains).

```bibtex
@mastersthesis{mattioli2017quadmesh,
  author       = {Mattioli, Dominik O.},
  title        = {{QuADMESH+}: A Quadrangular ADvanced Mesh Generator for Hydrodynamic Models},
  school       = {The Ohio State University},
  year         = {2017},
  url          = {http://rave.ohiolink.edu/etdc/view?acc_num=osu1500627779532088}
}
```

---

## References

- [Source of MATLAB implementation (Mattioli, 2017)](https://github.com/user-attachments/files/19727573/QuADMESH__Thesis_Doc.pdf)

---

## License

MIT License — See [LICENSE](LICENSE) for details.
