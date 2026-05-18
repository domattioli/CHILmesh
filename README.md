<h1 align="center">
  CHILmesh
</h1>

<p align="center">
  <strong>Fast 2D mesh generation and analysis for triangular, quadrilateral, and mixed-element meshes.</strong>
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
  <a href="https://github.com/domattioli/CHILmesh/blob/main/LICENSE"><img src="https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square" alt="License"></a>
</p>

> **Note for MATLAB users**: This Python implementation is the actively-developed successor to the original MATLAB QuADMesh+ codebase. It is **still in development** and the API may evolve. The original MATLAB code (no longer maintained) remains available for reference at [domattioli/QuADMesh-MATLAB](https://github.com/domattioli/QuADMesh-MATLAB).

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

### WNAT_Hagen — 52,774 vertices · 98,365 elements

![WNAT_Hagen quality plot and distribution](output/wnat_hagen_showcase.png?v=2)

Per-element quality (skew, `4√3·area / Σedge²`) and a 100-bin distribution. Median quality 0.797, full init + quality analysis in **~3.3 s** end-to-end. Reproduce: `python scripts/benchmark_wnat_hagen.py`.

### Mixed-element mesh — quads + ADMESH tri ring

![Mixed-element mesh: wireframe, layers, quality](output/mixed_mesh_showcase.png?v=2)

466 corner-graded triangles surrounding 60 quads after FEM smoothing (symmetric quad stiffness, boundary pinned, median quality 0.760). 16×12 quad core, distmesh1d-graded perimeter (`h(p) = 0.05 + 0.45·(1 − exp(−(d/0.5)²))`), Delaunay-stitched gap band, then FEM smoother. Reproduce: `python scripts/generate_mixed_truss_demo.py`.

### Skeletonization + quality plotting (3 smoothing states)

![CHILmesh skeletonization layers and quality plot across three smoothing states](output/annulus_quickstart.png?v=6)

The two flagship visualisations — layer-based skeletonization (centre) and per-element quality (right) — on three states of the same 580-triangle annulus: raw (median q ≈ 0.71), ADMESH warm-start truss (≈ 0.92), and FEM smoother (≈ 0.93). Reproduce: `python scripts/generate_3row_admesh.py`.

---

## Features

- **Fast** — 4,000×+ workflow speedup vs v0.1.1 via hash-mapped adjacencies and vectorised core ops
- **Mixed-element** — triangles, quads, and mixed meshes share one API
- **Smoothing** — angle-based FEM smoother for quality improvement (Zhou & Shimada 2000)
- **Analysis** — element quality, interior angles, layer-based skeletonization
- **I/O** — ADCIRC `.fort.14` and SMS `.2dm` read/write
- **Spatial queries** — point-in-element, k-nearest vertices, radius search (v0.3.0)
- **ADMESH-Domains integration** — `from_admesh_domain()` adapter for catalog meshes

---

## Installation

From PyPI:
```bash
pip install chilmesh
```

From source:
```bash
git clone https://github.com/domattioli/CHILmesh && cd CHILmesh
pip install -e .
```

---

## Performance (v0.3.0)

WNAT_Hagen workflow (52,774 vertices · 98,365 elements):

| Stage | v0.1.1 | v0.3.0 | Speedup |
|---|---:|---:|---:|
| **Total workflow** | 13,400 s | **3.33 s** | **4,027×** |

Full breakdown (init, quality, per-query latency) and methodology in [`docs/BENCHMARK.md`](docs/BENCHMARK.md). Reproduce: `python scripts/benchmark_wnat_hagen.py --json results.json`.

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

Angle-based FEM smoother (Zhou & Shimada 2000) for triangles, quads, and mixed meshes. One API across element types; boundary nodes are pinned, topology preserved, aspect ratio favoured.

```python
mesh.smooth_mesh(method='fem', acknowledge_change=True)        # any element type
new_points = mesh.direct_smoother(kinf=1e12)                   # boundary stiffness
```

Parameters, stiffness assembly details, and the angle-based fallback (`mesh.smooth_mesh(method='angle-based', ...)`) for mixed meshes are documented in [`docs/API.md`](docs/API.md).

**Reference:** Zhou, M., & Shimada, K. (2000). "An angle-based approach to two-dimensional mesh smoothing." *Proc. 9th International Meshing Roundtable*, 373–384.

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

## Downstream Projects

[**MADMESHR**](https://github.com/domattioli/MADMESHR) — Advancing-front mesh adaptation built on CHILmesh
[**ADMESH**](https://github.com/domattioli/ADMESH) — Optimized mesh generation and smoothing
[**ADMESH-Domains**](https://github.com/domattioli/ADMESH-Domains) — Mesh catalog for hydrodynamic applications

---

## Contributing

Issues and pull requests welcome at [github.com/domattioli/CHILmesh](https://github.com/domattioli/CHILmesh). Run `pytest -v` before opening a PR — see [`TESTING.md`](TESTING.md) for the test-suite guide.

---

## Citation

```bibtex
@mastersthesis{mattioli2017quadmesh,
  author       = {Mattioli, Dominik O.},
  title        = {{QuADMESH+}: A Quadrangular ADvanced Mesh Generator for Hydrodynamic Models},
  school       = {The Ohio State University},
  year         = {2017},
  url          = {http://rave.ohiolink.edu/etdc/view?acc_num=osu1500627779532088}
}
```

[Read thesis PDF](https://github.com/user-attachments/files/19727573/QuADMESH__Thesis_Doc.pdf)

---

## References

- [FEM Smoother (Zhou & Shimada, 2000)](https://api.semanticscholar.org/CorpusID:34335417)
- [Angle-Based Smoothing](https://www.andrew.cmu.edu/user/shimada/papers/00-imr-zhou.pdf)
- [ADMESH Paper](https://doi.org/10.1007/s10236-012-0574-0)
- Original MATLAB implementation funded by [Aquaveo](https://aquaveo.com/)

---

## License

MIT License — See [LICENSE](LICENSE) for details.
