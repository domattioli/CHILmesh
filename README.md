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

### WNAT_Hagen — 52,774 vertices · 98,365 elements

![WNAT_Hagen quality plot and distribution](output/wnat_hagen_showcase.png?v=2)

Per-element quality (skew, `4√3·area / Σedge²`) and a 100-bin distribution. Median quality 0.797, full init + quality analysis in **~3.3 s** end-to-end. Reproduce: `python scripts/benchmark_wnat_hagen.py`.

### Mixed-element mesh — quads + ADMESH tri ring

![Mixed-element mesh: wireframe, layers, quality](output/mixed_mesh_showcase.png?v=2)

Demonstrates **mixed-element support** end-to-end: a quad core stitched to an ADMESH triangle ring runs through wireframe rendering, layer-based skeletonization, and per-element quality analysis on a single mesh. Reproduce: `python scripts/generate_mixed_truss_demo.py`.

### Skeletonization + quality plotting (3 smoothing states)

![CHILmesh skeletonization layers and quality plot across three smoothing states](output/annulus_quickstart.png?v=6)

Shows the **two flagship visualisations** — `plot_layer()` (centre) and `plot_quality()` (right) — applied to the same annulus at three smoothing states (raw, ADMESH warm-start truss, FEM smoother). Same API, same fixture, different inputs: how skeletonization and quality reads track mesh evolution. Reproduce: `python scripts/generate_3row_admesh.py`.

---

## Features

- **Fast** — hash-mapped adjacencies and vectorised core ops; large meshes (~100k elements) initialise in seconds, not hours
- **Mixed-element** — triangles, quads, and mixed meshes share one API
- **Smoothing** — angle-based FEM smoother for quality improvement (Zhou & Shimada 2000)
- **Analysis** — element quality, interior angles, layer-based skeletonization
- **I/O** — ADCIRC `.fort.14` and SMS `.2dm` read/write
- **Spatial queries** — point-in-element, k-nearest vertices, radius search (v0.3.0)
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

Reference workload: WNAT_Hagen (52,774 vertices · 98,365 elements) initialises end-to-end in a few seconds on commodity hardware. Full numbers, per-stage breakdown, and reproducibility scripts in [`docs/BENCHMARK.md`](docs/BENCHMARK.md). Reproduce locally: `python scripts/benchmark_wnat_hagen.py --json results.json`.

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
