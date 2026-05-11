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
  <a href="https://github.com/domattioli/CHILmesh/actions/workflows/test.yml"><img src="https://img.shields.io/github/actions/workflow/status/domattioli/CHILmesh/test.yml?label=Tests&logo=github" alt="Tests"></a>
  <a href="https://github.com/domattioli/CHILmesh/blob/main/LICENSE"><img src="https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square" alt="License"></a>
</p>

---

## Table of Contents

- [Quick Start](#quick-start)
- [Features](#features)
- [Installation](#installation)
- [Performance](#performance-v030)
- [API Overview](#api-overview)
- [Downstream Projects](#downstream-projects)
- [Citation](#citation)
- [References](#references)
- [License](#license)

> **MATLAB users**: Python successor to original QuADMesh+ codebase. Still in development; API may evolve. Original MATLAB (unmaintained) at [domattioli/QuADMesh-MATLAB](https://github.com/domattioli/QuADMesh-MATLAB) — canonical algorithms (e.g., `meshLayers`) in `00_CHILMesh_Class/@CHILmesh/CHILmesh.m`.

---

## Quick Start

Generate, smooth, and analyze 2D meshes in seconds:

```python
import chilmesh
import matplotlib.pyplot as plt

# Load example mesh
mesh = chilmesh.examples.annulus()

# Smooth with FEM formulation
mesh.smooth_mesh(method='fem', acknowledge_change=True)

# Analyze quality
quality, angles, stats = mesh.elem_quality()
mesh.plot_quality()
plt.show()
```

### Showcase: WNAT_Hagen (52,774 verts · 98,365 elems)

![WNAT_Hagen quality plot and distribution](output/wnat_hagen_showcase.png?v=2)

Median q=0.797, mean 0.786. Full init + analysis: **~3.3s** end-to-end. Reproduce: `python scripts/benchmark_wnat_hagen.py`.

### Showcase: Mixed-Element Mesh

466 tris + 60 quads after FEM smoothing (symmetric quad stiffness, boundary pinned, q=0.760):

![Mixed-element mesh: wireframe, layers, quality](output/mixed_mesh_showcase.png?v=2)

Structured quad core (192 quads) → skeletonized outer ring (ADMESH 394 tris) + gap band Delaunay fill (72 tris) → FEM smooth. Reproduce: `python scripts/generate_mixed_truss_demo.py`.

### Showcase: Skeletonization & Mesh Plotting

Layer-based skeletonization (center, viridis) and quality plotting (right, `4√3·area / Σedge²`) on three smoothing states (same boundary, 580 triangles):

![CHILmesh skeletonization layers and quality plot across three smoothing states](output/annulus_quickstart.png?v=6)

Rows: (1) raw Delaunay (q≈0.71), (2) ADMESH warm-start (q≈0.92), (3) FEM smooth (q≈0.93).

**Regenerate**: `python scripts/generate_3row_admesh.py` → `output/annulus_quickstart.png`. See [`src/chilmesh/admesh_warmstart.py`](src/chilmesh/admesh_warmstart.py) and [`specs/005-admesh-warm-start-truss/`](specs/005-admesh-warm-start-truss/) for details.

---

## Features

- **Fast**: 937× speedup vs v0.1.1
- **Mixed-Element**: Triangles, quads, mixed meshes
- **Smoothing**: FEM and geometric smoothing
- **Analysis**: Quality, angles, skeletonization layers
- **I/O**: ADCIRC `.fort.14` and SMS `.2dm`
- **Catalog**: ADMESH-Domains integration

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

**4,000×+ faster** than v0.1.1 through systematic optimization. Reference mesh: **WNAT_Hagen — 52,774 vertices · 98,365 elements · 151,248 edges · 30 layers** (see [showcase image](#showcase-wnat_hagen-52774-vertices--98365-elements) above).

### Initialization

| Operation | v0.1.1 (est.) | v0.3.0 (measured) | Speedup |
|-----------|--------------:|------------------:|--------:|
| Fast init (no layers) | 3,200s | **0.44s** | **7,307×** |
| Full init (with layers) | 5,400s | **3.26s** | **1,658×** |
| Quality analysis | 4,800s | **0.07s** | **68,175×** |
| **Total workflow** | **13,400s** | **3.33s** | **4,027×** |

### Query latency (per call)

| Operation | v0.1.1 (est.) | v0.3.0 (measured) | Speedup |
|-----------|--------------:|------------------:|--------:|
| `elem2edge` (5k samples) | 2,000μs | **2.08μs** | **963×** |
| `Vert2Edge` lookup (5k samples) | 3,500μs | **0.17μs** | **21,092×** |
| `Elem2Edge` bulk (1k samples) | 4,500μs | **0.14μs** | **32,766×** |

Reproduce: `python scripts/benchmark_wnat_hagen.py --json results.json`. Historical April 2026 release numbers and methodology in [docs/BENCHMARK.md](docs/BENCHMARK.md).

---

## API Overview

```python
import chilmesh

# Load: built-in examples or fort.14 / 2dm
mesh = chilmesh.examples.annulus()
mesh = chilmesh.CHILmesh.read_from_fort14('mesh.14')
mesh = chilmesh.CHILmesh.read_from_2dm('mesh.2dm')

# Smooth (FEM or geometric)
mesh.smooth_mesh(method='fem', acknowledge_change=True)

# Analyze
quality, angles, stats = mesh.elem_quality()
interior_angles = mesh.interior_angles()

# Visualize
mesh.plot()                    # wireframe
mesh.plot_quality()            # per-element quality colormap
mesh.plot_layer()              # skeletonization layers
mesh.plot_boundary()           # boundary edges highlighted
mesh.plot_interior_edges()     # interior edges only

# Skeletonization output
layers = mesh.layers     # {'OE', 'IE', 'OV', 'IV'} per layer

# Topology
edges = mesh.boundary_edges()
boundary_nodes = mesh.boundary_node_indices()

# Optional: warm-start through ADMESH truss (boundary pinned bit-exact)
import numpy as np
sdf = lambda p: np.maximum(np.linalg.norm(p, axis=1) - 1.0,
                            0.3 - np.linalg.norm(p, axis=1))
mesh = chilmesh.optimize_with_admesh_truss(mesh, sdf, niter=500, Fscale=0.5)
```

---

## Downstream Projects

[**MADMESHR**](https://github.com/domattioli/MADMESHR) — Advancing-front mesh adaptation built on CHILmesh  
[**ADMESH**](https://github.com/domattioli/ADMESH) — Optimized mesh generation and smoothing  
[**ADMESH-Domains**](https://github.com/domattioli/ADMESH-Domains) — Mesh catalog for hydrodynamic applications  

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
