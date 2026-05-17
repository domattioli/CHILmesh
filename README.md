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
- [Performance](#performance-v020)
- [API Overview](#api-overview)
- [Mesh Element Types](#mesh-element-types)
- [Downstream Projects](#downstream-projects)
- [Citation](#citation)
- [References](#references)
- [License](#license)

> **MATLAB users**: Python successor to the original QuADMesh+ codebase. Still in development; API may evolve. Original MATLAB code (no longer maintained) at [domattioli/QuADMesh-MATLAB](https://github.com/domattioli/QuADMesh-MATLAB) — `00_CHILMesh_Class/@CHILmesh/CHILmesh.m` has canonical algorithms (e.g., `meshLayers` skeletonization).

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

![CHILmesh quickstart: raw Delaunay → ADMESH warm-start truss → FEM smoother](output/annulus_quickstart.png?v=5)

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

## Performance (v0.2.0)

| Operation | v0.1.1 | v0.2.0 | Speedup |
|-----------|--------|--------|---------|
| Fast init (52.7k verts) | 3,200s | 3.9s | **822×** |
| Full init (with layers) | 5,400s | 7.7s | **701×** |
| Quality analysis | 4,800s | 6.6s | **727×** |
| **Total workflow** | 13,400s | 14.3s | **937×** |

See [BENCHMARK.md](BENCHMARK.md) for detailed methodology.

---

## API Overview

```python
import chilmesh

# Load examples or from file
mesh = chilmesh.examples.annulus()
mesh = chilmesh.CHILmesh.read_from_fort14('mesh.14')
mesh = chilmesh.CHILmesh.read_from_2dm('mesh.2dm')

# Smooth mesh
mesh.smooth_mesh(method='fem', acknowledge_change=True)

# Warm-start an existing triangulation through ADMESH's distmesh truss loop.
# Boundary points are pinned bit-exactly; interior is relaxed to equilibrium.
import numpy as np
sdf = lambda p: np.maximum(np.linalg.norm(p, axis=1) - 1.0,
                            0.3 - np.linalg.norm(p, axis=1))
mesh = chilmesh.optimize_with_admesh_truss(
    mesh, sdf, niter=500, deltat=0.02, Fscale=0.5
)

# Analyze
quality, angles, stats = mesh.elem_quality()
interior_angles = mesh.interior_angles()

# Visualize
mesh.plot()
mesh.plot_quality()
mesh.plot_layer()

# Layer structure (skeletonization)
layers = mesh.layers  # {'OE': [...], 'IE': [...], 'OV': [...], 'IV': [...]}

# Access topology
edges = mesh.boundary_edges()
boundary_nodes = mesh.boundary_node_indices()
```

---

## Mesh Element Types

- **Triangles**: 3-vertex
- **Quads**: 4-vertex
- **Mixed**: both types (triangles padded to 4 columns)

---

## Downstream Projects

[**MADMESHR**](https://github.com/domattioli/MADMESHR) — Advancing-front mesh adaptation
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
