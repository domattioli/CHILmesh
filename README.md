<h1 align="center">
  CHILmesh
</h1>

<p align="center">
  <strong>Fast 2D mesh representation and analysis for triangular, quadrilateral, and mixed-element meshes of hydrodynamic domains.</strong>
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

> **Note for MATLAB users**: Python successor to the original [MATLAB QuADMesh+](https://github.com/domattioli/QuADMesh-MATLAB). API may evolve until v1.0.

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

Median quality 0.797. Full init + analysis: **~3.3 seconds** end-to-end. Reproduce: `python scripts/benchmark_wnat_hagen.py`.

### Showcase: Mixed-Element Mesh

466 corner-graded triangles + 60 quad core, smoothed with Zhou & Shimada angle-based method (boundary pinned, median quality 0.733):

![Mixed-element mesh: wireframe, layers, quality](output/mixed_mesh_showcase.png?v=1)

Structured quad core (192 quads) → skeletonized outer ring (ADMESH 394 tris) + gap band Delaunay fill (72 tris) → angle-based smooth. Reproduce with `python scripts/generate_mixed_truss_demo.py`.

### Showcase: Skeletonization & Mesh Plotting

CHILmesh's two flagship visualizations — **layer-based skeletonization** (center, viridis) and **per-element quality plotting** (right, cool, `4√3·area / Σedge²`) — rendered on three states of the same triangular annulus:

![CHILmesh skeletonization layers and quality plot across three smoothing states](output/annulus_quickstart.png?v=6)

**Rows**: (1) raw annulus (q ≈ 0.71), (2) ADMESH warm-start (q ≈ 0.92), (3) FEM smooth (q ≈ 0.93). Same boundary + 580 triangles.

**Regenerate**: `python scripts/generate_3row_admesh.py` → `output/annulus_quickstart.png`. See [`src/chilmesh/admesh_warmstart.py`](src/chilmesh/admesh_warmstart.py) for warm-start details.

---

## Features

- **Fast**: 937× speedup vs v0.1.1 through optimized data structures
- **Mixed-Element**: Triangles, quads, and mixed meshes with unified API
- **Smoothing**: FEM and geometric smoothing for quality improvement
- **Analysis**: Element quality metrics, interior angles, layer-based skeletonization
- **I/O**: Read/write ADCIRC `.fort.14` and SMS `.2dm` formats
- **Catalog**: ADMESH-Domains integration via `from_admesh_domain()` adapter

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

**4,000×+ faster** than v0.1.1. Benchmarks on WNAT_Hagen (52.7k verts · 98.4k elems):

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

# Load
mesh = chilmesh.examples.annulus()
mesh = chilmesh.CHILmesh.read_from_fort14('mesh.14')

# Smooth & Analyze
mesh.smooth_mesh(method='fem', acknowledge_change=True)
quality, angles, stats = mesh.elem_quality()

# Visualize
mesh.plot()           # wireframe
mesh.plot_quality()   # quality colormap
mesh.plot_layer()     # skeletonization layers

# Topology
boundary_nodes = mesh.boundary_node_indices()
layers = mesh.layers  # {'OE', 'IE', 'OV', 'IV'} per layer

# ADMESH warm-start (optional)
import numpy as np
sdf = lambda p: np.linalg.norm(p, axis=1) - 1.0
mesh = chilmesh.optimize_with_admesh_truss(mesh, sdf, niter=500)
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
