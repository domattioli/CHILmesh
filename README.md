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

> **Note for MATLAB users**: This Python implementation is the actively-developed successor to the original MATLAB QuADMesh+ codebase. It is **still in development** and the API may evolve. The original MATLAB code (no longer maintained) remains available for reference at [domattioli/QuADMesh-MATLAB](https://github.com/domattioli/QuADMesh-MATLAB) — see `00_CHILMesh_Class/@CHILmesh/CHILmesh.m` for the canonical algorithms (e.g., `meshLayers` skeletonization).

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

### Visual Comparison: Warm-Start Truss Optimization Pipeline

See how a raw Delaunay triangulation is progressively refined through warm-start truss equilibration and FEM smoothing — all three rows share the same boundary and element count (580 triangles):

![CHILmesh quickstart: raw Delaunay → ADMESH warm-start truss → FEM smoother](tests/output/annulus_quickstart.png?v=3)

**Row 1 — Raw Delaunay:** Unsmoothed input mesh from `chilmesh.examples.annulus()` (median quality ≈ 0.71).
**Row 2 — Row 1 + ADMESH Truss (warm-start):** Vendored `distmesh2d` truss loop, started from Row 1's points (boundary pinned bit-exactly via `pfix`). The vendored loop tracks the best-quality state across iterations and early-stops if median quality drops more than 10% from peak. Median quality jumps to ≈ 0.94.
**Row 3 — Row 2 + FEM Smoother:** CHILmesh's FEM relaxation applied to Row 2.

Columns: **left** = mesh wireframe · **center** = skeletonization layers (viridis) · **right** = per-element quality (cool, 4√3·area / Σedge²).

#### How to regenerate this figure

```bash
python generate_3row_admesh.py
```

The script lives at the repo root and writes `tests/output/annulus_quickstart.png`. It enforces fail-loud assertions for boundary preservation (V_BND), degeneracy (V_DEGENERACY), sibling chain (V_CHAIN), and positive-area connectivity (V_CONN) before saving. See [`src/chilmesh/admesh_warmstart.py`](src/chilmesh/admesh_warmstart.py) for the warm-start adapter API and [`specs/005-admesh-warm-start-truss/`](specs/005-admesh-warm-start-truss/) for the full design contract.

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

## Performance (v0.2.0)

**937× faster** than v0.1.1 through systematic optimization.

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

- **Triangles**: 3-vertex elements
- **Quads**: 4-vertex elements  
- **Mixed**: Both types in same mesh (triangles padded to 4 columns)

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
