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

## Table of Contents

- [Quick Start](#quick-start)
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

### Showcase: WNAT_Hagen (52,774 vertices · 98,365 elements)

Reference benchmark mesh — per-element quality (skew, `4√3·area / Σedge²`) and granular distribution histogram (100 bins):

![WNAT_Hagen quality plot and distribution](output/wnat_hagen_showcase.png?v=2)

Median quality 0.797, mean 0.786 across all 98k elements. Full init + quality analysis: **~3.3 seconds** end-to-end (see [Performance](#performance-v020) below). Reproduce with `python scripts/benchmark_wnat_hagen.py`.

### Showcase: Mixed-Element Mesh from Quad Core + ADMESH Tri Ring

466 corner-graded triangles surrounding 60 quads after FEM smoothing (symmetric quad stiffness, boundary pinned, median quality 0.760) — wireframe, skeletonization layers, and per-element quality:

![Mixed-element mesh: wireframe, layers, quality](output/mixed_mesh_showcase.png?v=2)

**Pipeline.** Start with a 16×12 structured quad rectangle (192 quads, 6 skeleton layers). Strip layers 0–1 as the ADMESH domain, drop layer 2, retain layers 3+ as the quad core. Apply `distmesh1d` to the outer rectangle perimeter with a corner-dense edge-length field `h(p) = 0.05 + 0.45·(1 − exp(−(d/0.5)²))` (d = distance to nearest corner) — perimeter edge lengths now span 0.035–0.499 vs. uniform 0.250. Grid-sample interior points within the ring SDF and run ADMESH full distmesh with both the outer rectangle perimeter and the layer-1/2 seam pinned — produces 394 quality-graded triangles densest at the four 90° corners. Delaunay-triangulate the layer-2 gap band from boundary nodes only (72 tris), stitch with the 60 quads into a combined mesh, then run the FEM smoother (symmetric quad decomposition, boundary pinned). Reproduce: `python scripts/generate_mixed_truss_demo.py` (also writes a 4-panel pipeline diagram to `output/mixed_truss_fem_demo.png`).

### Showcase: Skeletonization & Mesh Plotting

CHILmesh's two flagship visualizations — **layer-based skeletonization** (center, viridis) and **per-element quality plotting** (right, cool, `4√3·area / Σedge²`) — rendered on three states of the same triangular annulus:

![CHILmesh skeletonization layers and quality plot across three smoothing states](output/annulus_quickstart.png?v=6)

**Rows.** Row 1: raw `chilmesh.examples.annulus()` (median quality ≈ 0.71). Row 2: ADMESH warm-start truss applied to Row 1 (≈ 0.92). Row 3: CHILmesh FEM smoother applied to Row 2 (≈ 0.93). All three rows share the same boundary and 580 triangles; the smoothing passes are shown only to give the skeletonization & quality plots distinct inputs.

**Mixed-Element Showcase:** CHILmesh now supports mixed-element meshes (triangles + quads) throughout the full pipeline: load → adjacencies → skeletonization → FEM smoothing → analysis. The FEM smoother uses element-type-specific stiffness matrices (Zhou & Shimada for triangles, quad extension via analogy) with angle-based RHS forces that preserve ideal aspect ratios. Reproduce: `python scripts/generate_mixed_truss_demo.py`.

#### Regenerate

```bash
python scripts/generate_3row_admesh.py
```

Writes `output/annulus_quickstart.png`. Fail-loud assertions: boundary preservation (V_BND, V_BND_PROP), positive-area connectivity (V_CONN), sibling chain (V_CHAIN). See [`src/chilmesh/admesh_warmstart.py`](src/chilmesh/admesh_warmstart.py) for the warm-start adapter and [`specs/005-admesh-warm-start-truss/`](specs/005-admesh-warm-start-truss/) for the full contract.

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

## Mesh Smoothing

CHILmesh provides **angle-based FEM smoothing** for quality improvement across all element types: triangles, quads, and mixed-element meshes.

### Algorithm: Zhou & Shimada Angle-Based Smoothing

The smoother minimizes angular distortion via the Zhou & Shimada formulation:

- **Energy:** Element stiffness matrices computed from local angle and edge-length geometry
- **Objective:** Interior nodes driven toward configurations that maximize angles (60° ideal for triangles, 90° for quads)
- **Boundary:** Fixed pinned constraints; boundary nodes remain unchanged
- **Solver:** Direct (non-iterative) sparse linear solve via SciPy

### Element-Type Support

| Element Type | Stiffness Matrix | RHS Forces | Status |
|--------------|-----------------|-----------|--------|
| Triangles | Zhou & Shimada (D, T matrices) | Angle-based cotangent weighting | ✓ Production |
| Quads | Extended via analogy (2×2 integration) | Angle-based cotangent weighting | ✓ Production |
| Mixed | Per-element-type assembly | Per-element-type angle forces | ✓ Production |

**Reference:** Zhou, M., & Shimada, K. (2000). "An angle-based approach to two-dimensional mesh smoothing". *Proceedings of the 9th International Meshing Roundtable*, 373–384.

### Usage

```python
# Triangular mesh
mesh.smooth_mesh(method='fem', acknowledge_change=True)

# Quad or mixed-element mesh (same API)
mesh.smooth_mesh(method='fem', acknowledge_change=True)

# Direct access to smoother (returns new point array)
new_points = mesh.direct_smoother(kinf=1e12)  # kinf = boundary stiffness
```

### Properties

- **Preserves topology:** Element-node connectivity unchanged
- **Preserves boundary:** Boundary nodes remain fixed at input positions
- **Preserves aspect ratio:** Perfect-grid meshes (square quads, equilateral triangles) are minimally perturbed
- **Improves quality:** Interior nodes repositioned to maximize angles

### Limitations

- One-shot smoothing (not iterative); convergence is direct solve
- Assumes 2D planar meshes; z-coordinate preserved but not optimized
- No mesh mutation (add/remove elements); topology fixed

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
