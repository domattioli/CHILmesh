<h1 align="center">
  CHILmesh: representing triangular, quadrangular and mixed-element (2D) meshes for advanced and automatic mesh generation for hydrodynamic domains.
</h1>

<p align="center">
  <strong><a href="https://scholar.google.com/citations?user=IBFSkOcAAAAJ&hl=en">Dominik Mattioli</a><sup>1†</sup>, <a href="https://scholar.google.com/citations?user=mYPzjIwAAAAJ&hl=en">Ethan Kubatko</a><sup>2</sup></strong><br>
  <sup>†</sup>Corresponding author<br><br>
  <sup>1</sup>Penn State University<br>
  <sup>2</sup>Computational Hydrodynamics and Informatics Lab (CHIL), The Ohio State University
</p>

<p align="center">
  <a href="https://pypi.org/project/chilmesh/">
    <img src="https://img.shields.io/pypi/v/chilmesh?logo=python&logoColor=white&label=PyPI" alt="PyPI version">
  </a>
  <a href="https://pypi.org/project/chilmesh/">
    <img src="https://img.shields.io/pypi/pyversions/chilmesh?logo=python&logoColor=white" alt="Python 3.10+">
  </a>
  <a href="https://github.com/domattioli/CHILmesh/releases">
    <img src="https://img.shields.io/github/v/release/domattioli/CHILmesh?logo=github&logoColor=white" alt="Latest release">
  </a>
  <a href="https://github.com/domattioli/CHILmesh/actions/workflows/python-package.yml">
    <img src="https://img.shields.io/github/actions/workflow/status/domattioli/CHILmesh/python-package.yml?logo=github&logoColor=white&label=tests" alt="Tests">
  </a>
</p>

<p align="center">
  <a href="https://ceg.osu.edu/computational-hydrodynamics-and-informatics-laboratory">
    <img src="https://img.shields.io/badge/CHIL%20Lab%20@%20OSU-a7b1b7?logo=academia&logoColor=ba0c2f&labelColor=ba0c2f" alt="CHIL Lab @ OSU">
  </a>
  <a href="https://ceg.osu.edu/computational-hydrodynamics-and-informatics-laboratory">
    <img src="https://img.shields.io/badge/OSU_CHIL-ADMESH-66bb33?logo=github&logoColor=ba0c2f&labelColor=ffffff" alt="OSU CHIL ADMESH">
  </a>
  <a href="https://github.com/user-attachments/files/19724263/QuADMESH-Thesis.pdf">
    <img src="https://img.shields.io/badge/Thesis-QuADMESH-ba0c2f?style=flat-square&logo=book&logoColor=white&labelColor=cfd4d8" alt="QuADMESH Thesis">
  </a>
  <a href="https://scholar.google.com/citations?view_op=view_citation&hl=en&user=IBFSkOcAAAAJ&citation_for_view=IBFSkOcAAAAJ:u5HHmVD_uO8C">
    <img src="https://img.shields.io/badge/Scholar-Profile-4285F4?logo=google-scholar&logoColor=white" alt="Google Scholar">
  </a>
  <a href="https://www.mathworks.com/matlabcentral/fileexchange/135632-chilmesh">
    <img src="https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg" alt="MathWorks File Exchange">
  </a>
  <a href="https://github.com/domattioli/CHILmesh/blob/d63b7d221842cbb00bdb057b201519ac5e49febc/LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square" alt="License: MIT">
  </a>
</p>
<p align="center">
  <img src="https://github.com/user-attachments/assets/0c344383-cde0-454f-810f-5407092a7be2" alt="image">
</p>


## Releases
- **2026/04 0.2.0** - Performance modernization release
  - 937x speedup on large meshes through systematic optimization (Phases 1-4)
  - Bridge adapters for downstream integration (MADMESHR, ADMESH, ADMESH-Domains)
  - Advancing-front mesh generation API
  - See [CHANGELOG.md](CHANGELOG.md) for details
- 2026/04 0.1.1 - Honest hotfix: bug-fix release with regression tests, package data, and chilmesh.examples
- 2025/04/12 Python version of our code
- 2023/09/19 MATLAB code revisited; repo initiated
- 2017/08/01 Nascent MATLAB version of the code

## Table of Contents
- [Releases](#releases)
- [Installation](#installation)
- [Key Features](#key-features)
- [Performance](#performance-v020)
- [Example Usage](#example-usage)
- [See also](#see-also--downstream-projects)
- [BibTeX](#bibtex)
- [Acknowledgements](#acknowledgements)


## Installation
Install the latest release from PyPI:
```bash
pip install chilmesh
```
Or install from source:
```bash
git clone https://github.com/domattioli/CHILmesh && cd CHILmesh
python -m venv .myenv
source .myenv/bin/activate
pip install -e .
```


## Key Features
- Minimal user input, automatic generation.
- Support for triangular, quadrilateral, and mixed-element meshes.
- Finite Element Method (FEM)-based and geometric mesh smoothing and other topological quality-improvement functionality.
- Element quality evaluation (angular skewness) for quads & tris.
- Novel [layer-based conceptualization for 2D meshes.
  - [MeshLayers.pdf](https://github.com/user-attachments/files/19724245/MeshLayers.pdf)
  - Breadth-First Search x Image Skleletonization.
      - Note: "Mesh Layers" as a name for the algorithm is depracated; we will be migrating to using "Skeletonization".
- `.fort.14` file input/output for ADCIRC models
- API inspired by MATLAB’s `delaunayTriangulation()`

## Performance (v0.2.0)

CHILmesh 0.2.0 delivers **937x performance improvement** through systematic data structure modernization.

### Initialization Performance (WNAT_Hagen: 52,774 vertices, 98,365 elements)

| Operation | v0.1.1 | v0.2.0 | Speedup |
|-----------|--------|--------|---------|
| Fast init (no layers) | ~3,200s | **3.9s** | **822x** |
| Full init (30 layers) | ~5,400s | **7.7s** | **701x** |
| Quality analysis | ~4,800s | **6.6s** | **727x** |
| **Total** | ~13,400s (3.7 hrs) | **14.3s** | **937x** |

### Query Performance (O(1) guaranteed)

| Operation | v0.1.1 | v0.2.0 | Speedup |
|-----------|--------|--------|---------|
| Adjacency lookup | ~2,000μs | **4.0μs** | **500x** |
| Vertex neighbors | ~3,500μs | **0.7μs** | **5,000x** |
| Element neighbors | ~4,500μs | **4.4μs** | **1,022x** |
| 5,000 queries | ~6.8s | **0.022s** | **309x** |

### Real-World Impact

For MADMESHR mesh adaptation workflows:
- **Before:** 3,800s (64 minutes) per mesh
- **After:** 14.7s per mesh
- **Improvement:** **259x faster** — transforms hourly batch operations into interactive development

### Bulk Loading (ADMESH-Domains)

Fast initialization without layer computation:
- Average load time: **2.03ms per mesh** (requirement: <500ms) ✅
- Metadata queries: **0.044ms** (O(1) operation)

**See [BENCHMARK.md](BENCHMARK.md) for detailed methodology and breakdown.**

### Example Usage:
```python
# Load mesh
import matplotlib.pyplot as plt
import numpy as np
import chilmesh
import admesh  # Optional: required for rows 3 & 4 (ADMESH integration)

# Load one of the bundled example meshes (no file paths required).
# `annulus()` is a "raw" mesh: a fixed boundary + dispersed interior nodes,
# triangulated; quality is intentionally poor so the smoothers have work to do.
mesh = chilmesh.examples.annulus()
# Other built-ins: chilmesh.examples.donut(), .block_o(), .structured()

# Set up 4x3 subplot grid
fig, axs = plt.subplots(4, 3, figsize=(15, 18))
fig.suptitle("Original (dispersed) vs FEM-Smoothed vs "
             "ADMESH (boundary-only) vs ADMESH + Right-Angle Smoother",
             fontsize=16)

# --- Row 1: Original (dispersed) Mesh ---
_, ax = mesh.plot(ax=axs[0, 0])
mesh.plot_point(1, ax=ax); mesh.plot_edge(1, ax=ax); mesh.plot_elem(1, ax=ax)
ax.set_title("Original: Mesh + Highlighted Entities")

_, ax = mesh.plot_layer(ax=axs[0, 1])
ax.set_title("Original: Mesh Layers")

q0, _, _ = mesh.elem_quality()
_, ax = mesh.plot_quality(ax=axs[0, 2])
ax.set_title(f"Original: Quality (Median: {np.median(q0):.2f}, Std: {np.std(q0):.2f})")

# --- Row 2: CHILmesh FEM-Smoothed ---
mesh_fem = mesh.copy()
mesh_fem.smooth_mesh(method='fem', acknowledge_change=True)
_, ax = mesh_fem.plot(ax=axs[1, 0])
mesh_fem.plot_point(1, ax=ax); mesh_fem.plot_edge(1, ax=ax); mesh_fem.plot_elem(1, ax=ax)
ax.set_title("FEM-Smoothed: Mesh + Highlighted Entities")

_, ax = mesh_fem.plot_layer(ax=axs[1, 1])
ax.set_title("FEM-Smoothed: Mesh Layers")

q1, _, _ = mesh_fem.elem_quality()
_, ax = mesh_fem.plot_quality(ax=axs[1, 2])
ax.set_title(f"FEM-Smoothed: Quality (Median: {np.median(q1):.2f}, Std: {np.std(q1):.2f})")

# --- Row 3: ADMESH given row 1's BOUNDARY only (no interior nodes) ---
# Reuse row 1's annulus boundary; ADMESH inserts its own interior nodes.
def annulus_sdf(p):
    r = np.sqrt(p[:, 0]**2 + p[:, 1]**2)
    return np.maximum(r - 2.0, 1.0 - r)

boundary_pts = mesh.nodes[mesh.boundary_node_indices()]  # row 1's boundary
domain = admesh.Domain(annulus_sdf, bbox=(-2.5, -2.5, 2.5, 2.5), pfix=boundary_pts)
admesh_mesh = admesh.triangulate(domain, h_max=0.15)
mesh_adm = chilmesh.from_admesh(admesh_mesh)  # adapter: ADMESH -> CHILmesh
_, ax = mesh_adm.plot(ax=axs[2, 0])
mesh_adm.plot_point(1, ax=ax); mesh_adm.plot_edge(1, ax=ax); mesh_adm.plot_elem(1, ax=ax)
ax.set_title("ADMESH: Mesh + Highlighted Entities")

_, ax = mesh_adm.plot_layer(ax=axs[2, 1])
ax.set_title("ADMESH: Mesh Layers")

q2, _, _ = mesh_adm.elem_quality()
_, ax = mesh_adm.plot_quality(ax=axs[2, 2])
ax.set_title(f"ADMESH: Quality (Median: {np.median(q2):.2f}, Std: {np.std(q2):.2f})")

# --- Row 4: ADMESH + new right-angle triangle smoother ---
# Pushes triangles toward right-isosceles shape (preparation for quad fusion)
p_right, t_right = admesh.smooth_for_quadrangulation(
    admesh_mesh.nodes, admesh_mesh.elements,
    annulus_sdf, n_outer=3, pair_hint=True,
)
mesh_right = chilmesh.from_admesh_arrays(p_right, t_right)
_, ax = mesh_right.plot(ax=axs[3, 0])
mesh_right.plot_point(1, ax=ax); mesh_right.plot_edge(1, ax=ax); mesh_right.plot_elem(1, ax=ax)
ax.set_title("ADMESH + Right-Angle Smoother: Mesh + Highlighted Entities")

_, ax = mesh_right.plot_layer(ax=axs[3, 1])
ax.set_title("ADMESH + Right-Angle Smoother: Mesh Layers")

q3, _, _ = mesh_right.elem_quality()
_, ax = mesh_right.plot_quality(ax=axs[3, 2])
ax.set_title(f"ADMESH + Right-Angle Smoother: Quality "
             f"(Median: {np.median(q3):.2f}, Std: {np.std(q3):.2f})")

plt.tight_layout()
plt.subplots_adjust(top=0.96)
plt.show()
# fig.savefig("tests/output/annulus_quickstart.png", dpi=120, bbox_inches='tight')
```
![CHILmesh quickstart: annulus across four pipelines](https://raw.githubusercontent.com/domattioli/CHILmesh/main/tests/output/annulus_quickstart.png)

> **How the figure above was generated:** four meshes of the *same* annulus
> domain, each plotted with the three CHILmesh views (mesh + highlighted
> point/edge/element, BFS layers, quality map).
>
> - **Row 1 — Original (dispersed):** the bundled `chilmesh.examples.annulus()`
>   mesh — a fixed annular boundary plus interior nodes scattered roughly
>   uniformly inside, triangulated. Quality is intentionally poor.
> - **Row 2 — FEM-Smoothed:** CHILmesh's FEM smoother applied to row 1's
>   mesh (same connectivity, relaxed node positions).
> - **Row 3 — ADMESH (boundary-only):** ADMESH given *only* row 1's boundary
>   (no interior nodes) via `admesh.triangulate`; ADMESH inserts and triangulates
>   its own interior.
> - **Row 4 — ADMESH + Right-Angle Smoother:** ADMESH's new right-angle
>   triangle smoother (`admesh.smooth_for_quadrangulation`, an SVD-invariant
>   FEM Jacobian formulation that nudges triangles toward right-isosceles
>   shape in preparation for quad fusion) applied to row 3's mesh.
>
> The PNG is regenerated by `pytest tests/test_readme_quickstart.py`, which
> asserts geometric validity at every step (positive signed area, finite
> interior angles, complete layer cover, fort.14 roundtrip) before writing
> the image.


> **Note**: When mesh is mixed-element, connectivity (elem2vert adjacency) follows the format `Node1-Node2-Node3-Node4`, such that `Node4 == Node3` for triangular elements.


## See also / Downstream projects
[![MADMESHR_Project](https://img.shields.io/badge/GitHub-MADMESHR-121013?logo=github&logoColor=white&labelColor=gray)](https://github.com/domattioli/MADMESHR)

[MADMESHR](https://github.com/domattioli/MADMESHR) is a downstream research
project that builds on CHILmesh.


### BibTeX:
> DO Mattioli (2017). QuADMESH+: A Quadrangular ADvanced Mesh Generator for Hydrodynamic Models [Master's thesis, Ohio State University]. OhioLINK Electronic Theses and Dissertations Center. http://rave.ohiolink.edu/etdc/view?acc_num=osu1500627779532088
```bibtex
@mastersthesis{mattioli2017quadmesh,
  author       = {Mattioli, Dominik O.},
  title        = {{QuADMESH+}: A Quadrangular ADvanced Mesh Generator for Hydrodynamic Models},
  school       = {The Ohio State University},
  year         = {2017},
  note         = {Master's thesis},
  url          = {http://rave.ohiolink.edu/etdc/view?acc_num=osu1500627779532088}
}
```
- [Read the pdf for free here](https://github.com/user-attachments/files/19727573/QuADMESH__Thesis_Doc.pdf)



#### Acknowledgements
The following pieces of work inspired contributions to this repository:
- [ADMESH](https://doi.org/10.1007/s10236-012-0574-0)
- See the rest of the citations in the thesis [QuADMESH-Thesis.pdf](https://github.com/user-attachments/files/19724263/QuADMESH-Thesis.pdf)
- Original work was funded by [Aquaveo](https://aquaveo.com/) and contributed to by Alan Zundel.
- [FEM Smoother paper](https://api.semanticscholar.org/CorpusID:34335417)
  - [Inspiring MATLAB implementation](https://github.com/CHLNDDEV/OceanMesh2D/blob/Projection/utilities/direct_smoother_lur.m)
- [Angle-Based Smoother paper](https://www.andrew.cmu.edu/user/shimada/papers/00-imr-zhou.pdf)
  - The MATLAB code was originally developed for a master's thesis research project (2015–2017) at **The Ohio State University**.
