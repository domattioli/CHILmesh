<h1 align="center">
  CHILmesh: representing triangular, quadrangular and mixed-element (2D) meshes for advanced and automatic mesh generation for hydrodynamic domains.
</h1>

<p align="center">
  <strong><a href="#">Dominik Mattioli</a><sup>1†</sup>, <a href="#">Ethan Kubatko</a><sup>2</sup></strong><br>
  <sup>†</sup>Corresponding author<br><br>
  <sup>1</sup>Penn State University<br>
  <sup>2</sup>Computational Hydrodynamics and Informatics Lab (CHIL), The Ohio State University
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
- 2025/04/12 Python version of our code
- 2023/09/19 MATLAB code revisited; repo initiated
- 2017/08/01 Nascent MATLAB version of the code
#### Future Work
[![RL-QuADMESH_Project](https://img.shields.io/badge/GitHub-RL--QuADMESH-121013?logo=github&logoColor=white&labelColor=gray)](https://github.com/domattioli/RL-QuADMESH)


## Table of Contents
- [Releases](#releases)
- [Future Work](#future-work)
- [Installation](#installation)
- [Key Features](#key-features)
- [To-Do](#to-do)
- [Example Usage](#example-usage)
- [BibTeX](#bibtex)
- [Acknowledgements](#acknowledgements)


## Installation
```bash
git clone https://github.com/domattioli/CHILmesh && cd CHILmesh
python -m venv .myenv
source .myenv/bin/activate
.myenv/bin/pip install matplotlib numpy scipy
```


## Key Features
- Minimal user input, automatic generation.
- Support for triangular, quadrilateral, and mixed-element meshes.
- Finite Element Method (FEM)-based and geometric mesh smoothing and other topological quality-improvement functionality.
- Element quality evaluation (angular skewness) for quads & tris.
- Novel [layer-based conceptualization for 2D meshes.
  - [MeshLayers.pdf](https://github.com/user-attachments/files/19724245/MeshLayers.pdf)
- `.fort.14` file input/output for ADCIRC models
- API inspired by MATLAB’s `delaunayTriangulation()`

### To-Do📌
- Finish porting all functionality from original MATLAB code to python.
  - Add support for generating Delaunay meshes from scratch via zero-input CHILmesh().
  - Add support for [.gmsh](https://gmsh.info/doc/texinfo/gmsh.html) input/output.
  - Imperfections in the porting leading bugs, e.g., elements near 1st vertex being ignored in the plot() and plot_layer() methods but not plot_quality()
  - Extend .write_to_fort14() to support quadrilateral output

### Example Usage:
```python
# Load mesh
from CHILmesh import CHILmesh
mesh = CHILmesh.read_from_fort14("annulus_200pts.fort.14", grid_name="Annulus_200pts.fort.14")
# mesh = CHILmesh() # random delaunay to-do

# Set up 2x3 subplot grid
fig, axs = plt.subplots(2, 3, figsize=(18, 10))
axs = axs.flatten()
fig.suptitle("Original vs Smoothed Mesh Comparison", fontsize=16)

# --- Original Mesh Plots ---
# 0. Original: Mesh + point/edge/element
_, ax = mesh.plot(ax=axs[0])
mesh.plot_point(1, ax=ax)
mesh.plot_edge(1, ax=ax)
mesh.plot_elem(1, ax=ax)
ax.set_title("Original: Mesh + Highlighted Entities")

# 1. Original: Layers
_, ax = mesh.plot_layer(ax=axs[1])
ax.set_title("Original: Mesh Layers")

# 2. Original: Quality
q0, _, stats0 = mesh.elem_quality( )
print( stats0 )
_, ax = mesh.plot_quality(ax=axs[2])
ax.set_title(f"Original: Quality Map (Median: {np.median(q0):.2f}, Std: {np.std(q0):.2f})")

# --- Smoothed Mesh Plots ---
# 3. Smoothed: Mesh + point/edge/element
mesh_smoothed = mesh.copy()
mesh_smoothed.smooth_mesh( method='fem', acknowledge_change=True )
_, ax = mesh_smoothed.plot(ax=axs[3])
mesh_smoothed.plot_point(1, ax=ax)
mesh_smoothed.plot_edge(1, ax=ax)
mesh_smoothed.plot_elem(1, ax=ax)
ax.set_title("Smoothed: Mesh + Highlighted Entities")

# 4. Smoothed: Layers
_, ax = mesh_smoothed.plot_layer(ax=axs[4])
ax.set_title("Smoothed: Mesh Layers")

# 5. Smoothed: Quality
q, _, stats = mesh_smoothed.elem_quality( )
print( stats )
_, ax = mesh_smoothed.plot_quality(ax=axs[5])
ax.set_title(f"Smoothed: Quality Map (Median: {np.median(q):.2f}, Std: {np.std(q):.2f})")

# Layout tidy
plt.tight_layout()
plt.subplots_adjust(top=0.9)  # leave space for suptitle
plt.show()
```
![image](https://github.com/user-attachments/assets/7fdeda18-ad4a-4fd1-ad14-1f8cbc56e7da)

> **Note**: When mesh is mixed-element, connectivity (elem2vert adjacency) follows the format `Node1-Node2-Node3-Node4`, such that `Node4 == Node3` for triangular elements.


### BibTeX:
> Mattioli, D. D. (2017). QuADMESH+: A Quadrangular ADvanced Mesh Generator for Hydrodynamic Models [Master's thesis, Ohio State University]. OhioLINK Electronic Theses and Dissertations Center. http://rave.ohiolink.edu/etdc/view?acc_num=osu1500627779532088
```bibtex
@mastersthesis{mattioli2017quadmesh,
  author       = {Mattioli, Dominik D.},
  title        = {{QuADMESH+}: A Quadrangular ADvanced Mesh Generator for Hydrodynamic Models},
  school       = {The Ohio State University},
  year         = {2017},
  note         = {Master's thesis},
  url          = {http://rave.ohiolink.edu/etdc/view?acc_num=osu1500627779532088}
}
```
- [Read the pdf for free here](https://github.com/user-attachments/files/19724263/QuADMESH-Thesis.pdf)

#### Acknowledgements
The following pieces of work inspired contributions to this repository:
- [ADMESH](https://doi.org/10.1007/s10236-012-0574-0)
- See the rest of the citations in the thesis [QuADMESH-Thesis.pdf](https://github.com/user-attachments/files/19724263/QuADMESH-Thesis.pdf)
- Original work was funded by [Aquaveo](https://aquaveo.com/) and contributed to by Alan Zundel.
- FEM Smoother paper](https://api.semanticscholar.org/CorpusID:34335417)
  - [Inspiring MATLAB implementation](https://github.com/CHLNDDEV/OceanMesh2D/blob/Projection/utilities/direct_smoother_lur.m)
- [Angle-Based Smoother paper](https://www.andrew.cmu.edu/user/shimada/papers/00-imr-zhou.pdf)
  - The MATLAB code was originally developed for a master's thesis research project (2015–2017) at **The Ohio State University**.
