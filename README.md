This class is a bi-product of a project, funded by Aquaveo, at The Ohio State University in 2015-2017. The home page for the Computational Hydrodynamics and Informatics Laboratory at OSU, and further information on the works of Dr. Ethan Kubatko's CHIL lab may be found [here](https://ceg.osu.edu/projects-software-computational-hydrodynamics-and-informatics-lab).

The repository focuses on a class for representing triangular, quadrangular, and mixed-element polygonal meshes for hydrodynamic domains. One of the primary techniques that this work leverages, called 'Mesh Layering' is illustrated below. Research related to the original work is ongoing.

![donut_triangulation_layers](https://github.com/user-attachments/assets/5911a256-7e91-4634-91d3-b4ddc5054bbd)
![donut_triangulation_layers_verts](https://github.com/user-attachments/assets/1da016bc-ef71-48bc-a563-f7d5022b98ed)
 

Connectivity representation for elements (faces) is Node1-Node2-Node3-Node4, where Node4==Node3 for triangular elements. The functionality for all methods are inspired by MATLAB's built-in 'delaunayTriangulation()' class.

[![View CHILmesh on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/135632-chilmesh)

#### Citation:
> Mattioli, D. D. (2017). QuADMESH+: A Quadrangular ADvanced Mesh Generator for Hydrodynamic Models [Master's thesis, Ohio State University]. OhioLINK Electronic Theses and Dissertations Center. http://rave.ohiolink.edu/etdc/view?acc_num=osu1500627779532088
