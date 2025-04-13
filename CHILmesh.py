from plot_utils import CHILmeshPlotMixin

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.spatial import Delaunay
from typing import List, Tuple, Optional, Dict, Set, Union, Any

class CHILmesh(CHILmeshPlotMixin):
    """
    A class for triangular, quadrilateral, or mixed-element meshes in 2D.
    
    This Python implementation is based on the MATLAB CHILmesh class from the 
    Computational Hydrodynamics & Informatics Laboratory (CHIL) at The Ohio State University,
    focusing on the mesh layers approach described in Mattioli's thesis.
    """
    
    @property
    def grid_name(self):
        """Grid name property."""
        return self._grid_name
    
    @grid_name.setter
    def grid_name(self, value):
        """Set grid name."""
        self._grid_name = value
    
    @property
    def Layers(self):
        """
        Backwards compatibility property to access layers with uppercase name.
        
        Returns:
            The layers dictionary
        """
        return self.layers
    

    def __init__( self, connectivity: Optional[np.ndarray] = None, 
                 points: Optional[np.ndarray] = None, 
                 grid_name: Optional[str] = None ) -> None:
        """
        Initialize a CHILmesh object.
        
        Parameters:
            connectivity: Element connectivity list
            points: Vertex coordinates
            grid_name: Name of the mesh
        """
        # Public properties
        self.grid_name = grid_name
        self.points = points
        self.connectivity_list = connectivity
        self.boundary_condition = None
        
        # Hidden properties
        self.adjacencies: Dict[str, Any] = {}
        self.n_verts: int = 0
        self.n_elems: int = 0
        self.n_edges: int = 0
        self.n_layers: int = 0
        self.layers: Dict[str, List] = {"OE": [], "IE": [], "OV": [], "IV": [], "bEdgeIDs": []}
        self.type: Optional[str] = None
        
        # If no inputs are provided, create a random Delaunay triangulation
        if connectivity is None and points is None:
            self._create_random_triangulation()
        
        # Initialize the mesh
        self._initialize_mesh()
    
    def _create_random_triangulation( self ) -> None:
        """Create a random Delaunay triangulation for testing"""
        # Generate random points in a 10x10 domain
        points = np.random.rand( 20, 2 ) * 10
        # Create Delaunay triangulation
        tri = Delaunay( points )
        # Store points and connectivity
        self.points = np.column_stack( ( tri.points, np.zeros( tri.points.shape[0] ) ) )
        self.connectivity_list = tri.simplices
        self.grid_name = "Random Delaunay"

    def _initialize_mesh( self ) -> None:
        """Initialize the mesh properties"""
        if self.points is not None and self.connectivity_list is not None:
            self.n_verts = self.points.shape[0]
            self.n_elems = self.connectivity_list.shape[0]
            
            # Ensure points have z-coordinate
            if self.points.shape[1] == 2:
                self.points = np.column_stack( ( self.points, np.zeros( self.n_verts ) ) )
            
            # Check connectivity orientation and correct if needed
            self._ensure_ccw_orientation()
            
            # Build adjacency lists
            self._build_adjacencies()
            
            # Identify mesh layers
            self._mesh_layers()
    
    def _ensure_ccw_orientation( self ) -> None:
        """Ensure counter-clockwise orientation of elements"""
        # Calculate signed area of each element
        areas = self.signed_area()
        
        # Find elements with clockwise orientation (negative area)
        cw_elements = np.where( areas < 0 )[0]
        
        # Flip orientation of clockwise elements
        for elem_id in cw_elements:
            # For triangular elements (3 vertices)
            if self.connectivity_list.shape[1] == 3:
                self.connectivity_list[elem_id] = self.connectivity_list[elem_id, [0, 2, 1]]
            # For quadrilateral elements (4 vertices)
            elif self.connectivity_list.shape[1] == 4:
                self.connectivity_list[elem_id] = self.connectivity_list[elem_id, [0, 3, 2, 1]]
    
    def signed_area( self, elem_ids: Optional[Union[int, List[int], np.ndarray]] = None ) -> np.ndarray:
        """
        Calculate the signed area of elements.
        
        Parameters:
            elem_ids: Indices of elements to evaluate.
                If None, all elements are evaluated.
        
        Returns:
            Signed areas of elements
        """
        if elem_ids is None:
            elem_ids = np.arange( self.n_elems )
        
        if np.isscalar( elem_ids ):
            elem_ids = [elem_ids]
        
        areas = np.zeros( len( elem_ids ) )
        
        # Determine element types
        tri_elems, quad_elems = self._elem_type( elem_ids )
        
        # Calculate areas for triangular elements
        for i, elem_id in enumerate( elem_ids ):
            if elem_id in tri_elems:
                vertices = self.connectivity_list[elem_id][:3]  # First 3 vertices for triangles
                x = self.points[vertices, 0]
                y = self.points[vertices, 1]
                # Shoelace formula for triangle
                areas[i] = 0.5 * ((x[0]*(y[1]-y[2]) + x[1]*(y[2]-y[0]) + x[2]*(y[0]-y[1])))
            elif elem_id in quad_elems:
                vertices = self.connectivity_list[elem_id]
                x = self.points[vertices, 0]
                y = self.points[vertices, 1]
                # Shoelace formula for quadrilateral
                areas[i] = 0.5 * ((x[0]*(y[1]-y[3]) + x[1]*(y[2]-y[0]) + 
                                  x[2]*(y[3]-y[1]) + x[3]*(y[0]-y[2])))
        
        return areas
    
    def _elem_type( self, elem_ids: Optional[Union[int, List[int], np.ndarray]] = None 
                  ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Identify triangular and quadrilateral elements.
        
        Parameters:
            elem_ids: Indices of elements to check.
                If None, all elements are checked.
        
        Returns:
            Tuple of (tri_elems, quad_elems) arrays of element indices
        """
        if elem_ids is None:
            elem_ids = np.arange( self.n_elems )
        
        # Convert to numpy array if not already
        elem_ids = np.array( elem_ids )
        
        if self.connectivity_list.shape[1] == 3:
            # All triangular elements
            return elem_ids, np.array( [] )
        
        # Check for triangular elements in a quad/mixed-element mesh
        tri_mask = np.zeros( len( elem_ids ), dtype=bool )
        
        for i, elem_id in enumerate( elem_ids ):
            # Check for redundant vertices or zero-valued vertices
            vertices = self.connectivity_list[elem_id]
            if (vertices[0] == vertices[1] or
                vertices[1] == vertices[2] or
                vertices[2] == vertices[3] or
                vertices[3] == vertices[0] or
                vertices[3] == 0):
                tri_mask[i] = True
        
        tri_elems = elem_ids[tri_mask]
        quad_elems = elem_ids[~tri_mask]
        
        return tri_elems, quad_elems
    
    def _build_adjacencies( self ) -> None:
        """Build adjacency lists for the mesh"""
        # Identify triangular and quadrilateral elements
        tri_elems, quad_elems = self._elem_type()
        
        # Set mesh type based on element types
        if len( quad_elems ) == 0:
            self.type = "Triangular"
        elif len( tri_elems ) == 0:
            self.type = "Quadrilateral"
        else:
            self.type = "Mixed-Element"
            # For mixed-element mesh, triangles may need adjustment for consistency
            if self.connectivity_list.shape[1] == 4:  # Already have space for 4 vertices
                for elem_id in tri_elems:
                    self.connectivity_list[elem_id, 3] = self.connectivity_list[elem_id, 0]
        
        # Identify edges of the mesh
        edges = self._identify_edges()
        edge2vert = np.array( edges )
        self.n_edges = len( edge2vert )
        
        # Build Elem2Edge
        elem2edge = self._build_elem2edge( edge2vert )
        
        # Build Vert2Edge
        vert2edge = self._build_vert2edge( edge2vert )
        
        # Build Vert2Elem
        vert2elem = self._build_vert2elem()
        
        # Build Edge2Elem
        edge2elem = self._build_edge2elem( edge2vert )
        
        # Store adjacencies
        self.adjacencies = {
            "Elem2Vert": self.connectivity_list,
            "Edge2Vert": edge2vert,
            "Elem2Edge": elem2edge,
            "Vert2Edge": vert2edge,
            "Vert2Elem": vert2elem,
            "Edge2Elem": edge2elem
        }

    def _identify_edges( self ) -> List[Tuple[int, int]]:
        """
        Identify edges of the mesh.
        
        Returns:
            List of edges, where each edge is a tuple of two vertex indices
        """
        edges = set()
        
        for elem_id in range( self.n_elems ):
            vertices = self.connectivity_list[elem_id]
            n_vertices = 3 if self.type == "Triangular" else 4
            
            for i in range( n_vertices ):
                v1 = vertices[i]
                v2 = vertices[(i+1) % n_vertices]
                
                # Skip invalid edges (vertices with value 0)
                if v1 == 0 or v2 == 0:
                    continue
                
                # Store edge as a sorted tuple to avoid duplicates
                edge = tuple( sorted( [int(v1), int(v2)] ) )
                edges.add( edge )
        
        return list( edges )
    
    def _build_elem2edge( self, edge2vert: np.ndarray ) -> np.ndarray:
        """
        Build Elem2Edge adjacency.
        
        Parameters:
            edge2vert: Edge-to-vertex adjacency
        
        Returns:
            Element-to-edge adjacency
        """
        max_edges_per_elem = 4 if self.type != "Triangular" else 3
        elem2edge = np.zeros( ( self.n_elems, max_edges_per_elem ), dtype=int )
        
        # For each element
        for elem_id in range( self.n_elems ):
            vertices = self.connectivity_list[elem_id]
            n_vertices = 3 if self.type == "Triangular" else 4
            
            # For each edge of the element
            for i in range( n_vertices ):
                v1 = vertices[i]
                v2 = vertices[(i+1) % n_vertices]
                
                # Skip invalid edges (vertices with value 0)
                if v1 == 0 or v2 == 0:
                    continue
                
                # Find the edge index
                edge = tuple( sorted( [int(v1), int(v2)] ) )
                for j, e in enumerate( edge2vert ):
                    if set( e ) == set( edge ):
                        elem2edge[elem_id, i] = j
                        break
        
        return elem2edge

    def _build_vert2edge( self, edge2vert: np.ndarray ) -> List[List[int]]:
        """
        Build Vert2Edge adjacency.
        
        Parameters:
            edge2vert: Edge-to-vertex adjacency
        
        Returns:
            Vertex-to-edge adjacency as list of lists
        """
        # Initialize with empty lists for each vertex
        vert2edge = [[] for _ in range( self.n_verts )]
        
        # Populate the lists
        for edge_id, (v1, v2) in enumerate( edge2vert ):
            vert2edge[v1].append( edge_id )
            vert2edge[v2].append( edge_id )
        
        return vert2edge
    
    def _build_vert2elem( self ) -> List[List[int]]:
        """
        Build Vert2Elem adjacency.
        
        Returns:
            Vertex-to-element adjacency as list of lists
        """
        # Initialize with empty lists for each vertex
        vert2elem = [[] for _ in range( self.n_verts )]
        
        # Populate the lists
        for elem_id in range( self.n_elems ):
            vertices = self.connectivity_list[elem_id]
            for v in vertices:
                if v != 0:  # Skip zero vertices
                    vert2elem[v].append( elem_id )
        
        return vert2elem
    
    def _build_edge2elem( self, edge2vert: np.ndarray ) -> np.ndarray:
        """
        Build Edge2Elem adjacency.
        
        Parameters:
            edge2vert: Edge-to-vertex adjacency
        
        Returns:
            Edge-to-element adjacency
        """
        # Initialize with zeros - at most 2 elements per edge
        edge2elem = np.zeros( ( self.n_edges, 2 ), dtype=int )
        
        # Iterate through elements and their edges
        for elem_id in range( self.n_elems ):
            vertices = self.connectivity_list[elem_id]
            n_vertices = 3 if self.type == "Triangular" else 4
            
            # For each edge of the element
            for i in range( n_vertices ):
                v1 = vertices[i]
                v2 = vertices[(i+1) % n_vertices]
                
                # Skip invalid edges (vertices with value 0)
                if v1 == 0 or v2 == 0:
                    continue
                
                # Find the edge index
                edge = tuple( sorted( [int(v1), int(v2)] ) )
                for edge_id, e in enumerate( edge2vert ):
                    if set( e ) == set( edge ):
                        # Check if edge already has an element assigned
                        if edge2elem[edge_id, 0] == 0:
                            edge2elem[edge_id, 0] = elem_id + 1  # +1 to avoid 0
                        else:
                            edge2elem[edge_id, 1] = elem_id + 1
                        break
        
        # Adjust indices (remove the +1 offset)
        edge2elem[edge2elem > 0] -= 1
        
        return edge2elem
    
    def boundary_edges( self ) -> np.ndarray:
        """
        Identify boundary edges of the mesh.
        
        Returns:
            Indices of boundary edges
        """
        # Boundary edges have only one adjacent element
        edge2elem = self.adjacencies["Edge2Elem"]
        boundary_mask = (edge2elem[:, 1] == 0)  # Second element is zero
        return np.where( boundary_mask )[0]
    
    def edge2vert( self, edge_ids: Optional[Union[int, List[int], np.ndarray]] = None ) -> np.ndarray:
        """
        Get vertices that define specified edges.
        
        Parameters:
            edge_ids: Indices of edges to query.
                If None, all edges are queried.
        
        Returns:
            Array of vertex indices for each edge
        """
        if edge_ids is None:
            edge_ids = np.arange( self.n_edges )
        
        if np.isscalar( edge_ids ):
            edge_ids = [edge_ids]
        
        return self.adjacencies["Edge2Vert"][edge_ids]
    
    def elem2edge( self, elem_ids: Optional[Union[int, List[int], np.ndarray]] = None ) -> np.ndarray:
        """
        Get edges that define specified elements.
        
        Parameters:
            elem_ids: Indices of elements to query.
                If None, all elements are queried.
        
        Returns:
            Array of edge indices for each element
        """
        if elem_ids is None:
            elem_ids = np.arange( self.n_elems )
        
        if np.isscalar( elem_ids ):
            elem_ids = [elem_ids]
        
        return self.adjacencies["Elem2Edge"][elem_ids]
    
    def edge2elem( self, edge_ids: Optional[Union[int, List[int], np.ndarray]] = None ) -> np.ndarray:
        """
        Get elements adjacent to specified edges.
        
        Parameters:
            edge_ids: Indices of edges to query.
                If None, all edges are queried.
        
        Returns:
            Array of element indices for each edge
        """
        if edge_ids is None:
            edge_ids = np.arange( self.n_edges )
        
        if np.isscalar( edge_ids ):
            edge_ids = [edge_ids]
        
        return self.adjacencies["Edge2Elem"][edge_ids]
    
    def _mesh_layers(self) -> None:
        """
        Discretize the mesh into layers starting from the boundary.
        This implements the mesh layers approach described in Mattioli's thesis.
        """
        # Reset layers
        self.layers = {"OE": [], "IE": [], "OV": [], "IV": [], "bEdgeIDs": []}
        
        # Get boundary edges (edges with only one adjacent element)
        edge2elem = self.adjacencies["Edge2Elem"]
        boundary_mask = (edge2elem[:, 1] == 0)
        boundary_edges = np.where(boundary_mask)[0]
        
        # Keep track of which elements have been assigned to a layer
        remaining_elements = set(range(self.n_elems))
        
        # Get element-to-element connectivity using edge2elem
        elem2elem = [[] for _ in range(self.n_elems)]
        for edge_idx, (e1, e2) in enumerate(edge2elem):
            if e1 >= 0 and e2 >= 0:  # Both elements exist
                elem2elem[e1].append(e2)
                elem2elem[e2].append(e1)
        
        # Process layers from the boundary inward
        layer_idx = 0
        while remaining_elements and len(boundary_edges) > 0:
            # Get boundary vertices
            edge2vert = self.adjacencies["Edge2Vert"]
            outer_vertices = np.array(list(set(edge2vert[boundary_edges].flatten())))
            self.layers["OV"].append(outer_vertices)
            self.layers["bEdgeIDs"].append(boundary_edges)
            
            # Get outer elements (elements adjacent to boundary edges)
            outer_elems = []
            for edge_idx in boundary_edges:
                elems = edge2elem[edge_idx]
                for elem in elems:
                    if elem >= 0 and elem in remaining_elements:
                        outer_elems.append(elem)
            
            # Convert to numpy array of integers
            outer_elems = np.array(list(set(outer_elems)), dtype=int)
            
            # Skip if no outer elements found
            if len(outer_elems) == 0:
                break
                
            self.layers["OE"].append(outer_elems)
            for elem in outer_elems:
                remaining_elements.remove(elem)
            
            # Get inner elements (neighbors of outer elements that haven't been assigned yet)
            inner_elems = []
            for elem in outer_elems:
                for neighbor in elem2elem[elem]:
                    if neighbor in remaining_elements:
                        inner_elems.append(neighbor)
            
            # Convert to numpy array of integers and remove duplicates
            inner_elems = np.array(list(set(inner_elems)), dtype=int)
            
            # Store inner elements
            self.layers["IE"].append(inner_elems)
            for elem in inner_elems:
                if elem in remaining_elements:
                    remaining_elements.remove(elem)
            
            # Get inner vertices
            all_vertices = set()
            for elem in np.concatenate((outer_elems, inner_elems)):
                vertices = self.connectivity_list[elem]
                for v in vertices:
                    if v > 0:  # Skip zero vertices
                        all_vertices.add(v)
            
            inner_vertices = np.array(list(all_vertices - set(outer_vertices)), dtype=int)
            self.layers["IV"].append(inner_vertices)
            
            # Get new boundary by finding edges that have one element in the remaining set
            # and one element in the processed set
            boundary_edges = []
            for edge_idx, (e1, e2) in enumerate(edge2elem):
                # Skip boundary edges of the original mesh
                if e2 < 0:
                    continue
                    
                # An edge is a boundary if exactly one of its adjacent elements 
                # is in the remaining set
                if ((e1 in remaining_elements) != (e2 in remaining_elements)):
                    boundary_edges.append(edge_idx)
            
            boundary_edges = np.array(boundary_edges, dtype=int)
            
            # Move to next layer
            layer_idx += 1
        
        # Set number of layers
        self.n_layers = layer_idx
        
        # # Print summary
        # print(f"Created {self.n_layers} mesh layers")
        # for i in range(self.n_layers):
        #     print(f"  Layer {i}: {len(self.layers['OE'][i])} outer elements, {len(self.layers['IE'][i])} inner elements")

    def get_layer( self, layer_idx: int ) -> Dict[str, np.ndarray]:
        """
        Get the components of a specific mesh layer.
        
        Parameters:
            layer_idx: Index of the layer to retrieve
        
        Returns:
            Dictionary with outer elements (OE), inner elements (IE),
            outer vertices (OV), and inner vertices (IV) of the layer
        """
        if layer_idx < 0 or layer_idx >= self.n_layers:
            raise ValueError( f"Layer index {layer_idx} out of range [0, {self.n_layers-1}]" )
        
        return {
            "OE": self.layers["OE"][layer_idx],
            "IE": self.layers["IE"][layer_idx],
            "OV": self.layers["OV"][layer_idx],
            "IV": self.layers["IV"][layer_idx],
            "bEdgeIDs": self.layers["bEdgeIDs"][layer_idx]
        }
    
    
    @staticmethod
    def from_fort14(filepath: str, grid_name: str = None) -> "CHILmesh":
        """
        Load a mesh from a FORT.14 file.
        
        Parameters:
            filepath: Path to the FORT.14 file
            grid_name: Optional name for the grid
            
        Returns:
            A CHILmesh object
        """
        with open(filepath, 'r') as f:
            # Read header
            header = f.readline().strip()
            
            # Read element and node counts
            counts = f.readline().strip().split()
            n_elements = int(counts[0])
            n_nodes = int(counts[1])
            
            # Read nodes
            points = np.zeros((n_nodes, 3))  # x, y, z
            for i in range(n_nodes):
                line = f.readline().strip().split()
                # Fort14 format: node_idx x y z
                # Assuming the node index is 1-based
                points[i] = [float(line[1]), float(line[2]), float(line[3])]
            
            # Read elements
            elements = np.zeros((n_elements, 3), dtype=int)
            for i in range(n_elements):
                line = f.readline().strip().split()
                # Fort14 format: elem_idx num_nodes node1 node2 node3
                # Skip the element index, look at the number of nodes
                num_nodes = int(line[1])
                
                if num_nodes != 3:
                    raise ValueError(f"Only triangular elements supported, found element with {num_nodes} nodes")
                    
                # Convert 1-based node indices to 0-based
                node_indices = [int(line[j+2]) - 1 for j in range(num_nodes)]
                elements[i] = node_indices
        
        return CHILmesh(connectivity=elements, points=points, grid_name=grid_name or header)

    def interior_angles(self, elem_ids=None) -> np.ndarray:
        """
        Calculate interior angles of mesh elements.
        
        Parameters:
            elem_ids: Indices of elements to evaluate.
                If None, all elements are evaluated.
        
        Returns:
            Array of interior angles for each element
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        
        if np.isscalar(elem_ids):
            elem_ids = [elem_ids]
        
        # Determine element types
        tri_elems, quad_elems = self._elem_type(elem_ids)
        
        # Maximum number of angles per element
        max_angles = 4 if len(quad_elems) > 0 else 3
        
        # Initialize angles array
        angles = np.zeros((len(elem_ids), max_angles))
        
        # Calculate angles for each element
        for i, elem_id in enumerate(elem_ids):
            if elem_id in tri_elems:
                # Triangle angles
                vertices = self.connectivity_list[elem_id][:3]  # First 3 vertices for triangles
                coords = self.points[vertices, :2]  # Get x,y coordinates
                
                # Calculate angles at each vertex
                for j in range(3):
                    v1 = coords[(j+1)%3] - coords[j]
                    v2 = coords[(j-1)%3] - coords[j]
                    
                    # Normalize vectors
                    v1_norm = v1 / np.linalg.norm(v1)
                    v2_norm = v2 / np.linalg.norm(v2)
                    
                    # Calculate angle in degrees
                    dot_product = np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0)
                    angle = np.arccos(dot_product) * 180 / np.pi
                    angles[i, j] = angle
                    
            elif elem_id in quad_elems:
                # Quadrilateral angles
                vertices = self.connectivity_list[elem_id]  # All 4 vertices
                coords = self.points[vertices, :2]  # Get x,y coordinates
                
                # Calculate angles at each vertex
                for j in range(4):
                    v1 = coords[(j+1)%4] - coords[j]
                    v2 = coords[(j-1)%4] - coords[j]
                    
                    # Normalize vectors
                    v1_norm = v1 / np.linalg.norm(v1)
                    v2_norm = v2 / np.linalg.norm(v2)
                    
                    # Calculate angle in degrees
                    dot_product = np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0)
                    angle = np.arccos(dot_product) * 180 / np.pi
                    angles[i, j] = angle
        return angles

    def elem_quality(self, elem_ids=None, quality_type='skew') -> Tuple[np.ndarray, np.ndarray, dict]:
        """
        Calculate the quality of mesh elements.
        
        Parameters:
            elem_ids: Indices of elements to evaluate.
                If None, all elements are evaluated.
            quality_type: Type of quality metric to use.
                'skew', 'skewness', 'angular skewness': Measures deviation from ideal angles
        
        Returns:
            Tuple of (Quality, Angles) where:
            - Quality: Array of quality measurements for each element
            - Angles: Array of interior angles for each element
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        
        if np.isscalar(elem_ids):
            elem_ids = [elem_ids]
        
        # Determine element types
        tri_elems, quad_elems = self._elem_type(elem_ids)
        
        # Calculate interior angles
        angles = self.interior_angles(elem_ids)
        
        # Initialize quality array
        quality = np.zeros(len(elem_ids))
        
        # Compute quality based on the selected metric
        if quality_type in ['skew', 'skewness', 'angular skewness']:
            # Process triangular elements
            tri_mask = np.array([elem_id in tri_elems for elem_id in elem_ids])
            if np.any(tri_mask):
                # Get angles for triangular elements
                tri_angles = angles[tri_mask, :3]
                
                # Calculate max and min angles
                tri_max = np.max(tri_angles, axis=1)
                tri_min = np.min(tri_angles, axis=1)
                
                # Equiangular skew for triangles (ideal angle = 60°)
                quality[tri_mask] = 1 - np.maximum(
                    (tri_max - 60) / (180 - 60),
                    (60 - tri_min) / 60
                )
            
            # Process quadrilateral elements
            quad_mask = np.array([elem_id in quad_elems for elem_id in elem_ids])
            if np.any(quad_mask):
                # Get angles for quadrilateral elements
                quad_angles = angles[quad_mask, :]
                
                # Calculate max and min angles
                quad_max = np.max(quad_angles, axis=1)
                quad_min = np.min(quad_angles, axis=1)
                
                # Equiangular skew for quads (ideal angle = 90°)
                quality[quad_mask] = 1 - np.maximum(
                    (quad_max - 90) / (180 - 90),
                    (90 - quad_min) / 90
                )
            
            # Handle poor angle calculations (concave elements, etc.)
            # For triangles, sum of angles should be close to 180°
            tri_sum_mask = tri_mask & (np.sum(angles[:, :3], axis=1) <= 179.99)
            quality[tri_sum_mask] = 0
            
            # For quads, sum of angles should be close to 360°
            quad_sum_mask = quad_mask & (np.sum(angles, axis=1) <= 359.99)
            quality[quad_sum_mask] = 0
        
        else:
            raise ValueError(f"Unknown quality type: {quality_type}")
        
        # Calculate statistics for the computed quality
        stats = {
            'mean': np.mean(quality),
            'median': np.median(quality),
            'min': np.min(quality),
            'max': np.max(quality),
            'std': np.std(quality)
        }
        return quality, angles, stats
    
