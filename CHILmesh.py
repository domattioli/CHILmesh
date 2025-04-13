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
        
        # Get edges and edge neighbors
        edge2vert = self.adjacencies["Edge2Vert"]
        edge2elem = self.adjacencies["Edge2Elem"].copy()  # Work on a copy to modify
        
        # Keep track of which elements have been assigned to a layer
        assigned_elements = set()
        
        layer_idx = 0
        
        # Continue until all elements are assigned to a layer or we can't find any more
        while len(assigned_elements) < self.n_elems:
            # Identify boundary edges of current layer
            if layer_idx == 0:
                # First layer - use all mesh boundaries (external and internal)
                boundary_mask = (edge2elem[:, 1] == 0)  # Second element is zero
                boundary_edges = np.where(boundary_mask)[0]
                
                # If no boundary edges found, break the loop
                if len(boundary_edges) == 0:
                    break
                    
                outer_vertices = np.unique(edge2vert[boundary_edges].flatten())
                self.layers["OV"].append(outer_vertices)
            else:
                # Subsequent layers - use edges with only one neighbor
                boundary_mask = np.sum(edge2elem > 0, axis=1) == 1
                boundary_edges = np.where(boundary_mask)[0]
                
                # If no boundary edges found, but we still have unassigned elements,
                # we might have isolated components - break the loop
                if len(boundary_edges) == 0:
                    break
                    
                outer_vertices = np.unique(edge2vert[boundary_edges].flatten())
                self.layers["OV"].append(outer_vertices)
            
            self.layers["bEdgeIDs"].append(boundary_edges)
            
            # Identify outer elements of layer (adjacent to boundary edges)
            outer_elems = []
            for edge_idx in boundary_edges:
                for elem_idx in edge2elem[edge_idx]:
                    if elem_idx > 0:  # Skip 0 (no element) and negative (processed)
                        outer_elems.append(elem_idx)
            
            if not outer_elems:  # If no outer elements found, break loop
                break
                
            outer_elems = np.unique(outer_elems)
            self.layers["OE"].append(outer_elems)
            assigned_elements.update(outer_elems)
            
            # Mark used edges in edge2elem
            for elem_idx in outer_elems:
                edge2elem[edge2elem == elem_idx] = -1  # Mark as processed
            
            # Identify all edges connected to outer vertices
            layer_edges = []
            for v in outer_vertices:
                for i, (v1, v2) in enumerate(edge2vert):
                    if v == v1 or v == v2:
                        layer_edges.append(i)
                        
            layer_edges = np.unique(layer_edges)
            
            # Identify inner elements of layer (using these layer edges)
            inner_elems = []
            for edge_idx in layer_edges:
                for elem_idx in edge2elem[edge_idx]:
                    if elem_idx > 0:  # Skip 0 (no element) and negative (processed)
                        inner_elems.append(elem_idx)
            
            inner_elems = np.unique(inner_elems)
            self.layers["IE"].append(inner_elems)
            assigned_elements.update(inner_elems)
            
            # Mark used edges in edge2elem
            for elem_idx in inner_elems:
                edge2elem[edge2elem == elem_idx] = -1  # Mark as processed
            
            # Identify inner vertices (vertices of outer and inner elements not in outer vertices)
            all_vertices = []
            
            # Combine outer and inner elements for this layer
            layer_elems = np.concatenate((outer_elems, inner_elems))
            
            # Make sure we have valid element indices
            for elem_idx in layer_elems:
                # Safety check - make sure index is valid integer
                if not isinstance(elem_idx, (int, np.integer)) or elem_idx < 0 or elem_idx >= self.n_elems:
                    continue
                    
                vertices = self.connectivity_list[elem_idx]
                for v in vertices:
                    if v != 0:  # Skip zero vertices
                        all_vertices.append(v)
            
            all_vertices = np.unique(all_vertices)
            inner_vertices = np.setdiff1d(all_vertices, outer_vertices)
            self.layers["IV"].append(inner_vertices)
            
            # Move to next layer
            layer_idx += 1
            
            # Break if we've processed all elements or are stuck
            if len(assigned_elements) >= self.n_elems or len(outer_elems) + len(inner_elems) == 0:
                break
        
        # Check if all elements were assigned
        if len(assigned_elements) < self.n_elems:
            print(f"Warning: Only {len(assigned_elements)} out of {self.n_elems} elements were assigned to layers")
        
        # Set number of layers
        self.n_layers = layer_idx


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
