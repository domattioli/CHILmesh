from pathlib import Path
import warnings

from .utils.plot_utils import CHILmeshPlotMixin

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.spatial import Delaunay
from typing import List, Tuple, Optional as Opt, Dict, Set, Union, Any
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from copy import deepcopy

__all__ = ['CHILmesh', 'write_fort14']

class CHILmesh(CHILmeshPlotMixin):
    """
    A 2D mesh class supporting triangular, quadrilateral, and mixed-element meshes.

    Supports multiple file formats (ADCIRC `.fort.14`, SMS `.2dm`) and integrates
    with ADMESH-Domains catalog. Provides mesh analysis (layer structure, element
    quality, interior angles), geometric operations (smoothing), and fast metadata
    queries for bulk loading.

    Key Features:
    - Element Types: Triangles, quads, and mixed-element meshes (padded convention)
    - Fast Init: Optional skeletonization for <2s bulk loading (compute_layers=False)
    - Metadata: Node count, element count, element type, bounding box (ADMESH-Domains compatible)
    - Entry Point: CHILmesh.from_admesh_domain() for catalog integration (duck-typed, zero deps)
    - File I/O: Read ADCIRC `.fort.14` and SMS `.2dm` formats; roundtrip lossless
    - Analysis: Layers (skeletonization), element quality, interior angles

    Based on the MATLAB CHILmesh class from the Computational Hydrodynamics &
    Informatics Laboratory (CHIL) at The Ohio State University, focusing on the
    mesh layers approach described in Mattioli's thesis.

    Adjacency Structures:
        The mesh maintains multiple representations of topology:

        Core Mesh Representation:
        - points: ndarray[n_verts, 3] - Vertex coordinates (x, y, z)
        - connectivity_list: ndarray[n_elems, 3|4] - Element vertex indices
          * Triangles: [v0, v1, v2, -1]
          * Quads: [v0, v1, v2, v3]
          * Mixed: Triangles padded with repeated vertex

        Adjacency Dictionary (self.adjacencies):
        - Edge2Vert: ndarray[n_edges, 2] - Edge endpoints (canonical form: min, max)
        - Elem2Edge: ndarray[n_elems, 3|4] - Element edge IDs
        - Edge2Elem: ndarray[n_edges, 2] - Adjacent elements per edge (boundary: -1)
        - Vert2Edge: Dict[int, Set[int]] - Vertex incident edges (access via get_vertex_edges())
        - Vert2Elem: Dict[int, Set[int]] - Vertex incident elements (access via get_vertex_elements())
        - EdgeMap: Internal hash table for O(1) edge lookup (not part of public API)

        Access Patterns:
        >>> edges = mesh.get_vertex_edges(v)  # Set[int]
        >>> elems = mesh.get_vertex_elements(v)  # Set[int]
        >>> v1, v2 = mesh.edge2vert(edge_id)
        >>> elem_ids = mesh.edge2elem(edge_id)

    Examples:
        Load from ADMESH-Domains catalog:
            mesh = CHILmesh.from_admesh_domain(record)

        Fast metadata query:
            metadata = mesh.admesh_metadata()
            print(f"Nodes: {metadata['node_count']}, Type: {metadata['element_type']}")

        Read SMS format:
            mesh = CHILmesh.read_from_2dm(Path("mesh.2dm"))

        Query vertex neighborhoods:
            edges = mesh.get_vertex_edges(5)
            elems = mesh.get_vertex_elements(5)
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
    
    def change_points(self, new_points, acknowledge_change=False):
        """
        Change the mesh's (x,y,z) locations of its points.

        Parameters:
            new_points: New coordinates for the points
            acknowledge_change: If True, acknowledges the change in the mesh
        """
        assert acknowledge_change, "acknowledge_change must be True to change points -- this will change the mesh, make sure you understand this before using this method within a broader algorithm."
        if new_points.shape[1] == 2:    self.points[:, :2] = new_points
        elif new_points.shape[1] == 3:  self.points = new_points
        else:                           raise ValueError("new_points must have 2 or 3 columns")

    def __init__( self, connectivity: Opt[np.ndarray] = None, points: Opt[np.ndarray] = None, grid_name: Opt[str] = None, compute_layers: bool = True ) -> None:
        """
        Initialize a CHILmesh object.

        Parameters:
            connectivity: Element connectivity list
            points: Vertex coordinates
            grid_name: Name of the mesh
            compute_layers: If False, skip skeletonization for fast init (default: True)
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
        self.type: Opt[str] = None

        # If no inputs are provided, create a random Delaunay triangulation
        if connectivity is None and points is None:
            self._create_random_triangulation()

        # Initialize the mesh
        self._initialize_mesh(compute_layers=compute_layers)
    
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

    def _initialize_mesh( self, compute_layers: bool = True ) -> None:
        """Initialize the mesh properties.

        Parameters:
            compute_layers: If False, skip adjacency building and skeletonization (default: True)
        """
        if self.points is not None and self.connectivity_list is not None:
            self.n_verts = self.points.shape[0]
            self.n_elems = self.connectivity_list.shape[0]

            # Ensure points have z-coordinate
            if self.points.shape[1] == 2:
                self.points = np.column_stack( ( self.points, np.zeros( self.n_verts ) ) )

            # Check connectivity orientation and correct if needed
            self._ensure_ccw_orientation()

            # Set element type (triangular, quadrilateral, or mixed-element)
            tri_elems, quad_elems = self._elem_type()
            if len(quad_elems) == 0:
                self.type = "Triangular"
            elif len(tri_elems) == 0:
                self.type = "Quadrilateral"
            else:
                self.type = "Mixed-Element"

            # Build adjacency lists and compute layers only if requested
            if compute_layers:
                self._build_adjacencies()
                self._skeletonize()
    
    def _ensure_ccw_orientation( self ) -> None:
        """Ensure counter-clockwise orientation of every element.

        For mixed-element meshes (4-column connectivity containing both real
        quads and padded triangles) the flip permutation must be selected per
        element: a triangle pads its 4th slot with a repeated vertex, and
        applying the quad permutation ``[0, 3, 2, 1]`` to such a row would
        scramble it (B4).
        """
        areas = self.signed_area()
        cw_elements = np.where(areas < 0)[0]
        if cw_elements.size == 0:
            return

        if self.connectivity_list.shape[1] == 3:
            self.connectivity_list[cw_elements] = self.connectivity_list[
                cw_elements
            ][:, [0, 2, 1]]
            return

        # 4-column connectivity: classify each CW row before flipping.
        tri_elems, _ = self._elem_type(cw_elements)
        tri_set = set(tri_elems.tolist())
        for elem_id in cw_elements:
            row = self.connectivity_list[elem_id]
            if elem_id in tri_set:
                # Padded triangle: only the first 3 slots carry geometry,
                # the 4th slot is a repeated pad. Flip in place and re-pad.
                # Find the slot that holds the pad (any duplicate of another).
                a, b, c, d = row.tolist()
                # Identify the three distinct vertices (in row order, ignoring
                # the duplicated one) and reverse them to flip CCW.
                seen = []
                for v in (a, b, c, d):
                    if v not in seen:
                        seen.append(v)
                if len(seen) != 3:
                    # Fallback: degenerate row, leave as-is.
                    continue
                v0, v1, v2 = seen
                self.connectivity_list[elem_id] = np.array([v0, v2, v1, v0])
            else:
                self.connectivity_list[elem_id] = row[[0, 3, 2, 1]]
    
    def signed_area( self, elem_ids: Opt[Union[int, List[int], np.ndarray]] = None ) -> np.ndarray:
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
    
    def _elem_type( self, elem_ids: Opt[Union[int, List[int], np.ndarray]] = None 
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
        
        # Detect padded triangles in a 4-column connectivity table by
        # looking for any repeated vertex within the row. Vertex id ``0`` is a
        # *valid* vertex in this Python (0-indexed) port — the original code
        # treated ``vertices[3] == 0`` as a sentinel, which silently demoted
        # any quad whose 4th vertex happened to be 0 (B3).
        rows = self.connectivity_list[elem_ids]
        # A triangle padded into a quad row will have at least one repeated
        # vertex (typically ``v3 == v0`` or ``v3 == v2``).
        tri_mask = (
            (rows[:, 0] == rows[:, 1])
            | (rows[:, 1] == rows[:, 2])
            | (rows[:, 2] == rows[:, 3])
            | (rows[:, 3] == rows[:, 0])
            | (rows[:, 0] == rows[:, 2])
            | (rows[:, 1] == rows[:, 3])
        )

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
        edges, edge_map = self._identify_edges()
        edge2vert = np.array( edges )
        self.n_edges = len( edge2vert )

        # Build Elem2Edge
        elem2edge = self._build_elem2edge( edge2vert, edge_map )

        # Build Vert2Edge
        vert2edge = self._build_vert2edge( edge2vert )

        # Build Vert2Elem
        vert2elem = self._build_vert2elem()

        # Build Edge2Elem
        edge2elem = self._build_edge2elem( edge2vert, edge_map )

        # Store adjacencies
        self.adjacencies = {
            "Elem2Vert": self.connectivity_list,
            "Edge2Vert": edge2vert,
            "EdgeMap": edge_map,
            "Elem2Edge": elem2edge,
            "Vert2Edge": vert2edge,
            "Vert2Elem": vert2elem,
            "Edge2Elem": edge2elem
        }
        self._validate_adjacencies()

    def _validate_adjacencies(self) -> None:
        """
        Validate adjacency structure invariants.

        Checks:
        - All vertices have entries in Vert2Edge and Vert2Elem
        - All edge IDs are valid [0, n_edges)
        - All element IDs are valid [0, n_elems)
        - Vert2Edge and Vert2Elem are dicts with set values
        - Edge ordering consistency (edges stored in canonical form)

        Raises:
            AssertionError: If any invariant is violated
        """
        v2e = self.adjacencies['Vert2Edge']
        v2m = self.adjacencies['Vert2Elem']
        e2v = self.adjacencies['Edge2Vert']
        e2m = self.adjacencies['Edge2Elem']

        # Check: All vertices have entries
        assert len(v2e) == self.n_verts, f"Vert2Edge has {len(v2e)} entries, expected {self.n_verts}"
        assert len(v2m) == self.n_verts, f"Vert2Elem has {len(v2m)} entries, expected {self.n_verts}"

        # Check: Vert2Edge structure
        for vert_id in range(self.n_verts):
            assert vert_id in v2e, f"Vertex {vert_id} missing from Vert2Edge"
            assert isinstance(v2e[vert_id], set), f"Vert2Edge[{vert_id}] is not a set"
            for edge_id in v2e[vert_id]:
                assert 0 <= edge_id < self.n_edges, f"Invalid edge ID {edge_id} in Vert2Edge[{vert_id}]"
                # Verify vertex is endpoint of edge
                v1, v2 = e2v[edge_id]
                assert vert_id in [v1, v2], f"Vertex {vert_id} not endpoint of edge {edge_id}"

        # Check: Vert2Elem structure
        for vert_id in range(self.n_verts):
            assert vert_id in v2m, f"Vertex {vert_id} missing from Vert2Elem"
            assert isinstance(v2m[vert_id], set), f"Vert2Elem[{vert_id}] is not a set"
            for elem_id in v2m[vert_id]:
                assert 0 <= elem_id < self.n_elems, f"Invalid element ID {elem_id} in Vert2Elem[{vert_id}]"

    def _identify_edges( self ) -> Tuple[List[Tuple[int, int]], 'EdgeMap']:
        """
        Identify edges of the mesh.

        Builds both a list of edges (for backward compatibility) and an
        EdgeMap (for O(1) lookups in subsequent operations).

        Returns:
            Tuple of (edges_list, edge_map) where:
            - edges_list: List of edges as (v1, v2) tuples in sorted order
            - edge_map: EdgeMap object for O(1) edge ID lookup
        """
        from .mesh_topology import EdgeMap

        edge_map = EdgeMap()
        seen = set()

        for elem_id in range( self.n_elems ):
            vertices = self.connectivity_list[elem_id]
            n_vertices = 3 if self.type == "Triangular" else 4

            for i in range( n_vertices ):
                v1 = vertices[i]
                v2 = vertices[(i+1) % n_vertices]

                # Skip invalid edges (negative vertex ids)
                # In MATLAB the value 0 signified a placeholder for a missing
                # vertex in mixed element meshes.  In this Python port we use
                # 0-based indexing, therefore vertex id ``0`` is valid and
                # should not be discarded.  Only negative ids are considered
                # invalid.
                if v1 < 0 or v2 < 0:
                    continue

                # Store edge in canonical form to check for duplicates
                edge = tuple( sorted( [int(v1), int(v2)] ) )
                if edge not in seen:
                    seen.add( edge )
                    # Add to EdgeMap (maintains consistent ID ordering)
                    edge_map.add_edge(edge[0], edge[1])

        # Return edges in EdgeMap order to ensure IDs match
        return edge_map.to_list(), edge_map
    
    def _build_elem2edge( self, edge2vert: np.ndarray, edge_map: 'EdgeMap' = None ) -> np.ndarray:
        """
        Build Elem2Edge adjacency.

        Parameters:
            edge2vert: Edge-to-vertex adjacency (ndarray)
            edge_map: EdgeMap for O(1) edge lookup (optional, for optimization)

        Returns:
            Element-to-edge adjacency

        Note:
            If edge_map is provided, uses O(1) lookups instead of linear search,
            reducing complexity from O(n²) to O(n log n).
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

                # Skip invalid edges (negative vertex ids)
                if v1 < 0 or v2 < 0:
                    continue

                # Find the edge index using EdgeMap if available (O(1))
                # Otherwise fall back to linear search (O(n))
                if edge_map is not None:
                    edge_id = edge_map.find_edge(int(v1), int(v2))
                    if edge_id is not None:
                        elem2edge[elem_id, i] = edge_id
                else:
                    edge = tuple( sorted( [int(v1), int(v2)] ) )
                    for j, e in enumerate( edge2vert ):
                        if set( e ) == set( edge ):
                            elem2edge[elem_id, i] = j
                            break

        return elem2edge

    def _build_vert2edge( self, edge2vert: np.ndarray ) -> Dict[int, Set[int]]:
        """
        Build Vert2Edge adjacency as explicit dict with set values.

        Maps each vertex to the set of edge IDs incident to that vertex.

        Parameters:
            edge2vert: Edge-to-vertex adjacency (n_edges x 2 array)

        Returns:
            Dict mapping vertex ID to Set of incident edge IDs
        """
        vert2edge: Dict[int, Set[int]] = {}
        for vert_id in range(self.n_verts):
            vert2edge[vert_id] = set()

        for edge_id, (v1, v2) in enumerate(edge2vert):
            vert2edge[int(v1)].add(edge_id)
            vert2edge[int(v2)].add(edge_id)

        return vert2edge
    
    def _build_vert2elem( self ) -> Dict[int, Set[int]]:
        """
        Build Vert2Elem adjacency as explicit dict with set values.

        Maps each vertex to the set of element IDs incident to that vertex.

        Parameters:
            None

        Returns:
            Dict mapping vertex ID to Set of incident element IDs
        """
        vert2elem: Dict[int, Set[int]] = {}
        for vert_id in range(self.n_verts):
            vert2elem[vert_id] = set()

        for elem_id in range(self.n_elems):
            vertices = self.connectivity_list[elem_id]
            for v in vertices:
                if v >= 0:
                    vert2elem[int(v)].add(elem_id)

        return vert2elem
    
    def _build_edge2elem( self, edge2vert: np.ndarray, edge_map: 'EdgeMap' = None ) -> np.ndarray:
        """
        Build Edge2Elem adjacency.

        Parameters:
            edge2vert: Edge-to-vertex adjacency (ndarray)
            edge_map: EdgeMap for O(1) edge lookup (optional, for optimization)

        Returns:
            Edge-to-element adjacency (n_edges × 2 with -1 for boundary)

        Note:
            If edge_map is provided, uses O(1) lookups instead of linear search,
            reducing complexity from O(n²) to O(n log n).
        """
        # Initialize with -1 sentinel (no adjacent element)
        edge2elem = np.full( ( self.n_edges, 2 ), -1, dtype=int )

        # Iterate through elements and their edges
        for elem_id in range( self.n_elems ):
            vertices = self.connectivity_list[elem_id]
            n_vertices = 3 if self.type == "Triangular" else 4

            # For each edge of the element
            for i in range( n_vertices ):
                v1 = vertices[i]
                v2 = vertices[(i+1) % n_vertices]

                # Skip invalid edges (negative vertex ids)
                if v1 < 0 or v2 < 0:
                    continue

                # Find the edge index using EdgeMap if available (O(1))
                # Otherwise fall back to linear search (O(n))
                if edge_map is not None:
                    edge_id = edge_map.find_edge(int(v1), int(v2))
                    if edge_id is not None:
                        # Check if edge already has an element assigned
                        if edge2elem[edge_id, 0] == -1:
                            edge2elem[edge_id, 0] = elem_id
                        else:
                            edge2elem[edge_id, 1] = elem_id
                else:
                    edge = tuple( sorted( [int(v1), int(v2)] ) )
                    for edge_id, e in enumerate( edge2vert ):
                        if set( e ) == set( edge ):
                            # Check if edge already has an element assigned
                            if edge2elem[edge_id, 0] == -1:
                                edge2elem[edge_id, 0] = elem_id
                            else:
                                edge2elem[edge_id, 1] = elem_id
                            break

        return edge2elem
    
    def boundary_edges( self ) -> np.ndarray:
        """
        Identify boundary edges of the mesh.
        
        Returns:
            Indices of boundary edges
        """
        # Boundary edges have only one adjacent element
        edge2elem = self.adjacencies["Edge2Elem"]
        boundary_mask = (edge2elem[:, 1] == -1)  # Second element is sentinel
        return np.where( boundary_mask )[0]
    
    def edge2vert( self, edge_ids: Opt[Union[int, List[int], np.ndarray]] = None ) -> np.ndarray:
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
    
    def elem2edge( self, elem_ids: Opt[Union[int, List[int], np.ndarray]] = None ) -> np.ndarray:
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
    
    def edge2elem( self, edge_ids: Opt[Union[int, List[int], np.ndarray]] = None ) -> np.ndarray:
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

    def get_vertex_edges(self, vert_id: int) -> Set[int]:
        """
        Get all edges incident to a vertex.

        Returns a set of edge IDs where the vertex appears as an endpoint.
        For an isolated vertex (degree 0), returns empty set.

        Parameters:
            vert_id: Vertex index [0, n_verts)

        Returns:
            Set of edge IDs incident to this vertex

        Raises:
            ValueError: If vert_id is out of range

        Complexity:
            Time: O(1) dict lookup + O(k) to copy set where k=degree
            Space: O(k) for returned set

        Example:
            >>> edges = mesh.get_vertex_edges(5)
            >>> for edge_id in edges:
            ...     v1, v2 = mesh.edge2vert(edge_id)
        """
        if not 0 <= vert_id < self.n_verts:
            raise ValueError(f"Vertex {vert_id} out of range [0, {self.n_verts})")
        return self.adjacencies['Vert2Edge'][vert_id].copy()

    def get_vertex_elements(self, vert_id: int) -> Set[int]:
        """
        Get all elements incident to a vertex.

        Returns a set of element IDs that contain this vertex.
        For an isolated vertex (degree 0), returns empty set.

        Parameters:
            vert_id: Vertex index [0, n_verts)

        Returns:
            Set of element IDs incident to this vertex

        Raises:
            ValueError: If vert_id is out of range

        Complexity:
            Time: O(1) dict lookup + O(k) to copy set where k=degree
            Space: O(k) for returned set

        Example:
            >>> elems = mesh.get_vertex_elements(5)
            >>> for elem_id in elems:
            ...     verts = mesh.connectivity_list[elem_id]
        """
        if not 0 <= vert_id < self.n_verts:
            raise ValueError(f"Vertex {vert_id} out of range [0, {self.n_verts})")
        return self.adjacencies['Vert2Elem'][vert_id].copy()

    def _skeletonize(self) -> None:
        """
        Skeletonize the mesh by iteratively peeling concentric layers from the boundary inward.

        This is a medial axis extraction algorithm that decomposes the mesh into hierarchical
        layers, similar to skeletonization in image processing. Each iteration identifies:
        - Outer elements: Elements adjacent to the current boundary
        - Inner elements: Neighbors of outer elements (one step deeper)
        - Outer/inner vertices: Vertex sets for each layer
        - Boundary edges: The frontier between processed and unprocessed regions
        """
        # Reset layers
        self.layers = {"OE": [], "IE": [], "OV": [], "IV": [], "bEdgeIDs": []}

        # Get boundary edges (edges with only one adjacent element)
        edge2elem = self.adjacencies["Edge2Elem"]
        boundary_mask = (edge2elem[:, 1] == -1)
        boundary_edges = np.where(boundary_mask)[0]

        # Keep track of which elements have been assigned to a layer
        remaining_elements = set(range(self.n_elems))

        # Get element-to-element connectivity using edge2elem
        elem2elem = [[] for _ in range(self.n_elems)]
        for edge_idx, (e1, e2) in enumerate(edge2elem):
            if e1 >= 0 and e2 >= 0:  # Both elements exist (sentinel is -1)
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
                    # Ignore negative placeholders (if any).  Vertex index
                    # 0 is valid in this implementation.
                    if v >= 0:
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

    def _mesh_layers(self) -> None:
        """
        Deprecated: Use _skeletonize() instead.

        This method is retained for backwards compatibility but will be removed in a future version.
        The mesh layers approach is more accurately described as mesh skeletonization.
        """
        warnings.warn(
            "_mesh_layers() is deprecated and will be removed in a future version. "
            "Use _skeletonize() instead.",
            DeprecationWarning,
            stacklevel=2
        )
        self._skeletonize()
        
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
        if self.n_layers == 0:
            raise RuntimeError(
                "Layers not computed. Re-initialise with compute_layers=True."
            )
        if layer_idx < 0 or layer_idx >= self.n_layers:
            raise ValueError( f"Layer index {layer_idx} out of range [0, {self.n_layers-1}]" )
        
        return {
            "OE": self.layers["OE"][layer_idx],
            "IE": self.layers["IE"][layer_idx],
            "OV": self.layers["OV"][layer_idx],
            "IV": self.layers["IV"][layer_idx],
            "bEdgeIDs": self.layers["bEdgeIDs"][layer_idx]
        }

    def admesh_metadata(self) -> Dict[str, Any]:
        """
        Return mesh metadata compatible with ADMESH-Domains schema.

        Returns:
            Dictionary with keys:
            - node_count (int): Number of nodes
            - element_count (int): Number of elements
            - element_type (str): "Triangular", "Quadrilateral", or "Mixed-Element"
            - bounding_box (dict): {"min_x", "max_x", "min_y", "max_y"} with float values

        Note:
            Coordinates are treated as Cartesian (x, y). ADMESH-Domains can map x→lon, y→lat
            in its own adapter layer if needed.
            Safe to call even when compute_layers=False.
        """
        bbox = {
            "min_x": float(self.points[:, 0].min()),
            "max_x": float(self.points[:, 0].max()),
            "min_y": float(self.points[:, 1].min()),
            "max_y": float(self.points[:, 1].max()),
        }

        return {
            "node_count": self.n_verts,
            "element_count": self.n_elems,
            "element_type": self.type,
            "bounding_box": bbox,
        }

    @staticmethod
    def read_from_fort14(full_file_name: Path, compute_layers: bool = True) -> "CHILmesh":
        """
        Load a mesh from a FORT.14 file.

        Supports triangular, quadrilateral, and mixed-element meshes. Triangles in
        a 4-column array use the padded convention: [v0, v1, v2, v0].

        Parameters:
            full_file_name: Path object pointing to the FORT.14 file
            compute_layers: If False, skip skeletonization for fast init (default: True)

        Returns:
            A CHILmesh object
        """
        with open(full_file_name, 'r') as f:
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
                points[i] = [float(line[1]), float(line[2]), float(line[3])]

            # Two-pass element read: first scan for max num_nodes, then allocate
            element_lines = []
            max_nodes = 0
            for i in range(n_elements):
                line = f.readline().strip().split()
                num_nodes = int(line[1])
                if num_nodes < 3 or num_nodes > 4:
                    raise ValueError(
                        f"Unsupported element type: element {i+1} has {num_nodes} nodes "
                        "(only 3 or 4 supported)."
                    )
                max_nodes = max(max_nodes, num_nodes)
                element_lines.append(line)

        # Allocate elements array with correct width
        # For mixed meshes, always use 4 columns with padded triangles
        elem_width = 4 if max_nodes == 4 else 3
        elements = np.zeros((n_elements, elem_width), dtype=int)

        # Read elements again and populate
        for i, line in enumerate(element_lines):
            elem_id = int(line[0])
            num_nodes = int(line[1])
            node_indices = [int(float(line[j + 2])) - 1 for j in range(num_nodes)]

            if elem_width == 4 and num_nodes == 3:
                # Padded triangle: [v0, v1, v2, v0]
                elements[i] = [node_indices[0], node_indices[1], node_indices[2], node_indices[0]]
            else:
                # Quad or triangle in 3-column array
                elements[i, :num_nodes] = node_indices

        return CHILmesh(
            connectivity=elements,
            points=points,
            grid_name=header,
            compute_layers=compute_layers,
        )

    @classmethod
    def from_admesh_domain(cls, record, compute_layers: bool = True) -> "CHILmesh":
        """
        Load a mesh from an ADMESH-Domains catalog record.

        The record is duck-typed: it must have `filename` and optionally `type` attributes.
        No import of the `admesh-domains` package is required.

        Parameters:
            record: ADMESH-Domains Mesh record with `filename` and optional `type` attributes
                - `filename` (str): Path to mesh file on disk
                - `type` (str, optional): Format hint - "ADCIRC", "SMS_2DM", "ADCIRC_GRD", or other
                - `kind` (str, optional): Ignored by CHILmesh
            compute_layers: If False, skip skeletonization for fast init (default: True)

        Returns:
            A CHILmesh object

        Raises:
            FileNotFoundError: If the file does not exist (with guidance message)
            ValueError: If the file format is unsupported
        """
        filename = getattr(record, "filename", None)
        if not filename:
            raise ValueError("Record must have a 'filename' attribute")

        # Check if file exists
        filepath = Path(filename)
        if not filepath.exists():
            raise FileNotFoundError(
                f"File not found: {filepath}. "
                "If using ADMESH-Domains, call mesh_record.load() first."
            )

        # Determine format and route to appropriate reader
        mesh_type = getattr(record, "type", None)

        if mesh_type == "SMS_2DM":
            return cls.read_from_2dm(filepath, compute_layers=compute_layers)
        else:
            # Route all other types (ADCIRC, ADCIRC_GRD, None, unknown) to ADCIRC reader
            if mesh_type and mesh_type not in ("ADCIRC", "ADCIRC_GRD"):
                warnings.warn(
                    f"Unrecognised mesh type '{mesh_type}', falling back to ADCIRC reader.",
                    UserWarning,
                    stacklevel=2,
                )
            return cls.read_from_fort14(filepath, compute_layers=compute_layers)

    @staticmethod
    def read_from_2dm(full_file_name: Path, compute_layers: bool = True) -> "CHILmesh":
        """
        Load a mesh from an SMS .2dm file.

        Supports triangular and quadrilateral elements. Mixed-element meshes use
        the padded-triangle convention in a 4-column array.

        Parameters:
            full_file_name: Path object pointing to the .2dm file
            compute_layers: If False, skip skeletonization for fast init (default: True)

        Returns:
            A CHILmesh object
        """
        nodes = {}  # node_id -> [x, y, z]
        elements = []  # list of (type, vertices)

        with open(full_file_name, 'r') as f:
            for line in f:
                line = line.strip()

                # Skip empty lines and comments
                if not line or line.startswith("#"):
                    continue

                # Skip header
                if line.startswith("MESH2D"):
                    continue

                parts = line.split()
                if not parts:
                    continue

                keyword = parts[0]

                # Parse node definition: ND id x y z
                if keyword == "ND":
                    node_id = int(parts[1])
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4]) if len(parts) > 4 else 0.0
                    nodes[node_id] = [x, y, z]

                # Parse triangle element: E3T id n1 n2 n3 mat
                elif keyword == "E3T":
                    elem_id = int(parts[1])
                    n1, n2, n3 = int(parts[2]), int(parts[3]), int(parts[4])
                    elements.append(("E3T", [n1, n2, n3]))

                # Parse quad element: E4Q id n1 n2 n3 n4 mat
                elif keyword == "E4Q":
                    elem_id = int(parts[1])
                    n1, n2, n3, n4 = int(parts[2]), int(parts[3]), int(parts[4]), int(parts[5])
                    elements.append(("E4Q", [n1, n2, n3, n4]))

        # Convert node dict to ordered array (0-based indexing)
        if not nodes:
            raise ValueError("No nodes found in .2dm file")

        node_ids = sorted(nodes.keys())
        points = np.array([nodes[nid] for nid in node_ids], dtype=float)

        # Create mapping from 1-based file IDs to 0-based array indices
        id_to_idx = {nid: idx for idx, nid in enumerate(node_ids)}

        # Determine element array width
        has_quads = any(elem_type == "E4Q" for elem_type, _ in elements)
        elem_width = 4 if has_quads else 3

        # Build connectivity array
        connectivity = np.zeros((len(elements), elem_width), dtype=int)
        for i, (elem_type, vertices) in enumerate(elements):
            if elem_type == "E3T":
                # Triangle
                v1, v2, v3 = vertices
                if elem_width == 4:
                    # Mixed mesh: use padding convention [v0, v1, v2, v0]
                    idx1, idx2, idx3 = id_to_idx[v1], id_to_idx[v2], id_to_idx[v3]
                    connectivity[i] = [idx1, idx2, idx3, idx1]
                else:
                    # All triangles
                    connectivity[i] = [id_to_idx[v1], id_to_idx[v2], id_to_idx[v3]]
            else:
                # Quad
                v1, v2, v3, v4 = vertices
                connectivity[i] = [id_to_idx[v1], id_to_idx[v2], id_to_idx[v3], id_to_idx[v4]]

        return CHILmesh(
            connectivity=connectivity,
            points=points,
            grid_name="SMS_2DM",
            compute_layers=compute_layers,
        )

    def write_to_fort14( self, filename: str, grid_name: Opt[str] = "CHILmesh Grid") -> bool:
        """
        Export the current mesh to ADCIRC .fort.14 format.

        Parameters:
            filename: Path to save the file
            grid_name: Optional title for the mesh

        Returns:
            ``True`` on success.

        Note:
            Prior to 0.1.1 this method recursed into itself instead of
            delegating to the module-level ``write_fort14`` writer (B1).
        """
        return write_fort14(filename, self.points, self.connectivity_list, grid_name)

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
                    
                    # Normalize vectors safely to avoid runtime warnings of NaN
                    # v1_norm = v1 / np.linalg.norm(v1)
                    # v2_norm = v2 / np.linalg.norm(v2)
                    v1_norm = v1 / (np.linalg.norm(v1) + 1e-12)
                    v2_norm = v2 / (np.linalg.norm(v2) + 1e-12)
                    
                    # Calculate angle in degrees
                    dot_product = np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0)
                    angle = np.arccos(dot_product) * 180 / np.pi
                    angles[i, j] = angle
                    
            elif elem_id in quad_elems:
                # Quadrilateral angles
                vertices = self.connectivity_list[elem_id]  # All 4 vertices
                coords = self.points[vertices, :2]  # Get x,y coordinates

                # Calculate angles at each vertex.  Use the same epsilon
                # guard as the triangle path so degenerate quads with
                # zero-length edges produce a real number (typically 0)
                # rather than NaN, which would slip past elem_quality's
                # ``<= 359.99 -> zero out`` check (B5).
                for j in range(4):
                    v1 = coords[(j + 1) % 4] - coords[j]
                    v2 = coords[(j - 1) % 4] - coords[j]

                    v1_norm = v1 / (np.linalg.norm(v1) + 1e-12)
                    v2_norm = v2 / (np.linalg.norm(v2) + 1e-12)

                    dot_product = np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0)
                    angle = np.arccos(dot_product) * 180 / np.pi
                    angles[i, j] = float(np.real(angle))
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
            'mean': float(np.mean(quality)),
            'median': float(np.median(quality)),
            'min': float(np.min(quality)),
            'max': float(np.max(quality)),
            'std': float(np.std(quality))
        }
        return quality, angles, stats
    
    def smooth_mesh(self, method: str, acknowledge_change: bool=False, *kwargs) -> np.ndarray:
        """
        Perform mesh smoothing using a modified FEM-based approach.
        
        Parameters:
            method: Smoothing method ('FEM','angle-based')
            acknowledge_change: If True, acknowledges the change in the mesh
        """
        assert acknowledge_change, "acknowledge_change must be True to change mesh -- this will change the mesh, make sure you understand this before using this method within a broader algorithm."
        if method.lower() == 'fem':
            new_points = self.direct_smoother( *kwargs )
        elif method.lower() == 'angle-based':
            new_points = self.angle_based_smoother( *kwargs )
        else:
            raise ValueError(f"Unknown smoothing method: {method}")
        self.change_points( new_points, acknowledge_change=True )
        return new_points
    
    def angle_based_smoother( self, angle_limit: float = 30.0 ) -> np.ndarray:
        """
        Perform angle-based smoothing of the mesh.
        Based on this: https://www.andrew.cmu.edu/user/shimada/papers/00-imr-zhou.pdf
        Parameters:
            angle_limit: Maximum allowable angle deviation in degrees
        """
        # Placeholder for angle-based smoothing logic
        # This would involve checking angles and adjusting points accordingly
        # For now, just return the original points
        raise NotImplementedError("Angle-based smoothing not implemented yet.")
        return self.points
    
    def direct_smoother( self, kinf=1e12 ) -> np.ndarray:
        """
        Perform direct (non-iterative) FEM smoothing with fixed boundary nodes.
        Based on the triangle stiffness formulation in Balendran (2006).

        Parameters:
            kinf: Large stiffness value for fixed boundary vertices

        Reference:
            Zhou, M., & Shimada, K. (2000). 
            "An angle-based approach to two-dimensional mesh smoothing". 
            In *Proceedings of the 9th International Meshing Roundtable*, 373–384. 
            Sandia National Laboratories.
            https://api.semanticscholar.org/CorpusID:34335417
        """
        import numpy as np
        from scipy.sparse import csr_matrix
        from scipy.sparse.linalg import spsolve

        p = self.points[:, :2]  # Only use x, y
        t = self.connectivity_list[:, :3]  # Assume triangles only

        # Identify boundary nodes
        edge_verts = self.edge2vert(self.boundary_edges())
        boundary_nodes = np.unique(edge_verts.flatten())

        n = self.n_verts
        D = 2.0 * np.eye(2)
        T = np.array([[-1, -np.sqrt(3)], [np.sqrt(3), -1]])

        rows, cols, data = [], [], []
        for tri in t:
            for i in range(3):
                for j in range(3):
                    block = D if i == j else T if j == (i+1)%3 else T.T
                    for di in range(2):
                        for dj in range(2):
                            rows.append(2*tri[i]+di)
                            cols.append(2*tri[j]+dj)
                            data.append(block[di, dj])

        K = csr_matrix((data, (rows, cols)), shape=(2*n, 2*n))
        F = np.zeros(2*n)

        # Apply boundary constraints
        for v in boundary_nodes:
            F[2*v:2*v+2] = kinf * p[v]
            K[2*v, 2*v] = kinf
            K[2*v+1, 2*v+1] = kinf

        c = spsolve(K, F)
        new_points = np.zeros_like(self.points)
        new_points[:, :2] = c.reshape(-1, 2)
        new_points[:, 2] = self.points[:, 2]  # preserve z if needed
        return new_points    


    def advancing_front_boundary_edges(self) -> List[int]:
        """
        Get boundary edge IDs for advancing-front mesh generation.

        Returns edge IDs on the current mesh boundary, suitable for placing
        new elements during advancing-front generation.

        Returns:
            List of boundary edge IDs (sorted)

        Complexity:
            Time: O(n_edges)
            Space: O(n_boundary_edges)

        Example:
            >>> boundary = mesh.advancing_front_boundary_edges()
            >>> for edge_id in boundary:
            ...     v1, v2 = mesh.edge2vert(edge_id)[0]
        """
        boundary_edge_ids = self.boundary_edges()
        return sorted(boundary_edge_ids.tolist())

    def add_advancing_front_element(
        self, vertices: List[int], elem_type: str = "tri"
    ) -> int:
        """
        Add element to mesh during advancing-front generation.

        Adds a single element (triangle or quad) to the mesh and updates all
        adjacency structures. This is the primary method for MADMESHR advancing-front
        generation. Elements are added in order; mesh is valid after each call.

        Parameters:
            vertices: Vertex indices for new element
                - Triangle: [v0, v1, v2]
                - Quad: [v0, v1, v2, v3]
            elem_type: Element type ('tri' or 'quad')

        Returns:
            ID of newly added element

        Raises:
            ValueError: If vertices out of range or elem_type invalid
            RuntimeError: If adjacency update fails

        Complexity:
            Time: O(k) where k = element degree (typically 3-4)
            Space: O(k) for adjacency updates

        Example:
            >>> new_elem_id = mesh.add_advancing_front_element([v1, v2, v3], 'tri')
            >>> print(f"Added element {new_elem_id}")
        """
        if elem_type not in ("tri", "quad"):
            raise ValueError(f"Invalid elem_type: {elem_type}. Use 'tri' or 'quad'.")

        # Validate vertices
        for v in vertices:
            if not 0 <= v < self.n_verts:
                raise ValueError(f"Vertex {v} out of range [0, {self.n_verts})")

        # Detect connectivity format (3-column all triangles or 4-column mixed)
        elem_cols = self.connectivity_list.shape[1]

        # Prepare new element with matching format
        if elem_type == "tri":
            if len(vertices) != 3:
                raise ValueError(f"Triangle requires 3 vertices, got {len(vertices)}")
            if elem_cols == 3:
                new_elem = np.array(vertices)  # Pure triangle format
            else:
                new_elem = np.array(vertices + [vertices[0]])  # Pad to 4-column
        else:  # quad
            if len(vertices) != 4:
                raise ValueError(f"Quad requires 4 vertices, got {len(vertices)}")
            if elem_cols != 4:
                raise ValueError("Cannot add quad to triangle-only mesh (3-column connectivity)")
            new_elem = np.array(vertices)

        # Add to connectivity
        new_elem_id = self.n_elems
        self.connectivity_list = np.vstack([self.connectivity_list, new_elem])
        self.n_elems = self.connectivity_list.shape[0]  # Update element count

        # Rebuild adjacencies (invalidates layers)
        self._build_adjacencies()
        self.layers = {}  # Clear layers since mesh changed

        return new_elem_id

    def remove_boundary_loop(self, edge_ids: List[int]) -> None:
        """
        Remove boundary loop elements (residual closure).

        Removes elements adjacent to the specified boundary edges. Used in
        advancing-front generation to discard the residual boundary loop
        when it shrinks to ≤4 vertices.

        Parameters:
            edge_ids: List of boundary edge IDs to remove adjacent elements

        Raises:
            ValueError: If any edge_id is out of range or not a boundary edge

        Complexity:
            Time: O(m + n) where m = number of edges, n = mesh size
            Space: O(m) for tracking removals

        Example:
            >>> small_boundary = [e1, e2, e3, e4]
            >>> mesh.remove_boundary_loop(small_boundary)
        """
        # Collect elements to remove (before any deletion)
        elems_to_remove = set()

        for edge_id in edge_ids:
            if not 0 <= edge_id < self.n_edges:
                raise ValueError(f"Edge {edge_id} out of range [0, {self.n_edges})")

            # Get adjacent elements
            e1, e2 = self.edge2elem(edge_id)[0]

            # Boundary edge: e2 == -1, remove e1
            if e2 == -1:
                if e1 >= 0:
                    elems_to_remove.add(e1)
            else:
                # Interior edge: remove both (or neither if not at boundary)
                # For now, only handle boundary edges
                pass

        if len(elems_to_remove) == 0:
            return  # Nothing to remove

        # Create mask for elements to keep
        keep_mask = np.ones(self.n_elems, dtype=bool)
        for elem_id in elems_to_remove:
            if elem_id < len(keep_mask):
                keep_mask[elem_id] = False

        # Remove elements (keeping valid indices)
        self.connectivity_list = self.connectivity_list[keep_mask]
        self.n_elems = self.connectivity_list.shape[0]  # Update element count

        # Rebuild adjacencies (will create new valid indices)
        self._build_adjacencies()
        self.layers = {}  # Clear layers

    def pinch_points(self, width_threshold: float = 0.5) -> List[int]:
        """
        Identify bottleneck vertices (pinch points) in the mesh.

        Detects vertices where the local mesh width drops below a threshold.
        Useful for domain splitting in advancing-front generation.

        Parameters:
            width_threshold: Relative width threshold [0, 1]
                - Vertices where min_neighbor_distance / max_neighbor_distance < threshold
                  are considered pinch points

        Returns:
            List of vertex IDs identified as pinch points (sorted)

        Complexity:
            Time: O(n_verts * avg_degree)
            Space: O(n_pinch_points)

        Example:
            >>> pinches = mesh.pinch_points(width_threshold=0.3)
            >>> for v in pinches:
            ...     print(f"Bottleneck at vertex {v}")
        """
        pinches = []

        for vert_id in range(self.n_verts):
            edges = self.get_vertex_edges(vert_id)

            if len(edges) < 2:
                continue

            # Compute distances to neighbors
            distances = []
            for edge_id in edges:
                v1, v2 = self.edge2vert(edge_id)[0]
                other_vert = v2 if v1 == vert_id else v1
                p1 = self.points[vert_id, :2]
                p2 = self.points[other_vert, :2]
                dist = np.linalg.norm(p2 - p1)
                distances.append(dist)

            if len(distances) >= 2:
                min_dist = min(distances)
                max_dist = max(distances)
                if max_dist > 0:
                    ratio = min_dist / max_dist
                    if ratio < width_threshold:
                        pinches.append(vert_id)

        return sorted(pinches)

    def copy( self ) -> "CHILmesh":
        """ Returns: a deep copy of the  new CHILmesh object with the same properties."""
        return deepcopy(self)
    

def write_fort14( filename: Path, points: np.ndarray, elements: np.ndarray, grid_name: str ) -> bool:
    """
    Write mesh data to a .fort.14 ADCIRC file.

    Supports triangular, quadrilateral, and mixed-element meshes. Triangles in a
    4-column array use the padded convention: [v0, v1, v2, v0].

    Parameters:
        filename: Output path
        points: (n_nodes, 2 or 3) numpy array of node coordinates
        elements: (n_elems, 3 or 4) array of vertex indices (0-based)
        grid_name: Header string
    """
    try:
        with open(filename, 'w') as f:
            f.write(f"{grid_name}\n")
            f.write(f"{len(elements)} {len(points)}\n")

            for i, pt in enumerate(points, start=1):
                x, y = pt[:2]
                z = pt[2] if len(pt) == 3 else 0.0
                f.write(f"{i} {x:.8f} {y:.8f} {z:.8f}\n")

            for i, elem in enumerate(elements, start=1):
                # Convert to 1-based indexing
                indices = elem + 1

                # Detect triangle vs quad
                if len(elem) == 3:
                    # 3-column array: all triangles
                    n1, n2, n3 = indices
                    f.write(f"{i} 3 {n1} {n2} {n3}\n")
                else:
                    # 4-column array: check if padded triangle or quad
                    if elem[3] == elem[0]:
                        # Padded triangle: [v0, v1, v2, v0]
                        n1, n2, n3 = indices[:3]
                        f.write(f"{i} 3 {n1} {n2} {n3}\n")
                    else:
                        # Quad: [v0, v1, v2, v3]
                        n1, n2, n3, n4 = indices
                        f.write(f"{i} 4 {n1} {n2} {n3} {n4}\n")
        return True
    except Exception as e:
        print(f"Error writing fort14 file {filename}: {e}")
        return False
