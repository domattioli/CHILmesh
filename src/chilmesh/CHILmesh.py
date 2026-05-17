from __future__ import annotations

from pathlib import Path
import warnings

from .utils.plot_utils import CHILmeshPlotMixin

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.spatial import Delaunay, cKDTree
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

    Layer Structure (Skeletonization):
        After initialization, the mesh is decomposed into concentric layers via skeletonization.
        Access layers via the ``layers`` property (or ``Layers`` for backwards compatibility):

        - mesh.layers['OE'][i]: Array of **element indices** in the outer boundary of layer i
        - mesh.layers['IE'][i]: Array of **element indices** one step deeper than layer i's boundary
        - mesh.layers['OV'][i]: Array of **vertex indices** on the outer boundary of layer i
        - mesh.layers['IV'][i]: Array of **vertex indices** interior to layer i
        - mesh.layers['bEdgeIDs'][i]: Array of **edge indices** that form the boundary between layer i and i+1

        Layer 0 is the outermost boundary; layers progress inward toward the mesh interior.
        Use mesh.elements_in_layer(i) for a convenience method that returns all element indices
        (union of OE[i] and IE[i]) for a layer.

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

    def elements_in_layer(self, layer_idx: int) -> np.ndarray:
        """
        Get all element indices for a given skeletonization layer.

        This is a convenience method that returns the union of outer and inner
        elements for the specified layer.

        Parameters:
            layer_idx: Zero-based layer index (0 = outermost)

        Returns:
            Array of element indices (sorted) belonging to this layer

        Raises:
            IndexError: If layer_idx is outside [0, n_layers)

        Example:
            >>> mesh = chilmesh.examples.annulus()
            >>> elements = mesh.elements_in_layer(0)  # All elements in outermost layer
        """
        if layer_idx < 0 or layer_idx >= self.n_layers:
            raise IndexError(f"Layer index {layer_idx} out of range [0, {self.n_layers})")
        outer = self.layers["OE"][layer_idx]
        inner = self.layers["IE"][layer_idx]
        return np.sort(np.concatenate([outer, inner]).astype(int))

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
            connectivity: Element connectivity list (n_elems × 3 for triangles, n_elems × 4 for quads/mixed)
            points: Vertex coordinates (n_verts × 3, with z=0 for 2D meshes)
            grid_name: Name of the mesh
            compute_layers: If False, skip skeletonization for fast init (default: True)

        Attributes set during initialization:
            layers (dict): Skeletonization result with keys:
                - 'OE': List of arrays (outer elements per layer, **element indices**)
                - 'IE': List of arrays (inner elements per layer, **element indices**)
                - 'OV': List of arrays (outer vertices per layer, **vertex indices**)
                - 'IV': List of arrays (inner vertices per layer, **vertex indices**)
                - 'bEdgeIDs': List of arrays (boundary edges per layer, **edge indices**)
            n_layers (int): Number of skeletonization layers
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

            # Build spatial indices (always, regardless of compute_layers)
            self._build_spatial_indices()
    
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
            elem_ids = np.array( [elem_ids] )
        elem_ids = np.asarray(elem_ids)

        if self.connectivity_list.shape[1] == 3:
            # Pure triangular mesh — single vectorised path
            verts = self.connectivity_list[elem_ids]          # (n, 3)
            xy = self.points[verts, :2]                       # (n, 3, 2)
            x, y = xy[:, :, 0], xy[:, :, 1]
            return 0.5 * (
                x[:, 0] * (y[:, 1] - y[:, 2])
                + x[:, 1] * (y[:, 2] - y[:, 0])
                + x[:, 2] * (y[:, 0] - y[:, 1])
            )

        # 4-column connectivity: separate padded triangles from quads
        rows = self.connectivity_list[elem_ids]               # (n, 4)
        tri_mask = (
            (rows[:, 0] == rows[:, 1])
            | (rows[:, 1] == rows[:, 2])
            | (rows[:, 2] == rows[:, 3])
            | (rows[:, 3] == rows[:, 0])
            | (rows[:, 0] == rows[:, 2])
            | (rows[:, 1] == rows[:, 3])
        )
        areas = np.zeros(len(elem_ids))

        if tri_mask.any():
            verts = rows[tri_mask, :3]
            xy = self.points[verts, :2]
            x, y = xy[:, :, 0], xy[:, :, 1]
            areas[tri_mask] = 0.5 * (
                x[:, 0] * (y[:, 1] - y[:, 2])
                + x[:, 1] * (y[:, 2] - y[:, 0])
                + x[:, 2] * (y[:, 0] - y[:, 1])
            )

        quad_mask = ~tri_mask
        if quad_mask.any():
            verts = rows[quad_mask]
            xy = self.points[verts, :2]
            x, y = xy[:, :, 0], xy[:, :, 1]
            areas[quad_mask] = 0.5 * (
                x[:, 0] * (y[:, 1] - y[:, 3])
                + x[:, 1] * (y[:, 2] - y[:, 0])
                + x[:, 2] * (y[:, 3] - y[:, 1])
                + x[:, 3] * (y[:, 0] - y[:, 2])
            )

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

                # Skip invalid edges (negative vertex ids or padding edges v==v in padded triangles)
                # In MATLAB the value 0 signified a placeholder for a missing
                # vertex in mixed element meshes.  In this Python port we use
                # 0-based indexing, therefore vertex id ``0`` is valid and
                # should not be discarded.  Only negative ids are considered
                # invalid.
                if v1 < 0 or v2 < 0 or v1 == v2:
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

                # Skip invalid edges (negative vertex ids or padding edges v==v in padded triangles)
                if v1 < 0 or v2 < 0 or v1 == v2:
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
        if "Edge2Elem" in self.adjacencies:
            # Fast path: adjacency structure already built
            edge2elem = self.adjacencies["Edge2Elem"]
            boundary_mask = (edge2elem[:, 1] == -1)
            return np.where(boundary_mask)[0]

        # Fallback for meshes constructed with compute_layers=False:
        # count edge occurrences; boundary edges appear exactly once.
        edge_count: dict[tuple, int] = {}
        edge_list: list[tuple] = []
        for row in self.connectivity_list:
            verts = list(row[:3]) if (len(row) == 4 and row[3] == row[0]) else list(row)
            for i in range(len(verts)):
                a, b = int(verts[i]), int(verts[(i + 1) % len(verts)])
                key = (min(a, b), max(a, b))
                if key not in edge_count:
                    edge_list.append(key)
                edge_count[key] = edge_count.get(key, 0) + 1
        return np.array([i for i, key in enumerate(edge_list) if edge_count[key] == 1], dtype=int)

    def boundary_node_indices( self ) -> np.ndarray:
        """
        Return the indices of all boundary vertices.

        Returns:
            Sorted array of vertex indices that lie on at least one boundary edge.
        """
        edge_verts = self.edge2vert(self.boundary_edges())
        return np.unique(edge_verts.flatten())
    
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

    def _build_spatial_indices(self) -> None:
        """Build KD-trees for fast spatial queries on vertices and element centroids."""
        self._vertex_tree = cKDTree(self.points[:, :2])
        self._centroid_tree = cKDTree(self._get_centroids())

    def _get_centroids(self) -> np.ndarray:
        """Compute element centroids (mean of vertex coordinates)."""
        n_cols = self.connectivity_list.shape[1]
        centroids = np.zeros((self.n_elems, 2))
        for i, elem in enumerate(self.connectivity_list):
            verts = self.points[elem[:n_cols], :2]
            centroids[i] = np.mean(verts, axis=0)
        return centroids

    @property
    def centroids(self) -> np.ndarray:
        """Get element centroids (lazy computation)."""
        if not hasattr(self, '_centroids_cache'):
            self._centroids_cache = self._get_centroids()
        return self._centroids_cache

    def _point_in_triangle(self, point: np.ndarray, v0: np.ndarray, v1: np.ndarray, v2: np.ndarray) -> bool:
        """Check if point is inside triangle using barycentric coordinates."""
        denom = (v1[1] - v2[1]) * (v0[0] - v2[0]) + (v2[0] - v1[0]) * (v0[1] - v2[1])
        if abs(denom) < 1e-10:
            return False
        a = ((v1[1] - v2[1]) * (point[0] - v2[0]) + (v2[0] - v1[0]) * (point[1] - v2[1])) / denom
        b = ((v2[1] - v0[1]) * (point[0] - v2[0]) + (v0[0] - v2[0]) * (point[1] - v2[1])) / denom
        c = 1 - a - b
        return a >= -1e-10 and b >= -1e-10 and c >= -1e-10

    def _point_in_quad(self, point: np.ndarray, v0: np.ndarray, v1: np.ndarray, v2: np.ndarray, v3: np.ndarray) -> bool:
        """Check if point is inside quad by checking two triangles."""
        return (self._point_in_triangle(point, v0, v1, v2) or
                self._point_in_triangle(point, v0, v2, v3))

    def _point_in_element(self, point: np.ndarray, elem_id: int) -> bool:
        """Check if point is inside element (handles tri and quad)."""
        elem = self.connectivity_list[elem_id]
        verts = self.points[elem, :2]
        if elem.size == 3 or elem[3] == elem[2]:
            return self._point_in_triangle(point, verts[0], verts[1], verts[2])
        else:
            return self._point_in_quad(point, verts[0], verts[1], verts[2], verts[3])

    def find_element(self, point: np.ndarray) -> int:
        """
        Find element containing point using KD-tree + barycentric check.

        Parameters:
            point: 2D coordinate [x, y]

        Returns:
            Element ID if found, -1 if point is outside mesh

        Example:
            >>> elem_id = mesh.find_element(np.array([0.5, 0.5]))
        """
        point = np.asarray(point)
        _, candidates = self._centroid_tree.query(point, k=min(20, self.n_elems))
        if isinstance(candidates, np.int64):
            candidates = np.array([candidates])
        for elem_id in candidates:
            # Skip deleted elements (marked as all zeros)
            elem = self.connectivity_list[elem_id]
            if np.all(elem == 0):
                continue
            if self._point_in_element(point, elem_id):
                return int(elem_id)
        return -1

    def find_elements_in_radius(self, point: np.ndarray, radius: float) -> np.ndarray:
        """
        Find all elements within radius of point.

        Parameters:
            point: 2D coordinate [x, y]
            radius: Search radius (must be >= 0)

        Returns:
            Array of element IDs within radius

        Example:
            >>> elems = mesh.find_elements_in_radius(np.array([0.5, 0.5]), radius=0.1)
        """
        if radius < 0:
            raise ValueError(f"radius must be >= 0, got {radius}")
        point = np.asarray(point)
        elem_ids = self._centroid_tree.query_ball_point(point, radius)
        return np.array(elem_ids, dtype=int)

    def nearest_vertices(self, point: np.ndarray, k: int = 1) -> np.ndarray:
        """
        Find k nearest vertices to point.

        Parameters:
            point: 2D coordinate [x, y]
            k: Number of nearest vertices to return

        Returns:
            Array of k nearest vertex IDs

        Example:
            >>> verts = mesh.nearest_vertices(np.array([0.5, 0.5]), k=3)
        """
        point = np.asarray(point)
        k = min(k, self.n_verts)
        if k <= 0:
            return np.array([], dtype=int)
        _, vert_ids = self._vertex_tree.query(point, k=k)
        if k == 1:
            return np.array([vert_ids], dtype=int)
        return np.array(vert_ids, dtype=int)

    def nearest_elements(self, point: np.ndarray, k: int = 1) -> np.ndarray:
        """
        Find k nearest elements (by centroid) to point.

        Parameters:
            point: 2D coordinate [x, y]
            k: Number of nearest elements to return

        Returns:
            Array of k nearest element IDs

        Example:
            >>> elems = mesh.nearest_elements(np.array([0.5, 0.5]), k=5)
        """
        point = np.asarray(point)
        k = min(k, self.n_elems)
        if k <= 0:
            return np.array([], dtype=int)
        _, elem_ids = self._centroid_tree.query(point, k=k)
        if k == 1:
            return np.array([elem_ids], dtype=int)
        return np.array(elem_ids, dtype=int)

    def _skeletonize(self) -> None:
        """
        Skeletonize the mesh by iteratively peeling concentric layers inward.

        Faithful Python port of the MATLAB ``meshLayers`` function.  Each
        iteration of layer ``iL`` produces:

        - ``OV[iL]``: outer ring of vertices bounding the layer's outer side.
        - ``OE[iL]``: elements adjacent to the boundary edges of ``OV[iL]``.
        - ``IE[iL]``: active elements adjacent to ANY edge whose
          ``Edge2Vert`` entry references a vertex in ``OV[iL]``.
        - ``IV[iL]``: vertices of (``OE`` ∪ ``IE``) not in ``OV[iL]``.
        - ``bEdgeIDs[iL]``: boundary edges defining this layer's outer frontier.

        The vertex-adjacency definition of IE (broader than element-adjacency)
        guarantees the layer separation invariant: vertices from layer k cannot
        appear in any element of layer k+2 or beyond.

        Indexing: layers are 0-indexed; ``-1`` is the sentinel for consumed
        vertices/elements in the working copies of Edge2Vert/Edge2Elem.
        """
        self.layers = {"OE": [], "IE": [], "OV": [], "IV": [], "bEdgeIDs": []}

        edge2vert_work = self.adjacencies["Edge2Vert"].copy()
        edge2elem_work = self.adjacencies["Edge2Elem"].copy()

        iL = 0
        while np.any(edge2elem_work >= 0):
            # Step 1: boundary edges
            if iL == 0:
                iLbEdgeIDs = np.where(self.adjacencies["Edge2Elem"][:, 1] == -1)[0]
            else:
                active_count = np.sum(edge2elem_work >= 0, axis=1)
                iLbEdgeIDs = np.where(active_count == 1)[0]

            if len(iLbEdgeIDs) == 0:
                break

            # Step 2: OV = unique vertices on boundary edges (from working copy)
            ov_raw = edge2vert_work[iLbEdgeIDs].ravel()
            ov = np.unique(ov_raw[ov_raw >= 0]).astype(int)
            self.layers["OV"].append(ov)
            self.layers["bEdgeIDs"].append(iLbEdgeIDs)

            # Step 3: OE = active elements adjacent to those boundary edges
            oe_raw = edge2elem_work[iLbEdgeIDs].ravel()
            oe = np.unique(oe_raw[oe_raw >= 0]).astype(int)
            self.layers["OE"].append(oe)

            # Step 4: consume OE from edge2elem_work
            if len(oe) > 0:
                edge2elem_work[np.isin(edge2elem_work, oe)] = -1

            # Step 5: find edges touching ANY OV vertex (broader than just boundary edges)
            ov_edge_mask = np.any(np.isin(edge2vert_work, ov), axis=1)
            ov_edge_indices = np.where(ov_edge_mask)[0]

            # Step 6: IE = active elements adjacent to those edges
            if len(ov_edge_indices) > 0:
                ie_raw = edge2elem_work[ov_edge_indices].ravel()
                ie = np.unique(ie_raw[ie_raw >= 0]).astype(int)
            else:
                ie = np.empty(0, dtype=int)
            self.layers["IE"].append(ie)

            # Step 7: consume OV vertices and IE elements
            if len(ov) > 0:
                edge2vert_work[np.isin(edge2vert_work, ov)] = -1
            if len(ie) > 0:
                edge2elem_work[np.isin(edge2elem_work, ie)] = -1

            # Step 8: IV = vertices of (OE ∪ IE) connectivity, minus OV
            if len(oe) > 0 or len(ie) > 0:
                layer_elems = np.concatenate((oe, ie))
                lv = self.connectivity_list[layer_elems].ravel()
                lv = lv[lv >= 0]
                iv = np.setdiff1d(np.unique(lv), ov).astype(int)
            else:
                iv = np.empty(0, dtype=int)
            self.layers["IV"].append(iv)

            iL += 1

        self.n_layers = iL

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
            Array of interior angles (degrees) with shape (n_elems, 3) for
            triangular meshes or (n_elems, 4) for quad/mixed meshes.
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        if np.isscalar(elem_ids):
            elem_ids = np.array([elem_ids])
        elem_ids = np.asarray(elem_ids)

        if self.connectivity_list.shape[1] == 3:
            # Pure triangular mesh — vectorised over all elements,
            # inner loop only over 3 vertices (constant).
            verts = self.connectivity_list[elem_ids]   # (n, 3)
            xy = self.points[verts, :2]                # (n, 3, 2)
            angles = np.zeros((len(elem_ids), 3))
            for j in range(3):
                v1 = xy[:, (j + 1) % 3, :] - xy[:, j, :]  # (n, 2)
                v2 = xy[:, (j - 1) % 3, :] - xy[:, j, :]  # (n, 2)
                n1 = np.linalg.norm(v1, axis=1, keepdims=True) + 1e-12
                n2 = np.linalg.norm(v2, axis=1, keepdims=True) + 1e-12
                dot = np.clip(((v1 / n1) * (v2 / n2)).sum(axis=1), -1.0, 1.0)
                angles[:, j] = np.degrees(np.arccos(dot))
            return angles

        # 4-column connectivity
        rows = self.connectivity_list[elem_ids]        # (n, 4)
        tri_mask = (
            (rows[:, 0] == rows[:, 1])
            | (rows[:, 1] == rows[:, 2])
            | (rows[:, 2] == rows[:, 3])
            | (rows[:, 3] == rows[:, 0])
            | (rows[:, 0] == rows[:, 2])
            | (rows[:, 1] == rows[:, 3])
        )
        quad_mask = ~tri_mask
        n_angle_cols = 4 if quad_mask.any() else 3
        angles = np.zeros((len(elem_ids), n_angle_cols))

        if tri_mask.any():
            xy = self.points[rows[tri_mask, :3], :2]   # (n_tri, 3, 2)
            for j in range(3):
                v1 = xy[:, (j + 1) % 3, :] - xy[:, j, :]
                v2 = xy[:, (j - 1) % 3, :] - xy[:, j, :]
                n1 = np.linalg.norm(v1, axis=1, keepdims=True) + 1e-12
                n2 = np.linalg.norm(v2, axis=1, keepdims=True) + 1e-12
                dot = np.clip(((v1 / n1) * (v2 / n2)).sum(axis=1), -1.0, 1.0)
                angles[tri_mask, j] = np.degrees(np.arccos(dot))

        if quad_mask.any():
            xy = self.points[rows[quad_mask], :2]      # (n_quad, 4, 2)
            for j in range(4):
                v1 = xy[:, (j + 1) % 4, :] - xy[:, j, :]
                v2 = xy[:, (j - 1) % 4, :] - xy[:, j, :]
                n1 = np.linalg.norm(v1, axis=1, keepdims=True) + 1e-12
                n2 = np.linalg.norm(v2, axis=1, keepdims=True) + 1e-12
                dot = np.clip(((v1 / n1) * (v2 / n2)).sum(axis=1), -1.0, 1.0)
                angles[quad_mask, j] = np.degrees(np.arccos(dot))

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
            Tuple of (Quality, Angles, stats) where:
            - Quality: Array of quality measurements for each element
            - Angles: Array of interior angles for each element
            - stats: Dict with mean, median, min, max, std
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        if np.isscalar(elem_ids):
            elem_ids = np.array([elem_ids])
        elem_ids = np.asarray(elem_ids)

        # Compute boolean masks without O(n²) membership checks
        if self.connectivity_list.shape[1] == 3:
            tri_mask = np.ones(len(elem_ids), dtype=bool)
            quad_mask = np.zeros(len(elem_ids), dtype=bool)
        else:
            rows = self.connectivity_list[elem_ids]
            tri_mask = (
                (rows[:, 0] == rows[:, 1])
                | (rows[:, 1] == rows[:, 2])
                | (rows[:, 2] == rows[:, 3])
                | (rows[:, 3] == rows[:, 0])
                | (rows[:, 0] == rows[:, 2])
                | (rows[:, 1] == rows[:, 3])
            )
            quad_mask = ~tri_mask

        # Calculate interior angles
        angles = self.interior_angles(elem_ids)
        
        # Initialize quality array
        quality = np.zeros(len(elem_ids))
        
        if quality_type in ['skew', 'skewness', 'angular skewness']:
            if tri_mask.any():
                tri_angles = angles[tri_mask, :3]
                tri_max = np.max(tri_angles, axis=1)
                tri_min = np.min(tri_angles, axis=1)
                quality[tri_mask] = 1 - np.maximum(
                    (tri_max - 60) / (180 - 60),
                    (60 - tri_min) / 60
                )

            if quad_mask.any():
                quad_angles = angles[quad_mask, :]
                quad_max = np.max(quad_angles, axis=1)
                quad_min = np.min(quad_angles, axis=1)
                quality[quad_mask] = 1 - np.maximum(
                    (quad_max - 90) / (180 - 90),
                    (90 - quad_min) / 90
                )

            # Zero out elements with degenerate angle sums
            tri_sum_mask = tri_mask & (np.sum(angles[:, :3], axis=1) <= 179.99)
            quality[tri_sum_mask] = 0
            quad_sum_mask = quad_mask & (np.sum(angles, axis=1) <= 359.99)
            quality[quad_sum_mask] = 0
        
        else:
            raise ValueError(f"Unknown quality type: {quality_type}")
        
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
    
    def _ordered_vertex_ring(self, v_idx: int, elem_ids: list) -> list | None:
        """
        Return CCW-ordered ring of neighbor vertices around interior vertex v_idx.

        Uses element connectivity to chain pred→succ pairs. Returns None if
        non-manifold or the vertex is on the boundary (open ring).
        """
        succ_map: dict[int, int] = {}

        for eid in elem_ids:
            row = self.connectivity_list[eid]
            # Determine unique vertex sequence for this element
            if row.shape[0] == 4 and row[3] == row[0]:
                verts = row[:3].tolist()          # padded triangle [v0,v1,v2,v0]
            elif row.shape[0] == 4 and (
                row[0] == row[1] or row[1] == row[2]
                or row[2] == row[3] or row[1] == row[3]
                or row[0] == row[2]
            ):
                seen: set[int] = set()
                verts = []
                for x in row.tolist():
                    if x not in seen:
                        seen.add(x)
                        verts.append(x)
            else:
                verts = row.tolist()

            n_local = len(verts)
            try:
                i = verts.index(v_idx)
            except ValueError:
                continue

            pred = verts[(i - 1) % n_local]
            succ = verts[(i + 1) % n_local]
            succ_map[pred] = succ

        if len(succ_map) != len(elem_ids):
            return None  # non-manifold

        start = next(iter(succ_map))
        ring = [start]
        cur = start
        for _ in range(len(succ_map) - 1):
            nxt = succ_map.get(cur)
            if nxt is None or nxt == start:
                break
            ring.append(nxt)
            cur = nxt

        if succ_map.get(ring[-1]) != start or len(ring) != len(succ_map):
            return None  # open ring → boundary vertex

        return ring

    def angle_based_smoother(self, n_iter: int = 100, omega: float = 0.5,
                              tol: float = 1e-8) -> np.ndarray:
        """
        Iterative angle-based smoother (Zhou & Shimada 2000).

        For each interior vertex, computes a bisector-weighted correction that
        drives each sector angle toward the equiangular target 2π/m.  Deficit
        is clamped to ±π/3 to prevent overshoot in near-degenerate sectors
        (where even a tiny bisector displacement causes huge angle change).
        Updates are Gauss-Seidel and accepted only when the local minimum-quality
        metric strictly improves, guaranteeing monotone quality growth.

        Unlike Laplacian smoothing, corrections are driven by angle-deficit per
        sector and are accepted only when they actually improve mesh quality.

        Parameters:
            n_iter: Maximum number of passes over all interior vertices
            omega:  Initial relaxation factor (halved up to 6× in line search)
            tol:    Convergence threshold on max per-vertex displacement

        Reference:
            Zhou, M., & Shimada, K. (2000).
            "An angle-based approach to two-dimensional mesh smoothing".
            Proceedings of the 9th International Meshing Roundtable, 373–384.
        """
        p = self.points[:, :2].copy()
        n = self.n_verts

        edge_verts = self.edge2vert(self.boundary_edges())
        boundary_set = set(np.unique(edge_verts.flatten()).tolist())

        vert2elem = self.adjacencies['Vert2Elem']
        two_pi = 2.0 * np.pi
        deficit_cap = np.pi / 3.0   # 60° cap prevents wild corrections in thin sectors

        def _elem_verts(row):
            """Unique ordered vertices for a connectivity row (handles padded tris)."""
            if row.shape[0] == 4 and row[3] == row[0]:
                return row[:3].tolist()
            if row.shape[0] == 4 and (
                row[0] == row[1] or row[1] == row[2]
                or row[2] == row[3] or row[1] == row[3]
                or row[0] == row[2]
            ):
                seen: set[int] = set()
                out = []
                for x in row.tolist():
                    if x not in seen:
                        seen.add(x)
                        out.append(x)
                return out
            return row.tolist()

        import math as _math
        _R2D = 180.0 / _math.pi

        def _acos_deg(px0, py0, px1, py1, px2, py2):
            """Angle at vertex 0 in degrees, using fast scalar math."""
            ux = px1 - px0; uy = py1 - py0
            wx = px2 - px0; wy = py2 - py0
            lu2 = ux*ux + uy*uy
            lw2 = wx*wx + wy*wy
            if lu2 < 1e-28 or lw2 < 1e-28:
                return 0.0
            c = (ux*wx + uy*wy) / _math.sqrt(lu2 * lw2)
            return _math.acos(max(-1.0, min(1.0, c))) * _R2D

        def _local_min_quality(verts_lists, v_idx, vpos, p_cur):
            """Min angular-skewness quality over incident elements (matches elem_quality)."""
            vpx, vpy = float(vpos[0]), float(vpos[1])
            min_q = 1e9
            for verts in verts_lists:
                n_v = len(verts)
                px = [vpx if vi == v_idx else float(p_cur[vi][0]) for vi in verts]
                py = [vpy if vi == v_idx else float(p_cur[vi][1]) for vi in verts]

                if n_v == 3:
                    # Area check
                    dx1 = px[1]-px[0]; dy1 = py[1]-py[0]
                    dx2 = px[2]-px[0]; dy2 = py[2]-py[0]
                    if dx1*dy2 - dy1*dx2 <= 0.0:
                        return -1.0
                    a0 = _acos_deg(px[0],py[0], px[1],py[1], px[2],py[2])
                    a1 = _acos_deg(px[1],py[1], px[0],py[0], px[2],py[2])
                    a2 = 180.0 - a0 - a1
                    amax = max(a0, a1, a2); amin = min(a0, a1, a2)
                    q = 1.0 - max((amax - 60.0) / 120.0, (60.0 - amin) / 60.0)
                else:
                    # Quad: check two sub-triangle areas
                    dx1 = px[1]-px[0]; dy1 = py[1]-py[0]
                    dx2 = px[2]-px[0]; dy2 = py[2]-py[0]
                    if dx1*dy2 - dy1*dx2 <= 0.0:
                        return -1.0
                    dx1 = px[2]-px[0]; dy1 = py[2]-py[0]
                    dx2 = px[3]-px[0]; dy2 = py[3]-py[0]
                    if dx1*dy2 - dy1*dx2 <= 0.0:
                        return -1.0
                    a0 = _acos_deg(px[0],py[0], px[3],py[3], px[1],py[1])
                    a1 = _acos_deg(px[1],py[1], px[0],py[0], px[2],py[2])
                    a2 = _acos_deg(px[2],py[2], px[1],py[1], px[3],py[3])
                    a3 = 360.0 - a0 - a1 - a2
                    amax = max(a0, a1, a2, a3); amin = min(a0, a1, a2, a3)
                    q = 1.0 - max((amax - 90.0) / 90.0, (90.0 - amin) / 90.0)
                if q < min_q:
                    min_q = q

            return min_q

        for _iter in range(n_iter):
            max_move = 0.0

            for v in range(n):
                if v in boundary_set:
                    continue

                elem_ids = list(vert2elem[v])
                if not elem_ids:
                    continue

                ring = self._ordered_vertex_ring(v, elem_ids)
                if ring is None or len(ring) < 2:
                    continue

                m_ring = len(ring)
                theta_star = two_pi / m_ring
                v_pos = p[v]
                correction = np.zeros(2)

                for k in range(m_ring):
                    a = p[ring[k]]
                    b = p[ring[(k + 1) % m_ring]]

                    da = a - v_pos
                    db = b - v_pos
                    la = np.linalg.norm(da)
                    lb = np.linalg.norm(db)
                    if la < 1e-14 or lb < 1e-14:
                        continue

                    ua = da / la
                    ub = db / lb
                    cos_a = np.clip(ua @ ub, -1.0, 1.0)
                    alpha_k = np.arccos(cos_a)

                    bisector = ua + ub
                    bl = np.linalg.norm(bisector)
                    bisector = bisector / bl if bl > 1e-10 else np.array([-ua[1], ua[0]])

                    # Cap deficit: prevents wild corrections for near-degenerate sectors
                    deficit = np.clip(theta_star - alpha_k, -deficit_cap, deficit_cap)
                    avg_len = (la + lb) * 0.5
                    correction += deficit * avg_len * bisector

                if np.linalg.norm(correction) < 1e-14:
                    continue

                ev_lists = [_elem_verts(self.connectivity_list[eid]) for eid in elem_ids]
                current_q = _local_min_quality(ev_lists, v, v_pos, p)

                # Line search: accept only when local minimum quality strictly improves
                step = omega * correction / m_ring
                scale = 1.0
                accepted = False
                for _ in range(6):
                    candidate = v_pos + scale * step
                    new_q = _local_min_quality(ev_lists, v, candidate, p)
                    if new_q > current_q:
                        accepted = True
                        break
                    scale *= 0.5

                if accepted:
                    p[v] = candidate   # Gauss-Seidel: update immediately
                    move = np.linalg.norm(scale * step)
                    if move > max_move:
                        max_move = move

            if max_move < tol:
                break

        new_points = np.zeros_like(self.points)
        new_points[:, :2] = p
        new_points[:, 2] = self.points[:, 2]
        return new_points

    def _detect_element_types(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Detect element types from connectivity_list.

        Returns:
            (tris, quads): Arrays of triangle and quad element indices.
            - tris: Indices where elements have 3 unique vertices (padded or 3-column)
            - quads: Indices where elements have 4 distinct vertices
        """
        n_cols = self.connectivity_list.shape[1]

        if n_cols == 3:
            # Pure triangle mesh
            all_indices = np.arange(self.connectivity_list.shape[0])
            return all_indices, np.array([], dtype=int)
        elif n_cols == 4:
            # Mixed or pure quad mesh; detect padded triangles by repeated vertex.
            # Padded convention is [v0, v1, v2, v0], so rows[:, 3] == rows[:, 0]
            # for triangles. Use same vectorised logic as _elem_type() for consistency.
            rows = self.connectivity_list
            tri_mask = (
                (rows[:, 0] == rows[:, 1])
                | (rows[:, 1] == rows[:, 2])
                | (rows[:, 2] == rows[:, 3])
                | (rows[:, 3] == rows[:, 0])
                | (rows[:, 0] == rows[:, 2])
                | (rows[:, 1] == rows[:, 3])
            )
            all_indices = np.arange(len(rows))
            return all_indices[tri_mask], all_indices[~tri_mask]
        else:
            raise ValueError(f"Unexpected connectivity_list shape: {self.connectivity_list.shape}")

    def _tri_stiffness_assembly(self, tri_indices: np.ndarray, p: np.ndarray, n: int) -> tuple:
        """
        Assemble stiffness matrix contributions from triangle elements.

        Parameters:
            tri_indices: Array of triangle element indices
            p: (n_verts, 2) point array (x, y only)
            n: Total number of vertices

        Returns:
            (rows, cols, data): CSR matrix format lists for assembly
        """
        import numpy as np

        t = self.connectivity_list[tri_indices, :3]

        D = 2.0 * np.eye(2)
        T = np.array([[-1.0, -np.sqrt(3)], [np.sqrt(3), -1.0]])

        rows, cols, data = [], [], []
        for elem_idx, tri in enumerate(t):
            for i in range(3):
                for j in range(3):
                    block = D if i == j else T if j == (i+1)%3 else T.T
                    for di in range(2):
                        for dj in range(2):
                            rows.append(2*tri[i] + di)
                            cols.append(2*tri[j] + dj)
                            data.append(block[di, dj])

        return rows, cols, data

    def _quad_stiffness_assembly(self, quad_indices: np.ndarray, p: np.ndarray, n: int) -> tuple:
        """
        Assemble stiffness matrix contributions from quad elements.

        Averages both diagonal decompositions (0→2 and 1→3) so every quad vertex
        receives equal stiffness weighting (1.5×D each). Single-diagonal decomposition
        double-counts diagonal vertices (2×D vs 1×D), creating asymmetric forces at
        quad-triangle seams that cause element collapse in mixed meshes.

        Total stiffness magnitude unchanged: 4 triangles × 0.5 weight = 2 triangle
        contributions, same as the single-diagonal form.

        Parameters:
            quad_indices: Array of quad element indices
            p: (n_verts, 2) point array (x, y only)
            n: Total number of vertices

        Returns:
            (rows, cols, data): CSR matrix format lists for assembly
        """
        import numpy as np

        q = self.connectivity_list[quad_indices, :4]

        D = 2.0 * np.eye(2)
        T = np.array([[-1.0, -np.sqrt(3)], [np.sqrt(3), -1.0]])

        rows, cols, data = [], [], []
        for quad in q:
            # Average both diagonals so all 4 vertices get equal stiffness (1.5*D each).
            # Diagonal 0-2: (0,1,2) + (0,2,3); diagonal 1-3: (0,1,3) + (1,2,3).
            for diag_tris in [
                [(quad[0], quad[1], quad[2]), (quad[0], quad[2], quad[3])],
                [(quad[0], quad[1], quad[3]), (quad[1], quad[2], quad[3])],
            ]:
                for tri in diag_tris:
                    for i in range(3):
                        for j in range(3):
                            block = D if i == j else T if j == (i + 1) % 3 else T.T
                            for di in range(2):
                                for dj in range(2):
                                    rows.append(2 * tri[i] + di)
                                    cols.append(2 * tri[j] + dj)
                                    data.append(0.5 * block[di, dj])

        return rows, cols, data

    def _mixed_stiffness_assembly(self, tri_indices: np.ndarray, quad_indices: np.ndarray,
                                   p: np.ndarray, n: int) -> tuple:
        """
        Assemble stiffness matrix for mixed element mesh.
        Combines triangle and quad contributions into single stiffness matrix.

        Parameters:
            tri_indices: Array of triangle element indices
            quad_indices: Array of quad element indices
            p: (n_verts, 2) point array (x, y only)
            n: Total number of vertices

        Returns:
            (rows, cols, data): CSR matrix format lists for assembly
        """
        rows, cols, data = [], [], []

        # Add triangle contributions
        if len(tri_indices) > 0:
            tri_rows, tri_cols, tri_data = self._tri_stiffness_assembly(tri_indices, p, n)
            rows.extend(tri_rows)
            cols.extend(tri_cols)
            data.extend(tri_data)

        # Add quad contributions
        if len(quad_indices) > 0:
            quad_rows, quad_cols, quad_data = self._quad_stiffness_assembly(quad_indices, p, n)
            rows.extend(quad_rows)
            cols.extend(quad_cols)
            data.extend(quad_data)

        return rows, cols, data

    def _compute_angle_based_forces(self, p: np.ndarray, n: int) -> np.ndarray:
        """
        Compute angle-based RHS force vector for Zhou & Shimada smoother.

        For each interior vertex, computes forces that pull it toward angles that
        conform to ideal values (60° for triangles, 90° for quads). Uses cotangent
        weighting: F_i = -0.5 * Σ_{adjacent angles} cot(θ) * (edge_perp)

        Parameters:
            p: (n_verts, 2) point array (x, y only)
            n: Total number of vertices

        Returns:
            F: (2*n,) force vector with angle-based forces for interior vertices
        """
        import numpy as np

        F = np.zeros(2 * n)

        # Process each element to accumulate angle-based forces
        for elem_idx, elem in enumerate(self.connectivity_list):
            # Determine if triangle or quad
            is_tri = (elem[3] == elem[0]) if len(elem) == 4 else True

            if is_tri:
                verts = elem[:3]
                ideal_angle = np.pi / 3  # 60° for triangles
            else:
                verts = elem[:4]
                ideal_angle = np.pi / 2  # 90° for quads

            # Compute angles at each vertex
            for i in range(len(verts)):
                v0 = int(verts[i])
                v1 = int(verts[(i - 1) % len(verts)])
                v2 = int(verts[(i + 1) % len(verts)])

                # Edge vectors emanating from v0
                e1 = p[v1] - p[v0]  # to previous vertex
                e2 = p[v2] - p[v0]  # to next vertex

                # Angle at v0
                cos_angle = np.dot(e1, e2) / (np.linalg.norm(e1) * np.linalg.norm(e2) + 1e-14)
                cos_angle = np.clip(cos_angle, -1.0, 1.0)
                angle = np.arccos(cos_angle)

                # Cotangent weight (negative because we want to penalize deviations)
                cot_angle = 1.0 / (np.tan(angle) + 1e-14)

                # Force is perpendicular to edges, weighted by cotangent
                # Use perpendicular to e2: rotate 90° clockwise = (y, -x)
                e2_perp = np.array([e2[1], -e2[0]])
                force = -0.5 * cot_angle * e2_perp

                # Accumulate force (scaled by 0.1 to avoid over-correction)
                F[2*v0] += 0.1 * force[0]
                F[2*v0 + 1] += 0.1 * force[1]

        return F

    def direct_smoother( self, kinf=1e12 ) -> np.ndarray:
        """
        Perform direct (non-iterative) FEM smoothing with fixed boundary nodes.
        Supports triangle, quad, and mixed-element meshes.

        Based on the Zhou & Shimada formulation (triangles) extended to quads
        via analogy (maintaining energy-minimization principle).

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
        n = self.n_verts

        # Detect element types
        tri_indices, quad_indices = self._detect_element_types()

        # Assemble stiffness matrix based on element composition
        if len(quad_indices) == 0:
            # Pure triangle mesh
            rows, cols, data = self._tri_stiffness_assembly(tri_indices, p, n)
        elif len(tri_indices) == 0:
            # Pure quad mesh
            rows, cols, data = self._quad_stiffness_assembly(quad_indices, p, n)
        else:
            # Mixed-element mesh
            rows, cols, data = self._mixed_stiffness_assembly(tri_indices, quad_indices, p, n)

        K = csr_matrix((data, (rows, cols)), shape=(2*n, 2*n))

        # Compute angle-based RHS forces for interior nodes
        F = self._compute_angle_based_forces(p, n)

        # Identify boundary nodes and apply constraints.
        # Use adjacency structure when available; otherwise count edge occurrences
        # directly (supports compute_layers=False meshes used in tests).
        if "Edge2Elem" in self.adjacencies:
            edge_verts = self.edge2vert(self.boundary_edges())
            boundary_nodes = np.unique(edge_verts.flatten())
        else:
            edge_count: dict[tuple, int] = {}
            for row in self.connectivity_list:
                # strip padding for degenerate-quad triangles
                verts = list(row[:3]) if (len(row) == 4 and row[3] == row[0]) else list(row)
                for i in range(len(verts)):
                    a, b = int(verts[i]), int(verts[(i + 1) % len(verts)])
                    key = (min(a, b), max(a, b))
                    edge_count[key] = edge_count.get(key, 0) + 1
            boundary_nodes = np.unique([v for key, cnt in edge_count.items() if cnt == 1 for v in key])

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