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

    def submesh(
        self,
        elem_ids: Union[np.ndarray, List[int], Tuple[int, ...]],
        compute_layers: bool = True,
        compute_adjacencies: Opt[bool] = None,
    ) -> "CHILmesh":
        """
        Return a new :class:`CHILmesh` restricted to the given element indices (#138).

        The sub-mesh holds a compact copy of the selected elements: vertex indices
        are renumbered into [0, k) where k is the number of distinct vertices the
        selection references, and ``points`` is sliced accordingly. The original
        mesh is unchanged.

        Parameters:
            elem_ids: Indices of elements to include. May contain duplicates; the
                method deduplicates internally. Must be non-empty.
            compute_layers: Whether to skeletonize the sub-mesh (default: True).
                Forwarded to :class:`CHILmesh` constructor.
            compute_adjacencies: Whether to build adjacency dicts on the sub-mesh.
                ``None`` tracks ``compute_layers``. See :meth:`__init__`.

        Returns:
            A new :class:`CHILmesh` containing only the selected elements, with
            vertex indices remapped to a compact local numbering. ``grid_name`` is
            suffixed with ``"_submesh"`` when the parent grid had a name.

        Raises:
            ValueError: If ``elem_ids`` is empty or contains out-of-range indices.

        Notes:
            * Element ordering in the sub-mesh follows the **sorted** order of the
              unique input indices, not the input order. This keeps the operation
              deterministic regardless of caller iteration order.
            * For 4-column (quad or mixed) connectivity, padded triangle rows
              ``[v0, v1, v2, v0]`` remain padded after remapping (duplicates are
              preserved bit-for-bit through ``np.unique`` + dict lookup).
            * Boundary detection is recomputed from the sub-mesh topology: edges
              that were interior in the parent may become boundary edges in the
              sub-mesh when their partner element falls outside the selection.

        Example:
            >>> # Extract the outermost skeletonization layer as its own mesh.
            >>> outer = mesh.submesh(mesh.elements_in_layer(0))
            >>> outer.n_elems
            42
        """
        elem_ids = np.asarray(elem_ids, dtype=np.int64).ravel()
        if elem_ids.size == 0:
            raise ValueError("submesh requires at least one element id")
        if elem_ids.min() < 0 or elem_ids.max() >= self.n_elems:
            raise ValueError(
                f"elem_ids out of range [0, {self.n_elems}): "
                f"got min={int(elem_ids.min())}, max={int(elem_ids.max())}"
            )

        # Deduplicate and sort — np.unique gives both in one pass.
        unique_elems = np.unique(elem_ids)
        sub_conn = self.connectivity_list[unique_elems].copy()

        # -1 is the sentinel for triangle padding in 4-column mixed meshes.
        # Exclude it from the vertex remap; points[-1] is a real vertex and
        # must not be pulled into sub_points.
        valid_mask = sub_conn >= 0
        unique_verts = np.unique(sub_conn[valid_mask])
        remap_table = np.searchsorted(unique_verts, sub_conn)
        remap = np.where(valid_mask, remap_table, -1)

        sub_points = self.points[unique_verts].copy()

        sub_grid_name = (
            f"{self.grid_name}_submesh" if self.grid_name is not None else None
        )

        return CHILmesh(
            connectivity=remap,
            points=sub_points,
            grid_name=sub_grid_name,
            compute_layers=compute_layers,
            compute_adjacencies=compute_adjacencies,
        )

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

    def __init__( self, connectivity: Opt[np.ndarray] = None, points: Opt[np.ndarray] = None, grid_name: Opt[str] = None, compute_layers: bool = True, compute_adjacencies: Opt[bool] = None, seed_boundary_kinds: Opt[List[str]] = None, seed_ibtypes: Opt[List[int]] = None ) -> None:
        """
        Initialize a CHILmesh object.

        Parameters:
            connectivity: Element connectivity list (n_elems × 3 for triangles, n_elems × 4 for quads/mixed)
            points: Vertex coordinates (n_verts × 3, with z=0 for 2D meshes)
            grid_name: Name of the mesh
            compute_layers: If False, skip skeletonization for fast init (default: True).
                Skeletonization requires adjacencies, so when True, adjacencies are
                always built regardless of ``compute_adjacencies``.
            compute_adjacencies: Whether to build adjacency dicts (Vert2Edge, Vert2Elem,
                Edge2Vert, Edge2Elem, Elem2Edge, EdgeMap). If None (default), tracks
                ``compute_layers`` for backward compatibility. Set explicitly to True
                with ``compute_layers=False`` to obtain a mesh with usable adjacency
                queries (``get_vertex_edges``, ``edge2vert``, ``boundary_edges`` etc.)
                but no layer sweep — useful when downstream consumers need topology
                without paying the skeletonization cost (#134).
            seed_boundary_kinds: When provided, layer-0 peeling seeds only from boundary
                edges whose nodes belong to segments with a matching ``kind`` value
                (``'open'`` or ``'flow'``). ``None`` (default) uses all boundary edges
                — backward compatible. Ignored when ``compute_layers=False``. (#129)
            seed_ibtypes: When provided, layer-0 peeling seeds only from boundary edges
                whose nodes belong to segments with a matching ADCIRC IBTYPE value.
                Combined with ``seed_boundary_kinds`` via intersection. ``None``
                (default) applies no ibtype filter. Ignored when
                ``compute_layers=False``. (#129)

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
        # ADCIRC NOPE/NBOU boundary segments parsed from fort.14 (#129). Empty for
        # node+element-only meshes. Each entry: {"kind": "open"|"flow",
        # "ibtype": Optional[int], "nodes": np.ndarray of 0-based node indices}.
        self.boundary_segments: List[Dict[str, Any]] = []

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

        # Default compute_adjacencies follows compute_layers so old callers see
        # identical behavior; layer sweep cannot run without adjacencies, so
        # force-enable when layers are requested.
        if compute_adjacencies is None:
            compute_adjacencies = compute_layers
        if compute_layers:
            compute_adjacencies = True

        # Initialize the mesh
        self._initialize_mesh(
            compute_layers=compute_layers,
            compute_adjacencies=compute_adjacencies,
            seed_boundary_kinds=seed_boundary_kinds,
            seed_ibtypes=seed_ibtypes,
        )
    
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

    def _initialize_mesh( self, compute_layers: bool = True, compute_adjacencies: bool = True, seed_boundary_kinds: Opt[List[str]] = None, seed_ibtypes: Opt[List[int]] = None ) -> None:
        """Initialize the mesh properties.

        Parameters:
            compute_layers: If False, skip skeletonization (default: True)
            compute_adjacencies: If False, skip adjacency dict construction (default: True).
                Cannot be False when ``compute_layers`` is True — caller is responsible
                for that consistency (enforced in ``__init__``).
            seed_boundary_kinds: Forwarded to ``_skeletonize``; see ``__init__`` docs. (#129)
            seed_ibtypes: Forwarded to ``_skeletonize``; see ``__init__`` docs. (#129)
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

            if compute_adjacencies:
                self._build_adjacencies()
            if compute_layers:
                self._skeletonize(
                    seed_boundary_kinds=seed_boundary_kinds,
                    seed_ibtypes=seed_ibtypes,
                )

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
        elem_ids = np.asarray(elem_ids, dtype=np.intp)

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
        """Build adjacency lists for the mesh (vectorized)."""
        from .mesh_topology import EdgeMap

        tri_elems, quad_elems = self._elem_type()
        if len(quad_elems) == 0:
            self.type = "Triangular"
        elif len(tri_elems) == 0:
            self.type = "Quadrilateral"
        else:
            self.type = "Mixed-Element"
            if self.connectivity_list.shape[1] == 4:
                self.connectivity_list[tri_elems, 3] = self.connectivity_list[tri_elems, 0]

        conn = self.connectivity_list
        n_elems, n_cols = conn.shape[0], conn.shape[1]
        n_verts = self.n_verts

        # Build flat (elem, slot) → (v0, v1) arrays in element-major, slot-minor order
        # (same traversal order as the legacy Python loop) so first-encounter edge IDs
        # remain consistent with the C++ backend.
        v0 = conn.ravel(order='C').astype(np.int64)
        v1 = np.roll(conn, -1, axis=1).ravel(order='C').astype(np.int64)
        elem_flat = np.repeat(np.arange(n_elems, dtype=np.int64), n_cols)
        slot_flat = np.tile(np.arange(n_cols, dtype=np.int64), n_elems)

        valid = (v0 >= 0) & (v1 >= 0) & (v0 != v1)
        v0 = v0[valid]; v1 = v1[valid]
        elem_flat = elem_flat[valid]; slot_flat = slot_flat[valid]

        edges_canon = np.column_stack([np.minimum(v0, v1), np.maximum(v0, v1)])
        unique_lex, first_occur, inverse_lex = np.unique(
            edges_canon, axis=0, return_index=True, return_inverse=True
        )
        # Restore first-encounter ordering so edge IDs match the C++ backend.
        insert_order = np.argsort(first_occur, kind='stable')
        unique_edges = unique_lex[insert_order]
        lex_to_ins = np.empty(len(insert_order), dtype=np.int64)
        lex_to_ins[insert_order] = np.arange(len(insert_order), dtype=np.int64)
        inverse = lex_to_ins[inverse_lex]

        n_edges = len(unique_edges)
        self.n_edges = n_edges
        edge2vert = unique_edges

        elem2edge = np.zeros((n_elems, n_cols), dtype=int)
        elem2edge[elem_flat, slot_flat] = inverse

        edge2elem = np.full((n_edges, 2), -1, dtype=int)
        ord2 = np.argsort(inverse, kind='stable')
        si = inverse[ord2]; se = elem_flat[ord2]
        keep = np.ones(len(si), dtype=bool)
        keep[1:] = ~((si[1:] == si[:-1]) & (se[1:] == se[:-1]))
        si = si[keep]; se = se[keep]
        new_grp = np.concatenate([[True], si[1:] != si[:-1]])
        edge2elem[si, np.where(new_grp, 0, 1)] = se

        edge_map = EdgeMap()
        edge_map._map = {(int(r[0]), int(r[1])): i for i, r in enumerate(unique_edges)}
        edge_map._next_id = n_edges

        e_ids = np.arange(n_edges, dtype=np.int64)
        verts_e = np.concatenate([unique_edges[:, 0], unique_edges[:, 1]])
        edges_e = np.tile(e_ids, 2)
        oe = np.argsort(verts_e, kind='stable')
        sve = verts_e[oe]; see = edges_e[oe]
        bnde = np.concatenate([[0], np.where(sve[1:] != sve[:-1])[0] + 1, [len(sve)]])
        vert2edge: Dict[int, Set[int]] = {}
        for k in range(len(bnde) - 1):
            vert2edge[int(sve[bnde[k]])] = set(see[bnde[k]:bnde[k + 1]].tolist())
        for v in range(n_verts):
            if v not in vert2edge:
                vert2edge[v] = set()

        flat_v = conn.ravel().astype(np.int64)
        flat_e2 = np.repeat(np.arange(n_elems, dtype=np.int64), n_cols)
        ok = flat_v >= 0
        flat_v = flat_v[ok]; flat_e2 = flat_e2[ok]
        ov = np.argsort(flat_v, kind='stable')
        svv = flat_v[ov]; sev = flat_e2[ov]
        bndv = np.concatenate([[0], np.where(svv[1:] != svv[:-1])[0] + 1, [len(svv)]])
        vert2elem: Dict[int, Set[int]] = {}
        for k in range(len(bndv) - 1):
            vert2elem[int(svv[bndv[k]])] = set(sev[bndv[k]:bndv[k + 1]].tolist())
        for v in range(n_verts):
            if v not in vert2elem:
                vert2elem[v] = set()

        self.adjacencies = {
            "Elem2Vert": conn,
            "Edge2Vert": edge2vert,
            "EdgeMap": edge_map,
            "Elem2Edge": elem2edge,
            "Vert2Edge": vert2edge,
            "Vert2Elem": vert2elem,
            "Edge2Elem": edge2elem,
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

    def ccw_edges_around_vert(self, vert_id: int) -> List[int]:
        """
        Return the edge IDs incident to ``vert_id`` in counterclockwise order.

        Order is the geometric one in the XY plane: each incident edge is parameterised
        by ``atan2(dy, dx)`` from ``vert_id`` to its other endpoint, then sorted ascending
        (CCW from the positive x-axis). For a vertex with no incident edges this returns
        an empty list. Stable across runs given the same mesh (ties broken by edge ID).

        Parameters:
            vert_id: Vertex index [0, n_verts)

        Returns:
            List of edge IDs in CCW order around the vertex.

        Raises:
            ValueError: If vert_id is out of range
            RuntimeError: If adjacencies were not built (initialise with
                ``compute_adjacencies=True`` or ``compute_layers=True``)

        Complexity:
            Time: O(k log k) where k = degree of the vertex
            Space: O(k)

        Example:
            >>> # Walk neighbors of vertex 5 in CCW order
            >>> for edge_id in mesh.ccw_edges_around_vert(5):
            ...     v0, v1 = mesh.edge2vert(edge_id)
            ...     other = v1 if v0 == 5 else v0
        """
        if not 0 <= vert_id < self.n_verts:
            raise ValueError(f"Vertex {vert_id} out of range [0, {self.n_verts})")
        if "Vert2Edge" not in self.adjacencies or "Edge2Vert" not in self.adjacencies:
            raise RuntimeError(
                "Adjacencies not built. Re-initialise with compute_adjacencies=True "
                "(or compute_layers=True)."
            )

        edge_ids = list(self.adjacencies["Vert2Edge"][vert_id])
        if not edge_ids:
            return []

        edge2vert = self.adjacencies["Edge2Vert"]
        endpoints = edge2vert[edge_ids]
        # Pick the "other" endpoint per edge in one vectorised step.
        other = np.where(endpoints[:, 0] == vert_id, endpoints[:, 1], endpoints[:, 0])
        origin = self.points[vert_id, :2]
        deltas = self.points[other, :2] - origin
        angles = np.arctan2(deltas[:, 1], deltas[:, 0])
        order = np.lexsort((np.asarray(edge_ids), angles))
        return [int(edge_ids[i]) for i in order]

    def rebuild_adjacencies(self, rebuild_spatial_indices: bool = True) -> None:
        """Force a full rebuild of the adjacency cache from current connectivity.

        Use after a mid-sweep mutation of ``connectivity_list`` (or ``points``)
        where the cached adjacency dicts (``Edge2Vert``, ``Edge2Elem``,
        ``Vert2Edge``, ``Vert2Elem``, ``Elem2Edge``, ``EdgeMap``) are now stale.
        This is the public, non-private alternative to constructing a fresh
        ``CHILmesh(connectivity, points)`` and is the supported entry point for
        downstream consumers that mutate topology between adjacency queries.

        Parameters:
            rebuild_spatial_indices: When True (default), also rebuild the
                vertex / centroid KD-trees, which become stale when ``points``
                or ``connectivity_list`` change. Set False if only adjacency
                metadata changed without geometry shifts.

        Notes:
            - Does NOT recompute skeletonization layers. If the topology
              change altered the peel structure, also call ``_skeletonize()``
              or invalidate layers separately.
            - Preserves ``grid_name``, ``type``, and other metadata.
        """
        self._build_adjacencies()
        if rebuild_spatial_indices:
            self._build_spatial_indices()

    def invalidate_adjacencies(self) -> None:
        """Drop the cached adjacency dicts.

        Clears ``self.adjacencies`` so that subsequent adjacency-dependent
        accesses raise rather than return stale data. Call ``rebuild_adjacencies()``
        before any further adjacency queries. Spatial indices and layer caches
        are left untouched — invalidate them separately if needed.
        """
        self.adjacencies = {}
        self.n_edges = 0

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

    def _resolve_seed_nodes(self, seed_boundary_kinds: Opt[List[str]], seed_ibtypes: Opt[List[int]]) -> Opt[np.ndarray]:
        """Return the union of node arrays from boundary segments matching the filter (#129).

        Falls back to ``None`` (use all boundary edges) with a warning when
        ``boundary_segments`` is empty. Raises ``ValueError`` when segments exist but
        none match the filter.
        """
        if not self.boundary_segments:
            warnings.warn(
                "seed_boundary_kinds/seed_ibtypes specified but boundary_segments is "
                "empty; falling back to all boundary edges.",
                UserWarning,
                stacklevel=4,
            )
            return None
        matching = [
            seg for seg in self.boundary_segments
            if (seed_boundary_kinds is None or seg.get("kind") in seed_boundary_kinds)
            and (seed_ibtypes is None or seg.get("ibtype") in seed_ibtypes)
        ]
        if not matching:
            available_kinds = sorted({s.get("kind") for s in self.boundary_segments})
            available_ibtypes = sorted(
                {s.get("ibtype") for s in self.boundary_segments},
                key=lambda x: (x is None, x if x is not None else 0),
            )
            raise ValueError(
                f"No boundary segments match seed_boundary_kinds={seed_boundary_kinds!r}, "
                f"seed_ibtypes={seed_ibtypes!r}. "
                f"Available kinds: {available_kinds}, ibtypes: {available_ibtypes}"
            )
        return np.unique(np.concatenate([seg["nodes"] for seg in matching]))

    def _skeletonize(self, seed_boundary_kinds: Opt[List[str]] = None, seed_ibtypes: Opt[List[int]] = None) -> None:
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

        Parameters:
            seed_boundary_kinds: When not None, layer-0 seeds only from edges whose
                both vertices belong to boundary segments with matching ``kind``.
                ``None`` uses all boundary edges (default, MATLAB-parity). (#129)
            seed_ibtypes: When not None, layer-0 seeds only from edges whose both
                vertices belong to segments with matching IBTYPE. Combined with
                ``seed_boundary_kinds`` via intersection. (#129)
        """
        self.layers = {"OE": [], "IE": [], "OV": [], "IV": [], "bEdgeIDs": []}

        # Compute boundary-type seed node set for layer-0 filtering (#129).
        seed_nodes: Opt[np.ndarray] = None
        if seed_boundary_kinds is not None or seed_ibtypes is not None:
            seed_nodes = self._resolve_seed_nodes(seed_boundary_kinds, seed_ibtypes)

        edge2vert_work = self.adjacencies["Edge2Vert"].copy()
        edge2elem_work = self.adjacencies["Edge2Elem"].copy()

        iL = 0
        while np.any(edge2elem_work >= 0):
            # Step 1: boundary edges
            if iL == 0:
                iLbEdgeIDs = np.where(self.adjacencies["Edge2Elem"][:, 1] == -1)[0]
                if seed_nodes is not None and len(iLbEdgeIDs) > 0:
                    ev = self.adjacencies["Edge2Vert"][iLbEdgeIDs]  # (n, 2)
                    mask = np.isin(ev[:, 0], seed_nodes) & np.isin(ev[:, 1], seed_nodes)
                    iLbEdgeIDs = iLbEdgeIDs[mask]
                    if len(iLbEdgeIDs) == 0:
                        raise ValueError(
                            "seed_boundary_kinds/seed_ibtypes matched segments but no "
                            "boundary edges have both vertices in the seed node set. "
                            "Check that segment nodes lie on actual boundary edges."
                        )
            else:
                active_count = np.sum(edge2elem_work >= 0, axis=1)
                iLbEdgeIDs = np.where(active_count == 1)[0]

            if len(iLbEdgeIDs) == 0:
                break

            # Step 2: OV = unique vertices on boundary edges (from working copy)
            ov_raw = edge2vert_work[iLbEdgeIDs].ravel()
            iLOV = np.unique(ov_raw[ov_raw >= 0]).astype(int)
            self.layers["OV"].append(iLOV)
            self.layers["bEdgeIDs"].append(iLbEdgeIDs)

            # Step 3: OE = active elements adjacent to those boundary edges
            oe_raw = edge2elem_work[iLbEdgeIDs].ravel()
            iLOE = np.unique(oe_raw[oe_raw >= 0]).astype(int)
            self.layers["OE"].append(iLOE)

            # Step 4: consume OE from edge2elem_work
            if len(iLOE) > 0:
                edge2elem_work[np.isin(edge2elem_work, iLOE)] = -1

            # Step 5: find edges touching ANY OV vertex (broader than just boundary edges)
            ov_edge_mask = np.any(np.isin(edge2vert_work, iLOV), axis=1)
            ov_edge_indices = np.where(ov_edge_mask)[0]

            # Step 6: IE = active elements adjacent to those edges
            if len(ov_edge_indices) > 0:
                ie_raw = edge2elem_work[ov_edge_indices].ravel()
                iLIE = np.unique(ie_raw[ie_raw >= 0]).astype(int)
            else:
                iLIE = np.empty(0, dtype=int)
            self.layers["IE"].append(iLIE)

            # Step 7: consume OV vertices and IE elements
            if len(iLOV) > 0:
                edge2vert_work[np.isin(edge2vert_work, iLOV)] = -1
            if len(iLIE) > 0:
                edge2elem_work[np.isin(edge2elem_work, iLIE)] = -1

            # Step 8: IV = vertices of (OE ∪ IE) connectivity, minus OV
            if len(iLOE) > 0 or len(iLIE) > 0:
                layer_elems = np.concatenate((iLOE, iLIE))
                lv = self.connectivity_list[layer_elems].ravel()
                lv = lv[lv >= 0]
                iLIV = np.setdiff1d(np.unique(lv), iLOV).astype(int)
            else:
                iLIV = np.empty(0, dtype=int)
            self.layers["IV"].append(iLIV)

            iL += 1

        self.n_layers = iL

    def element_quality( self, metric: str = "aspect_ratio", elem_ids: Opt[Union[int, List[int], np.ndarray]] = None ) -> np.ndarray:
        """
        Compute element quality metric for elements in the mesh.

        Supported metrics:
            aspect_ratio (default): 2 * inradius / circumradius, range [0, 1].
                1 = equilateral, 0 = degenerate.
            min_angle: The minimum interior angle across all angles in the element (radians).
                Range [0, pi/3] for non-degenerate triangles.
            max_angle: The maximum interior angle across all angles in the element (radians).

        Parameters:
            metric: Quality metric to compute. One of 'aspect_ratio', 'min_angle', 'max_angle'.
            elem_ids: Indices of elements to evaluate.
                If None, all elements are evaluated.

        Returns:
            Array of quality values, one per element.

        Raises:
            ValueError: If metric is unknown.

        Example:
            >>> quality = mesh.element_quality()
            >>> print(f"Mean quality: {quality.mean():.3f}")
        """
        if metric not in ("aspect_ratio", "min_angle", "max_angle"):
            raise ValueError(f"Unknown metric: {metric!r}")

        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        if np.isscalar(elem_ids):
            elem_ids = np.array([elem_ids])
        elem_ids = np.asarray(elem_ids, dtype=int)

        qualities = np.zeros(len(elem_ids))
        for i, eid in enumerate(elem_ids):
            elem = self.connectivity_list[eid]
            verts = self.points[elem[:3], :2]  # triangles only for now
            a = verts[1] - verts[0]
            b = verts[2] - verts[1]
            c = verts[0] - verts[2]
            la, lb, lc = np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c)
            if metric == "aspect_ratio":
                s = (la + lb + lc) / 2
                area = abs(np.cross(verts[1] - verts[0], verts[2] - verts[0])) / 2
                if s <= 0 or la <= 0 or lb <= 0 or lc <= 0:
                    qualities[i] = 0.0
                else:
                    r_in = area / s
                    r_circ = (la * lb * lc) / (4 * area) if area > 0 else 0
                    qualities[i] = 2 * r_in / r_circ if r_circ > 0 else 0
            else:
                angles = np.array([
                    np.arccos(np.clip(np.dot(-c, a) / (lc * la + 1e-300), -1, 1)),
                    np.arccos(np.clip(np.dot(-a, b) / (la * lb + 1e-300), -1, 1)),
                    np.arccos(np.clip(np.dot(-b, c) / (lb * lc + 1e-300), -1, 1)),
                ])
                qualities[i] = angles.min() if metric == "min_angle" else angles.max()

        return qualities

    @staticmethod
    def read_from_msh(full_file_name: str, compute_layers: bool = False, compute_adjacencies: Opt[bool] = None) -> "CHILmesh":
        """
        Load a mesh from a Gmsh ASCII .msh file.

        Supports both format 2.2 and 4.1. Parses nodes and elements, supporting
        triangular and quadrilateral elements only. Mixed-element meshes use
        the padded-triangle convention in a 4-column array.

        Parameters:
            full_file_name: Path to the .msh file
            compute_layers: If False, skip skeletonization for fast init (default: False)
            compute_adjacencies: See ``CHILmesh.__init__``. If None (default), tracks
                ``compute_layers``.

        Returns:
            A CHILmesh object

        Raises:
            GmshParseError: If file format is unsupported or malformed.
        """
        from . import gmsh_io

        m = gmsh_io.read_msh(full_file_name)
        return CHILmesh(
            connectivity=m.connectivity_list,
            points=m.points,
            grid_name="Gmsh",
            compute_layers=compute_layers,
            compute_adjacencies=compute_adjacencies,
        )

    def write_to_msh(self, filename: str, grid_name: str = "CHILmesh Grid", version: str = "4.1") -> bool:
        """
        Export the current mesh to Gmsh ASCII .msh format.

        Supports format versions 2.2 and 4.1. Triangles in mixed-element meshes
        (padded convention) are written as 3-node triangles; quads as 4-node quads.

        Parameters:
            filename: Path to save the file
            grid_name: Optional title for the mesh
            version: Gmsh format version, "2.2" or "4.1" (default: "4.1")

        Returns:
            ``True`` on success.

        Raises:
            ValueError: If version is not "2.2" or "4.1".
        """
        from . import gmsh_io

        return gmsh_io.write_msh(
            filename,
            self.points,
            self.connectivity_list,
            grid_name,
            version,
        )

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

    def smooth(self, n_iter: int = 5, weight: float = 0.5, fixed_boundary: bool = True) -> None:
        """
        Apply Laplacian smoothing to the mesh (in-place).

        This method moves each interior vertex toward the average of its
        neighbors, optionally fixing boundary vertices. It modifies
        ``self.points`` in-place.

        Parameters:
            n_iter: Number of smoothing iterations (default: 5).
            weight: Relaxation weight in [0, 1] (default: 0.5).
                0 = no movement, 1 = full Laplacian step.
            fixed_boundary: If True (default), boundary vertices are not moved.

        Raises:
            ValueError: If weight is not in [0, 1].

        Note:
            Smoothing may reduce element quality for highly non-convex domains.
            Use small ``n_iter`` and ``weight`` values for such cases.

        Example:
            >>> mesh.smooth(n_iter=10, weight=0.3)
            >>> print(mesh.element_quality().mean())
        """
        if not 0 <= weight <= 1:
            raise ValueError(f"weight must be in [0, 1], got {weight}")

        boundary_nodes = set(self.boundary_node_indices().tolist()) if fixed_boundary else set()

        for _ in range(n_iter):
            new_pts = self.points.copy()
            for v in range(self.n_verts):
                if v in boundary_nodes:
                    continue
                neighbors = set()
                for e_id in self.adjacencies['Vert2Edge'][v]:
                    v1, v2 = self.adjacencies['Edge2Vert'][e_id]
                    neighbors.add(v1 if v2 == v else v2)
                if not neighbors:
                    continue
                neighbor_pts = self.points[list(neighbors), :2]
                avg = neighbor_pts.mean(axis=0)
                new_pts[v, :2] = (1 - weight) * self.points[v, :2] + weight * avg
            self.points = new_pts

        # Invalidate spatial indices after moving points
        self._build_spatial_indices()

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

        Uses element connectivity to chain pred->succ pairs. Returns None if
        non-manifold or the vertex is on the boundary (open ring).
        """
        succ_map: dict[int, int] = {}

        for eid in elem_ids:
            row = self.connectivity_list[eid]
            if row.shape[0] == 4 and row[3] == row[0]:
                verts = row[:3].tolist()
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
            return None

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
            return None

        return ring

    def angle_based_smoother(self, n_iter: int = 100, omega: float = 0.5,
                              tol: float = 1e-8) -> np.ndarray:
        """
        Iterative angle-based smoother (Zhou & Shimada 2000).

        For each interior vertex, computes a bisector-weighted correction that
        drives each sector angle toward the equiangular target 2pi/m. Deficit
        is clamped to +/-pi/3 to prevent overshoot in near-degenerate sectors.
        Updates are Gauss-Seidel and accepted only when the local minimum-quality
        metric strictly improves, guaranteeing monotone quality growth.

        Parameters:
            n_iter: Maximum number of passes over all interior vertices
            omega:  Initial relaxation factor (halved up to 6x in line search)
            tol:    Convergence threshold on max per-vertex displacement

        Reference:
            Zhou, M., & Shimada, K. (2000).
            "An angle-based approach to two-dimensional mesh smoothing".
            Proceedings of the 9th International Meshing Roundtable, 373-384.
        """
        p = self.points[:, :2].copy()
        n = self.n_verts

        edge_verts = self.edge2vert(self.boundary_edges())
        boundary_set = set(np.unique(edge_verts.flatten()).tolist())

        vert2elem = self.adjacencies['Vert2Elem']
        two_pi = 2.0 * np.pi
        deficit_cap = np.pi / 3.0

        def _elem_verts(row):
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
            ux = px1 - px0; uy = py1 - py0
            wx = px2 - px0; wy = py2 - py0
            lu2 = ux*ux + uy*uy
            lw2 = wx*wx + wy*wy
            if lu2 < 1e-28 or lw2 < 1e-28:
                return 0.0
            c = (ux*wx + uy*wy) / _math.sqrt(lu2 * lw2)
            return _math.acos(max(-1.0, min(1.0, c))) * _R2D

        def _local_min_quality_fast(v_idx, vpos, p_cur, cached_elem_verts):
            vpx, vpy = vpos[0], vpos[1]
            min_q = 1e9
            for ev in cached_elem_verts:
                n_v = len(ev)
                px = np.empty(n_v); py = np.empty(n_v)
                for i, vi in enumerate(ev):
                    if vi == v_idx:
                        px[i] = vpx; py[i] = vpy
                    else:
                        px[i] = p_cur[vi, 0]; py[i] = p_cur[vi, 1]

                if n_v == 3:
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
            elem_verts_cache = {}

            for v in range(n):
                if v in boundary_set:
                    continue

                elem_ids = list(vert2elem[v])
                if not elem_ids:
                    continue

                ring = self._ordered_vertex_ring(v, elem_ids)
                if ring is None or len(ring) < 2:
                    continue

                if v not in elem_verts_cache:
                    elem_verts_cache[v] = [_elem_verts(self.connectivity_list[eid]) for eid in elem_ids]
                cached_ev = elem_verts_cache[v]

                m_ring = len(ring)
                theta_star = two_pi / m_ring
                v_pos = p[v]

                ring_arr = np.array(ring)
                ring_next = np.roll(ring_arr, -1)

                a_vecs = p[ring_arr]
                b_vecs = p[ring_next]

                da = a_vecs - v_pos
                db = b_vecs - v_pos

                la = np.linalg.norm(da, axis=1)
                lb = np.linalg.norm(db, axis=1)

                valid = (la >= 1e-14) & (lb >= 1e-14)

                ua = np.zeros_like(da)
                ub = np.zeros_like(db)
                ua[valid] = da[valid] / la[valid, np.newaxis]
                ub[valid] = db[valid] / lb[valid, np.newaxis]

                cos_a = np.sum(ua * ub, axis=1)
                cos_a = np.clip(cos_a, -1.0, 1.0)
                alpha_k = np.arccos(cos_a)

                bisector = ua + ub
                bl = np.linalg.norm(bisector, axis=1)

                normalized_bisector = np.zeros_like(bisector)
                large_bl = bl > 1e-10
                normalized_bisector[large_bl] = bisector[large_bl] / bl[large_bl, np.newaxis]
                normalized_bisector[~large_bl] = np.column_stack([-ua[~large_bl, 1], ua[~large_bl, 0]])

                deficits = np.clip(theta_star - alpha_k, -deficit_cap, deficit_cap)
                avg_lens = (la + lb) * 0.5

                deficits[~valid] = 0.0
                avg_lens[~valid] = 0.0

                correction = np.sum(deficits[:, np.newaxis] * avg_lens[:, np.newaxis] * normalized_bisector, axis=0)

                if np.linalg.norm(correction) < 1e-14:
                    continue

                current_q = _local_min_quality_fast(v, v_pos, p, cached_ev)

                step = omega * correction / m_ring
                scale = 1.0
                accepted = False
                for _ in range(6):
                    candidate = v_pos + scale * step
                    new_q = _local_min_quality_fast(v, candidate, p, cached_ev)
                    if new_q > current_q:
                        accepted = True
                        break
                    scale *= 0.5

                if accepted:
                    p[v] = candidate
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
        """
        n_cols = self.connectivity_list.shape[1]

        if n_cols == 3:
            all_indices = np.arange(self.connectivity_list.shape[0])
            return all_indices, np.array([], dtype=int)
        elif n_cols == 4:
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
        """Assemble stiffness matrix contributions from triangle elements (Balendran 1999)."""
        t = self.connectivity_list[tri_indices, :3]

        D = 2.0 * np.eye(2)
        T = np.array([[-1.0, -np.sqrt(3)], [np.sqrt(3), -1.0]])

        rows, cols, data = [], [], []
        for tri in t:
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
        """Assemble stiffness from quad elements via the Q4 bilinear Laplacian."""
        q = self.connectivity_list[quad_indices, :4]

        KL = (1.0 / 6.0) * np.array(
            [[4.0, -1.0, -2.0, -1.0],
             [-1.0, 4.0, -1.0, -2.0],
             [-2.0, -1.0, 4.0, -1.0],
             [-1.0, -2.0, -1.0, 4.0]]
        )

        rows, cols, data = [], [], []
        for quad in q:
            for i in range(4):
                for j in range(4):
                    k = KL[i, j]
                    if k == 0.0:
                        continue
                    vi, vj = int(quad[i]), int(quad[j])
                    for d in range(2):
                        rows.append(2 * vi + d)
                        cols.append(2 * vj + d)
                        data.append(k)

        return rows, cols, data

    def _mixed_stiffness_assembly(self, tri_indices: np.ndarray, quad_indices: np.ndarray,
                                   p: np.ndarray, n: int) -> tuple:
        """Assemble stiffness matrix for mixed element mesh."""
        rows, cols, data = [], [], []

        if len(tri_indices) > 0:
            tri_rows, tri_cols, tri_data = self._tri_stiffness_assembly(tri_indices, p, n)
            rows.extend(tri_rows)
            cols.extend(tri_cols)
            data.extend(tri_data)

        if len(quad_indices) > 0:
            quad_rows, quad_cols, quad_data = self._quad_stiffness_assembly(quad_indices, p, n)
            rows.extend(quad_rows)
            cols.extend(quad_cols)
            data.extend(quad_data)

        return rows, cols, data

    def direct_smoother(self, kinf=1e12, freeze_quad_nodes: bool = False) -> np.ndarray:
        """
        Perform direct (non-iterative) FEM smoothing with fixed boundary nodes.
        Supports triangle, quad, and mixed-element meshes.

        Triangles use the Balendran rotation-based stiffness (equilateral target);
        quads use the Q4 bilinear Laplacian (square target).

        Interior RHS F = 0 per Balendran/MATLAB FEMSmooth.m — only boundary
        pinning terms are non-zero. (#173 fix: removed incorrect
        _compute_angle_based_forces call from interior RHS.)

        Note (size-field behavior, #168): this smoother is **isotropic**. The
        Balendran stiffness targets a uniform equilateral triangle (60 deg) /
        square quad (90 deg) and takes **no size-field input** — it equalizes
        element *shape*, not *size*. Applied to a graded mesh it grows fine
        (e.g. coastal) edges and shrinks coarse (offshore) edges, eroding the
        original sizing. For size-respecting smoothing, supply anisotropic
        targets (not yet implemented) or run a separate sizing pass afterward.

        Note (size-field behavior, #168): this smoother is **isotropic**. The
        Balendran stiffness targets a uniform equilateral triangle (60 deg) /
        square quad (90 deg) and takes **no size-field input** — it equalizes
        element *shape*, not *size*. Applied to a graded mesh it grows fine
        (e.g. coastal) edges and shrinks coarse (offshore) edges, eroding the
        original sizing. For size-respecting smoothing, supply anisotropic
        targets (not yet implemented) or run a separate sizing pass afterward.

        Parameters:
            kinf: Large stiffness value for fixed (pinned) vertices.
            freeze_quad_nodes: When True, pin every vertex that is a corner of any
                quad element in addition to the boundary.

        Reference:
            Balendran, B. (1999).
            "A direct smoothing method for surface meshes".
            Proceedings of the 8th International Meshing Roundtable, 189-193.
        """
        from scipy.sparse import csr_matrix
        from scipy.sparse.linalg import spsolve

        p = self.points[:, :2]
        n = self.n_verts

        tri_indices, quad_indices = self._detect_element_types()

        if len(quad_indices) == 0:
            rows, cols, data = self._tri_stiffness_assembly(tri_indices, p, n)
        elif len(tri_indices) == 0:
            rows, cols, data = self._quad_stiffness_assembly(quad_indices, p, n)
        else:
            rows, cols, data = self._mixed_stiffness_assembly(tri_indices, quad_indices, p, n)

        K = csr_matrix((data, (rows, cols)), shape=(2*n, 2*n))

        # Interior RHS = 0 per Balendran/MATLAB FEMSmooth.m (#173 fix)
        F = np.zeros(2 * n)

        if "Edge2Elem" in self.adjacencies:
            edge_verts = self.edge2vert(self.boundary_edges())
            boundary_nodes = np.unique(edge_verts.flatten())
        else:
            edge_count: dict[tuple, int] = {}
            for row in self.connectivity_list:
                verts = list(row[:3]) if (len(row) == 4 and row[3] == row[0]) else list(row)
                for i in range(len(verts)):
                    a, b = int(verts[i]), int(verts[(i + 1) % len(verts)])
                    key = (min(a, b), max(a, b))
                    edge_count[key] = edge_count.get(key, 0) + 1
            boundary_nodes = np.unique([v for key, cnt in edge_count.items() if cnt == 1 for v in key])

        pinned = set(int(v) for v in boundary_nodes)
        if freeze_quad_nodes and len(quad_indices) > 0:
            quad_corner_nodes = np.unique(self.connectivity_list[quad_indices, :4].flatten())
            pinned |= set(int(v) for v in quad_corner_nodes)

        for v in pinned:
            F[2*v:2*v+2] = kinf * p[v]
            K[2*v, 2*v] = kinf
            K[2*v+1, 2*v+1] = kinf

        c = spsolve(K, F)
        new_xy = c.reshape(-1, 2)
        domain_diag = float(np.linalg.norm(np.ptp(p, axis=0)))
        max_disp = float(np.max(np.linalg.norm(new_xy - p, axis=1)))
        if not np.isfinite(c).all() or max_disp > domain_diag:
            return self.points.copy()
        new_points = np.zeros_like(self.points)
        new_points[:, :2] = new_xy
        new_points[:, 2] = self.points[:, 2]
        return new_points

    def admesh_metadata(self) -> dict:
        """
        Return a metadata dictionary compatible with the ADMESH-Domains catalog schema.

        The returned dict contains all fields that ADMESH-Domains expects from a
        mesh record: node count, element count, element type, and bounding box.
        Designed to be callable on a ``compute_layers=False`` mesh for fast bulk
        loading.

        Returns:
            dict with keys:
                - 'node_count' (int)
                - 'element_count' (int)
                - 'element_type' (str): 'tri', 'quad', or 'mixed'
                - 'bbox' (list[float]): [xmin, ymin, xmax, ymax]

        Example:
            >>> mesh = CHILmesh.from_admesh_domain(record)
            >>> info = mesh.admesh_metadata()
            >>> print(info['node_count'])
        """
        type_map = {"Triangular": "tri", "Quadrilateral": "quad", "Mixed-Element": "mixed"}
        xy = self.points[:, :2]
        bbox = [float(xy[:, 0].min()), float(xy[:, 1].min()),
                float(xy[:, 0].max()), float(xy[:, 1].max())]
        return {
            "node_count": self.n_verts,
            "element_count": self.n_elems,
            "element_type": type_map.get(self.type, "tri"),
            "bbox": bbox,
        }

    @classmethod
    def from_admesh_domain(cls, record: object, compute_layers: bool = True, compute_adjacencies: Opt[bool] = None) -> "CHILmesh":
        """
        Construct a CHILmesh from an ADMESH-Domains catalog record.

        The catalog record is duck-typed: any object with ``connectivity``
        (ndarray, n_elems × 3|4) and ``points`` (ndarray, n_verts × 2|3)
        attributes qualifies. This avoids a hard dependency on the
        ``admesh_domains`` package — consumers can pass any compatible object.

        Alternatively, if the record has a ``type`` attribute, routes to the
        appropriate file reader:
        - type == "SMS_2DM" (or "2dm"): calls read_from_2dm(record.filename)
        - type == "fort14", "FORT14", "ADCIRC": calls read_from_fort14(record.filename)

        Parameters:
            record: An object with ``.connectivity`` and ``.points`` attributes,
                or a ``type`` and ``filename`` attribute for file-based loading.
                Optionally ``.grid_name`` (str) for naming the mesh.
            compute_layers: If False, skip skeletonization for fast init.
            compute_adjacencies: Forwarded to the file reader / constructor. (#192)

        Returns:
            A CHILmesh instance.

        Example:
            >>> from admesh_domains import get_mesh
            >>> record = get_mesh("WNAT/hagen@v1")
            >>> mesh = CHILmesh.from_admesh_domain(record)
        """
        # Check for type-based routing (file-based loading)
        record_type = getattr(record, "type", None)
        if record_type is not None:
            record_type_lower = str(record_type).lower()
            if record_type_lower in ("sms_2dm", "2dm"):
                filename = getattr(record, "filename", None)
                if filename is not None:
                    return cls.read_from_2dm(
                        Path(filename), compute_layers=compute_layers, compute_adjacencies=compute_adjacencies
                    )
            elif record_type_lower in ("fort14", "adcirc"):
                filename = getattr(record, "filename", None)
                if filename is not None:
                    return cls.read_from_fort14(
                        Path(filename), compute_layers=compute_layers, compute_adjacencies=compute_adjacencies
                    )

        # Fall through to duck-typed path (connectivity + points attributes)
        name = getattr(record, "grid_name", None)
        return cls(
            connectivity=np.asarray(record.connectivity),
            points=np.asarray(record.points),
            grid_name=name,
            compute_layers=compute_layers,
            compute_adjacencies=compute_adjacencies,
        )

    def save(self, filename: str) -> None:
        """
        Save the mesh to a file.

        Supported formats: ADCIRC fort.14 (`.14`, `.fort14`), SMS 2dm (`.2dm`).

        Parameters:
            filename: Output file path (format inferred from extension).

        Raises:
            ValueError: If the file format is not recognized.

        Example:
            >>> mesh.save("output.14")
            >>> mesh.save("output.2dm")
        """
        path = Path(filename)
        if path.suffix in (".14", ".fort14"):
            write_fort14(self, filename)
        elif path.suffix == ".2dm":
            self._write_2dm(filename)
        else:
            raise ValueError(f"Unrecognized file format: {path.suffix!r}. Supported: .14, .fort14, .2dm")

    @classmethod
    def load(cls, filename: str, compute_layers: bool = True, compute_adjacencies: Opt[bool] = None) -> "CHILmesh":
        """
        Load a mesh from a file.

        Supported formats: ADCIRC fort.14 (`.14`, `.fort14`), SMS 2dm (`.2dm`).

        Parameters:
            filename: Input file path.
            compute_layers: If False, skip skeletonization for fast init.
            compute_adjacencies: Forwarded to the reader. (#192)

        Returns:
            A CHILmesh instance.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file format is not recognized.

        Example:
            >>> mesh = CHILmesh.load("input.14")
        """
        path = Path(filename)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {filename}")
        if path.suffix in (".14", ".fort14"):
            return cls.read_from_fort14(path, compute_layers=compute_layers, compute_adjacencies=compute_adjacencies)
        elif path.suffix == ".2dm":
            return cls.read_from_2dm(path, compute_layers=compute_layers, compute_adjacencies=compute_adjacencies)
        else:
            raise ValueError(f"Unrecognized file format: {path.suffix!r}. Supported: .14, .fort14, .2dm")

    def _write_2dm(self, filename: str) -> None:
        """Write mesh to SMS 2dm format."""
        with open(filename, 'w') as f:
            f.write("MESH2D\n")
            for i, elem in enumerate(self.connectivity_list):
                if elem.size == 3 or (elem.size == 4 and elem[3] == elem[0]):
                    verts = elem[:3]
                    f.write(f"E3T {i+1} {' '.join(str(v+1) for v in verts)}\n")
                else:
                    verts = elem[:4]
                    f.write(f"E4Q {i+1} {' '.join(str(v+1) for v in verts)}\n")
            for i, pt in enumerate(self.points):
                f.write(f"ND {i+1} {pt[0]:.10f} {pt[1]:.10f} {pt[2]:.10f}\n")

    @classmethod
    def read_from_2dm(cls, filename: Path, compute_layers: bool = True, compute_adjacencies: Opt[bool] = None) -> "CHILmesh":
        """
        Read mesh from SMS 2dm format.

        Parameters:
            filename: Path to the 2dm file.
            compute_layers: If False, skip skeletonization for fast init.
            compute_adjacencies: See :meth:`read_from_fort14`. (#192)

        Returns:
            A CHILmesh instance.
        """
        filename = Path(filename)
        nodes = {}
        elems = []
        has_mixed = False
        has_tri = False
        has_quad = False

        with open(filename) as f:
            for line in f:
                parts = line.split()
                if not parts:
                    continue
                if parts[0] == 'E3T':
                    # E3T: elem_id n1 n2 n3 material_id
                    # Take only vertex indices (parts[2:5]), skip material_id
                    vertices = [int(x) - 1 for x in parts[2:5]]
                    elems.append(vertices)
                    has_tri = True
                elif parts[0] == 'E4Q':
                    # E4Q: elem_id n1 n2 n3 n4 material_id
                    # Take only vertex indices (parts[2:6]), skip material_id
                    vertices = [int(x) - 1 for x in parts[2:6]]
                    elems.append(vertices)
                    has_quad = True
                elif parts[0] == 'ND':
                    nid = int(parts[1]) - 1
                    nodes[nid] = [float(x) for x in parts[2:]]

        if not nodes:
            raise ValueError("No nodes found in 2dm file")

        # Check if mixed-element mesh
        has_mixed = has_tri and has_quad

        n_verts = max(nodes.keys()) + 1
        pts = np.zeros((n_verts, 3))
        for nid, coords in nodes.items():
            pts[nid] = coords if len(coords) == 3 else coords + [0.0]

        # Handle padding for mixed-element meshes
        if has_mixed:
            # Pad triangles to 4 columns by repeating first vertex
            padded_elems = []
            for elem in elems:
                if len(elem) == 3:
                    padded_elems.append([elem[0], elem[1], elem[2], elem[0]])
                else:
                    padded_elems.append(elem)
            conn = np.array(padded_elems, dtype=int)
        else:
            conn = np.array(elems, dtype=int)
        return cls(connectivity=conn, points=pts,
                   grid_name=filename.stem, compute_layers=compute_layers,
                   compute_adjacencies=compute_adjacencies)

    @classmethod
    def read_from_fort14(
        cls,
        filename: Path,
        compute_layers: bool = True,
        compute_adjacencies: Opt[bool] = None,
    ) -> "CHILmesh":
        """
        Read mesh from ADCIRC fort.14 format.

        Parameters:
            filename: Path to the fort.14 file.
            compute_layers: If False, skip skeletonization for fast init.
            compute_adjacencies: Whether to build adjacency dicts. ``None``
                (default) tracks ``compute_layers``; set True with
                ``compute_layers=False`` for a fast adjacency-only load, or
                False with ``compute_layers=False`` for a bare conn+pts mesh.
                Forced True when ``compute_layers=True``. (#192)

        Returns:
            A CHILmesh instance with boundary_segments populated from
            NOPE (open) and NBOU (flow) boundary records.
        """
        filename = Path(filename)
        with open(filename) as f:
            lines = f.readlines()

        i = 0
        grid_name = lines[i].strip(); i += 1
        parts = lines[i].split(); i += 1
        n_elems, n_verts = int(parts[0]), int(parts[1])

        pts = np.zeros((n_verts, 3))
        for j in range(n_verts):
            p = lines[i].split(); i += 1
            pts[j] = [float(p[1]), float(p[2]), float(p[3]) if len(p) > 3 else 0.0]

        elem_verts = []
        has_tri = False
        has_quad = False
        for j in range(n_elems):
            p = lines[i].split(); i += 1
            n_verts_elem = int(p[1])
            verts = [int(x) - 1 for x in p[2:2 + n_verts_elem]]
            elem_verts.append(verts)
            if n_verts_elem == 3:
                has_tri = True
            elif n_verts_elem == 4:
                has_quad = True

        if has_tri and has_quad:
            # Mixed mesh: pad triangles to 4 cols by repeating first vertex.
            conn = np.array(
                [v if len(v) == 4 else [v[0], v[1], v[2], v[0]] for v in elem_verts],
                dtype=int,
            )
        else:
            conn = np.array(elem_verts, dtype=int)

        mesh = cls(
            connectivity=conn,
            points=pts,
            grid_name=grid_name,
            compute_layers=False,
            compute_adjacencies=False,
        )

        # --- parse boundary segments (#129) ---
        boundary_segments = []
        try:
            # NOPE open boundaries
            nope = int(lines[i].split()[0]); i += 1
            total_nope = int(lines[i].split()[0]); i += 1
            for _ in range(nope):
                n_seg = int(lines[i].split()[0]); i += 1
                nodes = []
                for _ in range(n_seg):
                    nodes.append(int(lines[i].split()[0]) - 1)
                    i += 1
                boundary_segments.append(
                    {"kind": "open", "ibtype": None, "nodes": np.array(nodes, dtype=int)}
                )
            # NBOU flow boundaries
            nbou = int(lines[i].split()[0]); i += 1
            total_nbou = int(lines[i].split()[0]); i += 1
            for _ in range(nbou):
                hdr = lines[i].split(); i += 1
                n_seg = int(hdr[0])
                ibtype = int(hdr[1]) if len(hdr) > 1 else None
                nodes = []
                for _ in range(n_seg):
                    nodes.append(int(lines[i].split()[0]) - 1)
                    i += 1
                boundary_segments.append(
                    {"kind": "flow", "ibtype": ibtype, "nodes": np.array(nodes, dtype=int)}
                )
        except (IndexError, ValueError):
            pass  # Boundary section absent or malformed — leave segments empty

        mesh.boundary_segments = boundary_segments

        # Default compute_adjacencies follows compute_layers (mirrors __init__).
        if compute_adjacencies is None:
            compute_adjacencies = compute_layers
        if compute_layers:
            compute_adjacencies = True

        # Now run adjacencies + skeletonization with segments in place.
        if compute_layers:
            mesh._build_adjacencies()
            mesh._skeletonize()
        elif compute_adjacencies:
            mesh._build_adjacencies()
        mesh._build_spatial_indices()

        return mesh


def write_fort14(mesh, filename: str) -> None:
    """
    Write a CHILmesh to ADCIRC fort.14 format.

    Parameters:
        mesh: CHILmesh object to save.
        filename: Output file path.

    Notes:
        Only node and element data is written. Boundary condition records
        are not included in this basic implementation.

    Example:
        >>> write_fort14(mesh, "output.14")
    """
    with open(filename, 'w') as f:
        name = mesh.grid_name or "CHILmesh"
        f.write(f"{name}\n")
        f.write(f"{mesh.n_elems} {mesh.n_verts}\n")
        for i, pt in enumerate(mesh.points):
            f.write(f"{i+1} {pt[0]:.10f} {pt[1]:.10f} {pt[2]:.10f}\n")
        for i, elem in enumerate(mesh.connectivity_list):
            if elem.size == 3 or (elem.size == 4 and elem[3] == elem[0]):
                verts = elem[:3]
                f.write(f"{i+1} 3 {' '.join(str(v+1) for v in verts)}\n")
            else:
                verts = elem[:4]
                f.write(f"{i+1} 4 {' '.join(str(v+1) for v in verts)}\n")


def _check_fort14(filename) -> bool:
    """Quick structural-validity check for a fort.14 file.

    Returns True if the file parses correctly (node count, element count,
    and grid name header all match), False otherwise.
    """
    try:
        CHILmesh.read_from_fort14(Path(filename), compute_layers=False)
        return True
    except Exception as e:
        print(f"Error reading fort14 file {filename}: {e}")
        return False