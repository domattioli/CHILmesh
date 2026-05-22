"""Quad-edge topology backend for CHILmesh.

This module implements a quad-edge (4-connected) data structure as an alternative
to the current EdgeMap + dict-based adjacency model. The quad-edge representation
stores directed edges with 4 pointers: origin vertex, next-clockwise neighbor,
next-counter-clockwise neighbor, and opposite-edge index.

Schema (QuadEdgeTopology.edges):
    ndarray[n_edges, 4] with columns:
    - [0] origin_vertex: int — vertex ID at the start of this directed edge
    - [1] next_cw: int — index of next edge clockwise around the incident face
    - [2] next_ccw: int — index of next edge counter-clockwise around the opposite face
    - [3] opposite_idx: int — index of the opposite directed edge (-1 for boundary)

Sentinel conventions:
    - Boundary edges: next_ccw == -1 and opposite_idx == -1
    - No per-object Python overhead; all lookups are array indexing

Padding behavior (mixed-element support):
    - Triangles have elem2vert shape [n_elems, 3] (no padding in current CHILmesh)
    - Mixed-element support via element-type detection: if elem2vert[i, 3] == elem2vert[i, 0], treat as triangle
    - No spurious edges created for padded vertices

Design notes:
    - Directional: 2 edges per undirected pair (origin, destination vs. destination, origin)
    - Undirected edges extracted via canonical sorting in to_edge2vert()
    - 2-phase O(n) construction: create edges, pair opposites + assign next pointers
    - Based on Wikipedia quad-edge definition, adapted for 2D CHILmesh constraints

Example construction:
    >>> from chilmesh.examples import annulus
    >>> mesh = annulus()
    >>> elem2vert = mesh.adjacencies['Elem2Vert']
    >>> n_verts = elem2vert.max() + 1
    >>> qe_topo = build_quadegg_from_connectivity(elem2vert, n_verts)
    >>> edge2vert = qe_topo.to_edge2vert()  # Convert to standard adjacency format
"""

from __future__ import annotations

import numpy as np
from typing import Dict, List, Set, Tuple


class QuadEdgeTopology:
    """Container for all quad-edges and conversion methods.

    Equivalent in role to EdgeMap and HalfEdgeTopology; emits the same
    dict/ndarray adjacency outputs consumed by the rest of CHILmesh.
    """

    def __init__(self, edges: np.ndarray, n_verts: int, elem2vert: np.ndarray):
        """Initialize QuadEdgeTopology.

        Args:
            edges: ndarray[n_edges, 4] with (origin, next_cw, next_ccw, opposite_idx)
            n_verts: Total number of vertices
            elem2vert: Original element-to-vertex connectivity (for reference)
        """
        self.edges = edges
        self.n_verts = n_verts
        self.elem2vert = elem2vert
        self.n_elems = elem2vert.shape[0]

    def to_edge2vert(self) -> np.ndarray:
        """Extract Edge2Vert: sorted undirected edges.

        Returns unique undirected edges in canonical form (min_vert, max_vert),
        sorted lexicographically to match EdgeMap output.

        Built from element connectivity to avoid quad-edge traversal complexity.
        """
        edges_set = set()
        elem_cols = self.elem2vert.shape[1]
        for elem_idx in range(self.n_elems):
            elem_verts = self.elem2vert[elem_idx]
            # Detect element type: if elem_cols >= 4 and last vert is padding, triangle; else quad
            if elem_cols == 3:
                elem_type = 3
            else:
                elem_type = 3 if elem_verts[3] == elem_verts[0] else 4

            for i in range(elem_type):
                v0 = int(elem_verts[i])
                v1 = int(elem_verts[(i + 1) % elem_type])
                edge = tuple(sorted([v0, v1]))
                edges_set.add(edge)

        edges_list = sorted(edges_set)
        return np.array(edges_list, dtype=np.int32)

    def to_edge2elem(self) -> np.ndarray:
        """Extract Edge2Elem: edge to adjacent elements.

        Returns ndarray[n_edges, 2] with element indices (-1 for boundary).
        Convention: [elem1, elem2] where elem1 comes first in iteration order.

        Builds mapping from element connectivity directly, independent of
        quad-edge traversal (which may be expensive or error-prone).
        """
        # Build canonical edge list (same order as to_edge2vert)
        edges = self._to_canonical_edge_list()

        # Build undirected-edge-to-elements map
        # Use list to preserve order: first element encountered goes to position 0
        edge_to_elems_dict = {tuple(e): [] for e in edges}

        # Iterate elements to build mapping (preserves order)
        elem_cols = self.elem2vert.shape[1]
        for elem_idx in range(self.n_elems):
            elem_verts = self.elem2vert[elem_idx]
            elem_type = 3 if (elem_cols == 3 or elem_verts[3] == elem_verts[0]) else 4

            for i in range(elem_type):
                v1 = int(elem_verts[i])
                v2 = int(elem_verts[(i + 1) % elem_type])
                undirected_edge = tuple(sorted([v1, v2]))

                if undirected_edge in edge_to_elems_dict:
                    if len(edge_to_elems_dict[undirected_edge]) < 2:
                        edge_to_elems_dict[undirected_edge].append(elem_idx)

        # Convert to result array with padding for boundary edges
        result = []
        for edge in edges:
            elems = edge_to_elems_dict[edge]
            # Pad with -1 to get 2 elements
            while len(elems) < 2:
                elems.append(-1)
            result.append(elems[:2])

        return np.array(result, dtype=np.int32)

    def to_elem2edge(self) -> np.ndarray:
        """Extract Elem2Edge: element to incident edge indices.

        Returns ndarray[n_elems, 3|4] with edge IDs for each element.
        """
        edges = self._to_canonical_edge_list()
        edge_list_sorted = sorted(edges)
        edge_to_id = {tuple(e): i for i, e in enumerate(edge_list_sorted)}

        elem2edge = []
        for elem_idx in range(self.n_elems):
            elem_edges = []
            elem_verts = self.elem2vert[elem_idx]
            elem_cols = self.elem2vert.shape[1]
            elem_type = 3 if (elem_cols == 3 or elem_verts[3] == elem_verts[0]) else 4

            for j in range(elem_type):
                v1 = int(elem_verts[j])
                v2 = int(elem_verts[(j + 1) % elem_type])
                edge = tuple(sorted([v1, v2]))
                edge_id = edge_to_id.get(edge, -1)
                elem_edges.append(edge_id)

            elem2edge.append(elem_edges)

        max_edges = max(len(e) for e in elem2edge)
        result = np.full((self.n_elems, max_edges), -1, dtype=np.int32)
        for i, edges in enumerate(elem2edge):
            result[i, : len(edges)] = edges

        return result

    def to_vert2edge(self) -> Dict[int, Set[int]]:
        """Extract Vert2Edge: vertex to incident edge indices.

        Returns dict mapping vertex ID to set of edge IDs.
        """
        edges = self._to_canonical_edge_list()
        edge_list_sorted = sorted(edges)
        edge_to_id = {tuple(e): i for i, e in enumerate(edge_list_sorted)}

        vert2edge = {v: set() for v in range(self.n_verts)}
        for edge in edge_list_sorted:
            edge_id = edge_to_id[tuple(edge)]
            vert2edge[edge[0]].add(edge_id)
            vert2edge[edge[1]].add(edge_id)

        return vert2edge

    def to_vert2elem(self) -> Dict[int, Set[int]]:
        """Extract Vert2Elem: vertex to incident elements.

        Returns dict mapping vertex ID to set of element IDs.
        """
        vert2elem = {v: set() for v in range(self.n_verts)}
        for elem_idx in range(self.n_elems):
            elem_verts = self.elem2vert[elem_idx]
            elem_cols = self.elem2vert.shape[1]
            elem_type = 3 if (elem_cols == 3 or elem_verts[3] == elem_verts[0]) else 4

            for i in range(elem_type):
                v = int(elem_verts[i])
                if v >= 0:
                    vert2elem[v].add(elem_idx)

        return vert2elem

    def to_edgemap_list(self) -> List[Tuple[int, int]]:
        """Extract edge list in canonical form (sorted tuples).

        Returns list of unique undirected edges as (min_vert, max_vert) tuples,
        sorted lexicographically to match EdgeMap.to_list() output.
        """
        return self._to_canonical_edge_list()

    def _to_canonical_edge_list(self) -> List[Tuple[int, int]]:
        """Extract edge list in canonical form (sorted tuples).

        Returns list of unique undirected edges as (min_vert, max_vert) tuples,
        sorted lexicographically to match EdgeMap.to_list() output.
        """
        edges_set = set()
        elem_cols = self.elem2vert.shape[1]
        for elem_idx in range(self.n_elems):
            elem_verts = self.elem2vert[elem_idx]
            elem_type = 3 if (elem_cols == 3 or elem_verts[3] == elem_verts[0]) else 4

            for i in range(elem_type):
                v0 = int(elem_verts[i])
                v1 = int(elem_verts[(i + 1) % elem_type])
                edge = tuple(sorted([v0, v1]))
                edges_set.add(edge)

        return sorted(edges_set)


def build_quadegg_from_connectivity(elem2vert: np.ndarray, n_verts: int) -> QuadEdgeTopology:
    """Construct QuadEdgeTopology from element-to-vertex connectivity.

    O(n) algorithm:
    1. Create all directed edges per element in iteration order
    2. Hash-table lookup to pair opposite edges and assign next pointers
    3. Single data-structure pass (no separate phases)

    Algorithm detail:
    - Each directed edge stores: origin vertex, next_cw (next edge clockwise in same face),
      next_ccw (next edge counter-clockwise in opposite face), opposite_idx
    - During creation, record (v_origin, v_dest) → edge_idx mapping
    - During pairing: for each edge, look up reverse (v_dest, v_origin) to find opposite
    - Simultaneously assign next_cw within face and next_ccw in opposite face
    - Boundary edges: opposite_idx = -1, next_ccw = -1

    Args:
        elem2vert: ndarray[n_elems, 3|4] — element connectivity (CCW-oriented)
        n_verts: Total number of vertices

    Returns:
        QuadEdgeTopology instance with O(n) construction cost
    """
    elem2vert = elem2vert.copy()
    n_elems = elem2vert.shape[0]
    elem_cols = elem2vert.shape[1]

    # Determine element types (3 for triangle, 4 for quad)
    elem_type = np.zeros(n_elems, dtype=np.int32)
    for i in range(n_elems):
        if elem_cols == 3:
            elem_type[i] = 3
        else:
            elem_type[i] = 4 if elem2vert[i, 3] != elem2vert[i, 0] else 3

    # Phase 1: Create all directed edges in element order
    edges_list = []
    edge_map = {}  # (v_origin, v_dest) -> edge_idx
    elem_to_edge_indices = {}  # elem_idx -> [edge_idx, edge_idx, ...]

    edge_idx_counter = 0

    for elem_idx in range(n_elems):
        elem_verts = elem2vert[elem_idx]
        n_verts_in_elem = int(elem_type[elem_idx])

        elem_to_edge_indices[elem_idx] = []

        for i in range(n_verts_in_elem):
            v_origin = int(elem_verts[i])
            v_dest = int(elem_verts[(i + 1) % n_verts_in_elem])

            # Create edge: [origin, next_cw, next_ccw, opposite_idx]
            # next_cw will be assigned in same element loop below
            # next_ccw and opposite_idx assigned during Phase 2
            edge = [v_origin, -1, -1, -1]
            edges_list.append(edge)
            elem_to_edge_indices[elem_idx].append(edge_idx_counter)

            # Record directed edge for opposite-pairing lookup
            directed_edge = (v_origin, v_dest)
            edge_map[directed_edge] = edge_idx_counter

            edge_idx_counter += 1

    edges_array = np.array(edges_list, dtype=np.int32)

    # Assign next_cw within each element's edges
    for elem_idx in range(n_elems):
        edge_indices = elem_to_edge_indices[elem_idx]
        n_edges_in_elem = len(edge_indices)
        for i in range(n_edges_in_elem):
            curr_edge_idx = edge_indices[i]
            next_edge_idx = edge_indices[(i + 1) % n_edges_in_elem]
            edges_array[curr_edge_idx, 1] = next_edge_idx

    # Phase 2: Pair opposite edges and assign next_ccw
    for (v_origin, v_dest), edge_idx in edge_map.items():
        reverse_edge = (v_dest, v_origin)
        if reverse_edge in edge_map:
            # Found opposite edge
            opposite_edge_idx = edge_map[reverse_edge]
            edges_array[edge_idx, 3] = opposite_edge_idx
            edges_array[opposite_edge_idx, 3] = edge_idx

            # Assign next_ccw pointers (going around opposite face)
            # next_ccw of current edge = next_cw of opposite edge
            next_cw_of_opposite = int(edges_array[opposite_edge_idx, 1])
            edges_array[edge_idx, 2] = next_cw_of_opposite

            # Symmetrically: next_ccw of opposite edge = next_cw of current edge
            next_cw_of_current = int(edges_array[edge_idx, 1])
            edges_array[opposite_edge_idx, 2] = next_cw_of_current
        else:
            # Boundary edge: opposite_idx = -1, next_ccw = -1
            edges_array[edge_idx, 3] = -1
            edges_array[edge_idx, 2] = -1

    return QuadEdgeTopology(edges_array, n_verts, elem2vert)
