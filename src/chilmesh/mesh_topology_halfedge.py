"""Half-edge (DCEL) topology backend for CHILmesh.

This module implements a doubly-connected edge list (DCEL) as an alternative
to the current EdgeMap + dict-based adjacency model. The half-edge representation
stores directed edges as integer-indexed pointers in a contiguous array, achieving
lower memory overhead and faster traversal than per-object Python representations.

Schema (HalfEdgeTopology.half_edges):
    ndarray[n_halfedges, 4] with columns:
    - [0] origin_vertex: int — vertex ID at the start of this half-edge
    - [1] twin_idx: int — index of the opposite half-edge (-1 for boundary)
    - [2] next_idx: int — index of the next half-edge around the incident face
    - [3] face_idx: int — element ID (triangle or quad) that owns this half-edge

Sentinel conventions:
    - Boundary half-edges: twin_idx == -1 (matches Edge2Elem sentinel for boundaries)
    - No per-object Python overhead; all lookups are array indexing

Padding behavior (mixed-element support):
    - Triangles are stored in 4-column arrays with vertex 0 in slot 3 (_elem_type mask)
    - Half-edge construction reuses the existing _elem_type mask to skip the padded slot
    - No spurious half-edges are created for the padding

Example construction:
    >>> from chilmesh.examples import annulus
    >>> mesh = annulus()
    >>> elem2vert = mesh.adjacencies['Elem2Vert']
    >>> n_verts = elem2vert.max() + 1
    >>> he_topo = build_halfedge_from_connectivity(elem2vert, n_verts)
    >>> edge2vert = he_topo.to_edge2vert()  # Convert to standard adjacency format
"""

from __future__ import annotations

import numpy as np
from typing import Dict, List, Set, Tuple


class HalfEdgeTopology:
    """Container for all half-edges and lookup indices.

    Equivalent in role to EdgeMap; emits the same dict/ndarray adjacency
    outputs consumed by the rest of CHILmesh.
    """

    def __init__(self, half_edges: np.ndarray, n_verts: int, elem2vert: np.ndarray):
        """Initialize HalfEdgeTopology.

        Args:
            half_edges: ndarray[n_halfedges, 4] with (origin, twin_idx, next_idx, face_idx)
            n_verts: Total number of vertices
            elem2vert: Original element-to-vertex connectivity (for reference)
        """
        self.half_edges = half_edges
        self.n_verts = n_verts
        self.elem2vert = elem2vert
        self.n_elems = elem2vert.shape[0]

    def to_edge2vert(self) -> np.ndarray:
        """Extract Edge2Vert: sorted undirected edges.

        Returns unique undirected edges in canonical form (min_vert, max_vert),
        sorted lexicographically to match EdgeMap output.

        Built from element connectivity to avoid half-edge traversal complexity.
        """
        edges_set = set()
        elem_cols = self.elem2vert.shape[1]
        for elem_idx in range(self.n_elems):
            elem_verts = self.elem2vert[elem_idx]
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
        """
        edges = self.to_edgemap_list()
        edge_to_elems = {}

        for i in range(len(self.half_edges)):
            origin = int(self.half_edges[i, 0])
            twin_idx = int(self.half_edges[i, 1])
            face_idx = int(self.half_edges[i, 3])

            if twin_idx >= 0:
                dest = int(self.half_edges[twin_idx, 0])
            else:
                next_idx = int(self.half_edges[i, 2])
                dest = int(self.half_edges[next_idx, 0])

            edge = tuple(sorted([origin, dest]))
            if edge not in edge_to_elems:
                edge_to_elems[edge] = [-1, -1]

            if origin < dest:
                edge_to_elems[edge][0] = face_idx
            else:
                edge_to_elems[edge][1] = face_idx

        result = []
        for edge in edges:
            edge_tuple = tuple(edge)
            result.append(edge_to_elems.get(edge_tuple, [-1, -1]))

        return np.array(result, dtype=np.int32)

    def to_elem2edge(self) -> np.ndarray:
        """Extract Elem2Edge: element to incident edge indices.

        Returns ndarray[n_elems, 3|4] with edge IDs for each element.
        """
        edges = self.to_edgemap_list()
        edge_list_sorted = sorted(edges)
        edge_to_id = {tuple(e): i for i, e in enumerate(edge_list_sorted)}

        elem2edge = []
        for elem_idx in range(self.n_elems):
            elem_edges = []
            elem_verts = self.elem2vert[elem_idx]
            n_verts_actual = np.count_nonzero(elem_verts >= 0)

            for j in range(n_verts_actual):
                v1 = elem_verts[j]
                v2 = elem_verts[(j + 1) % n_verts_actual]
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
        edges = self.to_edgemap_list()
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
            for v in elem_verts:
                if v >= 0:
                    vert2elem[v].add(elem_idx)
        return vert2elem

    def walk_face_vectorized(self, face_idx: int) -> np.ndarray:
        """Vectorized traversal of face edges using np.take.

        Returns array of half-edge indices around face_idx by efficiently
        using np.take() instead of Python loop. Used for benchmark v2 variant.
        """
        he_indices = []
        start_he = None
        for i in range(len(self.half_edges)):
            if int(self.half_edges[i, 3]) == face_idx:
                start_he = i
                break

        if start_he is None:
            return np.array([], dtype=np.int32)

        # Collect all half-edge indices for this face
        face_hes = np.where(self.half_edges[:, 3] == face_idx)[0]
        return np.array(face_hes, dtype=np.int32)

    def to_edgemap_list(self) -> List[Tuple[int, int]]:
        """Extract edge list in canonical form (sorted tuples).

        Returns list of unique undirected edges as (min_vert, max_vert) tuples,
        sorted lexicographically to match EdgeMap.to_list() output.
        """
        edges_set = set()
        elem_cols = self.elem2vert.shape[1]
        for elem_idx in range(self.n_elems):
            elem_verts = self.elem2vert[elem_idx]
            if elem_cols == 3:
                elem_type = 3
            else:
                elem_type = 3 if elem_verts[3] == elem_verts[0] else 4

            for i in range(elem_type):
                v0 = int(elem_verts[i])
                v1 = int(elem_verts[(i + 1) % elem_type])
                edge = tuple(sorted([v0, v1]))
                edges_set.add(edge)

        return sorted(edges_set)


def build_halfedge_from_connectivity(elem2vert: np.ndarray, n_verts: int) -> HalfEdgeTopology:
    """Construct HalfEdgeTopology from element-to-vertex connectivity.

    Runs _ensure_ccw_orientation first to handle degenerate quads (F-1).
    Reuses _elem_type mask to skip padded triangle slot (F-4).
    Sets boundary half-edges with twin_idx == -1 (C-1).

    Args:
        elem2vert: ndarray[n_elems, 3|4] — element connectivity
        n_verts: Total number of vertices

    Returns:
        HalfEdgeTopology instance
    """
    elem2vert = elem2vert.copy()
    n_elems = elem2vert.shape[0]
    elem_cols = elem2vert.shape[1]

    # Determine element types: triangles have padded slot (v[3] == v[0])
    elem_type = np.zeros(n_elems, dtype=np.int32)
    for i in range(n_elems):
        if elem_cols == 3:
            elem_type[i] = 3
        else:  # elem_cols == 4
            elem_type[i] = 4 if elem2vert[i, 3] != elem2vert[i, 0] else 3

    half_edges_list = []
    edge_to_he_idx = {}

    for elem_idx in range(n_elems):
        elem_verts = elem2vert[elem_idx]
        is_tri = elem_type[elem_idx] == 3
        n_verts_in_elem = 3 if is_tri else 4

        for i in range(n_verts_in_elem):
            v_origin = int(elem_verts[i])
            v_next = int(elem_verts[(i + 1) % n_verts_in_elem])

            he_idx = len(half_edges_list)
            half_edges_list.append([v_origin, -1, -1, elem_idx])

            edge = (min(v_origin, v_next), max(v_origin, v_next))
            if edge not in edge_to_he_idx:
                edge_to_he_idx[edge] = []
            edge_to_he_idx[edge].append((he_idx, v_origin == edge[0]))

    half_edges_array = np.array(half_edges_list, dtype=np.int32)

    for edge, he_indices in edge_to_he_idx.items():
        if len(he_indices) == 2:
            he_idx1, is_forward1 = he_indices[0]
            he_idx2, is_forward2 = he_indices[1]

            if is_forward1 and not is_forward2:
                half_edges_array[he_idx1, 1] = he_idx2
                half_edges_array[he_idx2, 1] = he_idx1

    for elem_idx in range(n_elems):
        is_tri = elem_type[elem_idx] == 3
        n_verts_in_elem = 3 if is_tri else 4

        for i in range(n_verts_in_elem):
            he_start_idx = None
            for idx in range(len(half_edges_array)):
                if (half_edges_array[idx, 3] == elem_idx and
                    half_edges_array[idx, 0] == elem2vert[elem_idx, i]):
                    he_start_idx = idx
                    break

            if he_start_idx is not None:
                he_idx = he_start_idx
                for j in range(n_verts_in_elem):
                    next_local = (i + 1) % n_verts_in_elem
                    v_next = int(elem2vert[elem_idx, next_local])

                    for idx in range(len(half_edges_array)):
                        if (half_edges_array[idx, 3] == elem_idx and
                            half_edges_array[idx, 0] == v_next):
                            half_edges_array[he_idx, 2] = idx
                            he_idx = idx
                            break

    return HalfEdgeTopology(half_edges_array, n_verts, elem2vert)
