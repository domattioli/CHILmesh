"""Quad-edge (4-connected edge) topology backend for CHILmesh.

This module implements a quad-edge topology as an alternative to EdgeMap and half-edge.
Quad-edge stores directed edges (2 per undirected pair) as integer-indexed pointers in a
contiguous array. Each edge maintains four neighbors: origin, twin, next_clockwise, next_counter-clockwise.

Schema (QuadEdgeTopology.quad_edges):
    ndarray[n_edges, 4] with columns:
    - [0] origin: int — vertex ID at the start of this directed edge
    - [1] twin_idx: int — index of the opposite directed edge (-1 for boundary)
    - [2] next_cw: int — index of next edge clockwise around origin vertex
    - [3] next_ccw: int — index of next edge counter-clockwise around origin vertex

Directionality:
    DIRECTED (2 edges per undirected pair): Unlike Wikipedia's undirected 1-per-pair variant,
    this implementation uses 2 directed edges per edge pair. This matches half-edge (Phase 007)
    and is optimal for 2D meshes where each edge has at most 2 adjacent faces.

Sentinel conventions:
    - Boundary edges: twin_idx == -1 (no opposite face; matches Edge2Elem convention)
    - All pointers are 0-indexed; -1 indicates no valid neighbor (boundary only)

Padding behavior (mixed-element support):
    - Triangles are stored in 4-column arrays with vertex 0 in slot 3 (_elem_type mask)
    - Quad-edge construction respects elem_type and creates only 3 edges for triangles, 4 for quads
    - No spurious edges are created for the padding

Algorithm (O(n) complexity):
    Phase 1: Create all directed edges per element, record edge map for twin pairing
    Phase 2: Pair twins via hash-table lookup of reverse edges; boundary edges get -1
    Phase 3: Assign next_cw and next_ccw via face-walk logic per element

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
    """Container for all quad-edges and lookup indices.

    Equivalent in role to EdgeMap and HalfEdgeTopology; emits the same dict/ndarray adjacency
    outputs consumed by the rest of CHILmesh.
    """

    def __init__(self, quad_edges: np.ndarray, n_verts: int, elem2vert: np.ndarray):
        """Initialize QuadEdgeTopology.

        Args:
            quad_edges: ndarray[n_edges, 4] with (origin, twin_idx, next_cw, next_ccw)
            n_verts: Total number of vertices
            elem2vert: Original element-to-vertex connectivity (for reference)
        """
        self.quad_edges = quad_edges
        self.n_verts = n_verts
        self.elem2vert = elem2vert
        self.n_elems = elem2vert.shape[0]

    def to_edge2vert(self) -> np.ndarray:
        """Extract Edge2Vert: sorted undirected edges.

        Returns unique undirected edges in canonical form (min_vert, max_vert),
        sorted lexicographically to match EdgeMap output.

        Built from element connectivity to avoid edge traversal complexity.
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

        # Build element ownership map from quad-edges
        for i in range(len(self.quad_edges)):
            origin = int(self.quad_edges[i, 0])
            twin_idx = int(self.quad_edges[i, 1])

            # Infer face_idx by finding which element owns this edge
            # We scan elem2vert to locate the element containing this directed edge
            face_idx = None
            for elem_idx in range(self.n_elems):
                elem_verts = self.elem2vert[elem_idx]
                n_verts_actual = 3 if (elem_verts.shape[0] > 3 and elem_verts[3] == elem_verts[0]) else (4 if elem_verts.shape[0] == 4 else 3)

                for j in range(n_verts_actual):
                    if int(elem_verts[j]) == origin and int(elem_verts[(j + 1) % n_verts_actual]) == (int(self.quad_edges[i, 0]) if int(self.quad_edges[i, 0]) == origin else int(self.quad_edges[i, 0])):
                        face_idx = elem_idx
                        break
                if face_idx is not None:
                    break

            if face_idx is None:
                continue

            if twin_idx >= 0:
                dest = int(self.quad_edges[twin_idx, 0])
            else:
                next_cw = int(self.quad_edges[i, 2])
                if next_cw >= 0:
                    dest = int(self.quad_edges[next_cw, 0])
                else:
                    continue

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

            # Correct detection of actual element type
            if elem_verts.shape[0] > 3:
                n_verts_actual = 3 if elem_verts[3] == elem_verts[0] else 4

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


def build_quadegg_from_connectivity(elem2vert: np.ndarray, n_verts: int) -> QuadEdgeTopology:
    """Construct QuadEdgeTopology from element-to-vertex connectivity.

    O(n) algorithm:
    1. Create directed edge per element edge in iteration order
    2. Hash-table lookup to pair twins
    3. Single pass to assign next_cw and next_ccw within each element's face

    Args:
        elem2vert: ndarray[n_elems, 3|4] — element connectivity (CCW orientation assumed)
        n_verts: Total number of vertices

    Returns:
        QuadEdgeTopology instance

    Algorithm Complexity: O(n) — each edge processed 3 times across phases.

    Invariants:
    - All edges have valid origin and twin_idx (or -1 for boundary)
    - next_cw and next_ccw form a directed cycle around each vertex
    - Boundary edges have twin_idx == -1 and no opposite face neighbor
    - Mixed-element padding (triangles) creates only 3 edges per element
    """
    elem2vert = elem2vert.copy()
    n_elems = elem2vert.shape[0]
    elem_cols = elem2vert.shape[1]

    # Determine element types
    elem_type = np.zeros(n_elems, dtype=np.int32)
    for i in range(n_elems):
        if elem_cols == 3:
            elem_type[i] = 3
        else:
            elem_type[i] = 4 if elem2vert[i, 3] != elem2vert[i, 0] else 3

    # Phase 1: Create all directed edges in element order
    quad_edges_list = []
    qe_map = {}  # (origin, dest) -> qe_idx
    qe_idx_counter = 0
    elem_to_first_qe = {}  # elem_idx -> first QE index for this element

    for elem_idx in range(n_elems):
        elem_verts = elem2vert[elem_idx]
        n_verts_in_elem = int(elem_type[elem_idx])

        elem_to_first_qe[elem_idx] = qe_idx_counter

        for i in range(n_verts_in_elem):
            v_origin = int(elem_verts[i])
            v_dest = int(elem_verts[(i + 1) % n_verts_in_elem])

            qe = [v_origin, -1, -1, -1]  # origin, twin_idx, next_cw, next_ccw
            quad_edges_list.append(qe)

            # Record QE for twin pairing
            directed_edge = (v_origin, v_dest)
            if directed_edge not in qe_map:
                qe_map[directed_edge] = []
            qe_map[directed_edge].append(qe_idx_counter)

            qe_idx_counter += 1

    quad_edges_array = np.array(quad_edges_list, dtype=np.int32)

    # Phase 2: Pair twins using reverse-edge lookup
    for (v_origin, v_dest), qe_indices in qe_map.items():
        reverse_edge = (v_dest, v_origin)
        if reverse_edge in qe_map:
            reverse_indices = qe_map[reverse_edge]
            # Pair: QE(origin->dest) with QE(dest->origin)
            if len(qe_indices) == 1 and len(reverse_indices) == 1:
                qe_idx = qe_indices[0]
                rev_idx = reverse_indices[0]
                quad_edges_array[qe_idx, 1] = rev_idx
                quad_edges_array[rev_idx, 1] = qe_idx
        # Else: boundary edge; twin_idx remains -1

    # Phase 3: Assign next_cw and next_ccw
    # For a directed edge (v_origin, v_dest), next_cw and next_ccw are the adjacent
    # edges around v_origin in the same element.
    # In CCW-ordered elements:
    # - next_cw (clockwise) points to the edge that ends at v_origin (previous in CCW order)
    # - next_ccw (counter-clockwise) points to the edge that starts after v_dest (next in CCW order)

    for elem_idx in range(n_elems):
        n_verts_in_elem = int(elem_type[elem_idx])
        first_qe_idx = elem_to_first_qe[elem_idx]

        for i in range(n_verts_in_elem):
            curr_qe_idx = first_qe_idx + i
            # The next edge in CCW order around the face
            next_qe_idx = first_qe_idx + ((i + 1) % n_verts_in_elem)
            # The previous edge in CCW order (which is clockwise direction from current origin)
            prev_qe_idx = first_qe_idx + ((i - 1) % n_verts_in_elem)

            # next_ccw: the next edge in counter-clockwise direction around origin
            # (i.e., the next edge in the face's CCW order)
            quad_edges_array[curr_qe_idx, 3] = next_qe_idx

            # next_cw: the previous edge when viewed from current vertex
            # (i.e., it ends at the current origin)
            quad_edges_array[curr_qe_idx, 2] = prev_qe_idx

    return QuadEdgeTopology(quad_edges_array, n_verts, elem2vert)
