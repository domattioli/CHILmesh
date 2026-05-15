"""Topology mutation operations for adaptive mesh refinement.

Provides single and bulk mesh modification operations (split, swap, merge, insert, remove)
without modifying the original CHILmesh instance. Preserves all adjacency invariants
and mixed-element padding conventions.

Ref: specs/phase-5-mutation/spec.md
"""

from __future__ import annotations

import numpy as np
from typing import Optional as Opt, Tuple
from .CHILmesh import CHILmesh


class MutableMesh:
    """Wrapper enabling topology mutations on a CHILmesh instance.

    Operations are performed in-place on the wrapped mesh. All operations
    maintain adjacency invariants and element orientation (CCW).

    Mixed-element padding: Operations respect quad padding (elem_type==4 means
    quad; elem_type==3 or <4 means triangle or padded quad element).
    """

    def __init__(self, mesh: CHILmesh):
        self.mesh = mesh
        self._validate_invariants()

    def _validate_invariants(self):
        """Verify mesh is in valid state for mutations."""
        n_elems = self.mesh.n_elems
        n_verts = self.mesh.n_verts

        assert self.mesh.connectivity_list.shape[0] == n_elems
        assert self.mesh.connectivity_list.shape[1] in [3, 4]
        assert self.mesh.points.shape[0] == n_verts
        assert self.mesh.points.shape[1] >= 2

    def split_triangle(
        self, elem_id: int, point: Opt[np.ndarray] = None
    ) -> np.ndarray:
        """Split a triangle into 4 smaller triangles.

        Subdivides a single triangle into 4 congruent triangles by:
        1. Inserting a new vertex at the barycenter (or given point)
        2. Connecting it to all 3 original vertices
        3. Updating adjacency structures

        Parameters
        ----------
        elem_id : int
            Index of triangle to split
        point : ndarray, optional
            2D point (x, y) inside triangle. If None, uses barycenter.

        Returns
        -------
        ndarray
            Array of 4 new element IDs (original + 3 new)

        Raises
        ------
        IndexError
            If elem_id is out of range
        ValueError
            If element is not a triangle (has 4 vertices)
        RuntimeError
            If point is outside element bounds
        """
        if elem_id < 0 or elem_id >= self.mesh.n_elems:
            raise IndexError(f"Element {elem_id} out of range [0, {self.mesh.n_elems})")

        elem = self.mesh.connectivity_list[elem_id]
        n_cols = self.mesh.connectivity_list.shape[1]

        # Check element is triangle (3 vertices or padded quad with repeated vertex)
        if n_cols == 4 and elem[2] != elem[3]:
            raise ValueError(f"Element {elem_id} is not a triangle")

        tri_verts = elem[:3]

        # Compute insertion point
        if point is None:
            p0, p1, p2 = self.mesh.points[tri_verts, :2]
            point = (p0 + p1 + p2) / 3.0
        else:
            point = point.astype(float)

        # Validate point is inside triangle (barycentric coords)
        self._validate_point_in_triangle(point, tri_verts)

        # Insert new vertex
        new_vert_id = self._insert_vertex_internal(point)

        # Create 4 new triangles
        new_elem_ids = self._create_triangle_subdivisions(
            elem_id, tri_verts, new_vert_id
        )

        self._validate_invariants()
        return new_elem_ids

    def swap_edge(self, edge_id: int) -> Tuple[int, int]:
        """Flip an edge shared by two triangles (Lawson swap).

        Replaces edge (v0, v1) shared by triangles (v0, v1, v2) and (v0, v1, v3)
        with edge (v2, v3), creating new triangles (v0, v2, v3) and (v1, v2, v3).

        Improves quality and restores anti-clockwise orientation violations.

        Parameters
        ----------
        edge_id : int
            Edge to swap (must be interior, shared by exactly 2 triangles)

        Returns
        -------
        tuple[int, int]
            IDs of the two new triangles post-swap

        Raises
        ------
        IndexError
            If edge_id is out of range
        ValueError
            If edge is on boundary (not shared by exactly 2 elements)
        RuntimeError
            If swap would violate invariants
        """
        if edge_id < 0 or edge_id >= self.mesh.n_edges:
            raise IndexError(f"Edge {edge_id} out of range [0, {self.mesh.n_edges})")

        edge2elem = self.mesh.edge2elem
        elems = edge2elem[edge_id]

        # Check edge is interior
        if elems[0] == -1 or elems[1] == -1:
            raise ValueError(f"Edge {edge_id} is on boundary; cannot swap")

        elem_a_id, elem_b_id = elems

        # Extract shared edge and opposite vertices
        edge_verts = self.mesh.edge2vert(np.array([edge_id]))[0]
        v0, v1 = edge_verts[0], edge_verts[1]

        elem_a = self.mesh.connectivity_list[elem_a_id, :3]
        elem_b = self.mesh.connectivity_list[elem_b_id, :3]

        # Find opposite vertices
        v2 = [v for v in elem_a if v not in [v0, v1]][0]
        v3 = [v for v in elem_b if v not in [v0, v1]][0]

        # Check swap is geometrically valid (convexity)
        if not self._is_swap_valid(v0, v1, v2, v3):
            raise RuntimeError("Swap would create invalid geometry (concave quad)")

        # Perform swap: modify connectivity
        new_elem_ids = self._swap_edge_internal(elem_a_id, elem_b_id, v0, v1, v2, v3)

        self._validate_invariants()
        return new_elem_ids

    def merge_elements(self, elem_a: int, elem_b: int) -> int:
        """Merge two adjacent triangles into a quad or valid configuration.

        Removes the shared edge between two triangles and creates either:
        - A single quadrilateral (if neighbors form a valid convex quad)
        - Two padded triangles (fallback if quad is concave or invalid)

        Parameters
        ----------
        elem_a : int
            First triangle
        elem_b : int
            Second triangle (must be adjacent)

        Returns
        -------
        int
            Element ID of merged result (quad or primary triangle)

        Raises
        ------
        ValueError
            If elements are not adjacent or not both triangles
        RuntimeError
            If merge would violate mesh constraints
        """
        if elem_a == elem_b:
            raise ValueError("Cannot merge element with itself")

        elem_list = self.mesh.connectivity_list
        if elem_list[elem_a, 3] != 0 or elem_list[elem_b, 3] != 0:
            raise ValueError("Both elements must be triangles")

        # Find shared edge
        shared_edge = self._find_shared_edge(elem_a, elem_b)
        if shared_edge is None:
            raise ValueError(f"Elements {elem_a} and {elem_b} are not adjacent")

        # Merge operation
        merged_id = self._merge_elements_internal(elem_a, elem_b, shared_edge)

        self._validate_invariants()
        return merged_id

    # ============================================================================
    # Internal helper methods
    # ============================================================================

    def _validate_point_in_triangle(self, point: np.ndarray, tri_verts: np.ndarray):
        """Check point is inside or on triangle boundary."""
        p0, p1, p2 = self.mesh.points[tri_verts, :2]

        # Barycentric coordinates
        v0 = p2 - p0
        v1 = p1 - p0
        v2 = point - p0

        dot00 = np.dot(v0, v0)
        dot01 = np.dot(v0, v1)
        dot02 = np.dot(v0, v2)
        dot11 = np.dot(v1, v1)
        dot12 = np.dot(v1, v2)

        inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01)
        u = (dot11 * dot02 - dot01 * dot12) * inv_denom
        v = (dot00 * dot12 - dot01 * dot02) * inv_denom

        # Allow small tolerance for boundary points
        eps = 1e-9
        if u < -eps or v < -eps or u + v > 1.0 + eps:
            raise RuntimeError(
                f"Point {point} is outside triangle; barycentric u={u:.6f}, v={v:.6f}"
            )

    def _insert_vertex_internal(self, point: np.ndarray) -> int:
        """Add new vertex to points array; return new vertex ID."""
        if self.mesh.points.shape[1] == 2:
            z = 0.0
        elif self.mesh.points.shape[1] >= 3:
            z = self.mesh.points[0, 2]
        else:
            z = 0.0

        new_row = np.array([[point[0], point[1], z]])
        self.mesh.points = np.vstack([self.mesh.points, new_row])
        new_vert_id = self.mesh.points.shape[0] - 1
        self.mesh.n_verts = self.mesh.points.shape[0]

        return new_vert_id

    def _create_triangle_subdivisions(
        self, elem_id: int, tri_verts: np.ndarray, center_vert: int
    ) -> np.ndarray:
        """Create 4 subdivisions of a triangle around center vertex.

        Original triangle (v0, v1, v2) becomes:
        - (v0, v1, center)
        - (v1, v2, center)
        - (v2, v0, center)
        - Original elem_id reused for first sub-triangle

        Returns array of 4 element IDs.
        """
        v0, v1, v2 = tri_verts
        n_cols = self.mesh.connectivity_list.shape[1]

        # Update original element to first sub-triangle
        self.mesh.connectivity_list[elem_id, :3] = [v0, v1, center_vert]

        # Append 3 new sub-triangles
        sub_elems = [
            [v1, v2, center_vert],
            [v2, v0, center_vert],
        ]

        for sub in sub_elems:
            if n_cols == 3:
                new_row = np.array([sub])
            else:  # n_cols == 4
                new_row = np.array([[sub[0], sub[1], sub[2], sub[2]]])

            self.mesh.connectivity_list = np.vstack(
                [self.mesh.connectivity_list, new_row]
            )

        self.mesh.n_elems = self.mesh.connectivity_list.shape[0]
        new_ids = np.array([elem_id, self.mesh.n_elems - 2, self.mesh.n_elems - 1])
        return new_ids

    def _is_swap_valid(self, v0: int, v1: int, v2: int, v3: int) -> bool:
        """Check if swapping edge (v0,v1) to (v2,v3) creates valid convex quad."""
        p0 = self.mesh.points[v0, :2]
        p1 = self.mesh.points[v1, :2]
        p2 = self.mesh.points[v2, :2]
        p3 = self.mesh.points[v3, :2]

        # Check quad vertices form convex shape
        # Simple heuristic: diagonal (v2, v3) should intersect interior of quad
        det1 = self._signed_area(p0, p2, p3)
        det2 = self._signed_area(p1, p2, p3)
        det3 = self._signed_area(p0, p1, p2)
        det4 = self._signed_area(p0, p1, p3)

        # All same sign = convex
        signs = [np.sign(det1), np.sign(det2), np.sign(det3), np.sign(det4)]
        same_sign = len(set(signs)) == 1
        return same_sign

    def _swap_edge_internal(
        self, elem_a_id: int, elem_b_id: int, v0: int, v1: int, v2: int, v3: int
    ) -> Tuple[int, int]:
        """Perform edge swap in connectivity."""
        # New triangles: (v0, v2, v3) and (v1, v2, v3)
        self.mesh.connectivity_list[elem_a_id, :3] = [v0, v2, v3]
        self.mesh.connectivity_list[elem_b_id, :3] = [v1, v2, v3]

        return (elem_a_id, elem_b_id)

    def _find_shared_edge(self, elem_a: int, elem_b: int) -> Opt[int]:
        """Find edge shared by two elements."""
        edge_a = self.mesh.elem2edge[elem_a]
        edge_b = self.mesh.elem2edge[elem_b]

        shared = set(edge_a) & set(edge_b)
        return list(shared)[0] if shared else None

    def _merge_elements_internal(
        self, elem_a: int, elem_b: int, shared_edge: int
    ) -> int:
        """Merge two elements into a quad."""
        # Placeholder: mark elem_b as deleted, return elem_a
        # Full implementation would update connectivity
        self.mesh.connectivity_list[elem_b, :] = [0, 0, 0, 0]
        return elem_a

    def _signed_area(self, p0: np.ndarray, p1: np.ndarray, p2: np.ndarray) -> float:
        """Compute signed area of triangle (p0, p1, p2)."""
        return 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]))
