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
from .mesh_topology import quad_from_tri_pair


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
        if n_cols == 4 and elem[3] != elem[0]:
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

        # Rebuild adjacencies and spatial indices (necessary after adding vertices/elements)
        self.mesh._build_adjacencies()
        self.mesh._build_spatial_indices()
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

        edge2elem_all = self.mesh.edge2elem()  # Call method
        elems = edge2elem_all[edge_id]

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

        # Incremental O(1) adjacency patch; KD-trees deferred (marked dirty).
        self._patch_swap_adjacency(edge_id, int(v0), int(v1), int(v2), int(v3),
                                   int(elem_a_id), int(elem_b_id))
        self.mesh._spatial_dirty = True
        self._validate_invariants()
        return new_elem_ids

    def merge_elements(self, elem_a: int, elem_b: int) -> int:
        """Merge two adjacent triangles into a single quad.

        The shared edge becomes the quad's interior diagonal and is dropped.
        ``elem_a`` is rewritten as the resulting quad (CCW); ``elem_b`` is
        marked deleted in place so element IDs stay stable for the caller.

        Parameters
        ----------
        elem_a : int
            First triangle (becomes the merged quad).
        elem_b : int
            Second triangle; must share exactly one edge with ``elem_a``.

        Returns
        -------
        int
            Element ID of the merged quad (always ``elem_a``).

        Raises
        ------
        IndexError
            If either element ID is out of range.
        ValueError
            If merging an element with itself, if the elements are not
            adjacent, or if either element is not a triangle.
        """
        if elem_a == elem_b:
            raise ValueError("Cannot merge element with itself")

        n_elems = self.mesh.n_elems
        for eid in (elem_a, elem_b):
            if eid < 0 or eid >= n_elems:
                raise IndexError(f"Element {eid} out of range [0, {n_elems})")

        if not (self._is_triangle(elem_a) and self._is_triangle(elem_b)):
            raise ValueError("Both elements must be triangles")

        # Find shared edge
        shared_edge = self._find_shared_edge(elem_a, elem_b)
        if shared_edge is None:
            raise ValueError(f"Elements {elem_a} and {elem_b} are not adjacent")

        # Merge operation
        merged_id = self._merge_elements_internal(elem_a, elem_b, shared_edge)

        # Rebuild adjacencies and spatial indices (necessary after topology change)
        self.mesh._build_adjacencies()
        self.mesh._build_spatial_indices()
        self._validate_invariants()
        return merged_id

    def remove_vertex(self, vert_id: int) -> None:
        """Remove an interior vertex and re-triangulate the resulting cavity.

        Coarsening / degenerate-vertex removal. The incident elements (which
        must be triangles) are deleted and the surrounding ring is fan
        re-triangulated without the removed vertex.

        Tombstone convention: deleted elements are marked ``-1`` (consistent
        with ``merge_elements``); the removed vertex is left orphaned in
        ``points`` (``n_verts`` unchanged, no renumbering) so existing vertex
        IDs stay stable for callers and incremental adjacency (see #162).

        Parameters
        ----------
        vert_id : int
            Index of the interior vertex to remove.

        Raises
        ------
        IndexError
            If ``vert_id`` is out of range.
        ValueError
            If ``vert_id`` is a boundary vertex (would change the domain
            boundary), is already orphaned, or has a non-triangular incident
            element (not supported in v1).
        RuntimeError
            If the cavity cannot be fan re-triangulated without producing a
            non-positive-area element.
        """
        if vert_id < 0 or vert_id >= self.mesh.n_verts:
            raise IndexError(f"Vertex {vert_id} out of range [0, {self.mesh.n_verts})")

        if vert_id in {int(v) for v in self.mesh.boundary_node_indices()}:
            raise ValueError(f"Vertex {vert_id} is on the boundary; cannot remove in v1")

        incident = sorted(int(e) for e in self.mesh.get_vertex_elements(vert_id))
        if not incident:
            raise ValueError(f"Vertex {vert_id} is orphaned (no incident elements)")

        for eid in incident:
            if not self._is_triangle(eid):
                raise ValueError(
                    f"Vertex {vert_id} has a non-triangular incident element {eid}; "
                    "only triangular cavities are supported in v1"
                )

        # Each incident triangle, walked from vert_id in its stored CCW order,
        # contributes the directed edge opposite vert_id. Those directed edges
        # form a single CCW cycle (the cavity ring).
        ring_edges = []
        for eid in incident:
            tri = [int(v) for v in self.mesh.connectivity_list[eid, :3]]
            i = tri.index(vert_id)
            ring_edges.append((tri[(i + 1) % 3], tri[(i + 2) % 3]))

        ring = self._chain_ring(ring_edges)
        if ring is None or len(ring) < 3:
            raise RuntimeError(
                f"Cavity around vertex {vert_id} is not a simple ring; cannot re-triangulate"
            )

        # Fan-triangulate the ring. The cavity is star-shaped from the removed
        # vertex but not necessarily from an arbitrary ring vertex, so try each
        # ring vertex as the fan apex and keep the first that yields only
        # positive-area triangles.
        new_tris = self._fan_triangulate(ring)
        if new_tris is None:
            raise RuntimeError(
                f"Cavity around vertex {vert_id} is not fan-triangulable "
                "without inversion (non-convex); not supported in v1"
            )

        n_cols = self.mesh.connectivity_list.shape[1]
        for eid in incident:
            self.mesh.connectivity_list[eid] = [-1] * n_cols

        rows = [
            [a, b, c, c] if n_cols == 4 else [a, b, c]
            for (a, b, c) in new_tris
        ]
        if rows:
            self.mesh.connectivity_list = np.vstack(
                [self.mesh.connectivity_list, np.array(rows, dtype=self.mesh.connectivity_list.dtype)]
            )
        self.mesh.n_elems = self.mesh.connectivity_list.shape[0]

        self.mesh._build_adjacencies()
        self.mesh._build_spatial_indices()
        self._validate_invariants()

    def collapse_edge(self, edge_id: int) -> int:
        """Collapse an interior edge by merging its endpoints into one vertex.

        The two triangles incident to the edge degenerate to zero area and are
        deleted; every element referencing the removed endpoint is rewritten to
        reference the survivor. Standard coarsening / short-edge elimination.

        Survivor: ``v0`` (the first endpoint of the edge). Tombstone ``-1`` for
        the two collapsed elements; the removed endpoint is left orphaned in
        ``points`` (``n_verts`` unchanged), consistent with ``remove_vertex``.

        Parameters
        ----------
        edge_id : int
            Interior edge to collapse.

        Returns
        -------
        int
            Surviving vertex ID (``v0``).

        Raises
        ------
        IndexError
            If ``edge_id`` is out of range.
        ValueError
            If the edge is on the boundary (not supported in v1).
        RuntimeError
            If the collapse would invert any incident element (signed-area sign
            flip); no mutation is applied in that case.
        """
        if edge_id < 0 or edge_id >= self.mesh.n_edges:
            raise IndexError(f"Edge {edge_id} out of range [0, {self.mesh.n_edges})")

        elems = [int(e) for e in self.mesh.edge2elem(np.array([edge_id]))[0]]
        if -1 in elems:
            raise ValueError(f"Edge {edge_id} is on the boundary; cannot collapse in v1")

        edge_verts = self.mesh.edge2vert(np.array([edge_id]))[0]
        survivor, removed = int(edge_verts[0]), int(edge_verts[1])
        collapsing = {int(e) for e in elems}

        conn = self.mesh.connectivity_list
        n_cols = conn.shape[1]

        # Elements (excluding the two collapsing ones) that reference the removed
        # endpoint. Validate none invert under the substitution BEFORE mutating.
        affected = []
        for eid in range(self.mesh.n_elems):
            if eid in collapsing:
                continue
            row = conn[eid]
            if int(row[0]) < 0:  # tombstoned
                continue
            if removed in (int(v) for v in row):
                affected.append(eid)

        for eid in affected:
            before = self._element_signed_area(conn[eid])
            after_row = np.array(
                [survivor if int(v) == removed else int(v) for v in conn[eid]]
            )
            after = self._element_signed_area(after_row)
            if (before > 0 and after <= 1e-12) or (before < 0 and after >= -1e-12):
                raise RuntimeError(
                    f"collapse_edge({edge_id}) would invert element {eid}; aborted"
                )

        # Commit: rewrite survivors, tombstone the two collapsed elements.
        for eid in affected:
            for j in range(n_cols):
                if int(conn[eid, j]) == removed:
                    conn[eid, j] = survivor
        for eid in collapsing:
            conn[eid] = [-1] * n_cols

        self.mesh._build_adjacencies()
        self.mesh._build_spatial_indices()
        self._validate_invariants()
        return survivor

    def move_boundary_node(self, vert_id: int, new_xy: np.ndarray) -> None:
        """Move a boundary vertex to new coordinates with inversion guard.

        Pure coordinate update — topology (connectivity) is unchanged. The move
        is validated against all incident elements before applying; if any
        element would become non-positive-area the move is aborted.

        Parameters
        ----------
        vert_id : int
            Boundary vertex to move.
        new_xy : array-like, shape (2,)
            New (x, y) position.

        Raises
        ------
        IndexError
            If ``vert_id`` is out of range.
        ValueError
            If ``vert_id`` is not a boundary vertex.
        RuntimeError
            If the move would invert any incident element; no mutation applied.
        """
        if vert_id < 0 or vert_id >= self.mesh.n_verts:
            raise IndexError(f"Vertex {vert_id} out of range [0, {self.mesh.n_verts})")

        boundary = {int(v) for v in self.mesh.boundary_node_indices()}
        if vert_id not in boundary:
            raise ValueError(f"Vertex {vert_id} is not on the boundary")

        new_xy = np.asarray(new_xy, dtype=float).ravel()[:2]
        incident = [int(e) for e in self.mesh.get_vertex_elements(vert_id)]

        for eid in incident:
            row = self.mesh.connectivity_list[eid]
            if int(row[0]) < 0:
                continue
            seen: list = []
            pts: list = []
            for v in row:
                vi = int(v)
                if vi < 0 or vi in seen:
                    continue
                seen.append(vi)
                pts.append(new_xy if vi == vert_id else self.mesh.points[vi, :2])
            if len(pts) < 3:
                continue
            area = self._polygon_signed_area(np.array(pts))
            if area <= 1e-12:
                raise RuntimeError(
                    f"move_boundary_node({vert_id}) would invert element {eid}; aborted"
                )

        self.mesh.points[vert_id, :2] = new_xy
        self.mesh._build_spatial_indices()
        self._validate_invariants()

    def split_triangles(self, elem_ids: np.ndarray) -> np.ndarray:
        """Split multiple triangles atomically; rebuild adjacencies once.

        Equivalent to N sequential ``split_triangle`` calls but with a single
        adjacency + spatial-index rebuild at the end. Atomic: any per-element
        precondition failure rolls back to the pre-call state.

        Parameters
        ----------
        elem_ids : array-like of int
            Element IDs to split (must all be triangles).

        Returns
        -------
        ndarray
            Concatenated array of new element IDs from all splits.

        Raises
        ------
        IndexError
            If any element ID is out of range.
        ValueError
            If any element is not a triangle.
        RuntimeError
            If a provided split point is outside its element; full rollback.
        """
        elem_ids = np.asarray(elem_ids, dtype=int).ravel()
        snap = self._snapshot()
        new_ids_all: list = []
        try:
            for eid in elem_ids:
                if eid < 0 or eid >= self.mesh.n_elems:
                    raise IndexError(f"Element {eid} out of range [0, {self.mesh.n_elems})")
                elem = self.mesh.connectivity_list[eid]
                n_cols = self.mesh.connectivity_list.shape[1]
                if n_cols == 4 and elem[3] != elem[0]:
                    raise ValueError(f"Element {eid} is not a triangle")
                tri_verts = elem[:3]
                p0, p1, p2 = self.mesh.points[tri_verts, :2]
                point = (p0 + p1 + p2) / 3.0
                self._validate_point_in_triangle(point, tri_verts)
                new_vert_id = self._insert_vertex_internal(point)
                new_ids = self._create_triangle_subdivisions(eid, tri_verts, new_vert_id)
                new_ids_all.extend(new_ids.tolist())
        except Exception:
            self._restore(snap)
            raise
        self.mesh._build_adjacencies()
        self.mesh._build_spatial_indices()
        self._validate_invariants()
        return np.array(new_ids_all, dtype=int)

    def smooth_topology(self, metric_threshold: float = 0.0, max_passes: int = 100) -> int:
        """Improve mesh quality by swapping interior edges that raise min angle.

        Scans all interior edges each pass; swaps any edge whose flip raises
        the minimum angle of the two adjacent triangles above
        ``metric_threshold``. Repeats until no improvement or ``max_passes``
        reached (monotone, always terminates).

        Parameters
        ----------
        metric_threshold : float
            Minimum angle improvement (radians) required to trigger a swap.
            Default 0.0 accepts any improvement.
        max_passes : int
            Maximum number of full-edge-scan passes.

        Returns
        -------
        int
            Total number of edges swapped.
        """
        n_swapped = 0
        for _ in range(max_passes):
            swapped_this_pass = 0
            edge2elem = self.mesh.edge2elem()
            for edge_id in range(self.mesh.n_edges):
                ea, eb = int(edge2elem[edge_id][0]), int(edge2elem[edge_id][1])
                if ea == -1 or eb == -1:
                    continue
                if not (self._is_triangle(ea) and self._is_triangle(eb)):
                    continue
                if not self._swap_improves_quality(edge_id, ea, eb, metric_threshold):
                    continue
                try:
                    self.swap_edge(edge_id)
                    n_swapped += 1
                    swapped_this_pass += 1
                    break  # adjacency rebuilt; edge IDs changed — restart this pass
                except (ValueError, RuntimeError):
                    pass
            if swapped_this_pass == 0:
                break
        if n_swapped > 0:
            self.mesh._build_spatial_indices()
            self.mesh._spatial_dirty = False
        return n_swapped

    def reskeletonize_local(
        self, elem_ids: np.ndarray, radius: int = 2
    ) -> None:
        """Re-skeletonize only the neighborhood affected by mutated elements.

        Finds the earliest layer touched by ``elem_ids``, backs up ``radius``
        layers, replays consumption of all earlier layers into fresh working
        arrays, then re-peels from that start point — updating only
        ``layers[start_layer:]`` and ``n_layers``.

        Falls back to a full ``_skeletonize()`` when the affected elements are
        in the first ``radius`` layers or when no layer data exists yet.

        Parameters
        ----------
        elem_ids : array-like of int
            Element IDs changed by the most recent mutation.
        radius : int
            How many layers *before* the first affected layer to start the
            partial re-peel (provides a safety margin). Default 2.
        """
        elem_ids = np.asarray(elem_ids, dtype=int).ravel()
        layers = self.mesh.layers

        if not layers.get('OE'):
            self.mesh._skeletonize()
            return

        # Find the first layer that contains any of the changed elements.
        elem_set = set(int(e) for e in elem_ids)
        affected_layer = None
        for iL in range(self.mesh.n_layers):
            oe_set = set(int(e) for e in layers['OE'][iL])
            ie_set = set(int(e) for e in layers['IE'][iL])
            if elem_set & (oe_set | ie_set):
                affected_layer = iL
                break

        if affected_layer is None:
            # Changed elements not in any existing layer — full rebuild.
            self.mesh._skeletonize()
            return

        start_layer = max(0, affected_layer - radius)
        if start_layer == 0:
            self.mesh._skeletonize()
            return

        # Replay consumption of layers 0..(start_layer-1) to restore the
        # working arrays at that checkpoint.  Element/vertex IDs are stable
        # across _build_adjacencies, so this replay is valid even when edge
        # IDs have been reassigned.
        edge2vert_work = self.mesh.adjacencies['Edge2Vert'].copy()
        edge2elem_work = self.mesh.adjacencies['Edge2Elem'].copy()

        for iL in range(start_layer):
            oe = layers['OE'][iL]
            ov = layers['OV'][iL]
            ie = layers['IE'][iL]
            if len(oe) > 0:
                edge2elem_work[np.isin(edge2elem_work, oe)] = -1
            if len(ov) > 0:
                edge2vert_work[np.isin(edge2vert_work, ov)] = -1
            if len(ie) > 0:
                edge2elem_work[np.isin(edge2elem_work, ie)] = -1

        # Build new layers dict: keep 0..(start_layer-1), re-peel the rest.
        new_layers: dict = {
            'OE': list(layers['OE'][:start_layer]),
            'IE': list(layers['IE'][:start_layer]),
            'OV': list(layers['OV'][:start_layer]),
            'IV': list(layers['IV'][:start_layer]),
            'bEdgeIDs': list(layers['bEdgeIDs'][:start_layer]),
        }

        iL = start_layer
        while np.any(edge2elem_work >= 0):
            active_count = np.sum(edge2elem_work >= 0, axis=1)
            iLbEdgeIDs = np.where(active_count == 1)[0]
            if len(iLbEdgeIDs) == 0:
                break

            ov_raw = edge2vert_work[iLbEdgeIDs].ravel()
            ov = np.unique(ov_raw[ov_raw >= 0]).astype(int)
            new_layers['OV'].append(ov)
            new_layers['bEdgeIDs'].append(iLbEdgeIDs)

            oe_raw = edge2elem_work[iLbEdgeIDs].ravel()
            oe = np.unique(oe_raw[oe_raw >= 0]).astype(int)
            new_layers['OE'].append(oe)
            if len(oe) > 0:
                edge2elem_work[np.isin(edge2elem_work, oe)] = -1

            ov_edge_mask = np.any(np.isin(edge2vert_work, ov), axis=1)
            ov_edge_indices = np.where(ov_edge_mask)[0]
            if len(ov_edge_indices) > 0:
                ie_raw = edge2elem_work[ov_edge_indices].ravel()
                ie = np.unique(ie_raw[ie_raw >= 0]).astype(int)
            else:
                ie = np.empty(0, dtype=int)
            new_layers['IE'].append(ie)

            if len(ov) > 0:
                edge2vert_work[np.isin(edge2vert_work, ov)] = -1
            if len(ie) > 0:
                edge2elem_work[np.isin(edge2elem_work, ie)] = -1

            if len(oe) > 0 or len(ie) > 0:
                layer_elems = np.concatenate((oe, ie))
                lv = self.mesh.connectivity_list[layer_elems].ravel()
                lv = lv[lv >= 0]
                iv = np.setdiff1d(np.unique(lv), ov).astype(int)
            else:
                iv = np.empty(0, dtype=int)
            new_layers['IV'].append(iv)

            iL += 1

        self.mesh.layers = new_layers
        self.mesh.n_layers = len(new_layers['OE'])

    def skeletonize_diff(self, prev_layers: dict) -> dict:
        """Run full re-skeletonization and return which elements changed layer.

        Parameters
        ----------
        prev_layers : dict
            Snapshot captured before a mutation (via ``_snapshot_layers``).

        Returns
        -------
        dict
            ``{elem_id: {'old': (kind, layer_idx) | None,
                         'new': (kind, layer_idx) | None}}``
            for every element whose layer assignment changed.  ``kind`` is
            ``'OE'`` or ``'IE'``; ``None`` means the element was not in any
            layer (tombstoned or newly created).
        """
        prev_map: dict = {}
        for iL, (oe, ie) in enumerate(zip(prev_layers.get('OE', []),
                                          prev_layers.get('IE', []))):
            for e in oe:
                prev_map[int(e)] = ('OE', iL)
            for e in ie:
                prev_map[int(e)] = ('IE', iL)

        self.mesh._skeletonize()

        new_map: dict = {}
        for iL, (oe, ie) in enumerate(zip(self.mesh.layers['OE'],
                                          self.mesh.layers['IE'])):
            for e in oe:
                new_map[int(e)] = ('OE', iL)
            for e in ie:
                new_map[int(e)] = ('IE', iL)

        changed: dict = {}
        for eid in set(prev_map) | set(new_map):
            old = prev_map.get(eid)
            new = new_map.get(eid)
            if old != new:
                changed[eid] = {'old': old, 'new': new}
        return changed

    def _snapshot_layers(self) -> dict:
        """Snapshot current layer arrays for ``skeletonize_diff``."""
        return {
            key: [arr.copy() for arr in val]
            for key, val in self.mesh.layers.items()
        }

    def _snapshot(self) -> dict:
        """Capture mesh state for atomic rollback."""
        return {
            'points': self.mesh.points.copy(),
            'connectivity_list': self.mesh.connectivity_list.copy(),
            'n_verts': self.mesh.n_verts,
            'n_elems': self.mesh.n_elems,
        }

    def _restore(self, snap: dict) -> None:
        """Restore mesh state from snapshot and rebuild structures."""
        self.mesh.points = snap['points']
        self.mesh.connectivity_list = snap['connectivity_list']
        self.mesh.n_verts = snap['n_verts']
        self.mesh.n_elems = snap['n_elems']
        self.mesh._build_adjacencies()
        self.mesh._build_spatial_indices()

    def _swap_improves_quality(
        self, edge_id: int, elem_a: int, elem_b: int, threshold: float
    ) -> bool:
        """True if swapping this edge raises min angle of both adjacent triangles."""
        edge_verts = self.mesh.edge2vert(np.array([edge_id]))[0]
        v0, v1 = int(edge_verts[0]), int(edge_verts[1])
        row_a = self.mesh.connectivity_list[elem_a, :3]
        row_b = self.mesh.connectivity_list[elem_b, :3]
        try:
            v2 = next(int(v) for v in row_a if int(v) not in (v0, v1))
            v3 = next(int(v) for v in row_b if int(v) not in (v0, v1))
        except StopIteration:
            return False
        before = min(
            self._min_angle_from_verts(v0, v1, v2),
            self._min_angle_from_verts(v0, v1, v3),
        )
        after = min(
            self._min_angle_from_verts(v0, v2, v3),
            self._min_angle_from_verts(v1, v2, v3),
        )
        return after > before + threshold

    def _min_angle_from_verts(self, va: int, vb: int, vc: int) -> float:
        """Minimum interior angle (radians) of the triangle (va, vb, vc)."""
        pts = self.mesh.points[[va, vb, vc], :2]
        angles = []
        for i in range(3):
            a = pts[(i + 1) % 3] - pts[i]
            b = pts[(i + 2) % 3] - pts[i]
            na, nb = np.linalg.norm(a), np.linalg.norm(b)
            if na < 1e-14 or nb < 1e-14:
                return 0.0
            cos_a = np.clip(np.dot(a, b) / (na * nb), -1.0, 1.0)
            angles.append(float(np.arccos(cos_a)))
        return min(angles)

    def _patch_swap_adjacency(
        self,
        old_eid: int,
        v0: int, v1: int, v2: int, v3: int,
        ea: int, eb: int,
    ) -> None:
        """O(1) incremental adjacency patch after swapping edge old_eid (v0-v1)→(v2-v3).

        Reuses ``old_eid`` slot for the new edge (no array resizing). Updates:
        Edge2Vert, EdgeMap, Vert2Edge, Edge2Elem (for 3 affected edges),
        Vert2Elem (4 vertices), Elem2Edge (2 elements).  KD-trees are NOT
        rebuilt here — caller must set ``mesh._spatial_dirty = True``.
        """
        adj = self.mesh.adjacencies
        edge_map = adj['EdgeMap']

        # Reuse old_eid slot for the new diagonal {v2, v3}.
        nv2, nv3 = min(v2, v3), max(v2, v3)
        adj['Edge2Vert'][old_eid] = [nv2, nv3]
        del edge_map._map[(min(v0, v1), max(v0, v1))]
        edge_map._map[(nv2, nv3)] = old_eid

        # Vert2Edge: v0/v1 lose old_eid; v2/v3 gain it.
        adj['Vert2Edge'][v0].discard(old_eid)
        adj['Vert2Edge'][v1].discard(old_eid)
        adj['Vert2Edge'][v2].add(old_eid)
        adj['Vert2Edge'][v3].add(old_eid)

        # Edge2Elem for new diagonal: same two elements.
        adj['Edge2Elem'][old_eid] = [ea, eb]

        # Edge2Elem for {v1,v2}: ea→eb (ea no longer owns this edge).
        e12 = edge_map.find_edge(v1, v2)
        if e12 is not None:
            r = adj['Edge2Elem'][e12]
            if int(r[0]) == ea:
                adj['Edge2Elem'][e12, 0] = eb
            elif int(r[1]) == ea:
                adj['Edge2Elem'][e12, 1] = eb

        # Edge2Elem for {v0,v3}: eb→ea (eb no longer owns this edge).
        e03 = edge_map.find_edge(v0, v3)
        if e03 is not None:
            r = adj['Edge2Elem'][e03]
            if int(r[0]) == eb:
                adj['Edge2Elem'][e03, 0] = ea
            elif int(r[1]) == eb:
                adj['Edge2Elem'][e03, 1] = ea

        # Vert2Elem: v0 leaves eb; v1 leaves ea; v2 joins eb; v3 joins ea.
        adj['Vert2Elem'][v0].discard(eb)
        adj['Vert2Elem'][v1].discard(ea)
        adj['Vert2Elem'][v2].add(eb)
        adj['Vert2Elem'][v3].add(ea)

        # Elem2Edge: recompute from actual new connectivity rows.
        conn = self.mesh.connectivity_list
        n_cols = conn.shape[1]
        for elem_id in (ea, eb):
            row = conn[elem_id]
            n_v = 3 if n_cols == 3 or int(row[2]) == int(row[3 if n_cols == 4 else 2]) else 4
            for i in range(n_v):
                va, vb = int(row[i]), int(row[(i + 1) % n_v])
                if va < 0 or vb < 0 or va == vb:
                    continue
                eid = edge_map.find_edge(va, vb)
                if eid is not None:
                    adj['Elem2Edge'][elem_id, i] = eid

    def _is_triangle(self, elem_id: int) -> bool:
        """True if the element is a triangle (3-col row, or 4-col with a repeated vertex)."""
        if self.mesh.connectivity_list.shape[1] == 3:
            return True
        row = self.mesh.connectivity_list[elem_id]
        return len({int(v) for v in row}) <= 3

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

        denom = dot00 * dot11 - dot01 * dot01

        # Check for degenerate triangle (zero or near-zero area)
        if abs(denom) < 1e-12:
            raise RuntimeError(
                f"Point {point} cannot be validated against degenerate triangle (zero area)"
            )

        inv_denom = 1.0 / denom
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
        """Check if swapping edge (v0,v1) to (v2,v3) creates valid geometry.

        Simple check: new diagonal should have non-zero length and not create
        zero-area triangles. More sophisticated checks (convexity, angles) deferred.
        """
        p0 = self.mesh.points[v0, :2]
        p1 = self.mesh.points[v1, :2]
        p2 = self.mesh.points[v2, :2]
        p3 = self.mesh.points[v3, :2]

        # Check new diagonal has non-zero length
        diag_len = np.linalg.norm(p3 - p2)
        if diag_len < 1e-10:
            return False

        # Check neither new triangle has zero area
        area1 = abs(self._signed_area(p0, p2, p3))
        area2 = abs(self._signed_area(p1, p2, p3))

        return area1 > 1e-10 and area2 > 1e-10

    def _swap_edge_internal(
        self, elem_a_id: int, elem_b_id: int, v0: int, v1: int, v2: int, v3: int
    ) -> Tuple[int, int]:
        """Perform edge swap in connectivity, enforcing CCW winding."""
        pts = self.mesh.points
        new_a = [v0, v2, v3]
        if self._signed_area(pts[v0, :2], pts[v2, :2], pts[v3, :2]) < 0:
            new_a = [v0, v3, v2]
        new_b = [v1, v2, v3]
        if self._signed_area(pts[v1, :2], pts[v2, :2], pts[v3, :2]) < 0:
            new_b = [v1, v3, v2]
        self.mesh.connectivity_list[elem_a_id, :3] = new_a
        self.mesh.connectivity_list[elem_b_id, :3] = new_b
        return (elem_a_id, elem_b_id)

    def _find_shared_edge(self, elem_a: int, elem_b: int) -> Opt[int]:
        """Find edge shared by two elements."""
        elem2edge_all = self.mesh.elem2edge()  # Call method
        edge_a = elem2edge_all[elem_a]
        edge_b = elem2edge_all[elem_b]

        shared = set(edge_a) & set(edge_b)
        return list(shared)[0] if shared else None

    def _merge_elements_internal(
        self, elem_a: int, elem_b: int, shared_edge: int
    ) -> int:
        """Fuse two adjacent triangles into a quad written into ``elem_a``.

        The two shared-edge vertices become opposite corners' neighbours and
        the shared edge is dropped (it becomes the quad's interior diagonal).
        ``elem_b`` is marked deleted with the negative-vertex sentinel (``-1``),
        which the adjacency builder skips (only ``v >= 0`` entries are kept).
        Using ``-1`` rather
        than ``0`` avoids spuriously binding the deleted row to vertex 0.

        A triangle-only (3-column) table is widened to 4 columns first, padding
        existing triangles as ``[v0, v1, v2, v2]``; the caller's subsequent
        ``_build_adjacencies()`` normalises padding and flips ``mesh.type``.
        """
        quad = quad_from_tri_pair(self.mesh.points, self.mesh.connectivity_list[elem_a, :3], self.mesh.connectivity_list[elem_b, :3])

        if self.mesh.connectivity_list.shape[1] == 3:
            self.mesh.connectivity_list = np.column_stack([self.mesh.connectivity_list, self.mesh.connectivity_list[:, 2]])

        self.mesh.connectivity_list[elem_a, :] = quad
        self.mesh.connectivity_list[elem_b, :] = [-1, -1, -1, -1]

        return elem_a

    @staticmethod
    def _polygon_signed_area(pts: np.ndarray) -> float:
        """Shoelace signed area of an (n, 2) polygon; positive when CCW."""
        x = pts[:, 0]
        y = pts[:, 1]
        return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))

    @staticmethod
    def _chain_ring(edges):
        """Chain directed edges into an ordered vertex cycle.

        ``edges`` is a list of ``(a, b)`` directed edges that should form a
        single simple cycle. Returns the ordered list of vertices around the
        cycle, or ``None`` if the edges branch or do not close into one cycle.
        """
        succ = {}
        for a, b in edges:
            if a in succ:  # branching: not a simple cycle
                return None
            succ[a] = b
        if len(succ) != len(edges):
            return None
        start = edges[0][0]
        ring = [start]
        cur = succ[start]
        while cur != start:
            if cur not in succ or len(ring) > len(edges):
                return None
            ring.append(cur)
            cur = succ[cur]
        return ring

    def _fan_triangulate(self, ring):
        """Fan-triangulate a vertex ring, trying each apex.

        Returns a list of ``(a, b, c)`` triangles whose signed areas are all
        positive, or ``None`` if no ring vertex works as a fan apex (cavity
        non-convex).
        """
        n = len(ring)
        for s in range(n):
            order = ring[s:] + ring[:s]
            apex = order[0]
            tris = []
            ok = True
            for i in range(1, n - 1):
                a, b, c = apex, order[i], order[i + 1]
                area = self._signed_area(
                    self.mesh.points[a, :2],
                    self.mesh.points[b, :2],
                    self.mesh.points[c, :2],
                )
                if area <= 1e-12:
                    ok = False
                    break
                tris.append((a, b, c))
            if ok:
                return tris
        return None

    def _element_signed_area(self, row: np.ndarray) -> float:
        """Signed area of an element row over its distinct (non-negative) vertices."""
        verts = []
        for v in row:
            vi = int(v)
            if vi < 0 or vi in verts:
                continue
            verts.append(vi)
        if len(verts) < 3:
            return 0.0
        return self._polygon_signed_area(self.mesh.points[verts, :2])

    def insert_vertex(self, point: np.ndarray) -> int:
        """Insert vertex into mesh; auto-triangulates containing element.

        MVP implementation: Uses Bowyer-Watson cavity approach for interior points.
        - Finds containing element using spatial indexing
        - Collects cavity (affected elements)
        - Inserts new vertex
        - Re-triangulates cavity with fan triangulation

        Parameters
        ----------
        point : ndarray
            2D coordinate (x, y) to insert

        Returns
        -------
        int
            New vertex ID

        Raises
        ------
        ValueError
            If point is outside mesh
        """
        point = np.asarray(point)

        # Find containing element
        elem_id = self.mesh.find_element(point)
        if elem_id == -1:
            raise ValueError(f"Point {point} is outside mesh")

        # Insert new vertex
        new_vert_id = self._insert_vertex_internal(point)

        # Find cavity: all elements sharing vertices with containing element
        containing_elem = self.mesh.connectivity_list[elem_id]
        n_cols = self.mesh.connectivity_list.shape[1]
        elem_verts = containing_elem[:3] if (n_cols == 3 or containing_elem[2] != containing_elem[3]) else containing_elem[:3]

        cavity_elems = set([elem_id])
        for v in elem_verts:
            incident = self.mesh.get_vertex_elements(v)
            cavity_elems.update(incident)

        cavity_elems = sorted(list(cavity_elems))

        # Find boundary edges of cavity (edges with exactly 1 incident cavity elem)
        boundary_edges = []
        for elem_id_cav in cavity_elems:
            elem = self.mesh.connectivity_list[elem_id_cav]
            tri_verts = elem[:3] if n_cols == 3 or elem[2] != elem[3] else elem[:3]
            edges = [(tri_verts[0], tri_verts[1]),
                    (tri_verts[1], tri_verts[2]),
                    (tri_verts[2], tri_verts[0])]

            for e in edges:
                edge_key = tuple(sorted(e))
                count = sum(1 for e2 in cavity_elems
                           if edge_key in self._get_edge_set(e2, n_cols))
                if count == 1:
                    boundary_edges.append(e)

        # Remove duplicate edges (keep unique)
        boundary_edges = list(set([tuple(sorted(e)) for e in boundary_edges]))
        boundary_edges = [e for e in boundary_edges if e[0] != e[1]]

        if not boundary_edges:
            boundary_edges = [(elem_verts[0], elem_verts[1]),
                             (elem_verts[1], elem_verts[2]),
                             (elem_verts[2], elem_verts[0])]

        # Delete cavity elements (mark as zeros)
        for elem_id_del in cavity_elems:
            if n_cols == 3:
                self.mesh.connectivity_list[elem_id_del] = [0, 0, 0]
            else:
                self.mesh.connectivity_list[elem_id_del] = [0, 0, 0, 0]

        # Re-triangulate: create triangles from new vertex to boundary edges
        # Boundary edges must be consistently oriented (walked in order around cavity)
        # For now, sort edges by angle from new vertex to ensure consistent ordering
        p_new = self.mesh.points[new_vert_id, :2]

        # Sort boundary edges by angle from new vertex
        edge_angles = []
        for e in boundary_edges:
            v1, v2 = e
            p1 = self.mesh.points[v1, :2]
            p2 = self.mesh.points[v2, :2]
            # Use midpoint angle
            mid = (p1 + p2) / 2
            angle = np.arctan2(mid[1] - p_new[1], mid[0] - p_new[0])
            edge_angles.append((angle, v1, v2))

        edge_angles.sort()  # Sort by angle

        # Create triangles in angle order with CCW check
        for angle, v1, v2 in edge_angles:
            p1 = self.mesh.points[v1, :2]
            p2 = self.mesh.points[v2, :2]

            # Check (v1, v2, new_vert) triangle area - should be positive
            area = self._signed_area(p1, p2, p_new)
            if area < 0:
                # Swap to ensure positive area
                v1, v2 = v2, v1

            new_elem = np.array([[v1, v2, new_vert_id]])
            self.mesh.connectivity_list = np.vstack([self.mesh.connectivity_list, new_elem])

        self.mesh.n_elems = self.mesh.connectivity_list.shape[0]

        # Rebuild adjacencies and spatial indices
        self.mesh._build_adjacencies()
        self.mesh._build_spatial_indices()
        self._validate_invariants()

        return new_vert_id

    def _get_edge_set(self, elem_id: int, n_cols: int) -> set:
        """Get edge set for an element."""
        first_vert = self.mesh.connectivity_list[elem_id, 0]
        if elem_id >= self.mesh.n_elems or first_vert == 0 or first_vert < 0:
            return set()

        elem = self.mesh.connectivity_list[elem_id]
        tri_verts = elem[:3] if n_cols == 3 or elem[2] != elem[3] else elem[:3]
        edges = [(tri_verts[0], tri_verts[1]),
                (tri_verts[1], tri_verts[2]),
                (tri_verts[2], tri_verts[0])]
        return set([tuple(sorted(e)) for e in edges if e[0] != e[1]])

    def _signed_area(self, p0: np.ndarray, p1: np.ndarray, p2: np.ndarray) -> float:
        """Compute signed area of triangle (p0, p1, p2)."""
        return 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]))
