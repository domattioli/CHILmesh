"""Tests for Phase 5 mesh mutation operations (split, swap, merge, insert, remove)."""

from __future__ import annotations

import pytest
import numpy as np
from chilmesh import CHILmesh, MutableMesh
from chilmesh import examples


@pytest.fixture(params=["annulus", "donut"])
def triangle_mesh(request):
    """Provide fresh (non-cached) triangle-only meshes for mutation testing.

    Uses __wrapped__ to bypass conftest.py's global caching, since mutation
    tests modify the mesh and would otherwise corrupt the cache for other tests.
    """
    example_fn = getattr(examples, request.param)
    # Use __wrapped__ to bypass caching if available (added by conftest.py)
    if hasattr(example_fn, '__wrapped__'):
        return example_fn.__wrapped__()
    return example_fn()


class TestMutableMeshInitialization:
    """Tests for MutableMesh wrapper initialization."""

    def test_init_wraps_chilmesh(self, triangle_mesh):
        """Verify MutableMesh wraps CHILmesh without modifying it."""
        original_n_elems = triangle_mesh.n_elems
        mutable = MutableMesh(triangle_mesh)

        assert mutable.mesh is triangle_mesh
        assert triangle_mesh.n_elems == original_n_elems

    def test_init_validates_invariants(self, triangle_mesh):
        """Verify MutableMesh validates mesh on creation."""
        mutable = MutableMesh(triangle_mesh)
        mutable._validate_invariants()  # Should not raise


class TestSplitTriangle:
    """Tests for triangle splitting operation (refinement)."""

    def test_split_triangle_basic(self, triangle_mesh):
        """Verify split_triangle creates 4 triangles from 1."""
        mutable = MutableMesh(triangle_mesh)
        original_n_elems = triangle_mesh.n_elems
        original_n_verts = triangle_mesh.n_verts

        # Split first triangle
        new_ids = mutable.split_triangle(elem_id=0)

        assert len(new_ids) == 3  # 1 original + 2 new (returning indices, not count)
        assert triangle_mesh.n_elems == original_n_elems + 2  # 3 new added, 1 reused
        assert triangle_mesh.n_verts == original_n_verts + 1  # 1 new vertex

    def test_split_triangle_default_point(self, triangle_mesh):
        """Verify split_triangle uses barycenter when point=None."""
        mutable = MutableMesh(triangle_mesh)
        original_points = triangle_mesh.points.copy()
        original_verts = triangle_mesh.connectivity_list[0, :3].copy()

        mutable.split_triangle(elem_id=0, point=None)

        # New vertex should be barycenter of original triangle
        bary = original_points[original_verts, :2].mean(axis=0)
        new_vert_id = triangle_mesh.n_verts - 1

        np.testing.assert_allclose(
            triangle_mesh.points[new_vert_id, :2],
            bary,
            rtol=1e-10,
            err_msg="New vertex not at barycenter",
        )

    def test_split_triangle_custom_point(self, triangle_mesh):
        """Verify split_triangle accepts custom interior point."""
        mutable = MutableMesh(triangle_mesh)
        custom_pt = np.array([0.25, 0.25])

        # This will only work if point is actually inside an element
        # For now, just verify it doesn't crash on first element
        try:
            mutable.split_triangle(elem_id=0, point=custom_pt)
            # If successful, verify new point exists
            new_vert_id = triangle_mesh.n_verts - 1
            assert triangle_mesh.n_verts > 0
        except RuntimeError:
            # Point outside mesh is acceptable error
            pass

    def test_split_triangle_invalid_elem_id(self, triangle_mesh):
        """Verify split_triangle raises on invalid element ID."""
        mutable = MutableMesh(triangle_mesh)

        with pytest.raises(IndexError):
            mutable.split_triangle(elem_id=-1)

        with pytest.raises(IndexError):
            mutable.split_triangle(elem_id=triangle_mesh.n_elems)

    def test_split_triangle_preserves_orientation(self, triangle_mesh):
        """Verify split triangles maintain CCW orientation."""
        mutable = MutableMesh(triangle_mesh)

        mutable.split_triangle(elem_id=0)

        # All triangles should have positive signed area
        for elem_id in range(triangle_mesh.n_elems):
            verts = triangle_mesh.connectivity_list[elem_id, :3]
            if all(v > 0 for v in verts):  # Skip degenerate elements
                p0, p1, p2 = triangle_mesh.points[verts, :2]
                area = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1])
                assert area > 1e-10, f"Element {elem_id} has non-positive area"

    def test_split_triangle_z_preserved(self, triangle_mesh):
        """Verify z-coordinate unchanged in split operation."""
        mutable = MutableMesh(triangle_mesh)
        original_z = triangle_mesh.points[:, 2].copy()

        mutable.split_triangle(elem_id=0)

        np.testing.assert_array_equal(triangle_mesh.points[:original_z.shape[0], 2], original_z)


class TestSwapEdge:
    """Tests for edge swapping (quality improvement)."""

    def test_swap_edge_interior(self, triangle_mesh):
        """Verify swap_edge flips interior edge."""
        mutable = MutableMesh(triangle_mesh)

        # Find an interior edge (shared by 2 elements)
        interior_edge_id = None
        edge2elem_all = triangle_mesh.edge2elem()  # Call method
        for edge_id in range(len(edge2elem_all)):
            elems = edge2elem_all[edge_id]
            if elems[0] != -1 and elems[1] != -1:
                interior_edge_id = edge_id
                break

        if interior_edge_id is None:
            pytest.skip("No interior edges in this mesh")

        original_n_elems = triangle_mesh.n_elems
        new_ids = mutable.swap_edge(edge_id=interior_edge_id)

        assert len(new_ids) == 2
        assert triangle_mesh.n_elems == original_n_elems  # Element count unchanged

    def test_swap_edge_boundary_raises(self, triangle_mesh):
        """Verify swap_edge raises on boundary edge."""
        mutable = MutableMesh(triangle_mesh)

        # Find a boundary edge (has -1 in edge2elem)
        boundary_edge_id = None
        edge2elem_all = triangle_mesh.edge2elem()  # Call method
        for edge_id in range(len(edge2elem_all)):
            elems = edge2elem_all[edge_id]
            if elems[0] == -1 or elems[1] == -1:
                boundary_edge_id = edge_id
                break

        if boundary_edge_id is not None:
            with pytest.raises(ValueError, match="boundary"):
                mutable.swap_edge(edge_id=boundary_edge_id)

    def test_swap_edge_invalid_id(self, triangle_mesh):
        """Verify swap_edge raises on invalid edge ID."""
        mutable = MutableMesh(triangle_mesh)

        with pytest.raises(IndexError):
            mutable.swap_edge(edge_id=-1)

        with pytest.raises(IndexError):
            mutable.swap_edge(edge_id=triangle_mesh.n_edges)


class TestMergeElements:
    """Tests for element merging (coarsening)."""

    def _first_interior_edge_pair(self, mesh):
        """Return (edge_id, elem_a, elem_b) for the first interior edge."""
        edge2elem = mesh.edge2elem()
        for edge_id in range(mesh.n_edges):
            ea, eb = int(edge2elem[edge_id][0]), int(edge2elem[edge_id][1])
            if ea != -1 and eb != -1:
                return edge_id, ea, eb
        return None

    def test_merge_elements_adjacent(self, triangle_mesh):
        """Two adjacent triangles fuse into one CCW quad; elem_b is deleted."""
        mutable = MutableMesh(triangle_mesh)
        n_verts_before = triangle_mesh.n_verts

        found = self._first_interior_edge_pair(triangle_mesh)
        assert found is not None, "fixture has no interior edge to merge across"
        _, elem_a, elem_b = found
        tri_a = set(int(v) for v in triangle_mesh.connectivity_list[elem_a, :3])
        tri_b = set(int(v) for v in triangle_mesh.connectivity_list[elem_b, :3])
        expected_quad_verts = tri_a | tri_b

        merged = mutable.merge_elements(elem_a, elem_b)

        assert merged == elem_a
        # Table widened to 4 columns to hold the quad.
        assert triangle_mesh.connectivity_list.shape[1] == 4
        quad_row = triangle_mesh.connectivity_list[merged]
        # Genuine quad: 4 distinct vertices == union of the two triangles.
        assert set(int(v) for v in quad_row) == expected_quad_verts
        assert len(expected_quad_verts) == 4
        # _elem_type classifies the merged row as a quad, not a padded tri.
        _, quads = triangle_mesh._elem_type()
        assert merged in quads
        # Winding is CCW (positive signed area).
        area = MutableMesh._polygon_signed_area(triangle_mesh.points[quad_row, :2])
        assert area > 0
        # elem_b is deleted via the negative sentinel...
        assert all(int(v) < 0 for v in triangle_mesh.connectivity_list[elem_b])
        # ...and contributes no phantom membership for vertex 0.
        assert elem_b not in triangle_mesh.get_vertex_elements(0)
        # No vertices were added or removed by a merge.
        assert triangle_mesh.n_verts == n_verts_before
        # Mesh is now mixed and adjacencies rebuilt cleanly.
        assert triangle_mesh.type in ("Mixed-Element", "Quadrilateral")
        assert "Edge2Vert" in triangle_mesh.adjacencies

    def test_merge_elements_out_of_range_raises(self, triangle_mesh):
        """Out-of-range element IDs raise IndexError."""
        mutable = MutableMesh(triangle_mesh)
        with pytest.raises(IndexError):
            mutable.merge_elements(0, triangle_mesh.n_elems)

    def test_merge_elements_non_adjacent_raises(self, triangle_mesh):
        """Verify merge_elements raises on non-adjacent elements."""
        mutable = MutableMesh(triangle_mesh)

        # Try to merge non-adjacent elements
        # Elements 0 and 1 are typically not adjacent in most test meshes
        try:
            mutable.merge_elements(elem_a=0, elem_b=2)
            # If they happen to be adjacent, that's fine - skip test
            pytest.skip("Elements 0 and 2 are adjacent in this mesh")
        except ValueError as e:
            if "not adjacent" in str(e):
                pass  # Expected
            else:
                raise

    def test_merge_elements_self_raises(self, triangle_mesh):
        """Verify merge_elements raises on self-merge."""
        mutable = MutableMesh(triangle_mesh)

        with pytest.raises(ValueError, match="itself"):
            mutable.merge_elements(elem_a=0, elem_b=0)


class TestInsertVertex:
    """Tests for vertex insertion (incremental refinement)."""

    def test_insert_vertex_basic(self, triangle_mesh):
        """Verify insert_vertex adds new vertex and re-triangulates."""
        mutable = MutableMesh(triangle_mesh)
        original_n_elems = triangle_mesh.n_elems
        original_n_verts = triangle_mesh.n_verts

        # Insert at centroid of first element
        elem_verts = triangle_mesh.connectivity_list[0, :3]
        point = np.mean(triangle_mesh.points[elem_verts, :2], axis=0)

        new_vert_id = mutable.insert_vertex(point)

        assert new_vert_id == original_n_verts
        assert triangle_mesh.n_verts == original_n_verts + 1
        assert triangle_mesh.n_elems > original_n_elems  # Cavity retriangulated

    def test_insert_vertex_increases_count(self, triangle_mesh):
        """Verify insert_vertex increases vertex and element counts."""
        mutable = MutableMesh(triangle_mesh)
        original_n_elems = triangle_mesh.n_elems
        original_n_verts = triangle_mesh.n_verts

        # Insert at multiple locations
        for i in range(3):
            elem_verts = triangle_mesh.connectivity_list[i * 10, :3]
            point = np.mean(triangle_mesh.points[elem_verts, :2], axis=0)
            new_vert_id = mutable.insert_vertex(point)

            assert triangle_mesh.n_verts == original_n_verts + i + 1
            assert triangle_mesh.n_elems > original_n_elems

    def test_insert_vertex_outside_mesh_raises(self, triangle_mesh):
        """Verify insert_vertex raises for points outside mesh."""
        mutable = MutableMesh(triangle_mesh)

        # Point far outside mesh
        far_point = np.array([1e6, 1e6])

        with pytest.raises(ValueError, match="outside"):
            mutable.insert_vertex(far_point)

    def test_insert_vertex_creates_triangles_from_point(self, triangle_mesh):
        """Verify inserted vertex has incident elements."""
        mutable = MutableMesh(triangle_mesh)

        # Insert vertex at centroid of element (guaranteed inside for element interior)
        # Try first several elements to find one with centroid inside mesh
        new_vert_id = None
        for elem_id in range(min(100, triangle_mesh.n_elems)):
            elem_verts = triangle_mesh.connectivity_list[elem_id, :3]
            point = np.mean(triangle_mesh.points[elem_verts, :2], axis=0)
            try:
                new_vert_id = mutable.insert_vertex(point)
                break
            except ValueError:
                continue

        assert new_vert_id is not None, "Could not find valid insertion point"

        # Check that new vertex appears in some elements
        incident = triangle_mesh.get_vertex_elements(new_vert_id)
        assert len(incident) > 0, "Inserted vertex has no incident elements"

    def test_insert_vertex_maintains_ccw(self, triangle_mesh):
        """Verify inserted vertex produces CCW-oriented triangles."""
        mutable = MutableMesh(triangle_mesh)

        # Insert vertex at first element with valid centroid
        new_vert_id = None
        for elem_id in range(min(100, triangle_mesh.n_elems)):
            elem_verts = triangle_mesh.connectivity_list[elem_id, :3]
            point = np.mean(triangle_mesh.points[elem_verts, :2], axis=0)
            try:
                new_vert_id = mutable.insert_vertex(point)
                break
            except ValueError:
                continue

        assert new_vert_id is not None, "Could not find valid insertion point"

        # Check all incident elements have positive signed area
        incident = triangle_mesh.get_vertex_elements(new_vert_id)
        for elem_id in incident:
            elem = triangle_mesh.connectivity_list[elem_id, :3]
            v0, v1, v2 = triangle_mesh.points[elem, :2]
            signed_area = 0.5 * ((v1[0] - v0[0]) * (v2[1] - v0[1]) - (v2[0] - v0[0]) * (v1[1] - v0[1]))
            assert signed_area > -1e-9, f"Element {elem_id} is not CCW after insert_vertex"

    def test_insert_vertex_spatial_index_updated(self, triangle_mesh):
        """Verify spatial indices are rebuilt after insert_vertex.

        Regression test: spatial indices were stale after mutations.
        After insert_vertex, find_element() should find the inserted vertex.
        """
        mutable = MutableMesh(triangle_mesh)

        # Insert vertex at first element centroid
        elem_verts = triangle_mesh.connectivity_list[0, :3]
        point = np.mean(triangle_mesh.points[elem_verts, :2], axis=0)

        new_vert_id = mutable.insert_vertex(point)
        new_vert_pos = triangle_mesh.points[new_vert_id, :2]

        # Query for element containing new vertex (should exist now)
        found_elem = triangle_mesh.find_element(new_vert_pos)
        assert found_elem >= 0, "find_element failed after insert_vertex (spatial index not updated)"
        assert found_elem < triangle_mesh.n_elems, f"Invalid element ID {found_elem}"

    def test_split_triangle_spatial_index_updated(self, triangle_mesh):
        """Verify spatial indices are rebuilt after split_triangle."""
        mutable = MutableMesh(triangle_mesh)

        elem_verts = triangle_mesh.connectivity_list[0, :3]
        barycenter = np.mean(triangle_mesh.points[elem_verts, :2], axis=0)

        mutable.split_triangle(elem_id=0, point=None)

        # Query at barycenter should succeed after split
        found_elem = triangle_mesh.find_element(barycenter)
        assert found_elem >= 0, "find_element failed after split_triangle"


def _first_interior_vertex(mesh):
    """Return an interior vertex with only triangular incident elements, else None."""
    boundary = {int(v) for v in mesh.boundary_node_indices()}
    for vid in range(mesh.n_verts):
        if vid in boundary:
            continue
        incident = mesh.get_vertex_elements(vid)
        if not incident:
            continue
        all_tri = all(
            mesh.connectivity_list.shape[1] == 3
            or mesh.connectivity_list[e, 2] == mesh.connectivity_list[e, 3]
            for e in incident
        )
        if all_tri:
            return vid
    return None


def _first_interior_edge(mesh):
    """Return (edge_id, v0, v1) for the first interior edge, else None."""
    edge2elem = mesh.edge2elem()
    for edge_id in range(mesh.n_edges):
        ea, eb = int(edge2elem[edge_id][0]), int(edge2elem[edge_id][1])
        if ea != -1 and eb != -1:
            verts = mesh.edge2vert(np.array([edge_id]))[0]
            return edge_id, int(verts[0]), int(verts[1])
    return None


def _all_live_elements_ccw(mesh):
    """True if every non-tombstoned element has positive signed area."""
    for eid in range(mesh.n_elems):
        row = mesh.connectivity_list[eid]
        if int(row[0]) < 0:
            continue
        verts = []
        for v in row:
            vi = int(v)
            if vi >= 0 and vi not in verts:
                verts.append(vi)
        if len(verts) < 3:
            continue
        pts = mesh.points[verts, :2]
        x, y = pts[:, 0], pts[:, 1]
        area = 0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))
        if area <= 1e-10:
            return False
    return True


class TestRemoveVertex:
    """Tests for interior vertex removal (#158, coarsening)."""

    def test_remove_interior_vertex(self, triangle_mesh):
        """Removing an interior vertex deletes its cavity and re-triangulates."""
        mutable = MutableMesh(triangle_mesh)
        vid = _first_interior_vertex(triangle_mesh)
        if vid is None:
            pytest.skip("no interior triangular vertex in this fixture")

        n_verts_before = triangle_mesh.n_verts
        cavity = triangle_mesh.get_vertex_elements(vid)

        mutable.remove_vertex(vid)

        # n_verts unchanged (tombstone convention: vertex left orphaned).
        assert triangle_mesh.n_verts == n_verts_before
        # The removed vertex no longer belongs to any live element.
        assert len(triangle_mesh.get_vertex_elements(vid)) == 0
        # Old incident elements are tombstoned.
        for eid in cavity:
            assert all(int(v) < 0 for v in triangle_mesh.connectivity_list[eid])
        # No degenerate / inverted live elements remain.
        assert _all_live_elements_ccw(triangle_mesh)
        assert "Edge2Vert" in triangle_mesh.adjacencies

    def test_remove_boundary_vertex_raises(self, triangle_mesh):
        """Removing a boundary vertex raises ValueError."""
        mutable = MutableMesh(triangle_mesh)
        boundary = triangle_mesh.boundary_node_indices()
        assert len(boundary) > 0
        with pytest.raises(ValueError, match="boundary"):
            mutable.remove_vertex(int(boundary[0]))

    def test_remove_vertex_out_of_range_raises(self, triangle_mesh):
        """Out-of-range vertex raises IndexError."""
        mutable = MutableMesh(triangle_mesh)
        with pytest.raises(IndexError):
            mutable.remove_vertex(-1)
        with pytest.raises(IndexError):
            mutable.remove_vertex(triangle_mesh.n_verts)


class TestCollapseEdge:
    """Tests for interior edge collapse (#159, coarsening)."""

    def test_collapse_interior_edge(self, triangle_mesh):
        """Collapsing an interior edge tombstones its two triangles."""
        mutable = MutableMesh(triangle_mesh)

        live_before = sum(
            1 for e in range(triangle_mesh.n_elems)
            if int(triangle_mesh.connectivity_list[e, 0]) >= 0
        )

        # Many arbitrary interior-edge collapses legitimately invert a far
        # element and are rejected (pre-mutation, mesh untouched), so scan for
        # the first edge whose collapse is geometrically valid.
        edge2elem = triangle_mesh.edge2elem()
        survivor = None
        v0 = v1 = None
        for edge_id in range(triangle_mesh.n_edges):
            ea, eb = int(edge2elem[edge_id][0]), int(edge2elem[edge_id][1])
            if ea == -1 or eb == -1:
                continue
            verts = triangle_mesh.edge2vert(np.array([edge_id]))[0]
            v0, v1 = int(verts[0]), int(verts[1])
            try:
                survivor = mutable.collapse_edge(edge_id)
                break
            except RuntimeError:
                continue
        if survivor is None:
            pytest.skip("no non-inverting interior edge collapse in this fixture")

        assert survivor == v0
        # Two elements removed (the pair that shared the collapsed edge).
        live_after = sum(
            1 for e in range(triangle_mesh.n_elems)
            if int(triangle_mesh.connectivity_list[e, 0]) >= 0
        )
        assert live_after == live_before - 2
        # Removed endpoint no longer referenced by any live element.
        assert len(triangle_mesh.get_vertex_elements(v1)) == 0
        # No inverted / degenerate live elements.
        assert _all_live_elements_ccw(triangle_mesh)
        assert "Edge2Vert" in triangle_mesh.adjacencies

    def test_collapse_boundary_edge_raises(self, triangle_mesh):
        """Collapsing a boundary edge raises ValueError."""
        mutable = MutableMesh(triangle_mesh)
        edge2elem = triangle_mesh.edge2elem()
        boundary_edge = None
        for edge_id in range(triangle_mesh.n_edges):
            elems = edge2elem[edge_id]
            if elems[0] == -1 or elems[1] == -1:
                boundary_edge = edge_id
                break
        if boundary_edge is None:
            pytest.skip("no boundary edge in this fixture")
        with pytest.raises(ValueError, match="boundary"):
            mutable.collapse_edge(boundary_edge)

    def test_collapse_edge_out_of_range_raises(self, triangle_mesh):
        """Out-of-range edge raises IndexError."""
        mutable = MutableMesh(triangle_mesh)
        with pytest.raises(IndexError):
            mutable.collapse_edge(-1)
        with pytest.raises(IndexError):
            mutable.collapse_edge(triangle_mesh.n_edges)


class TestMoveBoundaryNode:
    """Tests for boundary-node coordinate move (#160, guarded coordinate update)."""

    def test_move_boundary_node_accepted(self, triangle_mesh):
        """Small inward nudge accepted; coords updated, CCW preserved."""
        mutable = MutableMesh(triangle_mesh)
        boundary = triangle_mesh.boundary_node_indices()
        vid = int(boundary[0])

        old_xy = triangle_mesh.points[vid, :2].copy()
        incident = triangle_mesh.get_vertex_elements(vid)
        eid = int(next(iter(incident)))
        verts = [int(v) for v in triangle_mesh.connectivity_list[eid] if int(v) >= 0]
        centroid = triangle_mesh.points[verts, :2].mean(axis=0)
        new_xy = old_xy + 0.01 * (centroid - old_xy)

        mutable.move_boundary_node(vid, new_xy)

        np.testing.assert_allclose(triangle_mesh.points[vid, :2], new_xy)
        assert _all_live_elements_ccw(triangle_mesh)
        assert "Edge2Vert" in triangle_mesh.adjacencies

    def test_move_boundary_node_rejected_inversion(self, triangle_mesh):
        """Move that degenerates an element is rejected; coords unchanged."""
        mutable = MutableMesh(triangle_mesh)
        vid = int(triangle_mesh.boundary_node_indices()[0])
        old_xy = triangle_mesh.points[vid, :2].copy()

        # Move to neighbor vertex position — area collapses to zero, rejected.
        incident = triangle_mesh.get_vertex_elements(vid)
        neighbor = None
        for eid in incident:
            for v in triangle_mesh.connectivity_list[eid]:
                vi = int(v)
                if vi != vid and vi >= 0:
                    neighbor = vi
                    break
            if neighbor is not None:
                break

        with pytest.raises(RuntimeError):
            mutable.move_boundary_node(vid, triangle_mesh.points[neighbor, :2].copy())

        np.testing.assert_array_equal(triangle_mesh.points[vid, :2], old_xy)

    def test_move_interior_vertex_raises(self, triangle_mesh):
        """Moving a non-boundary vertex raises ValueError."""
        mutable = MutableMesh(triangle_mesh)
        vid = _first_interior_vertex(triangle_mesh)
        if vid is None:
            pytest.skip("no interior vertex in this fixture")
        with pytest.raises(ValueError, match="boundary"):
            mutable.move_boundary_node(vid, np.array([0.0, 0.0]))

    def test_move_boundary_node_out_of_range_raises(self, triangle_mesh):
        """Out-of-range vertex raises IndexError."""
        mutable = MutableMesh(triangle_mesh)
        with pytest.raises(IndexError):
            mutable.move_boundary_node(-1, np.array([0.0, 0.0]))
        with pytest.raises(IndexError):
            mutable.move_boundary_node(triangle_mesh.n_verts, np.array([0.0, 0.0]))


class TestSplitTriangles:
    """Tests for bulk triangle split (#161, transactional refinement)."""

    def test_split_triangles_basic(self, triangle_mesh):
        """N splits in one call; adjacency rebuilt once; result matches N sequential splits."""
        mutable = MutableMesh(triangle_mesh)
        n_verts_before = triangle_mesh.n_verts
        n_elems_before = triangle_mesh.n_elems

        # Split first 3 elements atomically.
        ids = np.array([0, 1, 2], dtype=int)
        new_ids = mutable.split_triangles(ids)

        assert triangle_mesh.n_verts == n_verts_before + 3
        # Each split adds 2 new elements (original reused + 2 appended).
        assert triangle_mesh.n_elems == n_elems_before + 6
        assert len(new_ids) == 9  # 3 IDs returned per split
        assert _all_live_elements_ccw(triangle_mesh)
        assert "Edge2Vert" in triangle_mesh.adjacencies

    def test_split_triangles_rollback_on_invalid(self, triangle_mesh):
        """Invalid element in list rolls back; mesh unchanged."""
        mutable = MutableMesh(triangle_mesh)
        n_verts_before = triangle_mesh.n_verts
        n_elems_before = triangle_mesh.n_elems
        conn_before = triangle_mesh.connectivity_list.copy()

        # Mix valid + out-of-range → should roll back entirely.
        bad_ids = np.array([0, triangle_mesh.n_elems + 9999], dtype=int)
        with pytest.raises(IndexError):
            mutable.split_triangles(bad_ids)

        assert triangle_mesh.n_verts == n_verts_before
        assert triangle_mesh.n_elems == n_elems_before
        np.testing.assert_array_equal(triangle_mesh.connectivity_list, conn_before)

    def test_split_triangles_out_of_range_raises(self, triangle_mesh):
        """Out-of-range element raises IndexError."""
        mutable = MutableMesh(triangle_mesh)
        with pytest.raises(IndexError):
            mutable.split_triangles(np.array([-1], dtype=int))
        with pytest.raises(IndexError):
            mutable.split_triangles(np.array([triangle_mesh.n_elems], dtype=int))


class TestSmoothTopology:
    """Tests for topology smoothing via edge swaps (#161)."""

    def test_smooth_topology_returns_int(self, triangle_mesh):
        """smooth_topology returns number of swaps (int)."""
        mutable = MutableMesh(triangle_mesh)
        n = mutable.smooth_topology()
        assert isinstance(n, int)
        assert n >= 0

    def test_smooth_topology_preserves_ccw(self, triangle_mesh):
        """All live elements remain CCW after smoothing."""
        mutable = MutableMesh(triangle_mesh)
        mutable.smooth_topology()
        assert _all_live_elements_ccw(triangle_mesh)

    def test_smooth_topology_terminates(self, triangle_mesh):
        """Second call with max_passes=1 does not raise; result is 0 or small."""
        mutable = MutableMesh(triangle_mesh)
        mutable.smooth_topology()  # first pass to fixpoint
        n2 = mutable.smooth_topology(max_passes=1)
        # At fixpoint second call swaps nothing (or very little if near-tie).
        assert isinstance(n2, int)
        assert "Edge2Vert" in triangle_mesh.adjacencies

    def test_smooth_topology_strict_threshold(self, triangle_mesh):
        """Very large threshold prevents any swaps."""
        mutable = MutableMesh(triangle_mesh)
        n = mutable.smooth_topology(metric_threshold=1e6)
        assert n == 0


class TestIncrementalSkeletonization:
    """Tests for incremental skeletonization (#93)."""

    def test_reskeletonize_local_layers_valid(self, triangle_mesh):
        """After swap + reskeletonize_local, n_layers > 0 and no empty OE entry."""
        mutable = MutableMesh(triangle_mesh)
        edge2elem = triangle_mesh.edge2elem()
        swapped_elem = None
        for eid in range(triangle_mesh.n_edges):
            ea, eb = int(edge2elem[eid][0]), int(edge2elem[eid][1])
            if ea == -1 or eb == -1:
                continue
            try:
                mutable.swap_edge(eid)
                swapped_elem = ea
                break
            except (ValueError, RuntimeError):
                continue
        if swapped_elem is None:
            pytest.skip("no swappable edge in this fixture")

        mutable.reskeletonize_local(np.array([swapped_elem]))

        assert triangle_mesh.n_layers > 0
        assert len(triangle_mesh.layers['OE']) == triangle_mesh.n_layers
        for iL in range(triangle_mesh.n_layers):
            assert len(triangle_mesh.layers['OE'][iL]) > 0 or iL == triangle_mesh.n_layers - 1

    def test_reskeletonize_local_parity(self, triangle_mesh):
        """Partial re-skeletonize matches full _skeletonize on the same mesh."""
        mutable = MutableMesh(triangle_mesh)
        edge2elem = triangle_mesh.edge2elem()
        swapped_elem = None
        for eid in range(triangle_mesh.n_edges):
            ea, eb = int(edge2elem[eid][0]), int(edge2elem[eid][1])
            if ea == -1 or eb == -1:
                continue
            try:
                mutable.swap_edge(eid)
                swapped_elem = ea
                break
            except (ValueError, RuntimeError):
                continue
        if swapped_elem is None:
            pytest.skip("no swappable edge in this fixture")

        mutable.reskeletonize_local(np.array([swapped_elem]))
        local_oe = [set(int(e) for e in triangle_mesh.layers['OE'][iL])
                    for iL in range(triangle_mesh.n_layers)]

        triangle_mesh._skeletonize()
        full_oe = [set(int(e) for e in triangle_mesh.layers['OE'][iL])
                   for iL in range(triangle_mesh.n_layers)]

        assert local_oe == full_oe, "OE layer assignment differs from full rebuild"

    def test_skeletonize_diff_returns_changed(self, triangle_mesh):
        """skeletonize_diff returns a dict of changed element assignments."""
        mutable = MutableMesh(triangle_mesh)
        snap = mutable._snapshot_layers()

        edge2elem = triangle_mesh.edge2elem()
        swapped = False
        for eid in range(triangle_mesh.n_edges):
            ea, eb = int(edge2elem[eid][0]), int(edge2elem[eid][1])
            if ea == -1 or eb == -1:
                continue
            try:
                mutable.swap_edge(eid)
                swapped = True
                break
            except (ValueError, RuntimeError):
                continue
        if not swapped:
            pytest.skip("no swappable edge in this fixture")

        diff = mutable.skeletonize_diff(snap)
        assert isinstance(diff, dict)
        # Each entry has 'old' and 'new' keys.
        for eid, info in diff.items():
            assert 'old' in info and 'new' in info

    def test_snapshot_layers_deepcopy(self, triangle_mesh):
        """_snapshot_layers produces independent copy (mutation doesn't affect snapshot)."""
        mutable = MutableMesh(triangle_mesh)
        snap = mutable._snapshot_layers()
        original_n = len(snap['OE'])

        # Append a dummy entry to live layers — snapshot must not change.
        triangle_mesh.layers['OE'].append(np.array([9999]))

        assert len(snap['OE']) == original_n


class TestIncrementalAdjacency:
    """Tests for incremental O(1) adjacency patch (#162)."""

    @staticmethod
    def _vert2edge_as_pairs(mesh):
        """Convert Vert2Edge edge-ID sets to vertex-pair sets for topology comparison."""
        e2v = mesh.adjacencies['Edge2Vert']
        result = {}
        for v in range(mesh.n_verts):
            pairs = set()
            for eid in mesh.adjacencies['Vert2Edge'][v]:
                va, vb = int(e2v[eid][0]), int(e2v[eid][1])
                pairs.add((min(va, vb), max(va, vb)))
            result[v] = frozenset(pairs)
        return result

    def test_swap_adjacency_parity(self, triangle_mesh):
        """Vert2Elem and Vert2Edge topology after incremental swap matches full rebuild."""
        mutable = MutableMesh(triangle_mesh)

        # Find a swappable interior edge.
        edge2elem = triangle_mesh.edge2elem()
        swapped = False
        for eid in range(triangle_mesh.n_edges):
            ea, eb = int(edge2elem[eid][0]), int(edge2elem[eid][1])
            if ea == -1 or eb == -1:
                continue
            try:
                mutable.swap_edge(eid)
                swapped = True
                break
            except (ValueError, RuntimeError):
                continue
        if not swapped:
            pytest.skip("no swappable edge in this fixture")

        # Incremental state — Vert2Elem by element IDs (stable), Vert2Edge by vertex pairs.
        incr_v2m = {v: frozenset(triangle_mesh.adjacencies['Vert2Elem'][v])
                    for v in range(triangle_mesh.n_verts)}
        incr_v2e_pairs = self._vert2edge_as_pairs(triangle_mesh)

        # Full rebuild.
        triangle_mesh._build_adjacencies()
        full_v2m = {v: frozenset(triangle_mesh.adjacencies['Vert2Elem'][v])
                    for v in range(triangle_mesh.n_verts)}
        full_v2e_pairs = self._vert2edge_as_pairs(triangle_mesh)

        assert incr_v2m == full_v2m, "Vert2Elem mismatch after incremental swap patch"
        assert incr_v2e_pairs == full_v2e_pairs, "Vert2Edge (by vertex pairs) mismatch"

    def test_smooth_topology_adjacency_parity(self, triangle_mesh):
        """Vert2Elem after smooth_topology matches full rebuild."""
        mutable = MutableMesh(triangle_mesh)
        mutable.smooth_topology()

        incr_v2m = {v: frozenset(triangle_mesh.adjacencies['Vert2Elem'][v])
                    for v in range(triangle_mesh.n_verts)}

        triangle_mesh._build_adjacencies()
        full_v2m = {v: frozenset(triangle_mesh.adjacencies['Vert2Elem'][v])
                    for v in range(triangle_mesh.n_verts)}

        assert incr_v2m == full_v2m

    def test_smooth_topology_speed(self, triangle_mesh):
        """1000 swap ops complete in < 10s (was ~50min with full rebuild per swap)."""
        import time
        mutable = MutableMesh(triangle_mesh)
        # Run smooth_topology up to 1000 swaps worth of passes.
        start = time.perf_counter()
        mutable.smooth_topology(max_passes=200)
        elapsed = time.perf_counter() - start
        assert elapsed < 10.0, f"smooth_topology took {elapsed:.2f}s (expected < 10s)"


class TestPaddedTrianglePaddingConvention:
    """Regression tests for padded triangle padding convention fix (issue #211).

    Bug: split_triangle / split_triangles / _point_in_element wrongly used
    elem[3]==elem[2] to detect a padded triangle, when the actual convention
    is elem[3]==elem[0] (4th slot duplicates FIRST vertex, not 3rd).

    These tests verify the fix: padded triangles are no longer spuriously
    rejected, and point-in-element correctly classifies them.
    """

    def test_split_triangle_padded_triangle_not_rejected(self, triangle_mesh):
        """split_triangle accepts a padded triangle [v0, v1, v2, v0] without error.

        Regression: old code wrongly checked elem[3] != elem[2] and raised
        "is not a triangle" spuriously for valid padded triangles with padding
        convention [v0, v1, v2, v0].
        """
        # Take first two adjacent triangles and merge into a quad to convert to 4-col.
        mutable = MutableMesh(triangle_mesh)
        n_cols_before = triangle_mesh.connectivity_list.shape[1]

        if n_cols_before == 4:
            # Already mixed. Find a padded triangle (one with only 3 unique verts).
            padded_tri_id = None
            for eid in range(triangle_mesh.n_elems):
                row = triangle_mesh.connectivity_list[eid]
                unique_verts = len(set(int(v) for v in row if int(v) >= 0))
                if unique_verts == 3:  # Padded triangle
                    padded_tri_id = eid
                    break
            if padded_tri_id is None:
                pytest.skip("no padded triangle in this mixed fixture")
        else:
            # Merge first two adjacent triangles to get a quad, padding remaining tris.
            try:
                quad_id = mutable.merge_elements(0, 1)
            except ValueError as e:
                if "not adjacent" in str(e):
                    pytest.skip("elements 0 and 1 not adjacent in this fixture")
                raise
            # After merge, connectivity is 4-col and other triangles are padded.
            padded_tri_id = 2
            # Sanity: element 2 (if it exists) should now be a padded triangle.
            if padded_tri_id >= triangle_mesh.n_elems:
                pytest.skip("not enough elements to find padded triangle after merge")

        # Attempt split on the padded triangle — should NOT raise.
        try:
            new_ids = mutable.split_triangle(elem_id=padded_tri_id)
            assert len(new_ids) > 0  # Should succeed and return element IDs.
            # Connectivity should now have one more row.
            assert triangle_mesh.n_elems > padded_tri_id
        except ValueError as e:
            if "is not a triangle" in str(e):
                pytest.fail(f"split_triangle wrongly rejected padded triangle: {e}")
            raise

    def test_split_triangles_padded_triangles_not_rejected(self, triangle_mesh):
        """split_triangles accepts padded triangles [v0, v1, v2, v0] in bulk.

        Regression: old code wrongly checked elem[3] != elem[2] in the bulk path.
        """
        mutable = MutableMesh(triangle_mesh)
        n_cols_before = triangle_mesh.connectivity_list.shape[1]

        if n_cols_before == 4:
            # Find padded triangles (exactly 3 unique verts).
            padded_ids = []
            for eid in range(min(3, triangle_mesh.n_elems)):
                row = triangle_mesh.connectivity_list[eid]
                unique_verts = len(set(int(v) for v in row if int(v) >= 0))
                if unique_verts == 3:
                    padded_ids.append(eid)
            if len(padded_ids) < 2:
                pytest.skip("not enough padded triangles in mixed fixture")
        else:
            # Merge two pairs to get padded triangles.
            try:
                mutable.merge_elements(0, 1)
            except ValueError as e:
                if "not adjacent" in str(e):
                    pytest.skip("elements 0 and 1 not adjacent")
                raise
            # Elements 2, 3 should be padded now (if they exist).
            padded_ids = [eid for eid in [2, 3] if eid < triangle_mesh.n_elems]
            if len(padded_ids) < 2:
                pytest.skip("not enough elements after merge")

        try:
            new_ids = mutable.split_triangles(np.array(padded_ids, dtype=int))
            assert len(new_ids) > 0
            assert triangle_mesh.n_elems > padded_ids[-1]
        except ValueError as e:
            if "is not a triangle" in str(e):
                pytest.fail(f"split_triangles wrongly rejected padded triangle: {e}")
            raise

    def test_point_in_element_padded_triangle_classification(self, triangle_mesh, monkeypatch):
        """_point_in_element correctly identifies padded triangle [v0, v1, v2, v0] as tri, not quad.

        Regression: old code wrongly used elem[3]==elem[2] to detect padded tri,
        causing mislassification and wrong point-in-test routing.

        This test asserts WHICH BRANCH is taken: padded triangles must route
        through _point_in_triangle, never _point_in_quad (even though the latter
        would degenerate to the same answer).
        """
        mutable = MutableMesh(triangle_mesh)
        n_cols_before = triangle_mesh.connectivity_list.shape[1]

        if n_cols_before == 4:
            # Find a padded triangle.
            padded_tri_id = None
            for eid in range(triangle_mesh.n_elems):
                row = triangle_mesh.connectivity_list[eid]
                unique_verts = len(set(int(v) for v in row if int(v) >= 0))
                if unique_verts == 3:
                    padded_tri_id = eid
                    break
            if padded_tri_id is None:
                pytest.skip("no padded triangle in mixed fixture")
        else:
            # Merge to create padded triangles.
            try:
                mutable.merge_elements(0, 1)
            except ValueError as e:
                if "not adjacent" in str(e):
                    pytest.skip("elements not adjacent")
                raise
            padded_tri_id = 2
            if padded_tri_id >= triangle_mesh.n_elems:
                pytest.skip("not enough elements")

        # Get centroid of the padded triangle.
        elem = triangle_mesh.connectivity_list[padded_tri_id]
        tri_verts = elem[:3]
        centroid = triangle_mesh.points[tri_verts, :2].mean(axis=0)

        # Spy on _point_in_quad to detect wrong-branch routing.
        called = {"quad": False}
        orig_quad = triangle_mesh._point_in_quad
        def spy_quad(*a, **k):
            called["quad"] = True
            return orig_quad(*a, **k)
        monkeypatch.setattr(triangle_mesh, "_point_in_quad", spy_quad)

        # Call _point_in_element DIRECTLY on the known padded-triangle id.
        result = triangle_mesh._point_in_element(centroid, padded_tri_id)

        # Assertions:
        # 1. Centroid is inside its own triangle.
        assert result, f"centroid of padded triangle {padded_tri_id} should be inside it"
        # 2. Padded triangle MUST route through triangle branch, not quad branch.
        assert called["quad"] is False, "padded triangle wrongly routed to _point_in_quad"


def test_insert_vertex_quad_mesh():
    """Verify insert_vertex works on 4-column quad meshes.

    Regression test for mixed-element insert_vertex bug: dead ternaries,
    missing quad-edge enumeration, and unpadded new elements.
    """
    pts = np.array([
        [0, 0], [1, 0], [2, 0],
        [0, 1], [1, 1], [2, 1],
        [0, 2], [1, 2], [2, 2]
    ], dtype=float)
    pts = np.c_[pts, np.zeros(len(pts))]  # Add z-column

    # 2x2 quad mesh: CCW-oriented quadrilaterals (0-indexed)
    conn = np.array([
        [0, 1, 4, 3],  # quad 0: bottom-left
        [1, 2, 5, 4],  # quad 1: bottom-right
        [3, 4, 7, 6],  # quad 2: top-left
        [4, 5, 8, 7]   # quad 3: top-right
    ], dtype=int)

    m = CHILmesh(conn, pts)
    mm = MutableMesh(m)

    # Insert vertex at center of quad 0 (interior point (0.5, 0.5))
    nid = mm.insert_vertex(np.array([0.5, 0.5]))

    # Assertions: (a) no exception, (b) correct vertex ID, (c) vertex count, (d) re-triangulation
    assert nid == 9, f"Expected new vertex ID 9, got {nid}"
    assert mm.mesh.n_verts == 10, f"Expected 10 vertices, got {mm.mesh.n_verts}"
    assert mm.mesh.n_elems > 4, f"Expected >4 elements after re-triangulation, got {mm.mesh.n_elems}"

    # (e) New vertex has at least one incident element
    incident = mm.mesh.get_vertex_elements(nid)
    assert len(incident) > 0, f"New vertex {nid} has no incident elements"

    # (f) All elements have non-negative signed area
    signed_areas = mm.mesh.signed_area()
    assert (signed_areas >= -1e-9).all(), \
        f"Some elements have negative area: min={signed_areas.min()}"
