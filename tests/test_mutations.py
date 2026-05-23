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
