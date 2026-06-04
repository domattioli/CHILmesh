"""Tests for standalone element_quality() function.

Tests that the raw-array quality function matches CHILmesh.element_quality()
on existing fixtures and handles edge cases correctly.
"""
from __future__ import annotations

import numpy as np
import pytest
import chilmesh
from chilmesh import CHILmesh


class TestStandaloneQualityMatchesInstance:
    """Verify standalone function matches instance method on real fixtures."""

    def test_matches_annulus(self):
        """Test on annulus fixture (all triangles)."""
        mesh = CHILmesh.read_from_fort14("src/chilmesh/data/annulus_200pts.fort.14")
        expected = mesh.element_quality(metric="aspect_ratio")

        # Convert connectivity to list of lists for standalone function
        conn = [list(row) for row in mesh.connectivity_list]
        actual = chilmesh.element_quality(mesh.points, conn, metric="aspect_ratio")

        np.testing.assert_allclose(actual, expected, atol=1e-10)

    def test_matches_structured_mesh(self):
        """Test on structured quad mesh."""
        mesh = CHILmesh.read_from_fort14("src/chilmesh/data/structuredMesh1.14")
        expected = mesh.element_quality(metric="aspect_ratio")

        conn = [list(row) for row in mesh.connectivity_list]
        actual = chilmesh.element_quality(mesh.points, conn, metric="aspect_ratio")

        np.testing.assert_allclose(actual, expected, atol=1e-10)

    def test_matches_donut(self):
        """Test on donut fixture (all triangles)."""
        mesh = CHILmesh.read_from_fort14("src/chilmesh/data/donut_domain.fort.14")
        expected = mesh.element_quality(metric="aspect_ratio")

        conn = [list(row) for row in mesh.connectivity_list]
        actual = chilmesh.element_quality(mesh.points, conn, metric="aspect_ratio")

        np.testing.assert_allclose(actual, expected, atol=1e-10)


class TestQualityEdgeCases:
    """Test edge cases and degenerate elements."""

    def test_ideal_equilateral_triangle(self):
        """Equilateral triangle should have quality ~1.0."""
        h = np.sqrt(3) / 2
        verts = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, h]])
        conn = [[0, 1, 2]]

        q = chilmesh.element_quality(verts, conn, metric="aspect_ratio")

        assert len(q) == 1
        assert abs(q[0] - 1.0) < 1e-6

    def test_collinear_points_degenerate(self):
        """Collinear points (zero area) should return quality 0.0."""
        verts = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
        conn = [[0, 1, 2]]

        q = chilmesh.element_quality(verts, conn, metric="aspect_ratio")

        assert q[0] == 0.0

    def test_right_isoceles_triangle(self):
        """Right isoceles triangle should have aspect ratio < 1."""
        verts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        conn = [[0, 1, 2]]

        q = chilmesh.element_quality(verts, conn, metric="aspect_ratio")

        assert 0 < q[0] < 1.0

    def test_unit_square_quad(self):
        """Unit square quad (split as two triangles) should have decent quality."""
        verts = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
        conn = [[0, 1, 2, 3]]

        q = chilmesh.element_quality(verts, conn, metric="aspect_ratio")

        assert len(q) == 1
        assert 0 < q[0] <= 1.0

    def test_degenerate_quad_with_padding(self):
        """Quad padded with -1 (only 3 valid vertices) should return 0."""
        verts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 0.0]])  # degenerate
        conn = [[0, 1, 2, -1]]

        q = chilmesh.element_quality(verts, conn, metric="aspect_ratio")

        assert q[0] == 0.0

    def test_aspect_ratio_range(self):
        """All aspect ratios should be in [0, 1]."""
        # Use donut fixture for realistic data
        mesh = CHILmesh.read_from_fort14("src/chilmesh/data/donut_domain.fort.14")
        conn = [list(row) for row in mesh.connectivity_list]

        q = chilmesh.element_quality(mesh.points, conn, metric="aspect_ratio")

        assert np.all(q >= 0) and np.all(q <= 1.0)


class TestQualityMetrics:
    """Test different quality metrics."""

    def test_min_angle_metric(self):
        """Min angle metric should return values in [0, pi]."""
        mesh = CHILmesh.read_from_fort14("src/chilmesh/data/annulus_200pts.fort.14")
        conn = [list(row) for row in mesh.connectivity_list]

        q = chilmesh.element_quality(mesh.points, conn, metric="min_angle")

        assert np.all(q > 0)  # All non-degenerate triangles should have positive angles
        assert np.all(q <= np.pi)

    def test_max_angle_metric(self):
        """Max angle metric should return values in [0, pi]."""
        mesh = CHILmesh.read_from_fort14("src/chilmesh/data/annulus_200pts.fort.14")
        conn = [list(row) for row in mesh.connectivity_list]

        q = chilmesh.element_quality(mesh.points, conn, metric="max_angle")

        assert np.all(q > 0)
        assert np.all(q <= np.pi)

    def test_aspect_ratio_default(self):
        """Default metric should be aspect_ratio."""
        verts = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3) / 2]])
        conn = [[0, 1, 2]]

        # Call without metric argument
        q1 = chilmesh.element_quality(verts, conn)
        # Call with explicit metric
        q2 = chilmesh.element_quality(verts, conn, metric="aspect_ratio")

        np.testing.assert_allclose(q1, q2)

    def test_invalid_metric_raises(self):
        """Invalid metric should raise ValueError."""
        verts = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 0.5]])
        conn = [[0, 1, 2]]

        with pytest.raises(ValueError, match="Unknown metric"):
            chilmesh.element_quality(verts, conn, metric="invalid_metric")


class TestInputFormats:
    """Test various input format combinations."""

    def test_2d_vertices(self):
        """Standard 2D vertex input."""
        verts = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 0.5]])
        conn = [[0, 1, 2]]

        q = chilmesh.element_quality(verts, conn)

        assert len(q) == 1
        assert 0 <= q[0] <= 1.0

    def test_3d_vertices_uses_xy_only(self):
        """3D vertices should use only x, y coordinates."""
        verts = np.array([[0.0, 0.0, 10.0], [1.0, 0.0, 20.0], [0.5, 0.5, 15.0]])
        conn = [[0, 1, 2]]

        q = chilmesh.element_quality(verts, conn)

        # Should be same as 2D version (z ignored)
        verts_2d = verts[:, :2]
        q_2d = chilmesh.element_quality(verts_2d, conn)

        np.testing.assert_allclose(q, q_2d)

    def test_connectivity_as_numpy_array(self):
        """Connectivity can be numpy array."""
        verts = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 0.5]])
        conn_list = [[0, 1, 2]]
        conn_array = np.array(conn_list)

        q_list = chilmesh.element_quality(verts, conn_list)
        q_array = chilmesh.element_quality(verts, conn_array)

        np.testing.assert_allclose(q_list, q_array)

    def test_multiple_elements(self):
        """Multiple elements in connectivity."""
        verts = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 0.5],
            [2.0, 0.0],
            [1.5, 0.5],
        ])
        conn = [
            [0, 1, 2],  # Triangle 1
            [1, 3, 4],  # Triangle 2
        ]

        q = chilmesh.element_quality(verts, conn)

        assert len(q) == 2
        assert np.all(q >= 0) and np.all(q <= 1.0)


class TestMixedTriQuadConnectivity:
    """Test mixed triangular and quadrilateral elements."""

    def test_mixed_tri_and_quad(self):
        """Mesh with both triangles and quads."""
        verts = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
        ])
        conn = [
            [0, 1, 3, 2],  # Quad
            [1, 4, 5],      # Triangle
        ]

        q = chilmesh.element_quality(verts, conn)

        assert len(q) == 2
        assert np.all(q >= 0) and np.all(q <= 1.0)

    def test_quad_quality_is_min_of_triangles(self):
        """Quad quality should be min of split-triangle qualities."""
        # Create a unit square (will be split diagonally)
        verts = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
        conn = [[0, 1, 2, 3]]

        q_quad = chilmesh.element_quality(verts, conn)

        # Manually compute the two triangle qualities
        # Tri 0-1-2
        q_tri_1 = chilmesh.element_quality(verts, [[0, 1, 2]])
        # Tri 0-2-3
        q_tri_2 = chilmesh.element_quality(verts, [[0, 2, 3]])

        # Quad quality should equal min of the two triangles
        np.testing.assert_allclose(q_quad[0], min(q_tri_1[0], q_tri_2[0]))


class TestExportedInInit:
    """Test that element_quality is properly exported."""

    def test_importable_from_root(self):
        """element_quality should be importable from chilmesh."""
        from chilmesh import element_quality
        assert callable(element_quality)

    def test_in_all_exports(self):
        """element_quality should be in __all__."""
        import chilmesh
        assert "element_quality" in chilmesh.__all__
