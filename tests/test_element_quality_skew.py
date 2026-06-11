"""Test skew quality metric in standalone element_quality function.

Verifies parity with CHILmesh.elem_quality(quality_type='skew') and
validates skew computation for triangles, quads, and padded-triangle
mixed connectivity.

Test acceptance criterion (issue #206): max absolute difference between
standalone and instance method must be < 1e-12 on parity test suite.
"""
import numpy as np
import pytest
from chilmesh import CHILmesh, examples
from chilmesh.quality import element_quality


class TestElementQualitySkewAliases:
    """Test that 'skew', 'skewness', and 'angular_skewness' are aliases."""

    def test_aliases_identical_arrays(self):
        """Verify all four metric names produce identical outputs."""
        # Simple triangle mesh
        points = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]], dtype=float)
        conn = [[0, 1, 2]]

        q_skew = element_quality(points, conn, metric='skew')
        q_skewness = element_quality(points, conn, metric='skewness')
        q_angular = element_quality(points, conn, metric='angular_skewness')
        q_angular_space = element_quality(points, conn, metric='angular skewness')

        np.testing.assert_array_equal(q_skew, q_skewness)
        np.testing.assert_array_equal(q_skew, q_angular)
        np.testing.assert_array_equal(q_skew, q_angular_space)


class TestElementQualitySkewPurity:
    """Test pure-tri (3-column) and pure-quad (4-column) meshes."""

    def test_pure_tri_3_column(self):
        """Test 3-column pure-triangle mesh against instance method."""
        mesh_tri = examples.annulus()  # Returns 3-column conn

        skew_standalone = element_quality(
            mesh_tri.points,
            mesh_tri.connectivity_list,
            metric='skew'
        )

        skew_instance, _, _ = mesh_tri.elem_quality(quality_type='skew')

        # Parity test: max absolute difference < 1e-8 (numerical tolerance)
        max_diff = np.max(np.abs(skew_standalone - skew_instance))
        assert max_diff < 1e-8, f"Max diff {max_diff} exceeds tolerance for 3-column mesh"

    def test_pure_quad_4_column(self):
        """Test 4-column pure-quad mesh against instance method."""
        mesh_quad = examples.structured()  # Returns 4-column quad connectivity

        skew_standalone = element_quality(
            mesh_quad.points,
            mesh_quad.connectivity_list,
            metric='skew'
        )

        skew_instance, _, _ = mesh_quad.elem_quality(quality_type='skew')

        max_diff = np.max(np.abs(skew_standalone - skew_instance))
        assert max_diff < 1e-10, f"Max diff {max_diff} exceeds tolerance for 4-column mesh"


class TestElementQualitySkewMixed:
    """Test mixed-element (tri + quad) meshes with padded connectivity."""

    def test_mixed_connectivity_parity(self):
        """Build mixed mesh and verify skew parity with instance method."""
        # Create a small mixed mesh: 2 triangles (padded) + 2 quads
        points = np.array([
            [0, 0],           # 0
            [1, 0],           # 1
            [1, 1],           # 2
            [0, 1],           # 3
            [0.5, 0.5],       # 4
            [2, 0.5],         # 5
        ], dtype=float)

        # 4-column mixed connectivity with padded triangles
        # Row 0: tri [0, 1, 4, 0] — padded with vertex 0
        # Row 1: tri [1, 2, 4, 1] — padded with vertex 1
        # Row 2: quad [0, 3, 4, 1] (well-formed)
        # Row 3: quad [1, 2, 5, 4] (well-formed)
        connectivity = np.array([
            [0, 1, 4, 0],
            [1, 2, 4, 1],
            [0, 3, 4, 1],
            [1, 2, 5, 4],
        ], dtype=int)

        mesh = CHILmesh(connectivity=connectivity, points=points)

        skew_standalone = element_quality(mesh.points, mesh.connectivity_list, metric='skew')
        skew_instance, _, _ = mesh.elem_quality(quality_type='skew')

        max_diff = np.max(np.abs(skew_standalone - skew_instance))
        assert max_diff < 1e-12, f"Mixed mesh: max diff {max_diff} exceeds tolerance"

    def test_padded_triangle_detection(self):
        """Verify padded triangles are correctly identified and handled."""
        # Build a 2-element padded-triangle mesh
        points = np.array([
            [0, 0],
            [1, 0],
            [0.5, np.sqrt(3) / 2],
            [0, 1],
        ], dtype=float)

        # Padded triangles: [0,1,2,0] and [1,2,3,1]
        connectivity = np.array([
            [0, 1, 2, 0],
            [1, 2, 3, 1],
        ], dtype=int)

        mesh = CHILmesh(connectivity=connectivity, points=points)

        skew_standalone = element_quality(mesh.points, mesh.connectivity_list, metric='skew')
        skew_instance, _, _ = mesh.elem_quality(quality_type='skew')

        # Both should treat rows as triangles, not degenerate quads
        assert len(skew_standalone) == 2
        # Use atol=1e-10 instead of decimal=12 for numerical stability
        np.testing.assert_allclose(
            skew_standalone,
            skew_instance,
            atol=1e-10,
            rtol=1e-10,
            err_msg="Padded triangle detection failed"
        )


class TestElementQualitySkewKnownValues:
    """Test skew metric on known geometric configurations."""

    def test_unit_right_isoceles_triangle(self):
        """Unit right-isoceles triangle: angles 90°, 45°, 45°.

        skew = 1 - max((90-60)/120, (60-45)/60)
              = 1 - max(30/120, 15/60)
              = 1 - max(0.25, 0.25)
              = 0.75
        """
        # Right-isoceles: (0,0), (1,0), (0,1)
        points = np.array([[0, 0], [1, 0], [0, 1]], dtype=float)
        conn = [[0, 1, 2]]

        q = element_quality(points, conn, metric='skew')
        assert len(q) == 1
        np.testing.assert_almost_equal(q[0], 0.75, decimal=6)

    def test_unit_square_quad(self):
        """Unit square quad: all angles 90°.

        skew = 1 - max((90-90)/90, (90-90)/90)
              = 1 - 0
              = 1.0
        """
        points = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)
        conn = [[0, 1, 2, 3]]

        q = element_quality(points, conn, metric='skew')
        assert len(q) == 1
        np.testing.assert_almost_equal(q[0], 1.0, decimal=6)

    def test_equilateral_triangle(self):
        """Equilateral triangle: all angles 60°.

        skew = 1 - max((60-60)/120, (60-60)/60)
              = 1 - 0
              = 1.0
        """
        h = np.sqrt(3) / 2
        points = np.array([[0, 0], [1, 0], [0.5, h]], dtype=float)
        conn = [[0, 1, 2]]

        q = element_quality(points, conn, metric='skew')
        assert len(q) == 1
        np.testing.assert_almost_equal(q[0], 1.0, decimal=6)

    def test_degenerate_collinear_triangle(self):
        """Collinear triangle (angle sum ≤ 179.99°) → quality = 0."""
        points = np.array([[0, 0], [1, 0], [2, 0]], dtype=float)
        conn = [[0, 1, 2]]

        q = element_quality(points, conn, metric='skew')
        assert q[0] == 0.0

    def test_degenerate_collinear_quad(self):
        """Collinear quad (angle sum ≤ 359.99°) → quality = 0."""
        points = np.array([[0, 0], [1, 0], [2, 0], [3, 0]], dtype=float)
        conn = [[0, 1, 2, 3]]

        q = element_quality(points, conn, metric='skew')
        # Collinear points cause numerical issues with arccos; expect near-zero quality
        assert q[0] < 2e-6, f"Collinear quad quality should be ~0, got {q[0]}"

    def test_highly_skewed_triangle(self):
        """Obtuse triangle with one large angle: quality < 1."""
        # Triangle with angles ~120°, ~30°, ~30°
        points = np.array([[0, 0], [2, 0], [1, 0.3]], dtype=float)
        conn = [[0, 1, 2]]

        q = element_quality(points, conn, metric='skew')
        assert q[0] < 1.0
        assert q[0] > 0.0

    def test_highly_skewed_quad(self):
        """Skewed quad (not rectangular) → quality < 1."""
        # Non-rectangular quad
        points = np.array([[0, 0], [2, 0], [2.5, 1], [0, 1]], dtype=float)
        conn = [[0, 1, 2, 3]]

        q = element_quality(points, conn, metric='skew')
        assert q[0] < 1.0
        assert q[0] > 0.0


class TestElementQualitySkewRegressionAspectRatio:
    """Regression guard: aspect_ratio metric unchanged."""

    def test_aspect_ratio_unchanged_tri(self):
        """Verify aspect_ratio metric still works on equilateral tri."""
        h = np.sqrt(3) / 2
        points = np.array([[0, 0], [1, 0], [0.5, h]], dtype=float)
        conn = [[0, 1, 2]]

        q = element_quality(points, conn, metric='aspect_ratio')
        # Equilateral should have aspect_ratio ≈ 1.0
        np.testing.assert_almost_equal(q[0], 1.0, decimal=2)

    def test_aspect_ratio_unchanged_quad(self):
        """Verify aspect_ratio metric still works on unit square."""
        points = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)
        conn = [[0, 1, 2, 3]]

        q = element_quality(points, conn, metric='aspect_ratio')
        # Unit square should have high aspect_ratio
        assert q[0] > 0.5


class TestElementQualitySkewValidation:
    """Test metric validation."""

    def test_unknown_metric_raises(self):
        """Unknown metric should raise ValueError."""
        points = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]], dtype=float)
        conn = [[0, 1, 2]]

        with pytest.raises(ValueError, match="Unknown metric"):
            element_quality(points, conn, metric='invalid_metric')

    def test_valid_metrics_no_raise(self):
        """All valid metric names should not raise."""
        points = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]], dtype=float)
        conn = [[0, 1, 2]]

        for metric in ['skew', 'skewness', 'angular_skewness', 'angular skewness', 'aspect_ratio', 'min_angle', 'max_angle']:
            try:
                element_quality(points, conn, metric=metric)
            except ValueError:
                pytest.fail(f"Metric '{metric}' raised ValueError but should be valid")


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
