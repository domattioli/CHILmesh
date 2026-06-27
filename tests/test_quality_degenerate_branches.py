"""Tests for untested degenerate branches in quality.py.

Targets previously-uncovered code paths:
1. Lines 147-149: _quad_quality non-skew metric split into two triangles
2. Line 269: element_quality with fewer than 3 valid vertices

Plus behavioral pins for degenerate (collinear) triangles under the skew
metric. NOTE: the ``angle_sum <= 179.99`` guard in ``_triangle_quality``
(line 70) is an unreachable defensive floor for valid float geometry -- a
2M-sample near-collinear brute force never produced a computed angle sum
below 179.99999, because the three independent arccos values sum to ~180
for any real degenerate config, and the skew *formula* itself already
returns 0.0 for those. These tests assert that degenerate-returns-0.0
contract (via the formula), not the dead guard.
"""
from __future__ import annotations

import numpy as np
import chilmesh
from chilmesh.quality import _triangle_quality, _quad_quality


class TestTriangleQualitySkewDegenerate:
    """Degenerate collinear triangles return 0.0 under the skew metric."""

    def test_collinear_skew_metric(self):
        """Collinear triangle with skew metric returns 0.0."""
        # Three collinear points: all at y=0
        v0 = np.array([0.0, 0.0])
        v1 = np.array([1.0, 0.0])
        v2 = np.array([2.0, 0.0])

        q = _triangle_quality(v0, v1, v2, metric="skew")

        assert q == 0.0

    def test_collinear_angular_skewness_metric(self):
        """Collinear triangle with angular_skewness metric returns 0.0."""
        # Three collinear points
        v0 = np.array([0.0, 0.0])
        v1 = np.array([1.0, 0.0])
        v2 = np.array([2.0, 0.0])

        q = _triangle_quality(v0, v1, v2, metric="angular_skewness")

        assert q == 0.0

    def test_public_api_collinear_skew_via_element_quality(self):
        """Via public API: element_quality with collinear triangle, skew metric."""
        verts = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
        conn = [[0, 1, 2]]

        q = chilmesh.element_quality(verts, conn, metric="skew")

        assert q[0] == 0.0

    def test_public_api_collinear_angular_skewness_via_element_quality(self):
        """Via public API: element_quality with collinear triangle, angular_skewness metric."""
        verts = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
        conn = [[0, 1, 2]]

        q = chilmesh.element_quality(verts, conn, metric="angular_skewness")

        assert q[0] == 0.0


class TestQuadQualityNonSkewSplit:
    """Test lines 147-149: _quad_quality else branch for non-skew metrics."""

    def test_unit_square_min_angle_split(self):
        """Unit square split into two triangles, min_angle metric.

        Quad split as (v0,v1,v2) and (v0,v2,v3).
        Result should equal min of the two triangle qualities.
        """
        v0 = np.array([0.0, 0.0])
        v1 = np.array([1.0, 0.0])
        v2 = np.array([1.0, 1.0])
        v3 = np.array([0.0, 1.0])

        q_quad = _quad_quality(v0, v1, v2, v3, metric="min_angle")

        # Manually compute the two triangle qualities
        q_tri_1 = _triangle_quality(v0, v1, v2, metric="min_angle")
        q_tri_2 = _triangle_quality(v0, v2, v3, metric="min_angle")
        expected = min(q_tri_1, q_tri_2)

        np.testing.assert_allclose(q_quad, expected)

    def test_unit_square_aspect_ratio_split(self):
        """Unit square split into two triangles, aspect_ratio metric."""
        v0 = np.array([0.0, 0.0])
        v1 = np.array([1.0, 0.0])
        v2 = np.array([1.0, 1.0])
        v3 = np.array([0.0, 1.0])

        q_quad = _quad_quality(v0, v1, v2, v3, metric="aspect_ratio")

        # Manually compute the two triangle qualities
        q_tri_1 = _triangle_quality(v0, v1, v2, metric="aspect_ratio")
        q_tri_2 = _triangle_quality(v0, v2, v3, metric="aspect_ratio")
        expected = min(q_tri_1, q_tri_2)

        np.testing.assert_allclose(q_quad, expected)

    def test_quad_quality_returns_float(self):
        """_quad_quality always returns a float type."""
        v0 = np.array([0.0, 0.0])
        v1 = np.array([1.0, 0.0])
        v2 = np.array([1.0, 1.0])
        v3 = np.array([0.0, 1.0])

        q = _quad_quality(v0, v1, v2, v3, metric="min_angle")

        assert isinstance(q, float)

    def test_unit_square_max_angle_split(self):
        """Unit square with max_angle metric."""
        v0 = np.array([0.0, 0.0])
        v1 = np.array([1.0, 0.0])
        v2 = np.array([1.0, 1.0])
        v3 = np.array([0.0, 1.0])

        q_quad = _quad_quality(v0, v1, v2, v3, metric="max_angle")

        q_tri_1 = _triangle_quality(v0, v1, v2, metric="max_angle")
        q_tri_2 = _triangle_quality(v0, v2, v3, metric="max_angle")
        expected = min(q_tri_1, q_tri_2)

        np.testing.assert_allclose(q_quad, expected)


class TestElementQualityDegenerateFewVertices:
    """Test line 269: element_quality returns 0.0 for < 3 valid vertices."""

    def test_two_vertex_element(self):
        """Element with only 2 vertices (degenerate) returns quality 0.0."""
        verts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        # Only 2 valid vertices (indices 0, 1); index 2 used as third to form an element
        conn = [[0, 1]]  # Only 2 vertices in this element

        q = chilmesh.element_quality(verts, conn)

        assert len(q) == 1
        assert q[0] == 0.0

    def test_padded_tri_one_valid_vertex_in_quad_row(self):
        """Quad row with -1 padding but only 1 valid vertex: quality 0.0.

        A 4-row element [0, -1, -1, -1] after -1 filtering has only 1 valid vertex.
        """
        verts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        # Quad row with only vertex 0 valid; rest are -1 padding
        conn = [[0, -1, -1, -1]]

        q = chilmesh.element_quality(verts, conn)

        assert len(q) == 1
        assert q[0] == 0.0

    def test_single_vertex_element(self):
        """Single-vertex element returns quality 0.0."""
        verts = np.array([[0.0, 0.0], [1.0, 0.0]])
        conn = [[0]]

        q = chilmesh.element_quality(verts, conn)

        assert len(q) == 1
        assert q[0] == 0.0

    def test_empty_element(self):
        """Empty element (all -1) returns quality 0.0."""
        verts = np.array([[0.0, 0.0], [1.0, 0.0]])
        conn = [[-1, -1]]

        q = chilmesh.element_quality(verts, conn)

        assert len(q) == 1
        assert q[0] == 0.0
