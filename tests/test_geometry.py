"""Tests for chilmesh.geometry module (CRS-aware edge length computation)."""
from __future__ import annotations

import numpy as np
import pytest

from chilmesh.geometry import (
    haversine_m,
    edge_lengths,
    EARTH_RADIUS_M,
    convex_hull,
    is_antimeridian_wrapping,
    split_antimeridian_bbox,
    bbox_iou,
    hausdorff_distance,
)


class TestHaversine:
    """Tests for the haversine_m function."""

    def test_zero_distance(self):
        """Test that identical points return zero distance."""
        result = haversine_m(0.0, 0.0, 0.0, 0.0)
        assert result == 0.0

    def test_one_degree_longitude_at_equator(self):
        """Test 1 degree of longitude at equator ≈ 111,195 m."""
        # At equator, 1° of longitude ≈ 111195 m (pi/180 * EARTH_RADIUS_M)
        result = haversine_m(0.0, 0.0, 1.0, 0.0)
        expected = np.pi / 180.0 * EARTH_RADIUS_M
        assert result == pytest.approx(expected, rel=1e-3)

    def test_antipodal_points(self):
        """Test antipodal points (0,0) vs (180,0) ≈ pi * EARTH_RADIUS_M."""
        # Half the Earth circumference (clamped by min in the formula)
        result = haversine_m(0.0, 0.0, 180.0, 0.0)
        expected = np.pi * EARTH_RADIUS_M
        assert result == pytest.approx(expected, rel=1e-6)

    def test_vectorized_scalar_inputs(self):
        """Test vectorized computation with array inputs."""
        lons1 = np.array([0.0, 0.0, 45.0])
        lats1 = np.array([0.0, 0.0, 0.0])
        lons2 = np.array([1.0, 0.0, 46.0])
        lats2 = np.array([0.0, 1.0, 0.0])
        result = haversine_m(lons1, lats1, lons2, lats2)
        assert result.shape == (3,)
        assert isinstance(result, np.ndarray)

    def test_scalar_return_for_scalar_input(self):
        """Test that scalar inputs return a float, not an ndarray."""
        result = haversine_m(0.0, 0.0, 1.0, 0.0)
        assert isinstance(result, float)
        assert not isinstance(result, np.ndarray)


class TestEdgeLengths:
    """Tests for the edge_lengths function."""

    def test_cartesian_3_4_5_triangle(self):
        """Test that a 3-4-5 right triangle has distance 5 in Cartesian coords."""
        result = edge_lengths((0.0, 0.0), (3.0, 4.0), crs="cartesian")
        assert result == 5.0

    def test_cartesian_scalar_return(self):
        """Test that scalar points return float from Cartesian mode."""
        result = edge_lengths((0.0, 0.0), (3.0, 4.0), crs="cartesian")
        assert isinstance(result, float)

    def test_spherical_scalar_matches_haversine(self):
        """Test that spherical edge_lengths matches direct haversine_m call."""
        p1 = (10.0, 20.0)
        p2 = (11.0, 21.0)
        result_edge = edge_lengths(p1, p2, crs="spherical")
        result_hav = haversine_m(p1[0], p1[1], p2[0], p2[1])
        assert result_edge == pytest.approx(result_hav)

    def test_invalid_crs_raises_valueerror(self):
        """Test that an invalid CRS string raises ValueError."""
        with pytest.raises(ValueError, match="unsupported crs"):
            edge_lengths((0.0, 0.0), (1.0, 1.0), crs="bogus")

    def test_mismatched_shapes_raises_valueerror(self):
        """Test that mismatched point array shapes raise ValueError."""
        p1 = np.array([[0.0, 0.0]])  # shape (1, 2)
        p2 = np.array([[1.0, 1.0], [2.0, 2.0]])  # shape (2, 2)
        with pytest.raises(ValueError, match="must be matching"):
            edge_lengths(p1, p2, crs="cartesian")

    def test_vectorized_cartesian(self):
        """Test vectorized Cartesian distance computation."""
        p1 = np.array([[0.0, 0.0], [1.0, 1.0]])
        p2 = np.array([[3.0, 4.0], [4.0, 5.0]])
        result = edge_lengths(p1, p2, crs="cartesian")
        # Distance from (0,0) to (3,4) = sqrt(9+16) = 5.0
        # Distance from (1,1) to (4,5) = sqrt(9+16) = 5.0
        expected = np.array([5.0, 5.0])
        np.testing.assert_allclose(result, expected)

    def test_vectorized_spherical(self):
        """Test vectorized spherical distance computation."""
        lons1 = np.array([0.0, 45.0])
        lats1 = np.array([0.0, 0.0])
        lons2 = np.array([1.0, 46.0])
        lats2 = np.array([0.0, 0.0])
        p1 = np.column_stack([lons1, lats1])
        p2 = np.column_stack([lons2, lats2])
        result = edge_lengths(p1, p2, crs="spherical")
        assert result.shape == (2,)
        assert isinstance(result, np.ndarray)


class TestConvexHull:
    """Tests for the convex_hull function."""

    def test_unit_square(self):
        """Test that a unit square with interior point produces 4 hull vertices."""
        pts = np.array([[0, 0], [1, 1], [1, 0], [0, 1], [0.5, 0.5]])
        hull = convex_hull(pts)
        assert hull.shape == (4, 2)

    def test_ccw_order(self):
        """Test that hull vertices are in counter-clockwise order."""
        pts = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
        hull = convex_hull(pts)
        # First vertex should be lexicographically smallest: (0, 0)
        assert np.allclose(hull[0], [0, 0])

    def test_collinear_points_dropped(self):
        """Test that collinear points are excluded from the hull."""
        pts = np.array([[0, 0], [0.5, 0.5], [1, 1], [2, 2]])
        hull = convex_hull(pts)
        # All points are collinear, so hull should be the two endpoints
        assert hull.shape == (2, 2)
        assert np.allclose(hull[0], [0, 0])
        assert np.allclose(hull[1], [2, 2])

    def test_two_point_degenerate(self):
        """Test that a 2-point degenerate case returns both points sorted."""
        pts = np.array([[1, 1], [0, 0]])
        hull = convex_hull(pts)
        assert hull.shape == (2, 2)
        assert np.allclose(hull[0], [0, 0])
        assert np.allclose(hull[1], [1, 1])

    def test_single_point_degenerate(self):
        """Test that a single point is returned as-is."""
        pts = np.array([[5, 5]])
        hull = convex_hull(pts)
        assert hull.shape == (1, 2)
        assert np.allclose(hull[0], [5, 5])


class TestIsAntimeridianWrapping:
    """Tests for the is_antimeridian_wrapping function."""

    def test_wrapping_case(self):
        """Test that min_lon > max_lon returns True."""
        bbox = (170, -10, -170, 10)
        assert is_antimeridian_wrapping(bbox) is True

    def test_non_wrapping_case(self):
        """Test that min_lon <= max_lon returns False."""
        bbox = (-10, -10, 10, 10)
        assert is_antimeridian_wrapping(bbox) is False

    def test_edge_case_zero_crossing(self):
        """Test edge case where box straddles zero."""
        bbox = (-5, -5, 5, 5)
        assert is_antimeridian_wrapping(bbox) is False

    def test_edge_case_at_dateline(self):
        """Test edge case at exactly 180 / -180."""
        bbox = (179, -10, -179, 10)
        assert is_antimeridian_wrapping(bbox) is True


class TestSplitAntimeridianBbox:
    """Tests for the split_antimeridian_bbox function."""

    def test_wrapping_split(self):
        """Test that wrapping bbox is split into two non-wrapping parts."""
        bbox = (170, -10, -170, 10)
        result = split_antimeridian_bbox(bbox)
        assert len(result) == 2
        # Eastern part: (170, -10, 180, 10)
        assert np.allclose(result[0], (170, -10, 180.0, 10))
        # Western part: (-180, -10, -170, 10)
        assert np.allclose(result[1], (-180.0, -10, -170, 10))

    def test_non_wrapping_no_split(self):
        """Test that non-wrapping bbox is returned as single-element list."""
        bbox = (-10, -10, 10, 10)
        result = split_antimeridian_bbox(bbox)
        assert len(result) == 1
        assert result[0] == bbox


class TestBboxIOU:
    """Tests for the bbox_iou function."""

    def test_identical_boxes(self):
        """Test that identical boxes have IoU = 1.0."""
        bbox = (0, 0, 2, 2)
        assert bbox_iou(bbox, bbox) == pytest.approx(1.0)

    def test_disjoint_boxes(self):
        """Test that disjoint boxes have IoU = 0.0."""
        bbox_a = (0, 0, 2, 2)
        bbox_b = (3, 3, 5, 5)
        assert bbox_iou(bbox_a, bbox_b) == pytest.approx(0.0)

    def test_half_overlap(self):
        """Test half-overlapping boxes: (0,0,2,2) vs (1,0,3,2)."""
        bbox_a = (0, 0, 2, 2)
        bbox_b = (1, 0, 3, 2)
        # Intersection: 1 x 2 = 2
        # Union: 4 + 4 - 2 = 6
        # IoU = 2 / 6 ≈ 0.333
        expected = 2.0 / 6.0
        assert bbox_iou(bbox_a, bbox_b) == pytest.approx(expected, rel=1e-5)

    def test_wrapping_iou(self):
        """Test IoU with antimeridian-wrapping box."""
        # A wrapping box should compute IoU correctly
        bbox_wrapping = (170, -10, -170, 10)
        # Split parts: [(170, -10, 180, 10), (-180, -10, -170, 10)]
        # IoU with itself should be 1.0
        assert bbox_iou(bbox_wrapping, bbox_wrapping) == pytest.approx(1.0)


class TestHausdorffDistance:
    """Tests for the hausdorff_distance function."""

    def test_identical_sets_cartesian(self):
        """Test that identical sets have Hausdorff distance 0."""
        pts = np.array([[0, 0], [1, 1], [2, 0]])
        assert hausdorff_distance(pts, pts, crs="cartesian") == pytest.approx(0.0)

    def test_two_point_sets_cartesian(self):
        """Test Hausdorff distance for two-point sets: (0,0)/(3,4) vs (0,0)."""
        pts_a = np.array([[0, 0], [3, 4]])
        pts_b = np.array([[0, 0]])
        # Distance from (0,0) to (0,0) = 0
        # Distance from (3,4) to (0,0) = 5
        # max(0, 5) = 5, and max(5, 0) = 5
        assert hausdorff_distance(pts_a, pts_b, crs="cartesian") == pytest.approx(5.0)

    def test_spherical_path(self):
        """Test Hausdorff distance in spherical CRS."""
        pts_a = np.array([[0, 0]])
        pts_b = np.array([[1, 0]])
        result = hausdorff_distance(pts_a, pts_b, crs="spherical")
        # Should return the haversine distance
        expected = haversine_m(0, 0, 1, 0)
        assert result == pytest.approx(expected)

    def test_invalid_crs_raises_error(self):
        """Test that invalid crs raises ValueError."""
        pts_a = np.array([[0, 0]])
        pts_b = np.array([[1, 1]])
        with pytest.raises(ValueError, match="unsupported crs"):
            hausdorff_distance(pts_a, pts_b, crs="bogus")

    def test_empty_input_raises_error(self):
        """Test that empty point sets raise ValueError."""
        pts_a = np.array([[0, 0]])
        pts_b = np.array([]).reshape(0, 2)
        with pytest.raises(ValueError, match="non-empty"):
            hausdorff_distance(pts_a, pts_b, crs="cartesian")


class TestImportSurface:
    """Tests for public API exports."""

    def test_haversine_importable(self):
        """Test that haversine_m is importable from chilmesh package."""
        from chilmesh import haversine_m
        assert callable(haversine_m)

    def test_edge_lengths_importable(self):
        """Test that edge_lengths is importable from chilmesh package."""
        from chilmesh import edge_lengths
        assert callable(edge_lengths)

    def test_earth_radius_importable(self):
        """Test that EARTH_RADIUS_M is importable from chilmesh package."""
        from chilmesh import EARTH_RADIUS_M
        assert isinstance(EARTH_RADIUS_M, (int, float))
        assert EARTH_RADIUS_M > 0

    def test_convex_hull_importable(self):
        """Test that convex_hull is importable from chilmesh package."""
        from chilmesh import convex_hull
        assert callable(convex_hull)

    def test_is_antimeridian_wrapping_importable(self):
        """Test that is_antimeridian_wrapping is importable from chilmesh package."""
        from chilmesh import is_antimeridian_wrapping
        assert callable(is_antimeridian_wrapping)

    def test_split_antimeridian_bbox_importable(self):
        """Test that split_antimeridian_bbox is importable from chilmesh package."""
        from chilmesh import split_antimeridian_bbox
        assert callable(split_antimeridian_bbox)

    def test_bbox_iou_importable(self):
        """Test that bbox_iou is importable from chilmesh package."""
        from chilmesh import bbox_iou
        assert callable(bbox_iou)

    def test_hausdorff_distance_importable(self):
        """Test that hausdorff_distance is importable from chilmesh package."""
        from chilmesh import hausdorff_distance
        assert callable(hausdorff_distance)
