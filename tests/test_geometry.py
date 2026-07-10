"""Tests for chilmesh.geometry module (CRS-aware edge length computation)."""
from __future__ import annotations

import numpy as np
import pytest

from chilmesh.geometry import haversine_m, edge_lengths, EARTH_RADIUS_M


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
