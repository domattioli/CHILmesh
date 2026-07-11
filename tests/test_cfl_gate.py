"""Tests for CFL/Courant quality gate (courant_number and cfl_gate)."""
from __future__ import annotations

import numpy as np
import pytest
from chilmesh import courant_number, cfl_gate, haversine_m


class TestCourantNumber:
    """Test courant_number() function."""

    def test_known_value_cartesian(self):
        """Single edge (0,0)-(100,0) with depths [10,10], dt=1.
        Expected: C = sqrt(9.81*10)/100 ≈ 0.09905
        """
        points = np.array([[0, 0], [100, 0]], dtype=float)
        edges = np.array([[0, 1]], dtype=int)
        depths = np.array([10, 10], dtype=float)
        dt = 1.0

        c = courant_number(points, edges, depths, dt, crs="cartesian")

        expected = np.sqrt(9.81 * 10) / 100
        np.testing.assert_allclose(c[0], expected, rtol=1e-6)

    def test_deeper_endpoint_wins(self):
        """Depths [4, 9]: h=9 should be used in formula.
        Edge length 100, dt=1.
        Expected: C = sqrt(9.81*9)/100
        """
        points = np.array([[0, 0], [100, 0]], dtype=float)
        edges = np.array([[0, 1]], dtype=int)
        depths = np.array([4, 9], dtype=float)
        dt = 1.0

        c = courant_number(points, edges, depths, dt)

        expected = np.sqrt(9.81 * 9) / 100
        np.testing.assert_allclose(c[0], expected, rtol=1e-6)

    def test_dry_edge_nans(self):
        """Both depths <= 0: should produce NaN."""
        points = np.array([[0, 0], [100, 0]], dtype=float)
        edges = np.array([[0, 1]], dtype=int)
        depths = np.array([-5, 0], dtype=float)
        dt = 1.0

        c = courant_number(points, edges, depths, dt)

        assert np.isnan(c[0])

    def test_degenerate_zero_length_edge(self):
        """Edge with zero length: should produce NaN."""
        points = np.array([[0, 0], [0, 0]], dtype=float)
        edges = np.array([[0, 1]], dtype=int)
        depths = np.array([10, 10], dtype=float)
        dt = 1.0

        c = courant_number(points, edges, depths, dt)

        assert np.isnan(c[0])

    def test_dt_invalid_raises_valueerror(self):
        """dt <= 0 should raise ValueError."""
        points = np.array([[0, 0], [100, 0]], dtype=float)
        edges = np.array([[0, 1]], dtype=int)
        depths = np.array([10, 10], dtype=float)

        with pytest.raises(ValueError, match="dt must be positive"):
            courant_number(points, edges, depths, dt=0)

        with pytest.raises(ValueError, match="dt must be positive"):
            courant_number(points, edges, depths, dt=-1.0)

    def test_mismatched_shapes_raises_valueerror(self):
        """Mismatched points/edges/depths shapes should raise ValueError."""
        points = np.array([[0, 0], [100, 0]], dtype=float)
        edges = np.array([[0, 1]], dtype=int)
        depths = np.array([10], dtype=float)  # Only 1 depth, but need 2

        with pytest.raises(ValueError):
            courant_number(points, edges, depths, dt=1.0)

    def test_spherical_crs(self):
        """Spherical CRS (lon/lat degrees).
        Edge (0,0)-(1,0) degrees, depth 10, dt=1.
        Expected dx from haversine.
        """
        points = np.array([[0, 0], [1, 0]], dtype=float)
        edges = np.array([[0, 1]], dtype=int)
        depths = np.array([10, 10], dtype=float)
        dt = 1.0

        c = courant_number(points, edges, depths, dt, crs="spherical")

        # Compute expected dx via haversine
        dx_expected = haversine_m(0, 0, 1, 0)
        c_expected = np.sqrt(9.81 * 10) / dx_expected
        np.testing.assert_allclose(c[0], c_expected, rtol=1e-9)

    def test_empty_edges_array(self):
        """Empty edges array should return empty array."""
        points = np.array([[0, 0]], dtype=float)
        edges = np.array([], dtype=int).reshape(0, 2)
        depths = np.array([10], dtype=float)
        dt = 1.0

        c = courant_number(points, edges, depths, dt)

        assert len(c) == 0
        assert c.dtype == float


class TestCflGate:
    """Test cfl_gate() function."""

    def test_pass_all_edges_below_threshold(self):
        """All edges have C < courant_max: status='pass', n_over=0."""
        points = np.array([[0, 0], [100, 0], [200, 0]], dtype=float)
        edges = np.array([[0, 1], [1, 2]], dtype=int)
        depths = np.array([10, 10, 10], dtype=float)
        dt = 1.0
        courant_max = 1.0

        result = cfl_gate(points, edges, depths, dt, courant_max=courant_max)

        assert result["status"] == "pass"
        assert result["n_over"] == 0
        assert result["n_wet"] == 2
        assert result["n_skipped"] == 0
        assert len(result["offenders"]) == 0

    def test_fail_edges_exceed_threshold(self):
        """Mix edges where some exceed courant_max."""
        # Short edges with high depth -> high C values
        points = np.array([[0, 0], [5, 0], [25, 0]], dtype=float)
        edges = np.array([[0, 1], [1, 2]], dtype=int)
        depths = np.array([20, 20, 20], dtype=float)
        dt = 1.0
        courant_max = 0.5

        result = cfl_gate(points, edges, depths, dt, courant_max=courant_max)

        assert result["status"] == "fail"
        assert result["n_over"] >= 1  # Edge 0-1 is short, high C
        assert result["n_wet"] >= 1
        assert len(result["offenders"]) >= 1
        # First offender should be the one with highest C
        assert result["offenders"][0]["courant"] > courant_max

    def test_offenders_sorted_descending(self):
        """Offenders should be sorted by courant (descending)."""
        points = np.array([[0, 0], [5, 0], [15, 0], [30, 0]], dtype=float)
        edges = np.array([[0, 1], [1, 2], [2, 3]], dtype=int)
        depths = np.array([10, 10, 10, 10], dtype=float)
        dt = 1.0
        courant_max = 0.5

        result = cfl_gate(points, edges, depths, dt, courant_max=courant_max)

        assert result["n_over"] >= 1
        for i in range(len(result["offenders"]) - 1):
            assert result["offenders"][i]["courant"] >= result["offenders"][i+1]["courant"]

    def test_offenders_capped_at_max_offenders(self):
        """Number of offenders should not exceed max_offenders."""
        n_nodes = 52  # Create 51 edges, all with high C
        points = np.column_stack([np.arange(n_nodes) * 5, np.zeros(n_nodes)])
        edges = np.array([[i, i+1] for i in range(n_nodes - 1)], dtype=int)
        depths = np.full(n_nodes, 10, dtype=float)
        dt = 1.0
        courant_max = 0.5
        max_offenders = 25

        result = cfl_gate(points, edges, depths, dt, courant_max=courant_max,
                         max_offenders=max_offenders)

        assert len(result["offenders"]) <= max_offenders

    def test_dry_edges_excluded_from_stats(self):
        """Dry edges should be counted in n_skipped, excluded from stats.
        An edge is dry if max(h[endpoint0], h[endpoint1]) <= 0.
        """
        points = np.array([[0, 0], [100, 0], [200, 0]], dtype=float)
        edges = np.array([[0, 1], [1, 2]], dtype=int)
        depths = np.array([10, -5, -10], dtype=float)  # Edge 1-2 has both endpoints dry
        dt = 1.0

        result = cfl_gate(points, edges, depths, dt, courant_max=1.0)

        assert result["n_skipped"] >= 1  # Edge 1-2 is dry (max(-5, -10) = -5 <= 0)
        # The dry edge should not be in offenders
        for off in result["offenders"]:
            assert not np.isnan(off["courant"])

    def test_stats_present_and_valid(self):
        """Result dict should contain all expected keys with valid values."""
        points = np.array([[0, 0], [100, 0]], dtype=float)
        edges = np.array([[0, 1]], dtype=int)
        depths = np.array([10, 10], dtype=float)
        dt = 1.0

        result = cfl_gate(points, edges, depths, dt)

        # Check all expected keys
        expected_keys = {
            "status", "n_over", "n_wet", "n_skipped",
            "max", "median", "p95", "offenders", "thresholds"
        }
        assert set(result.keys()) == expected_keys

        # Check types
        assert isinstance(result["status"], str)
        assert result["status"] in ("pass", "fail")
        assert isinstance(result["n_over"], int)
        assert isinstance(result["n_wet"], int)
        assert isinstance(result["n_skipped"], int)
        assert isinstance(result["max"], float)
        assert isinstance(result["median"], float)
        assert isinstance(result["p95"], float)
        assert isinstance(result["offenders"], list)
        assert isinstance(result["thresholds"], dict)

    def test_offender_dict_structure(self):
        """Each offender dict should have required keys."""
        # Create scenario with offenders
        points = np.array([[0, 0], [5, 0], [20, 0]], dtype=float)
        edges = np.array([[0, 1], [1, 2]], dtype=int)
        depths = np.array([10, 10, 10], dtype=float)
        dt = 1.0
        courant_max = 1.0

        result = cfl_gate(points, edges, depths, dt, courant_max=courant_max)

        if result["offenders"]:
            off = result["offenders"][0]
            expected_keys = {"edge_index", "nodes", "courant", "edge_length", "depth"}
            assert set(off.keys()) == expected_keys
            assert isinstance(off["edge_index"], (int, np.integer))
            assert isinstance(off["nodes"], list) and len(off["nodes"]) == 2
            assert isinstance(off["courant"], float)
            assert isinstance(off["edge_length"], float)
            assert isinstance(off["depth"], float)

    def test_empty_edges_array(self):
        """Empty edges array should return all-zero stats."""
        points = np.array([[0, 0]], dtype=float)
        edges = np.array([], dtype=int).reshape(0, 2)
        depths = np.array([10], dtype=float)
        dt = 1.0

        result = cfl_gate(points, edges, depths, dt)

        assert result["status"] == "pass"
        assert result["n_over"] == 0
        assert result["n_wet"] == 0
        assert result["n_skipped"] == 0
        assert result["max"] == 0.0
        assert result["median"] == 0.0
        assert result["p95"] == 0.0
        assert len(result["offenders"]) == 0

    def test_thresholds_recorded(self):
        """Result should record the courant_max threshold."""
        courant_max = 0.8
        points = np.array([[0, 0], [100, 0]], dtype=float)
        edges = np.array([[0, 1]], dtype=int)
        depths = np.array([10, 10], dtype=float)
        dt = 1.0

        result = cfl_gate(points, edges, depths, dt, courant_max=courant_max)

        assert result["thresholds"]["courant_max"] == courant_max


class TestApiSurface:
    """Test public API import and availability."""

    def test_import_from_chilmesh(self):
        """Both functions should be importable from chilmesh top level."""
        from chilmesh import courant_number, cfl_gate
        assert callable(courant_number)
        assert callable(cfl_gate)

    def test_in_all_export(self):
        """Both should be in __all__ list."""
        import chilmesh
        assert "courant_number" in chilmesh.__all__
        assert "cfl_gate" in chilmesh.__all__
