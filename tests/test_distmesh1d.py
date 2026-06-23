"""Test suite for distmesh1d function."""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from chilmesh.admesh_warmstart import distmesh1d


def point_to_segment_distance(point, seg_start, seg_end):
    """Min distance from point to line segment."""
    seg_vec = seg_end - seg_start
    seg_len_sq = np.dot(seg_vec, seg_vec)
    if seg_len_sq < 1e-16:
        return np.linalg.norm(point - seg_start)
    p_vec = point - seg_start
    t = np.dot(p_vec, seg_vec) / seg_len_sq
    t = np.clip(t, 0.0, 1.0)
    closest = seg_start + t * seg_vec
    return np.linalg.norm(point - closest)


def min_distance_to_curve(point, curve_xy):
    """Min distance from point to piecewise-linear curve."""
    M = len(curve_xy)
    min_dist = float('inf')
    for i in range(M):
        seg_start = curve_xy[i]
        seg_end = curve_xy[(i + 1) % M]
        dist = point_to_segment_distance(point, seg_start, seg_end)
        min_dist = min(min_dist, dist)
    return min_dist


def create_clustered_square_boundary(n_points=20):
    """Unit square with nodes clustered near (0, 0) corner."""
    perimeter = 4.0
    t = np.linspace(0, 1, n_points, endpoint=False)
    cluster_mask = t < 0.2
    t_clustered = t[cluster_mask]
    t_rest = t[~cluster_mask]
    t_new = np.concatenate([
        t_clustered ** 2 * 0.2,
        0.2 + (t_rest - 0.2) * (1.0 - 0.2) / (1.0 - 0.2)
    ])
    t_new = np.sort(t_new) % 1.0
    curve_xy = []
    for ti in t_new:
        s = ti * perimeter
        if s <= 1.0:
            curve_xy.append([s, 0.0])
        elif s <= 2.0:
            curve_xy.append([1.0, s - 1.0])
        elif s <= 3.0:
            curve_xy.append([1.0 - (s - 2.0), 1.0])
        else:
            curve_xy.append([0.0, 1.0 - (s - 3.0)])
    return np.array(curve_xy)


class TestDistmesh1D:
    """Test suite for distmesh1d function."""

    def test_shape_count_dtype(self):
        """Output shape (M, 2), dtype float64, same M as input."""
        curve = np.array([
            [0.0, 0.0], [0.33, 0.0], [0.67, 0.0], [1.0, 0.0],
            [1.0, 0.33], [1.0, 0.67], [1.0, 1.0],
            [0.67, 1.0], [0.33, 1.0], [0.0, 1.0],
            [0.0, 0.67], [0.0, 0.33]
        ])
        def uniform_h(pts):
            return np.full(len(pts), 0.2)
        result = distmesh1d(curve, uniform_h, closed=True, niter=10)
        assert result.shape == (12, 2), f"Expected (12, 2), got {result.shape}"
        assert result.dtype == np.float64, f"Expected float64, got {result.dtype}"

    def test_short_curve_passthrough(self):
        """M < 3 returns input unchanged."""
        curve_2 = np.array([[0.0, 0.0], [1.0, 1.0]])
        def dummy_h(pts):
            return np.full(len(pts), 0.1)
        result = distmesh1d(curve_2, dummy_h, closed=True)
        np.testing.assert_allclose(result, curve_2, atol=1e-15)
        curve_1 = np.array([[0.5, 0.5]])
        result_1 = distmesh1d(curve_1, dummy_h, closed=True)
        np.testing.assert_allclose(result_1, curve_1, atol=1e-15)

    def test_input_not_mutated(self):
        """Original input array unchanged after call."""
        curve = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
        curve_copy = curve.copy()
        def uniform_h(pts):
            return np.full(len(pts), 0.2)
        _ = distmesh1d(curve, uniform_h, closed=True, niter=50)
        np.testing.assert_array_equal(curve, curve_copy)

    def test_nodes_stay_on_curve_closed(self):
        """KEY: All output nodes lie on original piecewise-linear curve."""
        angles = np.linspace(0, 2 * np.pi, 13)[:-1]
        curve = np.column_stack([np.cos(angles), np.sin(angles)])
        def uniform_h(pts):
            return np.full(len(pts), 0.3)
        result = distmesh1d(curve, uniform_h, closed=True, niter=200)
        for i, pt in enumerate(result):
            min_dist = min_distance_to_curve(pt, curve)
            assert min_dist < 1e-6, (
                f"Point {i} is {min_dist:.2e} from curve (should be ~0)"
            )

    def test_open_curve_endpoints_fixed(self):
        """With closed=False, first/last points equal input endpoints."""
        curve = np.array([
            [0.0, 0.0], [0.2, 0.2], [0.4, 0.4],
            [0.6, 0.6], [0.8, 0.8], [1.0, 1.0]
        ])
        def uniform_h(pts):
            return np.full(len(pts), 0.15)
        result = distmesh1d(curve, uniform_h, closed=False, niter=100)
        np.testing.assert_allclose(result[0], curve[0], atol=1e-9)
        np.testing.assert_allclose(result[-1], curve[-1], atol=1e-9)

    def test_determinism(self):
        """Two identical calls produce identical output."""
        curve = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
        def uniform_h(pts):
            return np.full(len(pts), 0.2)
        result1 = distmesh1d(curve, uniform_h, closed=True, niter=50)
        result2 = distmesh1d(curve, uniform_h, closed=True, niter=50)
        np.testing.assert_array_equal(result1, result2)

    def test_uniformization_reduces_gap_variance(self):
        """Uniform h_fn reduces variance of inter-node gaps."""
        curve = create_clustered_square_boundary(n_points=20)
        def uniform_h(pts):
            return np.full(len(pts), 0.15)
        # Before
        diffs_before = np.diff(curve, axis=0, append=curve[:1])
        gaps_before = np.linalg.norm(diffs_before, axis=1)
        var_before = np.var(gaps_before)
        # After
        result = distmesh1d(curve, uniform_h, closed=True, niter=200)
        diffs_after = np.diff(result, axis=0, append=result[:1])
        gaps_after = np.linalg.norm(diffs_after, axis=1)
        var_after = np.var(gaps_after)
        assert var_after < var_before, (
            f"Variance not reduced: {var_before:.6f} → {var_after:.6f}"
        )
