"""equiangle_skewness metric (#217): raw ANSYS/Fluent Qeas (0 = ideal,
1 = degenerate), the exact complement of the 'skew' quality metric on
non-degenerate elements."""
import numpy as np

from chilmesh import element_quality


def test_equilateral_triangle_zero_skewness():
    h = np.sqrt(3) / 2
    pts = np.array([[0, 0], [1, 0], [0.5, h]], dtype=float)
    q = element_quality(pts, [[0, 1, 2]], metric="equiangle_skewness")
    np.testing.assert_almost_equal(q[0], 0.0, decimal=6)


def test_right_isoceles_triangle_skewness():
    # angles 90/45/45 -> Qeas = max((90-60)/120, (60-45)/60) = 0.25
    pts = np.array([[0, 0], [1, 0], [0, 1]], dtype=float)
    q = element_quality(pts, [[0, 1, 2]], metric="equiangle_skewness")
    np.testing.assert_almost_equal(q[0], 0.25, decimal=6)


def test_unit_square_zero_skewness():
    pts = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)
    q = element_quality(pts, [[0, 1, 2, 3]], metric="equiangle_skewness")
    np.testing.assert_almost_equal(q[0], 0.0, decimal=6)


def test_degenerate_triangle_max_skewness():
    pts = np.array([[0, 0], [1, 0], [2, 0]], dtype=float)
    q = element_quality(pts, [[0, 1, 2]], metric="equiangle_skewness")
    np.testing.assert_almost_equal(q[0], 1.0, decimal=6)


def test_complement_of_skew_nondegenerate():
    # skew == 1 - Qeas on non-degenerate elements (tri + quad).
    pts_t = np.array([[0, 0], [2, 0], [1, 0.3]], dtype=float)
    conn_t = [[0, 1, 2]]
    pts_q = np.array([[0, 0], [2, 0], [2.5, 1], [0, 1]], dtype=float)
    conn_q = [[0, 1, 2, 3]]
    for pts, conn in ((pts_t, conn_t), (pts_q, conn_q)):
        skew = element_quality(pts, conn, metric="skew")
        eas = element_quality(pts, conn, metric="equiangle_skewness")
        np.testing.assert_allclose(skew + eas, 1.0, atol=1e-9)


def test_aliases_identical():
    pts = np.array([[0, 0], [1, 0], [0, 1]], dtype=float)
    conn = [[0, 1, 2]]
    a = element_quality(pts, conn, metric="equiangle_skewness")
    b = element_quality(pts, conn, metric="equiangle skewness")
    c = element_quality(pts, conn, metric="eas")
    np.testing.assert_array_equal(a, b)
    np.testing.assert_array_equal(a, c)
