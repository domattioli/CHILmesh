"""``signed_area`` regression tests (B2 + B4).

Triangle and quadrilateral signed-area paths must:

- match the analytic shoelace formula on synthetic geometry,
- be strictly positive once ``_ensure_ccw_orientation`` has run, and
- be parametrised over every shipped fixture.
"""
from __future__ import annotations

import numpy as np
import pytest

import chilmesh
from chilmesh import CHILmesh

FIXTURES = ["annulus", "donut", "block_o", "structured"]


@pytest.mark.parametrize("name", FIXTURES)
def test_signed_area_positive_after_ccw(name):
    mesh = getattr(chilmesh.examples, name)()
    areas = mesh.signed_area()
    assert areas.min() > 0, (
        f"{name}: min signed area {areas.min()} is non-positive after CCW reorient"
    )


def test_signed_area_unit_square_quad():
    """Unit square defined CCW must have area exactly 1.0."""
    pts = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    conn = np.array([[0, 1, 2, 3]])
    mesh = CHILmesh(connectivity=conn, points=pts, grid_name="unit_square")
    area = mesh.signed_area()
    np.testing.assert_allclose(area, [1.0], atol=1e-12)


def test_signed_area_quad_shoelace_matches_numpy():
    """Random non-degenerate convex quads must agree with numpy's shoelace."""
    rng = np.random.default_rng(0)
    # Build a small quad mesh from 4 random convex quads.
    quads = []
    pts = []
    for k in range(4):
        cx, cy = rng.uniform(-5, 5, size=2)
        # Sort vertices around the centroid to guarantee a convex CCW quad.
        raw = rng.uniform(-1, 1, size=(4, 2))
        ang = np.arctan2(raw[:, 1] - raw[:, 1].mean(), raw[:, 0] - raw[:, 0].mean())
        ordered = raw[np.argsort(ang)] + np.array([cx, cy])
        base = len(pts)
        pts.extend(ordered.tolist())
        quads.append([base, base + 1, base + 2, base + 3])

    pts_arr = np.array(pts)
    conn = np.array(quads)
    mesh = CHILmesh(connectivity=conn, points=pts_arr, grid_name="random_quads")

    areas = mesh.signed_area()

    # Reference using numpy's shoelace formula.
    expected = np.zeros(len(quads))
    for i, q in enumerate(mesh.connectivity_list):  # post-CCW
        x = mesh.points[q, 0]
        y = mesh.points[q, 1]
        expected[i] = 0.5 * np.abs(
            x[0] * (y[1] - y[3])
            + x[1] * (y[2] - y[0])
            + x[2] * (y[3] - y[1])
            + x[3] * (y[0] - y[2])
        )

    np.testing.assert_allclose(np.abs(areas), expected, rtol=1e-10, atol=1e-12)
    assert (areas > 0).all()


def test_signed_area_unit_triangle():
    pts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    conn = np.array([[0, 1, 2]])
    mesh = CHILmesh(connectivity=conn, points=pts, grid_name="unit_tri")
    np.testing.assert_allclose(mesh.signed_area(), [0.5], atol=1e-12)
