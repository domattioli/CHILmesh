"""``interior_angles`` regression tests (B5).

For triangles, the three interior angles sum to 180 deg. For quads, the
four interior angles sum to 360 deg. Degenerate quads (e.g. concave or
near-collinear) must not emit NaN — they must produce real numbers so
``elem_quality``'s "<= 359.99 -> zero out" guard does its job.
"""
from __future__ import annotations

import numpy as np
import pytest

import chilmesh
from chilmesh import CHILmesh

FIXTURES = ["annulus", "donut", "block_o", "structured"]


@pytest.mark.parametrize("name", FIXTURES)
def test_triangle_angles_sum_to_180(name):
    mesh = getattr(chilmesh.examples, name)()
    angles = mesh.interior_angles()
    # All shipped fixtures are pure-triangle.
    sums = angles[:, :3].sum(axis=1)
    assert not np.isnan(sums).any(), f"{name}: NaN in interior_angles"
    np.testing.assert_allclose(sums, 180.0, atol=1e-6)


def test_quad_angles_sum_to_360_unit_square():
    pts = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    conn = np.array([[0, 1, 2, 3]])
    mesh = CHILmesh(connectivity=conn, points=pts, grid_name="unit_square")
    angles = mesh.interior_angles()
    np.testing.assert_allclose(angles.sum(axis=1), 360.0, atol=1e-9)
    np.testing.assert_allclose(angles[0], [90.0, 90.0, 90.0, 90.0], atol=1e-9)


def test_quad_angles_no_nan_on_degenerate_quad():
    """B5: a degenerate quad with a zero-length edge must not produce NaN
    angles. Returning NaN bypasses ``elem_quality``'s ``<= 359.99`` zero-out
    because every NaN comparison is False."""
    pts = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    conn = np.array([[0, 1, 2, 3]])
    mesh = CHILmesh(connectivity=conn, points=pts, grid_name="zero_edge_quad")
    angles = mesh.interior_angles()
    assert not np.isnan(angles).any(), f"angles contain NaN: {angles}"

    # And the quality computation must still zero this element out.
    quality, _, _ = mesh.elem_quality()
    assert quality[0] == 0.0, f"elem_quality should reject degenerate quad, got {quality[0]}"
