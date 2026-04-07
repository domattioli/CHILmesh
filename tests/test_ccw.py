"""``_ensure_ccw_orientation`` regression tests (B4 + idempotency).

The legacy implementation switched on ``connectivity_list.shape[1]`` and
applied the quad permutation ``[0, 3, 2, 1]`` to padded-triangle rows in a
mixed-element mesh, scrambling them. The test below builds such a mesh
with deliberately CW elements and pins the corrected behaviour.
"""
from __future__ import annotations

import numpy as np
import pytest

import chilmesh
from chilmesh import CHILmesh

FIXTURES = ["annulus", "donut", "block_o", "structured"]


@pytest.mark.parametrize("name", FIXTURES)
def test_ccw_idempotent(name):
    mesh = getattr(chilmesh.examples, name)()
    before = mesh.connectivity_list.copy()
    mesh._ensure_ccw_orientation()
    after = mesh.connectivity_list.copy()
    np.testing.assert_array_equal(before, after)


def test_ccw_flips_clockwise_triangle():
    pts = np.array([[0.0, 0.0], [0.0, 1.0], [1.0, 0.0]])  # CW order
    conn = np.array([[0, 1, 2]])
    mesh = CHILmesh(connectivity=conn, points=pts, grid_name="cw_tri")
    assert mesh.signed_area().min() > 0  # constructor flipped it


def test_ccw_flips_clockwise_quad():
    pts = np.array([[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]])  # CW order
    conn = np.array([[0, 1, 2, 3]])
    mesh = CHILmesh(connectivity=conn, points=pts, grid_name="cw_quad")
    assert mesh.signed_area().min() > 0


def test_ccw_mixed_element_does_not_scramble_triangle_row():
    """B4: padded-triangle row in a 4-column mixed mesh must not be flipped
    with the quad permutation."""
    pts = np.array(
        [
            [0.0, 0.0],  # 0
            [1.0, 0.0],  # 1
            [0.0, 1.0],  # 2  triangle (0,1,2) padded with 0
            [2.0, 0.0],  # 3
            [2.0, 1.0],  # 4
            [1.0, 1.0],  # 5
        ]
    )
    # Triangle (0,1,2) is CCW; pad with vertex 0. Quad (1,3,4,5) is CCW.
    conn = np.array([[0, 1, 2, 0], [1, 3, 4, 5]])
    mesh = CHILmesh(connectivity=conn, points=pts, grid_name="mixed_ccw")

    # The padded triangle's first three vertices must remain {0,1,2}; the
    # quad-permutation [0,3,2,1] would have produced [0,0,2,1].
    tri_row = mesh.connectivity_list[0]
    assert set(tri_row[:3].tolist()) == {0, 1, 2}
    assert mesh.signed_area().min() > 0


def test_ccw_mixed_element_flips_padded_cw_triangle_correctly():
    """B4: a CW padded-triangle row in a 4-column mesh must be flipped using
    the *triangle* permutation [0,2,1,3], not the quad one [0,3,2,1] which
    would produce a row containing the pad twice."""
    pts = np.array(
        [
            [0.0, 0.0],  # 0
            [0.0, 1.0],  # 1
            [1.0, 0.0],  # 2  -- triangle (0,1,2) is CW
            [2.0, 0.0],  # 3
            [3.0, 0.0],  # 4
            [3.0, 1.0],  # 5
            [2.0, 1.0],  # 6
        ]
    )
    conn = np.array([[0, 1, 2, 0], [3, 4, 5, 6]])
    mesh = CHILmesh(connectivity=conn, points=pts, grid_name="cw_padded_tri")

    areas = mesh.signed_area()
    assert (areas > 0).all(), f"areas={areas}"

    # First three slots of the (now-flipped) triangle row must still be a
    # permutation of {0,1,2} — never duplicate the pad value.
    tri_row = mesh.connectivity_list[0]
    assert set(tri_row[:3].tolist()) == {0, 1, 2}, (
        f"B4 still scrambling padded triangle: row={tri_row.tolist()}"
    )


def test_ccw_mixed_element_flips_only_cw_quad():
    """A CW quad in a mixed mesh must be flipped, while the padded triangle
    next to it stays untouched."""
    pts = np.array(
        [
            [0.0, 0.0],  # 0
            [1.0, 0.0],  # 1
            [0.0, 1.0],  # 2  CCW triangle padded with 0
            [2.0, 0.0],  # 3
            [2.0, 1.0],  # 4
            [3.0, 1.0],  # 5
        ]
    )
    # Quad (1,3,5,4) traverses clockwise on purpose.
    conn = np.array([[0, 1, 2, 0], [1, 3, 5, 4]])
    mesh = CHILmesh(connectivity=conn, points=pts, grid_name="mixed_cw_quad")
    areas = mesh.signed_area()
    assert (areas > 0).all()
    # Padded triangle row's first 3 vertices unchanged.
    assert set(mesh.connectivity_list[0, :3].tolist()) == {0, 1, 2}
