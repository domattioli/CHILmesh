"""``_elem_type`` regression tests (B3).

In MATLAB, the value ``0`` was a placeholder for "missing 4th vertex" in
mixed-element connectivity. In this 0-indexed Python port, vertex id ``0``
is a perfectly valid vertex — but the original ``_elem_type`` reused the
MATLAB sentinel and silently demoted any quad whose 4th vertex happens
to be ``0`` to a triangle. This file pins the correct behaviour.
"""
from __future__ import annotations

import numpy as np

from chilmesh import CHILmesh


def _two_quads_with_vertex_zero():
    # Vertex 0 deliberately appears in the *4th* slot of one quad to trigger
    # the historical bug. Two unit squares stacked vertically:
    pts = np.array(
        [
            [0.0, 0.0],  # 0
            [1.0, 0.0],  # 1
            [1.0, 1.0],  # 2
            [0.0, 1.0],  # 3
            [1.0, 2.0],  # 4
            [0.0, 2.0],  # 5
        ]
    )
    # Second quad's 4th vertex == 0  (rolled cyclically: 1,4,5,0->1)
    # We use [1, 4, 5, 0] which is a valid CCW square traversal.
    conn = np.array([[0, 1, 2, 3], [1, 4, 5, 0]])
    return pts, conn


def test_elem_type_classifies_vertex_zero_quad_as_quad():
    pts, conn = _two_quads_with_vertex_zero()
    mesh = CHILmesh(connectivity=conn, points=pts, grid_name="vertex0_quads")

    tri_elems, quad_elems = mesh._elem_type()
    assert len(tri_elems) == 0, (
        f"Expected zero triangles, got {len(tri_elems)} (B3 still present)"
    )
    assert sorted(quad_elems.tolist()) == [0, 1]
    assert mesh.type == "Quadrilateral"


def test_elem_type_real_triangle_in_padded_mesh():
    """A genuine triangle padded with a repeated vertex must still be tri."""
    pts = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
        ]
    )
    # First element: triangle (0,1,2) padded by repeating vertex 0.
    # Second element: real quad (1,3,4,5).
    conn = np.array([[0, 1, 2, 0], [1, 3, 4, 5]])
    mesh = CHILmesh(connectivity=conn, points=pts, grid_name="mixed")
    tri_elems, quad_elems = mesh._elem_type()
    assert tri_elems.tolist() == [0]
    assert quad_elems.tolist() == [1]
    assert mesh.type == "Mixed-Element"
