"""Square-preservation regression tests for the smoothers (#105).

Operator dispute on #105: "the quads of that mixed test case should be perfect
squares." These tests pin that expectation as an executable check.

Ground truth established empirically (deterministic across runs):

| mesh           | smoother              | interior node (start [1,1]) |
|----------------|-----------------------|-----------------------------|
| pure quad      | direct_smoother       | [1.0000, 1.0000]  exact     |
| mixed tri/quad | direct_smoother       | [0.9269, 0.9154]  distorted |
| mixed tri/quad | angle_based_smoother  | [0.9879, 0.9879]  ~0.017 off|

The pure-quad case is exact because every interior quad angle is 90deg and
cot(90deg)=0, so the FEM right-hand side vanishes and the symmetric quad
stiffness leaves the node fixed.

The mixed case is NOT preserved by ``direct_smoother``: once triangle
stiffness couples into a shared interior node, the *linear* FEM equilibrium
shifts off the square position and no RHS force recovers it. That shift is the
subject of #105 and is captured here as a strict xfail so that any future fix
which restores square preservation will trip this test and force the marker's
removal.
"""
from __future__ import annotations

import numpy as np
import pytest

from chilmesh import CHILmesh


def _grid_3x3_quads() -> tuple[np.ndarray, list[list[int]]]:
    """3x3 node grid -> four unit quads; interior node id 4 sits at (1, 1)."""
    points = np.array([[i, j] for j in range(3) for i in range(3)], dtype=float)

    def nid(i: int, j: int) -> int:
        return j * 3 + i

    quads = [
        [nid(0, 0), nid(1, 0), nid(1, 1), nid(0, 1)],
        [nid(1, 0), nid(2, 0), nid(2, 1), nid(1, 1)],
        [nid(0, 1), nid(1, 1), nid(1, 2), nid(0, 2)],
        [nid(1, 1), nid(2, 1), nid(2, 2), nid(1, 2)],
    ]
    return points, quads


def _pure_quad_mesh() -> CHILmesh:
    points, quads = _grid_3x3_quads()
    return CHILmesh(connectivity=np.array(quads, dtype=int), points=points.copy())


def _mixed_mesh() -> CHILmesh:
    """Same grid, but the bottom-left quad is split into two padded triangles."""
    points, quads = _grid_3x3_quads()
    q0 = quads[0]
    conn = [
        [q0[0], q0[1], q0[2], q0[0]],  # padded triangle
        [q0[0], q0[2], q0[3], q0[0]],  # padded triangle
    ] + [list(q) for q in quads[1:]]
    return CHILmesh(connectivity=np.array(conn, dtype=int), points=points.copy())


INTERIOR = 4
SQUARE_POS = np.array([1.0, 1.0])


def test_direct_smoother_preserves_pure_quad_square():
    """On a pure-quad mesh the interior node must not move (exact)."""
    mesh = _pure_quad_mesh()
    assert mesh.type == "Quadrilateral"
    out = mesh.direct_smoother()
    np.testing.assert_allclose(out[INTERIOR, :2], SQUARE_POS, atol=1e-9)


def test_angle_based_smoother_keeps_mixed_node_near_square():
    """angle_based_smoother is the quality-gated path: stays within 0.02 of square."""
    mesh = _mixed_mesh()
    assert mesh.type == "Mixed-Element"
    out = mesh.angle_based_smoother(n_iter=50)
    disp = float(np.linalg.norm(out[INTERIOR, :2] - SQUARE_POS))
    assert disp < 0.02, f"interior node drifted {disp:.4f} from square"


@pytest.mark.xfail(
    strict=True,
    reason="#105: direct_smoother (linear FEM solve) does not preserve squares "
    "in mixed tri/quad meshes; interior node drifts ~0.11. Remove this marker "
    "when a square-preserving mixed-mesh smoother lands.",
)
def test_direct_smoother_preserves_mixed_quad_square():
    """Desired: a perfect square stays square even when a neighbour is triangulated."""
    mesh = _mixed_mesh()
    out = mesh.direct_smoother()
    disp = float(np.linalg.norm(out[INTERIOR, :2] - SQUARE_POS))
    assert disp < 1e-3, f"interior node drifted {disp:.4f} from square"
