"""Iterative (MINRES) FEM solver parity tests (#168).

direct_smoother gained an optional low-memory Krylov solver path. The
iterative result must match the direct spsolve result to a tight tolerance on
small meshes; the direct path stays the default (unchanged behavior).
"""
from __future__ import annotations

import numpy as np
import pytest

from chilmesh import CHILmesh


def _pure_quad_mesh() -> CHILmesh:
    points = np.array([[i, j] for j in range(3) for i in range(3)], dtype=float)
    def nid(i, j): return j * 3 + i
    quads = [
        [nid(0, 0), nid(1, 0), nid(1, 1), nid(0, 1)],
        [nid(1, 0), nid(2, 0), nid(2, 1), nid(1, 1)],
        [nid(0, 1), nid(1, 1), nid(1, 2), nid(0, 2)],
        [nid(1, 1), nid(2, 1), nid(2, 2), nid(1, 2)],
    ]
    return CHILmesh(connectivity=np.array(quads, dtype=int), points=points.copy())


def test_iterative_matches_direct_pure_quad():
    mesh = _pure_quad_mesh()
    direct = mesh.direct_smoother(solver="direct")
    iterative = mesh.direct_smoother(solver="iterative")
    np.testing.assert_allclose(iterative[:, :2], direct[:, :2], atol=1e-5)


def test_iterative_matches_direct_on_fixture():
    from chilmesh.examples import donut
    mesh = donut()
    direct = mesh.direct_smoother(solver="direct")
    # kinf boundary pinning inflates ||F||, so MINRES needs a tight relative
    # rtol to drive the absolute solution residual down (#168).
    iterative = mesh.direct_smoother(solver="iterative", tol=1e-12)
    np.testing.assert_allclose(iterative[:, :2], direct[:, :2], atol=1e-4)


def test_minres_alias_accepted():
    mesh = _pure_quad_mesh()
    out = mesh.direct_smoother(solver="minres")
    assert out.shape == mesh.points.shape


def test_unknown_solver_raises():
    mesh = _pure_quad_mesh()
    with pytest.raises(ValueError):
        mesh.direct_smoother(solver="bogus")
