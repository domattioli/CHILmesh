"""CHILmesh routes the layer peel through the C++ backend when it is selected.

The C++ path must give the same ``layers`` (bit-identical, sorted int64
arrays) and ``n_layers`` as the pure-Python peel. ``CHILMESH_BACKEND=python``
forces the Python peel; boundary-type seeding always uses Python.
"""
from __future__ import annotations

import numpy as np
import pytest

from chilmesh import CHILmesh, examples
from chilmesh.backends.cpp_backend import CPP_AVAILABLE

from _fixture_names import FIXTURE_NAMES

pytestmark = pytest.mark.skipif(not CPP_AVAILABLE, reason="chilmesh_cpp not built")

KEYS = ("OE", "IE", "OV", "IV", "bEdgeIDs")


def _mixed():
    # 3x2 vertex grid: one quad plus two triangles (padded with first vertex).
    pts = np.array([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], dtype=float)
    conn = np.array([[0, 1, 4, 3], [1, 2, 5, 1], [1, 5, 4, 1]])
    return pts, conn


def _source(name):
    if name == "mixed":
        return _mixed()
    mesh = getattr(examples, name)()
    return mesh.points.copy(), mesh.connectivity_list.copy()


def _build(monkeypatch, backend, pts, conn):
    monkeypatch.setenv("CHILMESH_BACKEND", backend)
    return CHILmesh(connectivity=conn.copy(), points=pts.copy())


@pytest.mark.parametrize("name", FIXTURE_NAMES + ["mixed"])
def test_cpp_peel_matches_python(monkeypatch, name):
    pts, conn = _source(name)
    py = _build(monkeypatch, "python", pts, conn)
    cpp = _build(monkeypatch, "cpp", pts, conn)
    assert py._peel_backend == py._adjacency_backend == "python"
    assert cpp._peel_backend == cpp._adjacency_backend == "cpp"
    assert cpp._cpp_init is None  # C++ mesh is dropped after init
    for key in ("Elem2Vert", "Edge2Vert", "Elem2Edge", "Edge2Elem"):
        a, b = np.asarray(py.adjacencies[key]), np.asarray(cpp.adjacencies[key])
        assert a.dtype == b.dtype, key
        np.testing.assert_array_equal(a, b, err_msg=f"{name} {key}")
    for key in ("Vert2Edge", "Vert2Elem"):
        assert py.adjacencies[key] == cpp.adjacencies[key], key
    assert py.adjacencies["EdgeMap"]._map == cpp.adjacencies["EdgeMap"]._map
    assert cpp.n_layers == py.n_layers > 0
    for key in KEYS:
        assert len(cpp.layers[key]) == py.n_layers
        for i in range(py.n_layers):
            a, b = py.layers[key][i], cpp.layers[key][i]
            assert a.dtype == b.dtype, (key, i)
            np.testing.assert_array_equal(a, b, err_msg=f"{name} {key}[{i}]")


def test_non_manifold_mesh_stays_python(monkeypatch):
    # Three triangles share edge (0, 1): the C++ half-edge build must not see it.
    pts = np.array([[0, 0], [1, 0], [0.5, 1], [0.5, -1], [0.5, 2]], dtype=float)
    conn = np.array([[0, 1, 2], [1, 0, 3], [0, 1, 4]])
    monkeypatch.setenv("CHILMESH_BACKEND", "cpp")
    mesh = CHILmesh(connectivity=conn, points=pts, compute_layers=False)
    assert mesh._peel_cpp() is False
    assert mesh._peel_backend == "python"


def test_seeded_peel_stays_python(monkeypatch):
    pts, conn = _source("annulus")
    mesh = _build(monkeypatch, "cpp", pts, conn)
    mesh.boundary_segments = [{"kind": "flow", "ibtype": 0,
                               "nodes": np.unique(mesh.adjacencies["Edge2Vert"][
                                   mesh.adjacencies["Edge2Elem"][:, 1] == -1])}]
    mesh.peel_layers(seed_ibtypes=[0])
    assert mesh._peel_backend == "python"
    assert mesh.n_layers > 0
