"""Tests for element coloring (spec 005 US2)."""
import numpy as np
import pytest

import chilmesh
from chilmesh import CHILmesh, check_coloring, color_elements

BUNDLED = ["annulus", "donut", "block_o", "structured", "quad_2x2"]


def _neighbors(mesh):
    v2e = mesh.adjacencies["Vert2Elem"]
    out = []
    for e, row in enumerate(mesh.connectivity_list):
        s = set()
        for v in set(int(x) for x in row):
            s |= v2e[v]
        s.discard(e)
        out.append(s)
    return out


def _assert_valid(mesh):
    col = color_elements(mesh)
    nb = _neighbors(mesh)
    for e, s in enumerate(nb):
        assert all(col.colors[m] != col.colors[e] for m in s)
    assert col.n_colors <= max(len(s) for s in nb) + 1
    assert check_coloring(mesh, col)
    return col


def test_single_element():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], float)
    m = CHILmesh(connectivity=np.array([[0, 1, 2]]), points=pts)
    col = color_elements(m)
    assert col.n_colors == 1 and col.colors.tolist() == [0]


def test_mixed_mesh():
    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0]], float)
    m = CHILmesh(connectivity=np.array([[0, 1, 2, 3], [1, 4, 2, 1]]), points=pts)
    col = _assert_valid(m)
    assert col.n_colors == 2


def test_deterministic():
    m = chilmesh.examples.donut()
    assert np.array_equal(color_elements(m).colors, color_elements(m).colors)


def test_check_coloring_detects_bad():
    m = chilmesh.examples.donut()
    assert not check_coloring(m, np.zeros(m.n_elems, dtype=int))
    assert not check_coloring(m, np.zeros(3, dtype=int))


@pytest.mark.parametrize("name", BUNDLED)
def test_all_bundled_meshes(name):  # SC-003
    _assert_valid(getattr(chilmesh.examples, name)())
