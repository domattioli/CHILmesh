"""Tests for MutableMesh.split_edge — valence-aware even (2->4) edge split (#237)."""

from __future__ import annotations

import numpy as np
import pytest
from chilmesh import CHILmesh, MutableMesh
from chilmesh import examples


@pytest.fixture(params=["annulus", "donut"])
def triangle_mesh(request):
    fn = getattr(examples, request.param)
    return fn.__wrapped__() if hasattr(fn, "__wrapped__") else fn()


def _first_interior_edge(mesh):
    e2e = mesh.edge2elem()
    for eid in range(mesh.n_edges):
        a, b = int(e2e[eid, 0]), int(e2e[eid, 1])
        if a != -1 and b != -1:
            verts = mesh.edge2vert(np.array([eid]))[0]
            return eid, int(verts[0]), int(verts[1])
    raise AssertionError("no interior edge found")


def _valence(mesh, vid):
    e2v = mesh.edge2vert()
    return int(np.sum(np.any(e2v == vid, axis=1)))


def _positive_area(mesh, eid):
    tri = mesh.connectivity_list[int(eid), :3]
    p0, p1, p2 = mesh.points[tri[0], :2], mesh.points[tri[1], :2], mesh.points[tri[2], :2]
    return 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]))


def test_split_edge_2_to_4(triangle_mesh):
    mutable = MutableMesh(triangle_mesh)
    n_elems0, n_verts0 = triangle_mesh.n_elems, triangle_mesh.n_verts
    _, u, v = _first_interior_edge(triangle_mesh)

    new_ids = mutable.split_edge(u, v)

    assert len(new_ids) == 4
    assert triangle_mesh.n_elems == n_elems0 + 2   # 2 -> 4
    assert triangle_mesh.n_verts == n_verts0 + 1


def test_split_edge_midpoint_location(triangle_mesh):
    mutable = MutableMesh(triangle_mesh)
    _, u, v = _first_interior_edge(triangle_mesh)
    pu = triangle_mesh.points[u, :2].copy()
    pv = triangle_mesh.points[v, :2].copy()

    mutable.split_edge(u, v)

    m = triangle_mesh.n_verts - 1
    np.testing.assert_allclose(triangle_mesh.points[m, :2], (pu + pv) / 2.0, rtol=1e-10)


def test_split_edge_new_vertex_valence_4(triangle_mesh):
    mutable = MutableMesh(triangle_mesh)
    _, u, v = _first_interior_edge(triangle_mesh)

    mutable.split_edge(u, v)

    m = triangle_mesh.n_verts - 1
    assert _valence(triangle_mesh, m) == 4


def test_split_edge_watertight_positive_area(triangle_mesh):
    mutable = MutableMesh(triangle_mesh)
    _, u, v = _first_interior_edge(triangle_mesh)

    new_ids = mutable.split_edge(u, v)

    for eid in new_ids:
        assert _positive_area(triangle_mesh, eid) > 1e-12, f"element {eid} non-positive area"


def test_split_edge_at_param(triangle_mesh):
    mutable = MutableMesh(triangle_mesh)
    _, u, v = _first_interior_edge(triangle_mesh)
    pu = triangle_mesh.points[u, :2].copy()
    pv = triangle_mesh.points[v, :2].copy()

    mutable.split_edge(u, v, at=0.25)

    m = triangle_mesh.n_verts - 1
    np.testing.assert_allclose(triangle_mesh.points[m, :2], pu + 0.25 * (pv - pu), rtol=1e-10)


def test_split_edge_boundary_raises(triangle_mesh):
    mutable = MutableMesh(triangle_mesh)
    e2e = triangle_mesh.edge2elem()
    bid = None
    for eid in range(triangle_mesh.n_edges):
        a, b = int(e2e[eid, 0]), int(e2e[eid, 1])
        if a == -1 or b == -1:
            bid = eid
            break
    if bid is None:
        pytest.skip("mesh has no boundary edge")
    verts = triangle_mesh.edge2vert(np.array([bid]))[0]
    with pytest.raises(ValueError):
        mutable.split_edge(int(verts[0]), int(verts[1]))


def test_split_edge_nonexistent_edge_raises(triangle_mesh):
    mutable = MutableMesh(triangle_mesh)
    e2v = triangle_mesh.edge2vert()
    existing = {frozenset((int(a), int(b))) for a, b in e2v}
    pair = None
    nv = triangle_mesh.n_verts
    for i in range(nv):
        for j in range(i + 1, nv):
            if frozenset((i, j)) not in existing:
                pair = (i, j)
                break
        if pair:
            break
    assert pair is not None
    with pytest.raises(ValueError):
        mutable.split_edge(*pair)


def test_split_edge_bad_at_raises(triangle_mesh):
    mutable = MutableMesh(triangle_mesh)
    _, u, v = _first_interior_edge(triangle_mesh)
    with pytest.raises(ValueError):
        mutable.split_edge(u, v, at=1.5)
