"""Tests for .2dm NS (nodestring) write/read (spec 005 US1)."""
import numpy as np

import chilmesh
from chilmesh import CHILmesh


def _seg(nodes):
    return {"kind": "nodestring", "ibtype": None, "nodes": np.array(nodes, dtype=int)}


def _tri_mesh():
    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], float)
    return CHILmesh(connectivity=np.array([[0, 1, 2], [0, 2, 3]]), points=pts,
                    compute_layers=False)


def _roundtrip(mesh, tmp_path):
    f = tmp_path / "m.2dm"
    mesh.write_to_2dm(str(f))
    return f, CHILmesh.read_from_2dm(f, compute_layers=False)


def test_tri_nodestrings_roundtrip(tmp_path):
    m = _tri_mesh()
    m.boundary_segments = [_seg([0, 1, 2]), _seg([2, 3, 0])]
    f, r = _roundtrip(m, tmp_path)
    assert [s["nodes"].tolist() for s in r.boundary_segments] == [[0, 1, 2], [2, 3, 0]]
    assert "NS 1 2 -3" in f.read_text()


def test_long_string_continuation(tmp_path):
    m = chilmesh.examples.annulus.__wrapped__()
    nodes = list(range(25))
    m.boundary_segments = [_seg(nodes)]
    _, r = _roundtrip(m, tmp_path)
    assert r.boundary_segments[0]["nodes"].tolist() == nodes


def test_quad_only(tmp_path):
    m = chilmesh.examples.quad_2x2.__wrapped__()
    m.boundary_segments = [_seg([0, 1, 2])]
    _, r = _roundtrip(m, tmp_path)
    assert r.boundary_segments[0]["nodes"].tolist() == [0, 1, 2]


def test_mixed(tmp_path):
    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0]], float)
    conn = np.array([[0, 1, 2, 3], [1, 4, 2, 1]])
    m = CHILmesh(connectivity=conn, points=pts, compute_layers=False)
    m.boundary_segments = [_seg([0, 1, 4])]
    _, r = _roundtrip(m, tmp_path)
    assert r.boundary_segments[0]["nodes"].tolist() == [0, 1, 4]


def test_no_boundary_writes_no_ns(tmp_path):
    m = _tri_mesh()
    m.boundary_segments = []
    f, r = _roundtrip(m, tmp_path)
    assert "NS" not in f.read_text()
    assert r.boundary_segments == []


def test_nodestrings_false(tmp_path):
    m = _tri_mesh()
    m.boundary_segments = [_seg([0, 1])]
    f = tmp_path / "x.2dm"
    m.write_to_2dm(str(f), nodestrings=False)
    assert "NS" not in f.read_text()


def test_alias_still_works(tmp_path):
    m = _tri_mesh()
    f = tmp_path / "a.2dm"
    m._write_2dm(str(f))
    assert CHILmesh.read_from_2dm(f, compute_layers=False).n_elems == 2


def test_unknown_sms_cards_ignored(tmp_path):
    f = tmp_path / "s.2dm"
    f.write_text("MESH2D\nMESHNAME \"x\"\nNUM_MATERIALS_PER_ELEM 1\nE3T 1 1 2 3 1\n"
                 "ND 1 0 0 0\nND 2 1 0 0\nND 3 0 1 0\nBEG\nGM \"m\"\nBCP 1 2\nEND\n"
                 "NS 1 2 -3\n")
    r = CHILmesh.read_from_2dm(f, compute_layers=False)
    assert r.n_elems == 1 and r.boundary_segments[0]["nodes"].tolist() == [0, 1, 2]
