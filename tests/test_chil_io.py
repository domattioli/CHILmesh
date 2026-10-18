"""Tests for the .chil container (C1): lossless coordinates and metadata."""
import json
import zipfile

import numpy as np
import pytest

from chilmesh import CHILmesh


def _mixed_mesh():
    pts = np.array(
        [[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1], [1, 2]], dtype=np.float64
    )
    pts = pts * np.pi * 1e5 + 0.1
    pts += 1e-13 * np.arange(pts.size).reshape(pts.shape)
    pts = np.column_stack([pts, np.zeros(len(pts))])
    conn = np.array([[0, 1, 4, 3], [1, 2, 5, 4], [3, 4, 6, 3]])
    return CHILmesh(connectivity=conn, points=pts, grid_name="mix")


def test_chil_lossless_roundtrip(tmp_path):
    m = _mixed_mesh()
    p = tmp_path / "m.chil"
    m.save(str(p))
    r = CHILmesh.load(str(p))
    assert np.array_equal(m.points, r.points)
    assert np.array_equal(m.connectivity_list, r.connectivity_list)
    assert np.array_equal(m.boundary_node_indices(), r.boundary_node_indices())


def test_chil_metadata_roundtrip(tmp_path):
    md = {"arm_or_baseline": "standin", "seeds": [1, 2], "variant": "raw"}
    m = _mixed_mesh()
    p = tmp_path / "m.chil"
    m.save(str(p), metadata=md)
    r = CHILmesh.load(str(p))
    assert r.metadata == md
    assert r.chil_c1_version == 1


def test_chil_metadata_bad_version(tmp_path):
    p = tmp_path / "m.chil"
    _mixed_mesh().save(str(p))
    q = tmp_path / "bad.chil"
    with zipfile.ZipFile(p) as zi, zipfile.ZipFile(q, "w") as zo:
        for n in zi.namelist():
            d = zi.read(n)
            if n == "manifest.json":
                mf = json.loads(d)
                mf["c1_version"] = 999
                d = json.dumps(mf).encode()
            zo.writestr(n, d)
    with pytest.raises(ValueError):
        CHILmesh.load(str(q))


def test_metadata_fort14_unchanged(tmp_path):
    m = _mixed_mesh()
    p = tmp_path / "m.14"
    m.save(str(p))
    assert p.exists()
    with pytest.raises(ValueError):
        m.save(str(p), metadata={"a": 1})
