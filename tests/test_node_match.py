"""Tests for cross-mesh node matching + nodal-field deltas (#239)."""
from __future__ import annotations

import numpy as np

from chilmesh import match_nodes, derive_tolerance, nodal_field_delta


def _grid(n=4, spacing=1.0):
    xs, ys = np.meshgrid(np.arange(n) * spacing, np.arange(n) * spacing)
    return np.column_stack([xs.ravel(), ys.ravel()]).astype(float)


def test_identity_match_same_points():
    pts = _grid()
    m = match_nodes(pts, pts.copy(), tolerance=0.1)
    assert m.method == "identity"
    assert m.matched == {i: i for i in range(len(pts))}
    assert m.added == [] and m.removed == []


def test_spatial_match_survives_renumbering():
    pts = _grid()
    perm = np.random.default_rng(0).permutation(len(pts))
    b = pts[perm]
    m = match_nodes(pts, b, tolerance=0.4)
    assert m.method == "spatial"
    assert len(m.matched) == len(pts)
    assert m.added == [] and m.removed == []
    # every a index maps to the b row holding the same coordinate
    for aid, bid in m.matched.items():
        assert np.allclose(pts[aid], b[bid])


def test_added_and_removed():
    a = _grid()  # 16 pts
    # b drops the last a-point's location and adds one far away
    b = np.vstack([a[:-1], [[100.0, 100.0]]])
    m = match_nodes(a, b, tolerance=0.4)
    assert m.removed == [15]           # a's last point has no b within tol
    assert m.added == [15]             # b's far point matched nothing


def test_nodal_field_delta_bathymetry():
    a = _grid()
    b = a.copy()
    fa = np.zeros(len(a))
    fb = fa.copy()
    fb[0] = 2.0    # deepened
    fb[1] = -1.0   # shallowed
    m = match_nodes(a, b, tolerance=0.1)
    d = nodal_field_delta(m, fa, fb)
    assert d["matched_nodes"] == len(a)
    assert d["changed_nodes"] == 2
    assert d["delta"][0] == 2.0
    assert d["delta"][1] == -1.0
    assert d["increased"]["count"] == 1
    assert d["increased"]["max"] == 2.0
    assert d["decreased"]["count"] == 1
    assert d["decreased"]["max"] == 1.0


def test_derive_tolerance_positive():
    a = _grid(spacing=2.0)
    b = _grid(spacing=1.0)
    tol = derive_tolerance(a, b)
    # finer mesh spacing = 1.0, fraction 0.5 => ~0.5
    assert 0.0 < tol <= 1.0
