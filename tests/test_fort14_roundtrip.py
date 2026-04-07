"""Roundtrip ``.fort.14`` IO regression test (B1).

Loads each shipped fixture, writes it to a temp directory via the
instance method ``CHILmesh.write_to_fort14`` (the historically-broken
recursive entry point), reloads, and asserts that vertex coordinates
and connectivity survive verbatim.

Catches B1 (recursive ``write_to_fort14``).
"""
from __future__ import annotations

import numpy as np

import chilmesh
from chilmesh import CHILmesh

FIXTURES = ["annulus", "donut", "block_o", "structured"]


import pytest


@pytest.mark.parametrize("name", FIXTURES)
def test_fort14_roundtrip_identity(name, tmp_path):
    mesh = getattr(chilmesh.examples, name)()

    out = tmp_path / f"{name}.fort.14"
    ok = mesh.write_to_fort14(str(out), grid_name=mesh.grid_name or name)
    assert ok is True, "write_to_fort14 returned a falsy value"
    assert out.exists() and out.stat().st_size > 0

    reloaded = CHILmesh.read_from_fort14(out)

    assert reloaded.n_verts == mesh.n_verts
    assert reloaded.n_elems == mesh.n_elems
    np.testing.assert_allclose(reloaded.points[:, :2], mesh.points[:, :2], atol=1e-7)

    # Connectivity is allowed to differ in vertex *order within an element*
    # because ``_ensure_ccw_orientation`` may flip CW elements on reload, but
    # the *set* of vertices defining each element must match.
    orig_sets = [frozenset(row[:3]) for row in mesh.connectivity_list]
    rel_sets = [frozenset(row[:3]) for row in reloaded.connectivity_list]
    assert sorted(orig_sets) == sorted(rel_sets)
