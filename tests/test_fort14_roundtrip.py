"""Roundtrip ``.fort.14`` IO regression test (B1).

Loads each shipped fixture, writes it to a temp directory via the
instance method ``CHILmesh.write_to_fort14`` (the historically-broken
recursive entry point), reloads, and asserts that vertex coordinates
and connectivity survive verbatim.

Catches B1 (recursive ``write_to_fort14``).
"""
from __future__ import annotations

import textwrap

import numpy as np

import pytest

import chilmesh
from chilmesh import CHILmesh
from conftest import TRI_FIXTURE_NAMES as FIXTURES


@pytest.mark.parametrize("name", FIXTURES)
def test_fort14_roundtrip_identity(name, tmp_path):
    mesh = getattr(chilmesh.examples, name)()

    out = tmp_path / f"{name}.fort.14"
    ok = mesh.write_to_fort14(str(out), grid_name=mesh.grid_name or name)
    assert ok is True, f"{name}: write_to_fort14 returned a falsy value"
    assert out.exists() and out.stat().st_size > 0, (
        f"{name}: fort14 output missing or empty at {out}"
    )

    reloaded = CHILmesh.read_from_fort14(out)

    assert reloaded.n_verts == mesh.n_verts, (
        f"{name}: n_verts drift {mesh.n_verts} -> {reloaded.n_verts}"
    )
    assert reloaded.n_elems == mesh.n_elems, (
        f"{name}: n_elems drift {mesh.n_elems} -> {reloaded.n_elems}"
    )
    np.testing.assert_allclose(reloaded.points[:, :2], mesh.points[:, :2], atol=1e-7)

    # Connectivity is allowed to differ in vertex *order within an element*
    # because ``_ensure_ccw_orientation`` may flip CW elements on reload, but
    # the *set* of vertices defining each element must match.
    orig_sets = [frozenset(row[:3]) for row in mesh.connectivity_list]
    rel_sets = [frozenset(row[:3]) for row in reloaded.connectivity_list]
    assert sorted(orig_sets) == sorted(rel_sets), (
        f"{name}: connectivity vertex-sets differ after roundtrip"
    )


def test_chilmesh_boundaries_present_distinguishes_empty_from_absent(tmp_path):
    """#259: CHILmesh.read_from_fort14 distinguishes present-but-empty from absent."""
    # Minimal 1-element triangle mesh WITH physical 0/0/0/0 boundary block.
    with_block = """
        test mesh with boundaries
        1 3
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 0.0 1.0 0.0
        1 3 1 2 3
        0
        0
        0
        0
        """
    p = tmp_path / "with_block.14"
    p.write_text(textwrap.dedent(with_block).strip() + "\n", encoding="utf-8")
    mesh_with = CHILmesh.read_from_fort14(p, compute_layers=False, compute_adjacencies=False)
    assert mesh_with.boundaries_present is True

    # Same mesh WITHOUT any boundary block at all.
    without_block = """
        test mesh no boundaries
        1 3
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 0.0 1.0 0.0
        1 3 1 2 3
        """
    p = tmp_path / "without_block.14"
    p.write_text(textwrap.dedent(without_block).strip() + "\n", encoding="utf-8")
    mesh_without = CHILmesh.read_from_fort14(p, compute_layers=False, compute_adjacencies=False)
    assert mesh_without.boundaries_present is False
