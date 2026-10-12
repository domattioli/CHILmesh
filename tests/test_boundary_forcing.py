"""fort.19, fort.20, fort.24 read/write (spec 005 US8, FR-017, FR-018)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent / "fixtures"))
import synthetic_adcirc as syn  # noqa: E402

from chilmesh import (  # noqa: E402
    read_fort19, write_fort19, read_fort20, write_fort20,
    read_fort24, write_fort24, examples,
)


@pytest.fixture
def mesh():
    m = examples.annulus.__wrapped__()
    m.boundary_segments = [
        {"kind": "open", "ibtype": None, "nodes": np.array([0, 1, 2])},
        {"kind": "open", "ibtype": None, "nodes": np.array([3, 4])},
        {"kind": "flow", "ibtype": 22, "nodes": np.array([5, 6])},
        {"kind": "flow", "ibtype": 2, "nodes": np.array([7, 8, 9])},
        {"kind": "flow", "ibtype": 0, "nodes": np.array([10, 11, 12, 13])},   # zero flux: excluded
        {"kind": "flow", "ibtype": 24, "nodes": np.array([14, 15])},           # barrier: excluded
    ]
    return m


# ---- fort.19 ----
def test_fort19_with_mesh_counts_open_nodes(tmp_path, mesh):
    vals = np.arange(15, dtype=float).reshape(3, 5) / 10
    p = syn.write_forcing_synth(tmp_path / "fort.19", 900.0, vals.ravel())
    s = read_fort19(p, mesh=mesh)
    assert s.time_increment == 900.0
    assert s.values.shape == (3, 5)
    np.testing.assert_array_equal(s.values, vals)


def test_fort19_n_nodes_kwarg(tmp_path):
    p = syn.write_forcing_synth(tmp_path / "fort.19", 60.0, np.arange(8.0))
    assert read_fort19(p, n_nodes=4).values.shape == (2, 4)


def test_fort19_without_count_returns_flat_series(tmp_path):
    p = syn.write_forcing_synth(tmp_path / "fort.19", 60.0, np.arange(7.0))
    s = read_fort19(p)
    assert s.time_increment == 60.0 and s.values.shape == (7,)


def test_fort19_count_mismatch_names_both(tmp_path, mesh):
    p = syn.write_forcing_synth(tmp_path / "fort.19", 900.0, np.arange(12.0))
    with pytest.raises(ValueError) as ei:
        read_fort19(p, mesh=mesh)
    assert "12" in str(ei.value) and "5" in str(ei.value)
    with pytest.raises(ValueError) as ei:
        read_fort19(p, mesh=mesh, n_nodes=4)
    assert "4" in str(ei.value) and "5" in str(ei.value)


def test_fort19_round_trip(tmp_path, mesh):
    vals = np.random.default_rng(1).normal(size=(4, 5))
    out = tmp_path / "a"
    write_fort19(read_fort19(syn.write_forcing_synth(tmp_path / "in", 3600.0, vals.ravel()), mesh=mesh), out)
    s = read_fort19(out, mesh=mesh)
    np.testing.assert_array_equal(s.values, vals)
    assert s.time_increment == 3600.0


# ---- fort.20 ----
def test_fort20_matching_nodes_use_only_flux_types(tmp_path, mesh):
    # IBTYPE 2 and 22 segments: 2 + 3 = 5 nodes; types 0 and 24 excluded
    vals = np.arange(10, dtype=float).reshape(2, 5)
    p = syn.write_forcing_synth(tmp_path / "fort.20", 3600.0, vals.ravel())
    s = read_fort20(p, mesh=mesh)
    assert s.values.shape == (2, 5)


def test_fort20_ignores_trailing_comments_and_parses_without_mesh(tmp_path):
    p = syn.write_forcing_synth(tmp_path / "fort.20", 31536000.0, [26.0, 30.3, -1.5], comments=True)
    s = read_fort20(p)
    assert s.time_increment == 31536000.0
    np.testing.assert_array_equal(s.values, [26.0, 30.3, -1.5])


def test_fort20_count_mismatch_names_both(tmp_path, mesh):
    p = syn.write_forcing_synth(tmp_path / "fort.20", 1.0, np.arange(7.0))
    with pytest.raises(ValueError) as ei:
        read_fort20(p, mesh=mesh)
    assert "7" in str(ei.value) and "5" in str(ei.value)


def test_fort20_round_trip(tmp_path, mesh):
    vals = np.random.default_rng(2).normal(size=(3, 5))
    p = syn.write_forcing_synth(tmp_path / "in", 60.0, vals.ravel())
    out = tmp_path / "out"
    write_fort20(read_fort20(p, mesh=mesh), out)
    np.testing.assert_array_equal(read_fort20(out, mesh=mesh).values, vals)


def test_mesh_without_matching_boundary_nodes_errors(tmp_path):
    m = examples.annulus.__wrapped__()
    m.boundary_segments = []
    p = syn.write_forcing_synth(tmp_path / "f", 1.0, [1.0, 2.0])
    with pytest.raises(ValueError, match="0"):
        read_fort19(p, mesh=m)


# ---- fort.24 ----
def _consts(n):
    rng = np.random.default_rng(3)
    return [("M2", 0.000140518902761, rng.random(n).round(8), (rng.random(n) * 360).round(6)),
            ("S2", 0.000145444104333, rng.random(n).round(8), (rng.random(n) * 360).round(6))]


def test_fort24_fields(tmp_path):
    cs = _consts(6)
    p = syn.write_fort24_synth(tmp_path / "fort.24", cs)
    got = read_fort24(p)
    assert [c.name for c in got] == ["M2", "S2"]
    assert got[0].frequency == pytest.approx(0.000140518902761)
    assert got[0].nodal_factor == 1.0
    np.testing.assert_allclose(got[1].amp, cs[1][2])
    np.testing.assert_allclose(got[1].phase, cs[1][3])
    assert got[0].amp.shape == (6,)


def test_fort24_mesh_count(tmp_path, mesh):
    p = syn.write_fort24_synth(tmp_path / "fort.24", _consts(mesh.n_verts))
    assert read_fort24(p, mesh=mesh)[0].amp.size == mesh.n_verts
    p2 = syn.write_fort24_synth(tmp_path / "bad", _consts(6))
    with pytest.raises(ValueError) as ei:
        read_fort24(p2, mesh=mesh)
    assert "6" in str(ei.value) and str(mesh.n_verts) in str(ei.value)


def test_fort24_round_trip(tmp_path):
    p = syn.write_fort24_synth(tmp_path / "in", _consts(9))
    first = read_fort24(p)
    out = tmp_path / "out"
    write_fort24(first, out)
    again = read_fort24(out)
    assert len(again) == 2
    for a, b in zip(first, again):
        assert (a.name, a.alpha, a.frequency, a.nodal_factor) == (b.name, b.alpha, b.frequency, b.nodal_factor)
        np.testing.assert_array_equal(a.amp, b.amp)
        np.testing.assert_array_equal(a.phase, b.phase)
