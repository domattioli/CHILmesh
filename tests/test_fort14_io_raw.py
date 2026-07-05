"""Faithful fort.14 reader/writer (#214): preserves id columns, winding, true
tri/quad arity, and barrier/weir boundary columns — the gaps in
CHILmesh.read_from_fort14 that break Valence's byte-parity."""
import textwrap

import pytest

from chilmesh.fort14_io import (
    Fort14ParseError,
    read_fort14_raw,
    write_fort14_raw,
)


def _write(tmp_path, text, name="mesh.14"):
    p = tmp_path / name
    p.write_text(textwrap.dedent(text).strip() + "\n", encoding="utf-8")
    return p


CW = """
    cw mesh
    1 3
    1 0.0 0.0 0.0
    2 0.0 1.0 0.0
    3 1.0 0.0 0.0
    1 3 1 2 3
    0
    0
    0
    0
    """

MIXED = """
    mixed mesh
    2 5
    10 0.0 0.0 0.0
    20 1.0 0.0 0.0
    30 1.0 1.0 5.0
    40 0.0 1.0 0.0
    50 2.0 0.5 0.0
    5 3 10 20 30
    6 4 10 30 40 50
    0
    0
    0
    0
    """

BARRIER = """
    barrier mesh
    1 4
    1 0.0 0.0 0.0
    2 1.0 0.0 0.0
    3 1.0 1.0 0.0
    4 0.0 1.0 0.0
    1 3 1 2 3
    1
    2
    2 1
    1
    2
    1
    2
    2 24
    1 3 2.5 1.0 1.0
    2 4 2.5 1.0 1.0
    """

LEGACY = """
    legacy mesh
    1 3
    1 0.0 0.0 0.0
    2 1.0 0.0 0.0
    3 0.0 1.0 0.0
    1 3 1 2 3
    """


def test_winding_preserved(tmp_path):
    # A CW element (signed area < 0) is kept in file order — NOT reoriented.
    raw = read_fort14_raw(_write(tmp_path, CW))
    assert raw.elements[1] == (1, 2, 3)


def test_mixed_arity_and_ids_preserved(tmp_path):
    raw = read_fort14_raw(_write(tmp_path, MIXED))
    assert raw.node_ids == [10, 20, 30, 40, 50]        # non-contiguous ids kept
    assert raw.coords[30] == (1.0, 1.0, 5.0)           # depth preserved
    assert raw.elem_ids == [5, 6]                      # element id column kept
    assert raw.elements[5] == (10, 20, 30)             # true tri arity (3), no pad
    assert raw.elements[6] == (10, 30, 40, 50)         # true quad arity (4)


def test_barrier_weir_columns_parsed(tmp_path):
    raw = read_fort14_raw(_write(tmp_path, BARRIER))
    assert len(raw.open_boundaries) == 1
    assert raw.open_boundaries[0].ibtype == 1
    assert raw.open_boundaries[0].nodes == [1, 2]
    assert len(raw.flow_boundaries) == 1
    fb = raw.flow_boundaries[0]
    assert fb.ibtype == 24
    assert fb.nodes == [1, 2]
    assert fb.back_nodes == [3, 4]
    assert fb.heights == [2.5, 2.5]
    assert fb.coeffs == [[1.0, 1.0], [1.0, 1.0]]


def test_legacy_no_boundary_block(tmp_path):
    raw = read_fort14_raw(_write(tmp_path, LEGACY))
    assert raw.open_boundaries == []
    assert raw.flow_boundaries == []


def test_roundtrip_structural_equality(tmp_path):
    r1 = read_fort14_raw(_write(tmp_path, BARRIER))
    dst = tmp_path / "rt.14"
    assert write_fort14_raw(r1, dst) is True
    r2 = read_fort14_raw(dst)
    assert r1.node_ids == r2.node_ids
    assert r1.elem_ids == r2.elem_ids
    assert r1.elements == r2.elements
    assert r1.coords == r2.coords
    assert [b.nodes for b in r1.flow_boundaries] == [b.nodes for b in r2.flow_boundaries]
    assert [b.back_nodes for b in r1.flow_boundaries] == [b.back_nodes for b in r2.flow_boundaries]
    assert [b.heights for b in r1.flow_boundaries] == [b.heights for b in r2.flow_boundaries]


def test_malformed_boundary_raises(tmp_path):
    bad = """
        bad mesh
        1 3
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 0.0 1.0 0.0
        1 3 1 2 3
        1
        2
        2 1
        not-a-node
        2
        0
        0
        """
    with pytest.raises(Fort14ParseError):
        read_fort14_raw(_write(tmp_path, bad))
