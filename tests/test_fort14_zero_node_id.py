"""Tests for zero/invalid node id handling in fort.14 parsing (#282).

CHILmesh.read_from_fort14 should raise ValueError when encountering
node ids < 1 or > n_verts, instead of silently wrapping to vertex -1.
"""
import textwrap

import numpy as np
import pytest

from chilmesh import CHILmesh
from chilmesh.fort14_io import read_fort14_raw


def _write(tmp_path, text, name="mesh.14"):
    """Helper to write a fort.14 fixture to tmp_path."""
    p = tmp_path / name
    p.write_text(textwrap.dedent(text).strip() + "\n", encoding="utf-8")
    return p


def test_zero_node_id_raises(tmp_path):
    """Case a: quad + 0-padded triangle (element 2 has node id 0)."""
    fixture = """
        pad0 repro
        2 5
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        5 2.0 0.0 0.0
        1 4 1 2 3 4
        2 4 2 5 3 0
        0
        0
        0
        0
        """
    p = _write(tmp_path, fixture)
    with pytest.raises(ValueError) as exc_info:
        CHILmesh.read_from_fort14(p)

    err = str(exc_info.value)
    assert "element 2" in err
    assert "0" in err


def test_two_zeros_raises(tmp_path):
    """Case b: row with two zero node ids."""
    fixture = """
        two zeros
        2 5
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        5 2.0 0.0 0.0
        1 3 1 2 3
        2 4 2 3 0 0
        0
        0
        0
        0
        """
    p = _write(tmp_path, fixture)
    with pytest.raises(ValueError) as exc_info:
        CHILmesh.read_from_fort14(p)

    err = str(exc_info.value)
    assert "element 2" in err


def test_node_id_above_count_raises(tmp_path):
    """Case c: node id (9) above node count (5)."""
    fixture = """
        id above count
        2 5
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        5 2.0 0.0 0.0
        1 3 1 2 3
        2 3 2 5 9
        0
        0
        0
        0
        """
    p = _write(tmp_path, fixture)
    with pytest.raises(ValueError) as exc_info:
        CHILmesh.read_from_fort14(p)

    err = str(exc_info.value)
    assert "9" in err


def test_valid_mixed_tri_quad_loads(tmp_path):
    """Case d: valid mixed tri/quad file loads successfully.

    Connectivity has no negative entry and triangle row is padded [v0, v1, v2, v0].
    """
    fixture = """
        valid mixed
        2 5
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        5 2.0 0.0 0.0
        1 4 1 2 3 4
        2 3 2 5 3
        0
        0
        0
        0
        """
    p = _write(tmp_path, fixture)
    mesh = CHILmesh.read_from_fort14(p, compute_layers=False, compute_adjacencies=False)

    # Mesh loaded successfully
    assert mesh.n_verts == 5
    assert mesh.n_elems == 2

    # Connectivity is 0-based (nodes 0..4)
    assert mesh.connectivity_list.shape == (2, 4)

    # Element 1 is quad: [0, 1, 2, 3] (original 1-based [1, 2, 3, 4])
    np.testing.assert_array_equal(mesh.connectivity_list[0], [0, 1, 2, 3])

    # Element 2 is triangle padded: [1, 4, 2, 1] (original 1-based [2, 5, 3, padded])
    np.testing.assert_array_equal(mesh.connectivity_list[1], [1, 4, 2, 1])

    # No negative indices
    assert np.all(mesh.connectivity_list >= 0)


def test_read_fort14_raw_preserves_zeros(tmp_path):
    """Case e: chilmesh.fort14_io.read_fort14_raw on the zero-id fixture
    does NOT raise and returns element 2 with raw ids (2, 5, 3, 0).
    """
    fixture = """
        pad0 repro
        2 5
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        5 2.0 0.0 0.0
        1 4 1 2 3 4
        2 4 2 5 3 0
        0
        0
        0
        0
        """
    p = _write(tmp_path, fixture)

    # read_fort14_raw should NOT raise
    raw = read_fort14_raw(p)

    # Element 2 should have raw ids preserved
    assert raw.elements[2] == (2, 5, 3, 0)
