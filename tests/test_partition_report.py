"""Tests for partition_report (spec 005 US3)."""
import numpy as np
import pytest

from chilmesh import CHILmesh, partition_report


def _strip(n):
    """Strip of n quads in a row: elements e and e+1 share an edge."""
    pts = [[i, j, 0] for i in range(n + 1) for j in (0, 1)]
    conn = [[2 * i, 2 * i + 2, 2 * i + 3, 2 * i + 1] for i in range(n)]
    return CHILmesh(connectivity=np.array(conn), points=np.array(pts, float))


def test_sizes_imbalance_cut_pieces():
    m = _strip(4)
    r = partition_report(m, np.array([0, 0, 0, 1]))
    assert r.sizes.tolist() == [3, 1]
    assert r.imbalance == pytest.approx(3 / 2)
    assert r.cut_edges == 1
    assert r.pieces.tolist() == [1, 1]


def test_disconnected_part():
    m = _strip(4)
    r = partition_report(m, np.array([0, 1, 0, 1]))
    assert r.cut_edges == 3
    assert r.pieces.tolist() == [2, 2]


def test_k_bounds():
    m = _strip(4)
    with pytest.raises(ValueError):
        partition_report(m, np.zeros(4, dtype=int))  # k = 1
    with pytest.raises(ValueError):
        partition_report(m, np.array([0, 1]))  # wrong length
    with pytest.raises(ValueError):
        partition_report(m, np.array([0, 1, 2, 9]))  # k=10 > n_elems
