"""Shared quality and Courant metrics on hand-checked tri, quad and mixed meshes."""
from __future__ import annotations

import math

import numpy as np
import pytest

from _grids import mixed_grid
from benchmarks.competitors.metrics import courant as C
from benchmarks.competitors.metrics import quality as Q

SQ = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]


def test_equilateral_triangle_quality_is_one():
    pts = np.array([[0, 0], [1, 0], [0.5, math.sqrt(3) / 2]])
    q = Q.element_quality(pts, np.array([[0, 1, 2]]))
    assert q[0] == pytest.approx(1.0)


def test_right_isosceles_triangle_quality():
    pts = np.array([[0, 0], [1, 0], [0, 1.0]])
    q = Q.element_quality(pts, np.array([[0, 1, 2]]))
    assert q[0] == pytest.approx((math.sqrt(2) / 2) * 2 / math.sqrt(3))  # 0.8165


def test_unit_square_and_rectangle_quality_is_one():
    pts = np.array(SQ)
    assert Q.element_quality(pts, np.array([[0, 1, 2, 3]]))[0] == pytest.approx(1.0)
    rect = np.array([[0, 0], [2, 0], [2, 1], [0, 1.0]])
    assert Q.element_quality(rect, np.array([[0, 1, 2, 3]]))[0] == pytest.approx(1.0)


def test_clockwise_orientation_gives_same_quality():
    pts = np.array(SQ)
    assert Q.element_quality(pts, np.array([[0, 3, 2, 1]]))[0] == pytest.approx(1.0)


def test_padded_triangle_in_four_column_array_is_a_triangle():
    pts = np.array([[0, 0], [1, 0], [0, 1.0]])
    q3 = Q.element_quality(pts, np.array([[0, 1, 2]]))
    q4 = Q.element_quality(pts, np.array([[0, 1, 2, 0]]))
    assert q4[0] == pytest.approx(q3[0])


def test_sliver_has_two_small_angles_of_three():
    pts = np.array([[0, 0], [1, 0], [0.5, 0.05]])
    s = Q.quality_summary(pts, np.array([[0, 1, 2]]))
    assert s["pct_angles_lt_10"] == pytest.approx(200 / 3)


def test_summary_fields_and_bar(grid):
    kind, pts, conn = grid
    s = Q.quality_summary(pts, conn)
    for key in ("q_min", "q_mean", "q_mean_minus_3sigma", "meets_0_75_bar",
                "pct_angles_lt_10", "pct_irregular_valence", "n_elems"):
        assert key in s
    assert s["n_elems"] == len(conn)
    assert s["q_mean_minus_3sigma"] == pytest.approx(s["q_mean"] - 3 * np.std(Q.element_quality(pts, conn)))
    assert s["meets_0_75_bar"] == (s["q_mean_minus_3sigma"] >= 0.75)


def test_quad_grid_is_perfect():
    pts = np.array([[i, j] for j in range(5) for i in range(5)], dtype=float)
    conn = np.array([[j * 5 + i, j * 5 + i + 1, (j + 1) * 5 + i + 1, (j + 1) * 5 + i]
                     for j in range(4) for i in range(4)])
    s = Q.quality_summary(pts, conn)
    assert s["q_min"] == pytest.approx(1.0)
    assert s["pct_angles_lt_10"] == 0.0
    assert s["pct_irregular_valence"] == 0.0  # interior nodes all valence 4


def test_hexagon_fan_has_regular_center():
    ang = np.linspace(0, 2 * np.pi, 7)[:-1]
    pts = np.vstack([[0, 0], np.column_stack([np.cos(ang), np.sin(ang)])])
    conn = np.array([[0, 1 + i, 1 + (i + 1) % 6] for i in range(6)])
    s = Q.quality_summary(pts, conn)
    assert s["q_min"] == pytest.approx(1.0)
    assert s["pct_irregular_valence"] == 0.0


def test_pentagon_fan_center_is_irregular():
    ang = np.linspace(0, 2 * np.pi, 6)[:-1]
    pts = np.vstack([[0, 0], np.column_stack([np.cos(ang), np.sin(ang)])])
    conn = np.array([[0, 1 + i, 1 + (i + 1) % 5] for i in range(5)])
    s = Q.quality_summary(pts, conn)
    assert s["pct_irregular_valence"] == 100.0


def test_no_interior_nodes_gives_none():
    pts = np.array([[0, 0], [1, 0], [0, 1.0]])
    s = Q.quality_summary(pts, np.array([[0, 1, 2]]))
    assert s["pct_irregular_valence"] is None


def test_mixed_grid_summary_runs_and_counts_both_types():
    pts, conn = mixed_grid()
    s = Q.quality_summary(pts, conn)
    assert s["n_elems"] == len(conn)
    assert 0 < s["q_min"] <= 1.0


# --- Courant -----------------------------------------------------------------

def test_dt_max_two_nodes_hand_check():
    pts = np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]])
    conn = np.array([[0, 1, 2]])
    depths = np.full(3, 10.0)
    # shortest edge 100 m, h = 10 m: C(dt=1) = sqrt(9.81*10)/100
    c1 = math.sqrt(9.81 * 10) / 100
    assert C.max_courant(pts, conn, depths, dt=1.0) == pytest.approx(c1)
    assert C.dt_max(pts, conn, depths, courant_max=1.0) == pytest.approx(1.0 / c1)
    assert C.dt_max(pts, conn, depths, courant_max=0.5) == pytest.approx(0.5 / c1)


def test_courant_is_linear_in_dt():
    pts = np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]])
    conn = np.array([[0, 1, 2]])
    depths = np.full(3, 10.0)
    assert C.max_courant(pts, conn, depths, dt=4.0) == pytest.approx(4 * C.max_courant(pts, conn, depths, dt=1.0))


def test_quad_edges_are_the_four_sides_not_diagonals():
    pts = np.array([[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]])
    edges = C.edges_from_connectivity(np.array([[0, 1, 2, 3]]))
    assert len(edges) == 4


def test_padded_triangle_gives_three_edges():
    edges = C.edges_from_connectivity(np.array([[0, 1, 2, 0]]))
    assert len(edges) == 3


def test_shared_edges_counted_once(grid):
    kind, pts, conn = grid
    edges = C.edges_from_connectivity(conn)
    assert len(edges) == len({tuple(sorted(e)) for e in edges.tolist()})


def test_all_dry_gives_none():
    pts = np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]])
    conn = np.array([[0, 1, 2]])
    assert C.dt_max(pts, conn, np.zeros(3)) is None
    assert C.max_courant(pts, conn, np.zeros(3), dt=1.0) is None


def test_mixed_mesh_dt_max_is_finite_positive():
    pts, conn = mixed_grid()
    assert C.dt_max(pts, conn, np.full(len(pts), 5.0)) > 0
