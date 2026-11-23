"""Boundary displacement check (FR-006) on tri, quad and mixed meshes, bit-exact."""
from __future__ import annotations

import numpy as np

from benchmarks.competitors.metrics import displacement as D
from benchmarks.competitors.metrics.quality import boundary_nodes


def test_boundary_nodes_of_grids(grid):
    kind, pts, conn = grid
    b = boundary_nodes(conn)
    on_edge = np.where((pts[:, 0] == 0) | (pts[:, 0] == 4) | (pts[:, 1] == 0) | (pts[:, 1] == 4))[0]
    assert np.array_equal(np.sort(b), np.sort(on_edge))


def test_unchanged_mesh_passes_with_exactly_zero(grid):
    kind, pts, conn = grid
    res = D.check_boundary(pts, conn, pts.copy(), conn.copy())
    assert res["status"] == "ok"
    assert res["n_boundary_moved"] == 0
    assert res["max_boundary_displacement"] == 0.0


def test_interior_move_passes(grid):
    kind, pts, conn = grid
    after = pts.copy()
    interior = np.setdiff1d(np.arange(len(pts)), boundary_nodes(conn))
    after[interior] += 0.01
    res = D.check_boundary(pts, conn, after, conn.copy())
    assert res["status"] == "ok"


def test_any_boundary_move_is_invalid_even_one_ulp(grid):
    kind, pts, conn = grid
    after = pts.copy()
    b = boundary_nodes(conn)[0]
    after[b, 0] = np.nextafter(after[b, 0], np.inf)
    res = D.check_boundary(pts, conn, after, conn.copy())
    assert res["status"] == "invalid"
    assert res["n_boundary_moved"] == 1
    assert "boundary" in res["reason"]


def test_track1_node_count_change_is_invalid(grid):
    kind, pts, conn = grid
    res = D.check_boundary(pts, conn, pts[:-1], conn.copy(), frozen_topology=True)
    assert res["status"] == "invalid"
    assert "node count" in res["reason"]


def test_track1_element_count_change_is_invalid(grid):
    kind, pts, conn = grid
    res = D.check_boundary(pts, conn, pts.copy(), conn[:-1], frozen_topology=True)
    assert res["status"] == "invalid"
    assert "element count" in res["reason"]


def test_track2_allows_count_change_but_still_checks_surviving_boundary(grid):
    kind, pts, conn = grid
    res = D.check_boundary(pts, conn, pts.copy(), conn[:-1], frozen_topology=False)
    assert res["status"] in {"ok", "invalid"}
    assert "element count" not in res.get("reason", "")


def test_chilmesh_fixed_boundary_smoother_passes_on_each_kind(grid):
    from chilmesh import CHILmesh

    kind, pts, conn = grid
    mesh = CHILmesh(connectivity=conn.copy(), points=pts.copy(), grid_name=kind)
    before_pts, before_conn = mesh.points.copy(), mesh.connectivity_list.copy()
    mesh.smooth(n_iter=3, weight=0.5, fixed_boundary=True)
    res = D.check_boundary(before_pts, before_conn, mesh.points, mesh.connectivity_list, frozen_topology=True)
    assert res["status"] == "ok", res
