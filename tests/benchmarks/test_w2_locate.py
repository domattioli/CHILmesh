"""Point-location harness: tree build timed apart from queries, outside count, real mesh."""
from __future__ import annotations

import numpy as np

from benchmarks.competitors.w2 import measures


class _Clock:
    def __init__(self):
        self.t = 0.0

    def __call__(self):
        return self.t


class _FakeMesh:
    """build costs 5 s; each query costs 1 s; points with x > 10 are outside."""

    def __init__(self, clock):
        self.clock = clock
        self.builds = 0
        self.queries = 0

    def _build_spatial_indices(self):
        self.builds += 1
        self.clock.t += 5.0

    def find_element(self, p):
        self.queries += 1
        self.clock.t += 1.0
        return -1 if p[0] > 10 else 0


def test_build_and_queries_timed_separately_without_double_counting():
    clock = _Clock()
    mesh = _FakeMesh(clock)
    pts = np.array([[1.0, 1.0], [2.0, 2.0], [20.0, 0.0], [3.0, 3.0]])
    res = measures.time_point_location(mesh, pts, clock=clock)
    assert res["tree_build_s"] == 5.0
    assert res["query_total_s"] == 4.0
    assert res["query_per_s"] == 1.0
    assert res["n_queries"] == 4
    assert res["n_outside"] == 1 and res["n_inside"] == 3
    assert mesh.builds == 1 and mesh.queries == 4
    assert clock.t == 9.0  # 5 + 4: nothing counted twice


def test_real_mesh_centroid_queries_are_all_inside(grid):
    from chilmesh import CHILmesh

    kind, pts, conn = grid
    mesh = CHILmesh(connectivity=conn, points=pts, grid_name=kind)
    q = measures.sample_centroids(mesh.points, mesh.connectivity_list, 20, seed=1)
    res = measures.time_point_location(mesh, q)
    assert res["n_queries"] == 20
    assert res["n_outside"] == 0
    assert res["tree_build_s"] >= 0 and res["query_total_s"] >= 0


def test_outside_points_are_reported_not_hidden():
    from chilmesh import CHILmesh
    from benchmarks.competitors.w2.measures import time_point_location

    q = np.array([[0.2, 0.2], [100.0, 100.0], [-5.0, 0.0]])
    # One element (k=1 KD-tree query) and two elements.
    single = CHILmesh(connectivity=np.array([[0, 1, 2]]),
                      points=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]), grid_name="t1")
    double = CHILmesh(connectivity=np.array([[0, 1, 2], [0, 2, 3]]),
                      points=np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]), grid_name="t2")
    for mesh in (single, double):
        res = time_point_location(mesh, q)
        assert res["n_outside"] == 2
        assert res["n_inside"] == 1


def test_sample_centroids_deterministic():
    pts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    conn = np.array([[0, 1, 2], [1, 3, 2]])
    a = measures.sample_centroids(pts, conn, 5, seed=3)
    b = measures.sample_centroids(pts, conn, 5, seed=3)
    assert np.array_equal(a, b) and a.shape == (5, 2)
