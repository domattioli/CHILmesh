"""Tests for public ``CHILmesh.ccw_edges_around_vert`` helper (#133).

Returns the edges incident to a vertex in counterclockwise order (sorted by
``atan2(dy, dx)`` from the vertex to the other endpoint of each edge).
"""
from __future__ import annotations

import math

import numpy as np
import pytest

import chilmesh


def _other_endpoint(mesh: chilmesh.CHILmesh, vert_id: int, edge_id: int) -> int:
    a, b = mesh.adjacencies["Edge2Vert"][edge_id]
    return int(b) if a == vert_id else int(a)


def _edge_angle(mesh: chilmesh.CHILmesh, vert_id: int, edge_id: int) -> float:
    other = _other_endpoint(mesh, vert_id, edge_id)
    dx = mesh.points[other, 0] - mesh.points[vert_id, 0]
    dy = mesh.points[other, 1] - mesh.points[vert_id, 1]
    return math.atan2(dy, dx)


class TestCCWEdgesAroundVert:
    def test_returns_all_incident_edges(self, mesh):
        """Result must contain exactly the edges in Vert2Edge[vert_id]."""
        vert_id = int(mesh.boundary_node_indices()[0])
        ccw = mesh.ccw_edges_around_vert(vert_id)
        assert set(ccw) == set(mesh.adjacencies["Vert2Edge"][vert_id])

    def test_returns_list_of_ints(self, mesh):
        ccw = mesh.ccw_edges_around_vert(0)
        assert isinstance(ccw, list)
        for edge_id in ccw:
            assert isinstance(edge_id, int)

    def test_order_is_monotone_in_angle(self, mesh):
        """Angles from the vertex to each neighbor must be monotone non-decreasing."""
        for vert_id in [0, mesh.n_verts // 2, mesh.n_verts - 1]:
            ccw = mesh.ccw_edges_around_vert(vert_id)
            if len(ccw) < 2:
                continue
            angles = [_edge_angle(mesh, vert_id, e) for e in ccw]
            for a, b in zip(angles, angles[1:]):
                assert a <= b + 1e-12, (
                    f"Angles not monotone for vert {vert_id} in {mesh.grid_name}: {angles}"
                )

    def test_synthetic_square_interior_vertex(self):
        """Hand-checkable case: vertex at origin with 4 neighbors at +x, +y, -x, -y."""
        points = np.array(
            [
                [0.0, 0.0, 0.0],   # 0 (center)
                [1.0, 0.0, 0.0],   # 1 (east)
                [0.0, 1.0, 0.0],   # 2 (north)
                [-1.0, 0.0, 0.0],  # 3 (west)
                [0.0, -1.0, 0.0],  # 4 (south)
            ]
        )
        connectivity = np.array(
            [
                [0, 1, 2],
                [0, 2, 3],
                [0, 3, 4],
                [0, 4, 1],
            ],
            dtype=int,
        )
        mesh = chilmesh.CHILmesh(
            connectivity=connectivity, points=points, compute_layers=True
        )
        ccw = mesh.ccw_edges_around_vert(0)
        # Expect 4 incident edges, ordered E → N → W → S (atan2 ascending from -pi to pi).
        assert len(ccw) == 4
        ordered_neighbors = [_other_endpoint(mesh, 0, e) for e in ccw]
        # atan2 returns -pi for west; ordering by ascending atan2 is: south(-pi/2), east(0), north(pi/2), west(pi).
        assert ordered_neighbors == [4, 1, 2, 3]

    def test_empty_for_isolated_vertex_index(self):
        """Out-of-range vert_id raises rather than returning empty."""
        fixture_path = chilmesh.examples.fixture_path("annulus_200pts.fort.14")
        mesh = chilmesh.CHILmesh.read_from_fort14(fixture_path)
        with pytest.raises(ValueError, match="out of range"):
            mesh.ccw_edges_around_vert(mesh.n_verts)
        with pytest.raises(ValueError, match="out of range"):
            mesh.ccw_edges_around_vert(-1)

    def test_raises_without_adjacencies(self):
        """Without adjacency build, helper must raise a clear error."""
        fixture_path = chilmesh.examples.fixture_path("annulus_200pts.fort.14")
        mesh = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=False)
        with pytest.raises(RuntimeError, match="Adjacencies not built"):
            mesh.ccw_edges_around_vert(0)

    def test_works_with_compute_adjacencies_only(self):
        """Works without skeletonization as long as adjacencies were built (#133 + #134)."""
        fixture_path = chilmesh.examples.fixture_path("annulus_200pts.fort.14")
        mesh = chilmesh.CHILmesh.read_from_fort14(
            fixture_path, compute_layers=False, compute_adjacencies=True
        )
        ccw = mesh.ccw_edges_around_vert(0)
        assert isinstance(ccw, list)
        assert len(ccw) > 0

    def test_deterministic_across_calls(self, mesh):
        """Repeated calls return identical order."""
        vert_id = int(mesh.boundary_node_indices()[0])
        ccw_a = mesh.ccw_edges_around_vert(vert_id)
        ccw_b = mesh.ccw_edges_around_vert(vert_id)
        assert ccw_a == ccw_b
