"""Synthetic negative fixtures for spec 007.

Each helper builds a mesh that plants exactly one violation. Because
``CHILmesh.__init__`` triggers a degeneracy fallback that would silently
repair some planted defects, fixtures construct a valid mesh first and
then post-mutate ``connectivity_list``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

import chilmesh


def _invalidate_caches(mesh: Any) -> None:
    if hasattr(mesh, "adjacencies"):
        mesh.adjacencies = {}
    if hasattr(mesh, "layers"):
        mesh.layers = {"OE": [], "IE": [], "OV": [], "IV": [], "bEdgeIDs": []}
    if hasattr(mesh, "n_layers"):
        mesh.n_layers = 0


def _minimal_quad_grid() -> Any:
    """3x3 vertex grid → 2x2 quad mesh (4 elements)."""
    pts = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
            [0.0, 2.0, 0.0],
            [1.0, 2.0, 0.0],
            [2.0, 2.0, 0.0],
        ]
    )
    conn = np.array(
        [
            [0, 1, 4, 3],
            [1, 2, 5, 4],
            [3, 4, 7, 6],
            [4, 5, 8, 7],
        ]
    )
    return chilmesh.CHILmesh(connectivity=conn, points=pts, compute_layers=True)


def bowtie_quad_mesh() -> Any:
    """Plant a bowtie quad by swapping v2 and v3 of element 0."""
    mesh = _minimal_quad_grid()
    new_row = np.array([0, 1, 3, 4])
    mesh.connectivity_list[0] = new_row
    _invalidate_caches(mesh)
    return mesh


def interior_triangle_mesh() -> Any:
    """Plant an interior triangle: take the center element of a 3x3 quad grid and turn it into a tri whose verts are all interior."""
    pts = np.array(
        [
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0],
            [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0], [3.0, 1.0, 0.0],
            [0.0, 2.0, 0.0], [1.0, 2.0, 0.0], [2.0, 2.0, 0.0], [3.0, 2.0, 0.0],
            [0.0, 3.0, 0.0], [1.0, 3.0, 0.0], [2.0, 3.0, 0.0], [3.0, 3.0, 0.0],
        ]
    )
    conn = np.array(
        [
            [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6],
            [4, 5, 9, 8], [5, 6, 10, 9], [6, 7, 11, 10],
            [8, 9, 13, 12], [9, 10, 14, 13], [10, 11, 15, 14],
        ]
    )
    mesh = chilmesh.CHILmesh(connectivity=conn, points=pts, compute_layers=True)
    mesh.connectivity_list[4] = np.array([5, 6, 10, 5])
    _invalidate_caches(mesh)
    return mesh


@dataclass
class _FakeMesh:
    """Test double for pentagon_mesh — CHILmesh connectivity is 4-column, so a real
    5-vertex element cannot be constructed via the real constructor.
    """

    connectivity_list: np.ndarray
    points: np.ndarray
    n_elems: int
    layers: dict
    n_verts: int = 0
    n_layers: int = 1

    def boundary_node_indices(self) -> np.ndarray:
        return np.unique(self.connectivity_list[self.connectivity_list != -1])


def pentagon_mesh() -> Any:
    pts = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.5, 0.5, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    conn = np.array([[0, 1, 2, 3, 4]])
    return _FakeMesh(
        connectivity_list=conn,
        points=pts,
        n_elems=1,
        layers={"OE": [np.array([0])], "IE": [np.array([])], "OV": [], "IV": [], "bEdgeIDs": []},
        n_verts=5,
    )


def overlapping_quads_mesh() -> Any:
    """Two disjoint-vertex quads whose interiors overlap.

    Two 4-vertex squares share no vertex IDs but the second sits inside the first.
    """
    pts = np.array(
        [
            [0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [4.0, 4.0, 0.0], [0.0, 4.0, 0.0],
            [1.0, 1.0, 0.0], [3.0, 1.0, 0.0], [3.0, 3.0, 0.0], [1.0, 3.0, 0.0],
        ]
    )
    conn = np.array(
        [
            [0, 1, 2, 3],
            [4, 5, 6, 7],
        ]
    )
    mesh = chilmesh.CHILmesh(connectivity=conn, points=pts, compute_layers=False)
    _invalidate_caches(mesh)
    return mesh


def edge_crossing_mesh() -> Any:
    """Two disjoint quads with edges that cross geometrically (no shared vertex)."""
    pts = np.array(
        [
            [0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [4.0, 2.0, 0.0], [0.0, 2.0, 0.0],
            [1.0, -1.0, 0.0], [3.0, -1.0, 0.0], [3.0, 4.0, 0.0], [1.0, 4.0, 0.0],
        ]
    )
    conn = np.array(
        [
            [0, 1, 2, 3],
            [4, 5, 6, 7],
        ]
    )
    mesh = chilmesh.CHILmesh(connectivity=conn, points=pts, compute_layers=False)
    _invalidate_caches(mesh)
    return mesh
