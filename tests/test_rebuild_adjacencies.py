"""Tests for public ``CHILmesh.rebuild_adjacencies`` and ``invalidate_adjacencies`` (#143).

Downstream consumers (QuADMesh aggressive routing, MADMESHR adapt loops) mutate
``connectivity_list`` between adjacency queries. These two methods give them a
public, non-private way to refresh / drop the cached adjacency dicts without
reconstructing a fresh ``CHILmesh`` (which loses ``grid_name`` and is the slow
path).
"""
from __future__ import annotations

import numpy as np
import pytest

import chilmesh


_ADJ_KEYS = {"Elem2Vert", "Edge2Vert", "Elem2Edge", "Vert2Edge", "Vert2Elem", "Edge2Elem", "EdgeMap"}


def _fresh_mesh() -> chilmesh.CHILmesh:
    fixture_path = chilmesh.examples.fixture_path("annulus_200pts.fort.14")
    return chilmesh.CHILmesh.read_from_fort14(fixture_path)


def test_invalidate_clears_adjacencies():
    mesh = _fresh_mesh()
    assert _ADJ_KEYS.issubset(mesh.adjacencies.keys())
    mesh.invalidate_adjacencies()
    assert mesh.adjacencies == {}
    assert mesh.n_edges == 0


def test_invalidate_preserves_connectivity_and_points():
    mesh = _fresh_mesh()
    conn_before = mesh.connectivity_list.copy()
    pts_before = mesh.points.copy()
    grid_name_before = mesh.grid_name

    mesh.invalidate_adjacencies()

    assert np.array_equal(mesh.connectivity_list, conn_before)
    assert np.array_equal(mesh.points, pts_before)
    assert mesh.grid_name == grid_name_before


def test_rebuild_restores_full_adjacency_bundle():
    mesh = _fresh_mesh()
    edges_before = mesh.n_edges
    edge2vert_before = mesh.adjacencies["Edge2Vert"].copy()

    mesh.invalidate_adjacencies()
    mesh.rebuild_adjacencies()

    assert _ADJ_KEYS.issubset(mesh.adjacencies.keys())
    assert mesh.n_edges == edges_before
    # Edge enumeration is canonical (sorted (min, max)) — bit-identical after rebuild.
    np.testing.assert_array_equal(mesh.adjacencies["Edge2Vert"], edge2vert_before)


def test_rebuild_after_connectivity_swap_picks_up_new_edges():
    """Simulate mid-sweep mutation: swap an edge's vertices. Rebuild must reflect it."""
    connectivity = np.array([[0, 1, 2], [0, 2, 3]], dtype=int)
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    mesh = chilmesh.CHILmesh(
        connectivity=connectivity,
        points=points,
        compute_layers=False,
        compute_adjacencies=True,
    )
    edges_before = int(mesh.n_edges)

    # Mutate: add a new triangle that introduces fresh edges.
    new_points = np.vstack([mesh.points, np.array([[2.0, 0.5, 0.0]])])
    new_conn = np.vstack([mesh.connectivity_list, np.array([[1, 4, 2]], dtype=int)])
    mesh.points = new_points
    mesh.connectivity_list = new_conn
    mesh.n_verts = new_points.shape[0]
    mesh.n_elems = new_conn.shape[0]

    mesh.rebuild_adjacencies()

    assert mesh.n_edges > edges_before
    # Each row referenced in connectivity must now appear in Vert2Elem.
    assert 4 in mesh.adjacencies["Vert2Elem"]
    assert 2 in mesh.adjacencies["Vert2Elem"][4]


def test_rebuild_preserves_grid_name():
    mesh = _fresh_mesh()
    name = mesh.grid_name
    mesh.rebuild_adjacencies()
    assert mesh.grid_name == name


def test_rebuild_works_from_compute_adjacencies_false_initial_state():
    fixture_path = chilmesh.examples.fixture_path("annulus_200pts.fort.14")
    mesh = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=False)
    assert mesh.adjacencies == {}

    mesh.rebuild_adjacencies()

    assert _ADJ_KEYS.issubset(mesh.adjacencies.keys())
    assert mesh.n_edges > 0


def test_rebuild_no_spatial_skip_keeps_existing_kdtree():
    mesh = _fresh_mesh()
    original_tree = mesh._vertex_tree

    mesh.rebuild_adjacencies(rebuild_spatial_indices=False)

    # KD-tree object identity preserved when skip flag is set.
    assert mesh._vertex_tree is original_tree


def test_rebuild_default_refreshes_spatial_indices():
    mesh = _fresh_mesh()
    original_tree = mesh._vertex_tree

    mesh.rebuild_adjacencies()

    assert mesh._vertex_tree is not original_tree


def test_idempotent_rebuild_is_stable():
    mesh = _fresh_mesh()
    mesh.rebuild_adjacencies()
    snapshot = mesh.adjacencies["Edge2Vert"].copy()
    mesh.rebuild_adjacencies()
    np.testing.assert_array_equal(mesh.adjacencies["Edge2Vert"], snapshot)


def test_post_invalidate_public_query_raises():
    mesh = _fresh_mesh()
    mesh.invalidate_adjacencies()
    with pytest.raises(RuntimeError):
        mesh.ccw_edges_around_vert(0)
