"""Tests for the compute_adjacencies flag (#134).

When set to True with compute_layers=False, the constructor must build adjacency
dicts (Vert2Edge, Vert2Elem, Edge2Vert, Edge2Elem, Elem2Edge, EdgeMap) without
running the skeletonization pass. This lets downstream consumers
(quadmesh-matlab, MADMESHR) use the public adjacency API without paying the
layer-sweep cost or calling private methods.
"""
from __future__ import annotations

import numpy as np
import pytest

import chilmesh


_ADJ_KEYS = {"Elem2Vert", "Edge2Vert", "Elem2Edge", "Vert2Edge", "Vert2Elem", "Edge2Elem", "EdgeMap"}


def test_layers_false_skips_adjacencies_by_default():
    """Backward-compat: compute_layers=False with no compute_adjacencies → no adjacencies."""
    fixture_path = chilmesh.examples.fixture_path("annulus_200pts.fort.14")
    mesh = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=False)
    assert mesh.adjacencies == {}
    assert mesh.n_layers == 0


def test_layers_false_adjacencies_true_builds_adjacencies():
    """New behavior: compute_layers=False + compute_adjacencies=True → adjacencies built, no layers."""
    fixture_path = chilmesh.examples.fixture_path("annulus_200pts.fort.14")
    mesh = chilmesh.CHILmesh.read_from_fort14(
        fixture_path, compute_layers=False, compute_adjacencies=True
    )
    assert mesh.adjacencies != {}
    assert _ADJ_KEYS.issubset(mesh.adjacencies.keys())
    assert mesh.n_edges > 0
    assert mesh.n_layers == 0
    assert mesh.layers["OE"] == []


def test_layers_true_forces_adjacencies_even_if_user_says_no():
    """compute_layers=True + compute_adjacencies=False → adjacencies still built (layers depend on them)."""
    fixture_path = chilmesh.examples.fixture_path("annulus_200pts.fort.14")
    mesh = chilmesh.CHILmesh.read_from_fort14(
        fixture_path, compute_layers=True, compute_adjacencies=False
    )
    assert mesh.adjacencies != {}
    assert mesh.n_layers > 0


def test_public_adjacency_api_works_with_compute_adjacencies_only():
    """Verifies the downstream use case from #134: public API works without compute_layers."""
    fixture_path = chilmesh.examples.fixture_path("annulus_200pts.fort.14")
    mesh = chilmesh.CHILmesh.read_from_fort14(
        fixture_path, compute_layers=False, compute_adjacencies=True
    )

    # Public adjacency queries should not require calling private _build_adjacencies()
    edges_around_v0 = mesh.get_vertex_edges(0)
    assert isinstance(edges_around_v0, set)
    assert len(edges_around_v0) > 0

    elems_around_v0 = mesh.get_vertex_elements(0)
    assert isinstance(elems_around_v0, set)
    assert len(elems_around_v0) > 0

    # boundary_edges should hit the fast (adjacency) path
    boundary = mesh.boundary_edges()
    assert boundary.dtype.kind in "iu"
    assert len(boundary) > 0


def test_direct_constructor_supports_compute_adjacencies():
    """Constructor accepts compute_adjacencies directly (not just file readers)."""
    connectivity = np.array([[0, 1, 2], [0, 2, 3]], dtype=int)
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    mesh = chilmesh.CHILmesh(
        connectivity=connectivity,
        points=points,
        compute_layers=False,
        compute_adjacencies=True,
    )
    assert _ADJ_KEYS.issubset(mesh.adjacencies.keys())
    assert mesh.n_layers == 0


def test_from_admesh_domain_threads_compute_adjacencies():
    """from_admesh_domain forwards compute_adjacencies to the reader."""
    from types import SimpleNamespace

    fixture_path = chilmesh.examples.fixture_path("quad_2x2.fort.14")
    record = SimpleNamespace(filename=str(fixture_path), type="ADCIRC")

    mesh = chilmesh.CHILmesh.from_admesh_domain(
        record, compute_layers=False, compute_adjacencies=True
    )
    assert mesh.adjacencies != {}
    assert mesh.n_layers == 0


def test_geometry_matches_full_init():
    """compute_adjacencies-only mesh has same topology as full init."""
    fixture_path = chilmesh.examples.fixture_path("annulus_200pts.fort.14")
    mesh_partial = chilmesh.CHILmesh.read_from_fort14(
        fixture_path, compute_layers=False, compute_adjacencies=True
    )
    mesh_full = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=True)

    assert mesh_partial.n_verts == mesh_full.n_verts
    assert mesh_partial.n_elems == mesh_full.n_elems
    assert mesh_partial.n_edges == mesh_full.n_edges
    assert mesh_partial.type == mesh_full.type
