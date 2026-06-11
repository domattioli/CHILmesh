"""Tests for CHILmesh constructor skip flags (#204)."""
import numpy as np
import pytest
from chilmesh.CHILmesh import CHILmesh


def _simple():
    """Create a simple two-triangle mesh forming a unit square."""
    pts = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)
    conn = np.array([[0, 1, 2], [0, 2, 3]], dtype=int)
    return conn, pts


def test_defaults_build_spatial_index():
    """Verify default behavior builds spatial indices."""
    conn, pts = _simple()
    m = CHILmesh(conn, pts)
    assert hasattr(m, '_vertex_tree'), "Default should build _vertex_tree"
    assert hasattr(m, '_centroid_tree'), "Default should build _centroid_tree"


def test_skip_spatial_index():
    """Verify build_spatial_indices=False skips KD-tree construction."""
    conn, pts = _simple()
    m = CHILmesh(
        conn, pts,
        compute_layers=False,
        compute_adjacencies=True,
        build_spatial_indices=False
    )
    assert not hasattr(m, '_vertex_tree'), "Should skip _vertex_tree"
    assert not hasattr(m, '_centroid_tree'), "Should skip _centroid_tree"


def test_skip_validate_keeps_topology():
    """Verify validate=False skips validation but preserves adjacency data."""
    conn, pts = _simple()
    m = CHILmesh(
        conn, pts,
        compute_layers=False,
        compute_adjacencies=True,
        build_spatial_indices=False,
        validate=False
    )
    # Check that adjacencies are built despite skipping validation
    assert 'Edge2Vert' in m.adjacencies, "Edge2Vert should be present"
    assert 'Edge2Elem' in m.adjacencies, "Edge2Elem should be present"
    assert m.adjacencies['Edge2Vert'] is not None
    assert m.adjacencies['Edge2Elem'] is not None
    # Verify shape is reasonable
    assert m.adjacencies['Edge2Vert'].shape[0] > 0, "Should have edges"
    assert m.adjacencies['Edge2Elem'].shape[0] > 0, "Should have edges"


def test_defaults_unchanged_validate_runs():
    """Verify default behavior (validation enabled) works without error."""
    conn, pts = _simple()
    m = CHILmesh(conn, pts)
    # Should construct without raising from validation
    assert m.n_elems == 2
    assert m.n_verts == 4
    assert m.n_edges > 0
