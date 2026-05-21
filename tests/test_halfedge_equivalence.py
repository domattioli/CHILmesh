"""Verify bit-identical adjacency between EdgeMap and half-edge backends.

Implements canonical-form comparator per clarify C-4:
- Edge2Vert: sorted lexicographically
- Edge2Elem, Elem2Edge: row-by-row set comparison (order-agnostic)
- Vert2Edge, Vert2Elem: as dicts, set values compared
"""

import pytest
import numpy as np
from chilmesh.examples import annulus, donut, block_o, structured


@pytest.mark.parametrize(
    "fixture_fn",
    [annulus, donut, block_o, structured],
    ids=["annulus", "donut", "block_o", "structured"]
)
def test_edge2vert_equivalence(fixture_fn):
    """Edge2Vert canonical form matches between backends."""
    mesh_edgemap = fixture_fn()
    mesh_halfedge = fixture_fn(topology_backend="halfedge")

    e2v_edgemap = mesh_edgemap.adjacencies["Edge2Vert"]
    e2v_halfedge = mesh_halfedge.adjacencies["Edge2Vert"]

    assert len(e2v_edgemap) == len(e2v_halfedge), "Edge count mismatch"

    e2v_edgemap_sorted = sorted(map(tuple, e2v_edgemap))
    e2v_halfedge_sorted = sorted(map(tuple, e2v_halfedge))

    assert e2v_edgemap_sorted == e2v_halfedge_sorted, "Edge2Vert mismatch"


@pytest.mark.parametrize(
    "fixture_fn",
    [annulus, donut, block_o, structured],
    ids=["annulus", "donut", "block_o", "structured"]
)
def test_edge2elem_equivalence(fixture_fn):
    """Edge2Elem matches between backends (same adjacency list)."""
    mesh_edgemap = fixture_fn()
    mesh_halfedge = fixture_fn(topology_backend="halfedge")

    e2m_edgemap = mesh_edgemap.adjacencies["Edge2Elem"]
    e2m_halfedge = mesh_halfedge.adjacencies["Edge2Elem"]

    assert e2m_edgemap.shape == e2m_halfedge.shape, "Edge2Elem shape mismatch"

    edgemap_list = sorted(map(tuple, mesh_edgemap.adjacencies["Edge2Vert"]))
    halfedge_list = sorted(map(tuple, mesh_halfedge.adjacencies["Edge2Vert"]))

    for i, (edge_em, edge_he) in enumerate(zip(edgemap_list, halfedge_list)):
        assert edge_em == edge_he, f"Edge {i} mismatch"
        em_adj = set(e2m_edgemap[i])
        he_adj = set(e2m_halfedge[i])
        assert em_adj == he_adj, f"Edge2Elem adjacency mismatch for edge {edge_em}"


@pytest.mark.parametrize(
    "fixture_fn",
    [annulus, donut, block_o, structured],
    ids=["annulus", "donut", "block_o", "structured"]
)
def test_elem2edge_equivalence(fixture_fn):
    """Elem2Edge matches between backends (set comparison per element)."""
    mesh_edgemap = fixture_fn()
    mesh_halfedge = fixture_fn(topology_backend="halfedge")

    e2e_edgemap = mesh_edgemap.adjacencies["Elem2Edge"]
    e2e_halfedge = mesh_halfedge.adjacencies["Elem2Edge"]

    assert e2e_edgemap.shape[0] == e2e_halfedge.shape[0], "Element count mismatch"

    for elem_idx in range(len(e2e_edgemap)):
        em_edges = set(e2e_edgemap[elem_idx])
        em_edges.discard(-1)
        he_edges = set(e2e_halfedge[elem_idx])
        he_edges.discard(-1)
        assert em_edges == he_edges, f"Elem2Edge mismatch for element {elem_idx}"


@pytest.mark.parametrize(
    "fixture_fn",
    [annulus, donut, block_o, structured],
    ids=["annulus", "donut", "block_o", "structured"]
)
def test_vert2edge_equivalence(fixture_fn):
    """Vert2Edge matches between backends (set of edge IDs per vertex)."""
    mesh_edgemap = fixture_fn()
    mesh_halfedge = fixture_fn(topology_backend="halfedge")

    v2e_edgemap = mesh_edgemap.adjacencies["Vert2Edge"]
    v2e_halfedge = mesh_halfedge.adjacencies["Vert2Edge"]

    assert len(v2e_edgemap) == len(v2e_halfedge), "Vertex count mismatch"

    for vert_idx in range(len(v2e_edgemap)):
        em_edges = v2e_edgemap[vert_idx]
        he_edges = v2e_halfedge[vert_idx]
        assert em_edges == he_edges, f"Vert2Edge mismatch for vertex {vert_idx}"


@pytest.mark.parametrize(
    "fixture_fn",
    [annulus, donut, block_o, structured],
    ids=["annulus", "donut", "block_o", "structured"]
)
def test_vert2elem_equivalence(fixture_fn):
    """Vert2Elem matches between backends (set of element IDs per vertex)."""
    mesh_edgemap = fixture_fn()
    mesh_halfedge = fixture_fn(topology_backend="halfedge")

    v2m_edgemap = mesh_edgemap.adjacencies["Vert2Elem"]
    v2m_halfedge = mesh_halfedge.adjacencies["Vert2Elem"]

    assert len(v2m_edgemap) == len(v2m_halfedge), "Vertex count mismatch"

    for vert_idx in range(len(v2m_edgemap)):
        em_elems = v2m_edgemap[vert_idx]
        he_elems = v2m_halfedge[vert_idx]
        assert em_elems == he_elems, f"Vert2Elem mismatch for vertex {vert_idx}"


@pytest.mark.parametrize(
    "fixture_fn",
    [annulus, donut, block_o, structured],
    ids=["annulus", "donut", "block_o", "structured"]
)
def test_elem2vert_unchanged(fixture_fn):
    """Elem2Vert is identical in both backends (not reordered)."""
    mesh_edgemap = fixture_fn()
    mesh_halfedge = fixture_fn(topology_backend="halfedge")

    e2v_edgemap = mesh_edgemap.adjacencies["Elem2Vert"]
    e2v_halfedge = mesh_halfedge.adjacencies["Elem2Vert"]

    np.testing.assert_array_equal(e2v_edgemap, e2v_halfedge)
