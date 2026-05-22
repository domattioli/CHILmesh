"""Verify bit-identical adjacency between EdgeMap and half-edge backends.

Implements canonical-form comparator per clarify C-4:
- Edge2Vert: sorted lexicographically
- Edge2Elem, Elem2Edge: row-by-row set comparison (order-agnostic)
- Vert2Edge, Vert2Elem: as dicts, set values compared (by edge tuples, not IDs)
"""

import pytest
import numpy as np
from chilmesh.examples import annulus, donut, block_o, structured


def _edge_sort_order(e2v):
    """Return argsort indices mapping sorted-edge-index → original edge index."""
    return sorted(range(len(e2v)), key=lambda i: tuple(sorted(e2v[i])))


def _normalized_edge_set(e2v, edge_ids):
    """Convert a set of edge IDs to a set of normalized (min,max) edge tuples."""
    return {tuple(sorted(e2v[eid])) for eid in edge_ids}


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
    """Edge2Elem matches between backends (same adjacency list).

    Edges are identified by their normalized vertex-pair tuple, not their
    local edge-ID (which differs between backends).
    """
    mesh_edgemap = fixture_fn()
    mesh_halfedge = fixture_fn(topology_backend="halfedge")

    e2v_em = mesh_edgemap.adjacencies["Edge2Vert"]
    e2v_he = mesh_halfedge.adjacencies["Edge2Vert"]
    e2m_em = mesh_edgemap.adjacencies["Edge2Elem"]
    e2m_he = mesh_halfedge.adjacencies["Edge2Elem"]

    assert e2m_em.shape == e2m_he.shape, "Edge2Elem shape mismatch"

    # Build mapping from normalized edge tuple → element set for each backend
    em_map = {tuple(sorted(e2v_em[i])): set(e2m_em[i]) for i in range(len(e2v_em))}
    he_map = {tuple(sorted(e2v_he[i])): set(e2m_he[i]) for i in range(len(e2v_he))}

    assert set(em_map.keys()) == set(he_map.keys()), "Edge sets differ"

    for edge_key in em_map:
        assert em_map[edge_key] == he_map[edge_key], (
            f"Edge2Elem adjacency mismatch for edge {edge_key}: "
            f"edgemap={em_map[edge_key]} halfedge={he_map[edge_key]}"
        )


@pytest.mark.parametrize(
    "fixture_fn",
    [annulus, donut, block_o, structured],
    ids=["annulus", "donut", "block_o", "structured"]
)
def test_elem2edge_equivalence(fixture_fn):
    """Elem2Edge matches between backends (set comparison per element).

    Edges are identified by their normalized vertex-pair tuple so that
    differing edge-ID numbering schemes do not cause false failures.
    """
    mesh_edgemap = fixture_fn()
    mesh_halfedge = fixture_fn(topology_backend="halfedge")

    e2v_em = mesh_edgemap.adjacencies["Edge2Vert"]
    e2v_he = mesh_halfedge.adjacencies["Edge2Vert"]
    e2e_em = mesh_edgemap.adjacencies["Elem2Edge"]
    e2e_he = mesh_halfedge.adjacencies["Elem2Edge"]

    assert e2e_em.shape[0] == e2e_he.shape[0], "Element count mismatch"

    for elem_idx in range(len(e2e_em)):
        # Convert edge IDs to normalized edge tuples for comparison
        em_edge_ids = [eid for eid in e2e_em[elem_idx] if eid != -1]
        he_edge_ids = [eid for eid in e2e_he[elem_idx] if eid != -1]

        em_edges = _normalized_edge_set(e2v_em, em_edge_ids)
        he_edges = _normalized_edge_set(e2v_he, he_edge_ids)

        assert em_edges == he_edges, (
            f"Elem2Edge mismatch for element {elem_idx}: "
            f"edgemap={em_edges} halfedge={he_edges}"
        )


@pytest.mark.parametrize(
    "fixture_fn",
    [annulus, donut, block_o, structured],
    ids=["annulus", "donut", "block_o", "structured"]
)
def test_vert2edge_equivalence(fixture_fn):
    """Vert2Edge matches between backends.

    Comparison is by normalized edge tuples (min_v, max_v) rather than
    raw edge IDs, which differ between backends.
    """
    mesh_edgemap = fixture_fn()
    mesh_halfedge = fixture_fn(topology_backend="halfedge")

    e2v_em = mesh_edgemap.adjacencies["Edge2Vert"]
    e2v_he = mesh_halfedge.adjacencies["Edge2Vert"]
    v2e_em = mesh_edgemap.adjacencies["Vert2Edge"]
    v2e_he = mesh_halfedge.adjacencies["Vert2Edge"]

    assert len(v2e_em) == len(v2e_he), "Vertex count mismatch"

    for vert_idx in range(len(v2e_em)):
        em_edges = _normalized_edge_set(e2v_em, v2e_em[vert_idx])
        he_edges = _normalized_edge_set(e2v_he, v2e_he[vert_idx])
        assert em_edges == he_edges, (
            f"Vert2Edge mismatch for vertex {vert_idx}: "
            f"edgemap={em_edges} halfedge={he_edges}"
        )


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
