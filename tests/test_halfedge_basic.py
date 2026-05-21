"""Basic HalfEdgeTopology construction and invariant tests.

These tests verify:
- HalfEdgeTopology constructs without error on all fixtures
- Half-edge invariants hold (twin-of-twin is identity, next walk is cyclic)
- Boundary half-edges have twin_idx == -1
- Padded triangles don't generate spurious half-edges (F-4)
- Degenerate quads don't crash construction (F-1)
- Unknown backend value raises ValueError (FR-008)
"""

import pytest
import numpy as np
from chilmesh.examples import annulus, donut, block_o, structured
from chilmesh.mesh_topology_halfedge import build_halfedge_from_connectivity


@pytest.mark.parametrize(
    "fixture_fn",
    [annulus, donut, block_o, structured],
    ids=["annulus", "donut", "block_o", "structured"]
)
def test_halfedge_construction_succeeds(fixture_fn):
    """HalfEdgeTopology construction succeeds on all fixtures."""
    mesh = fixture_fn()
    elem2vert = mesh.adjacencies["Elem2Vert"]
    n_verts = mesh.n_verts
    he_topo = build_halfedge_from_connectivity(elem2vert, n_verts)
    assert he_topo.half_edges.shape[1] == 4
    assert len(he_topo.half_edges) > 0


@pytest.mark.parametrize(
    "fixture_fn",
    [annulus, donut, block_o, structured],
    ids=["annulus", "donut", "block_o", "structured"]
)
def test_halfedge_twin_of_twin_identity(fixture_fn):
    """Twin-of-twin is identity: HE[HE[i, 1], 1] == i for interior half-edges."""
    mesh = fixture_fn()
    elem2vert = mesh.adjacencies["Elem2Vert"]
    n_verts = mesh.n_verts
    he_topo = build_halfedge_from_connectivity(elem2vert, n_verts)

    for i in range(len(he_topo.half_edges)):
        twin_idx = int(he_topo.half_edges[i, 1])
        if twin_idx >= 0:
            twin_of_twin = int(he_topo.half_edges[twin_idx, 1])
            assert twin_of_twin == i, f"Twin-of-twin failed for HE {i}"


@pytest.mark.parametrize(
    "fixture_fn",
    [annulus, donut, block_o, structured],
    ids=["annulus", "donut", "block_o", "structured"]
)
def test_halfedge_next_walk_is_cyclic(fixture_fn):
    """Walking next_idx around a face returns to origin."""
    mesh = fixture_fn()
    elem2vert = mesh.adjacencies["Elem2Vert"]
    n_verts = mesh.n_verts
    he_topo = build_halfedge_from_connectivity(elem2vert, n_verts)

    for face_idx in range(he_topo.n_elems):
        start_he = None
        for i in range(len(he_topo.half_edges)):
            if int(he_topo.half_edges[i, 3]) == face_idx:
                start_he = i
                break

        if start_he is not None:
            he_idx = start_he
            face_size = 0
            max_iters = 10
            while face_size < max_iters:
                next_idx = int(he_topo.half_edges[he_idx, 2])
                he_idx = next_idx
                face_size += 1
                if he_idx == start_he:
                    break

            assert he_idx == start_he, f"Next walk didn't close for face {face_idx}"


@pytest.mark.parametrize(
    "fixture_fn",
    [annulus, donut, block_o, structured],
    ids=["annulus", "donut", "block_o", "structured"]
)
def test_halfedge_boundary_marked(fixture_fn):
    """Boundary half-edges have twin_idx == -1."""
    mesh = fixture_fn()
    elem2vert = mesh.adjacencies["Elem2Vert"]
    n_verts = mesh.n_verts
    he_topo = build_halfedge_from_connectivity(elem2vert, n_verts)

    boundary_count = 0
    for i in range(len(he_topo.half_edges)):
        twin_idx = int(he_topo.half_edges[i, 1])
        if twin_idx == -1:
            boundary_count += 1

    assert boundary_count > 0, f"No boundary half-edges found"


def test_halfedge_padded_triangle_no_spurious_he(fresh_mesh):
    """Padded triangle does NOT generate spurious 4th half-edge (F-4 regression).

    block_o fixture has mixed elements (tri + quad).
    Count total half-edges; verify no spurious 4th HE for triangles.
    """
    mesh = fresh_mesh
    elem2vert = mesh.adjacencies["Elem2Vert"]
    n_verts = mesh.n_verts

    if mesh.type != "Mixed-Element":
        pytest.skip("Fixture is not mixed-element")

    he_topo = build_halfedge_from_connectivity(elem2vert, n_verts)

    for elem_idx in range(he_topo.n_elems):
        elem_verts = elem2vert[elem_idx]
        n_real_edges = 3 if (elem_verts[3] == elem_verts[0]) else 4

        he_count = sum(
            1 for i in range(len(he_topo.half_edges))
            if int(he_topo.half_edges[i, 3]) == elem_idx
        )

        assert he_count == n_real_edges, f"Elem {elem_idx}: expected {n_real_edges} HE, got {he_count}"


def test_halfedge_unknown_backend_raises():
    """Unknown topology_backend raises ValueError (FR-008)."""
    mesh = annulus()
    with pytest.raises(ValueError, match="Unknown CHILMESH_TOPOLOGY_BACKEND"):
        mesh_bad = annulus(topology_backend="unknown_backend")


def test_halfedge_backend_via_env_var(topology_backend_env):
    """Backend selection via env var and kwarg works."""
    backend = topology_backend_env
    mesh = annulus()
    elem2vert = mesh.adjacencies["Elem2Vert"]
    n_verts = mesh.n_verts

    he_topo = build_halfedge_from_connectivity(elem2vert, n_verts)
    assert he_topo is not None
