"""Equivalence tests: Quad-edge adjacency outputs vs. EdgeMap (canonical form).

These tests validate that QuadEdgeTopology produces bit-identical adjacency
outputs compared to the EdgeMap backend when using canonical-form comparison
(edges sorted by min_vert, max_vert).

Strategy: Build the same mesh with both edgemap and quadegg backends, then
compare adjacency structures using canonical-form sorting for edge-ordering
independence.
"""

import pytest
import numpy as np
from chilmesh.examples import annulus, donut, block_o, structured


FIXTURES = [
    ('annulus', annulus),
    ('donut', donut),
    ('block_o', block_o),
    ('structured', structured),
]


def _canonical_edge_set(edge2vert):
    """Convert edge2vert array to canonical set of (min_v, max_v) tuples.

    Handles different edge orderings between backends by normalizing each
    edge to (min_vert, max_vert) form before comparison.
    """
    canonical = set()
    for edge in edge2vert:
        min_v, max_v = min(edge[0], edge[1]), max(edge[0], edge[1])
        canonical.add((int(min_v), int(max_v)))
    return canonical


class TestQuadEdgeEquivalenceEdge2Vert:
    """Validate Edge2Vert (edge endpoint list) equivalence."""

    @pytest.mark.parametrize('fixture_name,fixture_fn', FIXTURES)
    def test_edge2vert_equivalence(self, fixture_name, fixture_fn):
        """Quad-edge Edge2Vert matches EdgeMap in canonical form."""
        # Build meshes with both backends
        mesh_em = fixture_fn(topology_backend='edgemap')
        mesh_qe = fixture_fn(topology_backend='quadegg')

        # Extract Edge2Vert
        e2v_em = mesh_em.adjacencies['Edge2Vert']
        e2v_qe = mesh_qe.adjacencies['Edge2Vert']

        # Convert to canonical sets
        edges_em = _canonical_edge_set(e2v_em)
        edges_qe = _canonical_edge_set(e2v_qe)

        # Compare
        assert edges_em == edges_qe, (
            f"Edge2Vert mismatch on {fixture_name}: "
            f"EdgeMap has {len(edges_em)} edges, quad-edge has {len(edges_qe)} edges"
        )


class TestQuadEdgeEquivalenceElem2Edge:
    """Validate Elem2Edge (element edge list) equivalence."""

    @pytest.mark.parametrize('fixture_name,fixture_fn', FIXTURES)
    def test_elem2edge_edge_sets(self, fixture_name, fixture_fn):
        """Quad-edge Elem2Edge has same edges per element as EdgeMap (set comparison)."""
        mesh_em = fixture_fn(topology_backend='edgemap')
        mesh_qe = fixture_fn(topology_backend='quadegg')

        e2v_em = mesh_em.adjacencies['Edge2Vert']
        e2v_qe = mesh_qe.adjacencies['Edge2Vert']
        e2e_em = mesh_em.adjacencies['Elem2Edge']
        e2e_qe = mesh_qe.adjacencies['Elem2Edge']

        # Map from EdgeMap edge ID to canonical edge tuple
        edges_em_list = [tuple(sorted(e)) for e in e2v_em]
        edges_qe_list = [tuple(sorted(e)) for e in e2v_qe]

        # For each element, collect incident edges in canonical form
        for elem_idx in range(mesh_em.n_elems):
            # EdgeMap edges
            edge_ids_em = e2e_em[elem_idx]
            edge_ids_em = edge_ids_em[edge_ids_em >= 0]  # Filter out padding
            edges_in_elem_em = set(edges_em_list[eid] for eid in edge_ids_em)

            # Quad-edge edges
            edge_ids_qe = e2e_qe[elem_idx]
            edge_ids_qe = edge_ids_qe[edge_ids_qe >= 0]  # Filter out padding
            edges_in_elem_qe = set(edges_qe_list[eid] for eid in edge_ids_qe)

            assert edges_in_elem_em == edges_in_elem_qe, (
                f"Elem2Edge mismatch for {fixture_name} element {elem_idx}: "
                f"EdgeMap edges {edges_in_elem_em} vs quad-edge edges {edges_in_elem_qe}"
            )


class TestQuadEdgeEquivalenceVert2Edge:
    """Validate Vert2Edge (vertex incident edge list) equivalence."""

    @pytest.mark.parametrize('fixture_name,fixture_fn', FIXTURES)
    def test_vert2edge_incident_edges(self, fixture_name, fixture_fn):
        """Quad-edge Vert2Edge has same incident edges per vertex as EdgeMap."""
        mesh_em = fixture_fn(topology_backend='edgemap')
        mesh_qe = fixture_fn(topology_backend='quadegg')

        e2v_em = mesh_em.adjacencies['Edge2Vert']
        e2v_qe = mesh_qe.adjacencies['Edge2Vert']
        v2e_em = mesh_em.adjacencies['Vert2Edge']
        v2e_qe = mesh_qe.adjacencies['Vert2Edge']

        # Map from edge IDs to canonical edge tuples
        edges_em_list = [tuple(sorted(e)) for e in e2v_em]
        edges_qe_list = [tuple(sorted(e)) for e in e2v_qe]

        # For each vertex, check incident edges
        for vert_idx in range(mesh_em.n_verts):
            # EdgeMap incident edges (in canonical form)
            edge_ids_em = v2e_em.get(vert_idx, set())
            incident_em = set(edges_em_list[eid] for eid in edge_ids_em)

            # Quad-edge incident edges (in canonical form)
            edge_ids_qe = v2e_qe.get(vert_idx, set())
            incident_qe = set(edges_qe_list[eid] for eid in edge_ids_qe)

            assert incident_em == incident_qe, (
                f"Vert2Edge mismatch for {fixture_name} vertex {vert_idx}: "
                f"EdgeMap incident {incident_em} vs quad-edge incident {incident_qe}"
            )


class TestQuadEdgeEquivalenceVert2Elem:
    """Validate Vert2Elem (vertex incident element list) equivalence."""

    @pytest.mark.parametrize('fixture_name,fixture_fn', FIXTURES)
    def test_vert2elem_incident_elements(self, fixture_name, fixture_fn):
        """Quad-edge Vert2Elem has same incident elements per vertex as EdgeMap."""
        mesh_em = fixture_fn(topology_backend='edgemap')
        mesh_qe = fixture_fn(topology_backend='quadegg')

        v2m_em = mesh_em.adjacencies['Vert2Elem']
        v2m_qe = mesh_qe.adjacencies['Vert2Elem']

        # For each vertex, check incident elements
        for vert_idx in range(mesh_em.n_verts):
            elems_em = v2m_em.get(vert_idx, set())
            elems_qe = v2m_qe.get(vert_idx, set())

            assert elems_em == elems_qe, (
                f"Vert2Elem mismatch for {fixture_name} vertex {vert_idx}: "
                f"EdgeMap incident {elems_em} vs quad-edge incident {elems_qe}"
            )


class TestQuadEdgeEquivalenceOverallCount:
    """High-level equivalence checks (mesh size invariants)."""

    @pytest.mark.parametrize('fixture_name,fixture_fn', FIXTURES)
    def test_vertex_count_unchanged(self, fixture_name, fixture_fn):
        """Vertex count unchanged across backends."""
        mesh_em = fixture_fn(topology_backend='edgemap')
        mesh_qe = fixture_fn(topology_backend='quadegg')

        assert mesh_em.n_verts == mesh_qe.n_verts, (
            f"Vertex count mismatch for {fixture_name}"
        )

    @pytest.mark.parametrize('fixture_name,fixture_fn', FIXTURES)
    def test_element_count_unchanged(self, fixture_name, fixture_fn):
        """Element count unchanged across backends."""
        mesh_em = fixture_fn(topology_backend='edgemap')
        mesh_qe = fixture_fn(topology_backend='quadegg')

        assert mesh_em.n_elems == mesh_qe.n_elems, (
            f"Element count mismatch for {fixture_name}"
        )

    @pytest.mark.parametrize('fixture_name,fixture_fn', FIXTURES)
    def test_edge_count_unchanged(self, fixture_name, fixture_fn):
        """Edge count unchanged across backends."""
        mesh_em = fixture_fn(topology_backend='edgemap')
        mesh_qe = fixture_fn(topology_backend='quadegg')

        assert mesh_em.n_edges == mesh_qe.n_edges, (
            f"Edge count mismatch for {fixture_name}: "
            f"EdgeMap {mesh_em.n_edges} vs quad-edge {mesh_qe.n_edges}"
        )
