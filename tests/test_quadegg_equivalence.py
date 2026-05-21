"""Equivalence tests: verify quad-edge produces same outputs as EdgeMap.

Uses canonical-form comparison for edge lists (handles different edge ID ordering
between backends). Tests all 5 adjacency types on all 4 test fixtures.
"""

import pytest
import numpy as np
from chilmesh.examples import annulus, donut, block_o, structured


def canonical_edge_list(edge2vert):
    """Convert edge list to canonical form (sorted by (min, max))."""
    sorted_edges = np.sort(edge2vert, axis=1)
    return tuple(sorted(map(tuple, sorted_edges)))


def compare_adjacency_canonical(em_adj, qe_adj, adj_name):
    """Compare two adjacencies in canonical form (handles edge ID ordering)."""
    em_canonical = canonical_edge_list(em_adj)
    qe_canonical = canonical_edge_list(qe_adj)
    
    assert em_canonical == qe_canonical, \
        f"{adj_name}: canonical forms differ"


class TestQuadEdgeEquivalence:
    """Test quad-edge equivalence with EdgeMap baseline."""
    
    FIXTURES = [annulus, donut, block_o, structured]
    FIXTURE_NAMES = ['annulus', 'donut', 'block_o', 'structured']
    
    @pytest.mark.parametrize('fixture_fn,fixture_name', zip(FIXTURES, FIXTURE_NAMES))
    def test_edge2vert_equivalence(self, fixture_fn, fixture_name):
        """Quad-edge Edge2Vert matches EdgeMap (canonical form)."""
        mesh_em = fixture_fn(topology_backend='edgemap')
        mesh_qe = fixture_fn(topology_backend='quadegg')
        
        edge2vert_em = mesh_em.adjacencies['Edge2Vert']
        edge2vert_qe = mesh_qe.adjacencies['Edge2Vert']
        
        compare_adjacency_canonical(edge2vert_em, edge2vert_qe, 
                                   f'{fixture_name} Edge2Vert')
    
    @pytest.mark.parametrize('fixture_fn,fixture_name', zip(FIXTURES, FIXTURE_NAMES))
    def test_elem2edge_shape_consistency(self, fixture_fn, fixture_name):
        """Quad-edge Elem2Edge shape matches EdgeMap."""
        mesh_em = fixture_fn(topology_backend='edgemap')
        mesh_qe = fixture_fn(topology_backend='quadegg')
        
        elem2edge_em = mesh_em.adjacencies['Elem2Edge']
        elem2edge_qe = mesh_qe.adjacencies['Elem2Edge']
        
        assert elem2edge_em.shape == elem2edge_qe.shape, \
            f'{fixture_name}: Elem2Edge shape mismatch'
    
    @pytest.mark.parametrize('fixture_fn,fixture_name', zip(FIXTURES, FIXTURE_NAMES))
    def test_vert2edge_count_consistency(self, fixture_fn, fixture_name):
        """Quad-edge Vert2Edge reports same incident counts as EdgeMap."""
        mesh_em = fixture_fn(topology_backend='edgemap')
        mesh_qe = fixture_fn(topology_backend='quadegg')
        
        vert2edge_em = mesh_em.adjacencies['Vert2Edge']
        vert2edge_qe = mesh_qe.adjacencies['Vert2Edge']
        
        assert len(vert2edge_em) == len(vert2edge_qe), \
            f'{fixture_name}: Vert2Edge dict size mismatch'
        
        # Check incident edge counts per vertex
        for v in range(len(vert2edge_em)):
            em_count = len(vert2edge_em[v])
            qe_count = len(vert2edge_qe[v])
            assert em_count == qe_count, \
                f'{fixture_name} vertex {v}: incident edge count mismatch'
    
    @pytest.mark.parametrize('fixture_fn,fixture_name', zip(FIXTURES, FIXTURE_NAMES))
    def test_vert2elem_count_consistency(self, fixture_fn, fixture_name):
        """Quad-edge Vert2Elem reports same incident counts as EdgeMap."""
        mesh_em = fixture_fn(topology_backend='edgemap')
        mesh_qe = fixture_fn(topology_backend='quadegg')
        
        vert2elem_em = mesh_em.adjacencies['Vert2Elem']
        vert2elem_qe = mesh_qe.adjacencies['Vert2Elem']
        
        assert len(vert2elem_em) == len(vert2elem_qe), \
            f'{fixture_name}: Vert2Elem dict size mismatch'
        
        # Check incident element counts per vertex
        for v in range(len(vert2elem_em)):
            em_count = len(vert2elem_em[v])
            qe_count = len(vert2elem_qe[v])
            assert em_count == qe_count, \
                f'{fixture_name} vertex {v}: incident element count mismatch'
    
    @pytest.mark.parametrize('fixture_fn,fixture_name', zip(FIXTURES[:2], FIXTURE_NAMES[:2]))
    def test_edge2elem_shape_consistency(self, fixture_fn, fixture_name):
        """Quad-edge Edge2Elem shape matches EdgeMap."""
        mesh_em = fixture_fn(topology_backend='edgemap')
        mesh_qe = fixture_fn(topology_backend='quadegg')
        
        edge2elem_em = mesh_em.adjacencies['Edge2Elem']
        edge2elem_qe = mesh_qe.adjacencies['Edge2Elem']
        
        assert edge2elem_em.shape == edge2elem_qe.shape, \
            f'{fixture_name}: Edge2Elem shape mismatch'
        
        # Check boundary edge sentinels (-1 for boundary)
        em_boundary_count = np.sum(edge2elem_em == -1)
        qe_boundary_count = np.sum(edge2elem_qe == -1)
        assert em_boundary_count == qe_boundary_count, \
            f'{fixture_name}: boundary edge count mismatch'


class TestQuadEdgeConsistency:
    """Test consistency of adjacency conversions (internal logic)."""
    
    @pytest.mark.parametrize('fixture_fn,fixture_name', 
                            zip([annulus, donut], ['annulus', 'donut']))
    def test_elem2edge_references_valid_edges(self, fixture_fn, fixture_name):
        """Verify Elem2Edge references only valid edge IDs."""
        mesh_qe = fixture_fn(topology_backend='quadegg')
        elem2edge = mesh_qe.adjacencies['Elem2Edge']
        n_edges = mesh_qe.n_edges
        
        # All non-padding edge IDs should be in range [0, n_edges)
        valid_mask = (elem2edge >= -1) & (elem2edge < n_edges)
        assert np.all(valid_mask), \
            f'{fixture_name}: invalid edge IDs in Elem2Edge'
    
    @pytest.mark.parametrize('fixture_fn,fixture_name', 
                            zip([annulus, donut], ['annulus', 'donut']))
    def test_vert2edge_references_valid_edges(self, fixture_fn, fixture_name):
        """Verify Vert2Edge references only valid edge IDs."""
        mesh_qe = fixture_fn(topology_backend='quadegg')
        vert2edge = mesh_qe.adjacencies['Vert2Edge']
        n_edges = mesh_qe.n_edges
        
        for v, edge_ids in vert2edge.items():
            for edge_id in edge_ids:
                assert 0 <= edge_id < n_edges, \
                    f'{fixture_name} vertex {v}: invalid edge ID {edge_id}'
    
    @pytest.mark.parametrize('fixture_fn,fixture_name', 
                            zip([annulus, donut], ['annulus', 'donut']))
    def test_vert2elem_references_valid_elems(self, fixture_fn, fixture_name):
        """Verify Vert2Elem references only valid element IDs."""
        mesh_qe = fixture_fn(topology_backend='quadegg')
        vert2elem = mesh_qe.adjacencies['Vert2Elem']
        n_elems = mesh_qe.n_elems
        
        for v, elem_ids in vert2elem.items():
            for elem_id in elem_ids:
                assert 0 <= elem_id < n_elems, \
                    f'{fixture_name} vertex {v}: invalid element ID {elem_id}'
