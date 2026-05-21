"""Basic construction and invariant tests for QuadEdgeTopology.

These tests validate the quad-edge data structure without relying on
equivalence with other backends. Focuses on construction algorithm
invariants and boundary handling.
"""

import pytest
import numpy as np
from chilmesh.mesh_topology_quadegg import QuadEdgeTopology, build_quadegg_from_connectivity
from chilmesh.examples import annulus, donut, block_o, structured


@pytest.fixture(params=['annulus', 'donut', 'block_o', 'structured'])
def quad_edge_fixture(request):
    """Parametrized fixture: build quad-edge for each test mesh."""
    mesh_name = request.param
    if mesh_name == 'annulus':
        mesh = annulus()
    elif mesh_name == 'donut':
        mesh = donut()
    elif mesh_name == 'block_o':
        mesh = block_o()
    else:  # structured
        mesh = structured()
    
    elem2vert = mesh.adjacencies['Elem2Vert']
    n_verts = elem2vert.max() + 1
    qe = build_quadegg_from_connectivity(elem2vert, n_verts)
    return mesh, qe, mesh_name


class TestQuadEdgeConstruction:
    """Test quad-edge construction algorithm."""
    
    def test_quadegg_shape(self, quad_edge_fixture):
        """Test that edges array has correct shape [n_edges, 4]."""
        mesh, qe, name = quad_edge_fixture
        assert qe.edges.ndim == 2, f"{name}: edges should be 2D"
        assert qe.edges.shape[1] == 4, f"{name}: edges should have 4 columns"
        assert qe.edges.shape[0] > 0, f"{name}: should have at least 1 edge"
    
    def test_quadegg_dtype(self, quad_edge_fixture):
        """Test that edges array has int32 dtype."""
        mesh, qe, name = quad_edge_fixture
        assert qe.edges.dtype == np.int32, f"{name}: edges should be int32"
    
    def test_boundary_sentinels(self, quad_edge_fixture):
        """Test that boundary edges have opposite_idx = -1."""
        mesh, qe, name = quad_edge_fixture
        # Some edges should be boundary (opposite_idx = -1)
        boundary_count = np.sum(qe.edges[:, 3] == -1)
        assert boundary_count > 0, f"{name}: should have boundary edges"
    
    def test_next_cw_validity(self, quad_edge_fixture):
        """Test that next_cw points to valid edges or -1."""
        mesh, qe, name = quad_edge_fixture
        next_cw = qe.edges[:, 1]
        valid_mask = (next_cw >= -1) & (next_cw < len(qe.edges))
        assert np.all(valid_mask), f"{name}: invalid next_cw indices found"
    
    def test_opposite_idx_validity(self, quad_edge_fixture):
        """Test that opposite_idx points to valid edges or -1."""
        mesh, qe, name = quad_edge_fixture
        opposite_idx = qe.edges[:, 3]
        valid_mask = (opposite_idx >= -1) & (opposite_idx < len(qe.edges))
        assert np.all(valid_mask), f"{name}: invalid opposite_idx indices found"
    
    def test_origin_vertex_validity(self, quad_edge_fixture):
        """Test that origin_vertex points to valid vertices."""
        mesh, qe, name = quad_edge_fixture
        origin = qe.edges[:, 0]
        valid_mask = (origin >= 0) & (origin < qe.n_verts)
        assert np.all(valid_mask), f"{name}: invalid origin vertices found"
    
    def test_opposite_symmetry(self, quad_edge_fixture):
        """Test that opposite edges point to each other."""
        mesh, qe, name = quad_edge_fixture
        for i in range(len(qe.edges)):
            opposite_idx = int(qe.edges[i, 3])
            if opposite_idx >= 0:
                # If i has an opposite, that opposite should point back to i
                reverse_opposite = int(qe.edges[opposite_idx, 3])
                assert reverse_opposite == i, \
                    f"{name}: opposite edge symmetry broken at edge {i}"
    
    def test_edge_count_reasonable(self, quad_edge_fixture):
        """Test that edge count is reasonable (2x undirected edge count)."""
        mesh, qe, name = quad_edge_fixture
        # Quad-edge has 2 directed edges per undirected edge
        undirected_count = len(qe.to_edge2vert())
        directed_count = len(qe.edges)
        assert directed_count >= undirected_count, \
            f"{name}: directed edges ({directed_count}) < undirected ({undirected_count})"
        assert directed_count <= 2 * undirected_count + 1, \
            f"{name}: directed edges ({directed_count}) too many vs undirected ({undirected_count})"
    
    def test_elem_count_preserved(self, quad_edge_fixture):
        """Test that element count is preserved."""
        mesh, qe, name = quad_edge_fixture
        elem2vert = mesh.adjacencies['Elem2Vert']
        assert qe.n_elems == elem2vert.shape[0], \
            f"{name}: element count mismatch"
    
    def test_conversion_methods_exist(self, quad_edge_fixture):
        """Test that all required conversion methods are available."""
        mesh, qe, name = quad_edge_fixture
        assert hasattr(qe, 'to_edge2vert'), f"{name}: missing to_edge2vert"
        assert hasattr(qe, 'to_elem2edge'), f"{name}: missing to_elem2edge"
        assert hasattr(qe, 'to_edge2elem'), f"{name}: missing to_edge2elem"
        assert hasattr(qe, 'to_vert2edge'), f"{name}: missing to_vert2edge"
        assert hasattr(qe, 'to_vert2elem'), f"{name}: missing to_vert2elem"
        assert hasattr(qe, 'to_edgemap_list'), f"{name}: missing to_edgemap_list"


class TestQuadEdgeConversions:
    """Test quad-edge adjacency conversions."""
    
    def test_edge2vert_shape(self, quad_edge_fixture):
        """Test Edge2Vert has correct shape [n_edges, 2]."""
        mesh, qe, name = quad_edge_fixture
        edge2vert = qe.to_edge2vert()
        assert edge2vert.ndim == 2, f"{name}: Edge2Vert should be 2D"
        assert edge2vert.shape[1] == 2, f"{name}: Edge2Vert should have 2 columns"
    
    def test_elem2edge_shape(self, quad_edge_fixture):
        """Test Elem2Edge has correct shape [n_elems, elem_type]."""
        mesh, qe, name = quad_edge_fixture
        elem2edge = qe.to_elem2edge()
        assert elem2edge.ndim == 2, f"{name}: Elem2Edge should be 2D"
        assert elem2edge.shape[0] == qe.n_elems, f"{name}: Elem2Edge row count mismatch"
    
    def test_edge2elem_shape(self, quad_edge_fixture):
        """Test Edge2Elem has correct shape [n_edges, 2]."""
        mesh, qe, name = quad_edge_fixture
        edge2elem = qe.to_edge2elem()
        assert edge2elem.ndim == 2, f"{name}: Edge2Elem should be 2D"
        assert edge2elem.shape[1] == 2, f"{name}: Edge2Elem should have 2 columns"
    
    def test_vert2edge_dict(self, quad_edge_fixture):
        """Test Vert2Edge is a dict with all vertices."""
        mesh, qe, name = quad_edge_fixture
        vert2edge = qe.to_vert2edge()
        assert isinstance(vert2edge, dict), f"{name}: Vert2Edge should be a dict"
        assert len(vert2edge) == qe.n_verts, f"{name}: Vert2Edge should have all vertices"
    
    def test_vert2elem_dict(self, quad_edge_fixture):
        """Test Vert2Elem is a dict with all vertices."""
        mesh, qe, name = quad_edge_fixture
        vert2elem = qe.to_vert2elem()
        assert isinstance(vert2elem, dict), f"{name}: Vert2Elem should be a dict"
        assert len(vert2elem) == qe.n_verts, f"{name}: Vert2Elem should have all vertices"
    
    def test_edgemap_list_format(self, quad_edge_fixture):
        """Test edgemap_list returns sorted canonical edges."""
        mesh, qe, name = quad_edge_fixture
        edge_list = qe.to_edgemap_list()
        assert isinstance(edge_list, list), f"{name}: edgemap_list should be a list"
        assert all(isinstance(e, tuple) for e in edge_list), \
            f"{name}: all edges should be tuples"
        assert all(len(e) == 2 for e in edge_list), \
            f"{name}: all edges should have 2 vertices"
        # Check sorted
        assert edge_list == sorted(edge_list), \
            f"{name}: edges should be sorted"


class TestQuadEdgeEquivalence:
    """Test quad-edge equivalence with EdgeMap on basic properties."""
    
    def test_edge_count_matches_edgemap(self, quad_edge_fixture):
        """Test that quad-edge produces same edge count as EdgeMap."""
        mesh_qe, qe, name = quad_edge_fixture
        mesh_em = annulus(topology_backend='edgemap') if name == 'annulus' else \
                  donut(topology_backend='edgemap') if name == 'donut' else \
                  block_o(topology_backend='edgemap') if name == 'block_o' else \
                  structured(topology_backend='edgemap')
        
        edge2vert_qe = qe.to_edge2vert()
        edge2vert_em = mesh_em.adjacencies['Edge2Vert']
        
        assert len(edge2vert_qe) == len(edge2vert_em), \
            f"{name}: edge count mismatch: QE={len(edge2vert_qe)} vs EM={len(edge2vert_em)}"
    
    def test_elem_count_matches_input(self, quad_edge_fixture):
        """Test that Elem2Edge reports same element count as input."""
        mesh_qe, qe, name = quad_edge_fixture
        elem2vert = mesh_qe.adjacencies['Elem2Vert']
        elem2edge = qe.to_elem2edge()
        
        assert elem2edge.shape[0] == elem2vert.shape[0], \
            f"{name}: element count mismatch"
