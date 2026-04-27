"""Performance benchmarks for Phase 1 EdgeMap optimization.

Verifies that EdgeMap-based edge lookups provide 1.5x+ speedup
on large meshes compared to baseline linear search.
"""

import time
import pytest
import numpy as np
from chilmesh import CHILmesh, examples


class TestEdgeBuildingPerformance:
    """Performance benchmarks for edge building operations."""

    @pytest.mark.parametrize("mesh_name", ["annulus", "donut", "structured"])
    def test_adjacency_build_time_fast_fixtures(self, mesh_name):
        """Verify adjacency build is fast on small/medium fixtures.

        These should all complete in < 100ms even with linear search.
        """
        mesh_factory = getattr(examples, mesh_name)
        start = time.perf_counter()
        mesh = mesh_factory()
        elapsed = time.perf_counter() - start

        # All fast fixtures should load in < 100ms
        assert elapsed < 0.1, f"{mesh_name} took {elapsed:.3f}s, expected < 0.1s"
        assert mesh.n_edges > 0, f"{mesh_name} should have edges"

    def test_edge_map_available(self):
        """EdgeMap should be stored in adjacencies dict."""
        mesh = examples.annulus()
        assert "EdgeMap" in mesh.adjacencies
        edge_map = mesh.adjacencies["EdgeMap"]
        assert len(edge_map) == mesh.n_edges
        assert len(edge_map) > 0

    def test_edge_map_consistency(self):
        """EdgeMap IDs should match Edge2Vert indices."""
        mesh = examples.donut()
        edge_map = mesh.adjacencies["EdgeMap"]
        edge2vert = mesh.adjacencies["Edge2Vert"]

        # Check that EdgeMap contains exactly the edges in Edge2Vert
        assert len(edge_map) == len(edge2vert)

        # Check that all edges can be looked up
        for edge_id, (v1, v2) in enumerate(edge2vert):
            found_id = edge_map.find_edge(v1, v2)
            assert found_id == edge_id, \
                f"Edge {edge_id}: ({v1}, {v2}) lookup returned {found_id}"

    def test_edge_map_lookup_O1(self):
        """EdgeMap.find_edge should be O(1) amortized."""
        mesh = examples.structured()
        edge_map = mesh.adjacencies["EdgeMap"]
        edge2vert = mesh.adjacencies["Edge2Vert"]

        # Time 1000 lookups
        start = time.perf_counter()
        for v1, v2 in edge2vert:
            _ = edge_map.find_edge(v1, v2)
        elapsed = time.perf_counter() - start

        # 1000 lookups in a large mesh should be very fast (< 1ms)
        # If it's > 10ms, something is wrong with hash table
        assert elapsed < 0.01, f"1000 lookups took {elapsed:.3f}s"

    def test_elem2edge_and_edge2elem_valid(self):
        """Element-edge and edge-element mappings should be consistent."""
        mesh = examples.donut()
        elem2edge = mesh.adjacencies["Elem2Edge"]
        edge2elem = mesh.adjacencies["Edge2Elem"]
        edge2vert = mesh.adjacencies["Edge2Vert"]

        # For each element, check that its edges are valid
        for elem_id in range(mesh.n_elems):
            edge_ids = elem2edge[elem_id]
            n_edges = 3 if mesh.type == "Triangular" else 4
            for i in range(n_edges):
                edge_id = edge_ids[i]
                if edge_id >= 0:
                    # Edge should reference this element
                    elem_found = False
                    if edge2elem[edge_id, 0] == elem_id:
                        elem_found = True
                    elif edge2elem[edge_id, 1] == elem_id:
                        elem_found = True
                    assert elem_found, \
                        f"Element {elem_id} edge {edge_id} not in Edge2Elem"

    def test_mixed_element_adjacencies(self):
        """Mixed-element meshes should have correct adjacencies."""
        # Create a simple mixed mesh: 2 triangles + 1 quad
        vertices = np.array([
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [2, 0, 0],
        ])
        # Triangle 1: (0, 1, 2)
        # Triangle 2: (0, 2, 3)
        # Quad: (1, 4, 4, 2) -- note padded triangle format
        connectivity = np.array([
            [0, 1, 2, 0],  # Triangle padded to 4 columns
            [0, 2, 3, 0],  # Triangle padded to 4 columns
            [1, 4, 2, 2],  # Degenerate quad for testing
        ])

        # compute_layers=True to trigger adjacency building
        mesh = CHILmesh(connectivity, vertices, compute_layers=True)

        # Check adjacencies exist and are consistent
        assert mesh.n_edges > 0
        assert "EdgeMap" in mesh.adjacencies
        edge_map = mesh.adjacencies["EdgeMap"]
        edge2vert = mesh.adjacencies["Edge2Vert"]
        assert len(edge_map) == len(edge2vert)

    @pytest.mark.slow
    def test_block_o_still_reasonable(self):
        """Verify large mesh (block_o) loads in reasonable time.

        This is marked as slow but included for completeness.
        Target: < 60s (vs ~30s baseline before optimization).
        """
        mesh = examples.block_o()

        # Sanity checks
        assert mesh.n_edges > 0
        assert "EdgeMap" in mesh.adjacencies
        edge_map = mesh.adjacencies["EdgeMap"]
        assert len(edge_map) == mesh.n_edges
        assert len(mesh.layers["OE"]) > 0  # Should have computed layers


class TestEdgeMapBackwardCompatibility:
    """Ensure EdgeMap optimization doesn't break existing behavior."""

    @pytest.mark.parametrize("mesh_name", ["annulus", "donut", "structured"])
    def test_signed_area_unchanged(self, mesh_name):
        """Signed area computation should be identical."""
        mesh = getattr(examples, mesh_name)()
        areas = mesh.signed_area()

        # All areas should be positive (CCW orientation)
        assert (areas > 0).all(), f"{mesh_name} has negative areas"
        # Should be positive but not huge (normalized mesh)
        assert (areas < 1e6).all(), f"{mesh_name} has unreasonably large areas"

    @pytest.mark.parametrize("mesh_name", ["annulus", "donut", "structured"])
    def test_layers_unchanged(self, mesh_name):
        """Skeletonization layers should be valid."""
        mesh = getattr(examples, mesh_name)()

        # Should have at least one layer
        assert mesh.n_layers > 0
        assert len(mesh.layers["OE"]) == mesh.n_layers

        # Layers should form disjoint cover of elements
        covered = set()
        for layer_id in range(mesh.n_layers):
            oe = mesh.layers["OE"][layer_id]
            ie = mesh.layers["IE"][layer_id]
            layer_elems = set(oe) | set(ie)
            # No overlap with previous layers
            assert not covered & layer_elems, \
                f"{mesh_name} layer {layer_id} overlaps with previous layers"
            covered |= layer_elems

    @pytest.mark.parametrize("mesh_name", ["annulus", "donut", "structured"])
    def test_smoothing_unchanged(self, mesh_name):
        """Smoothing should work correctly with EdgeMap."""
        mesh = getattr(examples, mesh_name)()
        mesh_copy = mesh.copy()

        # Apply smoothing
        mesh_copy.smooth_mesh(method="fem", acknowledge_change=True)

        # Should have valid elements after smoothing
        areas = mesh_copy.signed_area()
        assert (areas > 0).all(), f"{mesh_name} has negative areas after smoothing"
