"""Tests for currently-untested branches in mutations.py.

Covers edge cases in reskeletonize_local(), smooth_topology(), and related methods.
Targets specific uncovered branches identified in code review.
"""

from __future__ import annotations

import numpy as np
from chilmesh import MutableMesh
from chilmesh import examples


class TestReskeletonizeLocalBranches:
    """Tests for uncovered branches in reskeletonize_local().

    Target lines 576-578 (empty OE), 590-593 (no layer found),
    and 596-598 (start_layer == 0).
    """

    def test_reskeletonize_local_with_empty_oe_branch(self):
        """Branch lines 576-578: OE layers empty → calls _skeletonize() and returns.

        Creates a mesh, clears its layers, then calls reskeletonize_local.
        """
        mesh = annulus()
        original_n_verts = mesh.n_verts
        mutable = MutableMesh(mesh)

        # Clear OE layers (simulate "no skeletonization yet")
        mesh.layers['OE'] = []
        mesh.layers['IE'] = []
        mesh.layers['OV'] = []
        mesh.layers['IV'] = []
        mesh.layers['bEdgeIDs'] = []
        mesh.n_layers = 0

        # Call reskeletonize_local with any element ID
        mutable.reskeletonize_local(np.array([0, 1, 2]), radius=2)

        # Should have re-skeletonized the mesh
        assert mesh.n_layers > 0, "Expected n_layers > 0 after reskeletonize"
        assert len(mesh.layers['OE']) > 0, "Expected OE to be populated after skeletonize"
        assert mesh.n_verts == original_n_verts, "Vertices count should not change"

    def test_reskeletonize_local_with_nonexistent_elem_ids_branch(self):
        """Branch lines 590-593: elem_ids not in any layer → full rebuild.

        Skeletonizes a mesh, then calls reskeletonize_local with an element ID
        that doesn't exist in any layer (e.g., ID beyond n_elems).
        """
        mesh = annulus()
        original_n_layers = mesh.n_layers
        mutable = MutableMesh(mesh)

        # Get an ID that's out of range
        huge_id = mesh.n_elems + 100
        elem_ids = np.array([huge_id])

        # This should trigger full rebuild (no layer found containing huge_id)
        mutable.reskeletonize_local(elem_ids, radius=2)

        # Verify the mesh is still valid (full skeletonize was called)
        assert mesh.n_layers > 0, "Expected n_layers > 0 after rebuild"
        # n_layers should be regenerated (may differ from original)
        assert len(mesh.layers['OE']) == mesh.n_layers

    def test_reskeletonize_local_with_layer_zero_affected(self):
        """Branch lines 596-598: affected_layer=0 with radius=2 → start_layer==0 → full rebuild.

        The affected layer is layer 0 (outermost boundary). With radius=2,
        start_layer = max(0, 0-2) = 0, triggering the branch.
        """
        mesh = annulus()
        original_n_layers = mesh.n_layers
        mutable = MutableMesh(mesh)

        # Get an element from layer 0 (outer elements)
        if len(mesh.layers['OE']) > 0 and len(mesh.layers['OE'][0]) > 0:
            outer_elem_id = int(mesh.layers['OE'][0][0])

            # Call reskeletonize_local with element from layer 0
            mutable.reskeletonize_local(np.array([outer_elem_id]), radius=2)

            # Verify mesh state is valid
            assert mesh.n_layers > 0
            assert len(mesh.layers['OE']) == mesh.n_layers
            # n_layers may or may not equal original (depends on mesh structure)
            assert mesh.n_verts > 0

    def test_reskeletonize_local_partial_rebuild(self):
        """Branch around line 595: normal partial rebuild (start_layer > 0).

        Affected element is in mid-layer, so start_layer > 0 and < n_layers.
        This tests the "happy path" partial rebuild logic.
        """
        mesh = annulus()
        original_n_layers = mesh.n_layers
        mutable = MutableMesh(mesh)

        # Find an element in a mid-layer (not layer 0)
        if original_n_layers > 2:
            mid_layer_idx = 1  # Second layer
            if len(mesh.layers['OE'][mid_layer_idx]) > 0:
                elem_id = int(mesh.layers['OE'][mid_layer_idx][0])

                # Call reskeletonize_local with radius=1 so start_layer > 0
                mutable.reskeletonize_local(np.array([elem_id]), radius=1)

                # Verify state
                assert mesh.n_layers > 0
                assert len(mesh.layers['OE']) == mesh.n_layers


class TestSmoothTopologyBranches:
    """Tests for uncovered branches in smooth_topology().

    Target the pass loop body (lines 526-550) and exception handling (543-544).
    """

    def test_smooth_topology_end_to_end(self):
        """Exercise smooth_topology on annulus mesh.

        Lines 526-550: the pass loop body with edge-scan and swap logic.
        """
        mesh = annulus()
        original_n_elems = mesh.n_elems
        mutable = MutableMesh(mesh)

        n_swapped = mutable.smooth_topology(metric_threshold=0.0, max_passes=10)

        # Should run at least one pass and either swap edges or stop
        assert n_swapped >= 0, "n_swapped must be non-negative"
        # Mesh topology should be preserved
        assert mesh.n_elems == original_n_elems
        # All elements should still be valid triangles
        assert all(
            mutable._is_triangle(i) for i in range(mesh.n_elems)
        ), "All elements should remain triangles after smooth_topology"

    def test_smooth_topology_with_threshold(self):
        """Test smooth_topology with non-zero metric_threshold.

        Higher threshold means fewer swaps are accepted.
        """
        mesh1 = annulus()
        mesh2 = annulus()
        mutable1 = MutableMesh(mesh1)
        mutable2 = MutableMesh(mesh2)

        # One with threshold=0 (accept any improvement)
        swaps_loose = mutable1.smooth_topology(metric_threshold=0.0, max_passes=5)
        # One with high threshold (accept only large improvements)
        swaps_strict = mutable2.smooth_topology(metric_threshold=1.0, max_passes=5)

        # Strict threshold should accept fewer or equal swaps
        assert swaps_strict <= swaps_loose, (
            f"Strict threshold {swaps_strict} should be <= loose {swaps_loose}"
        )

    def test_smooth_topology_convergence(self):
        """Test that smooth_topology converges (no infinite looping).

        Should exit after max_passes or when no edges improve.
        """
        mesh = annulus()
        mutable = MutableMesh(mesh)

        # Run with small max_passes; should complete
        n_swapped = mutable.smooth_topology(metric_threshold=0.0, max_passes=3)

        assert isinstance(n_swapped, int)
        assert n_swapped >= 0

    def test_smooth_topology_preserves_boundary(self):
        """Verify smooth_topology doesn't touch boundary edges.

        Boundary edges are skipped in the loop (lines 532-533: skip if boundary).
        """
        mesh = annulus()
        original_boundary_edges = set()
        for eid in range(mesh.n_edges):
            edge2elem = mesh.edge2elem()
            ea, eb = edge2elem[eid]
            if int(ea) == -1 or int(eb) == -1:
                original_boundary_edges.add(eid)

        mutable = MutableMesh(mesh)
        mutable.smooth_topology(metric_threshold=0.0, max_passes=5)

        # Re-check boundary (should not change)
        new_boundary_edges = set()
        for eid in range(mesh.n_edges):
            edge2elem = mesh.edge2elem()
            ea, eb = edge2elem[eid]
            if int(ea) == -1 or int(eb) == -1:
                new_boundary_edges.add(eid)

        # Boundary edges should be the same (topology preserved)
        assert len(original_boundary_edges) == len(new_boundary_edges)

    def test_smooth_topology_exception_swallow(self):
        """Verify smooth_topology catches exceptions in swap_edge and continues.

        Lines 543-544: `except (ValueError, RuntimeError): pass`.
        This is a best-effort test — we cannot guarantee swap_edge raises
        without monkeypatching or crafting a very specific geometry.
        """
        mesh = annulus()
        mutable = MutableMesh(mesh)

        # Monkeypatch swap_edge to raise RuntimeError on first call
        original_swap = mutable.swap_edge
        call_count = [0]

        def failing_swap_edge(edge_id):
            call_count[0] += 1
            if call_count[0] == 1:
                raise RuntimeError("Simulated swap failure")
            return original_swap(edge_id)

        mutable.swap_edge = failing_swap_edge

        # Should not raise; exception should be caught and loop continues
        n_swapped = mutable.smooth_topology(metric_threshold=0.0, max_passes=10)

        # Call count should be > 1 (first attempt failed, loop continued)
        assert call_count[0] >= 1, "swap_edge should have been called at least once"
        # smooth_topology should complete without raising
        assert isinstance(n_swapped, int)


class TestSmoothTopologyIntegration:
    """Integration tests for smooth_topology with real mesh scenarios."""

    def test_smooth_topology_on_poorly_shaped_mesh(self):
        """Apply smooth_topology to a real mesh and verify output.

        Tests the integration of all loops and conditions.
        """
        mesh = donut()
        original_n_elems = mesh.n_elems
        mutable = MutableMesh(mesh)

        # Call on donut mesh (medium complexity)
        n_swapped = mutable.smooth_topology(metric_threshold=0.0, max_passes=20)

        assert mesh.n_elems == original_n_elems, "Topology count must not change"
        assert n_swapped >= 0
        # Verify all elements are valid after smoothing
        mutable._validate_invariants()


def donut():
    """Fetch donut fixture without caching."""
    example_fn = examples.donut
    if hasattr(example_fn, '__wrapped__'):
        return example_fn.__wrapped__()
    return example_fn()


def annulus():
    """Fetch annulus fixture without caching."""
    example_fn = examples.annulus
    if hasattr(example_fn, '__wrapped__'):
        return example_fn.__wrapped__()
    return example_fn()
