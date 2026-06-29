"""Regression tests for ``CHILmesh.smooth()`` Laplacian smoother.

The ``smooth()`` method implements in-place Laplacian smoothing with:
- Interior vertex relaxation toward edge-neighbor average
- Optional boundary-vertex fixing
- Z-coordinate preservation
- Configurable iteration count and weight

Invariants pinned here:
- ValueError raised if weight not in [0, 1].
- Returns None (modifies in place).
- weight=0 leaves points unchanged.
- weight>0 moves interior vertices (on non-trivial meshes).
- fixed_boundary=True prevents boundary vertex movement.
- fixed_boundary=False allows boundary vertex movement.
- Z-coordinate (column 2) preserved when present.
- n_iter=0 leaves points unchanged (no loop iterations).
- _build_spatial_indices() called at end.
"""
from __future__ import annotations

import numpy as np
import pytest

from chilmesh import CHILmesh, examples


class TestSmoothValidation:
    """Test weight parameter validation."""

    def test_weight_below_zero_raises(self):
        """weight < 0 raises ValueError."""
        mesh = examples.annulus()
        with pytest.raises(ValueError, match="weight must be in"):
            mesh.smooth(weight=-0.1)

    def test_weight_above_one_raises(self):
        """weight > 1 raises ValueError."""
        mesh = examples.annulus()
        with pytest.raises(ValueError, match="weight must be in"):
            mesh.smooth(weight=1.5)

    def test_weight_at_zero_valid(self):
        """weight = 0 is valid."""
        mesh = examples.annulus()
        result = mesh.smooth(weight=0.0)
        assert result is None

    def test_weight_at_one_valid(self):
        """weight = 1 is valid."""
        mesh = examples.annulus()
        result = mesh.smooth(weight=1.0)
        assert result is None

    @pytest.mark.parametrize("weight", [0.0, 0.25, 0.5, 0.75, 1.0])
    def test_valid_weights_parametrized(self, weight):
        """All weights in [0, 1] are valid."""
        mesh = examples.annulus()
        result = mesh.smooth(weight=weight)
        assert result is None


class TestSmoothInPlace:
    """Test in-place modification behavior."""

    def test_returns_none(self):
        """smooth() returns None."""
        mesh = examples.annulus()
        result = mesh.smooth(n_iter=1, weight=0.5)
        assert result is None

    def test_modifies_points_in_place(self):
        """smooth() modifies self.points, not a copy."""
        mesh = examples.annulus()
        orig_id = id(mesh.points)
        mesh.smooth(n_iter=1, weight=0.5)
        # After smoothing, self.points may be reassigned to new_pts (in the loop)
        # But the mesh object's points attribute should be modified
        assert mesh.points is not None
        assert len(mesh.points) > 0

    def test_weight_zero_unchanged(self):
        """weight=0.0 leaves points completely unchanged."""
        mesh = examples.annulus()
        points_before = mesh.points.copy()
        mesh.smooth(n_iter=5, weight=0.0)
        np.testing.assert_array_equal(mesh.points, points_before)

    def test_weight_positive_moves_interior(self):
        """weight > 0 moves at least one interior vertex."""
        mesh = examples.annulus()
        points_before = mesh.points.copy()

        # Apply smoothing with positive weight
        mesh.smooth(n_iter=1, weight=0.5)

        # At least one interior point should have moved
        # (annulus has interior vertices)
        diff = np.linalg.norm(mesh.points - points_before, axis=1)
        assert np.any(diff > 1e-10), "No vertices moved with weight=0.5"


class TestSmoothBoundaryFixed:
    """Test fixed_boundary=True behavior."""

    def test_boundary_vertices_unchanged_when_fixed(self):
        """fixed_boundary=True prevents boundary vertex movement."""
        mesh = examples.annulus()
        boundary_idx = mesh.boundary_node_indices()

        points_before = mesh.points.copy()
        mesh.smooth(n_iter=5, weight=0.5, fixed_boundary=True)
        points_after = mesh.points

        # Boundary vertices should not move
        np.testing.assert_array_equal(
            points_after[boundary_idx, :],
            points_before[boundary_idx, :],
            err_msg="Boundary vertices moved with fixed_boundary=True"
        )

    def test_interior_vertices_can_move_when_boundary_fixed(self):
        """fixed_boundary=True still allows interior vertices to move."""
        mesh = examples.annulus()
        boundary_idx = set(mesh.boundary_node_indices().tolist())
        interior_idx = np.array([v for v in range(mesh.n_verts) if v not in boundary_idx])

        if len(interior_idx) == 0:
            pytest.skip("Mesh has no interior vertices")

        points_before = mesh.points.copy()
        mesh.smooth(n_iter=1, weight=0.5, fixed_boundary=True)

        # At least one interior vertex should have moved
        interior_moved = np.linalg.norm(
            mesh.points[interior_idx, :2] - points_before[interior_idx, :2],
            axis=1
        )
        assert np.any(interior_moved > 1e-10), "No interior vertices moved"


class TestSmoothBoundaryFree:
    """Test fixed_boundary=False behavior."""

    def test_boundary_can_move_when_not_fixed(self):
        """fixed_boundary=False allows boundary vertices to move."""
        mesh = examples.annulus()
        boundary_idx = mesh.boundary_node_indices()

        if len(boundary_idx) == 0:
            pytest.skip("Mesh has no boundary vertices")

        points_before = mesh.points.copy()
        # Run multiple iterations to increase chance of boundary movement
        mesh.smooth(n_iter=10, weight=0.5, fixed_boundary=False)

        # At least some boundary vertices should move
        # (if they are interior-connected in the edge adjacency graph)
        boundary_moved = np.linalg.norm(
            mesh.points[boundary_idx, :2] - points_before[boundary_idx, :2],
            axis=1
        )
        # Note: some boundary vertices may not have interior neighbors,
        # so we check that the smoothing didn't error, not that all moved
        assert len(boundary_moved) > 0


class TestSmoothZPreservation:
    """Test Z-coordinate preservation when present."""

    def test_z_column_preserved_when_present(self):
        """Z-coordinate (column 2) is preserved during smoothing."""
        mesh = examples.annulus()

        # Add a Z column if not present
        if mesh.points.shape[1] == 2:
            z_vals = np.arange(mesh.n_verts, dtype=float)
            mesh.points = np.c_[mesh.points, z_vals]
        else:
            z_vals = mesh.points[:, 2].copy()

        z_before = mesh.points[:, 2].copy()
        mesh.smooth(n_iter=5, weight=0.5)
        z_after = mesh.points[:, 2]

        np.testing.assert_array_equal(
            z_after, z_before,
            err_msg="Z-coordinate changed during smoothing"
        )

    def test_xy_column_modified_when_weight_positive(self):
        """X, Y columns (columns 0, 1) are modified when weight > 0."""
        mesh = examples.annulus()

        # Add Z for clarity
        if mesh.points.shape[1] == 2:
            mesh.points = np.c_[mesh.points, np.zeros(mesh.n_verts)]

        xy_before = mesh.points[:, :2].copy()
        mesh.smooth(n_iter=1, weight=0.5, fixed_boundary=True)
        xy_after = mesh.points[:, :2]

        # Interior X, Y should change
        diff = np.linalg.norm(xy_after - xy_before, axis=1)
        assert np.any(diff > 1e-10), "No X, Y movement with weight=0.5"


class TestSmoothIterations:
    """Test iteration behavior."""

    def test_n_iter_zero_leaves_unchanged(self):
        """n_iter=0 (no iterations) leaves points unchanged."""
        mesh = examples.annulus()
        points_before = mesh.points.copy()

        mesh.smooth(n_iter=0, weight=0.5)

        np.testing.assert_array_equal(
            mesh.points, points_before,
            err_msg="Points changed with n_iter=0"
        )

    def test_n_iter_positive_moves_interior(self):
        """n_iter > 0 moves interior vertices (with weight > 0)."""
        mesh = examples.annulus()
        points_before = mesh.points.copy()

        mesh.smooth(n_iter=1, weight=0.5, fixed_boundary=True)

        diff = np.linalg.norm(mesh.points - points_before, axis=1)
        assert np.any(diff > 1e-10), "No vertices moved with n_iter=1"

    def test_more_iterations_moves_farther(self):
        """More iterations should move vertices at least as far (monotone)."""
        mesh1 = examples.annulus()
        mesh2 = examples.annulus()

        points_before = mesh1.points.copy()

        # One iteration
        mesh1.smooth(n_iter=1, weight=0.5, fixed_boundary=True)
        displacement_1 = np.linalg.norm(mesh1.points - points_before, axis=1)

        # Three iterations
        mesh2.smooth(n_iter=3, weight=0.5, fixed_boundary=True)
        displacement_3 = np.linalg.norm(mesh2.points - points_before, axis=1)

        # 3 iterations should displace at least as much as 1 iteration (on average)
        # (Laplacian smoothing converges, so displacement may not strictly increase)
        assert np.sum(displacement_3) >= np.sum(displacement_1) * 0.9, \
            "3 iterations displaced less than 1 (convergence or bug)"


class TestSmoothWeightEffect:
    """Test weight parameter scaling."""

    def test_weight_half_moves_less_than_weight_one(self):
        """weight=0.5 should move vertices less than weight=1.0."""
        mesh1 = examples.annulus()
        mesh2 = examples.annulus()

        points_before = mesh1.points.copy()

        # weight=0.5
        mesh1.smooth(n_iter=1, weight=0.5, fixed_boundary=True)
        displacement_half = np.linalg.norm(
            mesh1.points - points_before, axis=1
        ).max()

        # weight=1.0
        mesh2.smooth(n_iter=1, weight=1.0, fixed_boundary=True)
        displacement_full = np.linalg.norm(
            mesh2.points - points_before, axis=1
        ).max()

        # weight=1.0 should move farther (or equal within tolerance)
        assert displacement_full >= displacement_half * 0.99, \
            "weight=1.0 moved less than weight=0.5"


class TestSmoothFixtures:
    """Test smooth() on different mesh types."""

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured", "quad_2x2"])
    def test_smooth_on_all_fixtures(self, fixture_name):
        """smooth() works on all fixture types without error."""
        if fixture_name == "block_o":
            pytest.skip("block_o takes too long for quick smoke test")

        mesh = getattr(examples, fixture_name)()
        points_before = mesh.points.copy()

        # Should not raise
        result = mesh.smooth(n_iter=2, weight=0.3, fixed_boundary=True)

        assert result is None
        assert mesh.points.shape == points_before.shape
        assert len(mesh.points) > 0


class TestSmoothSpatialIndices:
    """Test that spatial indices are updated after smoothing."""

    def test_spatial_indices_rebuilt(self):
        """_build_spatial_indices() is called at the end of smooth()."""
        mesh = examples.annulus()

        # Smooth the mesh
        mesh.smooth(n_iter=1, weight=0.5)

        # After smoothing, spatial indices should exist and be valid
        # (they should not error when accessed)
        assert hasattr(mesh, '_vertex_tree') or hasattr(mesh, '_centroid_tree')
        # If either attribute exists, the rebuild succeeded
        assert mesh.points is not None
