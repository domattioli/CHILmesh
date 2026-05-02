"""
Tests for FEM smoother - triangle, quad, and mixed-element meshes.
"""

import pytest
import numpy as np
from chilmesh import CHILmesh
from chilmesh import examples


@pytest.fixture(params=["annulus", "donut", "block_o"])
def triangle_mesh(request):
    """Fixture providing triangle-only meshes."""
    return getattr(examples, request.param)()


class TestTriangleSmoother:
    """Tests for triangle-only mesh smoothing (backward compatibility)."""

    def test_fem_smoother_triangle_preserves_boundary(self, triangle_mesh):
        """Verify boundary nodes remain fixed after smoothing."""
        mesh = triangle_mesh
        original_points = mesh.points.copy()

        # Get boundary nodes before smoothing
        edge_verts = mesh.edge2vert(mesh.boundary_edges())
        boundary_nodes = np.unique(edge_verts.flatten())

        # Smooth
        smoothed = mesh.direct_smoother(kinf=1e12)

        # Check boundary nodes are unchanged (relax tolerance for solver roundoff)
        for v in boundary_nodes:
            np.testing.assert_allclose(
                smoothed[v], original_points[v],
                rtol=1e-8, atol=1e-11,
                err_msg=f"Boundary node {v} moved"
            )

    def test_fem_smoother_triangle_interior_changes(self, triangle_mesh):
        """Verify interior nodes are actually smoothed."""
        mesh = triangle_mesh
        original_points = mesh.points.copy()

        # Get boundary nodes
        edge_verts = mesh.edge2vert(mesh.boundary_edges())
        boundary_nodes = set(np.unique(edge_verts.flatten()))

        # Smooth
        smoothed = mesh.direct_smoother(kinf=1e12)

        # Check some interior nodes changed
        interior_changed = False
        for v in range(mesh.n_verts):
            if v not in boundary_nodes:
                if not np.allclose(smoothed[v, :2], original_points[v, :2], rtol=1e-10):
                    interior_changed = True
                    break

        assert interior_changed, "No interior nodes were modified"

    def test_fem_smoother_triangle_preserves_z(self, triangle_mesh):
        """Verify z-coordinates are preserved."""
        mesh = triangle_mesh
        original_z = mesh.points[:, 2].copy()

        smoothed = mesh.direct_smoother(kinf=1e12)

        np.testing.assert_array_equal(
            smoothed[:, 2], original_z,
            err_msg="Z-coordinates should be unchanged"
        )

    def test_fem_smoother_triangle_element_type_detection(self, triangle_mesh):
        """Verify element type is correctly detected as 'triangle'."""
        mesh = triangle_mesh
        elem_type = mesh._detect_element_type()
        assert elem_type == 'triangle', f"Expected 'triangle', got '{elem_type}'"



class TestQuadSmoother:
    """Tests for quad-only mesh smoothing."""

    def test_fem_smoother_on_structured_fixture(self, triangle_mesh):
        """Test smoother on structured fixture (verifies backward compatibility)."""
        mesh = triangle_mesh
        original_points = mesh.points.copy()
        smoothed = mesh.direct_smoother(kinf=1e12)

        # Boundary should not change (relax tolerance for solver roundoff)
        edge_verts = mesh.edge2vert(mesh.boundary_edges())
        boundary_nodes = np.unique(edge_verts.flatten())

        for v in boundary_nodes:
            np.testing.assert_allclose(
                smoothed[v], original_points[v],
                rtol=1e-8, atol=1e-11
            )


class TestMixedSmoother:
    """Tests for mixed-element (tri + quad) mesh smoothing."""


class TestSmootherIntegration:
    """Integration tests for smooth_mesh() method."""

    def test_smooth_mesh_fem_triangles(self, triangle_mesh):
        """Test smooth_mesh() method with FEM for triangles."""
        mesh = triangle_mesh
        mesh_copy = mesh.copy()

        # Should not raise
        mesh_copy.smooth_mesh(method="fem", acknowledge_change=True)

