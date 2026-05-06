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


@pytest.fixture
def quad_mesh():
    """Fixture providing a quad-only mesh."""
    return examples.quad_2x2()


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
        """Verify element types are correctly detected as triangles only."""
        mesh = triangle_mesh
        tri_indices, quad_indices = mesh._detect_element_types()

        assert len(tri_indices) == mesh.n_elems, "All elements should be detected as triangles"
        assert len(quad_indices) == 0, "No quad elements should be detected"



class TestQuadSmoother:
    """Tests for quad-only mesh smoothing."""

    def test_fem_smoother_quad_preserves_boundary(self, quad_mesh):
        """Verify boundary nodes remain fixed after smoothing on quad mesh."""
        mesh = quad_mesh
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

    def test_fem_smoother_quad_interior_changes(self, quad_mesh):
        """Verify interior nodes are actually smoothed on quad mesh."""
        mesh = quad_mesh
        original_points = mesh.points.copy()

        # Get boundary nodes
        edge_verts = mesh.edge2vert(mesh.boundary_edges())
        boundary_nodes = set(np.unique(edge_verts.flatten()))

        # Smooth
        smoothed = mesh.direct_smoother(kinf=1e12)

        # Check some interior nodes changed (if any exist)
        interior_changed = False
        for v in range(mesh.n_verts):
            if v not in boundary_nodes:
                if not np.allclose(smoothed[v, :2], original_points[v, :2], rtol=1e-10):
                    interior_changed = True
                    break

        # Note: quad_2x2 might not have interior nodes, so this is optional
        if len(boundary_nodes) < mesh.n_verts:
            assert interior_changed, "No interior nodes were modified"

    def test_fem_smoother_quad_preserves_z(self, quad_mesh):
        """Verify z-coordinates are preserved on quad mesh."""
        mesh = quad_mesh
        original_z = mesh.points[:, 2].copy()

        smoothed = mesh.direct_smoother(kinf=1e12)

        np.testing.assert_array_equal(
            smoothed[:, 2], original_z,
            err_msg="Z-coordinates should be unchanged"
        )

    def test_fem_smoother_quad_element_type_detection(self, quad_mesh):
        """Verify element type is correctly detected as 'quad'."""
        mesh = quad_mesh
        tri_indices, quad_indices = mesh._detect_element_types()

        assert len(quad_indices) == mesh.n_elems, "All elements should be detected as quads"
        assert len(tri_indices) == 0, "No triangle elements should be detected"


class TestQuadStiffnessAssembly:
    """Unit tests for quad stiffness matrix assembly."""

    def test_quad_stiffness_assembly_method_exists(self):
        """Verify _quad_stiffness_assembly method exists and is callable."""
        from chilmesh import CHILmesh
        mesh = CHILmesh()  # Random mesh
        assert hasattr(mesh, '_quad_stiffness_assembly'), "Method should exist"
        assert callable(mesh._quad_stiffness_assembly), "Method should be callable"

    def test_quad_stiffness_assembly_returns_correct_format(self, quad_mesh):
        """Verify _quad_stiffness_assembly returns correct format."""
        mesh = quad_mesh
        p = mesh.points[:, :2]
        n = mesh.n_verts
        quad_indices = np.arange(mesh.n_elems)

        # Call the method directly
        rows, cols, data = mesh._quad_stiffness_assembly(quad_indices, p, n)

        # Verify return types are lists
        assert isinstance(rows, list), "rows should be a list"
        assert isinstance(cols, list), "cols should be a list"
        assert isinstance(data, list), "data should be a list"

        # Verify they have the same length
        assert len(rows) == len(cols) == len(data), "Lists should have same length"

        # Verify data contains floats
        assert all(isinstance(x, (int, float, np.number)) for x in data), "data should contain numbers"


class TestMixedStiffnessAssembly:
    """Unit tests for mixed stiffness matrix assembly."""

    def test_mixed_stiffness_assembly_method_exists(self):
        """Verify _mixed_stiffness_assembly method exists and is callable."""
        from chilmesh import CHILmesh
        mesh = CHILmesh()  # Random mesh
        assert hasattr(mesh, '_mixed_stiffness_assembly'), "Method should exist"
        assert callable(mesh._mixed_stiffness_assembly), "Method should be callable"

    def test_mixed_stiffness_assembly_tri_only(self, triangle_mesh):
        """Verify _mixed_stiffness_assembly works on triangle-only mesh."""
        mesh = triangle_mesh
        p = mesh.points[:, :2]
        n = mesh.n_verts
        tri_indices, quad_indices = mesh._detect_element_types()

        # Call mixed stiffness assembly on pure triangle mesh
        rows, cols, data = mesh._mixed_stiffness_assembly(tri_indices, quad_indices, p, n)

        # Verify return types are lists
        assert isinstance(rows, list), "rows should be a list"
        assert isinstance(cols, list), "cols should be a list"
        assert isinstance(data, list), "data should be a list"

        # Verify they have the same length
        assert len(rows) == len(cols) == len(data), "Lists should have same length"

    def test_mixed_stiffness_assembly_quad_only(self, quad_mesh):
        """Verify _mixed_stiffness_assembly works on quad-only mesh."""
        mesh = quad_mesh
        p = mesh.points[:, :2]
        n = mesh.n_verts
        tri_indices, quad_indices = mesh._detect_element_types()

        # Call mixed stiffness assembly on pure quad mesh
        rows, cols, data = mesh._mixed_stiffness_assembly(tri_indices, quad_indices, p, n)

        # Verify return types are lists
        assert isinstance(rows, list), "rows should be a list"
        assert isinstance(cols, list), "cols should be a list"
        assert isinstance(data, list), "data should be a list"

        # Verify they have the same length
        assert len(rows) == len(cols) == len(data), "Lists should have same length"

    def test_mixed_stiffness_assembly_consistency_triangles(self, triangle_mesh):
        """Verify mixed assembly produces consistent results for triangle-only mesh."""
        mesh = triangle_mesh
        p = mesh.points[:, :2]
        n = mesh.n_verts
        tri_indices, quad_indices = mesh._detect_element_types()

        # Get stiffness from both methods
        tri_rows, tri_cols, tri_data = mesh._tri_stiffness_assembly(tri_indices, p, n)
        mixed_rows, mixed_cols, mixed_data = mesh._mixed_stiffness_assembly(tri_indices, quad_indices, p, n)

        # Should have same number of entries
        assert len(tri_data) == len(mixed_data), "Should have same number of entries"

    def test_mixed_stiffness_assembly_consistency_quads(self, quad_mesh):
        """Verify mixed assembly produces consistent results for quad-only mesh."""
        mesh = quad_mesh
        p = mesh.points[:, :2]
        n = mesh.n_verts
        tri_indices, quad_indices = mesh._detect_element_types()

        # Get stiffness from both methods
        quad_rows, quad_cols, quad_data = mesh._quad_stiffness_assembly(quad_indices, p, n)
        mixed_rows, mixed_cols, mixed_data = mesh._mixed_stiffness_assembly(tri_indices, quad_indices, p, n)

        # Should have same number of entries
        assert len(quad_data) == len(mixed_data), "Should have same number of entries"


class TestMixedSmoother:
    """Tests for mixed-element (tri + quad) mesh smoothing."""

    def test_fem_smoother_triangle_mesh_backward_compat(self, triangle_mesh):
        """Verify direct_smoother works on triangle-only meshes (backward compat)."""
        mesh = triangle_mesh
        original_points = mesh.points.copy()

        # Smooth using direct_smoother
        smoothed = mesh.direct_smoother(kinf=1e12)

        # Verify output structure
        assert smoothed.shape == original_points.shape, "Output should have same shape"

        # Verify boundary preservation
        edge_verts = mesh.edge2vert(mesh.boundary_edges())
        boundary_nodes = np.unique(edge_verts.flatten())

        for v in boundary_nodes:
            np.testing.assert_allclose(
                smoothed[v], original_points[v],
                rtol=1e-8, atol=1e-11
            )

    def test_fem_smoother_quad_mesh_support(self, quad_mesh):
        """Verify direct_smoother works on quad-only meshes."""
        mesh = quad_mesh
        original_points = mesh.points.copy()

        # Smooth using direct_smoother
        smoothed = mesh.direct_smoother(kinf=1e12)

        # Verify output structure
        assert smoothed.shape == original_points.shape, "Output should have same shape"

        # Verify boundary preservation
        edge_verts = mesh.edge2vert(mesh.boundary_edges())
        boundary_nodes = np.unique(edge_verts.flatten())

        for v in boundary_nodes:
            np.testing.assert_allclose(
                smoothed[v], original_points[v],
                rtol=1e-8, atol=1e-11
            )


class TestSmootherIntegration:
    """Integration tests for smooth_mesh() method."""

    def test_smooth_mesh_fem_triangles(self, triangle_mesh):
        """Test smooth_mesh() method with FEM for triangles."""
        mesh = triangle_mesh
        mesh_copy = mesh.copy()

        # Should not raise
        mesh_copy.smooth_mesh(method="fem", acknowledge_change=True)

