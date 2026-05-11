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
        """FEM smoother moves a perturbed interior node back toward equilibrium.

        A uniform grid is already at the FEM optimum — no movement is expected
        or required. Perturb the interior node first so the smoother has real
        work to do, then verify it returns closer to the symmetric center.
        """
        mesh = quad_mesh
        edge_verts = mesh.edge2vert(mesh.boundary_edges())
        boundary_nodes = set(np.unique(edge_verts.flatten()))
        interior_nodes = [v for v in range(mesh.n_verts) if v not in boundary_nodes]

        if not interior_nodes:
            return  # no interior nodes — nothing to test

        # Perturb the interior node away from center
        pts = mesh.points.copy()
        v = interior_nodes[0]
        original_xy = pts[v, :2].copy()
        pts[v, :2] += 0.3  # push off center
        mesh.change_points(pts, acknowledge_change=True)

        smoothed = mesh.direct_smoother(kinf=1e12)

        # After smoothing, interior node should be closer to original than the perturbed pos
        dist_after = np.linalg.norm(smoothed[v, :2] - original_xy)
        dist_before = np.linalg.norm(pts[v, :2] - original_xy)
        assert dist_after < dist_before, (
            f"Interior node {v} did not move toward equilibrium after smoothing; "
            f"dist_before={dist_before:.4f} dist_after={dist_after:.4f}"
        )

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


class TestMixedElementFEMSmoother:
    """Regression tests for FEM smoother on mixed tri+quad meshes (issues #95, #100).

    The symmetric quad decomposition (both diagonals, each weighted 0.5) must give
    equal stiffness to all four quad vertices and must not collapse interior nodes.
    """

    @pytest.fixture
    def mixed_mesh(self):
        """Mixed-element mesh: 8 padded triangles (top) + 4 quads (bottom), interior nodes at row y=1."""
        from chilmesh import CHILmesh

        # 3x3 point grid:  top row y=2, seam row y=1 (interior nodes 6,7,8), bottom row y=0
        pts = np.array(
            [[i * 1.0, j, 0.0] for j in [2.0, 1.0, 0.0] for i in range(5)],
            dtype=float,
        )
        # 8 CCW triangles covering top two rows, padded [v0,v1,v2,v0]
        def pad(a, b, c):
            return [a, b, c, a]

        top = np.array(
            [pad(0, 6, 5), pad(0, 1, 6), pad(1, 7, 6), pad(1, 2, 7),
             pad(2, 8, 7), pad(2, 3, 8), pad(3, 9, 8), pad(3, 4, 9)],
            dtype=int,
        )
        # 4 CCW quads covering bottom two rows
        bot = np.array(
            [[5, 6, 11, 10], [6, 7, 12, 11], [7, 8, 13, 12], [8, 9, 14, 13]],
            dtype=int,
        )
        conn = np.vstack([top, bot])

        mesh = CHILmesh()
        mesh.points = pts
        mesh.connectivity_list = conn
        mesh.n_verts = len(pts)
        mesh.n_elems = len(conn)
        mesh._build_adjacencies()
        return mesh

    def test_quad_stiffness_is_symmetric(self):
        """All four quad vertices must receive equal diagonal stiffness (1.5xD each)."""
        from chilmesh import examples
        from scipy.sparse import csr_matrix

        mesh = examples.quad_2x2()
        p = mesh.points[:, :2]
        n = mesh.n_verts
        _, quad_idx = mesh._detect_element_types()

        # Restrict to the first quad only to get element-local values
        rows, cols, data = mesh._quad_stiffness_assembly(quad_idx[:1], p, n)
        K = csr_matrix((data, (rows, cols)), shape=(2 * n, 2 * n)).toarray()

        first_quad = mesh.connectivity_list[quad_idx[0]]
        diag_vals = [K[2 * v, 2 * v] for v in first_quad]
        # All four vertices must have equal diagonal (1.5 x 2 = 3.0)
        assert all(abs(d - 3.0) < 1e-10 for d in diag_vals), (
            f"Quad diagonal stiffness is asymmetric: {diag_vals}"
        )

    def test_mixed_mesh_has_interior_nodes(self, mixed_mesh):
        """Fixture sanity: vertices 6, 7, 8 should be interior."""
        edge_verts = mixed_mesh.edge2vert(mixed_mesh.boundary_edges())
        boundary = set(np.unique(edge_verts.flatten()))
        interior = [v for v in range(mixed_mesh.n_verts) if v not in boundary]
        assert interior == [6, 7, 8], f"Expected interior nodes [6,7,8], got {interior}"

    def test_fem_smoother_mixed_no_nan(self, mixed_mesh):
        """FEM smoother must not produce NaN on a mixed-element mesh."""
        smoothed = mixed_mesh.direct_smoother(kinf=1e12)
        assert not np.any(np.isnan(smoothed)), "FEM smoother produced NaN on mixed mesh"

    def test_fem_smoother_mixed_no_element_collapse(self, mixed_mesh):
        """FEM smoother must not produce folded (zero-area) elements."""
        # Perturb interior nodes, smooth, check no element collapses
        pts = mixed_mesh.points.copy()
        for v in [6, 7, 8]:
            pts[v, :2] += 0.2
        mixed_mesh.change_points(pts, acknowledge_change=True)

        smoothed = mixed_mesh.direct_smoother(kinf=1e12)
        new_pts = mixed_mesh.points.copy()
        new_pts[:, :2] = smoothed[:, :2]
        mixed_mesh.change_points(new_pts, acknowledge_change=True)

        q, _, _ = mixed_mesh.elem_quality()
        assert q.min() > 0.0, f"FEM smoother folded element: min quality={q.min():.4f}"

    def test_fem_smoother_mixed_boundary_preserved(self, mixed_mesh):
        """Boundary nodes must not move after FEM smoothing on mixed mesh."""
        original_points = mixed_mesh.points.copy()
        edge_verts = mixed_mesh.edge2vert(mixed_mesh.boundary_edges())
        boundary_nodes = np.unique(edge_verts.flatten())

        smoothed = mixed_mesh.direct_smoother(kinf=1e12)

        for v in boundary_nodes:
            np.testing.assert_allclose(
                smoothed[v], original_points[v], rtol=1e-8, atol=1e-9,
                err_msg=f"Boundary node {v} moved after mixed-mesh FEM smooth"
            )

    def test_fem_smoother_mixed_interior_moves_to_equilibrium(self, mixed_mesh):
        """Perturbed interior nodes must move back toward their grid positions."""
        orig_pts = mixed_mesh.points.copy()

        # Perturb all 3 interior nodes
        pts = orig_pts.copy()
        for v in [6, 7, 8]:
            pts[v, :2] += 0.25
        mixed_mesh.change_points(pts, acknowledge_change=True)

        smoothed = mixed_mesh.direct_smoother(kinf=1e12)

        for v in [6, 7, 8]:
            dist_after = np.linalg.norm(smoothed[v, :2] - orig_pts[v, :2])
            dist_before = np.linalg.norm(pts[v, :2] - orig_pts[v, :2])
            assert dist_after < dist_before, (
                f"Interior node {v} did not move toward equilibrium: "
                f"dist_before={dist_before:.4f}, dist_after={dist_after:.4f}"
            )


class TestSmootherIntegration:
    """Integration tests for smooth_mesh() method."""

    def test_smooth_mesh_fem_triangles(self, triangle_mesh):
        """Test smooth_mesh() method with FEM for triangles."""
        mesh = triangle_mesh
        mesh_copy = mesh.copy()

        # Should not raise
        mesh_copy.smooth_mesh(method="fem", acknowledge_change=True)
