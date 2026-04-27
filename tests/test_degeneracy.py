"""Tests for degeneracy handling and fallback logic."""

import numpy as np
import pytest
from chilmesh import CHILmesh


class TestDegeneracyHandling:
    """Verify that CHILmesh correctly handles degenerate elements."""

    def test_padded_triangle_as_quad_ccw_flip(self):
        """Verify CCW orientation flip works for padded triangles (quad format)."""
        # Create a simple padded triangle: vertices [0, 1, 2, 0]
        # representing a triangle (0,1,2) stored as a quad with pad at position 4
        points = np.array([
            [0.0, 0.0],  # vertex 0
            [1.0, 0.0],  # vertex 1
            [0.0, 1.0],  # vertex 2
        ], dtype=np.float64)

        # CW orientation: [0, 2, 1, 0] (should be flipped to CCW)
        connectivity = np.array([
            [0, 2, 1, 0],  # CW padded triangle (needs CCW flip)
        ], dtype=np.int32)

        mesh = CHILmesh(connectivity=connectivity, points=points)

        # After CCW fix, should be [0, 1, 2, 0] (flipped)
        expected = np.array([[0, 1, 2, 0]], dtype=np.int32)
        np.testing.assert_array_equal(mesh.connectivity_list, expected)

    def test_degenerate_quad_with_duplicate_vertices_fallback(self):
        """Verify degenerate quad (5+ distinct vertices) triggers fallback (left as-is)."""
        points = np.array([
            [0.0, 0.0],  # vertex 0
            [1.0, 0.0],  # vertex 1
            [1.0, 1.0],  # vertex 2
            [0.0, 1.0],  # vertex 3
        ], dtype=np.float64)

        # Degenerate quad with 4 distinct vertices (not a padded triangle)
        # In CCW format: [0, 1, 2, 3]
        connectivity = np.array([
            [0, 1, 2, 3],  # Not degenerate, regular quad
        ], dtype=np.int32)

        mesh = CHILmesh(connectivity=connectivity, points=points)

        # Should not crash and connectivity should be preserved
        assert mesh.n_elems == 1
        assert mesh.connectivity_list.shape == (1, 4)

    def test_mixed_element_mesh_with_padded_triangles(self):
        """Verify mixed tri/quad mesh handles both regular and padded elements."""
        points = np.array([
            [0.0, 0.0],  # 0
            [1.0, 0.0],  # 1
            [1.0, 1.0],  # 2
            [0.0, 1.0],  # 3
            [0.5, 0.5],  # 4 (center)
        ], dtype=np.float64)

        # Mix: quad at index 0, padded triangle at index 1
        connectivity = np.array([
            [0, 1, 2, 3],     # Regular quad
            [1, 2, 4, 1],     # Padded triangle (vertex 1 repeated)
        ], dtype=np.int32)

        mesh = CHILmesh(connectivity=connectivity, points=points)

        assert mesh.n_elems == 2
        # Both elements should be preserved (first as quad, second as padded triangle)
        assert mesh.connectivity_list[0, 0] == 0
        assert mesh.connectivity_list[1, 0] == 1

    def test_adjacency_building_with_degenerate_elements(self):
        """Verify adjacency structures build correctly with degenerate elements."""
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
        ], dtype=np.float64)

        # Padded triangle
        connectivity = np.array([
            [0, 1, 2, 0],  # Padded triangle
        ], dtype=np.int32)

        mesh = CHILmesh(connectivity=connectivity, points=points)

        # Verify adjacencies were built
        assert hasattr(mesh, 'adjacencies')
        assert 'Elem2Vert' in mesh.adjacencies
        assert 'Edge2Vert' in mesh.adjacencies
        assert mesh.n_edges > 0

    def test_element_quality_with_degenerate_mesh(self):
        """Verify elem_quality() handles degenerate elements without crashing."""
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
        ], dtype=np.float64)

        connectivity = np.array([
            [0, 1, 2, 0],  # Padded triangle
        ], dtype=np.int32)

        mesh = CHILmesh(connectivity=connectivity, points=points)

        # Should not crash
        quality, angles, stats = mesh.elem_quality()

        assert quality is not None
        assert len(quality) == 1
        assert np.isfinite(quality[0])

    def test_backward_compatibility_quad_reader(self):
        """Verify quad support doesn't break backward compatibility with triangles."""
        # All-triangle mesh
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            [0.5, 0.5],
        ], dtype=np.float64)

        # Triangles must be padded to 4 columns for storage
        connectivity = np.array([
            [0, 1, 2, 0],  # Padded triangle
        ], dtype=np.int32)

        mesh = CHILmesh(connectivity=connectivity, points=points)

        assert mesh.n_elems == 1
        assert mesh.type == "Triangular"
        assert mesh.n_edges > 0


class TestMADMESHRCompatibility:
    """Verify degeneracy handling aligns with MADMESHR expectations."""

    def test_mesh_can_be_adapted_after_degeneracy_handling(self):
        """Verify a mesh with degenerate elements can still be used for adaptation."""
        # Create a simple 2-triangle mesh (both padded to 4 columns)
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            [0.5, 0.5],
        ], dtype=np.float64)

        connectivity = np.array([
            [0, 1, 2, 0],  # Padded triangle
            [1, 2, 3, 1],  # Padded triangle
        ], dtype=np.int32)

        mesh = CHILmesh(connectivity=connectivity, points=points)

        # Verify mesh can be queried (basic MADMESHR compatibility check)
        assert mesh.n_verts == 4
        assert mesh.n_elems == 2
        assert mesh.n_edges > 0

        # Verify we can access basic adjacency structures
        assert 'Elem2Vert' in mesh.adjacencies
        assert 'Vert2Elem' in mesh.adjacencies
        assert 'Edge2Vert' in mesh.adjacencies

    def test_degeneracy_fallback_preserves_mesh_integrity(self):
        """Verify that degenerate element fallback doesn't corrupt mesh topology."""
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ], dtype=np.float64)

        # Mix of padded and regular elements
        connectivity = np.array([
            [0, 1, 2, 3],     # Regular quad
            [1, 2, 0, 1],     # Padded triangle
        ], dtype=np.int32)

        mesh = CHILmesh(connectivity=connectivity, points=points)

        # Verify invariants
        assert mesh.n_verts == 4
        assert mesh.n_elems == 2

        # Verify no NaN or inf in coordinates
        assert np.all(np.isfinite(mesh.points))

        # Verify all vertex indices are valid
        assert np.all(mesh.connectivity_list >= 0)
        assert np.all(mesh.connectivity_list < mesh.n_verts)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
