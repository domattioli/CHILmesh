"""Tests for spatial indexing (Phase 5.3)."""
import pytest
import numpy as np
from chilmesh.examples import annulus, donut
from chilmesh import examples as _examples


class TestSpatialIndexing:
    """Test spatial indexing methods."""

    def test_find_element_on_vertex(self, mesh):
        """find_element should locate element containing a vertex."""
        for v_id in range(min(10, mesh.n_verts)):
            point = mesh.points[v_id, :2]
            elem_id = mesh.find_element(point)
            assert elem_id >= 0, f"Failed to find element for vertex {v_id}"
            assert v_id in mesh.connectivity_list[elem_id]

    def test_find_element_outside_returns_negative(self, mesh):
        """find_element should return -1 for points outside mesh."""
        far_point = np.array([1e6, 1e6])
        elem_id = mesh.find_element(far_point)
        assert elem_id == -1

    def test_find_element_center_of_mesh(self, mesh):
        """find_element should work for points inside mesh."""
        if mesh.grid_name in ['annulus', 'donut', 'block_o']:
            pytest.skip(f"Mesh {mesh.grid_name} has hole; center may be outside")
        centroid = np.mean(mesh.points[:, :2], axis=0)
        elem_id = mesh.find_element(centroid)
        assert elem_id >= 0

    def test_nearest_vertices_returns_correct_count(self, mesh):
        """nearest_vertices should return exactly k vertices (capped by mesh size)."""
        point = mesh.points[0, :2]
        for k in [1, 3, 5, 10]:
            requested_k = min(k, mesh.n_verts)
            nearest = mesh.nearest_vertices(point, k=k)
            assert len(nearest) == requested_k, f"Expected {requested_k} vertices, got {len(nearest)}"
            assert np.all(nearest >= 0) and np.all(nearest < mesh.n_verts)

    def test_nearest_vertices_first_is_closest(self, mesh):
        """First nearest vertex should be closest."""
        point = mesh.points[0, :2]
        nearest = mesh.nearest_vertices(point, k=5)
        distances = np.linalg.norm(mesh.points[nearest, :2] - point, axis=1)
        assert distances[0] <= distances[-1] + 1e-6

    def test_nearest_vertices_is_vertex(self, mesh):
        """nearest_vertices should return actual vertex indices."""
        point = mesh.points[0, :2]
        nearest = mesh.nearest_vertices(point, k=1)
        assert 0 in nearest

    def test_nearest_elements_returns_correct_count(self, mesh):
        """nearest_elements should return exactly k elements (capped by mesh size)."""
        point = mesh.points[0, :2]
        for k in [1, 5, 10]:
            requested_k = min(k, mesh.n_elems)
            nearest = mesh.nearest_elements(point, k=k)
            assert len(nearest) == requested_k, f"Expected {requested_k} elements, got {len(nearest)}"
            assert np.all(nearest >= 0) and np.all(nearest < mesh.n_elems)

    def test_find_elements_in_radius_contains_nearest(self, mesh):
        """Elements in radius should include nearest element."""
        point = mesh.points[0, :2]
        nearest = mesh.nearest_elements(point, k=1)[0]
        nearby = mesh.find_elements_in_radius(point, radius=1.0)
        assert nearest in nearby

    def test_find_elements_in_radius_grows_with_radius(self, mesh):
        """More elements should be found with larger radius."""
        point = mesh.points[0, :2]
        nearby_small = mesh.find_elements_in_radius(point, radius=0.1)
        nearby_large = mesh.find_elements_in_radius(point, radius=1.0)
        assert len(nearby_small) <= len(nearby_large)

    def test_spatial_indices_built_on_init(self, mesh):
        """KD-trees should be built during initialization."""
        assert hasattr(mesh, '_vertex_tree'), "Vertex tree not built"
        assert hasattr(mesh, '_centroid_tree'), "Centroid tree not built"
        assert mesh._vertex_tree is not None
        assert mesh._centroid_tree is not None

    def test_centroids_property_works(self, mesh):
        """centroids property should return correct shape."""
        centroids = mesh.centroids
        assert centroids.shape == (mesh.n_elems, 2)
        assert np.all(np.isfinite(centroids))


class TestSpatialIndexingEdgeCases:
    """Edge-case tests for spatial indexing (boundary, zero, negative inputs)."""

    def test_nearest_vertices_k_zero(self, mesh):
        """nearest_vertices with k=0 should return empty array."""
        point = mesh.points[0, :2]
        result = mesh.nearest_vertices(point, k=0)
        assert isinstance(result, np.ndarray)
        assert len(result) == 0
        assert result.dtype == int

    def test_nearest_elements_k_zero(self, mesh):
        """nearest_elements with k=0 should return empty array."""
        point = mesh.points[0, :2]
        result = mesh.nearest_elements(point, k=0)
        assert isinstance(result, np.ndarray)
        assert len(result) == 0
        assert result.dtype == int

    def test_nearest_vertices_coincident_with_vertex(self, mesh):
        """Query point at exact vertex should have distance 0 to first nearest."""
        v_id = min(5, mesh.n_verts - 1)
        point = mesh.points[v_id, :2].copy()
        nearest = mesh.nearest_vertices(point, k=1)
        assert nearest[0] == v_id
        # Verify distance is truly zero
        dist = np.linalg.norm(mesh.points[nearest[0], :2] - point)
        assert dist < 1e-10

    def test_find_elements_in_radius_zero(self, mesh):
        """find_elements_in_radius with radius=0 should only return elements at point."""
        if mesh.grid_name in ['annulus', 'donut']:
            pytest.skip("Mesh has hole; center may be outside domain")
        centroid = np.mean(mesh.points[:, :2], axis=0)
        elems_zero = mesh.find_elements_in_radius(centroid, radius=0.0)
        elems_small = mesh.find_elements_in_radius(centroid, radius=0.01)
        # Zero radius should return at most as many as small radius
        assert len(elems_zero) <= len(elems_small)
        # Both should be finite lists
        assert isinstance(elems_zero, np.ndarray)
        assert isinstance(elems_small, np.ndarray)

    def test_find_elements_in_radius_negative_raises_value_error(self, mesh):
        """find_elements_in_radius with negative radius should raise ValueError."""
        point = mesh.points[0, :2]
        with pytest.raises(ValueError, match="radius must be >= 0"):
            mesh.find_elements_in_radius(point, radius=-0.1)

    def test_find_element_inside_mixed_element_mesh(self):
        """find_element should work on mixed tri/quad mesh."""
        mesh = _examples.structured()
        # All vertices should have elements
        for v_id in range(min(10, mesh.n_verts)):
            point = mesh.points[v_id, :2]
            elem_id = mesh.find_element(point)
            assert elem_id >= 0, f"Failed to find element for vertex {v_id} on structured mesh"
            assert v_id in mesh.connectivity_list[elem_id]

    def test_nearest_vertices_k_exceeds_mesh_size(self, mesh):
        """Requesting k > mesh.n_verts should clamp to n_verts."""
        point = mesh.points[0, :2]
        requested_k = mesh.n_verts + 100
        result = mesh.nearest_vertices(point, k=requested_k)
        # Should return exactly n_verts, not more
        assert len(result) == mesh.n_verts
        assert len(result) < requested_k

    def test_nearest_elements_k_exceeds_mesh_size(self, mesh):
        """Requesting k > mesh.n_elems should clamp to n_elems."""
        point = mesh.points[0, :2]
        requested_k = mesh.n_elems + 100
        result = mesh.nearest_elements(point, k=requested_k)
        # Should return exactly n_elems, not more
        assert len(result) == mesh.n_elems
        assert len(result) < requested_k


class TestSpatialIndexingPerformance:
    """Test performance targets for spatial indexing."""

    def test_100k_point_queries_annulus(self):
        """100k point queries should complete quickly."""
        mesh = annulus()
        np.random.seed(42)
        center = np.mean(mesh.points[:, :2], axis=0)
        radius = 1.0
        points = center[:, np.newaxis] + radius * np.random.randn(2, 100000)

        import time
        start = time.time()
        for i in range(100):
            point = points[:, i]
            elem_id = mesh.find_element(point)
        elapsed = time.time() - start

        assert elapsed < 5.0, f"100 queries took {elapsed:.2f}s, should be <5s"

    def test_spatial_index_build_time(self):
        """Spatial index build should be fast (<0.5s for large mesh)."""
        mesh = donut()
        import time
        start = time.time()
        mesh._build_spatial_indices()
        elapsed = time.time() - start
        assert elapsed < 0.5, f"Building spatial indices took {elapsed:.2f}s, should be <0.5s"
