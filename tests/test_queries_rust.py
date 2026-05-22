"""
Equivalence tests for Rust quality analysis and mesh queries.

Tests that Rust-based signed area computation, orientation enforcement,
and vertex/element queries match the Python reference implementation.
"""

import pytest
import numpy as np
from pathlib import Path
from chilmesh_core import RustMesh
from chilmesh import examples


# Test fixtures
@pytest.fixture(params=["annulus", "donut", "block_o", "structured"])
def test_fixture(request):
    """Load a test fixture via the Python examples module."""
    fixture_name = request.param
    if fixture_name == "annulus":
        mesh = examples.annulus()
    elif fixture_name == "donut":
        mesh = examples.donut()
    elif fixture_name == "block_o":
        mesh = examples.block_o()
    elif fixture_name == "structured":
        mesh = examples.structured()
    else:
        raise ValueError(f"Unknown fixture: {fixture_name}")
    return mesh, fixture_name


@pytest.fixture(params=["annulus", "donut", "block_o", "structured"])
def rust_mesh_fixture(request):
    """Load a test fixture as a Rust mesh."""
    fixture_name = request.param
    if fixture_name == "annulus":
        fort14_path = str(Path(__file__).parent.parent / "src" / "chilmesh" / "data" / "annulus_200pts.fort.14")
    elif fixture_name == "donut":
        fort14_path = str(Path(__file__).parent.parent / "src" / "chilmesh" / "data" / "donut_domain.fort.14")
    elif fixture_name == "block_o":
        fort14_path = str(Path(__file__).parent.parent / "src" / "chilmesh" / "data" / "Block_O.14")
    elif fixture_name == "structured":
        fort14_path = str(Path(__file__).parent.parent / "src" / "chilmesh" / "data" / "structuredMesh1.14")
    else:
        raise ValueError(f"Unknown fixture: {fixture_name}")

    rust_mesh = RustMesh()
    rust_mesh.read_from_fort14(fort14_path)
    return rust_mesh, fixture_name


class TestSignedAreaComputation:
    """Test signed area computation on real fixtures."""

    def test_compute_quality_produces_areas(self, rust_mesh_fixture):
        """compute_quality should populate signed areas."""
        rust_mesh, fixture_name = rust_mesh_fixture

        # Should not have areas initially
        with pytest.raises(Exception):  # RuntimeError when not computed
            rust_mesh.get_signed_areas()

        # After compute_quality, should have areas
        rust_mesh.compute_quality()
        areas = rust_mesh.get_signed_areas()

        assert len(areas) == rust_mesh.n_elems, \
            f"Expected {rust_mesh.n_elems} areas, got {len(areas)} in {fixture_name}"

        # All areas should be finite
        assert np.all(np.isfinite(areas)), \
            f"Found non-finite areas in {fixture_name}"

    def test_signed_area_matches_python(self, test_fixture):
        """Signed areas from Rust should match Python reference within ±0.1%."""
        py_mesh, fixture_name = test_fixture

        # Load same mesh in Rust
        if fixture_name == "annulus":
            fort14_path = str(Path(__file__).parent.parent / "src" / "chilmesh" / "data" / "annulus_200pts.fort.14")
        elif fixture_name == "donut":
            fort14_path = str(Path(__file__).parent.parent / "src" / "chilmesh" / "data" / "donut_domain.fort.14")
        elif fixture_name == "block_o":
            fort14_path = str(Path(__file__).parent.parent / "src" / "chilmesh" / "data" / "Block_O.14")
        elif fixture_name == "structured":
            fort14_path = str(Path(__file__).parent.parent / "src" / "chilmesh" / "data" / "structuredMesh1.14")

        rust_mesh = RustMesh()
        rust_mesh.read_from_fort14(fort14_path)
        rust_mesh.compute_quality()
        rust_areas = rust_mesh.get_signed_areas()

        # Get Python reference areas
        py_areas = py_mesh.signed_area()

        # Check shapes match
        assert rust_areas.shape[0] == py_areas.shape[0], \
            f"Shape mismatch: Rust {rust_areas.shape} vs Python {py_areas.shape}"

        # Check all areas match within tolerance (±0.1%)
        # Use absolute tolerance for small areas
        relative_errors = np.abs((rust_areas - py_areas) / (np.abs(py_areas) + 1e-10))
        max_error = np.max(relative_errors)

        assert max_error < 0.001, \
            f"Max relative error {max_error:.6f} exceeds 0.1% tolerance for {fixture_name}"

        # Check a few specific elements match exactly (within floating point tolerance)
        sample_indices = [0, 1, min(10, len(py_areas) - 1), len(py_areas) - 1]
        for i in sample_indices:
            if i < len(py_areas):
                assert np.isclose(rust_areas[i], py_areas[i], rtol=1e-6, atol=1e-12), \
                    f"Element {i} mismatch: Rust {rust_areas[i]:.10f} vs Python {py_areas[i]:.10f}"


class TestCCWOrientation:
    """Test counter-clockwise orientation enforcement."""

    def test_ensure_ccw_all_positive_after(self, rust_mesh_fixture):
        """After ensure_ccw, all elements should have positive area."""
        rust_mesh, fixture_name = rust_mesh_fixture

        rust_mesh.ensure_ccw()
        areas = rust_mesh.get_signed_areas()

        # All areas should be positive
        assert np.all(areas > 0.0), \
            f"Found {np.sum(areas <= 0)} elements with non-positive area in {fixture_name}"

        # Check min and mean areas (sanity check)
        min_area = np.min(areas)
        mean_area = np.mean(areas)
        assert min_area > 0.0, f"Min area {min_area} is not positive"
        assert mean_area > 0.0, f"Mean area {mean_area} is not positive"

    def test_ensure_ccw_preserves_area_magnitude(self, rust_mesh_fixture):
        """ensure_ccw should only flip sign, not change magnitude."""
        rust_mesh, fixture_name = rust_mesh_fixture

        # Compute quality before CCW
        rust_mesh.compute_quality()
        areas_before = np.abs(rust_mesh.get_signed_areas()).copy()

        # Ensure CCW
        rust_mesh.ensure_ccw()
        areas_after = np.abs(rust_mesh.get_signed_areas())

        # Magnitudes should match
        assert np.allclose(areas_before, areas_after, rtol=1e-10), \
            f"Area magnitudes changed after CCW in {fixture_name}"


class TestVertexEdgeQueries:
    """Test get_vertex_edges query method."""

    def test_get_vertex_edges_returns_list(self, rust_mesh_fixture):
        """get_vertex_edges should return a list of edge IDs."""
        rust_mesh, fixture_name = rust_mesh_fixture
        rust_mesh.build_adjacencies()

        edges_v0 = rust_mesh.get_vertex_edges(0)
        assert isinstance(edges_v0, list), f"Should return list, got {type(edges_v0)}"
        assert all(isinstance(e, int) for e in edges_v0), "All edge IDs should be integers"

    def test_get_vertex_edges_deduplicates(self, rust_mesh_fixture):
        """get_vertex_edges should not return duplicate edge IDs."""
        rust_mesh, fixture_name = rust_mesh_fixture
        rust_mesh.build_adjacencies()

        # Query vertices (limit to first 50 for speed)
        for v in range(min(50, rust_mesh.n_verts)):
            edges = rust_mesh.get_vertex_edges(v)
            unique_edges = set(edges)
            assert len(edges) == len(unique_edges), \
                f"Vertex {v} in {fixture_name} has duplicate edges"

    def test_get_vertex_edges_sorted(self, rust_mesh_fixture):
        """get_vertex_edges should return edges in sorted order."""
        rust_mesh, fixture_name = rust_mesh_fixture
        rust_mesh.build_adjacencies()

        for v in range(min(50, rust_mesh.n_verts)):
            edges = rust_mesh.get_vertex_edges(v)
            if len(edges) > 1:
                assert edges == sorted(edges), \
                    f"Vertex {v} edges not sorted: {edges} in {fixture_name}"

    def test_get_vertex_edges_valid_range(self, rust_mesh_fixture):
        """All returned edge IDs should be in valid range."""
        rust_mesh, fixture_name = rust_mesh_fixture
        rust_mesh.build_adjacencies()

        edge2vert = rust_mesh.get_edge2vert()
        n_edges = edge2vert.shape[0]

        for v in range(min(50, rust_mesh.n_verts)):
            edges = rust_mesh.get_vertex_edges(v)
            for edge_id in edges:
                assert 0 <= edge_id < n_edges, \
                    f"Edge ID {edge_id} out of range [0, {n_edges}) for vertex {v} in {fixture_name}"


class TestVertexElementQueries:
    """Test get_vertex_elements query method."""

    def test_get_vertex_elements_returns_list(self, rust_mesh_fixture):
        """get_vertex_elements should return a list of element IDs."""
        rust_mesh, fixture_name = rust_mesh_fixture

        elems_v0 = rust_mesh.get_vertex_elements(0)
        assert isinstance(elems_v0, list), f"Should return list, got {type(elems_v0)}"
        assert all(isinstance(e, int) for e in elems_v0), "All element IDs should be integers"

    def test_get_vertex_elements_deduplicates(self, rust_mesh_fixture):
        """get_vertex_elements should not return duplicate element IDs."""
        rust_mesh, fixture_name = rust_mesh_fixture

        for v in range(min(50, rust_mesh.n_verts)):
            elems = rust_mesh.get_vertex_elements(v)
            unique_elems = set(elems)
            assert len(elems) == len(unique_elems), \
                f"Vertex {v} in {fixture_name} has duplicate elements"

    def test_get_vertex_elements_sorted(self, rust_mesh_fixture):
        """get_vertex_elements should return elements in sorted order."""
        rust_mesh, fixture_name = rust_mesh_fixture

        for v in range(min(50, rust_mesh.n_verts)):
            elems = rust_mesh.get_vertex_elements(v)
            if len(elems) > 1:
                assert elems == sorted(elems), \
                    f"Vertex {v} elements not sorted: {elems} in {fixture_name}"

    def test_get_vertex_elements_valid_ids(self, rust_mesh_fixture):
        """All returned element IDs should be in valid range."""
        rust_mesh, fixture_name = rust_mesh_fixture

        for v in range(min(50, rust_mesh.n_verts)):
            elems = rust_mesh.get_vertex_elements(v)
            for elem_id in elems:
                assert 0 <= elem_id < rust_mesh.n_elems, \
                    f"Element ID {elem_id} out of range [0, {rust_mesh.n_elems}) for vertex {v} in {fixture_name}"


class TestElementVertexQueries:
    """Test get_element_vertices query method."""

    def test_get_element_vertices_returns_array(self, rust_mesh_fixture):
        """get_element_vertices should return a 4-element array."""
        rust_mesh, fixture_name = rust_mesh_fixture

        verts = rust_mesh.get_element_vertices(0)
        assert len(verts) == 4, f"Expected 4 vertices, got {len(verts)}"

    def test_get_element_vertices_all_elements(self, rust_mesh_fixture):
        """Query all elements and verify returned vertices are valid."""
        rust_mesh, fixture_name = rust_mesh_fixture

        for elem_id in range(min(100, rust_mesh.n_elems)):
            verts = rust_mesh.get_element_vertices(elem_id)

            # Should be 4-element array
            assert len(verts) == 4, \
                f"Element {elem_id} in {fixture_name} should have 4 vertices, got {len(verts)}"

            # All should be valid vertex indices
            for v_id in verts:
                assert 0 <= v_id < rust_mesh.n_verts, \
                    f"Vertex {v_id} out of range [0, {rust_mesh.n_verts}) for element {elem_id} in {fixture_name}"

    def test_get_element_vertices_triangle_padding(self, rust_mesh_fixture):
        """For triangular elements, v3 should equal v2 (padding)."""
        rust_mesh, fixture_name = rust_mesh_fixture

        # Check first few elements
        for elem_id in range(min(20, rust_mesh.n_elems)):
            verts = rust_mesh.get_element_vertices(elem_id)
            v0, v1, v2, v3 = verts

            # Check if this is a padded triangle (any two adjacent or diagonal vertices equal)
            is_triangle = (v0 == v1) or (v1 == v2) or (v2 == v3) or (v3 == v0) or (v0 == v2) or (v1 == v3)

            if is_triangle:
                # For triangles, v3 should equal one of the first three
                assert v3 == v0 or v3 == v1 or v3 == v2, \
                    f"Triangle element {elem_id} has non-padded fourth vertex: {verts}"


class TestQueryConsistency:
    """Test consistency between different query methods."""

    def test_vertex_elements_contains_vertex(self, rust_mesh_fixture):
        """Every element from get_vertex_elements should contain the vertex."""
        rust_mesh, fixture_name = rust_mesh_fixture

        for v in range(min(50, rust_mesh.n_verts)):
            elem_list = rust_mesh.get_vertex_elements(v)

            for elem_id in elem_list:
                elem_verts = rust_mesh.get_element_vertices(elem_id)
                assert v in elem_verts, \
                    f"Vertex {v} returned in get_vertex_elements({v}) but not in element {elem_id} (verts={list(elem_verts)}) in {fixture_name}"

    def test_vertex_edges_exist_in_edge2vert(self, rust_mesh_fixture):
        """Every edge from get_vertex_edges should be in edge2vert."""
        rust_mesh, fixture_name = rust_mesh_fixture
        rust_mesh.build_adjacencies()

        # Build edge2vert for verification
        edge2vert = rust_mesh.get_edge2vert()

        for v in range(min(50, rust_mesh.n_verts)):
            edge_list = rust_mesh.get_vertex_edges(v)

            for edge_id in edge_list:
                assert edge_id < edge2vert.shape[0], \
                    f"Edge ID {edge_id} out of range for vertex {v} in {fixture_name}"

                v0, v1 = edge2vert[edge_id, 0], edge2vert[edge_id, 1]
                assert v == v0 or v == v1, \
                    f"Vertex {v} not in edge {edge_id} ({v0}, {v1}) in {fixture_name}"

    def test_element_vertices_consistency(self, rust_mesh_fixture):
        """get_element_vertices should match connectivity data."""
        rust_mesh, fixture_name = rust_mesh_fixture

        # Get connectivity from adjacency converters
        elem2vert = rust_mesh.get_elem2vert()

        for elem_id in range(min(50, rust_mesh.n_elems)):
            verts = rust_mesh.get_element_vertices(elem_id)

            # First 3 should match (or all 4 if quad)
            for i in range(min(4, elem2vert.shape[1])):
                expected_v = elem2vert[elem_id, i]
                actual_v = verts[i]
                assert expected_v == actual_v, \
                    f"Element {elem_id} vertex {i} mismatch: expected {expected_v}, got {actual_v} in {fixture_name}"


class TestQualityMetrics:
    """Test quality metrics derived from signed areas."""

    def test_area_statistics(self, rust_mesh_fixture):
        """Compute and verify area statistics."""
        rust_mesh, fixture_name = rust_mesh_fixture
        rust_mesh.compute_quality()
        areas = rust_mesh.get_signed_areas()

        # All areas should be finite
        assert np.all(np.isfinite(areas)), f"Found non-finite areas in {fixture_name}"

        # Compute statistics
        min_area = np.min(areas)
        max_area = np.max(areas)
        median_area = np.median(areas)
        mean_area = np.mean(areas)

        # Basic sanity checks
        assert min_area <= median_area <= max_area, \
            f"Median not between min and max in {fixture_name}"
        assert 1e-12 < abs(mean_area) < 1e6, \
            f"Mean area {mean_area} seems out of range for {fixture_name}"

    def test_positive_areas_after_ccw(self, rust_mesh_fixture):
        """After CCW enforcement, all areas should be positive."""
        rust_mesh, fixture_name = rust_mesh_fixture
        rust_mesh.ensure_ccw()
        areas = rust_mesh.get_signed_areas()

        # All should be positive
        assert np.all(areas > 0.0), \
            f"Found {np.sum(areas <= 0)} non-positive areas after CCW in {fixture_name}"

        # Check statistics
        min_area = np.min(areas)
        max_area = np.max(areas)
        assert 0 < min_area <= max_area, \
            f"Min/max area invalid in {fixture_name}: min={min_area}, max={max_area}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
