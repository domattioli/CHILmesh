"""Integration tests for MADMESHR advancing-front API.

Tests validate that advancing-front element addition maintains mesh consistency
and produces valid results suitable for MADMESHR mesh generation.
"""

import pytest
import numpy as np
from pathlib import Path
from chilmesh import CHILmesh


class TestAdvancingFrontBoundary:
    """Tests for advancing_front_boundary_edges() method."""

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_boundary_edges_match_direct_call(self, fixture_name):
        """Verify advancing_front_boundary_edges matches boundary_edges()."""
        mesh = self._load_fixture(fixture_name)

        direct_boundary = set(mesh.boundary_edges())
        af_boundary = set(mesh.advancing_front_boundary_edges())

        assert direct_boundary == af_boundary

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_boundary_edges_sorted(self, fixture_name):
        """Verify advancing_front_boundary_edges returns sorted list."""
        mesh = self._load_fixture(fixture_name)

        boundary = mesh.advancing_front_boundary_edges()

        assert isinstance(boundary, list)
        assert boundary == sorted(boundary)

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_boundary_edges_valid_range(self, fixture_name):
        """Verify all boundary edges are valid indices."""
        mesh = self._load_fixture(fixture_name)

        boundary = mesh.advancing_front_boundary_edges()

        for edge_id in boundary:
            assert 0 <= edge_id < mesh.n_edges

    def _load_fixture(self, name):
        """Load a test fixture mesh."""
        if name == "annulus":
            return CHILmesh.read_from_fort14(
                Path("src/chilmesh/data/annulus_200pts.fort.14")
            )
        elif name == "donut":
            return CHILmesh.read_from_fort14(
                Path("src/chilmesh/data/donut_domain.fort.14")
            )
        elif name == "structured":
            return CHILmesh.read_from_fort14(
                Path("src/chilmesh/data/quad_2x2.fort.14")
            )
        else:
            raise ValueError(f"Unknown fixture: {name}")


class TestAddAdvancingFrontElement:
    """Tests for add_advancing_front_element() method."""

    def test_add_triangle_element(self):
        """Verify adding a triangle element works correctly."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )
        initial_elems = mesh.n_elems
        elem_cols = mesh.connectivity_list.shape[1]

        # Add a triangle
        new_elem_id = mesh.add_advancing_front_element([0, 1, 2], elem_type="tri")

        assert new_elem_id == initial_elems
        assert mesh.n_elems == initial_elems + 1

        # Verify element was added
        elem = mesh.connectivity_list[new_elem_id]
        assert elem[0] == 0 and elem[1] == 1 and elem[2] == 2

        # Check padding based on mesh format
        if elem_cols == 4:
            assert elem[3] == 0  # Padded triangle in 4-column format

    def test_add_quad_element(self):
        """Verify adding a quad element works correctly."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/quad_2x2.fort.14"), compute_layers=False
        )
        initial_elems = mesh.n_elems
        elem_cols = mesh.connectivity_list.shape[1]

        # Verify mesh supports quads (4-column format)
        if elem_cols == 4:
            # Add a quad
            new_elem_id = mesh.add_advancing_front_element([0, 1, 2, 3], elem_type="quad")

            assert new_elem_id == initial_elems
            assert mesh.n_elems == initial_elems + 1

            # Verify element was added
            elem = mesh.connectivity_list[new_elem_id]
            assert list(elem) == [0, 1, 2, 3]
        else:
            # Mesh is triangles-only, skip quad test
            pytest.skip("Mesh does not support quad elements")

    def test_add_multiple_elements(self):
        """Verify adding multiple elements maintains consistency."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )
        initial_elems = mesh.n_elems

        # Add 5 triangles
        elem_ids = []
        for i in range(5):
            v1 = i % mesh.n_verts
            v2 = (i + 1) % mesh.n_verts
            v3 = (i + 2) % mesh.n_verts
            elem_id = mesh.add_advancing_front_element([v1, v2, v3], elem_type="tri")
            elem_ids.append(elem_id)

        # Verify all elements were added
        assert mesh.n_elems == initial_elems + 5

        # Verify element IDs are sequential
        assert elem_ids == list(range(initial_elems, initial_elems + 5))

        # Verify adjacencies are valid
        for elem_id in elem_ids:
            edges = mesh.elem2edge(elem_id)[0]  # Extract first row
            for edge_id in edges:
                assert 0 <= edge_id < mesh.n_edges

    def test_add_element_invalid_type(self):
        """Verify invalid element type raises error."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )

        with pytest.raises(ValueError):
            mesh.add_advancing_front_element([0, 1, 2], elem_type="invalid")

    def test_add_element_invalid_vertex(self):
        """Verify invalid vertex index raises error."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )

        with pytest.raises(ValueError):
            mesh.add_advancing_front_element([0, 1, mesh.n_verts], elem_type="tri")

    def test_add_element_wrong_vertex_count(self):
        """Verify wrong vertex count raises error."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )

        # Triangle with 4 vertices
        with pytest.raises(ValueError):
            mesh.add_advancing_front_element([0, 1, 2, 3], elem_type="tri")

        # Quad with 3 vertices
        with pytest.raises(ValueError):
            mesh.add_advancing_front_element([0, 1, 2], elem_type="quad")

    def test_add_element_invalidates_layers(self):
        """Verify adding elements clears cached layers."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14"), compute_layers=True
        )

        initial_layers = len(mesh.layers)
        if initial_layers == 0:
            pytest.skip("Mesh has no layers to test invalidation")

        # Add element
        mesh.add_advancing_front_element([0, 1, 2], elem_type="tri")

        # Layers should be cleared
        assert len(mesh.layers) == 0


class TestRemoveBoundaryLoop:
    """Tests for remove_boundary_loop() method."""

    def test_remove_single_boundary_element(self):
        """Verify removing a single boundary element works."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )
        initial_elems = mesh.n_elems

        # Get a boundary edge
        boundary = mesh.advancing_front_boundary_edges()
        if len(boundary) > 0:
            edge_id = boundary[0]

            # Remove the element adjacent to this edge
            mesh.remove_boundary_loop([edge_id])

            # Verify element count decreased
            assert mesh.n_elems == initial_elems - 1

    def test_remove_multiple_boundary_elements(self):
        """Verify removing multiple boundary elements works."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )
        initial_elems = mesh.n_elems

        # Get first 3 boundary edges
        boundary = mesh.advancing_front_boundary_edges()[:3]
        if len(boundary) > 0:
            mesh.remove_boundary_loop(boundary)

            # Verify element count decreased
            assert mesh.n_elems < initial_elems

    def test_remove_invalid_edge(self):
        """Verify invalid edge index raises error."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )

        with pytest.raises(ValueError):
            mesh.remove_boundary_loop([mesh.n_edges])

    def test_remove_boundary_invalidates_layers(self):
        """Verify removing elements clears cached layers."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )

        assert len(mesh.layers) > 0

        boundary = mesh.advancing_front_boundary_edges()
        if len(boundary) > 0:
            mesh.remove_boundary_loop([boundary[0]])

            # Layers should be cleared
            assert len(mesh.layers) == 0


class TestPinchPoints:
    """Tests for pinch_points() method."""

    def test_pinch_points_returns_list(self):
        """Verify pinch_points returns a list."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )

        pinches = mesh.pinch_points(width_threshold=0.5)

        assert isinstance(pinches, list)

    def test_pinch_points_sorted(self):
        """Verify pinch_points returns sorted list."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )

        pinches = mesh.pinch_points(width_threshold=0.5)

        assert pinches == sorted(pinches)

    def test_pinch_points_valid_range(self):
        """Verify all pinch points are valid vertex indices."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )

        pinches = mesh.pinch_points(width_threshold=0.5)

        for vert_id in pinches:
            assert 0 <= vert_id < mesh.n_verts

    def test_pinch_points_threshold_effect(self):
        """Verify threshold affects number of detected pinches."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )

        low_threshold = mesh.pinch_points(width_threshold=0.1)
        high_threshold = mesh.pinch_points(width_threshold=0.9)

        # Higher threshold should detect more pinch points
        assert len(high_threshold) >= len(low_threshold)


class TestAdvancingFrontWorkflow:
    """Test realistic advancing-front generation scenarios."""

    def test_small_fill_center_scenario(self):
        """Simulate advancing-front filling center of annulus."""
        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )

        initial_elems = mesh.n_elems
        initial_verts = mesh.n_verts

        # Add 10 elements
        for i in range(10):
            v1 = i % (initial_verts - 2)
            v2 = (i + 1) % (initial_verts - 2)
            v3 = (i + 2) % (initial_verts - 2)
            mesh.add_advancing_front_element([v1, v2, v3], elem_type="tri")

        # Verify state
        assert mesh.n_elems == initial_elems + 10

        # Verify mesh is still valid
        edges = mesh.elem2edge()
        assert edges.shape[0] == mesh.n_elems

    def test_export_after_modification(self):
        """Verify mesh can be exported after advancing-front modification."""
        import tempfile

        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )

        # Add an element
        mesh.add_advancing_front_element([0, 1, 2], elem_type="tri")

        # Export to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fort.14', delete=False) as f:
            temp_path = f.name

        try:
            success = mesh.write_to_fort14(temp_path)
            assert success

            # Re-import and verify
            mesh2 = CHILmesh.read_from_fort14(Path(temp_path), compute_layers=False)
            assert mesh2.n_elems == mesh.n_elems
            assert mesh2.n_verts == mesh.n_verts
        finally:
            import os
            if os.path.exists(temp_path):
                os.remove(temp_path)
