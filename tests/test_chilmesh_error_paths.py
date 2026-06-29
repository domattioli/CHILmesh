"""
Test error-path guards in CHILmesh.py.

Covers ValueError, AssertionError, and FileNotFoundError conditions in:
- get_vertex_edges / get_vertex_elements (range checks)
- find_elements_in_radius (negative radius)
- elem_quality (invalid quality_type)
- smooth_mesh (invalid method, missing sdf, missing acknowledge_change)
- add_advancing_front_element (invalid elem_type, out-of-range vertices, wrong counts)
- CHILmesh.load (nonexistent file)
"""

import pytest
import numpy as np
import tempfile
from pathlib import Path
import chilmesh
from chilmesh import examples


class TestGetVertexEdgesErrorPaths:
    """Test get_vertex_edges() error guards."""

    def test_vertex_edges_negative_id(self):
        """Negative vertex ID raises ValueError."""
        mesh = examples.annulus()
        with pytest.raises(ValueError, match=r"Vertex -1 out of range"):
            mesh.get_vertex_edges(-1)

    def test_vertex_edges_out_of_range(self):
        """Vertex ID >= n_verts raises ValueError."""
        mesh = examples.annulus()
        n_verts = mesh.n_verts
        with pytest.raises(ValueError, match=f"Vertex {n_verts} out of range"):
            mesh.get_vertex_edges(n_verts)


class TestGetVertexElementsErrorPaths:
    """Test get_vertex_elements() error guards."""

    def test_vertex_elements_negative_id(self):
        """Negative vertex ID raises ValueError."""
        mesh = examples.annulus()
        with pytest.raises(ValueError, match=r"Vertex -1 out of range"):
            mesh.get_vertex_elements(-1)

    def test_vertex_elements_out_of_range(self):
        """Vertex ID >= n_verts raises ValueError."""
        mesh = examples.annulus()
        n_verts = mesh.n_verts
        with pytest.raises(ValueError, match=f"Vertex {n_verts} out of range"):
            mesh.get_vertex_elements(n_verts)


class TestFindElementsInRadiusErrorPaths:
    """Test find_elements_in_radius() error guards."""

    def test_negative_radius(self):
        """Negative radius raises ValueError."""
        mesh = examples.annulus()
        with pytest.raises(ValueError, match=r"radius must be >= 0"):
            mesh.find_elements_in_radius(np.array([0.0, 0.0]), -1.0)


class TestElemQualityErrorPaths:
    """Test elem_quality() error guards."""

    def test_unknown_quality_type(self):
        """Invalid quality_type raises ValueError."""
        mesh = examples.annulus()
        with pytest.raises(ValueError, match=r"Unknown quality type"):
            mesh.elem_quality(quality_type='bogus')


class TestSmoothMeshErrorPaths:
    """Test smooth_mesh() error guards."""

    def test_unknown_smoothing_method(self):
        """Invalid smoothing method raises ValueError."""
        mesh = examples.annulus().copy()
        with pytest.raises(ValueError, match=r"Unknown smoothing method"):
            mesh.smooth_mesh(method='bogus', acknowledge_change=True)

    def test_sdf_method_missing_sdf(self):
        """sdf method without sdf= callable raises ValueError."""
        mesh = examples.annulus().copy()
        with pytest.raises(ValueError, match=r"method='sdf' requires sdf="):
            mesh.smooth_mesh(method='sdf', acknowledge_change=True)

    def test_smooth_mesh_missing_acknowledge_change(self):
        """smooth_mesh without acknowledge_change=True raises AssertionError."""
        mesh = examples.annulus().copy()
        # The assert at the top of smooth_mesh checks acknowledge_change
        with pytest.raises(AssertionError):
            mesh.smooth_mesh(method='fem')


class TestAddAdvancingFrontElementErrorPaths:
    """Test add_advancing_front_element() error guards."""

    def test_invalid_elem_type(self):
        """Invalid elem_type raises ValueError."""
        mesh = examples.annulus().copy()
        with pytest.raises(ValueError, match=r"Invalid elem_type"):
            mesh.add_advancing_front_element([0, 1, 2], 'pentagon')

    def test_vertex_out_of_range(self):
        """Vertex ID out of range raises ValueError."""
        mesh = examples.annulus().copy()
        with pytest.raises(ValueError, match=r"Vertex \d+ out of range"):
            mesh.add_advancing_front_element([0, 1, 999999], 'tri')

    def test_triangle_wrong_vertex_count(self):
        """Triangle with wrong vertex count raises ValueError."""
        mesh = examples.annulus().copy()
        # 4 vertices for a triangle should fail
        with pytest.raises(ValueError, match=r"Triangle requires 3 vertices"):
            mesh.add_advancing_front_element([0, 1, 2, 3], 'tri')

    def test_quad_wrong_vertex_count(self):
        """Quad with wrong vertex count raises ValueError."""
        mesh = examples.annulus().copy()
        # 3 vertices for a quad should fail
        with pytest.raises(ValueError, match=r"Quad requires 4 vertices"):
            mesh.add_advancing_front_element([0, 1, 2], 'quad')

    def test_quad_into_triangle_only_mesh(self):
        """Adding quad to triangle-only mesh raises ValueError."""
        # annulus has 3-column connectivity (triangle-only)
        mesh = examples.annulus().copy()
        assert mesh.connectivity_list.shape[1] == 3, "Expected triangle-only mesh"
        with pytest.raises(ValueError, match=r"Cannot add quad to triangle-only mesh"):
            mesh.add_advancing_front_element([0, 1, 2, 3], 'quad')

    def test_add_triangle_to_quad_mesh(self):
        """Adding triangle to quad mesh succeeds (pads to 4-column)."""
        # quad_2x2 is a quad-only (4-column) mesh
        mesh = examples.quad_2x2().copy()
        assert mesh.connectivity_list.shape[1] == 4, "Expected quad mesh"
        n_elems_before = mesh.n_elems
        # Triangle should be padded to 4-column
        new_elem_id = mesh.add_advancing_front_element([0, 1, 2], 'tri')
        assert mesh.n_elems == n_elems_before + 1
        assert new_elem_id == n_elems_before


class TestLoadFileErrorPaths:
    """Test CHILmesh.load() error guards."""

    def test_load_nonexistent_file(self):
        """Loading nonexistent file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match=r"File not found"):
            chilmesh.CHILmesh.load('/nonexistent/path/nope.14')

    def test_load_file_not_found_with_2dm_extension(self):
        """Loading nonexistent .2dm file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match=r"File not found"):
            chilmesh.CHILmesh.load('/nonexistent/path/nope.2dm')
