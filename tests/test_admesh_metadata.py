"""Tests for admesh_metadata() accuracy (Issue #42)."""
import pytest

import chilmesh


class TestADMESHMetadataAccuracy:
    """Test that admesh_metadata() returns accurate, ADMESH-Domains-compatible values."""

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured", "quad_2x2"])
    def test_metadata_completeness(self, fixture_name):
        """
        Test that metadata dict has all required keys.

        Success Criteria (SC-004): admesh_metadata() returns all 4 schema fields.
        """
        mesh = chilmesh.examples.__dict__[fixture_name]()
        metadata = mesh.admesh_metadata()

        required_keys = {"node_count", "element_count", "element_type", "bounding_box"}
        assert required_keys.issubset(metadata.keys()), \
            f"Missing keys in metadata. Expected {required_keys}, got {metadata.keys()}"

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured", "quad_2x2"])
    def test_metadata_node_count_accuracy(self, fixture_name):
        """Test that node_count matches actual mesh node count."""
        mesh = chilmesh.examples.__dict__[fixture_name]()
        metadata = mesh.admesh_metadata()

        assert metadata["node_count"] == mesh.n_verts, \
            f"node_count mismatch: {metadata['node_count']} != {mesh.n_verts}"

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured", "quad_2x2"])
    def test_metadata_element_count_accuracy(self, fixture_name):
        """Test that element_count matches actual mesh element count."""
        mesh = chilmesh.examples.__dict__[fixture_name]()
        metadata = mesh.admesh_metadata()

        assert metadata["element_count"] == mesh.n_elems, \
            f"element_count mismatch: {metadata['element_count']} != {mesh.n_elems}"

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured", "quad_2x2"])
    def test_metadata_element_type_canonical(self, fixture_name):
        """Test that element_type is one of the canonical values."""
        mesh = chilmesh.examples.__dict__[fixture_name]()
        metadata = mesh.admesh_metadata()

        canonical_types = {"Triangular", "Quadrilateral", "Mixed-Element"}
        assert metadata["element_type"] in canonical_types, \
            f"element_type '{metadata['element_type']}' not in {canonical_types}"

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured", "quad_2x2"])
    def test_metadata_element_type_matches_mesh_type(self, fixture_name):
        """Test that metadata element_type matches mesh.type."""
        mesh = chilmesh.examples.__dict__[fixture_name]()
        metadata = mesh.admesh_metadata()

        assert metadata["element_type"] == mesh.type, \
            f"element_type mismatch: {metadata['element_type']} != {mesh.type}"

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured", "quad_2x2"])
    def test_metadata_bounding_box_keys(self, fixture_name):
        """Test that bounding_box has correct keys."""
        mesh = chilmesh.examples.__dict__[fixture_name]()
        metadata = mesh.admesh_metadata()
        bbox = metadata["bounding_box"]

        required_keys = {"min_x", "max_x", "min_y", "max_y"}
        assert required_keys.issubset(bbox.keys()), \
            f"Missing keys in bounding_box. Expected {required_keys}, got {bbox.keys()}"

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured", "quad_2x2"])
    def test_metadata_bounding_box_validity(self, fixture_name):
        """Test that bounding_box values are sensible (min <= max)."""
        mesh = chilmesh.examples.__dict__[fixture_name]()
        metadata = mesh.admesh_metadata()
        bbox = metadata["bounding_box"]

        assert bbox["min_x"] <= bbox["max_x"], \
            f"Invalid bounding_box: min_x {bbox['min_x']} > max_x {bbox['max_x']}"
        assert bbox["min_y"] <= bbox["max_y"], \
            f"Invalid bounding_box: min_y {bbox['min_y']} > max_y {bbox['max_y']}"

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured", "quad_2x2"])
    def test_metadata_bounding_box_matches_points(self, fixture_name):
        """Test that bounding_box matches actual point extrema."""
        mesh = chilmesh.examples.__dict__[fixture_name]()
        metadata = mesh.admesh_metadata()
        bbox = metadata["bounding_box"]

        points = mesh.points
        expected_min_x = points[:, 0].min()
        expected_max_x = points[:, 0].max()
        expected_min_y = points[:, 1].min()
        expected_max_y = points[:, 1].max()

        assert bbox["min_x"] == pytest.approx(expected_min_x), \
            f"min_x mismatch: {bbox['min_x']} != {expected_min_x}"
        assert bbox["max_x"] == pytest.approx(expected_max_x), \
            f"max_x mismatch: {bbox['max_x']} != {expected_max_x}"
        assert bbox["min_y"] == pytest.approx(expected_min_y), \
            f"min_y mismatch: {bbox['min_y']} != {expected_min_y}"
        assert bbox["max_y"] == pytest.approx(expected_max_y), \
            f"max_y mismatch: {bbox['max_y']} != {expected_max_y}"

    def test_metadata_with_compute_layers_false(self):
        """Test that metadata is available even when compute_layers=False (SC-004)."""
        fixture_path = chilmesh.examples.fixture_path("quad_2x2.fort.14")
        mesh = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=False)

        # Metadata should still be available
        metadata = mesh.admesh_metadata()
        assert metadata["node_count"] == 9
        assert metadata["element_count"] == 4
        assert metadata["element_type"] == "Quadrilateral"
        assert "bounding_box" in metadata

    def test_metadata_with_from_admesh_domain(self):
        """Test metadata via from_admesh_domain entry point."""
        from types import SimpleNamespace

        fixture_path = chilmesh.examples.fixture_path("quad_2x2.fort.14")
        record = SimpleNamespace(filename=str(fixture_path), type="ADCIRC")

        mesh = chilmesh.CHILmesh.from_admesh_domain(record)
        metadata = mesh.admesh_metadata()

        assert metadata["node_count"] == 9
        assert metadata["element_count"] == 4
        assert metadata["element_type"] == "Quadrilateral"
        assert metadata["bounding_box"]["min_x"] == pytest.approx(0.0)
        assert metadata["bounding_box"]["max_x"] == pytest.approx(2.0)

    def test_metadata_triangular_vs_quadrilateral(self):
        """Test that triangular and quad meshes are correctly classified."""
        # Triangular fixture
        tri_mesh = chilmesh.examples.annulus()
        tri_metadata = tri_mesh.admesh_metadata()
        assert tri_metadata["element_type"] == "Triangular"

        # Quadrilateral fixture
        quad_mesh = chilmesh.examples.quad_2x2()
        quad_metadata = quad_mesh.admesh_metadata()
        assert quad_metadata["element_type"] == "Quadrilateral"
