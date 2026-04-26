"""Tests for metadata validation and contributor workflows (Issue #44)."""
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pytest

import chilmesh


class TestMetadataValidation:
    """Test metadata validation for contributors (User Story 3)."""

    def test_metadata_completeness_for_contribution(self):
        """Test that metadata dict is complete for contribution validation."""
        mesh = chilmesh.examples.quad_2x2()
        metadata = mesh.admesh_metadata()

        # All 4 fields required for contribution
        required_fields = {"node_count", "element_count", "element_type", "bounding_box"}
        assert required_fields.issubset(metadata.keys()), \
            "Metadata missing required fields for contribution"

    def test_metadata_values_verifiable(self):
        """
        Test that metadata values are verifiable against manifest.

        Contributor workflow: Load mesh, call admesh_metadata(), compare returned
        values against what they entered in the manifest.
        """
        mesh = chilmesh.examples.quad_2x2()
        metadata = mesh.admesh_metadata()

        # Simulated manifest entry (what contributor entered)
        manifest_entry = {
            "node_count": 9,
            "element_count": 4,
            "element_type": "Quadrilateral",
            "bounding_box": {"min_x": 0.0, "max_x": 2.0, "min_y": 0.0, "max_y": 2.0}
        }

        # Verify all match
        assert metadata["node_count"] == manifest_entry["node_count"], "node_count mismatch"
        assert metadata["element_count"] == manifest_entry["element_count"], "element_count mismatch"
        assert metadata["element_type"] == manifest_entry["element_type"], "element_type mismatch"
        for key in ["min_x", "max_x", "min_y", "max_y"]:
            assert metadata["bounding_box"][key] == manifest_entry["bounding_box"][key], \
                f"bounding_box[{key}] mismatch"

    def test_metadata_mismatch_detectable(self):
        """Test that mismatches between metadata and manifest are easily detectable."""
        mesh = chilmesh.examples.quad_2x2()
        metadata = mesh.admesh_metadata()

        # Simulated wrong manifest entry
        wrong_manifest = {
            "node_count": 10,  # Wrong!
            "element_count": 4,
            "element_type": "Quadrilateral",
        }

        # Mismatch should be obvious
        assert metadata["node_count"] != wrong_manifest["node_count"], \
            "Should detect node_count mismatch"

        # In a real scenario, contributor sees:
        # Manifest says: node_count = 10
        # CHILmesh says: node_count = 9
        # → Clear mismatch for correction

    def test_element_type_enum_validation(self):
        """Test that element_type is always one of the valid enum values."""
        for fixture_func in [chilmesh.examples.annulus, chilmesh.examples.donut,
                             chilmesh.examples.block_o, chilmesh.examples.structured,
                             chilmesh.examples.quad_2x2]:
            mesh = fixture_func()
            metadata = mesh.admesh_metadata()

            valid_types = {"Triangular", "Quadrilateral", "Mixed-Element"}
            assert metadata["element_type"] in valid_types, \
                f"Invalid element_type '{metadata['element_type']}' from {fixture_func.__name__}"

    def test_contributor_preparation_workflow(self):
        """
        Test full contributor preparation workflow:
        1. Load mesh
        2. Get metadata
        3. Compare with manifest
        """
        # Contributor loads their mesh file
        mesh_file = chilmesh.examples.fixture_path("quad_2x2.fort.14")
        mesh = chilmesh.CHILmesh.read_from_fort14(mesh_file)

        # Get CHILmesh-computed metadata
        computed = mesh.admesh_metadata()

        # What they plan to enter in the manifest
        planned_manifest = {
            "node_count": 9,
            "element_count": 4,
            "element_type": "Quadrilateral",
            "bounding_box": {"min_x": 0.0, "max_x": 2.0, "min_y": 0.0, "max_y": 2.0}
        }

        # Quick validation: do they match?
        assert computed["node_count"] == planned_manifest["node_count"]
        assert computed["element_count"] == planned_manifest["element_count"]
        assert computed["element_type"] == planned_manifest["element_type"]

    def test_wrong_element_type_caught(self):
        """Test that wrong element_type in metadata is easily caught."""
        mesh = chilmesh.examples.quad_2x2()
        metadata = mesh.admesh_metadata()

        # Contributor mistakenly writes "Mixed-Element" in manifest
        wrong_type = "Mixed-Element"

        # CHILmesh says "Quadrilateral"
        assert metadata["element_type"] != wrong_type, \
            "Should detect element_type error"

    def test_bounding_box_precision(self):
        """Test that bounding_box values are precise enough for validation."""
        mesh = chilmesh.examples.quad_2x2()
        metadata = mesh.admesh_metadata()
        bbox = metadata["bounding_box"]

        # Values should be exact floats, not rounded
        assert bbox["min_x"] == 0.0
        assert bbox["max_x"] == 2.0
        assert bbox["min_y"] == 0.0
        assert bbox["max_y"] == 2.0

    def test_metadata_via_from_admesh_domain_for_validation(self):
        """Test contributor validation via from_admesh_domain entry point."""
        # Contributor receives a mesh record from ADMESH-Domains
        fixture_path = chilmesh.examples.fixture_path("quad_2x2.fort.14")
        record = SimpleNamespace(filename=str(fixture_path), type="ADCIRC")

        # They want to validate the metadata before contribution
        mesh = chilmesh.CHILmesh.from_admesh_domain(record)
        metadata = mesh.admesh_metadata()

        # All fields should be present and accurate
        assert metadata["node_count"] == 9
        assert metadata["element_count"] == 4
        assert metadata["element_type"] == "Quadrilateral"
        assert len(metadata["bounding_box"]) == 4


class TestGetLayerGuard:
    """Test improved error messages for get_layer() when layers not computed."""

    def test_get_layer_error_message_clarity(self):
        """Test that get_layer() error is clear and actionable."""
        fixture_path = chilmesh.examples.fixture_path("quad_2x2.fort.14")
        mesh = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=False)

        with pytest.raises(RuntimeError) as exc_info:
            mesh.get_layer(0)

        error_msg = str(exc_info.value)
        assert "Layers not computed" in error_msg, "Error should mention layers not computed"
        assert "compute_layers=True" in error_msg, "Error should mention compute_layers=True"
        assert "Re-initialise" in error_msg or "re-initialise" in error_msg, \
            "Error should suggest re-initialization"

    def test_get_layer_error_before_range_check(self):
        """Test that no-layers error is raised before range validation."""
        fixture_path = chilmesh.examples.fixture_path("quad_2x2.fort.14")
        mesh = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=False)

        # Even with an invalid index, should get the "layers not computed" error first
        with pytest.raises(RuntimeError, match="Layers not computed"):
            mesh.get_layer(999)

    def test_get_layer_valid_index_with_layers(self):
        """Test that get_layer still works normally when layers are computed."""
        mesh = chilmesh.examples.annulus()  # Default: compute_layers=True

        # Should work fine
        layer = mesh.get_layer(0)
        assert "OE" in layer
        assert "IE" in layer
        assert "OV" in layer
        assert "IV" in layer
        assert "bEdgeIDs" in layer

    def test_get_layer_invalid_index_with_layers(self):
        """Test that range validation still works when layers are computed."""
        mesh = chilmesh.examples.annulus()

        with pytest.raises(ValueError, match="out of range"):
            mesh.get_layer(999)
