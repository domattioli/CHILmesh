"""Tests for from_admesh_domain entry point and metadata (Issues #42–43)."""
from pathlib import Path
from types import SimpleNamespace

import pytest

import chilmesh


class TestFromADMESHDomain:
    """Test CHILmesh.from_admesh_domain() entry point."""

    def test_adcirc_routing(self):
        """Test that type='ADCIRC' routes to read_from_fort14."""
        fixture_path = chilmesh.examples.fixture_path("quad_2x2.fort.14")
        record = SimpleNamespace(filename=str(fixture_path), type="ADCIRC")

        mesh = chilmesh.CHILmesh.from_admesh_domain(record)
        assert mesh.n_verts == 9
        assert mesh.n_elems == 4
        assert mesh.type == "Quadrilateral"

    def test_adcirc_grd_routing(self):
        """Test that type='ADCIRC_GRD' routes to read_from_fort14."""
        fixture_path = chilmesh.examples.fixture_path("quad_2x2.fort.14")
        record = SimpleNamespace(filename=str(fixture_path), type="ADCIRC_GRD")

        mesh = chilmesh.CHILmesh.from_admesh_domain(record)
        assert mesh.n_verts == 9
        assert mesh.n_elems == 4

    def test_unknown_type_fallback(self):
        """Test that unknown type falls back to ADCIRC reader with warning."""
        fixture_path = chilmesh.examples.fixture_path("quad_2x2.fort.14")
        record = SimpleNamespace(filename=str(fixture_path), type="UNKNOWN_FORMAT")

        with pytest.warns(UserWarning, match="Unrecognised mesh type"):
            mesh = chilmesh.CHILmesh.from_admesh_domain(record)

        assert mesh.n_verts == 9

    def test_none_type_defaults_to_adcirc(self):
        """Test that type=None defaults to ADCIRC reader."""
        fixture_path = chilmesh.examples.fixture_path("quad_2x2.fort.14")
        record = SimpleNamespace(filename=str(fixture_path), type=None)

        mesh = chilmesh.CHILmesh.from_admesh_domain(record)
        assert mesh.n_verts == 9

    def test_missing_file_error(self):
        """Test that missing file raises FileNotFoundError with guidance."""
        record = SimpleNamespace(filename="/nonexistent/path/mesh.fort.14", type="ADCIRC")

        with pytest.raises(FileNotFoundError) as exc_info:
            chilmesh.CHILmesh.from_admesh_domain(record)

        assert "File not found" in str(exc_info.value)
        assert "mesh_record.load()" in str(exc_info.value)

    def test_compute_layers_false_forwarded(self):
        """Test that compute_layers=False is forwarded to reader."""
        fixture_path = chilmesh.examples.fixture_path("quad_2x2.fort.14")
        record = SimpleNamespace(filename=str(fixture_path), type="ADCIRC")

        mesh = chilmesh.CHILmesh.from_admesh_domain(record, compute_layers=False)
        assert mesh.n_layers == 0
        assert mesh.adjacencies == {}


class TestADMESHMetadata:
    """Test CHILmesh.admesh_metadata() method."""

    def test_metadata_basic(self):
        """Test that metadata dict has required keys and correct values."""
        mesh = chilmesh.examples.quad_2x2()
        metadata = mesh.admesh_metadata()

        assert isinstance(metadata, dict)
        assert "node_count" in metadata
        assert "element_count" in metadata
        assert "element_type" in metadata
        assert "bounding_box" in metadata

    def test_metadata_values_quad_2x2(self):
        """Test metadata values for quad_2x2 fixture."""
        mesh = chilmesh.examples.quad_2x2()
        metadata = mesh.admesh_metadata()

        assert metadata["node_count"] == 9
        assert metadata["element_count"] == 4
        assert metadata["element_type"] == "Quadrilateral"

    def test_metadata_bounding_box(self):
        """Test bounding_box values."""
        mesh = chilmesh.examples.quad_2x2()
        metadata = mesh.admesh_metadata()
        bbox = metadata["bounding_box"]

        # quad_2x2 spans 0-2 in both x and y
        assert bbox["min_x"] == pytest.approx(0.0)
        assert bbox["max_x"] == pytest.approx(2.0)
        assert bbox["min_y"] == pytest.approx(0.0)
        assert bbox["max_y"] == pytest.approx(2.0)

    def test_metadata_with_no_layers(self):
        """Test that metadata is available even when compute_layers=False."""
        fixture_path = chilmesh.examples.fixture_path("quad_2x2.fort.14")
        mesh = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=False)

        metadata = mesh.admesh_metadata()
        assert metadata["node_count"] == 9
        assert metadata["element_count"] == 4
        assert metadata["element_type"] == "Quadrilateral"

    def test_metadata_triangular_mesh(self):
        """Test metadata for pure triangular mesh."""
        mesh = chilmesh.examples.annulus()
        metadata = mesh.admesh_metadata()

        assert metadata["element_type"] == "Triangular"
        assert metadata["node_count"] == mesh.n_verts
        assert metadata["element_count"] == mesh.n_elems
