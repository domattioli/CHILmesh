"""Tests for chilmesh.summary_io module.

Cover the summary() function for both file paths (fort14, 2dm) and CHILmesh
objects, with shallow and deep reading modes.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from chilmesh import CHILmesh, examples
from chilmesh.summary_io import summary, SummaryError


class TestSummaryFort14Shallow:
    """Test fort14 shallow summary (header-only reading)."""

    def test_fort14_quad_2x2_shallow(self):
        """Fort14 shallow summary of quad_2x2: basic fields present, no deep fields."""
        path = examples.fixture_path("quad_2x2.fort.14")
        result = summary(path, deep=False)

        # Shallow-mode required fields
        assert result['n_elems'] == 4, f"Expected n_elems=4, got {result['n_elems']}"
        assert result['n_nodes'] == 9, f"Expected n_nodes=9, got {result['n_nodes']}"
        assert result['format'] == 'fort14', f"Expected format='fort14', got {result['format']}"
        assert result['grid_name'] == '2x2 quad test mesh', (
            f"Expected grid_name='2x2 quad test mesh', got {result['grid_name']}"
        )
        assert result['file_bytes'] > 0, f"Expected file_bytes>0, got {result['file_bytes']}"
        assert str(path) in result['path'], (
            f"Expected path to contain {path}, got {result['path']}"
        )

        # Shallow mode: these should NOT be present
        assert 'element_type' not in result, (
            "element_type should not be in shallow summary"
        )
        assert 'bbox' not in result, (
            "bbox should not be in shallow summary"
        )

    def test_fort14_donut_shallow(self):
        """Fort14 shallow summary of donut fixture: count sanity."""
        path = examples.fixture_path("donut_domain.fort.14")
        result = summary(path, deep=False)

        assert result['n_elems'] == 276, f"Expected n_elems=276, got {result['n_elems']}"
        assert result['n_nodes'] == 188, f"Expected n_nodes=188, got {result['n_nodes']}"
        assert result['format'] == 'fort14', f"Expected format='fort14', got {result['format']}"


class TestSummaryFort14Deep:
    """Test fort14 deep summary (loads full mesh for element_type and bbox)."""

    def test_fort14_quad_2x2_deep(self):
        """Fort14 deep summary: includes element_type and bbox."""
        path = examples.fixture_path("quad_2x2.fort.14")
        result = summary(path, deep=True)

        # Deep mode adds element_type and bbox
        assert result['element_type'] == 'Quadrilateral', (
            f"Expected element_type='Quadrilateral', got {result['element_type']}"
        )

        # Bbox should be present with 4 keys
        assert 'bbox' in result, "bbox should be in deep summary"
        bbox = result['bbox']
        assert 'min_x' in bbox, f"bbox missing min_x: {bbox}"
        assert 'max_x' in bbox, f"bbox missing max_x: {bbox}"
        assert 'min_y' in bbox, f"bbox missing min_y: {bbox}"
        assert 'max_y' in bbox, f"bbox missing max_y: {bbox}"

        # Bbox bounds sanity
        assert bbox['min_x'] <= bbox['max_x'], (
            f"min_x={bbox['min_x']} > max_x={bbox['max_x']}"
        )
        assert bbox['min_y'] <= bbox['max_y'], (
            f"min_y={bbox['min_y']} > max_y={bbox['max_y']}"
        )

        # Expected values for quad_2x2
        assert bbox['min_x'] == 0.0, f"Expected min_x=0.0, got {bbox['min_x']}"
        assert bbox['max_x'] == 2.0, f"Expected max_x=2.0, got {bbox['max_x']}"
        assert bbox['min_y'] == 0.0, f"Expected min_y=0.0, got {bbox['min_y']}"
        assert bbox['max_y'] == 2.0, f"Expected max_y=2.0, got {bbox['max_y']}"


class TestSummaryMeshObject:
    """Test summary() on CHILmesh objects (always returns full info)."""

    def test_mesh_object_summary(self):
        """CHILmesh object summary: includes element_type and bbox, path is None."""
        path = examples.fixture_path("quad_2x2.fort.14")
        mesh = CHILmesh.read_from_fort14(path)
        result = summary(mesh)

        # Mesh-object mode fields
        assert result['format'] == 'mesh-object', (
            f"Expected format='mesh-object', got {result['format']}"
        )
        assert result['path'] is None, f"Expected path=None, got {result['path']}"
        assert result['n_nodes'] == 9, f"Expected n_nodes=9, got {result['n_nodes']}"
        assert result['n_elems'] == 4, f"Expected n_elems=4, got {result['n_elems']}"

        # Always includes element_type and bbox
        assert result['element_type'], f"element_type must be truthy, got {result['element_type']}"
        assert 'bbox' in result, "bbox must be present for mesh objects"

        bbox = result['bbox']
        assert len(bbox) == 4, f"bbox should have 4 keys, got {len(bbox)}"
        assert bbox['min_x'] <= bbox['max_x'], (
            f"min_x={bbox['min_x']} > max_x={bbox['max_x']}"
        )
        assert bbox['min_y'] <= bbox['max_y'], (
            f"min_y={bbox['min_y']} > max_y={bbox['max_y']}"
        )


class TestSummary2dm:
    """Test summary() on 2dm format files."""

    def test_2dm_roundtrip(self, tmp_path):
        """Create a 2dm file from quad_2x2, then test summary() on it."""
        # Load quad_2x2 fort14, save to 2dm
        fort14_path = examples.fixture_path("quad_2x2.fort.14")
        mesh = CHILmesh.read_from_fort14(fort14_path)

        out2dm = tmp_path / "quad_2x2.2dm"
        mesh.save(str(out2dm))
        assert out2dm.exists(), f"2dm file not created at {out2dm}"

        # Now test summary on the 2dm file
        result = summary(out2dm)

        assert result['format'] == '2dm', f"Expected format='2dm', got {result['format']}"
        assert result['n_nodes'] == 9, f"Expected n_nodes=9, got {result['n_nodes']}"
        assert result['n_elems'] == 4, f"Expected n_elems=4, got {result['n_elems']}"
        assert result['file_bytes'] > 0, f"Expected file_bytes>0, got {result['file_bytes']}"


class TestSummaryErrors:
    """Test error cases for summary()."""

    def test_unknown_format_error(self, tmp_path):
        """Unknown file suffix raises SummaryError with helpful message."""
        bad_path = tmp_path / "mesh.xyz"
        bad_path.write_text("not a mesh")

        with pytest.raises(SummaryError) as exc_info:
            summary(bad_path)

        assert "unknown mesh format" in str(exc_info.value).lower(), (
            f"SummaryError message should mention 'unknown mesh format', got: {exc_info.value}"
        )

    def test_malformed_fort14_header_error(self, tmp_path):
        """Malformed fort14 header (2nd line with <2 tokens) raises SummaryError."""
        bad_fort14 = tmp_path / "bad.fort.14"
        # First line OK, second line has only 1 token
        bad_fort14.write_text("valid grid name\n4\n")

        with pytest.raises(SummaryError) as exc_info:
            summary(bad_fort14)

        assert "malformed" in str(exc_info.value).lower(), (
            f"SummaryError should mention 'malformed', got: {exc_info.value}"
        )

    def test_single_line_fort14_error(self, tmp_path):
        """Fort14 with only 1 line raises SummaryError."""
        bad_fort14 = tmp_path / "single_line.fort.14"
        bad_fort14.write_text("only one line\n")

        with pytest.raises(SummaryError) as exc_info:
            summary(bad_fort14)

        assert "malformed" in str(exc_info.value).lower(), (
            f"SummaryError should mention 'malformed', got: {exc_info.value}"
        )


class TestSummaryFort14SuffixVariant:
    """Test fort14 suffix variant detection (lines 94-95)."""

    def test_fort14_suffix_with_fort14_extension(self, tmp_path):
        """File with .fort14 suffix detects as fort14 format."""
        mesh_file = tmp_path / "mesh.fort14"
        mesh_file.write_text("grid\n4 9\n")

        result = summary(mesh_file)

        assert result['format'] == 'fort14', (
            f"Expected format='fort14' for .fort14 suffix, got {result['format']}"
        )
        assert result['n_elems'] == 4, f"Expected n_elems=4, got {result['n_elems']}"
        assert result['n_nodes'] == 9, f"Expected n_nodes=9, got {result['n_nodes']}"


class TestSummaryFileStatErrors:
    """Test file stat errors (lines 104-105)."""

    def test_nonexistent_file_stat_error(self, tmp_path):
        """Calling summary() on nonexistent file raises FileNotFoundError (issue #235)."""
        nonexistent = tmp_path / "nonexistent.14"

        with pytest.raises(FileNotFoundError) as exc_info:
            summary(nonexistent)

        assert "no such file" in str(exc_info.value).lower(), (
            f"FileNotFoundError should mention 'no such file', got: {exc_info.value}"
        )


class TestSummaryFort14HeaderParsing:
    """Test fort14 header parsing errors (lines 154-155)."""

    def test_fort14_non_integer_tokens_error(self, tmp_path):
        """Fort14 with non-integer tokens in header line raises SummaryError."""
        bad_fort14 = tmp_path / "bad_header.fort.14"
        bad_fort14.write_text("grid name\nfoo bar\n")

        with pytest.raises(SummaryError) as exc_info:
            summary(bad_fort14)

        assert "parse error" in str(exc_info.value).lower(), (
            f"SummaryError should mention 'parse error', got: {exc_info.value}"
        )


class TestSummaryFort14IOErrors:
    """Test fort14 I/O errors (line 161)."""

    def test_fort14_is_directory_error(self, tmp_path):
        """Opening a directory as fort14 raises SummaryError (IsADirectoryError)."""
        mesh_dir = tmp_path / "dir.14"
        mesh_dir.mkdir()

        with pytest.raises(SummaryError) as exc_info:
            summary(mesh_dir)

        # IsADirectoryError is caught as IOError and wrapped in SummaryError
        assert exc_info.value is not None, "Expected SummaryError to be raised"


class TestSummaryMeshAttributeError:
    """Test CHILmesh object missing attribute (lines 83-84)."""

    def test_mesh_object_missing_points_attribute(self, tmp_path):
        """CHILmesh object without 'points' attribute raises SummaryError."""
        # Create a mock object that looks like CHILmesh (type().__name__ == 'CHILmesh')
        # but has no 'points' attribute
        class CHILmesh:
            """Mock CHILmesh with missing points attribute."""
            pass

        fake_mesh = CHILmesh()

        with pytest.raises(SummaryError) as exc_info:
            summary(fake_mesh)

        assert "failed to extract metadata" in str(exc_info.value).lower(), (
            f"SummaryError should mention 'failed to extract metadata', got: {exc_info.value}"
        )


class TestSummaryDeepLoadError:
    """Test deep-load failure (lines 133-134)."""

    def test_fort14_deep_load_with_invalid_body(self, tmp_path):
        """Fort14 with valid header but invalid body fails on deep=True."""
        bad_fort14 = tmp_path / "bad_body.fort.14"
        # Valid header, invalid body (not enough data to parse elements)
        bad_fort14.write_text("name\n1 1\n")

        with pytest.raises(SummaryError) as exc_info:
            summary(bad_fort14, deep=True)

        assert "deep summary" in str(exc_info.value).lower(), (
            f"SummaryError should mention 'deep summary', got: {exc_info.value}"
        )


class TestSummary2dmHeaderParsing:
    """Test 2dm header parsing (lines 175, 179-181, 194)."""

    def test_2dm_with_meshname_and_elements(self, tmp_path):
        """2dm file with MESHNAME, ND, and E3T records parses correctly."""
        mesh_2dm = tmp_path / "test_mesh.2dm"
        content = """MESHNAME my_test_mesh

ND 1 0.0 0.0 0.0
ND 2 1.0 0.0 0.0
ND 3 0.0 1.0 0.0
E3T 1 1 2 3 1
"""
        mesh_2dm.write_text(content)

        result = summary(mesh_2dm)

        assert result['format'] == '2dm', (
            f"Expected format='2dm', got {result['format']}"
        )
        assert result['n_nodes'] == 3, f"Expected n_nodes=3, got {result['n_nodes']}"
        assert result['n_elems'] == 1, f"Expected n_elems=1, got {result['n_elems']}"
        assert result['grid_name'] == 'my_test_mesh', (
            f"Expected grid_name='my_test_mesh', got {result['grid_name']}"
        )


class TestSummary2dmIOErrors:
    """Test 2dm I/O errors (lines 195-196)."""

    def test_2dm_is_directory_error(self, tmp_path):
        """Opening a directory as 2dm raises SummaryError (IsADirectoryError)."""
        mesh_dir = tmp_path / "dir.2dm"
        mesh_dir.mkdir()

        with pytest.raises(SummaryError) as exc_info:
            summary(mesh_dir)

        # IsADirectoryError is caught as IOError and wrapped in SummaryError
        assert exc_info.value is not None, "Expected SummaryError to be raised"
