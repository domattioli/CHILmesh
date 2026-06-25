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
