"""End-to-end tests for the ``chilmesh.summary`` API.

Covers lazy header-only reading (fort.14, 2dm) and deep full-mesh loading.
Tests run against bundled fixtures from ``chilmesh.examples.fixture_path()``.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from chilmesh import CHILmesh, examples
from chilmesh.summary_io import summary, SummaryError


# ============================================================================
# Unit Tests: Direct API
# ============================================================================


@pytest.fixture
def annulus_path():
    """Path to bundled annulus fixture."""
    return examples.fixture_path("annulus_200pts.fort.14")


@pytest.fixture
def donut_path():
    """Path to bundled donut fixture."""
    return examples.fixture_path("donut_domain.fort.14")


def test_summary_lazy_fort14_annulus(annulus_path):
    """Load annulus_200pts.fort.14 lazily; assert key metadata present."""
    result = summary(annulus_path)

    # Basic metadata always present
    assert result['format'] == 'fort14'
    assert result['n_nodes'] == 380
    assert result['n_elems'] == 580
    assert result['file_bytes'] > 0
    assert isinstance(result['path'], str)

    # Lazy read: no exception
    assert 'grid_name' in result


def test_summary_lazy_counts_match_mesh():
    """Lazy n_nodes/n_elems match authoritative CHILmesh.load() values.

    Parametrized over annulus and donut.
    """
    fixtures = [
        ('annulus_200pts.fort.14', 380, 580),
        ('donut_domain.fort.14', 188, 276),
    ]

    for fname, expected_nodes, expected_elems in fixtures:
        path = examples.fixture_path(fname)

        # Lazy summary
        lazy = summary(path)
        assert lazy['n_nodes'] == expected_nodes, \
            f"{fname}: lazy n_nodes mismatch"
        assert lazy['n_elems'] == expected_elems, \
            f"{fname}: lazy n_elems mismatch"

        # Full load to verify against authoritative source
        mesh = CHILmesh.read_from_fort14(path)
        assert mesh.n_verts == expected_nodes
        assert mesh.n_elems == expected_elems


def test_summary_deep_adds_bbox_and_type(annulus_path):
    """deep=True adds element_type and bbox dict with valid float bounds."""
    result = summary(annulus_path, deep=True)

    # Lazy keys still present
    assert result['format'] == 'fort14'
    assert result['n_nodes'] == 380
    assert result['n_elems'] == 580

    # Deep-only keys added
    assert 'element_type' in result
    assert result['element_type'] == 'Triangular'
    assert isinstance(result['element_type'], str)
    assert len(result['element_type']) > 0

    # Bbox validation
    assert 'bbox' in result
    bbox = result['bbox']
    assert isinstance(bbox, dict)
    assert set(bbox.keys()) == {'min_x', 'min_y', 'max_x', 'max_y'}

    # All bounds are floats
    for key in bbox:
        assert isinstance(bbox[key], float)

    # min <= max invariant
    assert bbox['min_x'] <= bbox['max_x']
    assert bbox['min_y'] <= bbox['max_y']


def test_summary_from_mesh_object(annulus_path):
    """Call summary() on loaded CHILmesh object; confirm all fields match object."""
    # Load the mesh
    mesh = CHILmesh.read_from_fort14(annulus_path)

    # Get summary from the object
    result = summary(mesh)

    # Mesh-object format
    assert result['format'] == 'mesh-object'
    assert result['path'] is None

    # Counts match object
    assert result['n_nodes'] == mesh.n_verts
    assert result['n_elems'] == mesh.n_elems

    # Element type matches
    assert result['element_type'] == mesh.type
    assert result['element_type'] == 'Triangular'

    # Bbox is present and valid
    assert 'bbox' in result
    bbox = result['bbox']

    # Bounds match actual mesh points
    points = mesh.points
    assert bbox['min_x'] == float(points[:, 0].min())
    assert bbox['max_x'] == float(points[:, 0].max())
    assert bbox['min_y'] == float(points[:, 1].min())
    assert bbox['max_y'] == float(points[:, 1].max())


def test_summary_unknown_suffix_raises(tmp_path):
    """Unsupported file suffix (.xyz) raises SummaryError."""
    bad_file = tmp_path / "x.xyz"
    bad_file.write_text("not a real mesh")

    with pytest.raises(SummaryError, match="Unknown mesh format"):
        summary(bad_file)


# ============================================================================
# CLI Integration Tests
# ============================================================================


def _run(args, **kw):
    """Invoke the CLI via ``python -m chilmesh`` for path-independence."""
    cmd = [sys.executable, "-m", "chilmesh", *args]
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def test_cli_summary_subcommand(annulus_path):
    """Subprocess `python -m chilmesh summary <annulus>` exits 0."""
    result = _run(["summary", str(annulus_path)])

    assert result.returncode == 0, result.stderr

    # Output contains expected metadata
    assert "Nodes:" in result.stdout
    assert "580" in result.stdout
    assert "Elements:" in result.stdout


def test_cli_summary_deep_flag(annulus_path):
    """With --deep flag, summary adds bbox to stdout."""
    result = _run(["summary", str(annulus_path), "--deep"])

    assert result.returncode == 0, result.stderr

    # Lazy output present
    assert "Nodes:" in result.stdout
    assert "580" in result.stdout

    # Deep output present
    assert "Bbox:" in result.stdout or "bbox" in result.stdout.lower()
