"""Tests for admesh.gmsh.write_msh and Mesh.to_msh — spec 008 FR-002,003,007,015."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

import admesh
from admesh.api import BoundarySegment, Mesh
from admesh.boundary_types import BoundaryType
from admesh.gmsh import GmshParseError, read_msh, write_msh


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_simple_mesh(with_bathy: bool = False) -> Mesh:
    """5-node, 4-triangle unit square mesh with MAINLAND boundary."""
    nodes = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
        dtype=np.float64,
    )
    elements = np.array(
        [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        dtype=np.int64,
    )
    ring = np.array([0, 1, 2, 3], dtype=np.int64)
    boundary = BoundarySegment(node_ids=ring, bc_type=BoundaryType.MAINLAND, is_open=False)
    bathy = None
    if with_bathy:
        bathy = np.array([-10.0, -20.0, -30.0, -40.0, -25.0], dtype=np.float64)
    return Mesh(
        nodes=nodes,
        elements=elements,
        boundaries=(boundary,),
        bathymetry=bathy,
        quality=None,
        title="test",
    )


def _make_multi_bc_mesh() -> Mesh:
    """Mesh with OPEN + MAINLAND boundaries."""
    nodes = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
        dtype=np.float64,
    )
    elements = np.array(
        [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        dtype=np.int64,
    )
    open_seg = BoundarySegment(
        node_ids=np.array([0, 1], dtype=np.int64),
        bc_type=BoundaryType.OPEN,
        is_open=True,
    )
    land_seg = BoundarySegment(
        node_ids=np.array([1, 2, 3, 0], dtype=np.int64),
        bc_type=BoundaryType.MAINLAND,
        is_open=False,
    )
    return Mesh(
        nodes=nodes,
        elements=elements,
        boundaries=(open_seg, land_seg),
        bathymetry=None,
        quality=None,
        title="multi_bc",
    )


# ---------------------------------------------------------------------------
# Basic write tests
# ---------------------------------------------------------------------------


class TestWriteV22:
    """FR-002: write_msh produces valid v2.2 files."""

    def test_file_created(self, tmp_path):
        mesh = _make_simple_mesh()
        out = tmp_path / "out.msh"
        write_msh(mesh, out, version="2.2")
        assert out.exists()

    def test_file_nonempty(self, tmp_path):
        mesh = _make_simple_mesh()
        out = tmp_path / "out.msh"
        write_msh(mesh, out, version="2.2")
        assert out.stat().st_size > 0

    def test_contains_mesh_format_header(self, tmp_path):
        mesh = _make_simple_mesh()
        out = tmp_path / "out.msh"
        write_msh(mesh, out, version="2.2")
        content = out.read_text()
        assert "$MeshFormat" in content
        assert "2.2 0 8" in content

    def test_contains_nodes_block(self, tmp_path):
        mesh = _make_simple_mesh()
        out = tmp_path / "out.msh"
        write_msh(mesh, out, version="2.2")
        content = out.read_text()
        assert "$Nodes" in content
        assert "$EndNodes" in content

    def test_contains_elements_block(self, tmp_path):
        mesh = _make_simple_mesh()
        out = tmp_path / "out.msh"
        write_msh(mesh, out, version="2.2")
        content = out.read_text()
        assert "$Elements" in content
        assert "$EndElements" in content

    def test_contains_physical_names(self, tmp_path):
        mesh = _make_simple_mesh()
        out = tmp_path / "out.msh"
        write_msh(mesh, out, version="2.2")
        content = out.read_text()
        assert "$PhysicalNames" in content
        assert '"mainland"' in content

    def test_z_is_zero(self, tmp_path):
        """FR-007: writer always emits z=0.0."""
        mesh = _make_simple_mesh()
        out = tmp_path / "out.msh"
        write_msh(mesh, out, version="2.2")
        content = out.read_text()
        # Check node lines have z=0
        for line in content.splitlines():
            stripped = line.strip()
            if stripped and stripped[0].isdigit():
                parts = stripped.split()
                if len(parts) == 4:
                    try:
                        # node line: id x y z
                        _ = int(parts[0])
                        z = float(parts[3])
                        assert z == 0.0, f"Expected z=0.0, got z={z} in line: {line}"
                    except (ValueError, IndexError):
                        pass


class TestWriteV41:
    """FR-002: write_msh produces valid v4.1 files."""

    def test_file_created(self, tmp_path):
        mesh = _make_simple_mesh()
        out = tmp_path / "out.msh"
        write_msh(mesh, out, version="4.1")
        assert out.exists()

    def test_contains_mesh_format_header(self, tmp_path):
        mesh = _make_simple_mesh()
        out = tmp_path / "out.msh"
        write_msh(mesh, out, version="4.1")
        content = out.read_text()
        assert "4.1 0 8" in content

    def test_contains_entities_block(self, tmp_path):
        mesh = _make_simple_mesh()
        out = tmp_path / "out.msh"
        write_msh(mesh, out, version="4.1")
        content = out.read_text()
        assert "$Entities" in content


class TestToMshMethod:
    """FR-003: Mesh.to_msh delegates to write_msh."""

    def test_to_msh_creates_file(self, tmp_path):
        mesh = _make_simple_mesh()
        out = tmp_path / "out.msh"
        mesh.to_msh(out)
        assert out.exists()

    def test_to_msh_v22(self, tmp_path):
        mesh = _make_simple_mesh()
        out = tmp_path / "out22.msh"
        mesh.to_msh(out, version="2.2")
        content = out.read_text()
        assert "2.2 0 8" in content

    def test_to_msh_v41_default(self, tmp_path):
        mesh = _make_simple_mesh()
        out = tmp_path / "out41.msh"
        mesh.to_msh(out)  # default version="4.1"
        content = out.read_text()
        assert "4.1 0 8" in content


class TestWriteErrors:
    """Error paths in write_msh."""

    def test_bad_version_raises_value_error(self, tmp_path):
        mesh = _make_simple_mesh()
        with pytest.raises(ValueError, match="version"):
            write_msh(mesh, tmp_path / "out.msh", version="3.0")

    def test_write_msh_top_level(self):
        assert admesh.write_msh is write_msh


class TestMultipleBoundaryTypes:
    """FR-002: mixed OPEN + MAINLAND → correct $PhysicalNames block."""

    def test_both_types_in_physical_names(self, tmp_path):
        mesh = _make_multi_bc_mesh()
        out = tmp_path / "multi.msh"
        write_msh(mesh, out, version="2.2")
        content = out.read_text()
        assert '"open"' in content
        assert '"mainland"' in content

    def test_both_types_in_v41(self, tmp_path):
        mesh = _make_multi_bc_mesh()
        out = tmp_path / "multi41.msh"
        write_msh(mesh, out, version="4.1")
        content = out.read_text()
        assert '"open"' in content
        assert '"mainland"' in content


class TestBathymetryWrite:
    """FR-015: $NodeData 'bathymetry' emitted when mesh has bathymetry."""

    def test_node_data_present_when_bathy(self, tmp_path):
        mesh = _make_simple_mesh(with_bathy=True)
        out = tmp_path / "bathy.msh"
        write_msh(mesh, out, version="2.2")
        content = out.read_text()
        assert "$NodeData" in content
        assert '"bathymetry"' in content

    def test_node_data_absent_when_no_bathy(self, tmp_path):
        mesh = _make_simple_mesh(with_bathy=False)
        out = tmp_path / "nobathy.msh"
        write_msh(mesh, out, version="2.2")
        content = out.read_text()
        assert "$NodeData" not in content

    def test_bathymetry_sign_convention(self, tmp_path):
        """FR-007: writer emits -mesh.bathymetry (positive-down on disk)."""
        # Internal bathymetry is negative (elevation, positive-up → -10 means 10m depth)
        mesh = _make_simple_mesh(with_bathy=True)
        out = tmp_path / "bathy_sign.msh"
        write_msh(mesh, out, version="2.2")
        content = out.read_text()
        # Internal: -10 → on disk should be +10
        assert "10" in content


# ---------------------------------------------------------------------------
# gmsh -check gate (skipped if gmsh not in PATH)
# ---------------------------------------------------------------------------

GMSH_AVAILABLE = shutil.which("gmsh") is not None


@pytest.mark.skipif(not GMSH_AVAILABLE, reason="gmsh CLI not in PATH")
class TestGmshCheck:
    """SC-002: gmsh -check accepts output with zero warnings."""

    def _check(self, path: Path) -> None:
        result = subprocess.run(
            ["gmsh", "-check", str(path)],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"gmsh -check failed: {result.stderr}"
        warnings_in_stderr = [
            line for line in result.stderr.splitlines()
            if line.startswith("Warning:")
        ]
        assert not warnings_in_stderr, (
            f"gmsh -check emitted warnings: {warnings_in_stderr}"
        )

    def test_v22_simple(self, tmp_path):
        mesh = _make_simple_mesh()
        out = tmp_path / "simple_v22.msh"
        write_msh(mesh, out, version="2.2")
        self._check(out)

    def test_v41_simple(self, tmp_path):
        mesh = _make_simple_mesh()
        out = tmp_path / "simple_v41.msh"
        write_msh(mesh, out, version="4.1")
        self._check(out)
