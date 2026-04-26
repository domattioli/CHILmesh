"""Tests for SMS .2dm format reader (Issue #44)."""
import tempfile
from pathlib import Path

import numpy as np
import pytest

import chilmesh


class Test2DMReader:
    """Test CHILmesh.read_from_2dm() SMS format support."""

    def test_read_2dm_triangular_mesh(self):
        """Test reading a pure triangular mesh from .2dm format."""
        # Create a minimal triangular mesh
        sms_content = """MESH2D
ND 1 0.0 0.0 0.0
ND 2 1.0 0.0 0.0
ND 3 1.0 1.0 0.0
ND 4 0.0 1.0 0.0
E3T 1 1 2 3 0
E3T 2 1 3 4 0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.2dm', delete=False) as f:
            f.write(sms_content)
            tmpfile = Path(f.name)

        try:
            mesh = chilmesh.CHILmesh.read_from_2dm(tmpfile)

            assert mesh.n_verts == 4, "Should have 4 nodes"
            assert mesh.n_elems == 2, "Should have 2 triangles"
            assert mesh.type == "Triangular", "Should be Triangular type"
            assert mesh.connectivity_list.shape[1] == 3, "Should be 3-column array"
        finally:
            tmpfile.unlink()

    def test_read_2dm_quad_mesh(self):
        """Test reading a pure quad mesh from .2dm format."""
        # Create a minimal quad mesh
        sms_content = """MESH2D
ND 1 0.0 0.0 0.0
ND 2 1.0 0.0 0.0
ND 3 2.0 0.0 0.0
ND 4 0.0 1.0 0.0
ND 5 1.0 1.0 0.0
ND 6 2.0 1.0 0.0
E4Q 1 1 2 5 4 0
E4Q 2 2 3 6 5 0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.2dm', delete=False) as f:
            f.write(sms_content)
            tmpfile = Path(f.name)

        try:
            mesh = chilmesh.CHILmesh.read_from_2dm(tmpfile)

            assert mesh.n_verts == 6, "Should have 6 nodes"
            assert mesh.n_elems == 2, "Should have 2 quads"
            assert mesh.type == "Quadrilateral", "Should be Quadrilateral type"
            assert mesh.connectivity_list.shape[1] == 4, "Should be 4-column array"
        finally:
            tmpfile.unlink()

    def test_read_2dm_mixed_mesh(self):
        """Test reading a mixed-element mesh from .2dm format."""
        # Create a mesh with both triangles and quads
        sms_content = """MESH2D
ND 1 0.0 0.0 0.0
ND 2 1.0 0.0 0.0
ND 3 2.0 0.0 0.0
ND 4 0.0 1.0 0.0
ND 5 1.0 1.0 0.0
ND 6 2.0 1.0 0.0
E3T 1 1 2 4 0
E4Q 2 2 3 6 5 0
E3T 3 2 5 4 0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.2dm', delete=False) as f:
            f.write(sms_content)
            tmpfile = Path(f.name)

        try:
            mesh = chilmesh.CHILmesh.read_from_2dm(tmpfile)

            assert mesh.n_verts == 6, "Should have 6 nodes"
            assert mesh.n_elems == 3, "Should have 3 elements"
            assert mesh.type == "Mixed-Element", "Should be Mixed-Element type"
            assert mesh.connectivity_list.shape[1] == 4, "Mixed mesh should be 4-column (padded)"
        finally:
            tmpfile.unlink()

    def test_read_2dm_with_comments(self):
        """Test that .2dm reader skips comment lines."""
        sms_content = """MESH2D
# This is a comment
ND 1 0.0 0.0 0.0
# Another comment
ND 2 1.0 0.0 0.0
ND 3 1.0 1.0 0.0
E3T 1 1 2 3 0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.2dm', delete=False) as f:
            f.write(sms_content)
            tmpfile = Path(f.name)

        try:
            mesh = chilmesh.CHILmesh.read_from_2dm(tmpfile)
            assert mesh.n_verts == 3
            assert mesh.n_elems == 1
        finally:
            tmpfile.unlink()

    def test_read_2dm_with_z_coordinates(self):
        """Test that .2dm reader handles z-coordinates."""
        sms_content = """MESH2D
ND 1 0.0 0.0 -1.0
ND 2 1.0 0.0 -0.5
ND 3 1.0 1.0 0.0
E3T 1 1 2 3 0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.2dm', delete=False) as f:
            f.write(sms_content)
            tmpfile = Path(f.name)

        try:
            mesh = chilmesh.CHILmesh.read_from_2dm(tmpfile)
            assert mesh.points[0, 2] == pytest.approx(-1.0), "Z-coordinate should be preserved"
            assert mesh.points[1, 2] == pytest.approx(-0.5)
            assert mesh.points[2, 2] == pytest.approx(0.0)
        finally:
            tmpfile.unlink()

    def test_read_2dm_missing_file(self):
        """Test that missing file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            chilmesh.CHILmesh.read_from_2dm(Path("/nonexistent/file.2dm"))

    def test_read_2dm_empty_file(self):
        """Test that empty file raises ValueError."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.2dm', delete=False) as f:
            f.write("")
            tmpfile = Path(f.name)

        try:
            with pytest.raises(ValueError, match="No nodes found"):
                chilmesh.CHILmesh.read_from_2dm(tmpfile)
        finally:
            tmpfile.unlink()

    def test_read_2dm_compute_layers_false(self):
        """Test that compute_layers=False is respected by .2dm reader."""
        sms_content = """MESH2D
ND 1 0.0 0.0 0.0
ND 2 1.0 0.0 0.0
ND 3 1.0 1.0 0.0
ND 4 0.0 1.0 0.0
E3T 1 1 2 3 0
E3T 2 1 3 4 0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.2dm', delete=False) as f:
            f.write(sms_content)
            tmpfile = Path(f.name)

        try:
            mesh = chilmesh.CHILmesh.read_from_2dm(tmpfile, compute_layers=False)
            assert mesh.n_layers == 0
            assert mesh.adjacencies == {}
        finally:
            tmpfile.unlink()

    def test_read_2dm_from_admesh_domain(self):
        """Test .2dm reader via from_admesh_domain routing."""
        from types import SimpleNamespace

        sms_content = """MESH2D
ND 1 0.0 0.0 0.0
ND 2 1.0 0.0 0.0
ND 3 1.0 1.0 0.0
E3T 1 1 2 3 0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.2dm', delete=False) as f:
            f.write(sms_content)
            tmpfile = Path(f.name)

        try:
            record = SimpleNamespace(filename=str(tmpfile), type="SMS_2DM")
            mesh = chilmesh.CHILmesh.from_admesh_domain(record)

            assert mesh.n_verts == 3
            assert mesh.n_elems == 1
            assert mesh.type == "Triangular"
        finally:
            tmpfile.unlink()

    def test_read_2dm_preserves_geometry(self):
        """Test that .2dm geometry is preserved correctly."""
        sms_content = """MESH2D
ND 1 0.5 0.5 0.0
ND 2 1.5 0.5 0.0
ND 3 1.5 1.5 0.0
E3T 1 1 2 3 0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.2dm', delete=False) as f:
            f.write(sms_content)
            tmpfile = Path(f.name)

        try:
            mesh = chilmesh.CHILmesh.read_from_2dm(tmpfile)
            assert mesh.points[0, 0] == pytest.approx(0.5)
            assert mesh.points[1, 1] == pytest.approx(0.5)
            assert mesh.points[2, 0] == pytest.approx(1.5)
        finally:
            tmpfile.unlink()
