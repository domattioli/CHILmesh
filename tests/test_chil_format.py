"""Test .chil binary format I/O."""

import numpy as np
import pytest
import tempfile
from pathlib import Path
from chilmesh._chil_format import write_chil, read_chil


def test_write_read_roundtrip():
    """Write mesh to .chil, read back, verify equality."""
    points = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [1, 1, 0],
    ], dtype=np.float64)
    connectivity = np.array([
        [0, 1, 2],
        [1, 3, 2],
    ], dtype=np.int32)
    grid_name = "test_quad"

    with tempfile.TemporaryDirectory() as tmpdir:
        filepath = Path(tmpdir) / "test.chil"

        # Write
        ok = write_chil(filepath, points, connectivity, grid_name)
        assert ok

        # Read back
        pts, conn, name = read_chil(filepath)

        assert np.allclose(pts, points)
        assert np.array_equal(conn, connectivity)
        assert name == grid_name


def test_quad_elements():
    """Test roundtrip with quad elements."""
    points = np.eye(5)  # 5 vertices
    # All quads: 2 four-node quads
    connectivity = np.array([
        [0, 1, 2, 3],
        [1, 2, 3, 4],
    ], dtype=np.int32)

    with tempfile.TemporaryDirectory() as tmpdir:
        filepath = Path(tmpdir) / "quads.chil"
        write_chil(filepath, points, connectivity)
        pts, conn, _ = read_chil(filepath)

        assert conn.shape[0] == 2
        assert conn.shape[1] == 4  # all quads


def test_read_nonexistent_file():
    """Attempting to read missing file raises ValueError."""
    with pytest.raises(ValueError):
        read_chil(Path("/nonexistent/file.chil"))


def test_write_creates_file():
    """write_chil creates a file on disk (as .chil.npz)."""
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=np.float64)
    connectivity = np.array([[0, 1, 2]], dtype=np.int32)  # triangle

    with tempfile.TemporaryDirectory() as tmpdir:
        filepath = Path(tmpdir) / "test.chil"
        assert not filepath.exists()

        ok = write_chil(filepath, points, connectivity)
        assert ok
        # savez_compressed appends .npz to the filename
        assert (filepath.parent / "test.chil.npz").exists()
