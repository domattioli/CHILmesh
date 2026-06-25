"""Tests for Gmsh .msh format I/O (Issue #184)."""
from pathlib import Path

import pytest
import numpy as np

import chilmesh


def _write_msh_v2_2(tmp_path: Path, content: str, name: str = "mesh.msh") -> Path:
    """Helper to write v2.2 .msh content to a temp file."""
    path = tmp_path / name
    path.write_text(content)
    return path


def _write_msh_v4_1(tmp_path: Path, content: str, name: str = "mesh.msh") -> Path:
    """Helper to write v4.1 .msh content to a temp file."""
    path = tmp_path / name
    path.write_text(content)
    return path


class TestGmshIOBasic:
    """Basic Gmsh I/O tests."""

    def test_read_msh_v2_2_triangular(self, tmp_path):
        """Test reading a pure triangular mesh from Gmsh v2.2 format."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
4 0.0 1.0 0.0
$EndNodes
$Elements
2
1 2 2 0 0 1 2 3
2 2 2 0 0 1 3 4
$EndElements
"""
        path = _write_msh_v2_2(tmp_path, content)
        mesh = chilmesh.read_msh(str(path))

        assert mesh.n_verts == 4
        assert mesh.n_elems == 2
        assert mesh.connectivity_list.shape == (2, 3)
        assert np.allclose(mesh.points[0], [0.0, 0.0, 0.0])
        assert np.allclose(mesh.points[2], [1.0, 1.0, 0.0])

    def test_read_msh_v4_1_triangular(self, tmp_path):
        """Test reading a pure triangular mesh from Gmsh v4.1 format."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
0 1 0 4
1
2
3
4
0.0 0.0 0.0
1.0 0.0 0.0
1.0 1.0 0.0
0.0 1.0 0.0
$EndNodes
$Elements
1 2 1 2
2 1 2 2
1 1 2 3
2 1 3 4
$EndElements
"""
        path = _write_msh_v4_1(tmp_path, content)
        mesh = chilmesh.read_msh(str(path))

        assert mesh.n_verts == 4
        assert mesh.n_elems == 2
        assert mesh.connectivity_list.shape == (2, 3)

    def test_read_msh_v2_2_quad(self, tmp_path):
        """Test reading a pure quad mesh from Gmsh v2.2 format."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
4 0.0 1.0 0.0
$EndNodes
$Elements
1
1 3 2 0 0 1 2 3 4
$EndElements
"""
        path = _write_msh_v2_2(tmp_path, content)
        mesh = chilmesh.read_msh(str(path))

        assert mesh.n_verts == 4
        assert mesh.n_elems == 1
        assert mesh.connectivity_list.shape == (1, 4)

    def test_read_msh_v4_1_quad(self, tmp_path):
        """Test reading a pure quad mesh from Gmsh v4.1 format."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
0 1 0 4
1
2
3
4
0.0 0.0 0.0
1.0 0.0 0.0
1.0 1.0 0.0
0.0 1.0 0.0
$EndNodes
$Elements
1 1 1 1
2 3 3 1
1 1 2 3 4
$EndElements
"""
        path = _write_msh_v4_1(tmp_path, content)
        mesh = chilmesh.read_msh(str(path))

        assert mesh.n_verts == 4
        assert mesh.n_elems == 1
        assert mesh.connectivity_list.shape == (1, 4)

    def test_read_msh_v2_2_mixed(self, tmp_path):
        """Test reading a mixed tri/quad mesh from Gmsh v2.2 format."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
6
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 2.0 0.0 0.0
4 0.0 1.0 0.0
5 1.0 1.0 0.0
6 2.0 1.0 0.0
$EndNodes
$Elements
3
1 2 2 0 0 1 2 4
2 3 2 0 0 2 3 5 6
3 3 2 0 0 3 6 4 1
$EndElements
"""
        path = _write_msh_v2_2(tmp_path, content)
        mesh = chilmesh.read_msh(str(path))

        assert mesh.n_verts == 6
        assert mesh.n_elems == 3
        # Mixed mesh should be 4-column (padded)
        assert mesh.connectivity_list.shape == (3, 4)


class TestGmshRoundtrip:
    """Round-trip read/write tests."""

    def test_roundtrip_v2_2_triangular(self, tmp_path):
        """Test round-trip v2.2 triangular mesh."""
        # Create a small mesh
        points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        conn = np.array([[0, 1, 2], [0, 2, 3]])

        # Write
        path1 = tmp_path / "mesh1.msh"
        assert chilmesh.write_msh(str(path1), points, conn, version="2.2")

        # Read back
        mesh = chilmesh.read_msh(str(path1))

        assert mesh.n_verts == 4
        assert mesh.n_elems == 2
        assert np.allclose(mesh.points, points)
        assert np.array_equal(mesh.connectivity_list, conn)

    def test_roundtrip_v4_1_triangular(self, tmp_path):
        """Test round-trip v4.1 triangular mesh."""
        points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        conn = np.array([[0, 1, 2], [0, 2, 3]])

        # Write
        path1 = tmp_path / "mesh1.msh"
        assert chilmesh.write_msh(str(path1), points, conn, version="4.1")

        # Read back
        mesh = chilmesh.read_msh(str(path1))

        assert mesh.n_verts == 4
        assert mesh.n_elems == 2
        assert np.allclose(mesh.points, points)
        assert np.array_equal(mesh.connectivity_list, conn)

    def test_roundtrip_v2_2_quad(self, tmp_path):
        """Test round-trip v2.2 quad mesh."""
        points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        conn = np.array([[0, 1, 2, 3]])

        path1 = tmp_path / "mesh1.msh"
        assert chilmesh.write_msh(str(path1), points, conn, version="2.2")

        mesh = chilmesh.read_msh(str(path1))

        assert mesh.n_verts == 4
        assert mesh.n_elems == 1
        assert np.allclose(mesh.points, points)
        assert np.array_equal(mesh.connectivity_list, conn)

    def test_roundtrip_v4_1_quad(self, tmp_path):
        """Test round-trip v4.1 quad mesh."""
        points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        conn = np.array([[0, 1, 2, 3]])

        path1 = tmp_path / "mesh1.msh"
        assert chilmesh.write_msh(str(path1), points, conn, version="4.1")

        mesh = chilmesh.read_msh(str(path1))

        assert mesh.n_verts == 4
        assert mesh.n_elems == 1
        assert np.allclose(mesh.points, points)
        assert np.array_equal(mesh.connectivity_list, conn)

    def test_roundtrip_v2_2_mixed(self, tmp_path):
        """Test round-trip v2.2 mixed mesh (padded triangles written as pure tris)."""
        points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        # Mixed: padded triangles. When written, they are pure triangles.
        # When read back, they stay as 3-column (pure triangles).
        conn = np.array([
            [0, 1, 2, 0],  # Padded triangle
            [1, 2, 3, 1],  # Padded triangle
        ])

        path1 = tmp_path / "mesh1.msh"
        assert chilmesh.write_msh(str(path1), points, conn, version="2.2")

        mesh = chilmesh.read_msh(str(path1))

        assert mesh.n_verts == 4
        assert mesh.n_elems == 2
        assert np.allclose(mesh.points, points)
        # After round-trip, padded triangles become pure triangles (3-column)
        expected_conn = np.array([[0, 1, 2], [1, 2, 3]])
        assert np.array_equal(mesh.connectivity_list, expected_conn)

    def test_roundtrip_v4_1_mixed(self, tmp_path):
        """Test round-trip v4.1 mixed mesh (padded triangles written as pure tris)."""
        points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        conn = np.array([
            [0, 1, 2, 0],  # Padded triangle
            [1, 2, 3, 1],  # Padded triangle
        ])

        path1 = tmp_path / "mesh1.msh"
        assert chilmesh.write_msh(str(path1), points, conn, version="4.1")

        mesh = chilmesh.read_msh(str(path1))

        assert mesh.n_verts == 4
        assert mesh.n_elems == 2
        assert np.allclose(mesh.points, points)
        # After round-trip, padded triangles become pure triangles (3-column)
        expected_conn = np.array([[0, 1, 2], [1, 2, 3]])
        assert np.array_equal(mesh.connectivity_list, expected_conn)


class TestGmshClassMethods:
    """Test CHILmesh class methods for Gmsh I/O."""

    def test_mesh_read_from_msh_static(self, tmp_path):
        """Test CHILmesh.read_from_msh() static method."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
$EndNodes
$Elements
1
1 2 2 0 0 1 2 3
$EndElements
"""
        path = _write_msh_v2_2(tmp_path, content)
        mesh = chilmesh.CHILmesh.read_from_msh(str(path))

        assert mesh.n_verts == 3
        assert mesh.n_elems == 1
        assert mesh.connectivity_list.shape == (1, 3)

    def test_mesh_write_to_msh_instance(self, tmp_path):
        """Test CHILmesh.write_to_msh() instance method."""
        points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
        conn = np.array([[0, 1, 2]])

        mesh = chilmesh.CHILmesh(conn, points, grid_name="Test", compute_layers=False, compute_adjacencies=False)
        path = tmp_path / "out.msh"

        assert mesh.write_to_msh(str(path), grid_name="MyMesh", version="2.2")
        assert path.exists()

        # Verify by reading back
        mesh2 = chilmesh.read_msh(str(path))
        assert mesh2.n_verts == 3
        assert mesh2.n_elems == 1

    def test_mesh_compute_layers_parameter(self, tmp_path):
        """Test that compute_layers parameter is honored."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
$EndNodes
$Elements
1
1 2 2 0 0 1 2 3
$EndElements
"""
        path = _write_msh_v2_2(tmp_path, content)

        # Fast init (compute_layers=False)
        mesh_fast = chilmesh.CHILmesh.read_from_msh(str(path), compute_layers=False)
        assert mesh_fast.n_layers == 0
        assert mesh_fast.adjacencies == {}

        # Normal init (compute_layers=True)
        mesh_layers = chilmesh.CHILmesh.read_from_msh(str(path), compute_layers=True)
        assert mesh_layers.n_layers > 0


class TestGmshErrors:
    """Test error handling."""

    def test_missing_meshformat(self, tmp_path):
        """Test error when $MeshFormat section missing."""
        content = "$Nodes\n3\n1 0 0 0\n"
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="Missing.*MeshFormat"):
            chilmesh.read_msh(str(path))

    def test_missing_nodes(self, tmp_path):
        """Test error when $Nodes section missing."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Elements
1
1 2 2 0 0 1 2 3
$EndElements
"""
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="Missing.*Nodes"):
            chilmesh.read_msh(str(path))

    def test_missing_elements(self, tmp_path):
        """Test error when $Elements section missing."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
$EndNodes
"""
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="Missing.*Elements"):
            chilmesh.read_msh(str(path))

    def test_unsupported_version(self, tmp_path):
        """Test error on unsupported Gmsh version."""
        content = """$MeshFormat
3.0 0 8
$EndMeshFormat
$Nodes
0
$EndNodes
$Elements
0
$EndElements
"""
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="Unsupported.*format"):
            chilmesh.read_msh(str(path))

    def test_invalid_write_version(self, tmp_path):
        """Test error on invalid version in write_msh."""
        points = np.array([[0.0, 0.0, 0.0]])
        conn = np.array([[0, 0, 0]])
        path = tmp_path / "out.msh"

        with pytest.raises(ValueError, match="Unsupported.*version"):
            chilmesh.write_msh(str(path), points, conn, version="3.0")


class TestGmshNodeOrdering:
    """Test that node ordering is preserved."""

    def test_v2_2_non_contiguous_node_ids(self, tmp_path):
        """Test v2.2 with non-contiguous node IDs."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
10 0.0 0.0 0.0
20 1.0 0.0 0.0
30 1.0 1.0 0.0
$EndNodes
$Elements
1
1 2 2 0 0 10 20 30
$EndElements
"""
        path = _write_msh_v2_2(tmp_path, content)
        mesh = chilmesh.read_msh(str(path))

        # Should map nodes in sorted order: 10→0, 20→1, 30→2
        assert mesh.n_verts == 3
        assert np.allclose(mesh.points[0], [0.0, 0.0, 0.0])  # Node 10
        assert np.allclose(mesh.points[1], [1.0, 0.0, 0.0])  # Node 20
        assert np.allclose(mesh.points[2], [1.0, 1.0, 0.0])  # Node 30

    def test_v4_1_non_contiguous_node_ids(self, tmp_path):
        """Test v4.1 with non-contiguous node IDs."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 3 10 30
0 1 0 3
10
20
30
0.0 0.0 0.0
1.0 0.0 0.0
1.0 1.0 0.0
$EndNodes
$Elements
1 1 1 1
2 1 2 1
1 10 20 30
$EndElements
"""
        path = _write_msh_v4_1(tmp_path, content)
        mesh = chilmesh.read_msh(str(path))

        assert mesh.n_verts == 3
        # Nodes sorted: 10, 20, 30 → indices 0, 1, 2
        assert np.allclose(mesh.points[0], [0.0, 0.0, 0.0])
        assert np.allclose(mesh.points[1], [1.0, 0.0, 0.0])
        assert np.allclose(mesh.points[2], [1.0, 1.0, 0.0])


class TestGmshCoordinates:
    """Test coordinate handling."""

    def test_v2_2_2d_coords(self, tmp_path):
        """Test v2.2 with implicit z=0 for 2D."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
1 0.5 0.5 0.0
2 1.5 0.5 0.0
3 1.5 1.5 0.0
$EndNodes
$Elements
1
1 2 2 0 0 1 2 3
$EndElements
"""
        path = _write_msh_v2_2(tmp_path, content)
        mesh = chilmesh.read_msh(str(path))

        assert np.allclose(mesh.points[:, 2], 0.0)
        assert np.allclose(mesh.points[0, :2], [0.5, 0.5])

    def test_v4_1_2d_coords(self, tmp_path):
        """Test v4.1 with explicit z values."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 3 1 3
0 1 0 3
1
2
3
0.5 0.5 -1.0
1.5 0.5 -0.5
1.5 1.5 0.0
$EndNodes
$Elements
1 1 1 1
2 1 2 1
1 1 2 3
$EndElements
"""
        path = _write_msh_v4_1(tmp_path, content)
        mesh = chilmesh.read_msh(str(path))

        assert mesh.points[0, 2] == pytest.approx(-1.0)
        assert mesh.points[1, 2] == pytest.approx(-0.5)
        assert mesh.points[2, 2] == pytest.approx(0.0)


class TestGmshErrorBranches:
    """Test error handling for uncovered branches in gmsh_io.py."""

    # v2.2 Reader Error Tests
    def test_v2_2_meshformat_last_line(self, tmp_path):
        """Test v2.2 error: $MeshFormat is last line (nothing after)."""
        content = "$MeshFormat\n"
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="incomplete"):
            chilmesh.read_msh(str(path))

    def test_v2_2_meshformat_blank_version_line(self, tmp_path):
        """Test v2.2 error: $MeshFormat followed by blank/empty version line."""
        content = "$MeshFormat\n\n$Nodes\n3\n1 0 0 0\n"
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="no version line"):
            chilmesh.read_msh(str(path))

    def test_v2_2_nodes_last_line(self, tmp_path):
        """Test v2.2 error: $Nodes is last line (incomplete)."""
        # $Elements must come first so section-scan succeeds; then $Nodes is last with no header.
        content = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Elements\n0\n$EndElements\n$Nodes"
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="incomplete"):
            chilmesh.read_msh(str(path))

    def test_v2_2_nodes_count_not_integer(self, tmp_path):
        """Test v2.2 error: $Nodes count line is not an integer."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
abc
$EndNodes
$Elements
0
"""
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="not an integer"):
            chilmesh.read_msh(str(path))

    def test_v2_2_nodes_fewer_than_declared(self, tmp_path):
        """Test v2.2 error: $Nodes count says N but fewer node lines present."""
        # $Elements comes first; then $Nodes with count 3 but only 1 node line, then file ends mid-parse.
        content = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Elements\n0\n$EndElements\n$Nodes\n3\n1 0.0 0.0 0.0"
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="expected 3 nodes"):
            chilmesh.read_msh(str(path))

    def test_v2_2_node_line_too_few_fields(self, tmp_path):
        """Test v2.2 error: node line has <4 whitespace fields."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
1
1 0.0 0.0
$EndNodes
$Elements
0
"""
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="Node line malformed"):
            chilmesh.read_msh(str(path))

    def test_v2_2_node_line_non_numeric_coords(self, tmp_path):
        """Test v2.2 error: node line 4 fields but non-numeric coords."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
1
1 x y z
$EndNodes
$Elements
0
"""
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="non-numeric"):
            chilmesh.read_msh(str(path))

    def test_v2_2_no_nodes_found(self, tmp_path):
        """Test v2.2 error: nodes count is 0 / no nodes parsed."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
0
$EndNodes
$Elements
1
1 2 2 0 0 1 2 3
$EndElements
"""
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="No nodes found"):
            chilmesh.read_msh(str(path))

    def test_v2_2_elements_last_line(self, tmp_path):
        """Test v2.2 error: $Elements is last line (incomplete)."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
1
1 0.0 0.0 0.0
$EndNodes
$Elements
"""
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="incomplete"):
            chilmesh.read_msh(str(path))

    def test_v2_2_elements_count_not_integer(self, tmp_path):
        """Test v2.2 error: $Elements count not an integer."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
1
1 0.0 0.0 0.0
$EndNodes
$Elements
xyz
1 2 2 0 0 1 1 1
"""
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="not an integer"):
            chilmesh.read_msh(str(path))

    def test_v2_2_elements_fewer_than_declared(self, tmp_path):
        """Test v2.2 error: $Elements count says N but fewer element lines."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
$EndNodes
$Elements
2
1 2 2 0 0 1 2 3
"""
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="expected 2 elements"):
            chilmesh.read_msh(str(path))

    def test_v2_2_element_line_too_few_fields(self, tmp_path):
        """Test v2.2 error: element line has <4 fields."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
$EndNodes
$Elements
1
1 2 2
"""
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="Element line malformed"):
            chilmesh.read_msh(str(path))

    def test_v2_2_triangle_element_line_too_short(self, tmp_path):
        """Test v2.2 error: triangle (type 2) element line too short for nodes."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
$EndNodes
$Elements
1
1 2 2 0 0 1 2
"""
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="too short"):
            chilmesh.read_msh(str(path))

    def test_v2_2_element_line_non_numeric(self, tmp_path):
        """Test v2.2 error: element line with non-numeric fields."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
$EndNodes
$Elements
1
1 abc 2 0 0 1 2 3
"""
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="malformed or non-numeric"):
            chilmesh.read_msh(str(path))

    def test_v2_2_no_supported_elements(self, tmp_path):
        """Test v2.2 error: only unsupported element types (e.g. lines/points)."""
        content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
2
1 0.0 0.0 0.0
2 1.0 0.0 0.0
$EndNodes
$Elements
1
1 1 2 0 0 1 2
"""
        path = _write_msh_v2_2(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="No supported elements"):
            chilmesh.read_msh(str(path))

    # v4.1 Reader Error Tests
    def test_v4_1_nodes_last_line(self, tmp_path):
        """Test v4.1 error: $Nodes is last line (incomplete)."""
        # $Elements first, then $Nodes as final line with no header.
        content = "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Elements\n0 0 1 0\n$EndElements\n$Nodes"
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="incomplete"):
            chilmesh.read_msh(str(path))

    def test_v4_1_nodes_header_too_few_fields(self, tmp_path):
        """Test v4.1 error: $Nodes header has <4 fields."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 2 3
$EndNodes
$Elements
0 0 0 0
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="header malformed"):
            chilmesh.read_msh(str(path))

    def test_v4_1_nodes_header_non_numeric(self, tmp_path):
        """Test v4.1 error: $Nodes header non-numeric."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
a b c d
$EndNodes
$Elements
0 0 0 0
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="non-numeric"):
            chilmesh.read_msh(str(path))

    def test_v4_1_node_block_truncated_before_header(self, tmp_path):
        """Test v4.1 error: node block truncated before its header."""
        # $Elements first; then $Nodes header saying 1 block exists, but file ends before block header.
        content = "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Elements\n0 0 1 0\n$EndElements\n$Nodes\n1 1 1 1"
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="incomplete"):
            chilmesh.read_msh(str(path))

    def test_v4_1_node_block_header_too_few_fields(self, tmp_path):
        """Test v4.1 error: node block header <4 fields."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 1 1 1
0 1 2
$EndNodes
$Elements
0 0 0 0
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="header malformed"):
            chilmesh.read_msh(str(path))

    def test_v4_1_node_block_header_non_numeric(self, tmp_path):
        """Test v4.1 error: node block header non-numeric."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 1 1 1
a b c d
$EndNodes
$Elements
0 0 0 0
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="non-numeric"):
            chilmesh.read_msh(str(path))

    def test_v4_1_node_tags_truncated(self, tmp_path):
        """Test v4.1 error: node-tags truncated before all read."""
        # $Elements first; block header says 2 tags but only 1 tag line present, file ends.
        content = "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Elements\n0 0 1 0\n$EndElements\n$Nodes\n1 2 1 2\n0 1 0 2\n1"
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="node-tags incomplete"):
            chilmesh.read_msh(str(path))

    def test_v4_1_node_tag_line_non_integer(self, tmp_path):
        """Test v4.1 error: node tag line non-integer."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 1 1 1
0 1 0 1
abc
$EndNodes
$Elements
0 0 0 0
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="Node tag line malformed"):
            chilmesh.read_msh(str(path))

    def test_v4_1_coordinates_truncated(self, tmp_path):
        """Test v4.1 error: coordinates truncated."""
        # $Elements first; $Nodes header + block header + 1 tag line but no coordinate line, file ends.
        content = "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Elements\n0 0 1 0\n$EndElements\n$Nodes\n1 1 1 1\n0 1 0 1\n1"
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="coordinates incomplete"):
            chilmesh.read_msh(str(path))

    def test_v4_1_coordinate_line_too_few_fields(self, tmp_path):
        """Test v4.1 error: coordinate line <3 fields."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 1 1 1
0 1 0 1
1
0.0 0.0
$EndNodes
$Elements
0 0 0 0
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="coordinate line malformed"):
            chilmesh.read_msh(str(path))

    def test_v4_1_coordinate_line_non_numeric(self, tmp_path):
        """Test v4.1 error: coordinate line non-numeric."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 1 1 1
0 1 0 1
1
x y z
$EndNodes
$Elements
0 0 0 0
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="non-numeric"):
            chilmesh.read_msh(str(path))

    def test_v4_1_no_nodes_found(self, tmp_path):
        """Test v4.1 error: no nodes parsed."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
0 0 0 0
$EndNodes
$Elements
1 1 1 1
2 1 2 1
1 1 2 3
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="No nodes found"):
            chilmesh.read_msh(str(path))

    def test_v4_1_elements_last_line(self, tmp_path):
        """Test v4.1 error: $Elements is last line (incomplete)."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 1 1 1
0 1 0 1
1
0.0 0.0 0.0
$EndNodes
$Elements
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="incomplete"):
            chilmesh.read_msh(str(path))

    def test_v4_1_elements_header_too_few_fields(self, tmp_path):
        """Test v4.1 error: $Elements header <4 fields."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 1 1 1
0 1 0 1
1
0.0 0.0 0.0
$EndNodes
$Elements
1 2 3
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="header malformed"):
            chilmesh.read_msh(str(path))

    def test_v4_1_elements_header_non_numeric(self, tmp_path):
        """Test v4.1 error: $Elements header non-numeric."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 1 1 1
0 1 0 1
1
0.0 0.0 0.0
$EndNodes
$Elements
a b c d
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="non-numeric"):
            chilmesh.read_msh(str(path))

    def test_v4_1_element_block_truncated_before_header(self, tmp_path):
        """Test v4.1 error: element block truncated before header."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 1 1 1
0 1 0 1
1
0.0 0.0 0.0
$EndNodes
$Elements
1 1 1 1
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="incomplete"):
            chilmesh.read_msh(str(path))

    def test_v4_1_element_block_header_too_few_fields(self, tmp_path):
        """Test v4.1 error: element block header <4 fields."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 1 1 1
0 1 0 1
1
0.0 0.0 0.0
$EndNodes
$Elements
1 1 1 1
2 1 2
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="header malformed"):
            chilmesh.read_msh(str(path))

    def test_v4_1_element_block_header_non_numeric(self, tmp_path):
        """Test v4.1 error: element block header non-numeric."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 1 1 1
0 1 0 1
1
0.0 0.0 0.0
$EndNodes
$Elements
1 1 1 1
a b c d
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="non-numeric"):
            chilmesh.read_msh(str(path))

    def test_v4_1_element_block_truncated_mid_block(self, tmp_path):
        """Test v4.1 error: element block truncated mid-block."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 1 1 1
0 1 0 1
1
0.0 0.0 0.0
$EndNodes
$Elements
1 2 1 1
2 1 2 2
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="incomplete"):
            chilmesh.read_msh(str(path))

    def test_v4_1_triangle_element_line_too_short(self, tmp_path):
        """Test v4.1 error: triangle (type 2) element line too short (<4 fields)."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 3 1 3
0 1 0 3
1
2
3
0.0 0.0 0.0
1.0 0.0 0.0
1.0 1.0 0.0
$EndNodes
$Elements
1 1 1 1
2 1 2 1
1 1 2
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="Triangle element line too short"):
            chilmesh.read_msh(str(path))

    def test_v4_1_quad_element_line_too_short(self, tmp_path):
        """Test v4.1 error: quad (type 3) element line too short (<5 fields)."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
0 1 0 4
1
2
3
4
0.0 0.0 0.0
1.0 0.0 0.0
1.0 1.0 0.0
0.0 1.0 0.0
$EndNodes
$Elements
1 1 1 1
2 2 3 1
1 1 2 3
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="Quad element line too short"):
            chilmesh.read_msh(str(path))

    def test_v4_1_element_line_non_numeric_tag(self, tmp_path):
        """Test v4.1 error: element line non-numeric tag."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 3 1 3
0 1 0 3
1
2
3
0.0 0.0 0.0
1.0 0.0 0.0
1.0 1.0 0.0
$EndNodes
$Elements
1 1 1 1
2 1 2 1
x 1 2 3
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="Element line malformed"):
            chilmesh.read_msh(str(path))

    def test_v4_1_no_supported_elements(self, tmp_path):
        """Test v4.1 error: only unsupported element types."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 2 1 2
0 1 0 2
1
2
0.0 0.0 0.0
1.0 0.0 0.0
$EndNodes
$Elements
1 1 1 1
2 1 1 1
1 1 2
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="No supported elements"):
            chilmesh.read_msh(str(path))

    # Writer Error Tests (unwritable path)
    def test_write_msh_v2_2_write_failure(self, tmp_path):
        """Test v2.2 write failure: unwritable path (nonexistent directory)."""
        points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
        conn = np.array([[0, 1, 2]])
        bad_path = str(tmp_path / "nonexistent_dir" / "out.msh")

        result = chilmesh.write_msh(bad_path, points, conn, version="2.2")
        assert result is False

    def test_write_msh_v4_1_write_failure(self, tmp_path):
        """Test v4.1 write failure: unwritable path (nonexistent directory)."""
        points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
        conn = np.array([[0, 1, 2]])
        bad_path = str(tmp_path / "nonexistent_dir" / "out.msh")

        result = chilmesh.write_msh(bad_path, points, conn, version="4.1")
        assert result is False

    def test_v4_1_missing_nodes_section(self, tmp_path):
        """Test v4.1 error: Missing $Nodes section."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Elements
0 0 1 0
$EndElements
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="Missing.*Nodes"):
            chilmesh.read_msh(str(path))

    def test_v4_1_missing_elements_section(self, tmp_path):
        """Test v4.1 error: Missing $Elements section."""
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 1 1 1
0 1 0 1
1
0.0 0.0 0.0
$EndNodes
"""
        path = _write_msh_v4_1(tmp_path, content, "bad.msh")

        with pytest.raises(chilmesh.GmshParseError, match="Missing.*Elements"):
            chilmesh.read_msh(str(path))
