"""Tests for admesh.gmsh.read_msh — spec 008 FR-001,004,006,010-012,016."""

from __future__ import annotations

import os
import warnings
from pathlib import Path

import numpy as np
import pytest

import admesh
from admesh.gmsh import GmshParseError, read_msh, PHYSICAL_GROUP_MAP
from admesh.boundary_types import BoundaryType
from admesh.api import Mesh

FIXTURES = Path(__file__).parent / "fixtures" / "gmsh"


# ---------------------------------------------------------------------------
# Happy-path: basic read
# ---------------------------------------------------------------------------


class TestReadV22:
    """FR-001: read Gmsh ASCII v2.2."""

    def test_node_count(self):
        mesh = read_msh(FIXTURES / "triangle_v22.msh")
        assert mesh.n_nodes == 5

    def test_element_count(self):
        mesh = read_msh(FIXTURES / "triangle_v22.msh")
        assert mesh.n_elements == 4

    def test_nodes_shape_and_dtype(self):
        mesh = read_msh(FIXTURES / "triangle_v22.msh")
        assert mesh.nodes.shape == (5, 2)
        assert mesh.nodes.dtype == np.float64

    def test_elements_shape_and_dtype(self):
        mesh = read_msh(FIXTURES / "triangle_v22.msh")
        assert mesh.elements.shape == (4, 3)
        assert mesh.elements.dtype == np.int64

    def test_elements_zero_based(self):
        """FR-006: no 1-based indices leak into Mesh."""
        mesh = read_msh(FIXTURES / "triangle_v22.msh")
        assert mesh.elements.min() >= 0
        assert mesh.elements.max() < mesh.n_nodes

    def test_returns_mesh_instance(self):
        mesh = read_msh(FIXTURES / "triangle_v22.msh")
        assert isinstance(mesh, Mesh)

    def test_mainland_boundary_type(self):
        """PhysicalNames 'mainland' → BoundaryType.MAINLAND."""
        mesh = read_msh(FIXTURES / "triangle_v22.msh")
        assert len(mesh.boundaries) >= 1
        bc_types = {int(s.bc_type) for s in mesh.boundaries}
        assert int(BoundaryType.MAINLAND) in bc_types


class TestReadV41:
    """FR-001: read Gmsh ASCII v4.1."""

    def test_node_count(self):
        mesh = read_msh(FIXTURES / "triangle_v41.msh")
        assert mesh.n_nodes == 5

    def test_element_count(self):
        mesh = read_msh(FIXTURES / "triangle_v41.msh")
        assert mesh.n_elements == 4

    def test_elements_zero_based(self):
        mesh = read_msh(FIXTURES / "triangle_v41.msh")
        assert mesh.elements.min() >= 0
        assert mesh.elements.max() < mesh.n_nodes

    def test_returns_mesh_instance(self):
        mesh = read_msh(FIXTURES / "triangle_v41.msh")
        assert isinstance(mesh, Mesh)


class TestVersionParity:
    """Both v2.2 and v4.1 on the same logical mesh produce equal results."""

    def test_same_node_count(self):
        m22 = read_msh(FIXTURES / "triangle_v22.msh")
        m41 = read_msh(FIXTURES / "triangle_v41.msh")
        assert m22.n_nodes == m41.n_nodes

    def test_same_element_count(self):
        m22 = read_msh(FIXTURES / "triangle_v22.msh")
        m41 = read_msh(FIXTURES / "triangle_v41.msh")
        assert m22.n_elements == m41.n_elements

    def test_same_node_coordinates(self):
        m22 = read_msh(FIXTURES / "triangle_v22.msh")
        m41 = read_msh(FIXTURES / "triangle_v41.msh")
        # Node order may differ; compare sets of coordinates
        coords22 = set(map(tuple, np.round(m22.nodes, 10).tolist()))
        coords41 = set(map(tuple, np.round(m41.nodes, 10).tolist()))
        assert coords22 == coords41


# ---------------------------------------------------------------------------
# Boundary / physical groups
# ---------------------------------------------------------------------------


class TestPhysicalGroups:
    """FR-008,009: PhysicalNames → BoundaryType mapping."""

    def test_open_boundary_type(self):
        mesh = read_msh(FIXTURES / "open_boundary_v22.msh")
        bc_types = {int(s.bc_type) for s in mesh.boundaries}
        assert int(BoundaryType.OPEN) in bc_types

    def test_mainland_boundary_type(self):
        mesh = read_msh(FIXTURES / "open_boundary_v22.msh")
        bc_types = {int(s.bc_type) for s in mesh.boundaries}
        assert int(BoundaryType.MAINLAND) in bc_types

    def test_physical_group_map_constant(self):
        """PHYSICAL_GROUP_MAP covers the 4 canonical names."""
        assert "open" in PHYSICAL_GROUP_MAP
        assert "mainland" in PHYSICAL_GROUP_MAP
        assert "island" in PHYSICAL_GROUP_MAP
        assert "mainland_flux" in PHYSICAL_GROUP_MAP

    def test_physical_group_map_values(self):
        assert PHYSICAL_GROUP_MAP["open"] == BoundaryType.OPEN
        assert PHYSICAL_GROUP_MAP["mainland"] == BoundaryType.MAINLAND
        assert PHYSICAL_GROUP_MAP["island"] == BoundaryType.ISLAND
        assert PHYSICAL_GROUP_MAP["mainland_flux"] == BoundaryType.MAINLAND_FLUX


class TestNoPhysicalGroups:
    """FR-012: missing $PhysicalNames → UserWarning + single MAINLAND ring."""

    def test_warns_on_no_physical_names(self):
        with pytest.warns(UserWarning, match="no \\$PhysicalNames"):
            mesh = read_msh(FIXTURES / "no_physical_groups.msh")

    def test_fallback_to_mainland(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mesh = read_msh(FIXTURES / "no_physical_groups.msh")
        assert len(mesh.boundaries) >= 1
        for seg in mesh.boundaries:
            assert int(seg.bc_type) == int(BoundaryType.MAINLAND)


# ---------------------------------------------------------------------------
# Bathymetry
# ---------------------------------------------------------------------------


class TestBathymetry:
    """FR-016: $NodeData 'bathymetry' → Mesh.bathymetry with sign flip."""

    def test_bathymetry_loaded(self):
        mesh = read_msh(FIXTURES / "with_bathymetry_v22.msh")
        assert mesh.bathymetry is not None

    def test_bathymetry_sign_flip(self):
        """Disk stores positive-down; internal is positive-up, so signs flipped."""
        mesh = read_msh(FIXTURES / "with_bathymetry_v22.msh")
        # File stores [10,20,30,40,25]; internal should be negated
        assert mesh.bathymetry is not None
        assert np.all(mesh.bathymetry <= 0)  # positive-down → negative elevation

    def test_no_bathymetry_when_absent(self):
        mesh = read_msh(FIXTURES / "triangle_v22.msh")
        assert mesh.bathymetry is None


# ---------------------------------------------------------------------------
# Error cases
# ---------------------------------------------------------------------------


class TestErrors:
    """FR-004,010,011: GmshParseError on bad inputs."""

    def test_nonplanar_raises(self):
        with pytest.raises(GmshParseError, match="non-planar"):
            read_msh(FIXTURES / "nonplanar.msh")

    def test_higher_order_raises(self):
        with pytest.raises(GmshParseError, match="higher-order"):
            read_msh(FIXTURES / "higher_order.msh")

    def test_unsupported_version_raises(self):
        with pytest.raises(GmshParseError, match="unsupported Gmsh version"):
            read_msh(FIXTURES / "malformed_header.msh")

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            read_msh("/nonexistent/path/mesh.msh")

    def test_gmsh_parse_error_is_value_error(self):
        """GmshParseError must inherit ValueError (contract)."""
        assert issubclass(GmshParseError, ValueError)

    def test_gmsh_parse_error_attributes(self):
        err = GmshParseError("test detail", line_number=42)
        assert err.detail == "test detail"
        assert err.line_number == 42

    def test_binary_mode_raises(self, tmp_path):
        msh = tmp_path / "binary.msh"
        msh.write_text("$MeshFormat\n2.2 1 8\n$EndMeshFormat\n")
        with pytest.raises(GmshParseError, match="binary"):
            read_msh(msh)

    def test_top_level_import(self):
        """FR-017: GmshParseError accessible from top-level admesh namespace."""
        assert admesh.GmshParseError is GmshParseError


# ---------------------------------------------------------------------------
# Exports check
# ---------------------------------------------------------------------------


def test_read_msh_top_level_import():
    assert admesh.read_msh is read_msh


def test_physical_group_map_importable():
    from admesh.gmsh import PHYSICAL_GROUP_MAP as pgm
    assert isinstance(pgm, dict)
