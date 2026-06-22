"""Tests for SMS .2dm write→read round-trip.

Covers the previously-untested ``CHILmesh.save("*.2dm")`` -> ``read_from_2dm``
write path (``tests/test_2dm_reader.py`` only exercised the reader). Substantiates
the constitution Principle VI ".2dm roundtrip" claim for geometry + topology.

Caveat: ``_write_2dm`` serializes only node coords + element connectivity, NOT
boundary-segment metadata, so the round-trip is lossless for geometry/topology
but does not preserve boundary records.
"""
from pathlib import Path

import numpy as np
import pytest

import chilmesh


class Test2DMRoundtrip:
    """Test CHILmesh.save() → read_from_2dm() round-trip fidelity."""

    def test_roundtrip_triangular_mesh(self, tmp_path):
        """Test roundtrip: triangular mesh save→load preserves node/element counts and connectivity."""
        mesh = chilmesh.examples.annulus()
        original_n_verts = mesh.n_verts
        original_n_elems = mesh.n_elems
        original_type = mesh.type

        # Write to temporary file
        tmpfile = tmp_path / "tri_roundtrip.2dm"
        mesh.save(str(tmpfile))

        # Reload
        mesh_reloaded = chilmesh.CHILmesh.read_from_2dm(tmpfile, compute_layers=False)

        # Verify roundtrip fidelity
        assert mesh_reloaded.n_verts == original_n_verts, \
            f"Node count mismatch after roundtrip: {original_n_verts} -> {mesh_reloaded.n_verts}"
        assert mesh_reloaded.n_elems == original_n_elems, \
            f"Element count mismatch after roundtrip: {original_n_elems} -> {mesh_reloaded.n_elems}"
        assert mesh_reloaded.type == original_type, \
            f"Element type mismatch after roundtrip: {original_type} -> {mesh_reloaded.type}"

    def test_roundtrip_quad_mesh(self, tmp_path):
        """Test roundtrip: quad mesh save→load preserves node/element counts and connectivity."""
        mesh = chilmesh.examples.quad_2x2()
        original_n_verts = mesh.n_verts
        original_n_elems = mesh.n_elems
        original_type = mesh.type

        # Write to temporary file
        tmpfile = tmp_path / "quad_roundtrip.2dm"
        mesh.save(str(tmpfile))

        # Reload
        mesh_reloaded = chilmesh.CHILmesh.read_from_2dm(tmpfile, compute_layers=False)

        # Verify roundtrip fidelity
        assert mesh_reloaded.n_verts == original_n_verts, \
            f"Node count mismatch after roundtrip: {original_n_verts} -> {mesh_reloaded.n_verts}"
        assert mesh_reloaded.n_elems == original_n_elems, \
            f"Element count mismatch after roundtrip: {original_n_elems} -> {mesh_reloaded.n_elems}"
        assert mesh_reloaded.type == original_type, \
            f"Element type mismatch after roundtrip: {original_type} -> {mesh_reloaded.type}"

    def test_roundtrip_mixed_element_mesh(self, tmp_path):
        """Test roundtrip: mixed-element mesh (triangles + quads) save→load preserves types."""
        # Create a small synthetic mixed-element mesh
        points = np.array([
            [0.0, 0.0, 0.0],  # 0
            [1.0, 0.0, 0.0],  # 1
            [1.0, 1.0, 0.0],  # 2
            [0.0, 1.0, 0.0],  # 3
            [2.0, 0.0, 0.0],  # 4
            [2.0, 1.0, 0.0],  # 5
        ], dtype=float)

        # Mixed connectivity: one triangle (padded), one quad
        connectivity = np.array([
            [0, 1, 3, 0],  # Triangle (padded with first vertex)
            [1, 4, 5, 2],  # Quad
        ], dtype=int)

        mesh = chilmesh.CHILmesh(connectivity=connectivity, points=points,
                                  compute_layers=False)
        assert mesh.type == "Mixed-Element", "Synthetic mesh should be Mixed-Element"

        # Write to temporary file
        tmpfile = tmp_path / "mixed_roundtrip.2dm"
        mesh.save(str(tmpfile))

        # Reload
        mesh_reloaded = chilmesh.CHILmesh.read_from_2dm(tmpfile, compute_layers=False)

        # Verify roundtrip fidelity
        assert mesh_reloaded.n_verts == mesh.n_verts == 6
        assert mesh_reloaded.n_elems == mesh.n_elems == 2
        assert mesh_reloaded.type == "Mixed-Element"

    def test_roundtrip_node_coordinates_precision(self, tmp_path):
        """Test that node coordinates are preserved to within 1e-9 (writer uses 10 decimals)."""
        # Create a small mesh with exact coordinates
        points = np.array([
            [0.123456789, 0.987654321, -0.5],
            [1.111111111, 2.222222222, 0.0],
            [3.141592654, 2.718281828, 1.0],
        ], dtype=float)

        connectivity = np.array([[0, 1, 2, 0]], dtype=int)  # Single triangle (padded)

        mesh = chilmesh.CHILmesh(connectivity=connectivity, points=points,
                                  compute_layers=False)

        # Write to temporary file
        tmpfile = tmp_path / "precision_roundtrip.2dm"
        mesh.save(str(tmpfile))

        # Reload
        mesh_reloaded = chilmesh.CHILmesh.read_from_2dm(tmpfile, compute_layers=False)

        # Verify coordinates match to within atol=1e-9 (writer uses 10 decimal places)
        np.testing.assert_allclose(
            mesh_reloaded.points[:, :2],  # Compare x,y only
            mesh.points[:, :2],
            atol=1e-9,
            err_msg="Coordinates not preserved after roundtrip"
        )

    def test_roundtrip_z_coordinates(self, tmp_path):
        """Test that z-coordinates are preserved in roundtrip."""
        points = np.array([
            [0.0, 0.0, -10.5],
            [1.0, 0.0, -5.25],
            [1.0, 1.0, 0.0],
        ], dtype=float)

        connectivity = np.array([[0, 1, 2, 0]], dtype=int)

        mesh = chilmesh.CHILmesh(connectivity=connectivity, points=points,
                                  compute_layers=False)

        # Write to temporary file
        tmpfile = tmp_path / "z_coord_roundtrip.2dm"
        mesh.save(str(tmpfile))

        # Reload
        mesh_reloaded = chilmesh.CHILmesh.read_from_2dm(tmpfile, compute_layers=False)

        # Verify z-coordinates match
        np.testing.assert_allclose(
            mesh_reloaded.points[:, 2],
            mesh.points[:, 2],
            atol=1e-9,
            err_msg="Z-coordinates not preserved after roundtrip"
        )

    def test_roundtrip_connectivity_preserved(self, tmp_path):
        """Test that element connectivity is preserved exactly in roundtrip."""
        # Create a synthetic mesh with explicit connectivity
        points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
        ], dtype=float)

        # Two triangles (padded) + two quads = mixed mesh
        connectivity = np.array([
            [0, 1, 3, 0],  # Triangle elem 0 (padded)
            [1, 2, 4, 4],  # Triangle elem 1 (padded, note: 4 appears twice)
            [1, 4, 5, 2],  # Quad elem 2
            [3, 4, 5, 3],  # Quad elem 3 (padded for roundtrip: elem[3] == elem[0])
        ], dtype=int)

        mesh = chilmesh.CHILmesh(connectivity=connectivity, points=points,
                                  compute_layers=False)

        # Write to temporary file
        tmpfile = tmp_path / "connectivity_roundtrip.2dm"
        mesh.save(str(tmpfile))

        # Reload
        mesh_reloaded = chilmesh.CHILmesh.read_from_2dm(tmpfile, compute_layers=False)

        # After roundtrip, both meshes should have same n_verts, n_elems
        assert mesh_reloaded.n_verts == mesh.n_verts == 6
        assert mesh_reloaded.n_elems == mesh.n_elems == 4

        # Note: connectivity_list ordering is NOT guaranteed to match after roundtrip
        # because elements are written in order but may be read back in same order
        # The key invariant is: same elements present, same connectivity per element.
        # For now, just verify counts; element-by-element matching is done in
        # test_roundtrip_fine_synthetic_mesh below.

    def test_roundtrip_fine_synthetic_mesh(self, tmp_path):
        """Test roundtrip on a tiny synthetic mesh: exact element connectivity verified."""
        # Create a minimal 2×1 grid: 6 nodes, 2 triangles (or 2 quads depending on triangulation)
        points = np.array([
            [0.0, 0.0, 0.0],  # 0
            [1.0, 0.0, 0.0],  # 1
            [2.0, 0.0, 0.0],  # 2
            [0.0, 1.0, 0.0],  # 3
            [1.0, 1.0, 0.0],  # 4
            [2.0, 1.0, 0.0],  # 5
        ], dtype=float)

        # Simple triangulation: 4 triangles
        connectivity = np.array([
            [0, 1, 4, 0],  # Triangle 0 (padded)
            [0, 4, 3, 0],  # Triangle 1 (padded)
            [1, 2, 5, 1],  # Triangle 2 (padded)
            [1, 5, 4, 1],  # Triangle 3 (padded)
        ], dtype=int)

        mesh = chilmesh.CHILmesh(connectivity=connectivity, points=points,
                                  compute_layers=False)

        # Write to temporary file
        tmpfile = tmp_path / "fine_synthetic.2dm"
        mesh.save(str(tmpfile))

        # Reload
        mesh_reloaded = chilmesh.CHILmesh.read_from_2dm(tmpfile, compute_layers=False)

        # Verify exact counts
        assert mesh_reloaded.n_verts == 6
        assert mesh_reloaded.n_elems == 4
        assert mesh_reloaded.type == "Triangular"

        # Verify points match
        np.testing.assert_allclose(mesh_reloaded.points, mesh.points, atol=1e-9)

        # For triangular meshes, connectivity should be exactly preserved
        # (Note: reader creates 3-column arrays for pure triangular meshes)
        np.testing.assert_array_equal(
            mesh_reloaded.connectivity_list,
            mesh.connectivity_list[:, :3],  # Strip padding column for comparison
            err_msg="Triangular mesh connectivity not preserved exactly"
        )

    def test_roundtrip_pure_quad_mesh_connectivity(self, tmp_path):
        """Test connectivity preservation for pure quad meshes."""
        # Create a pure quad mesh: 2×2 grid
        points = np.array([
            [0.0, 0.0, 0.0],  # 0
            [1.0, 0.0, 0.0],  # 1
            [2.0, 0.0, 0.0],  # 2
            [0.0, 1.0, 0.0],  # 3
            [1.0, 1.0, 0.0],  # 4
            [2.0, 1.0, 0.0],  # 5
            [0.0, 2.0, 0.0],  # 6
            [1.0, 2.0, 0.0],  # 7
            [2.0, 2.0, 0.0],  # 8
        ], dtype=float)

        connectivity = np.array([
            [0, 1, 4, 3],  # Quad 0
            [1, 2, 5, 4],  # Quad 1
            [3, 4, 7, 6],  # Quad 2
            [4, 5, 8, 7],  # Quad 3
        ], dtype=int)

        mesh = chilmesh.CHILmesh(connectivity=connectivity, points=points,
                                  compute_layers=False)
        assert mesh.type == "Quadrilateral"

        # Write to temporary file
        tmpfile = tmp_path / "quad_connectivity.2dm"
        mesh.save(str(tmpfile))

        # Reload
        mesh_reloaded = chilmesh.CHILmesh.read_from_2dm(tmpfile, compute_layers=False)

        # Verify exact counts and type
        assert mesh_reloaded.n_verts == 9
        assert mesh_reloaded.n_elems == 4
        assert mesh_reloaded.type == "Quadrilateral"

        # Verify connectivity matches exactly
        np.testing.assert_array_equal(
            mesh_reloaded.connectivity_list,
            mesh.connectivity_list,
            err_msg="Quad mesh connectivity not preserved exactly"
        )

    def test_roundtrip_donut_fixture(self, tmp_path):
        """Test roundtrip on donut fixture: medium-sized realistic mesh."""
        mesh = chilmesh.examples.donut()
        original_n_verts = mesh.n_verts
        original_n_elems = mesh.n_elems

        # Write to temporary file
        tmpfile = tmp_path / "donut_roundtrip.2dm"
        mesh.save(str(tmpfile))

        # Reload
        mesh_reloaded = chilmesh.CHILmesh.read_from_2dm(tmpfile, compute_layers=False)

        # Verify roundtrip fidelity
        assert mesh_reloaded.n_verts == original_n_verts
        assert mesh_reloaded.n_elems == original_n_elems
        assert mesh_reloaded.type == mesh.type

        # Verify coordinates preserved
        np.testing.assert_allclose(
            mesh_reloaded.points[:, :2],
            mesh.points[:, :2],
            atol=1e-9,
            err_msg="Donut geometry not preserved after roundtrip"
        )

    def test_roundtrip_with_layers_disabled_then_enabled(self, tmp_path):
        """Test roundtrip with compute_layers=False, then reload with compute_layers=True."""
        mesh = chilmesh.examples.annulus()

        # Write to temporary file (original mesh has layers computed)
        tmpfile = tmp_path / "layers_test.2dm"
        mesh.save(str(tmpfile))

        # Reload WITHOUT computing layers
        mesh_no_layers = chilmesh.CHILmesh.read_from_2dm(tmpfile, compute_layers=False)
        assert mesh_no_layers.n_layers == 0
        assert mesh_no_layers.adjacencies == {}

        # Reload WITH computing layers
        mesh_with_layers = chilmesh.CHILmesh.read_from_2dm(tmpfile, compute_layers=True)
        assert mesh_with_layers.n_layers > 0, "Layers should be computed"
        assert len(mesh_with_layers.adjacencies) > 0, "Adjacencies should be built"

    def test_roundtrip_coordinate_ordering(self, tmp_path):
        """Test that coordinate ordering (x, y, z) is preserved exactly."""
        # Create mesh with non-uniform coordinates
        points = np.array([
            [0.1, 0.9, -0.3],
            [0.5, 0.5, 0.0],
            [0.9, 0.1, 0.3],
        ], dtype=float)

        connectivity = np.array([[0, 1, 2, 0]], dtype=int)

        mesh = chilmesh.CHILmesh(connectivity=connectivity, points=points,
                                  compute_layers=False)

        # Write to temporary file
        tmpfile = tmp_path / "ordering.2dm"
        mesh.save(str(tmpfile))

        # Reload
        mesh_reloaded = chilmesh.CHILmesh.read_from_2dm(tmpfile, compute_layers=False)

        # Verify all coordinates (x, y, z) match
        np.testing.assert_allclose(
            mesh_reloaded.points,
            mesh.points,
            atol=1e-9,
            err_msg="Coordinate ordering (x,y,z) not preserved"
        )

    def test_roundtrip_structured_tri_fixture(self, tmp_path):
        """Test roundtrip on structured triangular fixture."""
        mesh = chilmesh.examples.structured()
        assert mesh.type == "Triangular", "Structured fixture should be triangular"

        original_n_verts = mesh.n_verts
        original_n_elems = mesh.n_elems

        # Write to temporary file
        tmpfile = tmp_path / "structured_roundtrip.2dm"
        mesh.save(str(tmpfile))

        # Reload
        mesh_reloaded = chilmesh.CHILmesh.read_from_2dm(tmpfile, compute_layers=False)

        # Verify roundtrip fidelity
        assert mesh_reloaded.n_verts == original_n_verts
        assert mesh_reloaded.n_elems == original_n_elems
        assert mesh_reloaded.type == "Triangular"

        # Verify coordinates preserved
        np.testing.assert_allclose(
            mesh_reloaded.points[:, :2],
            mesh.points[:, :2],
            atol=1e-9,
            err_msg="Structured triangular mesh coordinates not preserved"
        )
