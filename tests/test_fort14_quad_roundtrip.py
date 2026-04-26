"""Tests for quad and mixed-element mesh I/O (Issue #40)."""
import tempfile
from pathlib import Path

import pytest

import chilmesh


def test_quad_load_basic():
    """Test loading quad_2x2 fixture: verify node/element counts and element type."""
    mesh = chilmesh.examples.quad_2x2()

    assert mesh.n_verts == 9, f"Expected 9 nodes, got {mesh.n_verts}"
    assert mesh.n_elems == 4, f"Expected 4 elements, got {mesh.n_elems}"
    assert mesh.type == "Quadrilateral", f"Expected type='Quadrilateral', got {mesh.type}"
    assert mesh.connectivity_list.shape == (4, 4), f"Expected shape (4, 4), got {mesh.connectivity_list.shape}"


def test_quad_roundtrip():
    """Test roundtrip: load quad mesh, write to file, reload, verify counts."""
    mesh = chilmesh.examples.quad_2x2()
    original_n_verts = mesh.n_verts
    original_n_elems = mesh.n_elems
    original_type = mesh.type

    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fort.14', delete=False) as f:
        tmpfile = Path(f.name)

    try:
        # Write
        success = mesh.write_to_fort14(str(tmpfile))
        assert success, "Failed to write fort14 file"

        # Reload
        mesh_reloaded = chilmesh.CHILmesh.read_from_fort14(tmpfile)

        # Verify roundtrip fidelity
        assert mesh_reloaded.n_verts == original_n_verts, \
            f"Node count mismatch after roundtrip: {original_n_verts} -> {mesh_reloaded.n_verts}"
        assert mesh_reloaded.n_elems == original_n_elems, \
            f"Element count mismatch after roundtrip: {original_n_elems} -> {mesh_reloaded.n_elems}"
        assert mesh_reloaded.type == original_type, \
            f"Element type mismatch after roundtrip: {original_type} -> {mesh_reloaded.type}"
    finally:
        tmpfile.unlink()


def test_quad_compute_layers_false():
    """Test loading quad mesh with compute_layers=False for fast init."""
    import time

    # Load with layers disabled
    start = time.perf_counter()
    mesh = chilmesh.CHILmesh.read_from_fort14(
        chilmesh.examples.fixture_path("quad_2x2.fort.14"),
        compute_layers=False
    )
    elapsed = time.perf_counter() - start

    # Should be very fast
    assert elapsed < 1.0, f"Fast init took {elapsed:.2f}s, should be <1s"

    # Verify mesh loaded correctly
    assert mesh.n_verts == 9
    assert mesh.n_elems == 4
    assert mesh.type == "Quadrilateral"

    # Verify layers not computed
    assert mesh.n_layers == 0, f"Expected n_layers=0, got {mesh.n_layers}"
    assert mesh.adjacencies == {}, f"Expected empty adjacencies, got {mesh.adjacencies}"

    # Verify get_layer raises helpful error
    with pytest.raises(RuntimeError, match="Layers not computed"):
        mesh.get_layer(0)


def test_quad_padding_convention():
    """Test that quads and triangles use correct padding in 4-column arrays."""
    mesh = chilmesh.examples.quad_2x2()

    # All rows should have 4 columns (padded triangles or quads)
    assert mesh.connectivity_list.shape[1] == 4

    # All rows should be quads (no triangle padding needed in this fixture)
    for i, row in enumerate(mesh.connectivity_list):
        # For quads: all 4 values should be distinct
        unique_vals = len(set(row))
        assert unique_vals == 4, f"Row {i} should be a quad with 4 distinct vertices, got {row}"


def test_quad_bounding_box():
    """Test that bounding box is correct for quad mesh."""
    mesh = chilmesh.examples.quad_2x2()
    points = mesh.points

    # Quad_2x2 spans 0-2 in both x and y
    assert points[:, 0].min() == pytest.approx(0.0), "Min x should be 0.0"
    assert points[:, 0].max() == pytest.approx(2.0), "Max x should be 2.0"
    assert points[:, 1].min() == pytest.approx(0.0), "Min y should be 0.0"
    assert points[:, 1].max() == pytest.approx(2.0), "Max y should be 2.0"
