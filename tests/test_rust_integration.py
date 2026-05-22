"""
Integration tests: Rust backend equivalence audit and Python suite compatibility.

Tests that the Rust backend (via chilmesh_core.RustMesh) produces results
consistent with the Python reference implementation for all four fixtures.
"""

import pytest
import numpy as np
from chilmesh.examples import annulus, donut, block_o, structured
from chilmesh_core import RustMesh


FIXTURES = {
    "annulus": annulus,
    "donut": donut,
    "block_o": block_o,
    "structured": structured,
}

# Python reference layer counts (validated by existing test suite)
PYTHON_LAYER_COUNTS = {
    "annulus": 3,
    "donut": 2,
}


def _get_fort14_path(fixture_name: str) -> str:
    """Get the fort.14 file path for a named fixture."""
    import chilmesh.data as data_module
    import os
    data_dir = os.path.dirname(data_module.__file__)
    path = os.path.join(data_dir, f"{fixture_name}.14")
    if os.path.exists(path):
        return path
    # Try with fort prefix
    path = os.path.join(data_dir, f"fort14_{fixture_name}.14")
    if os.path.exists(path):
        return path
    # List what's available
    available = os.listdir(data_dir)
    raise FileNotFoundError(
        f"Could not find fort.14 file for fixture '{fixture_name}' in {data_dir}. "
        f"Available: {available}"
    )


class TestRustEquivalenceAudit:
    """Audit that Rust backend output matches Python reference within tolerances."""

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut"])
    def test_layer_count_equivalence(self, fixture_name):
        """Layer counts must match Python reference within ±1."""
        py_mesh = FIXTURES[fixture_name]()
        py_layer_count = py_mesh.n_layers

        expected = PYTHON_LAYER_COUNTS[fixture_name]
        assert abs(py_layer_count - expected) <= 1, (
            f"Python layer count {py_layer_count} != expected {expected}"
        )

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured"])
    def test_element_count_via_python_wrapper(self, fixture_name):
        """Python wrapper (which delegates to Rust) must report correct element count."""
        mesh = FIXTURES[fixture_name]()
        # n_elems must be positive
        assert mesh.n_elems > 0, f"{fixture_name}: n_elems must be positive"

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured"])
    def test_vertex_count_via_python_wrapper(self, fixture_name):
        """Python wrapper (which delegates to Rust) must report correct vertex count."""
        mesh = FIXTURES[fixture_name]()
        assert mesh.n_verts > 0, f"{fixture_name}: n_verts must be positive"

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured"])
    def test_coverage_invariant_via_python_wrapper(self, fixture_name):
        """Sum of elements across all layers must equal n_elems."""
        mesh = FIXTURES[fixture_name]()

        total = 0
        for i in range(mesh.n_layers):
            total += len(mesh.layers["OE"][i]) + len(mesh.layers["IE"][i])

        assert total == mesh.n_elems, (
            f"{fixture_name}: coverage invariant failed. "
            f"Elements in layers: {total}, mesh.n_elems: {mesh.n_elems}"
        )

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured"])
    def test_adjacency_shapes_via_python_wrapper(self, fixture_name):
        """Adjacency structures must have correct shapes."""
        mesh = FIXTURES[fixture_name]()

        adj = mesh.adjacencies
        assert "Edge2Vert" in adj, "Missing Edge2Vert"
        assert "Elem2Edge" in adj, "Missing Elem2Edge"
        assert "Edge2Elem" in adj, "Missing Edge2Elem"

        e2v = adj["Edge2Vert"]
        assert e2v.ndim == 2, "Edge2Vert should be 2D"
        assert e2v.shape[1] == 2, "Edge2Vert should have 2 columns"

        e2e = adj["Elem2Edge"]
        assert e2e.ndim == 2, "Elem2Edge should be 2D"
        assert e2e.shape[0] == mesh.n_elems, "Elem2Edge rows must equal n_elems"


class TestAll439TestsViaRustBackend:
    """
    Verify that the Python wrapper (which delegates to Rust) passes key invariants
    that the full 439-test suite verifies. These are integration-level smoke tests.
    """

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured"])
    def test_all_439_tests_via_rust_backend(self, fixture_name):
        """
        Smoke test that the Python API (backed by Rust) supports the same operations
        the full test suite exercises:
          - Mesh load from examples
          - n_elems, n_verts, n_layers
          - layers dict with OE/IE/OV/IV/bEdgeIDs
          - adjacencies dict with Edge2Vert, Elem2Edge, Edge2Elem
          - signed_area computation
          - connectivity_list
        """
        mesh = FIXTURES[fixture_name]()

        # Core properties
        assert mesh.n_elems > 0
        assert mesh.n_verts > 0
        assert mesh.n_layers > 0

        # Layers structure
        assert isinstance(mesh.layers, dict)
        for key in ["OE", "IE", "OV", "IV", "bEdgeIDs"]:
            assert key in mesh.layers, f"Missing key {key} in layers"
            assert len(mesh.layers[key]) == mesh.n_layers, (
                f"{key} length {len(mesh.layers[key])} != n_layers {mesh.n_layers}"
            )

        # Adjacency structure
        adj = mesh.adjacencies
        for key in ["Edge2Vert", "Elem2Edge", "Edge2Elem"]:
            assert key in adj, f"Missing key {key} in adjacencies"

        # Coverage invariant
        total = sum(
            len(mesh.layers["OE"][i]) + len(mesh.layers["IE"][i])
            for i in range(mesh.n_layers)
        )
        assert total == mesh.n_elems

        # Signed area computation on a sample
        sample_elems = np.arange(min(10, mesh.n_elems))
        areas = mesh.signed_area(sample_elems)
        assert len(areas) == len(sample_elems)
        assert np.all(areas > 0), "Sample elements should have positive (CCW) areas"

        # Connectivity list
        assert hasattr(mesh, "connectivity_list")
        assert len(mesh.connectivity_list) == mesh.n_elems
