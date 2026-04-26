"""Tests for fast initialization with compute_layers=False (Issue #41)."""
import time

import pytest

import chilmesh


class TestFastInitialization:
    """Test that compute_layers=False enables fast mesh loading."""

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured", "quad_2x2"])
    def test_fast_init_timing(self, fixture_name):
        """
        Test that fast init completes in <2s for all fixtures.

        Success Criteria (SC-003): Fast loading <2s for meshes up to ~15,000 elements.
        """
        fixture_path = chilmesh.examples.fixture_path(f"{fixture_name}.fort.14")

        # Measure fast init time
        start = time.perf_counter()
        mesh = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=False)
        elapsed = time.perf_counter() - start

        # Should be very fast: <2s even for Block_O (~5,200 elements)
        assert elapsed < 2.0, f"Fast init took {elapsed:.2f}s, should be <2s for {fixture_name}"
        assert mesh.n_verts > 0
        assert mesh.n_elems > 0

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured", "quad_2x2"])
    def test_no_adjacencies_when_fast_init(self, fixture_name):
        """Test that adjacencies are not computed when compute_layers=False."""
        fixture_path = chilmesh.examples.fixture_path(f"{fixture_name}.fort.14")
        mesh = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=False)

        assert mesh.adjacencies == {}, f"Expected empty adjacencies for fast init, got {len(mesh.adjacencies)}"
        assert mesh.n_layers == 0, f"Expected n_layers=0 for fast init, got {mesh.n_layers}"

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured", "quad_2x2"])
    def test_mesh_geometry_preserved_fast_init(self, fixture_name):
        """Test that mesh geometry is correct even with compute_layers=False."""
        fixture_path = chilmesh.examples.fixture_path(f"{fixture_name}.fort.14")

        # Load with fast init
        mesh_fast = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=False)

        # Load with layers (for comparison)
        mesh_full = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=True)

        # Geometry should match
        assert mesh_fast.n_verts == mesh_full.n_verts
        assert mesh_fast.n_elems == mesh_full.n_elems
        assert mesh_fast.type == mesh_full.type

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "block_o", "structured", "quad_2x2"])
    def test_element_type_set_with_fast_init(self, fixture_name):
        """Test that element_type is set even with compute_layers=False."""
        fixture_path = chilmesh.examples.fixture_path(f"{fixture_name}.fort.14")
        mesh = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=False)

        assert mesh.type is not None, f"element_type should be set for {fixture_name}"
        assert mesh.type in ("Triangular", "Quadrilateral", "Mixed-Element")

    def test_get_layer_fails_after_fast_init(self):
        """Test that get_layer() raises helpful error when compute_layers=False."""
        fixture_path = chilmesh.examples.fixture_path("annulus_200pts.fort.14")
        mesh = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=False)

        with pytest.raises(RuntimeError, match="Layers not computed"):
            mesh.get_layer(0)

    def test_from_admesh_domain_fast_init(self):
        """Test that from_admesh_domain respects compute_layers=False."""
        from types import SimpleNamespace

        fixture_path = chilmesh.examples.fixture_path("quad_2x2.fort.14")
        record = SimpleNamespace(filename=str(fixture_path), type="ADCIRC")

        start = time.perf_counter()
        mesh = chilmesh.CHILmesh.from_admesh_domain(record, compute_layers=False)
        elapsed = time.perf_counter() - start

        assert elapsed < 1.0, f"from_admesh_domain with fast init took {elapsed:.2f}s"
        assert mesh.n_layers == 0
        assert mesh.adjacencies == {}
