"""Integration tests for bridge adapters.

Tests verify that bridge adapters (MADMESHR, ADMESH, ADMESH-Domains)
work correctly with realistic downstream workflows and all mesh fixtures.
"""

import pytest
import numpy as np
from chilmesh import CHILmesh
from chilmesh.bridge import (
    MeshAdapterForMADMESHR,
    MeshAdapterForADMESH,
    MeshAdapterForADMESHDomains,
)


class TestMADMESHRBridge:
    """Integration tests for MADMESHR bridge adapter."""

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_adapter_initialization(self, fixture_name):
        """Verify adapter can be initialized on all fixtures."""
        mesh = self._load_fixture(fixture_name)
        adapter = MeshAdapterForMADMESHR(mesh)
        assert adapter.mesh is mesh

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_get_element_neighbors_valid_range(self, fixture_name):
        """Verify element neighbor queries work for valid indices."""
        mesh = self._load_fixture(fixture_name)
        adapter = MeshAdapterForMADMESHR(mesh)

        # Test a few different elements
        for elem_id in [0, mesh.n_elems // 2, mesh.n_elems - 1]:
            neighbors = adapter.get_element_neighbors(elem_id)
            assert isinstance(neighbors, set)
            # Neighbors should not include the element itself
            assert elem_id not in neighbors
            # All neighbor IDs should be valid
            for n in neighbors:
                assert 0 <= n < mesh.n_elems

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_get_element_neighbors_invalid_range(self, fixture_name):
        """Verify element neighbor queries raise for invalid indices."""
        mesh = self._load_fixture(fixture_name)
        adapter = MeshAdapterForMADMESHR(mesh)

        with pytest.raises(ValueError):
            adapter.get_element_neighbors(-1)

        with pytest.raises(ValueError):
            adapter.get_element_neighbors(mesh.n_elems)

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_get_element_quality_neighborhood(self, fixture_name):
        """Verify element quality neighborhood computation."""
        mesh = self._load_fixture(fixture_name)
        adapter = MeshAdapterForMADMESHR(mesh)

        metrics = adapter.get_element_quality_neighborhood(0)

        assert "center" in metrics
        assert "neighbors" in metrics
        assert "mean" in metrics
        assert "stats" in metrics

        assert isinstance(metrics["center"], (float, np.floating))
        assert isinstance(metrics["mean"], float)
        assert 0 <= metrics["center"] <= 1
        assert 0 <= metrics["mean"] <= 1

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_get_refinement_region(self, fixture_name):
        """Verify refinement region generation."""
        mesh = self._load_fixture(fixture_name)
        adapter = MeshAdapterForMADMESHR(mesh)

        # Get low quality elements
        quality, _, _ = mesh.elem_quality()
        poor = np.where(quality < 0.5)[0]

        if len(poor) > 0:
            # Get refinement region with neighbors
            region = adapter.get_refinement_region(list(poor[:3]))
            assert isinstance(region, set)
            assert len(region) >= 3  # At least the seed elements
            # Check all elements are valid
            for e in region:
                assert 0 <= e < mesh.n_elems

            # Region without neighbors should be smaller or equal
            region_no_neighbors = adapter.get_refinement_region(
                list(poor[:3]), include_neighbors=False
            )
            assert len(region_no_neighbors) <= len(region)

    def _load_fixture(self, name):
        """Load a test fixture mesh."""
        if name == "annulus":
            from pathlib import Path

            return CHILmesh.read_from_fort14(
                Path("src/chilmesh/data/annulus_200pts.fort.14")
            )
        elif name == "donut":
            from pathlib import Path

            return CHILmesh.read_from_fort14(
                Path("src/chilmesh/data/donut_domain.fort.14")
            )
        elif name == "structured":
            from pathlib import Path

            return CHILmesh.read_from_fort14(
                Path("src/chilmesh/data/quad_2x2.fort.14")
            )
        else:
            raise ValueError(f"Unknown fixture: {name}")


class TestADMESHBridge:
    """Integration tests for ADMESH bridge adapter."""

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_adapter_initialization(self, fixture_name):
        """Verify adapter can be initialized on all fixtures."""
        mesh = self._load_fixture(fixture_name)
        adapter = MeshAdapterForADMESH(mesh)
        assert adapter.mesh is mesh

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_get_mesh_quality_report(self, fixture_name):
        """Verify mesh quality report generation."""
        mesh = self._load_fixture(fixture_name)
        adapter = MeshAdapterForADMESH(mesh)

        report = adapter.get_mesh_quality_report()

        # Check required keys
        required_keys = [
            "mean",
            "min",
            "max",
            "std",
            "poor_count",
            "poor_fraction",
            "poor_elements",
        ]
        for key in required_keys:
            assert key in report, f"Missing key: {key}"

        # Check value ranges
        assert 0 <= report["min"] <= report["mean"] <= report["max"] <= 1
        assert 0 <= report["poor_fraction"] <= 1
        assert len(report["poor_elements"]) == report["poor_count"]

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_get_element_angles_summary(self, fixture_name):
        """Verify element angle summary computation."""
        mesh = self._load_fixture(fixture_name)
        adapter = MeshAdapterForADMESH(mesh)

        summary = adapter.get_element_angles_summary()

        required_keys = [
            "min_angle_degrees",
            "max_angle_degrees",
            "elements_with_acute",
            "elements_with_obtuse",
        ]
        for key in required_keys:
            assert key in summary

        # Check angle ranges (0 to 180 degrees)
        assert 0 < summary["min_angle_degrees"] < 180
        assert 0 < summary["max_angle_degrees"] < 180
        assert summary["min_angle_degrees"] <= summary["max_angle_degrees"]

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_quality_report_consistency(self, fixture_name):
        """Verify quality report matches direct elem_quality call."""
        mesh = self._load_fixture(fixture_name)
        adapter = MeshAdapterForADMESH(mesh)

        report = adapter.get_mesh_quality_report()
        quality, _, stats = mesh.elem_quality()

        assert np.isclose(report["mean"], stats["mean"], atol=1e-5)
        assert np.isclose(report["min"], stats["min"], atol=1e-5)
        assert np.isclose(report["max"], stats["max"], atol=1e-5)

    def _load_fixture(self, name):
        """Load a test fixture mesh."""
        if name == "annulus":
            from pathlib import Path

            return CHILmesh.read_from_fort14(
                Path("src/chilmesh/data/annulus_200pts.fort.14")
            )
        elif name == "donut":
            from pathlib import Path

            return CHILmesh.read_from_fort14(
                Path("src/chilmesh/data/donut_domain.fort.14")
            )
        elif name == "structured":
            from pathlib import Path

            return CHILmesh.read_from_fort14(
                Path("src/chilmesh/data/quad_2x2.fort.14")
            )
        else:
            raise ValueError(f"Unknown fixture: {name}")


class TestADMESHDomainsBridge:
    """Integration tests for ADMESH-Domains bridge adapter."""

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_adapter_initialization(self, fixture_name):
        """Verify adapter can be initialized on all fixtures."""
        mesh = self._load_fixture(fixture_name)
        adapter = MeshAdapterForADMESHDomains(mesh)
        assert adapter.mesh is mesh

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_get_domain_boundaries(self, fixture_name):
        """Verify domain boundary extraction."""
        mesh = self._load_fixture(fixture_name)
        adapter = MeshAdapterForADMESHDomains(mesh)

        boundaries = adapter.get_domain_boundaries()

        assert isinstance(boundaries, dict)
        assert 0 in boundaries  # Single domain
        assert isinstance(boundaries[0], set)

        # Boundary vertices should be valid
        for vert_id in boundaries[0]:
            assert 0 <= vert_id < mesh.n_verts

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_get_mesh_connectivity_info(self, fixture_name):
        """Verify mesh connectivity information."""
        mesh = self._load_fixture(fixture_name)
        adapter = MeshAdapterForADMESHDomains(mesh)

        info = adapter.get_mesh_connectivity_info()

        required_keys = [
            "n_verts",
            "n_elems",
            "n_edges",
            "avg_vertex_degree",
            "avg_element_degree",
            "boundary_edge_count",
            "interior_edge_count",
        ]
        for key in required_keys:
            assert key in info

        # Check values are reasonable
        assert info["n_verts"] == mesh.n_verts
        assert info["n_elems"] == mesh.n_elems
        assert info["n_edges"] == mesh.n_edges
        assert info["boundary_edge_count"] > 0
        assert (
            info["boundary_edge_count"] + info["interior_edge_count"]
            == info["n_edges"]
        )
        assert info["avg_vertex_degree"] > 0
        assert info["avg_element_degree"] > 0

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_boundary_consistency_with_direct_call(self, fixture_name):
        """Verify boundary vertices match direct boundary_edges call."""
        mesh = self._load_fixture(fixture_name)
        adapter = MeshAdapterForADMESHDomains(mesh)

        boundaries = adapter.get_domain_boundaries()
        boundary_verts = boundaries[0]

        # Get boundary directly
        boundary_edges = mesh.boundary_edges()
        direct_boundary_verts = set()
        for edge_id in boundary_edges:
            v1, v2 = mesh.edge2vert(edge_id)[0]
            direct_boundary_verts.add(v1)
            direct_boundary_verts.add(v2)

        # Should match
        assert boundary_verts == direct_boundary_verts

    def _load_fixture(self, name):
        """Load a test fixture mesh."""
        if name == "annulus":
            from pathlib import Path

            return CHILmesh.read_from_fort14(
                Path("src/chilmesh/data/annulus_200pts.fort.14")
            )
        elif name == "donut":
            from pathlib import Path

            return CHILmesh.read_from_fort14(
                Path("src/chilmesh/data/donut_domain.fort.14")
            )
        elif name == "structured":
            from pathlib import Path

            return CHILmesh.read_from_fort14(
                Path("src/chilmesh/data/quad_2x2.fort.14")
            )
        else:
            raise ValueError(f"Unknown fixture: {name}")


class TestBridgeAdaptersErrorHandling:
    """Test error handling in bridge adapters."""

    def test_madmeshr_invalid_element_id(self):
        """Verify MADMESHR adapter handles invalid element IDs."""
        from pathlib import Path

        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )
        adapter = MeshAdapterForMADMESHR(mesh)

        with pytest.raises(ValueError):
            adapter.get_element_neighbors(-1)

        with pytest.raises(ValueError):
            adapter.get_element_neighbors(mesh.n_elems + 10)

    def test_madmeshr_quality_neighborhood_invalid(self):
        """Verify quality neighborhood handles invalid IDs."""
        from pathlib import Path

        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )
        adapter = MeshAdapterForMADMESHR(mesh)

        with pytest.raises(ValueError):
            adapter.get_element_quality_neighborhood(mesh.n_elems)


class TestBridgeAdaptersCommonWorkflows:
    """Test realistic downstream workflows using bridge adapters."""

    def test_madmeshr_refinement_workflow(self):
        """Simulate typical MADMESHR mesh refinement workflow."""
        from pathlib import Path

        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )
        adapter = MeshAdapterForMADMESHR(mesh)

        # Identify poor quality elements
        quality, _, _ = mesh.elem_quality()
        poor_quality = np.where(quality < 0.4)[0]

        if len(poor_quality) > 0:
            # Get refinement region
            region = adapter.get_refinement_region(list(poor_quality))

            # Verify region is valid
            assert len(region) >= len(poor_quality)
            for e in region:
                assert 0 <= e < mesh.n_elems

    def test_admesh_quality_assessment_workflow(self):
        """Simulate typical ADMESH quality assessment workflow."""
        from pathlib import Path

        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )
        adapter = MeshAdapterForADMESH(mesh)

        # Get quality report
        report = adapter.get_mesh_quality_report()

        # Get angles
        angles = adapter.get_element_angles_summary()

        # Use these for adaptation decisions
        if report["poor_count"] > 0:
            assert angles["elements_with_acute"] >= 0

    def test_admesh_domains_multi_domain_setup(self):
        """Simulate ADMESH-Domains domain initialization."""
        from pathlib import Path

        mesh = CHILmesh.read_from_fort14(
            Path("src/chilmesh/data/annulus_200pts.fort.14")
        )
        adapter = MeshAdapterForADMESHDomains(mesh)

        # Get boundaries for domain setup
        boundaries = adapter.get_domain_boundaries()
        connectivity = adapter.get_mesh_connectivity_info()

        # Verify both are consistent
        assert len(boundaries) > 0
        assert connectivity["n_verts"] == mesh.n_verts
        assert connectivity["n_elems"] == mesh.n_elems
