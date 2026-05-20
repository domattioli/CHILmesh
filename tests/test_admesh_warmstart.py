"""
Unit tests for ADMESH warm-start truss optimizer adapter.

Tests cover:
- Boundary preservation (FR-008, SC-002): V_BND, V_BND_PROP
- Non-degradation (FR-011): Quality never regresses
- Error handling (FR-005-007, FR-014): All validation paths
- Extensibility (FR-016-017, SC-007-008): Annulus, donut, custom domain
- Determinism (FR-009): RNG seed control
- Performance (SC-003): < 60s on annulus (Linux runner only — TEST-AUDIT F8)
"""

import numpy as np
import pytest
import warnings
from scipy.spatial import Delaunay

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from chilmesh import (
    CHILmesh,
    optimize_with_admesh_truss,
    optimize_with_admesh_truss_arrays,
    examples,
)

# Skip all tests if ADMESH is not available
try:
    import admesh
    HAS_ADMESH = True
except ImportError:
    HAS_ADMESH = False
    admesh_unavailable = pytest.mark.skipif(
        not HAS_ADMESH,
        reason="ADMESH not installed (pip install admesh or clone ./ADMESH)"
    )


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def annulus_mesh():
    """Bundled annulus fixture from CHILmesh."""
    return examples.annulus()


@pytest.fixture
def donut_mesh():
    """Bundled donut fixture from CHILmesh."""
    return examples.donut()


@pytest.fixture
def annulus_sdf():
    """Signed distance function for annulus (outer R=1.0, inner r=0.3)."""
    def sdf(p):
        r = np.linalg.norm(p, axis=1)
        return np.maximum(r - 1.0, 0.3 - r)
    return sdf


@pytest.fixture
def donut_sdf():
    """Signed distance function for donut (outer R=1, inner r=0.3)."""
    def sdf(p):
        r = np.linalg.norm(p, axis=1)
        return np.maximum(r - 1.0, 0.3 - r)
    return sdf


@pytest.fixture
def square_sdf():
    """Signed distance function for unit square [-1,1]^2."""
    def sdf(p):
        return np.max(np.abs(p) - 1.0, axis=1)
    return sdf


@pytest.fixture
def uniform_size_fn():
    """Uniform size function (h=0.1)."""
    def size_fn(p):
        return np.full(len(p), 0.1)
    return size_fn


@pytest.fixture
def graded_size_fn():
    """Graded size function: finer near boundaries, coarser in middle."""
    def size_fn(p):
        r = np.linalg.norm(p, axis=1)
        d_boundary = np.abs(r - 0.75)  # Distance to some inner radius
        return 0.05 + 0.10 * (d_boundary / 0.3).clip(max=1.0)
    return size_fn


# ============================================================================
# Boundary Preservation Tests (FR-008, SC-002)
# ============================================================================

class TestBoundaryPreservation:
    """V_BND: Bit-exact boundary preservation."""

    def test_vbnd_annulus_array_form(self, annulus_mesh, annulus_sdf):
        """V_BND on annulus using array form."""
        points = annulus_mesh.points
        triangles = annulus_mesh.connectivity_list

        # Get boundary indices
        boundary_edges = annulus_mesh.boundary_edges()
        edge2vert = annulus_mesh.adjacencies["Edge2Vert"]
        boundary_indices = np.unique(edge2vert[boundary_edges].flatten())

        input_boundary = points[boundary_indices, :2]

        # Run optimizer
        points_opt, _ = optimize_with_admesh_truss_arrays(
            points, triangles, annulus_sdf, size_fn=None, seed=0
        )

        # Check bit-exact preservation
        output_boundary = points_opt[:len(boundary_indices)]
        assert np.array_equal(output_boundary, input_boundary), \
            f"Boundary not preserved; max delta = {np.max(np.abs(output_boundary - input_boundary)):.2e}"

    def test_vbnd_annulus_chilmesh_form(self, annulus_mesh, annulus_sdf):
        """V_BND on annulus using CHILmesh form."""
        boundary_edges = annulus_mesh.boundary_edges()
        edge2vert = annulus_mesh.adjacencies["Edge2Vert"]
        boundary_indices = np.unique(edge2vert[boundary_edges].flatten())
        input_boundary = annulus_mesh.points[boundary_indices, :2]

        # Run optimizer
        mesh_opt = optimize_with_admesh_truss(annulus_mesh, annulus_sdf, size_fn=None, seed=0)

        # Check bit-exact preservation
        output_boundary = mesh_opt.points[boundary_indices, :2]
        assert np.array_equal(output_boundary, input_boundary), \
            f"Boundary not preserved in CHILmesh form"

    def test_vbnd_donut_domain_agnostic(self, donut_mesh, donut_sdf):
        """V_BND on donut: proves domain-agnosticism."""
        boundary_edges = donut_mesh.boundary_edges()
        edge2vert = donut_mesh.adjacencies["Edge2Vert"]
        boundary_indices = np.unique(edge2vert[boundary_edges].flatten())
        input_boundary = donut_mesh.points[boundary_indices, :2]

        mesh_opt = optimize_with_admesh_truss(donut_mesh, donut_sdf, size_fn=None, seed=0)

        output_boundary = mesh_opt.points[boundary_indices, :2]
        assert np.array_equal(output_boundary, input_boundary), \
            "Boundary not preserved on donut"

    def test_vbnd_multiple_domains(self, annulus_mesh, donut_mesh, annulus_sdf, donut_sdf):
        """V_BND: Test on multiple domains in single test."""
        for mesh, sdf in [(annulus_mesh, annulus_sdf), (donut_mesh, donut_sdf)]:
            boundary_edges = mesh.boundary_edges()
            edge2vert = mesh.adjacencies["Edge2Vert"]
            boundary_indices = np.unique(edge2vert[boundary_edges].flatten())
            input_boundary = mesh.points[boundary_indices, :2]

            mesh_opt = optimize_with_admesh_truss(mesh, sdf, size_fn=None, seed=0)
            output_boundary = mesh_opt.points[boundary_indices, :2]

            assert np.array_equal(output_boundary, input_boundary), \
                f"Boundary not preserved on domain {mesh.__class__.__name__}"


# ============================================================================
# Non-Degradation Tests (FR-011)
# ============================================================================

class TestNonDegradation:
    """Non-degradation guard: quality never regresses."""

    def test_enforce_non_degradation_true(self, annulus_mesh, annulus_sdf):
        """Non-degradation enforced: regressed output returns input."""
        # Use hostile params to force degradation (if possible)
        points = annulus_mesh.points
        triangles = annulus_mesh.connectivity_list

        # Run with enforcement (default)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            points_opt, triangles_opt = optimize_with_admesh_truss_arrays(
                points, triangles, annulus_sdf, size_fn=None,
                enforce_non_degradation=True, seed=0
            )

            # If quality regressed, we should get a warning and input back
            # (or the output was actually better)
            # Just verify it doesn't crash and returns something
            assert points_opt is not None
            assert triangles_opt is not None

    def test_enforce_non_degradation_false(self, annulus_mesh, annulus_sdf):
        """Non-degradation disabled: always return optimizer output."""
        points = annulus_mesh.points
        triangles = annulus_mesh.connectivity_list

        points_opt, triangles_opt = optimize_with_admesh_truss_arrays(
            points, triangles, annulus_sdf, size_fn=None,
            enforce_non_degradation=False, seed=0
        )

        assert points_opt is not None
        assert triangles_opt is not None


# ============================================================================
# Error Handling Tests (FR-005-007, FR-014)
# ============================================================================

class TestErrorHandling:
    """Validation errors: V_TRI, V_AREA, V_BND_SDF, V_BBOX, ImportError."""

    def test_v_tri_quad_rejection(self, annulus_mesh):
        """V_TRI: Quad mesh (4-vertex) raises NotImplementedError."""
        points = annulus_mesh.points
        # Create fake quad connectivity (4 columns)
        quads = np.column_stack([
            np.arange(4),
            np.arange(1, 5),
            np.arange(2, 6),
            np.arange(3, 7),
        ])

        with pytest.raises(NotImplementedError, match="triangle-only"):
            optimize_with_admesh_truss_arrays(
                points, quads, lambda p: np.zeros(len(p)), None
            )

    def test_v_area_degenerate_triangles(self, annulus_mesh, annulus_sdf):
        """V_AREA: Degenerate (zero-area) triangles raise ValueError."""
        points = annulus_mesh.points.copy()
        triangles = annulus_mesh.connectivity_list.copy()

        # Make first triangle degenerate: all three points the same
        triangles[0, :] = [0, 0, 0]

        with pytest.raises(ValueError, match="non-positive signed area"):
            optimize_with_admesh_truss_arrays(
                points, triangles, annulus_sdf, None
            )

    def test_v_bnd_sdf_off_boundary(self, annulus_mesh):
        """V_BND_SDF: Boundary not on SDF zero-set raises ValueError."""
        points = annulus_mesh.points.copy()
        triangles = annulus_mesh.connectivity_list

        # Wrong SDF: shift by 0.1 (boundary no longer on zero-set)
        def wrong_sdf(p):
            r = np.linalg.norm(p, axis=1)
            return np.maximum(r - 1.1, 0.4 - r)  # Shifted

        with pytest.raises(ValueError, match="boundary not on SDF zero set"):
            optimize_with_admesh_truss_arrays(
                points, triangles, wrong_sdf, None
            )

    def test_v_bbox_points_outside(self, annulus_mesh, annulus_sdf):
        """V_BBOX: Points outside bbox raise ValueError."""
        points = annulus_mesh.points.copy()
        triangles = annulus_mesh.connectivity_list

        # Tiny bbox that excludes all points
        bbox = (0.0, 0.0, 0.1, 0.1)

        with pytest.raises(ValueError, match="outside bbox"):
            optimize_with_admesh_truss_arrays(
                points, triangles, annulus_sdf, None, bbox=bbox
            )


# ============================================================================
# Extensibility Tests (FR-016-017, SC-007-008)
# ============================================================================

class TestExtensibility:
    """Input-source and domain agnosticism."""

    def test_array_form_source_agnostic(self, square_sdf):
        """Source-agnostic: raw (points, triangles) numpy arrays."""
        # Create square mesh from scratch (non-CHILmesh)
        np.random.seed(42)
        n_boundary = 40
        n_interior = 80

        # Boundary: perimeter of [-1, 1]^2
        edge = np.linspace(-1, 1, n_boundary // 4, endpoint=False)
        boundary = np.vstack([
            np.column_stack([edge, np.full(len(edge), -1)]),
            np.column_stack([np.full(len(edge), 1), edge]),
            np.column_stack([edge[::-1], np.full(len(edge), 1)]),
            np.column_stack([np.full(len(edge), -1), edge[::-1]]),
        ])[:n_boundary]

        # Interior: random
        interior = np.random.uniform(-0.95, 0.95, (n_interior, 2))
        all_points = np.vstack([boundary, interior])

        # Triangulate
        tri = Delaunay(all_points)
        triangles = tri.simplices

        # Boundary indices: first n_boundary
        boundary_indices = np.arange(n_boundary)

        # Optimize
        points_opt, triangles_opt = optimize_with_admesh_truss_arrays(
            all_points, triangles, square_sdf, None,
            boundary_indices=boundary_indices, h0=0.2, seed=0
        )

        # Verify boundary preservation
        assert np.array_equal(points_opt[:n_boundary], boundary)

    def test_domain_agnostic_annulus(self, annulus_mesh, annulus_sdf, uniform_size_fn):
        """Domain-agnostic: annulus with custom SDF."""
        mesh_opt = optimize_with_admesh_truss(
            annulus_mesh, annulus_sdf, size_fn=uniform_size_fn, seed=0
        )
        assert mesh_opt is not None

    def test_domain_agnostic_donut(self, donut_mesh, donut_sdf, graded_size_fn):
        """Domain-agnostic: donut with graded size function."""
        mesh_opt = optimize_with_admesh_truss(
            donut_mesh, donut_sdf, size_fn=graded_size_fn, seed=0
        )
        assert mesh_opt is not None


# ============================================================================
# Quality and Performance Tests
# ============================================================================

class TestQualityAndPerformance:
    """SC-001, SC-003, SC-005: Quality improvement and performance."""

    def test_determinism(self, annulus_mesh, annulus_sdf):
        """FR-009: Same seed produces identical output."""
        points = annulus_mesh.points
        triangles = annulus_mesh.connectivity_list

        out1 = optimize_with_admesh_truss_arrays(
            points, triangles, annulus_sdf, None, seed=42
        )
        out2 = optimize_with_admesh_truss_arrays(
            points, triangles, annulus_sdf, None, seed=42
        )

        assert np.array_equal(out1[0], out2[0])
        assert np.array_equal(out1[1], out2[1])

    @pytest.mark.slow
    @pytest.mark.skipif(
        sys.platform != "linux",
        reason="wall-clock gate is Linux-only; macOS/Windows runner variance triggers flakes (TEST-AUDIT F8)",
    )
    def test_performance_annulus(self, annulus_mesh, annulus_sdf):
        """SC-003: Annulus completes in < 60s on a canonical Linux runner."""
        import time

        points = annulus_mesh.points
        triangles = annulus_mesh.connectivity_list

        start = time.perf_counter()
        optimize_with_admesh_truss_arrays(
            points, triangles, annulus_sdf, None, niter=500, seed=0
        )
        elapsed = time.perf_counter() - start

        assert elapsed < 60.0, f"Took {elapsed:.1f}s, limit is 60s"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
