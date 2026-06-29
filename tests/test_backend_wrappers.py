"""Tests for backend wrapper modules (cpp_backend, rust_backend).

Covers both unavailable (natural False) and available (via monkeypatch) paths,
exercising error handling and marshalling logic.
"""
from __future__ import annotations

import types
from typing import Any

import numpy as np
import pytest

import chilmesh.backends.cpp_backend as cpp_backend
import chilmesh.backends.rust_backend as rust_backend


# ============================================================================
# C++ BACKEND TESTS
# ============================================================================

class TestCppBackendUnavailable:
    """Test C++ backend error paths when CPP_AVAILABLE is False.

    Relies on natural False (no extension installed) but guards with skip
    if the extension is present on the test machine.
    """

    def test_require_cpp_raises_when_unavailable(self):
        """_require_cpp() raises ImportError with documented message."""
        if cpp_backend.CPP_AVAILABLE:
            pytest.skip("C++ extension present on this machine")

        with pytest.raises(ImportError) as exc_info:
            cpp_backend._require_cpp()

        assert "chilmesh_cpp C++ extension not available" in str(exc_info.value)

    def test_fast_init_raises_when_unavailable(self):
        """fast_init() raises ImportError when C++ not available."""
        if cpp_backend.CPP_AVAILABLE:
            pytest.skip("C++ extension present")

        pts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        conn = np.array([[0, 1, 2]], dtype=np.int64)

        with pytest.raises(ImportError) as exc_info:
            cpp_backend.fast_init(pts, conn)

        assert "chilmesh_cpp C++ extension not available" in str(exc_info.value)

    def test_full_init_raises_when_unavailable(self):
        """full_init() raises ImportError when C++ not available."""
        if cpp_backend.CPP_AVAILABLE:
            pytest.skip("C++ extension present")

        pts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        conn = np.array([[0, 1, 2]], dtype=np.int64)

        with pytest.raises(ImportError) as exc_info:
            cpp_backend.full_init(pts, conn)

        assert "chilmesh_cpp C++ extension not available" in str(exc_info.value)

    def test_quality_analysis_raises_when_unavailable(self):
        """quality_analysis() raises ImportError when C++ not available."""
        if cpp_backend.CPP_AVAILABLE:
            pytest.skip("C++ extension present")

        fake_mesh = object()
        with pytest.raises(ImportError) as exc_info:
            cpp_backend.quality_analysis(fake_mesh)

        assert "chilmesh_cpp C++ extension not available" in str(exc_info.value)

    def test_get_vertex_edges_raises_when_unavailable(self):
        """get_vertex_edges() raises ImportError when C++ not available."""
        if cpp_backend.CPP_AVAILABLE:
            pytest.skip("C++ extension present")

        fake_mesh = object()
        with pytest.raises(ImportError) as exc_info:
            cpp_backend.get_vertex_edges(fake_mesh, 0)

        assert "chilmesh_cpp C++ extension not available" in str(exc_info.value)


class TestCppBackendAvailableViaMonkeypatch:
    """Test C++ backend with simulated available backend (monkeypatch).

    Uses fake backend module to verify marshalling and forwarding logic.
    """

    def test_fast_init_marshalls_points_and_connectivity(self, monkeypatch):
        """fast_init() converts points to float64 C-order and conn to int32 C-order."""
        # Create a fake backend module
        class FakeCpp:
            call_record = {}

            @staticmethod
            def fast_init(points, connectivity):
                FakeCpp.call_record['points'] = points
                FakeCpp.call_record['connectivity'] = connectivity
                return "sentinel_mesh"

        # Monkeypatch
        monkeypatch.setattr(cpp_backend, "CPP_AVAILABLE", True)
        monkeypatch.setattr(cpp_backend, "_cpp", FakeCpp)

        # Call with mixed-type inputs
        pts = np.array([[0.0, 0.0, 9.0], [1.0, 0.0, 9.0], [0.0, 1.0, 9.0]])  # 3 cols
        conn = np.array([[0, 1, 2]], dtype=np.int64)  # int64

        result = cpp_backend.fast_init(pts, conn)

        # Verify result
        assert result == "sentinel_mesh"

        # Verify marshalled arrays
        forwarded_pts = FakeCpp.call_record['points']
        forwarded_conn = FakeCpp.call_record['connectivity']

        # Points should be sliced to [:, :2], float64, C-contiguous
        assert forwarded_pts.shape == (3, 2), f"Expected (3, 2), got {forwarded_pts.shape}"
        assert forwarded_pts.dtype == np.float64
        assert forwarded_pts.flags['C_CONTIGUOUS']
        assert not np.any(forwarded_pts == 9.0), "Third column (9.0) should not be present"

        # Connectivity should be int32, C-contiguous
        assert forwarded_conn.dtype == np.int32
        assert forwarded_conn.flags['C_CONTIGUOUS']
        assert forwarded_conn.shape == (1, 3)

    def test_full_init_marshalls_points_and_connectivity(self, monkeypatch):
        """full_init() converts points to float64 C-order and conn to int32 C-order."""
        class FakeCpp:
            call_record = {}

            @staticmethod
            def full_init(points, connectivity):
                FakeCpp.call_record['points'] = points
                FakeCpp.call_record['connectivity'] = connectivity
                return "sentinel_mesh"

        monkeypatch.setattr(cpp_backend, "CPP_AVAILABLE", True)
        monkeypatch.setattr(cpp_backend, "_cpp", FakeCpp)

        pts = np.array([[0.0, 0.0, 9.0], [1.0, 0.0, 9.0], [0.0, 1.0, 9.0]])
        conn = np.array([[0, 1, 2]], dtype=np.int64)

        result = cpp_backend.full_init(pts, conn)

        assert result == "sentinel_mesh"
        assert FakeCpp.call_record['points'].shape == (3, 2)
        assert FakeCpp.call_record['points'].dtype == np.float64
        assert FakeCpp.call_record['connectivity'].dtype == np.int32

    def test_quality_analysis_forwards_mesh_and_returns_result(self, monkeypatch):
        """quality_analysis() forwards mesh object and returns C++ result."""
        sentinel_result = np.array([0.5, 0.6, 0.7])

        class FakeCpp:
            call_record = {}

            @staticmethod
            def quality_analysis(mesh):
                FakeCpp.call_record['mesh'] = mesh
                return sentinel_result

        monkeypatch.setattr(cpp_backend, "CPP_AVAILABLE", True)
        monkeypatch.setattr(cpp_backend, "_cpp", FakeCpp)

        fake_mesh = object()
        result = cpp_backend.quality_analysis(fake_mesh)

        assert result is sentinel_result
        assert FakeCpp.call_record['mesh'] is fake_mesh

    def test_get_vertex_edges_forwards_mesh_and_vertex_id(self, monkeypatch):
        """get_vertex_edges() forwards mesh and vertex_id, returns C++ result."""
        sentinel_result = np.array([0, 1, 2], dtype=np.int32)

        class FakeCpp:
            call_record = {}

            @staticmethod
            def get_vertex_edges(mesh, vertex_id):
                FakeCpp.call_record['mesh'] = mesh
                FakeCpp.call_record['vertex_id'] = vertex_id
                return sentinel_result

        monkeypatch.setattr(cpp_backend, "CPP_AVAILABLE", True)
        monkeypatch.setattr(cpp_backend, "_cpp", FakeCpp)

        fake_mesh = object()
        result = cpp_backend.get_vertex_edges(fake_mesh, 42)

        assert result is sentinel_result
        assert FakeCpp.call_record['mesh'] is fake_mesh
        assert FakeCpp.call_record['vertex_id'] == 42


# ============================================================================
# RUST BACKEND TESTS
# ============================================================================

class TestRustBackendUnavailable:
    """Test Rust backend error paths when RUST_AVAILABLE is False.

    Relies on natural False (no extension installed) but guards with skip
    if the extension is present.
    """

    def test_require_rust_raises_when_unavailable(self):
        """_require_rust() raises ImportError with documented message."""
        if rust_backend.RUST_AVAILABLE:
            pytest.skip("Rust extension present on this machine")

        with pytest.raises(ImportError) as exc_info:
            rust_backend._require_rust()

        assert "chilmesh_core Rust extension not available" in str(exc_info.value)

    def test_full_init_raises_when_unavailable(self):
        """full_init() raises ImportError when Rust not available."""
        if rust_backend.RUST_AVAILABLE:
            pytest.skip("Rust extension present")

        pts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        conn = np.array([[0, 1, 2]], dtype=np.int64)

        with pytest.raises(ImportError) as exc_info:
            rust_backend.full_init(pts, conn)

        assert "chilmesh_core Rust extension not available" in str(exc_info.value)


class TestRustBackendAvailableViaMonkeypatch:
    """Test Rust backend with simulated available backend (monkeypatch).

    Uses fake backend module to verify marshalling and method-call logic.
    """

    def test_full_init_builds_mesh_and_calls_all_methods(self, monkeypatch):
        """full_init() creates RustMesh and calls set_points, set_connectivity,
        build_adjacencies, and skeletonize in the correct order.
        """
        call_order = []

        class FakeRustMesh:
            def set_points(self, pts):
                call_order.append(('set_points', pts))

            def set_connectivity(self, conn):
                call_order.append(('set_connectivity', conn))

            def build_adjacencies(self):
                call_order.append(('build_adjacencies',))

            def skeletonize(self):
                call_order.append(('skeletonize',))

        class FakeRust:
            @staticmethod
            def RustMesh():
                return FakeRustMesh()

        monkeypatch.setattr(rust_backend, "RUST_AVAILABLE", True)
        monkeypatch.setattr(rust_backend, "_rust", FakeRust)

        pts = np.array([[0.0, 0.0, 9.0], [1.0, 0.0, 9.0], [0.0, 1.0, 9.0]])
        conn = np.array([[0, 1, 2]], dtype=np.int64)

        result = rust_backend.full_init(pts, conn)

        # Verify result is the mesh object
        assert isinstance(result, FakeRustMesh)

        # Verify call order: set_points, set_connectivity, build_adjacencies, skeletonize
        assert len(call_order) == 4
        assert call_order[0][0] == 'set_points'
        assert call_order[1][0] == 'set_connectivity'
        assert call_order[2][0] == 'build_adjacencies'
        assert call_order[3][0] == 'skeletonize'

    def test_rust_full_init_marshalls_points_to_float64_c_order(self, monkeypatch):
        """Rust full_init() converts points to float64 C-order and slices to [:, :2]."""
        forwarded_pts = None

        class FakeRustMesh:
            def set_points(self, pts):
                nonlocal forwarded_pts
                forwarded_pts = pts

            def set_connectivity(self, conn):
                pass

            def build_adjacencies(self):
                pass

            def skeletonize(self):
                pass

        class FakeRust:
            @staticmethod
            def RustMesh():
                return FakeRustMesh()

        monkeypatch.setattr(rust_backend, "RUST_AVAILABLE", True)
        monkeypatch.setattr(rust_backend, "_rust", FakeRust)

        pts = np.array([[0.0, 0.0, 9.0], [1.0, 0.0, 9.0], [0.0, 1.0, 9.0]])
        conn = np.array([[0, 1, 2]], dtype=np.int64)

        rust_backend.full_init(pts, conn)

        # Verify marshalled points
        assert forwarded_pts.shape == (3, 2)
        assert forwarded_pts.dtype == np.float64
        assert forwarded_pts.flags['C_CONTIGUOUS']
        assert not np.any(forwarded_pts == 9.0)

    def test_rust_full_init_marshalls_connectivity_to_int32_c_order(self, monkeypatch):
        """Rust full_init() converts connectivity to int32 C-order."""
        forwarded_conn = None

        class FakeRustMesh:
            def set_points(self, pts):
                pass

            def set_connectivity(self, conn):
                nonlocal forwarded_conn
                forwarded_conn = conn

            def build_adjacencies(self):
                pass

            def skeletonize(self):
                pass

        class FakeRust:
            @staticmethod
            def RustMesh():
                return FakeRustMesh()

        monkeypatch.setattr(rust_backend, "RUST_AVAILABLE", True)
        monkeypatch.setattr(rust_backend, "_rust", FakeRust)

        pts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        conn = np.array([[0, 1, 2]], dtype=np.int64)  # int64 input

        rust_backend.full_init(pts, conn)

        # Verify marshalled connectivity
        assert forwarded_conn.dtype == np.int32
        assert forwarded_conn.flags['C_CONTIGUOUS']
