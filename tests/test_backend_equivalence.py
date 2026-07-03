"""Cross-backend equivalence tests for CHILmesh.

These tests ensure that the Python, Rust, and C++ backends produce
bit-identical mesh topology, skeletonization layers, and quality outputs
across all built-in fixtures plus the WNAT_Hagen reference mesh.

If a backend is unavailable in the test environment, the corresponding
tests are skipped (not failed) so users can run the suite with only
some backends compiled.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from chilmesh import CHILmesh, examples

try:
    import chilmesh_cpp
    # Namespace-package import can succeed without compiled functions;
    # require the actual API surface to consider C++ "available".
    CPP_AVAILABLE = hasattr(chilmesh_cpp, "fast_init") and hasattr(chilmesh_cpp, "full_init")
except ImportError:
    CPP_AVAILABLE = False

try:
    import chilmesh_core
    RUST_AVAILABLE = hasattr(chilmesh_core, "RustMesh")
except ImportError:
    RUST_AVAILABLE = False


FIXTURE_NAMES = ["annulus", "donut", "block_o", "structured"]


def _load_python_mesh(fixture_name: str) -> CHILmesh:
    """Load a built-in fixture as a CHILmesh."""
    if fixture_name == "wnat_hagen":
        path = Path("/tmp/WNAT_Hagen.14")
        if not path.exists():
            pytest.skip("WNAT_Hagen.14 not available in this environment")
        return CHILmesh.read_from_fort14(path, compute_layers=True)
    return getattr(examples, fixture_name)()


def _arrays_for_cpp(mesh: CHILmesh):
    """Extract points (2D) and connectivity arrays from a Python mesh
    in the format expected by chilmesh_cpp.fast_init / full_init."""
    points = mesh.points[:, :2].astype(np.float64)
    conn = np.array(
        [list(c[:3]) for c in mesh.connectivity_list],
        dtype=np.int32,
    )
    return points, conn


@pytest.mark.skipif(not CPP_AVAILABLE, reason="chilmesh_cpp not built")
@pytest.mark.parametrize("fixture", FIXTURE_NAMES)
def test_cpp_layer_count_matches_python(fixture):
    """C++ skeletonization must produce same number of layers as Python."""
    mesh = _load_python_mesh(fixture)
    points, conn = _arrays_for_cpp(mesh)
    cpp = chilmesh_cpp.full_init(points, conn)
    assert mesh.n_layers == cpp.n_layers, (
        f"layer count mismatch on {fixture}: "
        f"Python={mesh.n_layers}, C++={cpp.n_layers}"
    )


@pytest.mark.skipif(not CPP_AVAILABLE, reason="chilmesh_cpp not built")
@pytest.mark.parametrize("fixture", FIXTURE_NAMES)
@pytest.mark.parametrize("layer_key", ["OE", "IE", "OV", "IV", "bEdgeIDs"])
def test_cpp_layer_member_sets_match_python(fixture, layer_key):
    """For each layer, the set of OE/IE/OV/IV/bEdgeIDs must match Python."""
    mesh = _load_python_mesh(fixture)
    points, conn = _arrays_for_cpp(mesh)
    cpp = chilmesh_cpp.full_init(points, conn)

    assert mesh.n_layers == cpp.n_layers

    for i in range(mesh.n_layers):
        py_set = set(mesh.Layers[layer_key][i].tolist())
        cpp_set = set(cpp.layers[i][layer_key].tolist())
        assert py_set == cpp_set, (
            f"{fixture} layer {i} {layer_key} mismatch: "
            f"Python {len(py_set)}, C++ {len(cpp_set)}, "
            f"symmetric diff size {len(py_set ^ cpp_set)}"
        )


@pytest.mark.skipif(not CPP_AVAILABLE, reason="chilmesh_cpp not built")
@pytest.mark.parametrize("fixture", FIXTURE_NAMES)
def test_cpp_signed_area_matches_python(fixture):
    """quality_analysis (signed areas) must match Python to float64 tolerance."""
    mesh = _load_python_mesh(fixture)
    points, conn = _arrays_for_cpp(mesh)
    cpp = chilmesh_cpp.full_init(points, conn)
    cpp_areas = np.asarray(chilmesh_cpp.quality_analysis(cpp))

    # Python's signed_area is on CHILmesh
    py_areas = mesh.signed_area(np.arange(mesh.n_elems, dtype=np.intp))
    assert py_areas.shape == cpp_areas.shape
    np.testing.assert_allclose(
        py_areas, cpp_areas, rtol=1e-9, atol=1e-12,
        err_msg=f"signed areas diverge on {fixture}",
    )


@pytest.mark.skipif(not CPP_AVAILABLE, reason="chilmesh_cpp not built")
@pytest.mark.parametrize("fixture", FIXTURE_NAMES)
def test_cpp_vert2edge_consistent_with_python(fixture):
    """Vert2Edge edge ID sets per vertex must match Python."""
    mesh = _load_python_mesh(fixture)
    points, conn = _arrays_for_cpp(mesh)
    cpp = chilmesh_cpp.full_init(points, conn)

    py_v2e = mesh.adjacencies["Vert2Edge"]
    cpp_v2e = cpp.vert2edge

    n_check = min(mesh.n_verts, 200)  # spot-check first 200 verts
    for v in range(n_check):
        py_edges = set(py_v2e[v].tolist()) if hasattr(py_v2e[v], 'tolist') else set(py_v2e[v])
        cpp_edges = set(cpp_v2e[v])
        assert py_edges == cpp_edges, (
            f"{fixture} vert {v}: Vert2Edge mismatch "
            f"Python={py_edges}, C++={cpp_edges}"
        )


def test_backend_info_reports_available():
    """backend_info() returns a structured dict with expected keys."""
    from chilmesh import backend_info
    info = backend_info()
    assert "available" in info
    assert "selected" in info
    assert "versions" in info
    assert "python" in info["available"]
    assert info["selected"] in info["available"]


def test_backend_info_honors_env_override(monkeypatch):
    """CHILMESH_BACKEND env var forces a specific backend selection."""
    from chilmesh import backend_info
    monkeypatch.setenv("CHILMESH_BACKEND", "python")
    info = backend_info()
    assert info["selected"] == "python"


@pytest.mark.skipif(not CPP_AVAILABLE, reason="chilmesh_cpp not built")
def test_backend_info_reports_cpp_when_available():
    """When chilmesh_cpp is importable, backend_info lists 'cpp'."""
    from chilmesh import backend_info
    info = backend_info()
    assert "cpp" in info["available"]


def _build_rust_mesh(mesh):
    """Build a chilmesh_core.RustMesh from a Python CHILmesh via the array API."""
    points, conn = _arrays_for_cpp(mesh)
    r = chilmesh_core.RustMesh()
    r.set_points(points)
    r.set_connectivity(conn)
    r.build_adjacencies()
    r.skeletonize()
    return r


# NOTE: Rust parity is asserted on n_layers, the OE/IE/OV/IV layer member sets,
# per-layer bEdgeIDs, signed areas, full-mesh edge2vert ordering, and Vert2Edge
# (CHILmesh #163 — Rust edge IDs were aligned to Python's first-encounter
# ordering, and bEdgeIDs are now recorded as full-mesh edge IDs to match Python).
# Element-id (OE/IE), vertex-id (OV/IV), and full-mesh edge-id (bEdgeIDs)
# membership are stable identities and DO match.
@pytest.mark.skipif(not RUST_AVAILABLE, reason="chilmesh_core (Rust) not built")
@pytest.mark.parametrize("fixture", FIXTURE_NAMES)
def test_rust_layer_count_matches_python(fixture):
    """Rust skeletonization must produce same number of layers as Python."""
    mesh = _load_python_mesh(fixture)
    rust = _build_rust_mesh(mesh)
    assert mesh.n_layers == rust.get_num_layers(), (
        f"layer count mismatch on {fixture}: "
        f"Python={mesh.n_layers}, Rust={rust.get_num_layers()}"
    )


@pytest.mark.skipif(not RUST_AVAILABLE, reason="chilmesh_core (Rust) not built")
@pytest.mark.parametrize("fixture", FIXTURE_NAMES)
@pytest.mark.parametrize("layer_key", ["OE", "IE", "OV", "IV", "bEdgeIDs"])
def test_rust_layer_member_sets_match_python(fixture, layer_key):
    """For each layer, the OE/IE/OV/IV/bEdgeIDs member sets must match Python.

    bEdgeIDs are full-mesh edge IDs matching Python (#163).
    """
    mesh = _load_python_mesh(fixture)
    rust = _build_rust_mesh(mesh)
    assert mesh.n_layers == rust.get_num_layers()
    layers = rust.get_all_layers()
    for i in range(mesh.n_layers):
        py_set = set(mesh.Layers[layer_key][i].tolist())
        rust_set = set(layers[i][layer_key])
        assert py_set == rust_set, (
            f"{fixture} layer {i} {layer_key} mismatch: "
            f"Python {len(py_set)}, Rust {len(rust_set)}, "
            f"symmetric diff size {len(py_set ^ rust_set)}"
        )


@pytest.mark.skipif(not RUST_AVAILABLE, reason="chilmesh_core (Rust) not built")
@pytest.mark.parametrize("fixture", FIXTURE_NAMES)
def test_rust_bedge_ids_order_matches_python(fixture):
    """Per-layer bEdgeIDs must match Python as ascending full-mesh edge-ID
    arrays, not merely as sets (CHILmesh #163: Rust records full-mesh edge IDs)."""
    mesh = _load_python_mesh(fixture)
    rust = _build_rust_mesh(mesh)
    layers = rust.get_all_layers()
    assert mesh.n_layers == rust.get_num_layers()
    for i in range(mesh.n_layers):
        py_ids = sorted(int(x) for x in mesh.Layers["bEdgeIDs"][i].tolist())
        rust_ids = sorted(int(x) for x in layers[i]["bEdgeIDs"])
        assert py_ids == rust_ids, (
            f"{fixture} layer {i} bEdgeIDs mismatch: "
            f"Python {len(py_ids)} ids, Rust {len(rust_ids)} ids"
        )


@pytest.mark.skipif(not RUST_AVAILABLE, reason="chilmesh_core (Rust) not built")
@pytest.mark.parametrize("fixture", FIXTURE_NAMES)
def test_rust_signed_area_matches_python(fixture):
    """quality_analysis (signed areas) must match Python to float64 tolerance."""
    mesh = _load_python_mesh(fixture)
    rust = _build_rust_mesh(mesh)
    rust.compute_quality()
    rust_areas = np.asarray(rust.get_signed_areas())
    py_areas = mesh.signed_area(np.arange(mesh.n_elems, dtype=np.intp))
    assert py_areas.shape == rust_areas.shape
    np.testing.assert_allclose(
        py_areas, rust_areas, rtol=1e-9, atol=1e-12,
        err_msg=f"signed areas diverge on {fixture}",
    )


@pytest.mark.skipif(not RUST_AVAILABLE, reason="chilmesh_core (Rust) not built")
@pytest.mark.parametrize("fixture", FIXTURE_NAMES)
def test_rust_edge2vert_order_matches_python(fixture):
    """Full-mesh Edge2Vert must match Python in both content AND ordering.

    Edge IDs are the index into this array; matching order means Rust edge IDs
    equal Python's, so everything keyed on them (Vert2Edge) is comparable (#163).
    """
    mesh = _load_python_mesh(fixture)
    rust = _build_rust_mesh(mesh)
    py_e2v = np.asarray(mesh.adjacencies["Edge2Vert"])
    rust_e2v = np.asarray(rust.get_edge2vert())
    assert py_e2v.shape == rust_e2v.shape, (
        f"{fixture}: Edge2Vert shape mismatch "
        f"Python={py_e2v.shape}, Rust={rust_e2v.shape}"
    )
    assert np.array_equal(py_e2v, rust_e2v), (
        f"{fixture}: Edge2Vert ordering diverges between Python and Rust"
    )


@pytest.mark.skipif(not RUST_AVAILABLE, reason="chilmesh_core (Rust) not built")
@pytest.mark.parametrize("fixture", FIXTURE_NAMES)
def test_rust_vert2edge_consistent_with_python(fixture):
    """Vert2Edge edge-ID sets per vertex must match Python (#163)."""
    mesh = _load_python_mesh(fixture)
    rust = _build_rust_mesh(mesh)
    py_v2e = mesh.adjacencies["Vert2Edge"]
    rust_v2e = rust.get_vert2edge()
    n_check = min(mesh.n_verts, 200)  # spot-check first 200 verts
    for v in range(n_check):
        py_edges = set(py_v2e[v].tolist()) if hasattr(py_v2e[v], "tolist") else set(py_v2e[v])
        rust_edges = set(rust_v2e[v])
        assert py_edges == rust_edges, (
            f"{fixture} vert {v}: Vert2Edge mismatch "
            f"Python={py_edges}, Rust={rust_edges}"
        )


@pytest.mark.skipif(not RUST_AVAILABLE, reason="chilmesh_core (Rust) not built")
def test_backend_info_reports_rust_when_available():
    """When chilmesh_core is importable, backend_info lists 'rust'."""
    from chilmesh import backend_info
    info = backend_info()
    assert "rust" in info["available"]
