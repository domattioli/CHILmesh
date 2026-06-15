"""Tests for chilmesh.backend_info() — backend detection correctness."""
from __future__ import annotations

import chilmesh
from chilmesh.backends.cpp_backend import CPP_AVAILABLE
from chilmesh.backends.rust_backend import RUST_AVAILABLE


def test_backend_info_structure():
    """Test that backend_info() returns correct dict structure."""
    info = chilmesh.backend_info()
    assert isinstance(info, dict)
    assert "available" in info
    assert "selected" in info
    assert "versions" in info
    assert isinstance(info["available"], list)
    assert isinstance(info["selected"], str)
    assert isinstance(info["versions"], dict)


def test_python_always_available():
    """Test that 'python' backend is always in available list."""
    info = chilmesh.backend_info()
    assert "python" in info["available"]
    assert info["available"][-1] == "python"


def test_backend_availability_consistency():
    """Test that available backends match their *_AVAILABLE flags.

    If a backend is listed in available, its *_AVAILABLE flag must be True.
    If a backend's *_AVAILABLE flag is True, it must be in available.
    """
    info = chilmesh.backend_info()
    available = info["available"]

    # Check cpp consistency
    if "cpp" in available:
        assert CPP_AVAILABLE, "cpp listed in available but CPP_AVAILABLE is False"
    else:
        assert not CPP_AVAILABLE, "cpp not in available but CPP_AVAILABLE is True"

    # Check rust consistency
    if "rust" in available:
        assert RUST_AVAILABLE, "rust listed in available but RUST_AVAILABLE is False"
    else:
        assert not RUST_AVAILABLE, "rust not in available but RUST_AVAILABLE is True"


def test_selected_in_available():
    """Test that selected backend is in the available list."""
    info = chilmesh.backend_info()
    assert info["selected"] in info["available"]


def test_versions_match_available():
    """Test that every version entry corresponds to an available backend."""
    info = chilmesh.backend_info()
    available = set(info["available"])
    versions_keys = set(info["versions"].keys())
    assert versions_keys == available, (
        f"version keys {versions_keys} do not match available {available}"
    )


def test_cpp_priority_when_available():
    """Test that cpp is prioritized (first in list) when available.

    If cpp is available, it should be the first item in available list,
    and selected should be cpp (unless CHILMESH_BACKEND env override is set).
    """
    info = chilmesh.backend_info()
    if CPP_AVAILABLE:
        assert info["available"][0] == "cpp"


def test_rust_placement_when_available():
    """Test that rust is placed before python but after cpp.

    When rust is available:
    - If cpp also available: order is [cpp, rust, python]
    - If cpp not available: order is [rust, python]
    """
    info = chilmesh.backend_info()
    if RUST_AVAILABLE:
        assert info["available"][-1] == "python", "python should always be last"
        rust_idx = info["available"].index("rust")
        python_idx = info["available"].index("python")
        assert rust_idx < python_idx, "rust should come before python"
        if CPP_AVAILABLE:
            cpp_idx = info["available"].index("cpp")
            assert cpp_idx < rust_idx, "cpp should come before rust"


import importlib
import sys
import types

import numpy as np
import pytest


def test_namespace_stub_reports_cpp_unavailable():
    """#163/#202: importable-but-API-less chilmesh_cpp must report unavailable."""
    import chilmesh.backends.cpp_backend as cpp_backend
    stub = types.ModuleType("chilmesh_cpp")  # no full_init attribute
    saved = sys.modules.get("chilmesh_cpp")
    sys.modules["chilmesh_cpp"] = stub
    try:
        reloaded = importlib.reload(cpp_backend)
        assert reloaded.CPP_AVAILABLE is False
        assert reloaded._cpp is None
    finally:
        if saved is not None:
            sys.modules["chilmesh_cpp"] = saved
        else:
            sys.modules.pop("chilmesh_cpp", None)
        importlib.reload(cpp_backend)  # restore real detection


def test_namespace_stub_reports_rust_unavailable():
    """#163/#202: importable-but-API-less chilmesh_core must report unavailable."""
    import chilmesh.backends.rust_backend as rust_backend
    stub = types.ModuleType("chilmesh_core")  # no RustMesh attribute
    saved = sys.modules.get("chilmesh_core")
    sys.modules["chilmesh_core"] = stub
    try:
        reloaded = importlib.reload(rust_backend)
        assert reloaded.RUST_AVAILABLE is False
        assert reloaded._rust is None
    finally:
        if saved is not None:
            sys.modules["chilmesh_core"] = saved
        else:
            sys.modules.pop("chilmesh_core", None)
        importlib.reload(rust_backend)  # restore real detection


def test_stub_with_api_reports_cpp_available():
    """Guard is not over-eager: a module exposing full_init reports available."""
    import chilmesh.backends.cpp_backend as cpp_backend
    stub = types.ModuleType("chilmesh_cpp")
    stub.full_init = lambda *a, **k: None
    saved = sys.modules.get("chilmesh_cpp")
    sys.modules["chilmesh_cpp"] = stub
    try:
        reloaded = importlib.reload(cpp_backend)
        assert reloaded.CPP_AVAILABLE is True
    finally:
        if saved is not None:
            sys.modules["chilmesh_cpp"] = saved
        else:
            sys.modules.pop("chilmesh_cpp", None)
        importlib.reload(cpp_backend)  # restore real detection


def test_slow_path_warning_when_pure_python(monkeypatch):
    """#202: pure-Python skeletonization of a large mesh warns once (non-silent)."""
    if CPP_AVAILABLE or RUST_AVAILABLE:
        pytest.skip("compiled fast backend present; slow-path warning N/A")
    import chilmesh.CHILmesh  # ensure submodule imported
    cm = sys.modules["chilmesh.CHILmesh"]
    monkeypatch.setattr(cm, "_SLOW_PATH_ELEM_THRESHOLD", 1)
    monkeypatch.setattr(cm, "_SLOW_PATH_WARNED", False)
    pts = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    conn = np.array([[0, 1, 2], [0, 2, 3]])
    with pytest.warns(UserWarning, match="pure-Python backend"):
        cm.CHILmesh(connectivity=conn, points=pts, compute_layers=True)


def test_no_slow_path_warning_below_threshold(monkeypatch, recwarn):
    """Small meshes never trigger the slow-path warning."""
    if CPP_AVAILABLE or RUST_AVAILABLE:
        pytest.skip("compiled fast backend present; slow-path warning N/A")
    import chilmesh.CHILmesh  # ensure submodule imported
    cm = sys.modules["chilmesh.CHILmesh"]
    monkeypatch.setattr(cm, "_SLOW_PATH_ELEM_THRESHOLD", 10_000)
    monkeypatch.setattr(cm, "_SLOW_PATH_WARNED", False)
    pts = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    conn = np.array([[0, 1, 2], [0, 2, 3]])
    cm.CHILmesh(connectivity=conn, points=pts, compute_layers=True)
    slow = [w for w in recwarn.list if "pure-Python backend" in str(w.message)]
    assert slow == []
