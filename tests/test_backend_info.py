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
