"""C++ half-edge backend for CHILmesh.

This module wraps the compiled chilmesh_cpp extension (pybind11/C++17).
It provides the same logical interface as the Rust backend wrapper.

Usage::

    from chilmesh.backends.cpp_backend import CPP_AVAILABLE, fast_init, full_init
    if CPP_AVAILABLE:
        mesh = full_init(points, connectivity)
"""
from __future__ import annotations

import numpy as np
from typing import Any

try:
    import chilmesh_cpp as _cpp
    CPP_AVAILABLE = True
except ImportError:
    CPP_AVAILABLE = False
    _cpp = None  # type: ignore


def _require_cpp() -> None:
    if not CPP_AVAILABLE:
        raise ImportError(
            "chilmesh_cpp C++ extension not available. "
            "Build it with: cd src/chilmesh_cpp && pip install -e . --no-build-isolation"
        )


def fast_init(points: np.ndarray, connectivity: np.ndarray) -> Any:
    """Build half-edge mesh without adjacency or skeletonization.

    Parameters
    ----------
    points : ndarray [n_verts, 2], float64
        Vertex x, y coordinates.
    connectivity : ndarray [n_elems, 3|4], int32
        Element vertex indices. Use -1 padding for triangles in a quad mesh.

    Returns
    -------
    CppMesh
        C++ half-edge mesh object.
    """
    _require_cpp()
    pts = np.asarray(points[:, :2], dtype=np.float64, order='C')
    conn = np.asarray(connectivity, dtype=np.int32, order='C')
    return _cpp.fast_init(pts, conn)


def full_init(points: np.ndarray, connectivity: np.ndarray) -> Any:
    """Build half-edge mesh with adjacency and skeletonization.

    Parameters
    ----------
    points : ndarray [n_verts, 2], float64
        Vertex x, y coordinates.
    connectivity : ndarray [n_elems, 3|4], int32
        Element vertex indices.

    Returns
    -------
    CppMesh
        C++ half-edge mesh with edge2vert, elem2edge, vert2edge, and layers.
    """
    _require_cpp()
    pts = np.asarray(points[:, :2], dtype=np.float64, order='C')
    conn = np.asarray(connectivity, dtype=np.int32, order='C')
    return _cpp.full_init(pts, conn)


def quality_analysis(mesh: Any) -> np.ndarray:
    """Compute signed area for each element.

    Parameters
    ----------
    mesh : CppMesh
        Mesh returned by fast_init or full_init.

    Returns
    -------
    ndarray [n_elems], float64
        Signed area per element (positive = CCW).
    """
    _require_cpp()
    return _cpp.quality_analysis(mesh)


def get_vertex_edges(mesh: Any, vertex_id: int) -> np.ndarray:
    """Return edge IDs incident to a vertex.

    Parameters
    ----------
    mesh : CppMesh
        Mesh returned by full_init (adjacency required).
    vertex_id : int
        0-indexed vertex ID.

    Returns
    -------
    ndarray of int32
        Edge IDs incident to the given vertex.
    """
    _require_cpp()
    return _cpp.get_vertex_edges(mesh, vertex_id)


__all__ = [
    "CPP_AVAILABLE",
    "fast_init",
    "full_init",
    "quality_analysis",
    "get_vertex_edges",
]
