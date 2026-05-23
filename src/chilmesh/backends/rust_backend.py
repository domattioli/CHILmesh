"""Rust quad-edge backend for CHILmesh.

This module wraps the compiled chilmesh_core extension (PyO3/maturin).
It provides the same logical interface as the C++ backend wrapper.

Usage::

    from chilmesh.backends.rust_backend import RUST_AVAILABLE, full_init
    if RUST_AVAILABLE:
        mesh = full_init(points, connectivity)
"""
from __future__ import annotations

import numpy as np
from typing import Any

try:
    import chilmesh_core as _rust
    # Guard against an importable-but-empty namespace package (CHILmesh #163).
    RUST_AVAILABLE = hasattr(_rust, "RustMesh")
    if not RUST_AVAILABLE:
        _rust = None  # type: ignore
except ImportError:
    RUST_AVAILABLE = False
    _rust = None  # type: ignore


def _require_rust() -> None:
    if not RUST_AVAILABLE:
        raise ImportError(
            "chilmesh_core Rust extension not available. "
            "Build it with: maturin develop --release -C src/chilmesh_core"
        )


def full_init(points: np.ndarray, connectivity: np.ndarray) -> Any:
    """Build quad-edge mesh with full adjacency and skeletonization.

    Parameters
    ----------
    points : ndarray [n_verts, 2|3], float64
        Vertex coordinates (x, y, optional z).
    connectivity : ndarray [n_elems, 3|4], int32|int64
        Element vertex indices. Use -1 padding for triangles in a quad mesh.

    Returns
    -------
    RustMesh
        Rust quad-edge mesh object with adjacency and skeletonization.
    """
    _require_rust()
    pts = np.asarray(points[:, :2], dtype=np.float64, order='C')
    conn = np.asarray(connectivity, dtype=np.int32, order='C')
    mesh = _rust.RustMesh()
    mesh.set_points(pts)
    mesh.set_connectivity(conn)
    mesh.build_adjacencies()
    mesh.skeletonize()
    return mesh
