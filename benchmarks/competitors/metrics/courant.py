"""Courant number and the largest stable timestep from raw arrays.

``C`` is linear in ``dt`` (``C = sqrt(g h) dt / dx``), so
``dt_max = courant_max / max C(dt=1)``. The per-edge computation is
``chilmesh.quality.courant_number``; this module only derives the edge list from
connectivity and the ``dt_max`` helper (plan R7).
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from chilmesh.quality import courant_number

from .quality import _edge_table


def edges_from_connectivity(connectivity) -> np.ndarray:
    """Unique element-side edges, shape (E, 2). Quad diagonals are not edges."""
    edges, _ = _edge_table(connectivity)
    return edges


def _max_c(points, connectivity, depths, dt, crs) -> Optional[float]:
    edges = edges_from_connectivity(connectivity)
    c = courant_number(np.asarray(points, dtype=float)[:, :2], edges, np.asarray(depths, dtype=float), dt, crs=crs)
    if c.size == 0 or np.all(np.isnan(c)):
        return None
    return float(np.nanmax(c))


def max_courant(points, connectivity, depths, dt: float, crs: str = "cartesian") -> Optional[float]:
    """Maximum Courant number over wet edges at ``dt``; None if every edge is dry."""
    return _max_c(points, connectivity, depths, dt, crs)


def dt_max(points, connectivity, depths, courant_max: float = 1.0, crs: str = "cartesian") -> Optional[float]:
    """Largest stable timestep (s) for ``courant_max``; None if every edge is dry."""
    c1 = _max_c(points, connectivity, depths, 1.0, crs)
    if c1 is None or c1 <= 0:
        return None
    return courant_max / c1
