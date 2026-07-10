"""Geodesic and planar geometry helpers (CRS-aware edge length).

Raw-array functions independent of mesh object construction, mirroring the
style of :mod:`chilmesh.quality`. The CRS is caller-declared ("cartesian"
projected planar units vs "spherical" lon/lat degrees) — for ADCIRC inputs it
can be resolved from fort.15 ICS (1=cartesian, 2=spherical), see
:class:`chilmesh.fort15_io.Fort15`.
"""
from __future__ import annotations

import numpy as np

EARTH_RADIUS_M = 6_371_000.0  # mean Earth radius for the haversine metric


def haversine_m(lon1, lat1, lon2, lat2):
    """Great-circle distance in metres between lon/lat points (degrees).

    Accepts scalars or broadcastable numpy arrays; returns float for scalar
    input, ndarray otherwise.
    """
    lon1, lat1, lon2, lat2 = (np.asarray(a, dtype=float) for a in (lon1, lat1, lon2, lat2))
    p1 = np.radians(lat1)
    p2 = np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlam = np.radians(lon2 - lon1)
    h = np.sin(dphi / 2.0) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlam / 2.0) ** 2
    d = 2.0 * EARTH_RADIUS_M * np.arcsin(np.minimum(1.0, np.sqrt(h)))
    return float(d) if d.ndim == 0 else d


def edge_lengths(p1, p2, crs: str = "cartesian"):
    """Distance between paired points, honouring the declared CRS.

    Parameters
    ----------
    p1, p2 : array-like, shape (2,) or (N, 2)
        Point coordinates: (x, y) when crs="cartesian", (lon, lat) degrees
        when crs="spherical".
    crs : {"cartesian", "spherical"}
        "cartesian": planar Euclidean in native units. "spherical": haversine
        great-circle metres.

    Returns
    -------
    float or ndarray of shape (N,)
    """
    a = np.atleast_2d(np.asarray(p1, dtype=float))
    b = np.atleast_2d(np.asarray(p2, dtype=float))
    if a.shape != b.shape or a.shape[-1] < 2:
        raise ValueError(f"p1/p2 must be matching (N, 2) arrays, got {a.shape} vs {b.shape}")
    if crs == "cartesian":
        d = np.hypot(b[:, 0] - a[:, 0], b[:, 1] - a[:, 1])
    elif crs == "spherical":
        d = haversine_m(a[:, 0], a[:, 1], b[:, 0], b[:, 1])
    else:
        raise ValueError(f"unsupported crs={crs!r}: expected 'cartesian' or 'spherical'")
    d = np.atleast_1d(d)
    return float(d[0]) if np.asarray(p1).ndim == 1 else d


__all__ = ["EARTH_RADIUS_M", "haversine_m", "edge_lengths"]
