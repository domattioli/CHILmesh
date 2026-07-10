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


def convex_hull(points):
    """Andrew's monotone chain 2D convex hull algorithm.

    Parameters
    ----------
    points : array-like, shape (N, 2)
        Point cloud in (x, y) coordinates.

    Returns
    -------
    ndarray, shape (H, 2)
        Vertices of the convex hull in counter-clockwise order, starting at
        the lexicographically smallest point. Collinear points are excluded.

    Examples
    --------
    >>> import numpy as np
    >>> from chilmesh.geometry import convex_hull
    >>> pts = np.array([[0, 0], [1, 1], [1, 0], [0, 1], [0.5, 0.5]])
    >>> hull = convex_hull(pts)
    >>> hull.shape[0]
    4
    """
    pts = sorted(set(map(tuple, np.asarray(points, dtype=float).reshape(-1, 2))))
    if len(pts) <= 2:
        return np.asarray(pts, dtype=float)

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    # Build lower hull
    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    # Build upper hull
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    return np.asarray(lower[:-1] + upper[:-1], dtype=float)


def is_antimeridian_wrapping(bbox):
    """Test if a geographic bounding box wraps the antimeridian (date line).

    Parameters
    ----------
    bbox : tuple of float
        Bounding box as (min_lon, min_lat, max_lon, max_lat) in degrees.
        Longitude in range [-180, 180].

    Returns
    -------
    bool
        True if min_lon > max_lon (wrapping case), False otherwise.

    Examples
    --------
    >>> from chilmesh.geometry import is_antimeridian_wrapping
    >>> is_antimeridian_wrapping((170, -10, -170, 10))
    True
    >>> is_antimeridian_wrapping((-10, -10, 10, 10))
    False
    """
    min_lon, min_lat, max_lon, max_lat = bbox
    return min_lon > max_lon


def split_antimeridian_bbox(bbox):
    """Split a geographic bounding box that wraps the antimeridian.

    Parameters
    ----------
    bbox : tuple of float
        Bounding box as (min_lon, min_lat, max_lon, max_lat) in degrees.

    Returns
    -------
    list of tuple
        If wrapping (min_lon > max_lon), returns a 2-element list of
        non-wrapping bboxes: [(min_lon, min_lat, 180.0, max_lat),
        (-180.0, min_lat, max_lon, max_lat)]. Otherwise returns [bbox].

    Examples
    --------
    >>> from chilmesh.geometry import split_antimeridian_bbox
    >>> result = split_antimeridian_bbox((170, -10, -170, 10))
    >>> len(result)
    2
    >>> result[0]
    (170, -10, 180.0, 10)
    """
    min_lon, min_lat, max_lon, max_lat = bbox
    if min_lon > max_lon:
        return [
            (min_lon, min_lat, 180.0, max_lat),
            (-180.0, min_lat, max_lon, max_lat),
        ]
    return [bbox]


def bbox_iou(a, b):
    """Intersection-over-union of two geographic bounding boxes, antimeridian-aware.

    Parameters
    ----------
    a, b : tuple of float
        Bounding boxes as (min_lon, min_lat, max_lon, max_lat) in degrees.

    Returns
    -------
    float
        Intersection over union in [0, 1]. Returns 0.0 if bboxes do not
        overlap or union area is 0.

    Examples
    --------
    >>> from chilmesh.geometry import bbox_iou
    >>> bbox_iou((0, 0, 2, 2), (0, 0, 2, 2))  # identical
    1.0
    >>> bbox_iou((0, 0, 2, 2), (3, 3, 5, 5))  # disjoint
    0.0
    """
    def area(b):
        min_lon, min_lat, max_lon, max_lat = b
        return max(0, max_lon - min_lon) * max(0, max_lat - min_lat)

    def inter_area(b1, b2):
        min_lon1, min_lat1, max_lon1, max_lat1 = b1
        min_lon2, min_lat2, max_lon2, max_lat2 = b2
        inter_lon = max(0, min(max_lon1, max_lon2) - max(min_lon1, min_lon2))
        inter_lat = max(0, min(max_lat1, max_lat2) - max(min_lat1, min_lat2))
        return inter_lon * inter_lat

    # Split each bbox and compute all pairwise intersections
    parts_a = split_antimeridian_bbox(a)
    parts_b = split_antimeridian_bbox(b)

    total_inter = sum(inter_area(ba, bb) for ba in parts_a for bb in parts_b)
    area_a = sum(area(ba) for ba in parts_a)
    area_b = sum(area(bb) for bb in parts_b)
    union = area_a + area_b - total_inter

    if union <= 0:
        return 0.0
    return total_inter / union


def hausdorff_distance(pts_a, pts_b, crs: str = "cartesian"):
    """Discrete symmetric Hausdorff distance between two point sets.

    Parameters
    ----------
    pts_a, pts_b : array-like, shape (N, 2) or (M, 2)
        Point sets in (x, y) coordinates when crs="cartesian" or (lon, lat)
        degrees when crs="spherical".
    crs : {"cartesian", "spherical"}
        "cartesian": planar Euclidean distance via hypot. "spherical":
        great-circle distance via haversine_m.

    Returns
    -------
    float
        Symmetric Hausdorff distance: max(max_a, max_b) where max_a is the
        maximum distance from any point in pts_a to its nearest point in pts_b.

    Raises
    ------
    ValueError
        If crs is not "cartesian" or "spherical", or if either point set is empty.

    Examples
    --------
    >>> import numpy as np
    >>> from chilmesh.geometry import hausdorff_distance
    >>> pts_a = np.array([[0, 0], [1, 0]])
    >>> pts_b = np.array([[0, 0]])
    >>> hausdorff_distance(pts_a, pts_b, crs="cartesian")
    1.0
    """
    if crs not in ("cartesian", "spherical"):
        raise ValueError(f"unsupported crs={crs!r}: expected 'cartesian' or 'spherical'")

    a = np.atleast_2d(np.asarray(pts_a, dtype=float))
    b = np.atleast_2d(np.asarray(pts_b, dtype=float))

    if a.size == 0 or b.size == 0:
        raise ValueError("point sets must be non-empty")

    if crs == "cartesian":
        # Broadcast: (N, 1, 2) vs (1, M, 2) -> (N, M)
        dists = np.hypot(
            a[:, None, 0] - b[None, :, 0], a[:, None, 1] - b[None, :, 1]
        )
    else:  # spherical
        dists = haversine_m(
            a[:, None, 0], a[:, None, 1], b[None, :, 0], b[None, :, 1]
        )

    # max over all min-to-B distances in A
    max_a = np.max(np.min(dists, axis=1))
    # max over all min-to-A distances in B
    max_b = np.max(np.min(dists, axis=0))

    return float(np.max([max_a, max_b]))


__all__ = [
    "EARTH_RADIUS_M",
    "haversine_m",
    "edge_lengths",
    "convex_hull",
    "is_antimeridian_wrapping",
    "split_antimeridian_bbox",
    "bbox_iou",
    "hausdorff_distance",
]
