"""Shared element-quality metrics computed from raw ``points`` and ``connectivity``.

Operating on plain arrays keeps the comparison fair: each tool's smoothed mesh is
scored by this one implementation, never by the tool's own metric.

Element quality is the minimum scaled Jacobian over the element's corners.
For a corner with edge vectors ``a`` (to next vertex) and ``b`` (to previous
vertex) of a counter-clockwise element, ``sin(theta) = (a x b) / (|a||b|)``.
Quads use that value directly (1 for a rectangle). Triangles are scaled by
``2/sqrt(3)`` so an equilateral triangle scores 1.

Connectivity may have 3 columns (triangles) or 4 columns. In a 4-column array a
row whose fourth vertex repeats an earlier vertex is a padded triangle.
"""
from __future__ import annotations

import numpy as np

BAR = 0.75
_TRI_SCALE = 2.0 / np.sqrt(3.0)


def split_elements(connectivity) -> tuple[np.ndarray, np.ndarray]:
    """Return (triangles (m,3), quads (k,4)); padded triangles go to the triangle array."""
    conn = np.asarray(connectivity, dtype=int)
    if conn.ndim != 2 or conn.shape[1] not in (3, 4):
        raise ValueError(f"connectivity must be (n, 3) or (n, 4), got {conn.shape}")
    if conn.shape[1] == 3:
        return conn, np.empty((0, 4), dtype=int)
    padded = (conn[:, 3] == conn[:, 0]) | (conn[:, 3] == conn[:, 1]) | (conn[:, 3] == conn[:, 2])
    return conn[padded, :3], conn[~padded]


def _corner_terms(xy: np.ndarray):
    """(sin, angle_deg) per corner for polygons ``xy`` of shape (n, k, 2), CCW-normalised."""
    n, k, _ = xy.shape
    nxt = np.roll(xy, -1, axis=1) - xy
    prv = np.roll(xy, 1, axis=1) - xy
    cross = nxt[..., 0] * prv[..., 1] - nxt[..., 1] * prv[..., 0]
    area2 = np.sum(xy[..., 0] * np.roll(xy[..., 1], -1, axis=1) - np.roll(xy[..., 0], -1, axis=1) * xy[..., 1], axis=1)
    cross = np.where(area2[:, None] < 0, -cross, cross)
    lens = np.linalg.norm(nxt, axis=2) * np.linalg.norm(prv, axis=2)
    dot = np.sum(nxt * prv, axis=2)
    with np.errstate(divide="ignore", invalid="ignore"):
        sin = np.where(lens > 0, cross / lens, 0.0)
        ang = np.degrees(np.arctan2(cross, dot))
    return sin, ang


def _polygon_arrays(points, connectivity):
    pts = np.asarray(points, dtype=float)[:, :2]
    tris, quads = split_elements(connectivity)
    return pts, tris, quads


def element_quality(points, connectivity) -> np.ndarray:
    """Per-element quality in input order of (triangles, then quads)."""
    pts, tris, quads = _polygon_arrays(points, connectivity)
    parts = []
    if len(tris):
        sin, _ = _corner_terms(pts[tris])
        parts.append(sin.min(axis=1) * _TRI_SCALE)
    if len(quads):
        sin, _ = _corner_terms(pts[quads])
        parts.append(sin.min(axis=1))
    return np.concatenate(parts) if parts else np.empty(0)


def interior_angles_deg(points, connectivity) -> np.ndarray:
    pts, tris, quads = _polygon_arrays(points, connectivity)
    parts = []
    for polys in (tris, quads):
        if len(polys):
            parts.append(_corner_terms(pts[polys])[1].ravel())
    return np.concatenate(parts) if parts else np.empty(0)


def _edge_table(connectivity):
    """Unique undirected edges and how many elements share each."""
    tris, quads = split_elements(connectivity)
    pairs = []
    for polys in (tris, quads):
        if len(polys):
            pairs.append(np.stack([polys, np.roll(polys, -1, axis=1)], axis=2).reshape(-1, 2))
    allp = np.sort(np.concatenate(pairs), axis=1) if pairs else np.empty((0, 2), dtype=int)
    edges, counts = np.unique(allp, axis=0, return_counts=True)
    return edges, counts


def boundary_nodes(connectivity) -> np.ndarray:
    """Sorted node ids on edges used by exactly one element."""
    edges, counts = _edge_table(connectivity)
    return np.unique(edges[counts == 1])


def pct_irregular_valence(points, connectivity):
    """Percent of interior nodes whose valence differs from the ideal; None if no interior nodes.

    Ideal valence is 6 when every incident element is a triangle, 4 when every
    incident element is a quad, and 5 for a mix (midpoint heuristic).
    """
    edges, counts = _edge_table(connectivity)
    n = int(np.asarray(connectivity).max()) + 1
    valence = np.bincount(edges.ravel(), minlength=n)
    tris, quads = split_elements(connectivity)
    n_tri_inc = np.bincount(tris.ravel(), minlength=n) if len(tris) else np.zeros(n, dtype=int)
    n_quad_inc = np.bincount(quads.ravel(), minlength=n) if len(quads) else np.zeros(n, dtype=int)
    boundary = np.zeros(n, dtype=bool)
    boundary[boundary_nodes(connectivity)] = True
    used = (n_tri_inc + n_quad_inc) > 0
    interior = used & ~boundary
    if not interior.any():
        return None
    ideal = np.where(n_quad_inc == 0, 6, np.where(n_tri_inc == 0, 4, 5))
    return float(100.0 * np.mean(valence[interior] != ideal[interior]))


def quality_summary(points, connectivity) -> dict:
    q = element_quality(points, connectivity)
    ang = interior_angles_deg(points, connectivity)
    mean, sigma = float(q.mean()), float(q.std())
    return {
        "n_elems": int(len(q)),
        "q_min": float(q.min()),
        "q_mean": mean,
        "q_mean_minus_3sigma": mean - 3 * sigma,
        "meets_0_75_bar": bool(mean - 3 * sigma >= BAR),
        "pct_angles_lt_10": float(100.0 * np.mean(np.abs(ang) < 10.0)),
        "pct_irregular_valence": pct_irregular_valence(points, connectivity),
    }
