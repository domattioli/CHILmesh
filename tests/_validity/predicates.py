"""Geometric predicates for mesh-element validation (spec 007).

All predicates take an explicit ``tol`` and use a robust orientation sign.
"""
from __future__ import annotations

import numpy as np


def bbox_diag(points_xy: np.ndarray) -> float:
    """Return ||(max-min)||_2 of an (N, 2) point array. Returns 0 for empty/single point."""
    if points_xy.shape[0] == 0:
        return 0.0
    bbox_min = points_xy.min(axis=0)
    bbox_max = points_xy.max(axis=0)
    return float(np.linalg.norm(bbox_max - bbox_min))


def effective_tol(points_xy: np.ndarray, tol_override: float | None) -> float:
    """Spec FR-013: 1e-12 * bbox_diag (floored at 1e-15)."""
    if tol_override is not None:
        return float(tol_override)
    return max(1e-15, 1e-12 * bbox_diag(points_xy))


def classify_element(row: np.ndarray) -> str:
    """Classify a connectivity row.

    Returns one of: ``TRI``, ``QUAD``, ``DEGENERATE_QUAD_DUPLICATE_VERTEX``,
    ``UNSUPPORTED_ELEMENT_ARITY``.
    """
    cols = len(row)
    if cols < 3:
        return "UNSUPPORTED_ELEMENT_ARITY"
    if cols > 4:
        return "UNSUPPORTED_ELEMENT_ARITY"
    if cols == 3:
        a, b, c = (int(x) for x in row)
        if len({a, b, c}) < 3:
            return "UNSUPPORTED_ELEMENT_ARITY"
        return "TRI"
    a, b, c, d = (int(x) for x in row)
    if d == -1 or d == a or d == b or d == c:
        if len({a, b, c}) < 3:
            return "UNSUPPORTED_ELEMENT_ARITY"
        return "TRI"
    if len({a, b, c, d}) < 4:
        return "DEGENERATE_QUAD_DUPLICATE_VERTEX"
    return "QUAD"


def element_vertex_ids(row: np.ndarray) -> tuple[int, ...]:
    """Distinct vertex IDs in CCW order (drops -1 padding and trailing duplicate of v0)."""
    verts = [int(x) for x in row if int(x) != -1]
    if len(verts) == 4 and verts[3] == verts[0]:
        verts = verts[:3]
    return tuple(verts)


def _orient(ax: float, ay: float, bx: float, by: float, cx: float, cy: float, tol: float) -> int:
    cross = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
    if cross > tol:
        return 1
    if cross < -tol:
        return -1
    return 0


def segment_proper_cross(
    p: np.ndarray, q: np.ndarray, r: np.ndarray, s: np.ndarray, tol: float
) -> bool:
    """True iff open segments (p,q) and (r,s) intersect at strict interior point of both.

    Shared endpoints (any of p,q equal to any of r,s within tol) → False.
    Collinear overlap → False (treated as touching, not crossing).
    """
    if (
        _points_close(p, r, tol) or _points_close(p, s, tol)
        or _points_close(q, r, tol) or _points_close(q, s, tol)
    ):
        return False
    o1 = _orient(p[0], p[1], q[0], q[1], r[0], r[1], tol)
    o2 = _orient(p[0], p[1], q[0], q[1], s[0], s[1], tol)
    o3 = _orient(r[0], r[1], s[0], s[1], p[0], p[1], tol)
    o4 = _orient(r[0], r[1], s[0], s[1], q[0], q[1], tol)
    if o1 == 0 or o2 == 0 or o3 == 0 or o4 == 0:
        return False
    return (o1 != o2) and (o3 != o4)


def _points_close(a: np.ndarray, b: np.ndarray, tol: float) -> bool:
    return float(np.linalg.norm(a - b)) <= tol


def is_self_intersecting_quad(quad_xy: np.ndarray, tol: float) -> tuple[bool, tuple[int, int] | None]:
    """Quad vertices (4, 2) → (is_bowtie, crossing_edge_pair_indices) or (False, None).

    A quad ``v0,v1,v2,v3`` is bowtie iff edge (v0,v1) properly crosses (v2,v3),
    OR edge (v1,v2) properly crosses (v3,v0).
    """
    v0, v1, v2, v3 = quad_xy[0], quad_xy[1], quad_xy[2], quad_xy[3]
    if segment_proper_cross(v0, v1, v2, v3, tol):
        return True, (0, 2)
    if segment_proper_cross(v1, v2, v3, v0, tol):
        return True, (1, 3)
    return False, None


def signed_area_polygon(poly_xy: np.ndarray) -> float:
    """Signed area of an (N, 2) polygon (CCW positive)."""
    n = poly_xy.shape[0]
    if n < 3:
        return 0.0
    x = poly_xy[:, 0]
    y = poly_xy[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def point_strictly_in_polygon(p: np.ndarray, poly_xy: np.ndarray, tol: float) -> bool:
    """Ray-casting PIP. Strict interior — on-edge → False.

    Horizontal edges are skipped (they don't contribute to a horizontal-ray
    crossing-number test). Bbox pre-test cheaply rejects obvious misses.
    """
    n = poly_xy.shape[0]
    if n < 3:
        return False
    xmin = poly_xy[:, 0].min()
    xmax = poly_xy[:, 0].max()
    ymin = poly_xy[:, 1].min()
    ymax = poly_xy[:, 1].max()
    if p[0] < xmin - tol or p[0] > xmax + tol or p[1] < ymin - tol or p[1] > ymax + tol:
        return False
    for i in range(n):
        a = poly_xy[i]
        b = poly_xy[(i + 1) % n]
        if _point_on_segment(p, a, b, tol):
            return False
    inside = False
    for i in range(n):
        a = poly_xy[i]
        b = poly_xy[(i + 1) % n]
        dy = b[1] - a[1]
        if abs(dy) <= tol:
            continue
        if (a[1] > p[1]) != (b[1] > p[1]):
            x_int = a[0] + (p[1] - a[1]) * (b[0] - a[0]) / dy
            if x_int > p[0] + tol:
                inside = not inside
    return inside


def _point_on_segment(p: np.ndarray, a: np.ndarray, b: np.ndarray, tol: float) -> bool:
    if _orient(a[0], a[1], b[0], b[1], p[0], p[1], tol) != 0:
        return False
    dx = b - a
    L2 = float(dx[0] ** 2 + dx[1] ** 2)
    if L2 <= tol * tol:
        return _points_close(p, a, tol)
    t = ((p[0] - a[0]) * dx[0] + (p[1] - a[1]) * dx[1]) / L2
    return -1e-12 <= t <= 1.0 + 1e-12


def element_bboxes(poly_xy_per_elem: list[np.ndarray]) -> np.ndarray:
    """(n_elems, 4) array of [min_x, min_y, max_x, max_y]."""
    out = np.empty((len(poly_xy_per_elem), 4), dtype=float)
    for i, poly in enumerate(poly_xy_per_elem):
        out[i, 0] = poly[:, 0].min()
        out[i, 1] = poly[:, 1].min()
        out[i, 2] = poly[:, 0].max()
        out[i, 3] = poly[:, 1].max()
    return out
