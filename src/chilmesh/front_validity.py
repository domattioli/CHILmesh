"""Brute-force checks for single counter-clockwise loops without repeated endpoints.

Relative tolerances scale the shortest domain edge. Reason codes are returned
in the order documented by :func:`check_quad_against_front`.
"""

import math

import numpy as np


def _polygon(points):
    points = np.asarray(points, dtype=float).reshape(-1, 2)
    if len(points) > 1 and np.array_equal(points[0], points[-1]):
        points = points[:-1]
    return points


def _cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def _area(poly):
    if len(poly) < 3:
        return 0.0
    p = poly - poly[0]
    return sum(_cross(p[i], p[(i + 1) % len(p)]) for i in range(len(p))) / 2


def _edges(poly):
    return [(poly[i], poly[(i + 1) % len(poly)]) for i in range(len(poly))]


def _distance_segment(p, a, b):
    ab = b - a
    denom = np.dot(ab, ab)
    t = min(1.0, max(0.0, np.dot(p - a, ab) / denom)) if denom else 0
    delta = p - (a + t * ab)
    return math.hypot(delta[0], delta[1])


def _proper_cross(a, b, c, d, tol):
    ab, cd = b - a, d - c
    x, y = _cross(ab, c - a), _cross(ab, d - a)
    z, w = _cross(cd, a - c), _cross(cd, b - c)
    eps_ab, eps_cd = tol * math.hypot(*ab), tol * math.hypot(*cd)
    return ((x > eps_ab and y < -eps_ab or x < -eps_ab and y > eps_ab)
            and (z > eps_cd and w < -eps_cd or z < -eps_cd and w > eps_cd))


def _location(point, poly, tol):
    """Return zero on the boundary, one inside, and minus one outside."""
    starts = poly
    ends = np.roll(poly, -1, axis=0)
    edges = ends - starts
    offsets = point - starts
    denom = np.sum(edges * edges, axis=1)
    t = np.divide(np.sum(offsets * edges, axis=1), denom,
                  out=np.zeros_like(denom), where=denom != 0)
    t = np.maximum(0, np.minimum(1, t))
    delta = point - (starts + t[:, None] * edges)
    if np.any(np.hypot(delta[:, 0], delta[:, 1]) <= tol):
        return 0
    crossing = (starts[:, 1] > point[1]) != (ends[:, 1] > point[1])
    a, b = starts[crossing], ends[crossing]
    x = a[:, 0] + (point[1] - a[:, 1]) * (b[:, 0] - a[:, 0]) / (b[:, 1] - a[:, 1])
    return 1 if np.count_nonzero(x > point[0]) % 2 else -1


def _outside(quad, polygon, tol):
    samples = list(quad) + [np.mean(quad, axis=0)]
    samples += [(a + b) / 2 for a, b in _edges(quad)]
    return any(_location(p, polygon, tol) < 0 for p in samples)


def check_quad_against_front(front_xy, quad_xy, placed_quads=(), domain_xy=None, rel_tol=1e-4,
                             *, angle_margin_deg=0.0):
    """Check CCW polygons in a single loop; rel_tol scales the shortest domain edge.

    Reasons in check order: degenerate, self_intersecting, inverted, nonconvex, crosses_front,
    outside_region, contains_front_vertex, hanging_vertex, overlaps_element,
    outside_domain, ok. Loops have no repeated final point. A positive
    angle_margin_deg rejects interior angles outside the specified margin.
    """
    front, quad = _polygon(front_xy), _polygon(quad_xy)
    domain = _polygon(domain_xy) if domain_xy is not None else front
    tol = rel_tol * min(np.linalg.norm(b - a) for a, b in _edges(domain))
    edges = _edges(quad)
    # Unsigned fan area recognizes a bowtie as crossing rather than degenerate.
    extent_area = (abs(_cross(quad[1] - quad[0], quad[2] - quad[0]))
                   + abs(_cross(quad[2] - quad[0], quad[3] - quad[0]))) / 2
    if any(np.linalg.norm(b - a) < tol for a, b in edges) or extent_area < tol ** 2:
        return False, 'degenerate'
    if any(_proper_cross(*edges[i], *edges[j], tol) for i, j in ((0, 2), (1, 3))):
        return False, 'self_intersecting'
    if _area(quad) <= 0:
        return False, 'inverted'
    for i in range(len(quad)):
        previous = quad[i - 1] - quad[i]
        following = quad[(i + 1) % len(quad)] - quad[i]
        if _cross(following, previous) <= tol * np.linalg.norm(previous) * np.linalg.norm(following):
            return False, 'nonconvex'
        if angle_margin_deg:
            angle = math.degrees(math.atan2(_cross(following, previous), np.dot(following, previous)))
            if angle < angle_margin_deg or angle > 180 - angle_margin_deg:
                return False, 'nonconvex'
    quad_min, quad_max = quad.min(axis=0), quad.max(axis=0)
    front_ends = np.roll(front, -1, axis=0)
    relevant = np.all(np.maximum(front, front_ends) >= quad_min, axis=1) & np.all(
        np.minimum(front, front_ends) <= quad_max, axis=1)
    front_edges = zip(front[relevant], front_ends[relevant])
    if any(_proper_cross(*qe, *fe, tol) for fe in front_edges for qe in edges):
        return False, 'crosses_front'
    if _outside(quad, front, tol):
        return False, 'outside_region'
    if any(_location(p, quad, tol) == 1 for p in front):
        return False, 'contains_front_vertex'
    if any(all(np.linalg.norm(p - v) > tol for v in quad) and
           any(_distance_segment(p, a, b) <= tol for a, b in edges) for p in front):
        return False, 'hanging_vertex'
    for element in placed_quads:
        element = _polygon(element)
        if np.any(element.max(axis=0) + tol < quad_min) or np.any(element.min(axis=0) - tol > quad_max):
            continue
        if (any(_proper_cross(*qe, *ee, tol) for qe in edges for ee in _edges(element))
                or any(_location(p, element, tol) == 1 for p in quad)
                or any(_location(p, quad, tol) == 1 for p in element)
                or any(_location(p, element, tol) == 1 for p in
                       [np.mean(quad, axis=0)] + [(a + b) / 2 for a, b in edges])
                or _location(np.mean(element, axis=0), quad, tol) == 1
                or (len(element) == len(quad) and all(
                    any(np.linalg.norm(p - q) <= tol for q in element) for p in quad))):
            return False, 'overlaps_element'
    if domain_xy is not None and _outside(quad, domain, tol):
        return False, 'outside_domain'
    return True, 'ok'


def check_front_invariants(domain_xy, front_xy, elements):
    """Check one CCW loop, without a repeated endpoint, and CCW elements.

    Simplicity uses 1e-9 times the shortest domain edge as its tolerance.
    """
    domain, front = _polygon(domain_xy), _polygon(front_xy)
    elements = [_polygon(element) for element in elements]
    edges = _edges(front)
    tol = 1e-9 * min(np.linalg.norm(b - a) for a, b in _edges(domain))
    simple = len(front) == 0 or (len(front) >= 3 and all(
        not _segments_touch(*edges[i], *edges[j], tol)
        for i in range(len(edges)) for j in range(i + 1, len(edges))
        if j != i + 1 and (i, j) != (0, len(edges) - 1)))
    domain_area = _area(domain)
    result = {
        'simple_loop': bool(simple),
        'area_conserved': bool(abs(_area(front) + sum(_area(e) for e in elements) - domain_area)
                               <= 1e-6 * abs(domain_area)),
        'even_edge_count': len(front) % 2 == 0,
        'no_inverted_element': all(_area(e) > 0 for e in elements),
    }
    result['ok'] = all(result.values())
    return result


def _segments_touch(a, b, c, d, tol=0):
    """Detect any intersection, including endpoint and collinear contact."""
    if _proper_cross(a, b, c, d, tol):
        return True
    return any(_distance_segment(p, x, y) <= tol for p, x, y in
               ((a, c, d), (b, c, d), (c, a, b), (d, a, b)))
