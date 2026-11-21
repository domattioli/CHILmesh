import numpy as np
import pytest

from chilmesh import check_front_invariants, check_quad_against_front


SQUARE = np.array([(x, 0) for x in range(10)] + [(10, y) for y in range(10)]
                  + [(x, 10) for x in range(10, 0, -1)]
                  + [(0, y) for y in range(10, 0, -1)], dtype=float)
L_SHAPE = np.array([(0, 0), (10, 0), (10, 5), (5, 5), (5, 10), (0, 10)], dtype=float)
CELL = np.array([(0, 0), (1, 0), (1, 1), (0, 1)], dtype=float)
NOTCH = np.array([(9, 5), (10, 5), (5, 10), (5, 9)], dtype=float)


def transformed(points):
    angle = np.deg2rad(37)
    matrix = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    return np.asarray(points) @ matrix.T * 1000 + [123, -45]


@pytest.mark.parametrize('transform', [lambda x: x, transformed], ids=['original', 'transformed'])
@pytest.mark.parametrize('quad,expected', [
    (CELL, 'ok'),
    (CELL[::-1], 'inverted'),
    (np.array([(9.5, 2.5), (10.5, 2.5), (10.5, 3.5), (9.5, 3.5)]), 'crosses_front'),
    (np.array([(2, 2), (3, 3), (2, 3), (3, 2)]), 'self_intersecting'),
    (np.array([(2, 2), (3, 2), (4, 2), (5, 2)]), 'degenerate'),
])
def test_square(quad, expected, transform):
    assert check_quad_against_front(transform(SQUARE), transform(quad)) == (expected == 'ok', expected)


@pytest.mark.parametrize('transform', [lambda x: x, transformed], ids=['original', 'transformed'])
def test_notch(transform):
    assert check_quad_against_front(transform(L_SHAPE), transform(NOTCH)) == (False, 'outside_region')


@pytest.mark.parametrize('transform', [lambda x: x, transformed], ids=['original', 'transformed'])
def test_overlap(transform):
    quad = np.array([(2, 2), (4, 2), (4, 4), (2, 4)])
    placed = np.array([(3, 3), (5, 3), (5, 5), (3, 5)])
    assert check_quad_against_front(transform(SQUARE), transform(quad), [transform(placed)]) == (False, 'overlaps_element')


def test_grid_invariants():
    cells = [CELL + [x, y] for x in range(10) for y in range(10)]
    assert check_front_invariants(SQUARE, [], cells)['ok']
    assert not check_front_invariants(SQUARE, [], cells[:-1])['area_conserved']


@pytest.mark.parametrize('transform', [lambda x: x, transformed], ids=['original', 'transformed'])
def test_duplicate_element(transform):
    quad = CELL + [2, 2]
    assert check_quad_against_front(transform(SQUARE), transform(quad), [transform(quad)]) == (False, 'overlaps_element')


@pytest.mark.parametrize('transform', [lambda x: x, transformed], ids=['original', 'transformed'])
def test_hanging_vertex(transform):
    front = np.array([(0, 0), (1, 0), (2, 0), (2, 2), (0, 2)])
    quad = np.array([(0, 0), (2, 0), (2, 1), (0, 1)])
    assert check_quad_against_front(transform(front), transform(quad)) == (False, 'hanging_vertex')


@pytest.mark.parametrize('transform', [lambda x: x, transformed], ids=['original', 'transformed'])
@pytest.mark.parametrize('quad', [
    np.array([(2, 2), (4, 2), (3, 2.5), (2, 4)]),
    np.array([(10, 8), (10, 9), (10, 10), (9, 10)]),
], ids=['reflex', 'collinear'])
def test_nonconvex_quad(transform, quad):
    # Both polygons have positive area, so the degeneracy check does not fire.
    assert check_quad_against_front(transform(SQUARE), transform(quad)) == (False, 'nonconvex')


@pytest.mark.parametrize('quad', [
    np.array([(2, 2), (8, 2), (8, 4), (5, 3.05)]),
    np.array([(2, 2), (8, 2), (8.5, 2.01), (2.5, 2.01)]),
], ids=['large-angle', 'small-angle'])
def test_angle_margin(quad):
    assert check_quad_against_front(SQUARE, quad) == (True, 'ok')
    assert check_quad_against_front(SQUARE, quad, angle_margin_deg=2) == (False, 'nonconvex')


# Frozen reference from before the performance changes.
def _reference_polygon(points):
    points = np.asarray(points, dtype=float).reshape(-1, 2)
    if len(points) > 1 and np.array_equal(points[0], points[-1]):
        points = points[:-1]
    return points


def _reference_cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def _reference_area(poly):
    if len(poly) < 3:
        return 0.0
    p = poly - poly[0]
    return sum(_reference_cross(p[i], p[(i + 1) % len(p)]) for i in range(len(p))) / 2


def _reference_edges(poly):
    return [(poly[i], poly[(i + 1) % len(poly)]) for i in range(len(poly))]


def _reference_distance_segment(p, a, b):
    ab = b - a
    denom = np.dot(ab, ab)
    t = np.clip(np.dot(p - a, ab) / denom, 0, 1) if denom else 0
    return np.linalg.norm(p - (a + t * ab))


def _reference_proper_cross(a, b, c, d, tol):
    ab, cd = b - a, d - c
    x, y = _reference_cross(ab, c - a), _reference_cross(ab, d - a)
    z, w = _reference_cross(cd, a - c), _reference_cross(cd, b - c)
    eps_ab, eps_cd = tol * np.linalg.norm(ab), tol * np.linalg.norm(cd)
    return ((x > eps_ab and y < -eps_ab or x < -eps_ab and y > eps_ab)
            and (z > eps_cd and w < -eps_cd or z < -eps_cd and w > eps_cd))


def _reference_location(point, poly, tol):
    """Return zero on the boundary, one inside, and minus one outside."""
    edges = _reference_edges(poly)
    if any(_reference_distance_segment(point, a, b) <= tol for a, b in edges):
        return 0
    inside = False
    for a, b in edges:
        if (a[1] > point[1]) != (b[1] > point[1]):
            x = a[0] + (point[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
            if x > point[0]:
                inside = not inside
    return 1 if inside else -1


def _reference_outside(quad, polygon, tol):
    samples = list(quad) + [np.mean(quad, axis=0)]
    samples += [(a + b) / 2 for a, b in _reference_edges(quad)]
    return any(_reference_location(p, polygon, tol) < 0 for p in samples)


def _referencecheck_quad_against_front(front_xy, quad_xy, placed_quads=(), domain_xy=None, rel_tol=1e-4):
    """Check CCW polygons in a single loop; rel_tol scales the shortest domain edge.

    Reasons in check order: degenerate, self_intersecting, inverted, nonconvex, crosses_front,
    outside_region, contains_front_vertex, hanging_vertex, overlaps_element,
    outside_domain, ok. Loops have no repeated final point.
    """
    front, quad = _reference_polygon(front_xy), _reference_polygon(quad_xy)
    domain = _reference_polygon(domain_xy) if domain_xy is not None else front
    tol = rel_tol * min(np.linalg.norm(b - a) for a, b in _reference_edges(domain))
    edges = _reference_edges(quad)
    # Unsigned fan area recognizes a bowtie as crossing rather than degenerate.
    extent_area = (abs(_reference_cross(quad[1] - quad[0], quad[2] - quad[0]))
                   + abs(_reference_cross(quad[2] - quad[0], quad[3] - quad[0]))) / 2
    if any(np.linalg.norm(b - a) < tol for a, b in edges) or extent_area < tol ** 2:
        return False, 'degenerate'
    if any(_reference_proper_cross(*edges[i], *edges[j], tol) for i, j in ((0, 2), (1, 3))):
        return False, 'self_intersecting'
    if _reference_area(quad) <= 0:
        return False, 'inverted'
    for i in range(len(quad)):
        previous = quad[i - 1] - quad[i]
        following = quad[(i + 1) % len(quad)] - quad[i]
        if _reference_cross(following, previous) <= tol * np.linalg.norm(previous) * np.linalg.norm(following):
            return False, 'nonconvex'
    if any(_reference_proper_cross(*qe, *fe, tol) for qe in edges for fe in _reference_edges(front)):
        return False, 'crosses_front'
    if _reference_outside(quad, front, tol):
        return False, 'outside_region'
    if any(_reference_location(p, quad, tol) == 1 for p in front):
        return False, 'contains_front_vertex'
    if any(all(np.linalg.norm(p - v) > tol for v in quad) and
           any(_reference_distance_segment(p, a, b) <= tol for a, b in edges) for p in front):
        return False, 'hanging_vertex'
    for element in placed_quads:
        element = _reference_polygon(element)
        if (any(_reference_proper_cross(*qe, *ee, tol) for qe in edges for ee in _reference_edges(element))
                or any(_reference_location(p, element, tol) == 1 for p in quad)
                or any(_reference_location(p, quad, tol) == 1 for p in element)
                or any(_reference_location(p, element, tol) == 1 for p in
                       [np.mean(quad, axis=0)] + [(a + b) / 2 for a, b in edges])
                or _reference_location(np.mean(element, axis=0), quad, tol) == 1
                or (len(element) == len(quad) and all(
                    any(np.linalg.norm(p - q) <= tol for q in element) for p in quad))):
            return False, 'overlaps_element'
    if domain_xy is not None and _reference_outside(quad, domain, tol):
        return False, 'outside_domain'
    return True, 'ok'



def test_randomized_reference_agreement():
    rng = np.random.default_rng(2381)
    placed = [CELL + offset for offset in [(1, 1), (3, 3), (7, 2), (2, 7)]]
    for i in range(2000):
        if i % 2:
            origin = rng.uniform(-1, 10, size=2)
            size = rng.uniform(0.01, 3, size=2)
            quad = origin + CELL * size
        else:
            quad = rng.uniform(-1, 11, size=(4, 2))
        assert check_quad_against_front(L_SHAPE, quad, placed, L_SHAPE) == (
            _referencecheck_quad_against_front(L_SHAPE, quad, placed, L_SHAPE))
