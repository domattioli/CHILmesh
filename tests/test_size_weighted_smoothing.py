"""Optional ``size_fn`` for ``direct_smoother`` / ``angle_based_smoother`` (#197, option A).

Pinned here:
- ``size_fn=None`` is the unchanged algorithm (the byte-identity proof against
  stored hashes lives in the supervisor harness; here we check None == omitted).
- With a size field, every fixture (tri, quad) keeps a finite, boundary-fixed,
  non-inverted mesh for both FEM solvers and the angle smoother.
- On a radially graded mesh with jittered interior nodes, the size-weighted run
  ends closer to the graded target sizes than the size-blind base, which
  erodes the grading (direction only, not exact numbers).
- ``smooth_mesh`` forwards ``size_fn`` to the ``fem`` and ``angle-based`` paths.
"""
from __future__ import annotations

import numpy as np
import pytest

from chilmesh import CHILmesh
from chilmesh.CHILmesh import _shoelace, _unique_in_order
from conftest import FIXTURE_NAMES, _load


def _graded_x_size(mesh: CHILmesh):
    """Analytic size field growing linearly with x (factor 1 to 2.5 over the domain)."""
    xy = mesh.points[:, :2]
    lo, hi = xy[:, 0].min(), xy[:, 0].max()
    e = mesh.adjacencies["Edge2Vert"]
    h0 = float(np.linalg.norm(xy[e[:, 0]] - xy[e[:, 1]], axis=1).mean())
    return lambda q: h0 * (1.0 + 1.5 * (q[:, 0] - lo) / (hi - lo))


def _boundary_nodes(mesh: CHILmesh) -> np.ndarray:
    return np.unique(mesh.edge2vert(mesh.boundary_edges()).flatten())


def _min_signed_area(mesh: CHILmesh, pts: np.ndarray) -> float:
    return min(_shoelace(pts[_unique_in_order(r)][:, :2]) for r in mesh.connectivity_list)


def _graded_mesh(quad: bool, jitter: float = 0.0, n: int = 21, power: float = 2.0) -> CHILmesh:
    """Square grid mapped by ``z -> z |z|**(power-1)``: fine at the centre, coarse at the rim.

    The map is near-conformal, so the grading is roughly isotropic (a scalar size
    field describes it). Interior nodes can be jittered to make the base smoother
    work on something that is not already at its fixed point.
    """
    t = np.linspace(-1.0, 1.0, n)
    X, Y = np.meshgrid(t, t)
    z = np.column_stack([X.ravel(), Y.ravel()])
    r = np.linalg.norm(z, axis=1, keepdims=True)
    pts = z * np.where(r > 0, r, 1.0) ** (power - 1.0)
    idx = np.arange(n * n).reshape(n, n)
    inner = idx[1:-1, 1:-1].ravel()
    if jitter:
        h = 2.0 / (n - 1) * np.maximum(np.linalg.norm(pts[inner], axis=1), 0.1)
        pts[inner] += np.random.default_rng(0).uniform(-jitter, jitter, (len(inner), 2)) * h[:, None]
    conn = []
    for j in range(n - 1):
        for i in range(n - 1):
            a, b, c, d = idx[j, i], idx[j, i + 1], idx[j + 1, i + 1], idx[j + 1, i]
            conn += [[a, b, c, d]] if quad else [[a, b, c], [a, c, d]]
    return CHILmesh(connectivity=np.array(conn),
                    points=np.column_stack([pts, np.zeros(len(pts))]))


def _own_size(mesh: CHILmesh):
    """Per-vertex mean incident edge length of ``mesh``, linearly interpolated."""
    from scipy.interpolate import LinearNDInterpolator
    xy = mesh.points[:, :2]
    e = mesh.adjacencies["Edge2Vert"]
    length = np.linalg.norm(xy[e[:, 0]] - xy[e[:, 1]], axis=1)
    total, count = np.zeros(len(xy)), np.zeros(len(xy))
    for col in (0, 1):
        np.add.at(total, e[:, col], length)
        np.add.at(count, e[:, col], 1.0)
    return LinearNDInterpolator(xy, total / count, fill_value=float(length.mean()))


def _size_error(mesh: CHILmesh, pts: np.ndarray, fn) -> float:
    """Mean relative edge-length error ``|L - h(mid)| / h(mid)`` against the target field."""
    e = mesh.adjacencies["Edge2Vert"]
    p1, p2 = pts[e[:, 0], :2], pts[e[:, 1], :2]
    h = fn(0.5 * (p1 + p2))
    return float(np.mean(np.abs(np.linalg.norm(p1 - p2, axis=1) - h) / h))


def _mixed_mesh() -> CHILmesh:
    """3x3-cell grid: left column quads, the rest triangles padded to 4 columns."""
    g = np.arange(16).reshape(4, 4)
    xy = np.array([[i, j] for j in range(4) for i in range(4)], dtype=float)
    xy[g[1, 1]] += [0.15, -0.1]
    xy[g[2, 2]] += [-0.1, 0.12]
    conn = []
    for j in range(3):
        for i in range(3):
            a, b, c, d = g[j, i], g[j, i + 1], g[j + 1, i + 1], g[j + 1, i]
            if i == 0:
                conn.append([a, b, c, d])
            else:
                conn += [[a, b, c, a], [a, c, d, a]]
    return CHILmesh(connectivity=np.array(conn),
                    points=np.column_stack([xy, np.zeros(16)]))


SMOOTHERS = [
    ("fem_direct", lambda m, **kw: m.direct_smoother(**kw)),
    ("fem_iterative", lambda m, **kw: m.direct_smoother(solver="iterative", **kw)),
    ("angle", lambda m, **kw: m.angle_based_smoother(n_iter=10, **kw)),
]


# angle smoother on block_o (~5k elements, pure Python loop) is slow: slow-marked.
FIXTURE_CASES = [
    pytest.param(name, label, run,
                 marks=[pytest.mark.slow] if (name == "block_o" and label == "angle") else [],
                 id=f"{name}-{label}")
    for name in FIXTURE_NAMES for label, run in SMOOTHERS
]


@pytest.mark.parametrize("name,label,run", FIXTURE_CASES)
def test_size_fn_none_equals_omitted(name, label, run):
    mesh = _load(name).copy()
    np.testing.assert_array_equal(run(mesh), run(mesh, size_fn=None))


@pytest.mark.parametrize("name,label,run", FIXTURE_CASES)
def test_size_weighted_valid_on_fixtures(name, label, run):
    mesh = _load(name).copy()
    new = run(mesh, size_fn=_graded_x_size(mesh))
    assert new.shape == mesh.points.shape
    assert np.isfinite(new).all()
    bnd = _boundary_nodes(mesh)
    # penalty pinning (kinf) holds boundary to ~1e-11, as in the base path
    np.testing.assert_allclose(new[bnd], mesh.points[bnd], atol=1e-8, rtol=0)
    np.testing.assert_array_equal(new[:, 2], mesh.points[:, 2])
    assert _min_signed_area(mesh, new) > 0.0, f"{name}/{label}: inverted element"


@pytest.mark.parametrize("label,run", SMOOTHERS, ids=[s[0] for s in SMOOTHERS])
def test_size_weighted_mixed_mesh(label, run):
    mesh = _mixed_mesh()
    new = run(mesh, size_fn=lambda q: 1.0 + 0.2 * q[:, 0])
    assert np.isfinite(new).all()
    bnd = _boundary_nodes(mesh)
    # penalty pinning (kinf) holds boundary to ~1e-11, as in the base path
    np.testing.assert_allclose(new[bnd], mesh.points[bnd], atol=1e-8, rtol=0)
    assert _min_signed_area(mesh, new) > 0.0


# The angle smoother needs its default pass count for the size pull to accumulate.
GRADED_SMOOTHERS = SMOOTHERS[:2] + [("angle", lambda m, **kw: m.angle_based_smoother(**kw))]


@pytest.mark.parametrize("quad", [False, True], ids=["tri", "quad"])
@pytest.mark.parametrize("label,run", GRADED_SMOOTHERS, ids=[s[0] for s in GRADED_SMOOTHERS])
def test_graded_size_field_tracks_target_better_than_base(quad, label, run):
    mesh = _graded_mesh(quad, jitter=0.2)
    # Target = the sizes of the same grid before the interior was jittered.
    fn = _own_size(_graded_mesh(quad, jitter=0.0))
    base = _size_error(mesh, run(mesh), fn)
    sized = _size_error(mesh, run(mesh, size_fn=fn), fn)
    assert sized < base, f"size-weighted error {sized:.4f} not below base {base:.4f}"


@pytest.mark.parametrize("method", ["fem", "angle-based"])
def test_smooth_mesh_forwards_size_fn(method, monkeypatch):
    mesh = _load("annulus").copy()
    seen = {}
    target = "direct_smoother" if method == "fem" else "angle_based_smoother"

    def spy(self, **kw):
        seen.update(kw)
        return self.points.copy()

    monkeypatch.setattr(CHILmesh, target, spy)
    fn = _graded_x_size(mesh)
    mesh.smooth_mesh(method, acknowledge_change=True, size_fn=fn)
    assert seen.get("size_fn") is fn


def test_smooth_mesh_without_size_fn_passes_no_size_kwarg(monkeypatch):
    mesh = _load("annulus").copy()
    seen = {}

    def spy(self, **kw):
        seen.update(kw)
        return self.points.copy()

    monkeypatch.setattr(CHILmesh, "direct_smoother", spy)
    mesh.smooth_mesh("fem", acknowledge_change=True)
    assert "size_fn" not in seen


@pytest.mark.parametrize("bad", [lambda q: -np.ones(len(q)), lambda q: np.ones(3)])
def test_invalid_size_fn_raises(bad):
    mesh = _load("annulus").copy()
    with pytest.raises(ValueError):
        mesh.direct_smoother(size_fn=bad)
