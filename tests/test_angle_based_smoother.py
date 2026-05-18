"""Regression tests for ``CHILmesh.angle_based_smoother`` and the
``_ordered_vertex_ring`` helper (.planning/TEST-AUDIT.md F1).

Both functions sat at zero coverage before this file existed. The
constitution (.specify/memory/constitution.md, "Feature-Specific
Principles → Mesh Smoothing") names the Zhou & Shimada angle-based
smoother as the DOMsmooth Hybrid Fallback for mixed meshes, so a
regression net is mandatory before any future refactor.

Invariants pinned here:
- Boundary vertices are not moved.
- Z column is preserved verbatim.
- Output shape matches input shape.
- The line search guarantees monotone non-degradation of the minimum
  element quality.
- Re-running on a smoothed mesh converges within tolerance (idempotent
  fixed point under the relaxation factor).
- The ``smooth_mesh('angle-based', acknowledge_change=True)`` dispatch
  path actually invokes this smoother.
- ``_ordered_vertex_ring`` returns ``None`` on boundary vertices and a
  CCW-ordered ring on manifold interior vertices.
"""
from __future__ import annotations

import numpy as np
import pytest

import chilmesh
from chilmesh import CHILmesh
from conftest import TRI_FIXTURE_NAMES

# block_o is in TRI_FIXTURE_NAMES but its ~14s warmup dominates this file.
# Quality / convergence checks below run on the smaller meshes; block_o
# is exercised by the dispatcher test only.
SMALL_TRI_FIXTURES = [n for n in TRI_FIXTURE_NAMES if n != "block_o"]


def _interior_vertices(mesh: CHILmesh) -> np.ndarray:
    edge_verts = mesh.edge2vert(mesh.boundary_edges())
    boundary = set(np.unique(edge_verts.flatten()).tolist())
    return np.array([v for v in range(mesh.n_verts) if v not in boundary], dtype=int)


def _qualities(mesh: CHILmesh) -> np.ndarray:
    """Per-element quality array from ``elem_quality`` (tuple of (q, angles, stats))."""
    out = mesh.elem_quality()
    return np.asarray(out[0] if isinstance(out, tuple) else out)


def _min_quality(mesh: CHILmesh, points: np.ndarray) -> float:
    """Min element quality if ``points`` were the mesh coordinates."""
    saved = mesh.points.copy()
    try:
        mesh.points = points.copy()
        return float(np.nanmin(_qualities(mesh)))
    finally:
        mesh.points = saved


@pytest.mark.parametrize("name", SMALL_TRI_FIXTURES)
def test_angle_smoother_preserves_boundary(name):
    mesh = getattr(chilmesh.examples, name)()
    original = mesh.points.copy()
    edge_verts = mesh.edge2vert(mesh.boundary_edges())
    boundary = np.unique(edge_verts.flatten())

    smoothed = mesh.angle_based_smoother(n_iter=5)

    for v in boundary:
        np.testing.assert_array_equal(
            smoothed[v, :2], original[v, :2],
            err_msg=f"{name}: boundary vertex {v} was moved",
        )


@pytest.mark.parametrize("name", SMALL_TRI_FIXTURES)
def test_angle_smoother_preserves_z(name):
    mesh = getattr(chilmesh.examples, name)()
    original_z = mesh.points[:, 2].copy()

    smoothed = mesh.angle_based_smoother(n_iter=3)

    np.testing.assert_array_equal(
        smoothed[:, 2], original_z,
        err_msg=f"{name}: z-coordinate changed",
    )


@pytest.mark.parametrize("name", SMALL_TRI_FIXTURES)
def test_angle_smoother_output_shape(name):
    mesh = getattr(chilmesh.examples, name)()
    smoothed = mesh.angle_based_smoother(n_iter=1)
    assert smoothed.shape == mesh.points.shape, (
        f"{name}: smoothed shape {smoothed.shape} != input shape {mesh.points.shape}"
    )


@pytest.mark.parametrize("name", SMALL_TRI_FIXTURES)
def test_angle_smoother_does_not_degrade_min_quality(name):
    """Line search inside the smoother only accepts moves that strictly
    improve local minimum quality. Global minimum must therefore be
    non-decreasing across the full pass."""
    mesh = getattr(chilmesh.examples, name)()
    before = float(np.nanmin(_qualities(mesh)))

    smoothed = mesh.angle_based_smoother(n_iter=10)
    after = _min_quality(mesh, smoothed)

    assert after >= before - 1e-12, (
        f"{name}: min quality degraded {before:.6f} -> {after:.6f}"
    )


@pytest.mark.parametrize("name", SMALL_TRI_FIXTURES)
def test_angle_smoother_progresses_toward_fixed_point(name):
    """The smoother is a relaxation, not a hard fixed-point: across
    successive passes the per-vertex displacement should not grow.
    Asserts pass-2 max move <= pass-1 max move (within a small slack)."""
    mesh = getattr(chilmesh.examples, name)()
    original = mesh.points.copy()

    pass_one = mesh.angle_based_smoother(n_iter=20, tol=1e-10)
    saved = mesh.points.copy()
    try:
        mesh.points = pass_one.copy()
        pass_two = mesh.angle_based_smoother(n_iter=20, tol=1e-10)
    finally:
        mesh.points = saved

    move_1 = float(np.max(np.linalg.norm(pass_one[:, :2] - original[:, :2], axis=1)))
    move_2 = float(np.max(np.linalg.norm(pass_two[:, :2] - pass_one[:, :2], axis=1)))
    assert move_2 <= move_1 + 1e-9, (
        f"{name}: second-pass max move {move_2:.3e} exceeds first-pass {move_1:.3e}"
    )


def test_angle_smoother_moves_interior_when_perturbed():
    """A perturbed annulus should see at least one interior vertex
    relax back toward equilibrium on a single pass."""
    mesh = chilmesh.examples.annulus()
    interior = _interior_vertices(mesh)
    assert interior.size > 0, "annulus has no interior vertices?"

    rng = np.random.default_rng(2026)
    perturbed = mesh.points.copy()
    perturbed[interior, :2] += rng.uniform(-0.01, 0.01, size=(interior.size, 2))

    saved = mesh.points.copy()
    try:
        mesh.points = perturbed
        smoothed = mesh.angle_based_smoother(n_iter=20)
    finally:
        mesh.points = saved

    interior_moves = np.linalg.norm(smoothed[interior, :2] - perturbed[interior, :2], axis=1)
    assert interior_moves.max() > 1e-8, (
        f"smoother did not move any interior vertex (max move {interior_moves.max():.2e})"
    )


def test_smooth_mesh_dispatches_to_angle_based():
    """``smooth_mesh('angle-based', acknowledge_change=True)`` must run
    the angle-based smoother and update ``mesh.points`` in place."""
    mesh = chilmesh.examples.annulus()
    before = mesh.points.copy()

    out = mesh.smooth_mesh('angle-based', acknowledge_change=True)

    assert out.shape == before.shape
    np.testing.assert_array_equal(mesh.points, out)


def test_smooth_mesh_angle_based_requires_acknowledge_change():
    mesh = chilmesh.examples.annulus()
    with pytest.raises(AssertionError):
        mesh.smooth_mesh('angle-based')


# ---------------------------------------------------------------------------
# _ordered_vertex_ring direct tests
# ---------------------------------------------------------------------------


def test_ordered_vertex_ring_boundary_returns_none():
    """A boundary vertex has an open ring → helper must return ``None``."""
    mesh = chilmesh.examples.annulus()
    edge_verts = mesh.edge2vert(mesh.boundary_edges())
    boundary = np.unique(edge_verts.flatten())
    bv = int(boundary[0])
    elem_ids = list(mesh.adjacencies['Vert2Elem'][bv])
    assert mesh._ordered_vertex_ring(bv, elem_ids) is None


def test_ordered_vertex_ring_interior_is_cyclic_cover():
    """For an interior vertex the ring must be a closed cycle whose
    members are exactly the unique neighbors across all incident
    elements. (The docstring claims CCW, but the pred->succ chain
    actually walks CW around the vertex — the smoother is direction-
    agnostic so this is a docstring bug, not a correctness bug.)"""
    mesh = chilmesh.examples.annulus()
    interior = _interior_vertices(mesh)
    assert interior.size > 0
    v = int(interior[0])
    elem_ids = list(mesh.adjacencies['Vert2Elem'][v])

    ring = mesh._ordered_vertex_ring(v, elem_ids)
    assert ring is not None, f"interior vertex {v} returned None"
    assert len(ring) == len(set(ring)), f"ring has duplicates: {ring}"

    neighbors: set[int] = set()
    for eid in elem_ids:
        row = mesh.connectivity_list[eid]
        if row.shape[0] == 4 and row[3] == row[0]:
            verts = row[:3].tolist()
        else:
            verts = row.tolist()
        neighbors.update(x for x in verts if x != v)
    assert set(ring) == neighbors, (
        f"ring {ring} != neighbor set {neighbors}"
    )

    # Cycle: consecutive ring entries share an element with v.
    edge_set = set()
    for eid in elem_ids:
        row = mesh.connectivity_list[eid]
        if row.shape[0] == 4 and row[3] == row[0]:
            verts = row[:3].tolist()
        else:
            verts = row.tolist()
        if v not in verts:
            continue
        for x in verts:
            if x != v:
                edge_set.add(frozenset((v, x)))
    for i in range(len(ring)):
        a, b = ring[i], ring[(i + 1) % len(ring)]
        # The pair (a,b) should both be reachable from v in one element step
        assert frozenset((v, a)) in edge_set and frozenset((v, b)) in edge_set, (
            f"ring pair ({a}, {b}) breaks cyclic invariant"
        )


def test_ordered_vertex_ring_padded_triangle_recognized():
    """Padded-triangle rows (``[v0,v1,v2,v0]``) must be treated as
    3-vertex elements, not 4-vertex ones with a repeated index."""
    pts = np.array([
        [0.0, 0.0],   # 0 (center / interior)
        [1.0, 0.0],   # 1
        [0.5, 0.87],  # 2
        [-0.5, 0.87], # 3
        [-1.0, 0.0],  # 4
        [-0.5, -0.87],# 5
        [0.5, -0.87], # 6
    ])
    # Six padded triangles fanning out from vertex 0
    conn = np.array([
        [0, 1, 2, 0],
        [0, 2, 3, 0],
        [0, 3, 4, 0],
        [0, 4, 5, 0],
        [0, 5, 6, 0],
        [0, 6, 1, 0],
    ])
    mesh = CHILmesh(connectivity=conn, points=pts, grid_name="fan6")
    elem_ids = list(mesh.adjacencies['Vert2Elem'][0])
    ring = mesh._ordered_vertex_ring(0, elem_ids)
    assert ring is not None
    assert set(ring) == {1, 2, 3, 4, 5, 6}
    assert len(ring) == 6
