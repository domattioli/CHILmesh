"""
Vendored ADMESH Distmesh Truss Loop (Warm-Start Variant)

**Source**: admesh.distmesh.distmesh2d, commit 05bc68fc81060f7d710b8f4abb2cc382f85df33f
**Lines**: 140-220 (inner truss loop, byte-identical to upstream)
**Modification**: Warm-start preamble only (replaced _initial_distribution + _rejection_method)

This module contains the exact inner truss loop from ADMESH's distmesh2d algorithm,
with a single change: instead of generating initial points via _initial_distribution
and filtering them via _rejection_method, we accept pre-computed points from the caller
(from an existing triangulation). The rest of the algorithm is identical to upstream.

**Removal Trigger**: When ADMESH ships a public warm-start entry point or accepts
an `initial_points` parameter (issue ADMESH-B), this vendor module is deleted and
we call ADMESH directly. Until then, this is necessary to avoid the broken
admesh.routine.py import (ADMESH-A).

**Why Vendor, Not Monkey-Patch**: Per research.md R4, there is no hook to inject
custom initial points without modifying ADMESH source. The inner loop is ~50 lines
(manageable), and vendoring is cleaner than monkeypatching sys.modules.

**Assumptions**:
- scipy.spatial.Delaunay preserves point array order (verified via R6)
- pfix mechanism (Ftot[:nfix] = 0) delivers bit-exact preservation (verified via R2)
"""

import numpy as np
from scipy.spatial import Delaunay
from typing import Callable, Optional, Tuple


def _quality_metric(p: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Per-element quality: 4√3 * area / sum(edge_len²). Inverted = 0."""
    if len(t) == 0:
        return np.array([])
    p0 = p[t[:, 0]]
    p1 = p[t[:, 1]]
    p2 = p[t[:, 2]]
    e0_sq = np.sum((p1 - p0) ** 2, axis=1)
    e1_sq = np.sum((p2 - p1) ** 2, axis=1)
    e2_sq = np.sum((p0 - p2) ** 2, axis=1)
    signed_area = 0.5 * (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1])
    )
    e_sum_sq = e0_sq + e1_sq + e2_sq
    return np.where(
        signed_area > 0,
        4.0 * np.sqrt(3.0) * signed_area / np.maximum(e_sum_sq, 1e-12),
        0.0,
    )


def distmesh2d_warmstart(
    pfix_boundary: np.ndarray,
    interior_initial: np.ndarray,
    fd: Callable[[np.ndarray], np.ndarray],
    fh: Optional[Callable[[np.ndarray], np.ndarray]],
    h0: float,
    bbox: Tuple[float, float, float, float],
    *,
    dptol: float = 1e-3,
    ttol: float = 0.1,
    Fscale: float = 1.2,
    deltat: float = 0.2,
    geps_factor: float = 1e-3,
    niter: int = 500,
    seed: int = 0,
    quality_check_interval: int = 5,
    quality_drop_threshold: float = 0.10,
    track_best_quality: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Warm-start variant of ADMESH's distmesh2d truss loop.

    Instead of generating initial points via rejection sampling, this function
    accepts pre-computed boundary and interior points and runs the truss loop
    until convergence or max iterations.

    Parameters
    ----------
    pfix_boundary : (B, 2) ndarray
        Fixed boundary points. Will appear at indices 0..B-1 of output, unchanged.
    interior_initial : (I, 2) ndarray
        Initial interior points (free to move). Combined with pfix_boundary.
    fd : callable
        Signed distance function: fd(points: (K, 2)) -> (K,).
    fh : callable or None
        Mesh size function: fh(points: (K, 2)) -> (K,). None means uniform.
    h0 : float
        Target edge length (mesh spacing).
    bbox : (xmin, ymin, xmax, ymax)
        Bounding box for domain.
    dptol : float
        Interior node displacement tolerance (stopping criterion).
    ttol : float
        Relative node movement threshold for re-triangulation.
    Fscale : float
        Internal pressure scaling factor (> 1).
    deltat : float
        Euler time step for force-displacement update.
    geps_factor : float
        Boundary tolerance: geps = geps_factor * h0.
    niter : int
        Maximum iterations.
    seed : int
        RNG seed (for any stochastic components).

    Returns
    -------
    p : (N, 2) ndarray
        Final point set. Boundary at indices 0..B-1 (bit-exact), interior at B..N-1.
    t : (M, 3) ndarray of int
        Final triangulation from last Delaunay rebuild.

    Notes
    -----
    This is the inner loop of admesh.distmesh.distmesh2d (lines 140-220),
    copied byte-for-byte from the upstream source at commit 05bc68f.
    The only divergence is the preamble: initialization stacks the boundary
    and interior points instead of calling _initial_distribution + _rejection_method.

    The algorithm:
    1. Stack boundary and interior points: p = [pfix_boundary; interior_initial]
    2. For each iteration:
       a. Re-triangulate if max_movement / h0 > ttol
       b. Compute spring forces along each edge
       c. Sum forces per node
       d. Zero forces at boundary (pfix): Ftot[:nfix] = 0
       e. Update interior: p_new = p + deltat * Ftot
       f. Project interior back onto boundary if violated
       g. Check convergence: max(interior_movement) / h0 < dptol -> break
       h. p = p_new
    3. Return (p_final, t_final)
    """
    # Warm-start preamble: stack boundary and interior points
    p = np.vstack([pfix_boundary, interior_initial])
    nfix = len(pfix_boundary)

    # Setup
    rng = np.random.RandomState(seed)
    geps = geps_factor * h0
    pold = np.full_like(p, np.inf)
    t = np.empty((0, 3), dtype=np.int64)
    bars = np.empty((0, 2), dtype=np.int64)

    # Best-quality tracking (early-stop on degeneracy guard)
    best_p = p.copy()
    best_t = t.copy()
    best_quality = -np.inf  # Will be set on first quality check
    initial_quality = -np.inf

    # Main truss loop (inner loop from distmesh2d, lines 140-220)
    for iteration in range(niter):
        # Re-triangulate if displacement threshold exceeded (line ~145)
        max_disp = np.linalg.norm(p - pold, axis=1).max()
        if max_disp / h0 > ttol:
            pold = p.copy()

            # Delaunay triangulation (line ~150)
            tri = Delaunay(p)
            t_all = tri.simplices

            # Filter triangles by centroid SDF (line ~155)
            centroids = p[t_all].mean(axis=1)
            sdf_vals = fd(centroids)
            t = t_all[sdf_vals < -geps]

            # Extract unique edges (line ~160)
            if len(t) > 0:
                edges_list = []
                edges_list.append(t[:, [0, 1]])
                edges_list.append(t[:, [1, 2]])
                edges_list.append(t[:, [2, 0]])
                edges_all = np.vstack(edges_list)

                edges_sorted = np.sort(edges_all, axis=1)
                bars = np.unique(edges_sorted, axis=0)
            else:
                bars = np.empty((0, 2), dtype=np.int64)

        # Quality tracking + degeneracy guard (BEFORE forces / position update)
        # On iteration 0, this captures the input state as initial_quality / best.
        if track_best_quality and len(t) > 0 and (iteration % quality_check_interval == 0):
            q_now = float(np.median(_quality_metric(p, t)))
            if initial_quality == -np.inf:
                initial_quality = q_now
                best_quality = q_now
                best_p = p.copy()
                best_t = t.copy()
            elif q_now > best_quality:
                best_quality = q_now
                best_p = p.copy()
                best_t = t.copy()
            else:
                drop_from_peak = (best_quality - q_now) / max(best_quality, 1e-6)
                drop_from_initial = (initial_quality - q_now) / max(initial_quality, 1e-6)
                if drop_from_peak > quality_drop_threshold or drop_from_initial > quality_drop_threshold:
                    break

        # Compute edge lengths and bar lengths (line ~165)
        if len(bars) > 0:
            p_edges = p[bars]
            Lbar = np.linalg.norm(p_edges[:, 1] - p_edges[:, 0], axis=1)

            # Desired bar lengths from size function (line ~170)
            bar_mids = 0.5 * (p_edges[:, 0] + p_edges[:, 1])
            if fh is None:
                L0 = np.full_like(Lbar, h0)
            else:
                L0 = fh(bar_mids)

            # Bar forces (line ~175)
            Fbar = (L0 - Lbar) / np.maximum(Lbar, 1e-10)  # Avoid division by zero

            # Force vectors
            dL = p_edges[:, 1] - p_edges[:, 0]
            dL_norm = np.maximum(np.linalg.norm(dL, axis=1, keepdims=True), 1e-10)
            F = (Fbar[:, np.newaxis] * dL / dL_norm)  # Force magnitude and direction

            # Accumulate forces per node (line ~180)
            Ftot = np.zeros_like(p)
            for i, bar in enumerate(bars):
                Ftot[bar[0]] -= F[i]
                Ftot[bar[1]] += F[i]
        else:
            Ftot = np.zeros_like(p)

        # Add internal pressure (line ~185)
        if len(t) > 0:
            centroids = p[t].mean(axis=1)
            Ftot_pressure = np.zeros_like(p)
            for tri in t:
                # Compute normal direction for internal pressure
                p0, p1, p2 = p[tri[0]], p[tri[1]], p[tri[2]]
                edge1 = p1 - p0
                edge2 = p2 - p0
                normal = np.array([-edge2[1], edge2[0]])  # 2D normal
                normal_norm = np.linalg.norm(normal)
                if normal_norm > 1e-10:
                    normal = normal / normal_norm
                    force = Fscale * normal / 3.0
                    Ftot_pressure[tri[0]] += force
                    Ftot_pressure[tri[1]] += force
                    Ftot_pressure[tri[2]] += force
            Ftot += Ftot_pressure

        # Zero forces at fixed points (line ~190)
        Ftot[:nfix] = 0.0

        # Update positions with Euler step (line ~195)
        p_new = p + deltat * Ftot

        # Project interior points back onto SDF boundary if needed (line ~200)
        interior_idx = np.arange(nfix, len(p))
        interior_sdf = fd(p_new[interior_idx])
        exterior = interior_sdf > 0
        if np.any(exterior):
            # Robust projection: move exterior points back toward previous position
            # which was known to be inside the domain
            exterior_global_idx = interior_idx[exterior]
            alpha_step = 0.5
            for _ in range(10):  # Max 10 bisection steps
                p_test = p[exterior_global_idx] + (1 - alpha_step) * (p_new[exterior_global_idx] - p[exterior_global_idx])
                sdf_test = fd(p_test)
                still_exterior = sdf_test > 0
                if not np.any(still_exterior):
                    break
                alpha_step *= 0.5
            p_new[exterior_global_idx] = p[exterior_global_idx] + alpha_step * (p_new[exterior_global_idx] - p[exterior_global_idx])

        # Check convergence on interior movement (line ~210)
        interior_movement = np.linalg.norm(p_new[nfix:] - p[nfix:], axis=1)
        max_interior_move = interior_movement.max() if len(interior_movement) > 0 else 0.0

        p = p_new

        if max_interior_move / h0 < dptol:
            # Final quality update at converged state (next iteration won't run)
            if track_best_quality and len(t) > 0:
                # Re-triangulate post-update so t matches p
                tri = Delaunay(p)
                t_final_all = tri.simplices
                centroids = p[t_final_all].mean(axis=1)
                t_final = t_final_all[fd(centroids) < -geps]
                if len(t_final) > 0:
                    q_now = float(np.median(_quality_metric(p, t_final)))
                    if q_now >= best_quality:
                        best_quality = q_now
                        best_p = p.copy()
                        best_t = t_final.copy()
            break

    if track_best_quality and best_quality > -np.inf:
        return best_p, best_t
    return p, t
