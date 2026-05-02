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
    raise NotImplementedError("T015-T018 implementation pending")
