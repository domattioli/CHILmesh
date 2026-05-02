"""
ADMESH Warm-Start Truss Optimizer Adapter

This module implements a warm-start wrapper around ADMESH's distmesh truss algorithm,
enabling optimization of existing triangulations while preserving boundary points bit-exactly.

**Functional Requirements** (FR-001a, FR-001b, FR-002 through FR-018):
- Two layered entry points: array form (FR-001a) and CHILmesh wrapper (FR-001b)
- Deterministic boundary identification (FR-002)
- Boundary pinning via pfix mechanism (FR-003)
- Warm-start initialization from existing points (FR-004)
- Pre-truss validation: SDF, triangle-only, positive areas (FR-005-007)
- Bit-exact boundary preservation (FR-008)
- RNG seed control for determinism (FR-009)
- Input immutability (FR-010)
- Non-degradation graceful fallback (FR-011)
- 4-row demo restructure (FR-012)
- Fail-loud demo assertions (FR-013)
- ADMESH import fallback (FR-014)
- ADMESH version pinning to 05bc68f (FR-015)
- Input-source agnosticism (FR-016)
- Domain agnosticism (FR-017)
- Extensibility examples (FR-018)

See also:
- contracts/api-contract.md: Function signatures, validation order, error catalog
- contracts/visualization-output.md: 4-row demo structure and assertions
- data-model.md: Runtime entities and state transitions
- quickstart.md: Worked examples for extensibility contract
- research.md: Design rationale and cross-repo tracking (ADMESH-A, B, C)
"""

import numpy as np
from typing import Callable, Optional, Tuple
from . import CHILmesh
from ._vendor_admesh_truss import distmesh2d_warmstart


def optimize_with_admesh_truss_arrays(
    points: np.ndarray,
    triangles: np.ndarray,
    sdf: Callable[[np.ndarray], np.ndarray],
    size_fn: Optional[Callable[[np.ndarray], np.ndarray]] = None,
    *,
    boundary_indices: Optional[np.ndarray] = None,
    h0: Optional[float] = None,
    bbox: Optional[Tuple[float, float, float, float]] = None,
    dptol: float = 1e-3,
    ttol: float = 0.1,
    Fscale: float = 1.2,
    deltat: float = 0.2,
    geps_factor: float = 1e-3,
    niter: int = 500,
    seed: int = 0,
    sdf_tolerance: float = 1e-6,
    enforce_non_degradation: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Optimize an existing triangulation using ADMESH's distmesh truss algorithm,
    preserving boundary points bit-exactly.

    This is the low-level, source-agnostic entry point. It accepts raw numpy arrays
    and returns optimized (points, triangles) with no dependency on CHILmesh.

    Parameters
    ----------
    points : (N, 2) or (N, 3) ndarray
        Vertex coordinates. Z column (if present) is preserved as zero in the
        output but never used by the truss.
    triangles : (M, 3) ndarray of int
        Triangle vertex indices. Quad / mixed-element inputs raise NotImplementedError.
    sdf : callable
        Signed distance function: `sdf(points: (K,2)) -> (K,)`. Negative inside.
    size_fn : callable or None
        Element size function: `size_fn(points: (K,2)) -> (K,)`. None means uniform (h0).
    boundary_indices : (B,) ndarray of int or None
        Indices of points to pin (preserved bit-exactly). If None, inferred from
        triangulation (edges in exactly one triangle).
    h0 : float or None
        Target edge length. If None, computed as the mean input edge length.
    bbox : (xmin, ymin, xmax, ymax) or None
        Bounding box. If None, inferred from input points with 5% margin.
    dptol, ttol, Fscale, deltat, geps_factor, niter, seed : truss-loop tunables
        Forwarded verbatim to the vendored distmesh2d_warmstart. See research.md R3
        for descriptions.
    sdf_tolerance : float
        Tolerance for boundary-on-SDF check (default 1e-6).
    enforce_non_degradation : bool
        If True (default), output median quality < input median quality causes
        the input to be returned unchanged plus a RuntimeWarning. If False,
        always return the truss-loop output regardless of quality.

    Returns
    -------
    points_out : (Nout, 2) ndarray
        Optimized vertex coordinates. Boundary points (the ones at indices
        `boundary_indices` in the input) appear at indices 0..B-1 of the output
        in the same order, byte-equal to the input boundary subset.
    triangles_out : (Mout, 3) ndarray of int
        Final triangulation from the truss loop's last Delaunay rebuild.

    Raises
    ------
    NotImplementedError
        If `triangles.shape[1] != 3`.
    ValueError
        - Input triangles have non-positive signed area.
        - Boundary points fail the SDF-zero-set check (|sdf(p)| >= sdf_tolerance).
        - Input points fall outside the supplied/inferred bbox.
    ImportError
        If ADMESH cannot be imported.
    RuntimeWarning (emitted, not raised)
        - Truss loop hit niter cap without convergence (returns best-so-far).
        - Output quality regressed and enforce_non_degradation=True (returns input).
        - Interior point count below threshold (optimization may have minimal effect).
        - Interior points have sdf > 0 (outside domain).
    """
    raise NotImplementedError("T009-T013 implementation pending")


def optimize_with_admesh_truss(
    mesh: CHILmesh,
    sdf: Callable[[np.ndarray], np.ndarray],
    size_fn: Optional[Callable[[np.ndarray], np.ndarray]] = None,
    **kwargs,
) -> CHILmesh:
    """
    CHILmesh-native wrapper around `optimize_with_admesh_truss_arrays`.

    Extracts (points, triangles, boundary_indices) from `mesh`, calls the
    array form, and rewraps the result in a fresh CHILmesh.

    All kwargs are forwarded to the array form. See `optimize_with_admesh_truss_arrays`
    for the full parameter list and behavior.

    Parameters
    ----------
    mesh : CHILmesh
        Input mesh. Must be triangle-only. Must have a non-empty boundary
        (i.e., be an open / non-closed manifold).
    sdf, size_fn, **kwargs
        See `optimize_with_admesh_truss_arrays`.

    Returns
    -------
    CHILmesh
        Fresh CHILmesh instance. Adjacencies and layers are recomputed from
        scratch. The input mesh is not modified.

    Raises
    ------
    All exceptions from `optimize_with_admesh_truss_arrays`, plus:
    NotImplementedError
        If `mesh` is mixed-element (has both triangles and quads). The error
        message says "high-level wrapper requires triangle-only CHILmesh; use
        the array form with manual filtering for mixed meshes".
    """
    raise NotImplementedError("T019-T021 implementation pending")
