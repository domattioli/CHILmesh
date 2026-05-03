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
import warnings
from typing import Callable, Optional, Tuple
from . import CHILmesh
from ._vendor_admesh_truss import distmesh2d_warmstart


def _infer_boundary_from_triangles(triangles: np.ndarray) -> np.ndarray:
    """
    Infer boundary vertex indices from triangle connectivity.

    A boundary edge is one that appears in exactly one triangle.
    The boundary vertices are all vertices incident to boundary edges.

    Parameters
    ----------
    triangles : (M, 3) ndarray of int
        Triangle vertex indices.

    Returns
    -------
    boundary_indices : (B,) ndarray of int
        Unique boundary vertex indices, sorted deterministically via np.unique.

    Notes
    -----
    This uses a hash-based edge counting approach for O(M) time complexity
    and deterministic ordering (per FR-007, T005).
    """
    edges = np.vstack([
        triangles[:, [0, 1]],
        triangles[:, [1, 2]],
        triangles[:, [2, 0]],
    ])

    # Sort edge endpoints so (i, j) and (j, i) become the same canonical form
    edges_sorted = np.sort(edges, axis=1)
    edges_tuples = [tuple(e) for e in edges_sorted]

    edge_count = {}
    for e in edges_tuples:
        edge_count[e] = edge_count.get(e, 0) + 1

    boundary_vertices = set()
    for e, count in edge_count.items():
        if count == 1:  # Boundary edge appears in exactly one triangle
            boundary_vertices.add(e[0])
            boundary_vertices.add(e[1])

    return np.unique(np.array(sorted(boundary_vertices), dtype=np.int32))


def _compute_h0_from_edges(points: np.ndarray, triangles: np.ndarray) -> float:
    """
    Estimate target edge length as mean of all edge lengths.

    Parameters
    ----------
    points : (N, 2|3) ndarray
        Vertex coordinates.
    triangles : (M, 3) ndarray of int
        Triangle connectivity.

    Returns
    -------
    h0 : float
        Mean edge length across all unique edges.
    """
    edges = np.vstack([
        triangles[:, [0, 1]],
        triangles[:, [1, 2]],
        triangles[:, [2, 0]],
    ])

    edges_sorted = np.sort(edges, axis=1)
    unique_edges = np.unique(edges_sorted, axis=0)

    edge_lengths = np.linalg.norm(
        points[unique_edges[:, 1], :2] - points[unique_edges[:, 0], :2],
        axis=1
    )

    return float(np.mean(edge_lengths))


def _infer_bbox_from_points(points: np.ndarray, margin: float = 0.05) -> Tuple[float, float, float, float]:
    """
    Infer axis-aligned bounding box from point cloud with margin.

    Parameters
    ----------
    points : (N, 2|3) ndarray
        Vertex coordinates (uses first two columns).
    margin : float
        Margin as fraction of box size (default 5%).

    Returns
    -------
    bbox : (xmin, ymin, xmax, ymax)
        Bounding box.
    """
    xy = points[:, :2]
    xmin, ymin = xy.min(axis=0)
    xmax, ymax = xy.max(axis=0)

    dx = xmax - xmin
    dy = ymax - ymin

    xmin -= margin * dx
    xmax += margin * dx
    ymin -= margin * dy
    ymax += margin * dy

    return (float(xmin), float(ymin), float(xmax), float(ymax))


def _signed_area(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray) -> np.ndarray:
    """
    Compute signed area of triangles using cross product.

    Parameters
    ----------
    p0, p1, p2 : (N, 2) ndarray
        Three vertices of N triangles.

    Returns
    -------
    areas : (N,) ndarray
        Signed areas (positive = CCW, negative = CW).
    """
    return 0.5 * ((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) -
                  (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1]))


def _validate_input(
    points: np.ndarray,
    triangles: np.ndarray,
    sdf: Callable[[np.ndarray], np.ndarray],
    boundary_indices: Optional[np.ndarray],
    bbox: Optional[Tuple[float, float, float, float]],
    sdf_tolerance: float = 1e-6,
) -> Tuple[np.ndarray, Tuple[float, float, float, float]]:
    """
    Validate input and return normalized boundary indices and bbox.

    Validation order (deterministic, per contracts/api-contract.md):
    1. V_TRI: Triangle-only check
    2. V_AREA: Positive signed area check
    3. Infer boundary_indices if None
    4. V_BND_SDF: Boundary on SDF zero-set check
    5. Infer bbox if None; verify all points within bbox
    6. V_BBOX: Points within bbox check

    Parameters
    ----------
    points : (N, 2|3) ndarray
        Vertex coordinates.
    triangles : (M, 3) ndarray of int
        Triangle connectivity.
    sdf : callable
        Signed distance function.
    boundary_indices : (B,) ndarray of int or None
        Pre-computed boundary indices (optional).
    bbox : (xmin, ymin, xmax, ymax) or None
        Pre-computed bounding box (optional).
    sdf_tolerance : float
        Tolerance for SDF zero-set check.

    Returns
    -------
    boundary_indices : (B,) ndarray of int
        Validated boundary indices.
    bbox : (xmin, ymin, xmax, ymax)
        Validated bounding box.

    Raises
    ------
    NotImplementedError
        If triangles are not triangle-only (V_TRI).
    ValueError
        If triangles have non-positive area (V_AREA).
        If boundary points are not on SDF zero-set (V_BND_SDF).
        If points are outside bbox (V_BBOX).
    RuntimeWarning (emitted, not raised)
        If interior points have sdf > 0 (outside domain).
        If interior point count < 10.
    """
    # V_TRI: Triangle-only check (line 1)
    if triangles.shape[1] != 3:
        raise NotImplementedError(
            f"warm-start truss requires triangle-only mesh; got {triangles.shape[1]}-vertex elements"
        )

    # V_AREA: Positive area check (line 2)
    p0 = points[triangles[:, 0], :2]
    p1 = points[triangles[:, 1], :2]
    p2 = points[triangles[:, 2], :2]
    areas = _signed_area(p0, p1, p2)

    bad_area_indices = np.where(areas <= 0)[0]
    if len(bad_area_indices) > 0:
        indices_str = str(bad_area_indices[:10].tolist()) + (
            "..." if len(bad_area_indices) > 10 else ""
        )
        raise ValueError(
            f"input has {len(bad_area_indices)} triangles with non-positive signed area: "
            f"indices {indices_str}"
        )

    # Infer boundary if None (line 3)
    if boundary_indices is None:
        boundary_indices = _infer_boundary_from_triangles(triangles)

    # V_BND_SDF: Boundary on SDF zero-set (line 4)
    boundary_xy = points[boundary_indices, :2]
    boundary_sdf = np.abs(sdf(boundary_xy))
    bad_bnd_indices = np.where(boundary_sdf >= sdf_tolerance)[0]

    if len(bad_bnd_indices) > 0:
        worst_idx = bad_bnd_indices[np.argmax(boundary_sdf[bad_bnd_indices])]
        worst_sdf = boundary_sdf[worst_idx]
        raise ValueError(
            f"boundary not on SDF zero set: {len(bad_bnd_indices)} of {len(boundary_indices)} "
            f"boundary points fail; max |sdf| = {worst_sdf:.2e} at index {boundary_indices[worst_idx]}; "
            f"tolerance = {sdf_tolerance:.2e}"
        )

    # Infer bbox if None (line 5)
    if bbox is None:
        bbox = _infer_bbox_from_points(points)

    # V_BBOX: Points within bbox (line 5 cont.)
    xy = points[:, :2]
    xmin, ymin, xmax, ymax = bbox

    out_of_bounds = (
        (xy[:, 0] < xmin) | (xy[:, 0] > xmax) |
        (xy[:, 1] < ymin) | (xy[:, 1] > ymax)
    )

    bad_bbox_indices = np.where(out_of_bounds)[0]
    if len(bad_bbox_indices) > 0:
        indices_str = str(bad_bbox_indices[:10].tolist()) + (
            "..." if len(bad_bbox_indices) > 10 else ""
        )
        raise ValueError(
            f"input has {len(bad_bbox_indices)} points outside bbox {bbox}: indices {indices_str}"
        )

    # Warnings (line 7-8)
    interior_indices = np.setdiff1d(np.arange(len(points)), boundary_indices)
    interior_xy = points[interior_indices, :2]

    interior_sdf = sdf(interior_xy)
    interior_outside = interior_sdf > 0
    if np.any(interior_outside):
        n_outside = np.sum(interior_outside)
        warnings.warn(
            f"{n_outside} interior points have sdf > 0 (outside domain); "
            "truss may push them onto boundary or stall",
            RuntimeWarning,
            stacklevel=3
        )

    if len(interior_indices) < 10:
        warnings.warn(
            f"warm-start truss called with only {len(interior_indices)} interior points "
            "(threshold: 10); optimization may have minimal effect",
            RuntimeWarning,
            stacklevel=3
        )

    return boundary_indices, bbox


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
    # T006-T008: Validate and normalize inputs
    boundary_indices, bbox = _validate_input(
        points, triangles, sdf, boundary_indices, bbox, sdf_tolerance
    )

    # T007: Compute h0 if None
    if h0 is None:
        h0 = _compute_h0_from_edges(points, triangles)

    # T010: Prepare normalized ValidatedInput
    points_xy = points[:, :2].astype(np.float64, copy=True)
    interior_indices = np.setdiff1d(np.arange(len(points)), boundary_indices)
    boundary_xy = points_xy[boundary_indices]
    interior_xy = points_xy[interior_indices]

    # T011: Call vendored truss loop
    try:
        points_opt, triangles_opt = distmesh2d_warmstart(
            boundary_xy, interior_xy,
            sdf, size_fn, h0, bbox,
            dptol=dptol, ttol=ttol, Fscale=Fscale, deltat=deltat,
            geps_factor=geps_factor, niter=niter, seed=seed
        )
    except ImportError as e:
        raise ImportError(
            "ADMESH (admesh.distmesh) is required for warm-start truss; "
            "pinned to SHA 05bc68f. Install via 'pip install -e /path/to/admesh' "
            "or 'pip install admesh2D'."
        ) from e

    # T012: Non-degradation guard (FR-011)
    if enforce_non_degradation:
        # Compute input and output quality metrics
        p0_in = points[triangles[:, 0], :2]
        p1_in = points[triangles[:, 1], :2]
        p2_in = points[triangles[:, 2], :2]
        areas_in = _signed_area(p0_in, p1_in, p2_in)
        quality_in = areas_in / np.linalg.norm(
            p1_in - p0_in, axis=1
        ) ** 2  # Simplified quality metric

        p0_out = points_opt[triangles_opt[:, 0], :2]
        p1_out = points_opt[triangles_opt[:, 1], :2]
        p2_out = points_opt[triangles_opt[:, 2], :2]
        areas_out = _signed_area(p0_out, p1_out, p2_out)
        quality_out = areas_out / np.linalg.norm(
            p1_out - p0_out, axis=1
        ) ** 2

        median_q_in = np.median(quality_in)
        median_q_out = np.median(quality_out)

        if median_q_out < median_q_in:
            warnings.warn(
                f"warm-start truss did not improve quality "
                f"(input median={median_q_in:.4f}, output median={median_q_out:.4f}); "
                "returning input unchanged (set enforce_non_degradation=False to override)",
                RuntimeWarning,
                stacklevel=2
            )
            return points[:, :2].copy(), triangles.copy()

    # T013: Return boundary-preserving output with bit-exact guarantee
    # Verify boundary preservation (sanity check)
    boundary_out = points_opt[:len(boundary_indices)]
    assert np.array_equal(boundary_out, boundary_xy), \
        "Internal error: boundary not preserved by distmesh2d_warmstart"

    return points_opt, triangles_opt


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
    # T019: Check triangle-only (FR-001b, FR-006)
    if mesh.connectivity_list.shape[1] != 3:
        raise NotImplementedError(
            "high-level wrapper requires triangle-only CHILmesh; use "
            "the array form with manual filtering for mixed meshes"
        )

    # T019: Extract boundary indices (FR-002)
    boundary_edge_ids = mesh.boundary_edges()
    edge2vert = mesh.adjacencies["Edge2Vert"]
    boundary_indices = np.unique(edge2vert[boundary_edge_ids].flatten())

    # T019: Extract points and triangles
    points = mesh.points
    triangles = mesh.connectivity_list

    # T019: Call array form
    points_opt, triangles_opt = optimize_with_admesh_truss_arrays(
        points, triangles, sdf, size_fn,
        boundary_indices=boundary_indices,
        **kwargs
    )

    # T020: Rewrap in fresh CHILmesh with 3D points (z=0)
    points_3d = np.column_stack([points_opt, np.zeros(len(points_opt))])
    mesh_out = CHILmesh(connectivity=triangles_opt, points=points_3d)

    return mesh_out
