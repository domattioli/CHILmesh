"""Standalone element quality computation without mesh object construction.

This module provides raw-array functions to compute per-element quality metrics
without requiring a full CHILmesh object. Useful for quality assessment of
element arrays before mesh construction or in standalone analysis.
"""
from __future__ import annotations

import numpy as np
from typing import Union

from .geometry import edge_lengths

GRAVITY = 9.81  # Gravitational acceleration (m/s^2)


def _triangle_quality(
    v0: np.ndarray,
    v1: np.ndarray,
    v2: np.ndarray,
    metric: str = "aspect_ratio",
) -> float:
    """Compute quality of a single triangle.

    Parameters
    ----------
    v0, v1, v2 : np.ndarray
        2D vertex coordinates (length 2).
    metric : str
        Quality metric: 'aspect_ratio', 'min_angle', 'max_angle', 'skew', 'skewness',
        'angular_skewness', or 'angular skewness'. Both underscore and space variants accepted.

    Returns
    -------
    float
        Quality value. For aspect_ratio: [0, 1] where 1 is ideal.
        For angles: raw radian values.
        For skew: [0, 1] where 1 is ideal (equilateral).
    """
    a = v1 - v0
    b = v2 - v1
    c = v0 - v2

    la = np.linalg.norm(a)
    lb = np.linalg.norm(b)
    lc = np.linalg.norm(c)

    if metric == "aspect_ratio":
        s = (la + lb + lc) / 2
        d1 = v1 - v0
        d2 = v2 - v0
        area = abs(d1[0] * d2[1] - d1[1] * d2[0]) / 2

        # Degenerate cases: collinear, zero area, or invalid side length
        if s <= 0 or la <= 0 or lb <= 0 or lc <= 0 or area <= 0:
            return 0.0

        r_in = area / s
        r_circ = (la * lb * lc) / (4 * area) if area > 0 else 0

        return float(2 * r_in / r_circ) if r_circ > 0 else 0.0

    elif metric in ("skew", "skewness", "angular_skewness", "angular skewness"):
        # Compute angles in degrees (same as interior_angles method)
        angles_rad = np.array([
            np.arccos(np.clip(np.dot(-c, a) / (lc * la + 1e-300), -1, 1)),
            np.arccos(np.clip(np.dot(-a, b) / (la * lb + 1e-300), -1, 1)),
            np.arccos(np.clip(np.dot(-b, c) / (lb * lc + 1e-300), -1, 1)),
        ])
        angles_deg = np.degrees(angles_rad)

        # Check degenerate: angle sum <= 179.99
        if np.sum(angles_deg) <= 179.99:
            return 0.0

        angle_max = np.max(angles_deg)
        angle_min = np.min(angles_deg)

        # Skew formula for triangle: ideal angle is 60 degrees
        quality = 1.0 - np.maximum(
            (angle_max - 60.0) / 120.0,
            (60.0 - angle_min) / 60.0
        )
        return float(quality)

    elif metric in ("equiangle_skewness", "equiangle skewness", "eas"):
        # Equiangle skewness Qeas in [0,1]: 0 = ideal (equilateral, all 60deg),
        # 1 = fully degenerate. Standard ANSYS/Fluent definition; exact
        # complement of the "skew" quality metric (skew == 1 - Qeas) on
        # non-degenerate elements. Ideal angle theta_e = 60deg for a tri. (#217)
        angles_rad = np.array([
            np.arccos(np.clip(np.dot(-c, a) / (lc * la + 1e-300), -1, 1)),
            np.arccos(np.clip(np.dot(-a, b) / (la * lb + 1e-300), -1, 1)),
            np.arccos(np.clip(np.dot(-b, c) / (lb * lc + 1e-300), -1, 1)),
        ])
        angles_deg = np.degrees(angles_rad)
        if np.sum(angles_deg) <= 179.99:
            return 1.0
        angle_max = np.max(angles_deg)
        angle_min = np.min(angles_deg)
        skewness = np.maximum(
            (angle_max - 60.0) / 120.0,
            (60.0 - angle_min) / 60.0,
        )
        return float(skewness)

    else:  # min_angle or max_angle
        angles = np.array([
            np.arccos(np.clip(np.dot(-c, a) / (lc * la + 1e-300), -1, 1)),
            np.arccos(np.clip(np.dot(-a, b) / (la * lb + 1e-300), -1, 1)),
            np.arccos(np.clip(np.dot(-b, c) / (lb * lc + 1e-300), -1, 1)),
        ])
        return float(angles.min() if metric == "min_angle" else angles.max())


def _quad_quality(
    v0: np.ndarray,
    v1: np.ndarray,
    v2: np.ndarray,
    v3: np.ndarray,
    metric: str = "aspect_ratio",
) -> float:
    """Compute quality of a single quadrilateral.

    Parameters
    ----------
    v0, v1, v2, v3 : np.ndarray
        2D vertex coordinates (length 2), in order around the quad.
    metric : str
        Quality metric: 'aspect_ratio', 'min_angle', 'max_angle', 'skew',
        'skewness', 'angular_skewness', or 'angular skewness'. Both underscore and space variants accepted.

    Returns
    -------
    float
        Quality value.
    """
    if metric in ("skew", "skewness", "angular_skewness", "angular skewness"):
        # Compute interior angles in degrees for all 4 vertices
        angles_deg = np.zeros(4)
        verts = np.array([v0, v1, v2, v3])

        for j in range(4):
            v_curr = verts[j]
            v_next = verts[(j + 1) % 4]
            v_prev = verts[(j - 1) % 4]

            v_edge1 = v_next - v_curr
            v_edge2 = v_prev - v_curr

            n1 = np.linalg.norm(v_edge1) + 1e-12
            n2 = np.linalg.norm(v_edge2) + 1e-12

            dot = np.clip(np.dot(v_edge1 / n1, v_edge2 / n2), -1.0, 1.0)
            angles_deg[j] = np.degrees(np.arccos(dot))

        # Check degenerate: angle sum <= 359.99
        if np.sum(angles_deg) <= 359.99:
            return 0.0

        angle_max = np.max(angles_deg)
        angle_min = np.min(angles_deg)

        # Skew formula for quad: ideal angle is 90 degrees
        quality = 1.0 - np.maximum(
            (angle_max - 90.0) / 90.0,
            (90.0 - angle_min) / 90.0
        )
        return float(quality)
    elif metric in ("equiangle_skewness", "equiangle skewness", "eas"):
        # Equiangle skewness Qeas in [0,1] for a quad: 0 = ideal (all 90deg),
        # 1 = degenerate. theta_e = 90deg. Complement of "skew". (#217)
        angles_deg = np.zeros(4)
        verts = np.array([v0, v1, v2, v3])
        for j in range(4):
            v_curr = verts[j]
            v_next = verts[(j + 1) % 4]
            v_prev = verts[(j - 1) % 4]
            v_edge1 = v_next - v_curr
            v_edge2 = v_prev - v_curr
            n1 = np.linalg.norm(v_edge1) + 1e-12
            n2 = np.linalg.norm(v_edge2) + 1e-12
            dot = np.clip(np.dot(v_edge1 / n1, v_edge2 / n2), -1.0, 1.0)
            angles_deg[j] = np.degrees(np.arccos(dot))
        if np.sum(angles_deg) <= 359.99:
            return 1.0
        angle_max = np.max(angles_deg)
        angle_min = np.min(angles_deg)
        skewness = np.maximum(
            (angle_max - 90.0) / 90.0,
            (90.0 - angle_min) / 90.0,
        )
        return float(skewness)
    elif metric in ("min_angle", "max_angle"):
        # Quad angle metrics use the element's own four interior angles
        # (radians), same convention as the skew/eas branches above. The
        # prior triangle-split path returned min-of-sub-triangle-maxima,
        # under-reporting a genuinely obtuse quad (#260).
        verts = np.array([v0, v1, v2, v3])
        angles = np.zeros(4)
        for j in range(4):
            v_curr = verts[j]
            v_next = verts[(j + 1) % 4]
            v_prev = verts[(j - 1) % 4]
            v_edge1 = v_next - v_curr
            v_edge2 = v_prev - v_curr
            n1 = np.linalg.norm(v_edge1) + 1e-12
            n2 = np.linalg.norm(v_edge2) + 1e-12
            dot = np.clip(np.dot(v_edge1 / n1, v_edge2 / n2), -1.0, 1.0)
            angles[j] = np.arccos(dot)
        return float(angles.min() if metric == "min_angle" else angles.max())
    else:
        # For other metrics, split quad into two triangles and take min
        q1 = _triangle_quality(v0, v1, v2, metric)
        q2 = _triangle_quality(v0, v2, v3, metric)
        return float(min(q1, q2))


def element_quality(
    verts: np.ndarray,
    conn: Union[list, np.ndarray],
    metric: str = "aspect_ratio",
) -> np.ndarray:
    """Compute per-element quality scores without constructing a CHILmesh object.

    Computes quality metrics for triangular and quadrilateral elements given
    vertex coordinates and connectivity. Handles mixed tri/quad meshes with
    variable-length connectivity rows.

    Parameters
    ----------
    verts : np.ndarray
        Vertex coordinates, shape (n_verts, 2) or (n_verts, 3).
        If 3D, only the first two columns are used (z ignored).
    conn : list[list[int]] or np.ndarray
        Element connectivity. Either:
        - list of lists (variable-length, may contain padding like -1)
        - 2D array (fixed rows, rows may be padded)
        Each row is a list of vertex indices forming the element.
        Triangles: 3 indices. Quads: 4 indices.
        Padding (-1) is filtered out.
    metric : str
        Quality metric to compute. One of:
        - 'aspect_ratio' (default): 2 * inradius / circumradius, range [0, 1].
          1 = equilateral, 0 = degenerate.
        - 'min_angle': Minimum interior angle (radians).
        - 'max_angle': Maximum interior angle (radians).
        - 'skew', 'skewness', 'angular_skewness', 'angular skewness': Angular deviation metric.
          Both underscore and space variants accepted. Range [0, 1].
          1 = ideal angles (60° for tri, 90° for quad), 0 = degenerate.

    Returns
    -------
    np.ndarray
        Quality scores, shape (n_elems,). One value per element.

    Raises
    ------
    ValueError
        If metric is not one of the supported values.

    Examples
    --------
    >>> import numpy as np
    >>> from chilmesh import element_quality
    >>> # Equilateral triangle
    >>> h = np.sqrt(3) / 2
    >>> verts = np.array([[0, 0], [1, 0], [0.5, h]])
    >>> conn = [[0, 1, 2]]
    >>> q = element_quality(verts, conn)
    >>> print(q[0])  # Should be ~1.0
    0.9999...

    >>> # Collinear (degenerate)
    >>> verts = np.array([[0, 0], [1, 0], [2, 0]])
    >>> conn = [[0, 1, 2]]
    >>> q = element_quality(verts, conn)
    >>> print(q[0])  # Should be 0.0
    0.0
    """
    if metric not in ("aspect_ratio", "min_angle", "max_angle", "skew", "skewness", "angular_skewness", "angular skewness", "equiangle_skewness", "equiangle skewness", "eas"):
        raise ValueError(f"Unknown metric: {metric!r}")

    verts_array = np.asarray(verts, dtype=float)
    if verts_array.shape[1] > 2:
        verts_array = verts_array[:, :2]

    n_elems = len(conn)
    qualities = np.zeros(n_elems, dtype=float)

    for i, elem in enumerate(conn):
        # Filter out -1 padding
        elem_valid = [int(e) for e in elem if int(e) >= 0]

        # Detect if padded triangle: any repeated vertex among 6 pairs
        # Check rows[:,0]==rows[:,1] | [1]==[2] | [2]==[3] | [3]==[0] | [0]==[2] | [1]==[3]
        if len(elem) == 4:
            row = np.asarray(elem, dtype=int)
            is_padded_tri = (
                (row[0] == row[1])
                | (row[1] == row[2])
                | (row[2] == row[3])
                | (row[3] == row[0])
                | (row[0] == row[2])
                | (row[1] == row[3])
            )
            if is_padded_tri:
                # Treat as triangle using first 3 unique vertices
                elem_valid = [int(e) for e in elem[:3] if int(e) >= 0]

        if len(elem_valid) == 3:
            # Triangle
            v0 = verts_array[elem_valid[0]]
            v1 = verts_array[elem_valid[1]]
            v2 = verts_array[elem_valid[2]]
            qualities[i] = _triangle_quality(v0, v1, v2, metric)

        elif len(elem_valid) == 4:
            # Quad
            v0 = verts_array[elem_valid[0]]
            v1 = verts_array[elem_valid[1]]
            v2 = verts_array[elem_valid[2]]
            v3 = verts_array[elem_valid[3]]

            if metric in ("skew", "skewness", "angular_skewness", "angular skewness", "equiangle_skewness", "equiangle skewness", "eas", "min_angle", "max_angle"):
                # Dedicated quad function: skew/eas + raw-interior-angle
                # min/max (#260 — angle metrics must not triangle-split).
                qualities[i] = _quad_quality(v0, v1, v2, v3, metric)
            else:
                # For other metrics, split into two triangles and take minimum
                q1 = _triangle_quality(v0, v1, v2, metric)
                q2 = _triangle_quality(v0, v2, v3, metric)
                qualities[i] = min(q1, q2)

        else:
            # Degenerate: fewer than 3 vertices
            qualities[i] = 0.0

    return qualities


def courant_number(points, edges, depths, dt, crs="cartesian", gravity=GRAVITY):
    """Per-edge Courant number ``sqrt(g*h) * dt / dx`` (shallow-water celerity).

    Parameters
    ----------
    points : array-like, shape (N, 2)
        Node coordinates: (x, y) when crs="cartesian", (lon, lat) degrees
        when crs="spherical" (dx then in geodesic metres).
    edges : array-like of int, shape (E, 2)
        Node-index pairs (e.g. a CHILmesh ``Edge2Vert`` adjacency).
    depths : array-like, shape (N,)
        Positive-down nodal depths (fort.14 convention). Edge depth h is the
        deeper of the two endpoint depths.
    dt : float
        Timestep in seconds (e.g. fort.15 DTDP).
    crs : {"cartesian", "spherical"}
        Coordinate reference system. "cartesian": planar Euclidean in native
        units. "spherical": haversine great-circle metres.
    gravity : float
        Gravitational acceleration, m/s^2. Default 9.81.

    Returns
    -------
    ndarray, shape (E,)
        Courant number per edge. NaN where the edge is dry (h <= 0) or
        degenerate (dx <= 0).

    Raises
    ------
    ValueError
        If shapes are incompatible or dt <= 0.
    """
    pts = np.asarray(points, dtype=float)
    e = np.asarray(edges, dtype=int)
    d = np.asarray(depths, dtype=float)

    # Validate shapes
    if pts.ndim != 2 or pts.shape[1] < 2:
        raise ValueError(f"points must be (N, >=2) array, got shape {pts.shape}")
    if e.ndim != 2 or e.shape[1] != 2:
        raise ValueError(f"edges must be (E, 2) array, got shape {e.shape}")
    if d.ndim != 1:
        raise ValueError(f"depths must be 1-D array, got shape {d.shape}")
    if len(e) > 0 and (e.min() < 0 or e.max() >= len(d)):
        raise ValueError(f"edge indices out of bounds: min={e.min()}, max={e.max()}, n_depths={len(d)}")
    if dt <= 0:
        raise ValueError(f"dt must be positive, got {dt}")

    # Handle empty edges
    if len(e) == 0:
        return np.array([], dtype=float)

    # Compute edge depth: max of the two endpoints
    h = np.maximum(d[e[:, 0]], d[e[:, 1]])

    # Compute edge lengths using edge_lengths
    p1 = pts[e[:, 0], :2]
    p2 = pts[e[:, 1], :2]
    dx = np.atleast_1d(edge_lengths(p1, p2, crs=crs))

    # Initialize Courant array
    c = np.full(len(e), np.nan, dtype=float)

    # Compute Courant only for wet, non-degenerate edges
    wet = (h > 0) & (dx > 0)
    c[wet] = np.sqrt(gravity * h[wet]) * dt / dx[wet]

    return c


def cfl_gate(points, edges, depths, dt, courant_max=1.0, crs="cartesian", gravity=GRAVITY, max_offenders=25):
    """CFL flag (not a stability proof): which edges exceed ``courant_max``.

    Parameters
    ----------
    points : array-like, shape (N, 2)
        Node coordinates.
    edges : array-like of int, shape (E, 2)
        Edge connectivity.
    depths : array-like, shape (N,)
        Nodal depths.
    dt : float
        Timestep in seconds.
    courant_max : float
        Maximum allowable Courant number. Default 1.0.
    crs : {"cartesian", "spherical"}
        Coordinate reference system.
    gravity : float
        Gravitational acceleration, m/s^2.
    max_offenders : int
        Maximum number of offending edges to report. Default 25.

    Returns
    -------
    dict with keys:
        status : str
            "pass" if all wet edges satisfy C <= courant_max, else "fail".
        n_over : int
            Number of edges exceeding courant_max.
        n_wet : int
            Number of wet (non-dry, non-degenerate) edges.
        n_skipped : int
            Number of dry or degenerate edges.
        max : float
            Maximum Courant number over wet edges (0.0 if none).
        median : float
            Median Courant number over wet edges (0.0 if none).
        p95 : float
            95th percentile Courant number over wet edges (0.0 if none).
        offenders : list of dict
            Worst-first (by courant), capped at max_offenders. Each dict:
            {"edge_index": int, "nodes": [a, b], "courant": float,
             "edge_length": float, "depth": float}
        thresholds : dict
            {"courant_max": courant_max}
    """
    pts = np.asarray(points, dtype=float)
    e = np.asarray(edges, dtype=int)
    d = np.asarray(depths, dtype=float)

    # Compute Courant numbers
    c = courant_number(pts, e, d, dt, crs=crs, gravity=gravity)

    # Separate wet from skipped
    wet_mask = np.isfinite(c)
    c_wet = c[wet_mask]
    n_wet = int(np.sum(wet_mask))
    n_skipped = len(c) - n_wet

    # Compute stats over wet edges
    if n_wet > 0:
        c_max = float(np.max(c_wet))
        c_median = float(np.median(c_wet))
        c_p95 = float(np.percentile(c_wet, 95))
    else:
        c_max = c_median = c_p95 = 0.0

    # Find offenders
    offender_mask = c > courant_max
    n_over = int(np.sum(offender_mask))
    offender_indices = np.where(offender_mask)[0]

    # Sort offenders by Courant (descending)
    offender_indices = offender_indices[np.argsort(-c[offender_indices])][:max_offenders]

    # Build offender records
    offenders = []
    p1 = pts[e[:, 0], :2]
    p2 = pts[e[:, 1], :2]
    dx = np.atleast_1d(edge_lengths(p1, p2, crs=crs))
    h = np.maximum(d[e[:, 0]], d[e[:, 1]])

    for idx in offender_indices:
        offenders.append({
            "edge_index": int(idx),
            "nodes": [int(e[idx, 0]), int(e[idx, 1])],
            "courant": float(c[idx]),
            "edge_length": float(dx[idx]),
            "depth": float(h[idx])
        })

    return {
        "status": "pass" if n_over == 0 else "fail",
        "n_over": n_over,
        "n_wet": n_wet,
        "n_skipped": n_skipped,
        "max": c_max,
        "median": c_median,
        "p95": c_p95,
        "offenders": offenders,
        "thresholds": {"courant_max": courant_max}
    }
