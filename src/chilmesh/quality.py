"""Standalone element quality computation without mesh object construction.

This module provides raw-array functions to compute per-element quality metrics
without requiring a full CHILmesh object. Useful for quality assessment of
element arrays before mesh construction or in standalone analysis.
"""
from __future__ import annotations

import numpy as np
from typing import Union


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
    if metric not in ("aspect_ratio", "min_angle", "max_angle", "skew", "skewness", "angular_skewness", "angular skewness"):
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

            if metric in ("skew", "skewness", "angular_skewness", "angular skewness"):
                # Use dedicated quad skew quality function
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
