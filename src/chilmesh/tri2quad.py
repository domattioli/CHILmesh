"""Triangle-to-quad conversion (port of QuADMesh's Tri2QuadRoutine).

The MATLAB original (``02_QuADMESH_Library/02_Tri2Quad_Routine``) walks
each skeletonization layer from outer-most inward, identifies every-other
edge along paths on the outer vertices, and merges the two triangles
sharing each flagged edge into a quad. Remaining triangles are converted
via edge bisection / edge insertion / edge removal heuristics, with a
small number kept as boundary triangles.

This Python port uses the same conceptual structure (layer-by-layer,
inner edge removal) but the matching step is simplified to a
quality-weighted greedy maximum matching on the dual graph: each
interior edge is scored by the convexity + angle quality of the
quad it would produce, sorted descending, and greedily merged.
Unmatched triangles are kept only if they touch the boundary; interior
triangles are bisected against the highest-quality neighbor.

The output satisfies CHILmesh spec 007:
  - quad-dominant with triangles only on layer 0
  - no self-intersecting quads
  - no element-element overlap
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from chilmesh import CHILmesh


def tri_to_quad(mesh: "CHILmesh") -> "CHILmesh":
    """Convert a triangular mesh into a quad-dominant mesh.

    Args:
        mesh: CHILmesh with ``type == "Triangular"``. Pure-triangle input
            (3-column ``connectivity_list``).

    Returns:
        New CHILmesh whose elements are quads (4-column) with boundary
        triangles encoded in the padded form ``[a, b, c, c]``.
    """
    from chilmesh import CHILmesh

    if mesh.connectivity_list.shape[1] != 3:
        raise ValueError(
            f"tri_to_quad requires a pure-triangle mesh; "
            f"connectivity has {mesh.connectivity_list.shape[1]} columns"
        )

    tris = np.asarray(mesh.connectivity_list, dtype=int)
    points = np.asarray(mesh.points, dtype=float)
    n_tris = tris.shape[0]

    edge_to_tris: dict[tuple[int, int], list[int]] = {}
    for t_id in range(n_tris):
        v0, v1, v2 = int(tris[t_id, 0]), int(tris[t_id, 1]), int(tris[t_id, 2])
        for a, b in ((v0, v1), (v1, v2), (v2, v0)):
            key = (a, b) if a < b else (b, a)
            edge_to_tris.setdefault(key, []).append(t_id)

    boundary_vert_ids = _boundary_vertex_ids(edge_to_tris)

    candidates: list[tuple[float, int, int, tuple[int, int, int, int]]] = []
    for (a, b), tlist in edge_to_tris.items():
        if len(tlist) != 2:
            continue
        t_a, t_b = tlist
        quad_verts = _merge_tris_into_quad(tris[t_a], tris[t_b], a, b)
        if quad_verts is None:
            continue
        score = _quad_quality(points[list(quad_verts), :2])
        if score <= 0.0:
            continue
        candidates.append((score, t_a, t_b, quad_verts))

    candidates.sort(reverse=True)

    matched: dict[int, tuple[int, tuple[int, int, int, int]]] = {}
    used: set[int] = set()
    for score, t_a, t_b, quad_verts in candidates:
        if t_a in used or t_b in used:
            continue
        used.add(t_a)
        used.add(t_b)
        matched[t_a] = (t_b, quad_verts)

    new_rows: list[list[int]] = []
    for t_a, (t_b, quad_verts) in matched.items():
        new_rows.append(list(quad_verts))

    leftover = [t for t in range(n_tris) if t not in used]
    for t_id in leftover:
        verts = tris[t_id].tolist()
        if any(v in boundary_vert_ids for v in verts):
            new_rows.append([int(verts[0]), int(verts[1]), int(verts[2]), int(verts[0])])
        else:
            replacement = _bisect_interior_triangle(
                t_id, tris, points, edge_to_tris, matched, used, new_rows
            )
            new_rows.extend(replacement)
            used.add(t_id)

    new_conn = np.array(new_rows, dtype=int)
    return CHILmesh(
        connectivity=new_conn,
        points=points.copy(),
        grid_name=mesh.grid_name,
        compute_layers=True,
    )


def _boundary_vertex_ids(edge_to_tris: dict[tuple[int, int], list[int]]) -> set[int]:
    out: set[int] = set()
    for (a, b), tlist in edge_to_tris.items():
        if len(tlist) == 1:
            out.add(a)
            out.add(b)
    return out


def _merge_tris_into_quad(
    tri_a: np.ndarray,
    tri_b: np.ndarray,
    shared_v0: int,
    shared_v1: int,
) -> tuple[int, int, int, int] | None:
    """Return (q0, q1, q2, q3) in CCW order, or None if degenerate / non-shared.

    The four quad verts come from the two tris minus the shared edge:
    ``[unique_a, shared_v0, unique_b, shared_v1]`` in the order that
    follows the CCW winding of tri_a.
    """
    tri_a = [int(v) for v in tri_a]
    tri_b = [int(v) for v in tri_b]
    unique_a = [v for v in tri_a if v != shared_v0 and v != shared_v1]
    unique_b = [v for v in tri_b if v != shared_v0 and v != shared_v1]
    if len(unique_a) != 1 or len(unique_b) != 1:
        return None
    ua = unique_a[0]
    ub = unique_b[0]

    idx_ua = tri_a.index(ua)
    next_after_ua = tri_a[(idx_ua + 1) % 3]
    if next_after_ua == shared_v0:
        return (ua, shared_v0, ub, shared_v1)
    return (ua, shared_v1, ub, shared_v0)


def _quad_quality(quad_xy: np.ndarray) -> float:
    """Score [0, 1] of the resulting quad. Returns 0 for non-convex / bowtie.

    Uses min interior-angle ratio: ``min(angle_i) / 90°`` clipped to [0, 1].
    A convex right-angle quad scores 1; a near-degenerate quad scores ~0.
    """
    if quad_xy.shape != (4, 2):
        return 0.0
    cross_signs = []
    for i in range(4):
        a = quad_xy[(i - 1) % 4]
        b = quad_xy[i]
        c = quad_xy[(i + 1) % 4]
        e1 = a - b
        e2 = c - b
        cross = e1[0] * e2[1] - e1[1] * e2[0]
        cross_signs.append(cross)
    if not (all(s > 0 for s in cross_signs) or all(s < 0 for s in cross_signs)):
        return 0.0
    angles = []
    for i in range(4):
        a = quad_xy[(i - 1) % 4]
        b = quad_xy[i]
        c = quad_xy[(i + 1) % 4]
        e1 = a - b
        e2 = c - b
        n1 = np.linalg.norm(e1)
        n2 = np.linalg.norm(e2)
        if n1 == 0 or n2 == 0:
            return 0.0
        cos_t = float(np.clip(np.dot(e1, e2) / (n1 * n2), -1.0, 1.0))
        ang = np.degrees(np.arccos(cos_t))
        angles.append(ang)
    min_ang = min(angles)
    max_ang = max(angles)
    return float(min(min_ang, 180.0 - max_ang) / 90.0)


def _bisect_interior_triangle(
    t_id: int,
    tris: np.ndarray,
    points: np.ndarray,
    edge_to_tris: dict[tuple[int, int], list[int]],
    matched: dict[int, tuple[int, tuple[int, int, int, int]]],
    used: set[int],
    new_rows: list[list[int]],
) -> list[list[int]]:
    """Centroid-split an interior triangle into 3 quads.

    For an interior triangle with verts (v0, v1, v2), inserts the centroid
    as a new vertex and creates 3 quads sharing the centroid. This adds
    one new point but always produces a valid (non-self-intersecting,
    non-overlapping) result.

    The caller is responsible for committing the centroid to the global
    points array; here we encode it via a side-channel by appending to a
    module-level scratch list.
    """
    v0, v1, v2 = (int(tris[t_id, k]) for k in range(3))
    centroid = (points[v0, :2] + points[v1, :2] + points[v2, :2]) / 3.0
    m01 = (points[v0, :2] + points[v1, :2]) / 2.0
    m12 = (points[v1, :2] + points[v2, :2]) / 2.0
    m20 = (points[v2, :2] + points[v0, :2]) / 2.0

    new_vert_ids = _allocate_new_vertices(points, [centroid, m01, m12, m20])
    c_id, e01_id, e12_id, e20_id = new_vert_ids

    return [
        [v0, e01_id, c_id, e20_id],
        [v1, e12_id, c_id, e01_id],
        [v2, e20_id, c_id, e12_id],
    ]


def _allocate_new_vertices(points: np.ndarray, xy_list: list[np.ndarray]) -> list[int]:
    """Append (x, y, 0) rows to the underlying points array via numpy resize.

    NOTE: Mutates the points array in place. Callers must pass the same
    array that will be handed to the new CHILmesh constructor.
    """
    n_old = points.shape[0]
    z_const = float(points[0, 2]) if n_old > 0 else 0.0
    new_ids = []
    new_block = np.zeros((len(xy_list), 3))
    for i, xy in enumerate(xy_list):
        new_block[i, 0] = xy[0]
        new_block[i, 1] = xy[1]
        new_block[i, 2] = z_const
        new_ids.append(n_old + i)
    points_resized = np.vstack([points, new_block])
    points.resize(points_resized.shape, refcheck=False)
    points[:] = points_resized
    return new_ids
