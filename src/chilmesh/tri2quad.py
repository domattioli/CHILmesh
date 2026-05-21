"""Triangle-to-quad conversion (port of QuADMesh's Tri2QuadRoutine).

The MATLAB original walks skeletonization layers and pairs triangles
across every-other interior edge into quads. This Python port runs a
single global quality-weighted greedy matching on the dual graph: each
interior edge is scored by the convexity + angle quality of the quad
it would produce, sorted descending, and merged greedily.

No centroid bisection, no edge insertion — the output is **conforming**
(no T-junctions / hanging nodes). Boundary triangles (touching at least
one boundary vertex) that remain unmatched are kept as padded-tri rows
``[a, b, c, a]``. Interior unmatched triangles raise RuntimeError; this
input cannot be quadified without non-conforming operations.

The output satisfies CHILmesh spec 007:
  - quad-dominant with triangles only on the boundary layer
  - no self-intersecting quads
  - no element-element overlap
  - no hanging nodes (conforming)
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from chilmesh import CHILmesh


def tri_to_quad(mesh: "CHILmesh", *, strict: bool = True) -> "CHILmesh":
    """Convert a triangular mesh into a quad-dominant mesh.

    Args:
        mesh: CHILmesh with ``type == "Triangular"``. Pure-triangle input.
        strict: When True (default), raise RuntimeError if any interior
            triangle remains unmatched after maximum matching. When False,
            keep unmatched interior tris in the output as padded-tri rows;
            the result violates spec-007 (FR-006) but is still useful for
            inspection / visualization.

    Returns:
        New CHILmesh whose elements are quads (4-column) with boundary
        triangles encoded in the padded form ``[a, b, c, a]``.

    Raises:
        ValueError: If the input is not a pure-triangle mesh.
        RuntimeError: If ``strict`` and an interior triangle remains
            unmatched (would require non-conforming bisection to fix).
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

    edge_to_tris = _build_edge_map(tris)
    boundary_vert_ids = _boundary_vertex_ids(edge_to_tris)

    candidates: dict[tuple[int, int], tuple[float, tuple[int, int, int, int]]] = {}
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
        key = (t_a, t_b) if t_a < t_b else (t_b, t_a)
        candidates[key] = (score, quad_verts)

    try:
        import networkx as nx

        G = nx.Graph()
        for t_id in range(n_tris):
            G.add_node(t_id)
        for (a, b), (score, _) in candidates.items():
            G.add_edge(a, b, weight=score)
        matching = nx.max_weight_matching(G, maxcardinality=True)
        matched_partner: dict[int, int] = {}
        matched_quad: dict[int, tuple[int, int, int, int]] = {}
        used: set[int] = set()
        for a, b in matching:
            key = (a, b) if a < b else (b, a)
            score, quad_verts = candidates[key]
            t_a = a if a < b else b
            t_b = b if a < b else a
            matched_partner[t_a] = t_b
            matched_quad[t_a] = quad_verts
            used.add(a)
            used.add(b)
    except ImportError:
        ordered = sorted(
            candidates.items(),
            key=lambda kv: kv[1][0],
            reverse=True,
        )
        matched_partner = {}
        matched_quad = {}
        used = set()
        for (a, b), (_score, quad_verts) in ordered:
            if a in used or b in used:
                continue
            used.add(a)
            used.add(b)
            t_a = a if a < b else b
            t_b = b if a < b else a
            matched_partner[t_a] = t_b
            matched_quad[t_a] = quad_verts

    new_rows: list[list[int]] = []
    for t_a, quad_verts in matched_quad.items():
        new_rows.append(
            [int(quad_verts[0]), int(quad_verts[1]),
             int(quad_verts[2]), int(quad_verts[3])]
        )

    interior_unmatched: list[int] = []
    for t_id in range(n_tris):
        if t_id in used:
            continue
        verts = tris[t_id].tolist()
        if any(v in boundary_vert_ids for v in verts):
            new_rows.append(
                [int(verts[0]), int(verts[1]), int(verts[2]), int(verts[0])]
            )
        else:
            interior_unmatched.append(t_id)

    if interior_unmatched:
        if strict:
            raise RuntimeError(
                f"tri_to_quad left {len(interior_unmatched)} interior triangles "
                f"unmatched (first 10 elem IDs: {interior_unmatched[:10]}). "
                f"This input is not convertible by the current algorithm without "
                f"non-conforming bisection."
            )
        for t_id in interior_unmatched:
            verts = tris[t_id].tolist()
            new_rows.append(
                [int(verts[0]), int(verts[1]), int(verts[2]), int(verts[0])]
            )

    new_conn = np.array(new_rows, dtype=int)
    return CHILmesh(
        connectivity=new_conn,
        points=points.copy(),
        grid_name=mesh.grid_name,
        compute_layers=True,
    )



def _build_edge_map(tris: np.ndarray) -> dict[tuple[int, int], list[int]]:
    edge_to_tris: dict[tuple[int, int], list[int]] = {}
    for t_id in range(tris.shape[0]):
        v0, v1, v2 = int(tris[t_id, 0]), int(tris[t_id, 1]), int(tris[t_id, 2])
        for a, b in ((v0, v1), (v1, v2), (v2, v0)):
            key = (a, b) if a < b else (b, a)
            edge_to_tris.setdefault(key, []).append(t_id)
    return edge_to_tris


def _boundary_vertex_ids(
    edge_to_tris: dict[tuple[int, int], list[int]],
) -> set[int]:
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
    tri_a_l = [int(v) for v in tri_a]
    tri_b_l = [int(v) for v in tri_b]
    unique_a = [v for v in tri_a_l if v != shared_v0 and v != shared_v1]
    unique_b = [v for v in tri_b_l if v != shared_v0 and v != shared_v1]
    if len(unique_a) != 1 or len(unique_b) != 1:
        return None
    ua = unique_a[0]
    ub = unique_b[0]

    idx_ua = tri_a_l.index(ua)
    next_after_ua = tri_a_l[(idx_ua + 1) % 3]
    if next_after_ua == shared_v0:
        return (ua, shared_v0, ub, shared_v1)
    return (ua, shared_v1, ub, shared_v0)


def _quad_quality(quad_xy: np.ndarray) -> float:
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
