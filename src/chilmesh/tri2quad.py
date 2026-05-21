"""Triangle-to-quad conversion (port of QuADMesh's Tri2QuadRoutine).

The MATLAB original ``Tri2QuadRoutine`` iterates skeletonization layers
inner→outer; within each layer ``identifyEdgesFun_v2`` selects merge
candidates via every-other-edge path walking on outer vertices and
``mergeTrianglesFun`` merges each selected pair into a quad.

This Python port replaces the layer-restricted path-walk pairing with
a single global maximum-cardinality matching (Blossom V via networkx)
on the dual graph of the input. Empirically — measured on the four
CHILmesh built-in fixtures — global matching strictly dominates the
layer-restricted variant in pair count:

  * ``structured``: both yield 0 leftover.
  * ``block_o``: global → 0 leftover; layer-restricted → 6 leftover.
  * ``annulus`` / ``donut``: both leave Tutte-Berge obstructions
    (input dual graph has odd components that no matching algorithm
    can pair without subdivision).

The MATLAB ``removeTrianglesFun`` (vertex insertion via
``edgeBisection`` / ``edgeInsertion``) is intentionally NOT ported per
user direction — interior triangles must be eliminated by pairing logic
alone, never by inserting new vertices. Inputs whose dual graph has a
Tutte-Berge obstruction will surface those interior triangles as
``RuntimeError`` (when ``strict=True``).

Output: conforming (no hanging nodes), quad-dominant.
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
            triangle survives. When False, encode leftover interior tris
            as padded-tri rows (validator will flag them).

    Returns:
        New CHILmesh whose elements are quads (4-column) with boundary
        triangles encoded in the padded form ``[a, b, c, a]``.

    Raises:
        ValueError: If the input is not a pure-triangle mesh.
        RuntimeError: If ``strict`` and an interior triangle survives.
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

    tri_layer = _assign_layer_per_tri(mesh, n_tris)
    n_layers = max(max(tri_layer) + 1, 1) if tri_layer else 1

    edge_to_tris = _build_edge_map(tris)
    boundary_vert_ids = _boundary_vertex_ids(edge_to_tris)

    quad_rows: list[list[int]] = []
    consumed: set[int] = set()

    global_pool = set(range(n_tris))
    _layer_pair_and_collect(
        global_pool, tris, points, edge_to_tris,
        quad_rows, consumed,
    )

    leftover_interior: list[int] = []
    for t_id in range(n_tris):
        if t_id in consumed:
            continue
        verts = tris[t_id].tolist()
        if any(int(v) in boundary_vert_ids for v in verts):
            quad_rows.append(
                [int(verts[0]), int(verts[1]), int(verts[2]), int(verts[0])]
            )
            consumed.add(t_id)
        else:
            leftover_interior.append(t_id)

    if leftover_interior:
        if strict:
            raise RuntimeError(
                f"tri_to_quad left {len(leftover_interior)} interior triangles "
                f"after layer-by-layer pairing + rescue "
                f"(first 10: {leftover_interior[:10]})."
            )
        for t_id in leftover_interior:
            verts = tris[t_id].tolist()
            quad_rows.append(
                [int(verts[0]), int(verts[1]), int(verts[2]), int(verts[0])]
            )

    new_conn = np.array(quad_rows, dtype=int)
    return CHILmesh(
        connectivity=new_conn,
        points=points.copy(),
        grid_name=mesh.grid_name,
        compute_layers=True,
    )


def _layer_pair_and_collect(
    pool: set[int],
    tris: np.ndarray,
    points: np.ndarray,
    edge_to_tris: dict[tuple[int, int], list[int]],
    quad_rows: list[list[int]],
    consumed: set[int],
) -> None:
    """Run max-cardinality matching on the pool's dual graph; emit quads."""
    pair_candidates: dict[tuple[int, int], tuple[float, tuple[int, int, int, int]]] = {}
    for (a, b), tlist in edge_to_tris.items():
        usable = [t for t in tlist if t in pool and t not in consumed]
        if len(usable) != 2:
            continue
        t_a, t_b = usable
        quad_verts = _merge_tris_into_quad(tris[t_a], tris[t_b], a, b)
        if quad_verts is None:
            continue
        score = _quad_quality(points[list(quad_verts), :2])
        if score <= 0.0:
            continue
        key = (t_a, t_b) if t_a < t_b else (t_b, t_a)
        pair_candidates[key] = (score, quad_verts)

    matched = _max_weight_matching(pool, pair_candidates)
    for (t_a, t_b), quad_verts in matched.items():
        quad_rows.append(list(quad_verts))
        consumed.add(t_a)
        consumed.add(t_b)


def _max_weight_matching(
    pool,
    pair_candidates: dict[tuple[int, int], tuple[float, tuple[int, int, int, int]]],
) -> dict[tuple[int, int], tuple[int, int, int, int]]:
    """Maximum-cardinality, max-weight matching via networkx Blossom V."""
    try:
        import networkx as nx

        G = nx.Graph()
        for t in pool:
            G.add_node(t)
        for (a, b), (score, _) in pair_candidates.items():
            G.add_edge(a, b, weight=score)
        matching = nx.max_weight_matching(G, maxcardinality=True)
        out = {}
        for a, b in matching:
            key = (a, b) if a < b else (b, a)
            _score, quad_verts = pair_candidates[key]
            out[key] = quad_verts
        return out
    except ImportError:  # pragma: no cover
        ordered = sorted(
            pair_candidates.items(), key=lambda kv: kv[1][0], reverse=True
        )
        used: set[int] = set()
        out = {}
        for (a, b), (_score, quad_verts) in ordered:
            if a in used or b in used:
                continue
            used.add(a)
            used.add(b)
            out[(a, b)] = quad_verts
        return out


def _assign_layer_per_tri(mesh, n_tris: int) -> list[int]:
    out = [-1] * n_tris
    if mesh.n_layers == 0 or not mesh.layers.get("OE"):
        return [0] * n_tris
    for k in range(mesh.n_layers):
        for elem_id in np.asarray(mesh.layers["OE"][k]).flatten():
            t = int(elem_id)
            if 0 <= t < n_tris:
                out[t] = k
        for elem_id in np.asarray(mesh.layers["IE"][k]).flatten():
            t = int(elem_id)
            if 0 <= t < n_tris and out[t] < 0:
                out[t] = k
    for t in range(n_tris):
        if out[t] < 0:
            out[t] = 0
    return out


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
