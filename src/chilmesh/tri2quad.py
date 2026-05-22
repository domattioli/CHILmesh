"""Triangle-to-quad conversion — faithful port of QuADMesh's Tri2QuadRoutine.

Direct port of ``QuADMesh-MATLAB/02_QuADMESH_Library/02_Tri2Quad_Routine``:

  1. Skeletonize the input. Iterate layers inner → outer (k = nLayers-1 → 0).
  2. ``identifyEdgesFun_v2``: for the current layer's sub-domain
     (``OE[k] ∪ IE[k]``), walk each outer-vertex path
     (``paths_on_outer_vertices``). At each path vertex sort the
     incident sub-domain edges CCW. The "up-edge" is the sub-domain
     boundary edge to the previous path vertex; the "down-edge" is the
     boundary edge to the next. Rotate so up-edge is first; if
     down-edge is at position 1 (immediately after up), reverse the
     remainder. Pick interior edges between up and down using the
     every-other-edge rule (when an edge is selected, its two adjacent
     tris are marked consumed; future edges sharing a consumed tri are
     skipped — the every-other pattern emerges from this).
  3. ``mergeTrianglesFun``: each selected edge merges the two adjacent
     tris into a quad.
  4. Leftover tris in layer k pass to layer k-1's pool.
  5. After all layers, leftover tris touching the mesh boundary are
     kept as padded-tri rows ``[a, b, c, a]`` (this matches
     ``removeTrianglesFun``'s boundary-only role per user direction).

No Blossom, no greedy matching, no vertex insertion, no edge flips,
no midpoint bisection. Pure faithful port. Whatever interior tris
remain are reported as RuntimeError under ``strict=True``.
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
            as padded-tri rows.

    Returns:
        New CHILmesh whose elements are quads (4-column) with boundary
        triangles encoded in the padded form ``[a, b, c, a]``.

    Raises:
        ValueError: If the input is not a pure-triangle mesh.
        RuntimeError: If ``strict`` and an interior triangle survives.
    """
    from chilmesh import CHILmesh
    from chilmesh.layer_paths import paths_on_outer_vertices

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

    edge_to_tris_global = _build_edge_map(tris)
    boundary_vert_ids = _boundary_vertex_ids(edge_to_tris_global)

    quad_rows: list[list[int]] = []
    consumed: set[int] = set()
    deferred: list[int] = []

    for layer in range(n_layers - 1, -1, -1):
        layer_pool: set[int] = {
            t for t in range(n_tris)
            if t not in consumed and tri_layer[t] == layer
        }
        layer_pool.update(deferred)
        deferred = []
        if not layer_pool:
            continue

        sub_edges = _restrict_edge_map(edge_to_tris_global, layer_pool)

        ov_ids: set[int] = set()
        if layer < mesh.n_layers and mesh.layers.get("OV"):
            ov_ids = {int(v) for v in np.asarray(mesh.layers["OV"][layer]).flatten()}

        try:
            paths = paths_on_outer_vertices(mesh, layer) if layer < mesh.n_layers else []
        except IndexError:
            paths = []

        selected = _identify_edges_fun_v2(
            paths, sub_edges, points, ov_ids,
        )

        for (a, b) in selected:
            tlist = [t for t in sub_edges.get((a, b), []) if t not in consumed]
            if len(tlist) != 2:
                continue
            t_a, t_b = tlist
            quad_verts = _merge_tris_fun(tris[t_a], tris[t_b], a, b)
            if quad_verts is None:
                continue
            quad_rows.append(list(quad_verts))
            consumed.add(t_a)
            consumed.add(t_b)

        leftover = [t for t in layer_pool if t not in consumed]
        if layer > 0:
            deferred.extend(leftover)

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
                f"after layer-by-layer path-walk identifyEdgesFun_v2 + "
                f"mergeTrianglesFun (first 10: {leftover_interior[:10]})."
            )
        for t_id in leftover_interior:
            verts = tris[t_id].tolist()
            quad_rows.append(
                [int(verts[0]), int(verts[1]), int(verts[2]), int(verts[0])]
            )
            consumed.add(t_id)

    new_conn = np.array(quad_rows, dtype=int)
    return CHILmesh(
        connectivity=new_conn,
        points=points.copy(),
        grid_name=mesh.grid_name,
        compute_layers=True,
    )


def _identify_edges_fun_v2(
    paths,
    sub_edges: dict[tuple[int, int], list[int]],
    points: np.ndarray,
    ov_set: set[int],
) -> list[tuple[int, int]]:
    """Faithful port of MATLAB identifyEdgesFun_v2."""
    incident_by_vert: dict[int, list[tuple[int, int]]] = {}
    for (a, b) in sub_edges.keys():
        incident_by_vert.setdefault(a, []).append((a, b))
        incident_by_vert.setdefault(b, []).append((a, b))

    sub_boundary_edges_all: set[tuple[int, int]] = {
        key for key, tlist in sub_edges.items() if len(tlist) == 1
    }
    sub_boundary_edges = {
        (a, b) for (a, b) in sub_boundary_edges_all
        if a in ov_set or b in ov_set
    }

    edge_consumed: set[tuple[int, int]] = set()
    tri_consumed: set[int] = set()
    selected: list[tuple[int, int]] = []

    for path in paths:
        path = [int(v) for v in path]
        if len(path) < 2:
            continue
        if path[0] == path[-1]:
            path_seq = path[:-1]
            is_closed = True
        else:
            path_seq = path
            is_closed = False
        if not path_seq:
            continue

        # MATLAB: pathVertIDs = [end; path; start] for closed paths
        if is_closed:
            ext = [path_seq[-1]] + path_seq + [path_seq[0]]
        else:
            ext = path_seq

        down_edge = _initial_down_edge(ext, sub_boundary_edges)

        for i in range(1, len(ext) - 1):
            v_prev = ext[i - 1]
            v_cur = ext[i]
            v_next = ext[i + 1]

            up_edge = down_edge
            new_down = _make_key(v_cur, v_next)
            if new_down in sub_boundary_edges:
                down_edge = new_down
            else:
                cand = None
                for e in incident_by_vert.get(v_cur, []):
                    if e in sub_boundary_edges and v_next in e:
                        cand = e
                        break
                down_edge = cand

            if up_edge is None or down_edge is None:
                continue
            if v_cur not in ov_set:
                continue
            if v_cur not in incident_by_vert:
                continue

            ccw_all = _ccw_edges_around_vert(v_cur, incident_by_vert[v_cur], points)
            ccw_active = [e for e in ccw_all if e not in edge_consumed]
            if up_edge not in ccw_active:
                continue
            idx_up = ccw_active.index(up_edge)
            rotated = ccw_active[idx_up:] + ccw_active[:idx_up]
            if down_edge not in rotated:
                continue
            idx_down = rotated.index(down_edge)
            if idx_down == 1 and len(rotated) > 2:
                rotated = [rotated[0]] + list(reversed(rotated[1:]))
                idx_down = rotated.index(down_edge)
            if idx_down <= 1:
                continue

            for j in range(1, idx_down):
                edge = rotated[j]
                if edge in edge_consumed:
                    continue
                tlist = sub_edges.get(edge, [])
                if len(tlist) != 2:
                    continue
                t_a, t_b = tlist
                if t_a in tri_consumed or t_b in tri_consumed:
                    continue
                selected.append(edge)
                edge_consumed.add(edge)
                tri_consumed.add(t_a)
                tri_consumed.add(t_b)
    return selected


def _initial_down_edge(
    ext: list[int],
    sub_boundary_edges: set[tuple[int, int]],
) -> tuple[int, int] | None:
    if len(ext) < 2:
        return None
    k = _make_key(ext[0], ext[1])
    if k in sub_boundary_edges:
        return k
    for edge in sub_boundary_edges:
        if ext[0] in edge:
            return edge
    return None


def _ccw_edges_around_vert(
    vert_id: int,
    incident_edges: list[tuple[int, int]],
    points: np.ndarray,
) -> list[tuple[int, int]]:
    p0 = points[vert_id, :2]
    angled: list[tuple[float, tuple[int, int]]] = []
    for edge in incident_edges:
        other = edge[1] if edge[0] == vert_id else edge[0]
        d = points[other, :2] - p0
        angled.append((float(np.arctan2(d[1], d[0])), edge))
    angled.sort()
    return [e for _ang, e in angled]


def _make_key(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


def _restrict_edge_map(
    edge_to_tris: dict[tuple[int, int], list[int]],
    pool: set[int],
) -> dict[tuple[int, int], list[int]]:
    out: dict[tuple[int, int], list[int]] = {}
    for key, tlist in edge_to_tris.items():
        usable = [t for t in tlist if t in pool]
        if usable:
            out[key] = usable
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
            key = _make_key(a, b)
            edge_to_tris.setdefault(key, []).append(t_id)
    return edge_to_tris


def _boundary_vertex_ids(
    edge_to_tris: dict[tuple[int, int], list[int]],
) -> set[int]:
    out: set[int] = set()
    for key, tlist in edge_to_tris.items():
        if len(tlist) == 1:
            out.add(key[0])
            out.add(key[1])
    return out


def _merge_tris_fun(
    tri_a: np.ndarray,
    tri_b: np.ndarray,
    shared_v0: int,
    shared_v1: int,
) -> tuple[int, int, int, int] | None:
    """Faithful port of MATLAB mergeTrianglesFun: merge 2 tris into 1 quad."""
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


def _count_boundary_edges(
    tri_verts: list[int],
    boundary_vert_ids: set[int],
    edge_to_tris: dict[tuple[int, int], list[int]],
) -> tuple[int, list[int]]:
    """Count boundary edges of a triangle.

    Returns: (count, indices of boundary edges in [edge01, edge12, edge20])
    """
    v0, v1, v2 = tri_verts
    edges = [(v0, v1), (v1, v2), (v2, v0)]
    count = 0
    boundary_edge_indices = []
    for i, (a, b) in enumerate(edges):
        key = _make_key(a, b)
        if len(edge_to_tris.get(key, [])) == 1:
            count += 1
            boundary_edge_indices.append(i)
    return count, boundary_edge_indices


def _flip_shared_edge(
    tri_a: np.ndarray,
    tri_b: np.ndarray,
    shared_v0: int,
    shared_v1: int,
) -> tuple[tuple[int, int, int], tuple[int, int, int]] | None:
    """Flip shared edge between two adjacent triangles (no vertex insertion).

    Takes two tris sharing edge (shared_v0, shared_v1) and returns them with the
    diagonal flipped. This is a topological operation only; coordinates unchanged.

    For tris [a,b,c] and [a,b,d] sharing edge (a,b), returns:
      ([a,d,c], [b,c,d])
    where the new shared edge is (c,d).

    Per user direction: no midpoint insertion, no vertex addition.
    Result tris are deferred to next layer for processing.
    """
    tri_a_l = [int(v) for v in tri_a]
    tri_b_l = [int(v) for v in tri_b]

    unique_a = [v for v in tri_a_l if v != shared_v0 and v != shared_v1]
    unique_b = [v for v in tri_b_l if v != shared_v0 and v != shared_v1]

    if len(unique_a) != 1 or len(unique_b) != 1:
        return None

    ua = unique_a[0]
    ub = unique_b[0]

    return ((shared_v0, ub, ua), (shared_v1, ua, ub))
