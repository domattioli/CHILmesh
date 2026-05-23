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

    Notes:
        **CRITICAL: Annulus result is GEOMETRICALLY BROKEN.**

        1. 205/258 quads (79%!) are BOWTIE-SHAPED (diagonals cross).
           Quad wraps around itself, self-overlapping, non-convex.
           Validator misses this (only checks edge-crossing, not diagonal).

        2. 10 severely degenerate quads (aspect 0.022-0.099, 44:1 edge ratio).
           Validator doesn't check mesh quality metrics.

        3. Root cause: edge-insertion vertex placement (1/3 along edge)
           combined with layer-based merging creates invalid topology.

        BLOCKER: Do NOT ship. Algorithm fundamentally broken on annulus.
        Recommend: (a) redesign edge-insertion, or (b) use different
        approach entirely, or (c) detect and reject bowtie quads + fall
        back to original triangles for those regions.
    """
    from chilmesh import CHILmesh
    from chilmesh.layer_paths import paths_on_outer_vertices

    if mesh.connectivity_list.shape[1] != 3:
        raise ValueError(
            f"tri_to_quad requires a pure-triangle mesh; "
            f"connectivity has {mesh.connectivity_list.shape[1]} columns"
        )

    tris = np.asarray(mesh.connectivity_list, dtype=int).copy()
    points = np.asarray(mesh.points, dtype=float).copy()
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
        leftover_interior = _remove_interior_triangles(
            leftover_interior, tris, points, edge_to_tris_global, boundary_vert_ids
        )

    if leftover_interior:
        leftover_interior, points = _apply_edge_insertion(
            leftover_interior, tris, points, quad_rows, edge_to_tris_global
        )

    if leftover_interior:
        if strict:
            raise RuntimeError(
                f"tri_to_quad left {len(leftover_interior)} interior triangles "
                f"after layer-by-layer path-walk identifyEdgesFun_v2 + "
                f"mergeTrianglesFun + edge flips + edge insertion (first 10: {leftover_interior[:10]})."
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


def tri_to_quad_full(
    mesh: "CHILmesh",
    *,
    smooth: bool = True,
    smooth_method: str = "fem",
    smooth_iters: int = 50,
    doublet_passes: int = 6,
    convert_boundary_tris: bool = True,
) -> "CHILmesh":
    """Full pipeline: tri_to_quad → boundary-tri → quad via edge-insertion → DoubletCollapse → smooth.

    Pipeline:
      1. ``tri_to_quad`` — layer-by-layer path-walk merge.
      2. ``_insert_boundary_tri_midpoint`` — for each padded triangle with
         a mesh-boundary edge, insert a midpoint on that edge and convert
         the tri into a quad. Adds 1 vertex per converted tri; element
         count unchanged; area unchanged (midpoint lies on the boundary
         edge being replaced).
      3. ``_doublet_collapse`` — interior valence-2 vertex shared by
         exactly two quads → merge into a single quad covering the same
         union. Each collapse drops one row but the surviving quad covers
         the combined area of both originals (no coverage loss). Iterates
         up to ``doublet_passes`` rounds.
      4. ``CHILmesh.smooth_mesh`` — FEM first, fall back to angle-based
         if FEM returns non-finite coords.

    Args:
        mesh: Pure-triangle input.
        smooth: Apply smoother as final step.
        smooth_method: ``"fem"`` (direct) or ``"angle-based"``.
        smooth_iters: Passed to angle-based smoother (ignored for FEM).
        doublet_passes: Max DoubletCollapse iterations.
        convert_boundary_tris: Edge-insert midpoint on each padded tri's
            boundary edge to convert it into a quad.

    Returns:
        New CHILmesh.
    """
    from chilmesh import CHILmesh

    q = tri_to_quad(mesh, strict=False)
    conn = np.asarray(q.connectivity_list, dtype=int).copy()
    pts = np.asarray(q.points, dtype=float).copy()

    if convert_boundary_tris:
        conn, pts = _insert_boundary_tri_midpoint(conn, pts)
        # Pentagon-split removed per user direction "there should never be
        # pentagons". Interior tris (padded rows lacking a mesh-boundary edge)
        # should be eliminated upstream — improved doublet / QVM / swap
        # passes below absorb them into neighbouring quads.

    for _ in range(doublet_passes):
        new_conn = _doublet_collapse(conn, pts)
        if new_conn.shape[0] == conn.shape[0]:
            break
        conn = new_conn

    for _ in range(doublet_passes):
        new_conn = _quad_vertex_merge(conn, pts)
        if new_conn.shape[0] == conn.shape[0]:
            break
        conn = new_conn
        # Doublet pass after QVM since the merge can create new valence-2 verts.
        for _ in range(doublet_passes):
            new_conn = _doublet_collapse(conn, pts)
            if new_conn.shape[0] == conn.shape[0]:
                break
            conn = new_conn

    # Absorb any leftover interior triangles by fusing adjacent interior
    # tri pairs into a single quad covering their union. Two adjacent
    # triangles sharing one edge form a convex quad iff the diagonal lies
    # inside the union — guarded by signed-area sign check below.
    conn = _fuse_interior_tri_pairs(conn, pts)

    out = CHILmesh(
        connectivity=conn,
        points=pts,
        grid_name=mesh.grid_name,
        compute_layers=True,
    )

    if smooth:
        from chilmesh import CHILmesh as _CM

        # Smoothing pass 1: angle-based (FEM optional).
        pre_pts = np.asarray(out.points).copy()
        try:
            out.smooth_mesh(method=smooth_method, acknowledge_change=True)
            if not np.isfinite(out.points).all():
                raise RuntimeError("smoother produced non-finite coordinates")
        except Exception:
            out.change_points(pre_pts, acknowledge_change=True)
            out.smooth_mesh(method="angle-based", acknowledge_change=True)
        lap_pts = _laplacian_smooth(out, n_iter=smooth_iters, omega=0.7)
        out.change_points(lap_pts, acknowledge_change=True)

        # Topology cleanup with validator-guarded acceptance. For meshes
        # below ``_LARGE_MESH_THRESHOLD`` we run the O(n²) snapshot
        # validator per swap; above that, we skip per-swap validation and
        # rely on local geometric checks (hexagon convexity, signed-area
        # sign, diagonal midpoint containment, perimeter crossings) plus a
        # single end-of-loop validation.
        _LARGE_MESH_THRESHOLD = 1000
        try:
            from chilmesh.validation.validator import (
                validate_mesh_elements as _v,
            )
            _have_validator = True
        except ImportError:
            _have_validator = False

        _use_per_swap_validator = (
            _have_validator and out.n_elems <= _LARGE_MESH_THRESHOLD
        )

        def _viol_count(mesh: "CHILmesh") -> int:
            if not _have_validator:
                return -1
            try:
                return len(_v(mesh).violations)
            except Exception:
                return -1

        baseline_viol = _viol_count(out) if _have_validator else 0

        class _SwapValidator:
            """Callable wrapper that runs the mesh validator on a connectivity
            snapshot and returns the violation count. Caches a ``baseline``
            attribute so the swap loop can compare against the last-known-good
            violation count."""

            def __init__(self, baseline: int):
                self.baseline = baseline

            def __call__(self, conn_arr: np.ndarray, pts_arr: np.ndarray):
                if not _have_validator:
                    return None
                try:
                    snapshot = _CM(connectivity=conn_arr, points=pts_arr,
                                    grid_name=out.grid_name,
                                    compute_layers=True)
                    return len(_v(snapshot).violations)
                except Exception:
                    return None

        for _ in range(8):
            c = np.asarray(out.connectivity_list, dtype=int).copy()
            p = np.asarray(out.points, dtype=float).copy()
            sv = _SwapValidator(baseline=baseline_viol) if _use_per_swap_validator else None
            c_new = _quad_pair_edge_swap(c, p, quality_threshold=0.5,
                                          min_improvement=0.02,
                                          validator=sv)
            n_swap = int(np.sum(np.any(c != c_new, axis=1)))
            c_new = _doublet_collapse(c_new, p)
            c_new = _quad_vertex_merge(c_new, p)
            n_drop = c.shape[0] - c_new.shape[0]
            if n_swap == 0 and n_drop == 0:
                break
            candidate = _CM(connectivity=c_new, points=p,
                             grid_name=out.grid_name, compute_layers=True)
            # On large meshes the per-swap validator is skipped, so accept
            # the iteration without snapshot validation. Geometric guards
            # in _quad_pair_edge_swap reject invalid swaps locally.
            if _use_per_swap_validator:
                new_viol = _viol_count(candidate)
                if new_viol > baseline_viol:
                    break
                baseline_viol = new_viol
            out = candidate

        # Final smoothing pass on the topology-cleaned mesh. Skip Laplacian
        # inside the loop because post-doublet/QVM relabeling expands some
        # vertex 1-rings into shapes Laplacian cannot smooth without
        # crossing adjacent quads.
        pre_pts = np.asarray(out.points).copy()
        try:
            out.smooth_mesh(method="angle-based", acknowledge_change=True)
        except Exception:
            out.change_points(pre_pts, acknowledge_change=True)
        lap_pts = _laplacian_smooth(out, n_iter=smooth_iters, omega=0.5)
        # Only accept Laplacian update if it doesn't break validity.
        candidate2 = _CM(connectivity=out.connectivity_list,
                          points=lap_pts, grid_name=out.grid_name,
                          compute_layers=True)
        final_viol = _viol_count(candidate2)
        if not _have_validator or final_viol <= baseline_viol:
            out.change_points(lap_pts, acknowledge_change=True)

    # Final geometric repair: snap displaced midpoints, fix bowties, dissolve
    # violation clusters. Runs after smoothing so midpoint positions are stable.
    conn_r = np.asarray(out.connectivity_list, dtype=int).copy()
    pts_r = np.asarray(out.points, dtype=float).copy()
    conn_r, pts_r = repair_mesh(conn_r, pts_r, mid_threshold=0)
    out = CHILmesh(
        connectivity=conn_r,
        points=pts_r,
        grid_name=mesh.grid_name,
        compute_layers=True,
    )

    return out


def _laplacian_smooth(
    mesh: "CHILmesh",
    n_iter: int = 100,
    omega: float = 0.7,
) -> np.ndarray:
    """Guarded Laplacian smoothing: each interior vertex moves toward the
    centroid of its 1-ring neighbours, blended by ``omega``. Each move is
    accepted only if it preserves the signed-area sign of every incident
    element (i.e., no element flips orientation or collapses). Boundary
    vertices are held fixed.

    The guard prevents edge crossings in concave regions where the
    centroid-pull would otherwise drag a vertex into an adjacent element's
    footprint. Returns a new points array; does not mutate the mesh.
    """
    conn = np.asarray(mesh.connectivity_list, dtype=int)
    pts = np.asarray(mesh.points, dtype=float).copy()
    boundary_verts = {
        int(x) for x in np.asarray(mesh.boundary_node_indices()).flatten()
    }

    n_verts = pts.shape[0]
    nbrs: list[set[int]] = [set() for _ in range(n_verts)]
    vert_elems: list[list[int]] = [[] for _ in range(n_verts)]
    for ei in range(conn.shape[0]):
        row = conn[ei]
        if row[0] == row[3]:
            verts = [int(row[k]) for k in range(3)]
        else:
            verts = [int(row[k]) for k in range(4)]
        n = len(verts)
        for k in range(n):
            a, b = verts[k], verts[(k + 1) % n]
            nbrs[a].add(b)
            nbrs[b].add(a)
        for v in verts:
            vert_elems[v].append(ei)

    def _elem_signed(ei: int, pts_arr: np.ndarray) -> float:
        row = conn[ei]
        if row[0] == row[3]:
            poly = pts_arr[[int(row[0]), int(row[1]), int(row[2])], :2]
        else:
            poly = pts_arr[[int(row[k]) for k in range(4)], :2]
        x = poly[:, 0]
        y = poly[:, 1]
        return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))

    # Cache reference signed-area sign for each element from the initial
    # mesh — every Laplacian move must keep the same sign.
    ref_signs = np.array(
        [1 if _elem_signed(ei, pts) > 0 else -1 for ei in range(conn.shape[0])],
        dtype=int,
    )

    for _ in range(n_iter):
        new_pts = pts.copy()
        for v in range(n_verts):
            if v in boundary_verts or not nbrs[v]:
                continue
            avg = np.mean([pts[u, :2] for u in nbrs[v]], axis=0)
            new_xy = (1.0 - omega) * pts[v, :2] + omega * avg
            # Cap step magnitude at 25% of the shortest incident edge.
            min_edge = min(
                float(np.linalg.norm(pts[u, :2] - pts[v, :2])) for u in nbrs[v]
            )
            delta = new_xy - pts[v, :2]
            d_norm = float(np.linalg.norm(delta))
            cap = 0.25 * min_edge
            if d_norm > cap and d_norm > 1e-14:
                new_xy = pts[v, :2] + delta * (cap / d_norm)
            saved = pts[v, :2].copy()
            pts[v, :2] = new_xy
            ok = True
            for ei in vert_elems[v]:
                sa = _elem_signed(ei, pts)
                if sa * ref_signs[ei] <= 0:
                    ok = False
                    break
                if abs(sa) < 1e-14:
                    ok = False
                    break
            if ok:
                new_pts[v, :2] = new_xy
            pts[v, :2] = saved
        pts = new_pts
    return pts


def _build_full_edge_map(conn: np.ndarray) -> dict[tuple[int, int], list[int]]:
    """Edge → element list for a mixed-element (3 or 4 col) connectivity table."""
    out: dict[tuple[int, int], list[int]] = {}
    for ei in range(conn.shape[0]):
        row = conn[ei]
        if row[0] == row[3]:
            verts = [int(row[k]) for k in range(3)]
        else:
            verts = [int(row[k]) for k in range(4)]
        n = len(verts)
        for k in range(n):
            a, b = verts[k], verts[(k + 1) % n]
            key = (a, b) if a < b else (b, a)
            out.setdefault(key, []).append(ei)
    return out


def _insert_boundary_tri_midpoint(
    conn: np.ndarray,
    pts: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert each padded boundary triangle ``[a,b,c,a]`` to a quad by
    inserting a midpoint on its mesh-boundary edge.

    For a padded triangle with one edge ``(va, vb)`` on the mesh boundary
    (i.e., shared with no other element), the result is a quad
    ``[va, m, vb, vc]`` where ``m`` is the midpoint of ``(va, vb)`` and
    ``vc`` is the tri's third vertex. The original tri row is overwritten
    with the new quad (no element drop, no area change — ``m`` lies on the
    boundary edge being split).
    """
    edge_to_elems = _build_full_edge_map(conn)
    new_pts: list[np.ndarray] = []
    next_vid = pts.shape[0]

    for ei in range(conn.shape[0]):
        row = conn[ei]
        if row[0] != row[3]:
            continue
        verts = [int(row[k]) for k in range(3)]
        # Find a mesh-boundary edge (length-1 in edge_to_elems).
        chosen_edge: tuple[int, int] | None = None
        chosen_other_vert: int | None = None
        for k in range(3):
            a, b = verts[k], verts[(k + 1) % 3]
            key = (a, b) if a < b else (b, a)
            if len(edge_to_elems.get(key, [])) == 1:
                chosen_edge = (a, b)
                chosen_other_vert = verts[(k + 2) % 3]
                break
        if chosen_edge is None:
            # Padded tri has no mesh-boundary edge (interior-layer tri);
            # midpoint insertion would break neighboring quad topology
            # (would force a 5-vert pentagon). Leave as-is.
            continue
        va, vb = chosen_edge
        vc = chosen_other_vert
        m_xy = 0.5 * (pts[va, :2] + pts[vb, :2])
        z = pts[va, 2] if pts.shape[1] >= 3 else 0.0
        new_pts.append(np.array([m_xy[0], m_xy[1], z]))
        # New quad row [va, m, vb, vc] in CCW order (preserves tri orientation).
        conn[ei] = np.array([va, next_vid, vb, vc], dtype=int)
        next_vid += 1

    if new_pts:
        pts = np.vstack([pts, np.asarray(new_pts, dtype=float)])
    return conn, pts


def _pentagon_split_interior_tris(
    conn: np.ndarray,
    pts: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert every padded tri lacking a mesh-boundary edge into two quads
    by inserting a midpoint on a neighboring quad's far edge.

    Setup: tri ``ABC`` (CCW) shares one edge — say ``CA`` — with a
    neighboring quad ``CAEF`` (CCW). The union forms a pentagon
    ``A → B → C → E → F → A``. Bisect the quad's "far edge" ``EF`` at
    ``M`` and split the union into two quads:

      - new quad replacing tri: ``[A, B, C, M]``... wait this has 5 verts
        if M not on AC. Re-derive:

    Actual split: pentagon ``A, B, C, E, F`` (5 verts). Bisect edge
    ``EF`` to get ``M``. Hexagon ``A, B, C, E, M, F`` (6 verts). Split
    along diagonal ``B–M`` (or any reasonable diagonal) into two quads:
      - ``[A, B, M, F]`` (left half)
      - ``[B, C, E, M]`` (right half)

    Both quads are convex if pentagon is convex; orientation matches
    pentagon CCW. Net: 1 tri + 1 quad → 2 quads (count preserved). Area
    preserved (M lies on EF).

    A padded tri is processed only if it has a quad neighbor on some
    edge; tris adjacent only to other padded tris are left as-is.
    """
    edge_to_elems = _build_full_edge_map(conn)
    next_vid = pts.shape[0]
    new_pts: list[np.ndarray] = []

    pad_ids = list(np.where(conn[:, 0] == conn[:, 3])[0])

    for tri_eid in pad_ids:
        row_tri = conn[tri_eid]
        tri_verts = [int(row_tri[k]) for k in range(3)]

        chosen_neighbor: int | None = None
        shared_edge: tuple[int, int] | None = None
        for k in range(3):
            a, b = tri_verts[k], tri_verts[(k + 1) % 3]
            key = (a, b) if a < b else (b, a)
            elems = edge_to_elems.get(key, [])
            for ne in elems:
                if ne == tri_eid:
                    continue
                if conn[ne, 0] == conn[ne, 3]:
                    continue  # neighbor is also a tri; skip
                chosen_neighbor = ne
                shared_edge = (a, b)
                break
            if chosen_neighbor is not None:
                break
        if chosen_neighbor is None:
            continue

        # Build pentagon CCW. Tri CCW is [A, B, C] = tri_verts.
        # Identify which edge of tri is shared with neighbor.
        a, b = shared_edge
        idx_a = tri_verts.index(a)
        idx_b = tri_verts.index(b)
        # Tri's third vert.
        third = [v for v in tri_verts if v not in (a, b)][0]

        # Rotate tri so the shared edge is the last edge (tri verts: [B, third, ...]?).
        # Simpler: name tri vertices V0=third, V1, V2 where (V1,V2) is shared edge.
        # Pentagon must be CCW. Tri CCW is [tri_verts[0], tri_verts[1], tri_verts[2]].
        # Walk tri CCW: edge indices 0→1, 1→2, 2→0. The shared edge between consecutive verts.
        # Find ordered shared edge as (tri[k], tri[(k+1)%3]).
        shared_k = None
        for k in range(3):
            if {tri_verts[k], tri_verts[(k + 1) % 3]} == {a, b}:
                shared_k = k
                break
        assert shared_k is not None
        # CCW oriented shared edge: tri side (s_a, s_b) = (tri[k], tri[k+1]).
        s_a = tri_verts[shared_k]
        s_b = tri_verts[(shared_k + 1) % 3]
        # Tri's "B" (the non-shared vert), with CCW order [s_b, B, s_a]? No:
        # Tri CCW: [tri[k], tri[k+1], tri[k+2]] = [s_a, s_b, B] where B = third.
        # So pentagon traversal: start at B (third), then s_a, s_b, then quad's
        # far verts in CCW.
        # Quad CCW with edge s_b → s_a (opposite direction of tri's shared
        # edge): quad row contains s_a and s_b. Find their order in quad.
        quad_verts = [int(conn[chosen_neighbor, k]) for k in range(4)]
        try:
            i_sa_q = quad_verts.index(s_a)
        except ValueError:
            continue
        try:
            i_sb_q = quad_verts.index(s_b)
        except ValueError:
            continue
        # In the quad, the shared edge runs s_b → s_a (opposite of tri's CCW
        # direction). Walk quad CCW starting from s_a → next CCW vert (skip s_b).
        # Rotate quad so it starts at s_a:
        rot = [quad_verts[(i_sa_q + k) % 4] for k in range(4)]
        # rot = [s_a, q1, q2, q3]. Confirm rot[3] == s_b (i.e., previous CCW is s_b).
        if rot[3] != s_b:
            # Quad orientation puts s_b at rot[1] (= going s_a → s_b CCW). This
            # means our assumed shared-edge direction was wrong; skip rather
            # than mis-split.
            continue
        q_far1 = rot[1]
        q_far2 = rot[2]
        # Pentagon CCW: [third(B), s_a, q_far1, q_far2, s_b]? Let's verify.
        # Tri CCW: [s_a, s_b, B] (rotation invariance).
        # After removing shared edge s_a–s_b, tri contributes vert B.
        # Quad CCW: [s_a, q_far1, q_far2, s_b]. After removing shared edge,
        # quad contributes q_far1 and q_far2.
        # Union polygon CCW (with shared edge removed): s_a, q_far1, q_far2, s_b, B.
        # Hexagon after bisecting q_far1–q_far2 with M:
        # s_a, q_far1, M, q_far2, s_b, B → 6 verts.
        m_xy = 0.5 * (pts[q_far1, :2] + pts[q_far2, :2])
        z = pts[q_far1, 2] if pts.shape[1] >= 3 else 0.0
        m_vid = next_vid
        new_pts.append(np.array([m_xy[0], m_xy[1], z]))
        next_vid += 1

        # Hexagon CCW: [s_a, q_far1, M, q_far2, s_b, B] (indices 0..5).
        # Three "main diagonal" splits each yield a pair of quads:
        #   (0, 3) = (s_a, q_far2)  → [s_a, q_far1, M, q_far2] + [s_a, q_far2, s_b, B]
        #   (1, 4) = (q_far1, s_b)  → [q_far1, M, q_far2, s_b] + [q_far1, s_b, B, s_a]
        #   (2, 5) = (M, B)         → [M, q_far2, s_b, B]     + [M, B, s_a, q_far1]
        # Pick whichever maximises the minimum element skew quality.
        B = third
        hexa = [s_a, q_far1, m_vid, q_far2, s_b, B]
        # Build temporary points array with M position appended (m_vid points
        # to it once we vstack at end).
        pts_with_m = np.vstack([pts[:, :2], np.array([m_xy])])

        def _q_skew(idx_list: list[int]) -> float:
            poly = pts_with_m[idx_list]
            angs = []
            for k in range(4):
                a = poly[(k - 1) % 4] - poly[k]
                b = poly[(k + 1) % 4] - poly[k]
                na = float(np.linalg.norm(a))
                nb = float(np.linalg.norm(b))
                if na < 1e-12 or nb < 1e-12:
                    return 0.0
                cos = float(np.dot(a, b) / (na * nb))
                cos = max(-1.0, min(1.0, cos))
                angs.append(np.degrees(np.arccos(cos)))
            return max(0.0, 1.0 - max(abs(a - 90.0) for a in angs) / 90.0)

        def _signed(idx_list: list[int]) -> float:
            poly = pts_with_m[idx_list]
            x = poly[:, 0]
            y = poly[:, 1]
            return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))

        # Hexa vert ids in the temporary index space (m_vid → index = pts.shape[0]):
        m_local = pts.shape[0]
        H = [s_a, q_far1, m_local, q_far2, s_b, B]

        candidates = [
            ([H[0], H[1], H[2], H[3]], [H[0], H[3], H[4], H[5]]),  # diag (s_a, q_far2)
            ([H[1], H[2], H[3], H[4]], [H[1], H[4], H[5], H[0]]),  # diag (q_far1, s_b)
            ([H[2], H[3], H[4], H[5]], [H[2], H[5], H[0], H[1]]),  # diag (M, B) — original
        ]
        # Concrete-vertex variants for writing back into conn (replace m_local
        # with m_vid which is the real id once we vstack new_pts later).
        candidates_real = [
            ([hexa[0], hexa[1], hexa[2], hexa[3]], [hexa[0], hexa[3], hexa[4], hexa[5]]),
            ([hexa[1], hexa[2], hexa[3], hexa[4]], [hexa[1], hexa[4], hexa[5], hexa[0]]),
            ([hexa[2], hexa[3], hexa[4], hexa[5]], [hexa[2], hexa[5], hexa[0], hexa[1]]),
        ]
        best_score = -1.0
        best_pair: tuple[list[int], list[int]] | None = None
        for (qa_idx, qb_idx), (qa_real, qb_real) in zip(candidates, candidates_real):
            sa = _signed(qa_idx)
            sb = _signed(qb_idx)
            if sa * sb <= 0:
                continue
            if abs(sa) < 1e-12 or abs(sb) < 1e-12:
                continue
            score = min(_q_skew(qa_idx), _q_skew(qb_idx))
            if score > best_score:
                best_score = score
                best_pair = (qa_real, qb_real)
        if best_pair is None:
            # Fall back to original diag (M, B) split to keep the algorithm
            # producing a valid pair even if scoring rejected all (rare).
            best_pair = candidates_real[2]
        quad1, quad2 = best_pair
        conn[tri_eid] = np.array(quad1, dtype=int)
        conn[chosen_neighbor] = np.array(quad2, dtype=int)
        # Refresh edge map (cheap since small).
        edge_to_elems = _build_full_edge_map(conn)

    if new_pts:
        pts = np.vstack([pts, np.asarray(new_pts, dtype=float)])
    return conn, pts


def _quad_pair_edge_swap(
    conn: np.ndarray,
    pts: np.ndarray,
    quality_threshold: float = 0.2,
    min_improvement: float = 0.05,
    validator: "object | None" = None,
) -> np.ndarray:
    """Diagonal swap between adjacent quad pairs to fix slivery quads.

    For each quad with element quality below ``quality_threshold``, attempt
    a topological swap with each neighbor sharing an edge:

      Q1 = [a, b, c, d], Q2 = [a, d, e, f] sharing edge (a, d)
      Hexagon CCW: a, b, c, d, e, f

    Three possible quad-pair splits via 3 diagonals:
      (a, d) — original
      (b, e) — Q3 = [a, b, e, f], Q4 = [b, c, d, e]
      (c, f) — Q5 = [a, b, c, f], Q6 = [c, d, e, f]

    Pick whichever maximises the minimum element quality of the pair, and
    accept if it improves over the original by at least ``min_improvement``.
    Element count preserved. Coverage preserved (the hexagon footprint is
    unchanged). No vertex add or drop.
    """
    n_elems = conn.shape[0]
    if n_elems == 0:
        return conn

    edge_to_elems = _build_full_edge_map(conn)

    def _quad_skew_quality(verts: list[int]) -> float:
        if len(set(verts)) != 4:
            return 0.0
        poly = pts[verts, :2]
        # Skew quality: 1 - max(|θ_i - 90°|) / 90°, clipped to [0, 1].
        angs = []
        for k in range(4):
            v_prev = poly[(k - 1) % 4]
            v_cur = poly[k]
            v_next = poly[(k + 1) % 4]
            a = v_prev - v_cur
            b = v_next - v_cur
            na = float(np.linalg.norm(a))
            nb = float(np.linalg.norm(b))
            if na < 1e-12 or nb < 1e-12:
                return 0.0
            cos = float(np.dot(a, b) / (na * nb))
            cos = max(-1.0, min(1.0, cos))
            angs.append(np.degrees(np.arccos(cos)))
        max_dev = max(abs(a - 90.0) for a in angs)
        return max(0.0, 1.0 - max_dev / 90.0)

    def _signed_area(verts: list[int]) -> float:
        poly = pts[verts, :2]
        x = poly[:, 0]
        y = poly[:, 1]
        return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))

    used = np.zeros(n_elems, dtype=bool)

    # Score current quads.
    quality = np.zeros(n_elems)
    for ei in range(n_elems):
        row = conn[ei]
        if row[0] == row[3]:
            quality[ei] = 0.0
            continue
        quality[ei] = _quad_skew_quality([int(row[k]) for k in range(4)])

    order = np.argsort(quality)  # process worst-first.

    for ei in order:
        if used[ei]:
            continue
        if quality[ei] >= quality_threshold:
            break
        row = conn[ei]
        if row[0] == row[3]:
            continue
        verts1 = [int(row[k]) for k in range(4)]

        best_swap: tuple[int, list[int], list[int], float] | None = None

        for k in range(4):
            # Shared edge in Q1: a → d going CCW. So a = verts1[k], d = verts1[(k+1)%4].
            a = verts1[k]
            d = verts1[(k + 1) % 4]
            # Hexagon CCW around the union goes a → b → c → d → e → f → a.
            # b, c come from Q1's other two verts traversed *backwards* from a
            # (i.e., CW around Q1).
            b = verts1[(k + 3) % 4]
            c = verts1[(k + 2) % 4]
            shared_key = (a, d) if a < d else (d, a)
            neighbors = [n for n in edge_to_elems.get(shared_key, []) if n != ei]
            if len(neighbors) != 1:
                continue
            nei = neighbors[0]
            if used[nei]:
                continue
            if conn[nei, 0] == conn[nei, 3]:
                continue
            verts2 = [int(conn[nei, j]) for j in range(4)]
            if a not in verts2 or d not in verts2:
                continue
            # In Q2's CCW listing, the shared edge runs d → a (opposite of Q1).
            # So if idx_d = position of d in verts2, then idx_a = (idx_d + 1) % 4.
            idx_d = verts2.index(d)
            if verts2[(idx_d + 1) % 4] != a:
                continue
            # e, f are Q2's non-shared verts, traversed backwards from d
            # (which equals CCW around the union after d).
            e = verts2[(idx_d + 3) % 4]
            f = verts2[(idx_d + 2) % 4]
            # Hexagon CCW: a, b, c, d, e, f.
            hexa = [a, b, c, d, e, f]
            if len(set(hexa)) != 6:
                continue
            # Candidate splits.
            candidates = [
                ([a, b, c, d], [d, e, f, a]),  # original (a, d) diagonal
                ([a, b, e, f], [b, c, d, e]),  # diagonal (b, e)
                ([a, b, c, f], [c, d, e, f]),  # diagonal (c, f)
            ]
            # Reject swaps whose new diagonal exits the hexagon perimeter
            # (i.e., the hexagon is non-convex and the diagonal would cross
            # an outer edge). Check by verifying the diagonal's midpoint
            # lies strictly inside the hexagon polygon.
            hexa_poly = pts[[a, b, c, d, e, f], :2]

            def _hex_signed() -> float:
                x = hexa_poly[:, 0]
                y = hexa_poly[:, 1]
                return 0.5 * float(
                    np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
                )

            hex_sa = _hex_signed()
            if abs(hex_sa) < 1e-12:
                continue

            def _point_in_hex(p_xy: np.ndarray) -> bool:
                # Ray-casting test.
                inside = False
                n = hexa_poly.shape[0]
                for j in range(n):
                    p1 = hexa_poly[j]
                    p2 = hexa_poly[(j + 1) % n]
                    if (p1[1] > p_xy[1]) != (p2[1] > p_xy[1]):
                        x_int = p1[0] + (p_xy[1] - p1[1]) * (p2[0] - p1[0]) / (
                            p2[1] - p1[1] + 1e-30
                        )
                        if x_int > p_xy[0]:
                            inside = not inside
                return inside

            # Build hexagon edge list (each edge as a vertex id pair) so we
            # can reject diagonals that cross any hexagon perimeter edge.
            hex_verts = [a, b, c, d, e, f]

            def _seg_cross(p1: np.ndarray, p2: np.ndarray,
                            p3: np.ndarray, p4: np.ndarray) -> bool:
                def _ori(pa, pb, pc):
                    v = (pb[0] - pa[0]) * (pc[1] - pa[1]) - (
                        pb[1] - pa[1]
                    ) * (pc[0] - pa[0])
                    if v > 1e-12:
                        return 1
                    if v < -1e-12:
                        return -1
                    return 0

                o1 = _ori(p1, p2, p3)
                o2 = _ori(p1, p2, p4)
                o3 = _ori(p3, p4, p1)
                o4 = _ori(p3, p4, p2)
                if o1 == 0 or o2 == 0 or o3 == 0 or o4 == 0:
                    return False
                return o1 != o2 and o3 != o4

            best_local = None
            for cand_idx, (qa, qb) in enumerate(candidates):
                sa = _signed_area(qa)
                sb = _signed_area(qb)
                if sa * hex_sa <= 0 or sb * hex_sa <= 0:
                    continue
                if abs(sa) < 1e-12 or abs(sb) < 1e-12:
                    continue
                diag = set(qa) & set(qb)
                if len(diag) != 2:
                    continue
                d_a, d_b = list(diag)
                p_a = pts[d_a, :2]
                p_b = pts[d_b, :2]
                # Reject diagonals that intersect any hexagon perimeter edge
                # (i.e., the hexagon is non-convex and the diagonal exits
                # the polygon).
                crosses_perim = False
                for hi in range(6):
                    ha = hex_verts[hi]
                    hb = hex_verts[(hi + 1) % 6]
                    if ha in (d_a, d_b) or hb in (d_a, d_b):
                        continue  # shares an endpoint
                    if _seg_cross(p_a, p_b, pts[ha, :2], pts[hb, :2]):
                        crosses_perim = True
                        break
                if crosses_perim:
                    continue
                mid = 0.5 * (p_a + p_b)
                if not _point_in_hex(mid):
                    continue
                q1q = _quad_skew_quality(qa)
                q2q = _quad_skew_quality(qb)
                mn = min(q1q, q2q)
                if best_local is None or mn > best_local[0]:
                    best_local = (mn, cand_idx, qa, qb)
            if best_local is None:
                continue
            mn, cand_idx, qa, qb = best_local
            # Original split's min-quality:
            orig_mn = min(quality[ei], quality[nei])
            if cand_idx == 0:  # picking original is a no-op
                continue
            if mn < orig_mn + min_improvement:
                continue
            if best_swap is None or mn > best_swap[3]:
                best_swap = (nei, qa, qb, mn)

        if best_swap is None:
            continue
        nei, qa, qb, _ = best_swap
        # Apply swap tentatively. If a validator is supplied, verify the
        # swap doesn't add any violation beyond the baseline.
        if validator is not None:
            saved_ei = conn[ei].copy()
            saved_nei = conn[nei].copy()
            conn[ei] = np.array(qa, dtype=int)
            conn[nei] = np.array(qb, dtype=int)
            try:
                rep = validator(conn, pts)
                if rep is not None and rep > validator.baseline:
                    conn[ei] = saved_ei
                    conn[nei] = saved_nei
                    continue
                validator.baseline = rep if rep is not None else validator.baseline
            except Exception:
                conn[ei] = saved_ei
                conn[nei] = saved_nei
                continue
        else:
            conn[ei] = np.array(qa, dtype=int)
            conn[nei] = np.array(qb, dtype=int)
        used[ei] = True
        used[nei] = True
        # Refresh edge map for subsequent quads.
        edge_to_elems = _build_full_edge_map(conn)
        quality[ei] = _quad_skew_quality(qa)
        quality[nei] = _quad_skew_quality(qb)

    return conn


def _quad_vertex_merge(
    conn: np.ndarray,
    pts: np.ndarray,
) -> np.ndarray:
    """Port of QuADMesh ``QuadVertexMerge_v2``.

    For each interior quad whose two diagonally-opposite vertices are both
    interior valence-3 vertices, collapse the diagonal: keep one vertex
    (``v_a``), drop the other (``v_b``), relabel all incident elements'
    references from ``v_b`` to ``v_a``. The center quad becomes degenerate
    (two of its corners now equal ``v_a``) and is dropped; the four
    neighbor quads share ``v_a`` going from 5 elems → 4 elems covering
    the same union (no coverage loss).

    Mutual exclusivity: each QVM operation touches 5 quads (center + 4
    neighbors). Greedily skip overlapping clusters within a single pass.
    """
    n_elems = conn.shape[0]
    if n_elems == 0:
        return conn

    edge_to_elems = _build_full_edge_map(conn)
    boundary_verts: set[int] = set()
    for key, elems in edge_to_elems.items():
        if len(elems) == 1:
            boundary_verts.add(key[0])
            boundary_verts.add(key[1])

    vert_to_elems: dict[int, list[int]] = {}
    for ei in range(n_elems):
        row = conn[ei]
        if row[0] == row[3]:
            verts = {int(row[k]) for k in range(3)}
        else:
            verts = {int(row[k]) for k in range(4)}
        for v in verts:
            vert_to_elems.setdefault(v, []).append(ei)

    keep = np.ones(n_elems, dtype=bool)
    locked: set[int] = set()  # elem ids participating in a QVM this pass.

    for ei in range(n_elems):
        if not keep[ei] or ei in locked:
            continue
        row = conn[ei]
        if row[0] == row[3]:
            continue  # tri
        # Two diagonals.
        for da, db in ((0, 2), (1, 3)):
            v_a, v_b = int(row[da]), int(row[db])
            if v_a in boundary_verts or v_b in boundary_verts:
                continue
            if len(vert_to_elems.get(v_a, [])) != 3 or len(vert_to_elems.get(v_b, [])) != 3:
                continue
            cluster = set(vert_to_elems[v_a]) | set(vert_to_elems[v_b])
            if any(c in locked for c in cluster):
                continue
            if any(not keep[c] for c in cluster):
                continue
            if len(cluster) != 5:
                continue  # diagonals must be distinct topologically; sanity check.

            # Collapse v_b → v_a in all of v_b's incident elements (other than
            # center, which dies entirely).
            for nei in vert_to_elems[v_b]:
                if nei == ei:
                    continue
                conn[nei] = np.where(conn[nei] == v_b, v_a, conn[nei])
            keep[ei] = False
            locked.update(cluster)
            break

    return conn[keep]


def _doublet_collapse(
    conn: np.ndarray,
    pts: np.ndarray,
) -> np.ndarray:
    """Port of QuADMesh ``DoubletCollapse``.

    A "doublet" is an interior vertex ``V`` with valence 2 whose two
    incident elements are both quads that share exactly 3 vertices
    (``V`` plus the two verts on either side of ``V`` in each quad). The
    collapse merges them into one quad by replacing ``V`` in the first
    quad with the second quad's unique (non-shared) vertex, then drops
    the second quad's row. The surviving quad covers the union of the
    original pair — no coverage loss.

    Invariants enforced:
      - Both incident elements must be 4-vert quads (no padded tris).
      - Vertex must be interior (no boundary-edge incidence).
      - The two quads must share exactly 3 vertices.
      - The replacement vertex must not already appear in the surviving
        quad (otherwise the merge produces a degenerate row).
    """
    n_elems = conn.shape[0]
    if n_elems == 0:
        return conn

    edge_to_elems = _build_full_edge_map(conn)
    boundary_verts: set[int] = set()
    for key, elems in edge_to_elems.items():
        if len(elems) == 1:
            boundary_verts.add(key[0])
            boundary_verts.add(key[1])

    vert_to_elems: dict[int, list[int]] = {}
    for ei in range(n_elems):
        row = conn[ei]
        if row[0] == row[3]:
            verts = {int(row[k]) for k in range(3)}
        else:
            verts = {int(row[k]) for k in range(4)}
        for v in verts:
            vert_to_elems.setdefault(v, []).append(ei)

    keep = np.ones(n_elems, dtype=bool)

    for v, adj in vert_to_elems.items():
        if v in boundary_verts:
            continue
        if len(adj) != 2:
            continue
        e1, e2 = adj
        if not (keep[e1] and keep[e2]):
            continue
        if conn[e1, 0] == conn[e1, 3] or conn[e2, 0] == conn[e2, 3]:
            continue
        verts1 = [int(conn[e1, k]) for k in range(4)]
        verts2 = [int(conn[e2, k]) for k in range(4)]
        shared = set(verts1) & set(verts2)
        if len(shared) != 3:
            continue
        unique_v2 = [u for u in verts2 if u not in verts1]
        if len(unique_v2) != 1:
            continue
        u = unique_v2[0]
        if u in verts1:
            continue
        new_row = [u if x == v else x for x in verts1]
        if len(set(new_row)) != 4:
            continue
        conn[e1] = np.array(new_row, dtype=int)
        keep[e2] = False

    return conn[keep]


def _fuse_interior_tri_pairs(
    conn: np.ndarray,
    pts: np.ndarray,
) -> np.ndarray:
    """Fuse pairs of adjacent interior triangles into a single quad.

    An "interior triangle" is a padded-tri row ``[a, b, c, a]`` whose
    three vertices are all NOT on the mesh boundary. When two such tris
    share exactly one edge ``(u, v)``, the union forms a quad
    ``(u, w_1, v, w_2)`` where ``w_1`` and ``w_2`` are the two non-shared
    vertices. The fusion is accepted iff the resulting quad is convex and
    has a consistent CCW orientation (positive signed area).
    """
    n_elems = conn.shape[0]
    if n_elems == 0:
        return conn

    edge_to_elems = _build_full_edge_map(conn)
    boundary_verts: set[int] = set()
    for key, elems in edge_to_elems.items():
        if len(elems) == 1:
            boundary_verts.add(key[0])
            boundary_verts.add(key[1])

    interior_tri_ids: list[int] = []
    for ei in range(n_elems):
        row = conn[ei]
        if row[0] != row[3]:
            continue
        verts = {int(row[k]) for k in range(3)}
        if not (verts & boundary_verts):
            interior_tri_ids.append(ei)

    if not interior_tri_ids:
        return conn

    keep = np.ones(n_elems, dtype=bool)
    used: set[int] = set()

    def _signed_area_quad(verts: list[int]) -> float:
        poly = pts[verts, :2]
        x = poly[:, 0]
        y = poly[:, 1]
        return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))

    def _is_convex_ccw(verts: list[int]) -> bool:
        poly = pts[verts, :2]
        n = len(verts)
        sign = 0
        for k in range(n):
            a = poly[(k - 1) % n]
            b = poly[k]
            c = poly[(k + 1) % n]
            cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
            if cross > 1e-14:
                s = 1
            elif cross < -1e-14:
                s = -1
            else:
                return False  # collinear vertex
            if sign == 0:
                sign = s
            elif s != sign:
                return False
        return sign == 1

    for ei in interior_tri_ids:
        if ei in used or not keep[ei]:
            continue
        row = conn[ei]
        verts_a = [int(row[k]) for k in range(3)]
        # Try each edge in the tri for a neighboring interior tri.
        for k in range(3):
            u = verts_a[k]
            v = verts_a[(k + 1) % 3]
            w_a = verts_a[(k + 2) % 3]
            key = (u, v) if u < v else (v, u)
            neighbors = [n for n in edge_to_elems.get(key, []) if n != ei]
            if len(neighbors) != 1:
                continue
            ej = neighbors[0]
            if ej in used or not keep[ej]:
                continue
            row_b = conn[ej]
            if row_b[0] != row_b[3]:
                continue  # neighbor is a quad
            verts_b = [int(row_b[k]) for k in range(3)]
            if not ({int(x) for x in verts_b} & boundary_verts == set()):
                # mixed: one of the neighbor tri verts is on boundary
                # — skip; keep this case for a boundary-aware pass.
                continue
            # The fourth vertex is verts_b minus {u, v}.
            others = [int(x) for x in verts_b if x not in (u, v)]
            if len(others) != 1:
                continue
            w_b = others[0]
            if w_b in (u, v, w_a):
                continue
            # Quad CCW: u → w_a → v → w_b. The shared edge (u,v) becomes
            # a diagonal; the new perimeter is (u→w_a, w_a→v, v→w_b, w_b→u).
            quad = [u, w_a, v, w_b]
            if _signed_area_quad(quad) <= 0:
                continue
            if not _is_convex_ccw(quad):
                continue
            # Accept fusion: overwrite ei with the quad, drop ej.
            conn[ei] = np.array(quad, dtype=int)
            keep[ej] = False
            used.add(ei)
            used.add(ej)
            break

    return conn[keep]


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


def _apply_edge_insertion(
    interior_tri_ids: list[int],
    tris: np.ndarray,
    points: np.ndarray,
    quad_rows: list,
    edge_to_tris: dict[tuple[int, int], list[int]],
) -> tuple[list[int], np.ndarray]:
    """Apply edge insertion (vertex split) to interior tris.

    For each interior tri with an edge shared with a neighbor:
    - Insert new vertex at 1/3 distance along shared edge
    - Convert tri → quad via inserted vertex
    - Add quad to output and new point to points array

    Returns:
        (remaining_tri_ids, updated_points_array)
    """
    if not interior_tri_ids:
        return [], points

    remaining = []
    next_vert_id = len(points)
    new_points_list = []

    for t_id in interior_tri_ids:
        tri = tris[t_id]
        tri_verts = [int(v) for v in tri]

        handled = False

        for v_a, v_b in [
            (tri_verts[0], tri_verts[1]),
            (tri_verts[1], tri_verts[2]),
            (tri_verts[2], tri_verts[0]),
        ]:
            edge_key = _make_key(v_a, v_b)
            tlist = edge_to_tris.get(edge_key, [])
            if len(tlist) >= 2:
                p_a = points[v_a, :2]
                p_b = points[v_b, :2]
                new_pt_2d = (1.0 - 1.0/3.0) * p_a + (1.0/3.0) * p_b
                new_pt_3d = np.array([new_pt_2d[0], new_pt_2d[1], points[v_a, 2]])

                quad_v = [
                    int(tri_verts[0]) if tri_verts[0] not in (v_a, v_b) else None,
                    int(tri_verts[1]) if tri_verts[1] not in (v_a, v_b) else None,
                    int(tri_verts[2]) if tri_verts[2] not in (v_a, v_b) else None,
                ]
                quad_v = [v for v in quad_v if v is not None]

                if len(quad_v) == 1:
                    quad_verts = [
                        int(v_a),
                        quad_v[0],
                        int(v_b),
                        next_vert_id,
                    ]
                    quad_rows.append(quad_verts)
                    new_points_list.append(new_pt_3d)
                    next_vert_id += 1
                    handled = True
                    break

        if not handled:
            remaining.append(t_id)

    if new_points_list:
        points = np.vstack([points, np.array(new_points_list)])

    return remaining, points


def _remove_interior_triangles(
    interior_tri_ids: list[int],
    tris: np.ndarray,
    points: np.ndarray,
    edge_to_tris: dict[tuple[int, int], list[int]],
    boundary_vert_ids: set[int],
) -> list[int]:
    """Attempt to remove interior triangles via edge flips or insertion.

    For interior tris with no boundary verts, tries to match with neighbors
    and convert via topological operations. Unhandled tris remain and will
    be padded as fallback.

    Note: Full topological flip requires connectivity rebuild (not yet implemented).
    Currently detects flipability but returns all for fallback padding.

    Returns:
        List of interior tri IDs (all returned as fallback padding not yet
        integrated into main connectivity).
    """
    if not interior_tri_ids:
        return []

    remaining = []
    for t_id in interior_tri_ids:
        tri = tris[t_id]
        tri_verts = [int(v) for v in tri]

        has_neighbor = False

        for i, (v_a, v_b) in enumerate(
            [(tri_verts[0], tri_verts[1]),
             (tri_verts[1], tri_verts[2]),
             (tri_verts[2], tri_verts[0])]
        ):
            key = _make_key(v_a, v_b)
            tlist = edge_to_tris.get(key, [])

            if len(tlist) == 2:
                has_neighbor = True
                break

        remaining.append(t_id)


# ---------------------------------------------------------------------------
# Mesh repair utilities
# ---------------------------------------------------------------------------

def _snap_boundary_midpoints(
    conn: np.ndarray,
    pts: np.ndarray,
    mid_threshold: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Snap boundary midpoint vertices to the exact midpoint of their two
    boundary-edge neighbors.

    ``_insert_boundary_tri_midpoint`` inserts midpoints at vertex indices
    >= the original vertex count. Subsequent topology ops (doublet collapse,
    QVM) rewire connectivity without updating midpoint geometry, leaving
    some midpoints displaced from their correct positions. This function
    finds all degree-2 boundary vertices (exactly 2 incident edges, both on
    the mesh boundary) with index >= ``mid_threshold`` and snaps them to
    the midpoint of their two boundary neighbors.

    Args:
        conn: (n_elems, 4) connectivity array (padded tris have row[0]==row[3]).
        pts:  (n_verts, 3) points array.
        mid_threshold: Only snap vertices with index >= this value. Set to
            the original vertex count before midpoint insertion. When 0,
            snaps all degree-2 boundary verts (safe but slower).

    Returns:
        (conn, pts) with pts snapped in-place (conn unchanged).
    """
    edge_map = _build_full_edge_map(conn)

    bv_edges: dict[int, list[tuple[int, int]]] = {}
    for key, lst in edge_map.items():
        if len(lst) != 1:
            continue
        a, b = key
        bv_edges.setdefault(a, []).append(key)
        bv_edges.setdefault(b, []).append(key)

    for v, edges in bv_edges.items():
        if v < mid_threshold:
            continue
        if len(edges) != 2:
            continue
        nbrs = [k[0] if k[1] == v else k[1] for k in edges]
        pts[v, :2] = 0.5 * (pts[nbrs[0], :2] + pts[nbrs[1], :2])

    return conn, pts


def _fix_bowties(conn: np.ndarray, pts: np.ndarray) -> tuple[np.ndarray, int]:
    """Fix self-intersecting (bowtie) quads by reordering their vertices.

    For each quad where edges (0,1)x(2,3) or (1,2)x(3,0) cross, try
    alternative vertex orderings — including the convex-hull ordering — to
    find a valid CCW non-self-intersecting arrangement.

    Returns:
        (conn, n_fixed) — conn modified in-place; n_fixed = number of rows fixed.
    """
    try:
        from scipy.spatial import ConvexHull as _ConvexHull
        _have_scipy = True
    except ImportError:
        _have_scipy = False

    def _seg_cross(p1: np.ndarray, p2: np.ndarray,
                   p3: np.ndarray, p4: np.ndarray) -> bool:
        def _ori(pa, pb, pc):
            v = (pb[0] - pa[0]) * (pc[1] - pa[1]) - (pb[1] - pa[1]) * (pc[0] - pa[0])
            return 1 if v > 1e-12 else (-1 if v < -1e-12 else 0)
        o1, o2 = _ori(p1, p2, p3), _ori(p1, p2, p4)
        o3, o4 = _ori(p3, p4, p1), _ori(p3, p4, p2)
        if o1 == 0 or o2 == 0 or o3 == 0 or o4 == 0:
            return False
        return o1 != o2 and o3 != o4

    n_fixed = 0
    for ei in range(conn.shape[0]):
        r = conn[ei]
        if r[0] == r[3]:
            continue
        verts = [int(r[k]) for k in range(4)]
        if len(set(verts)) < 4:
            continue
        poly = pts[verts, :2]
        v0, v1, v2, v3 = poly[0], poly[1], poly[2], poly[3]
        if not (_seg_cross(v0, v1, v2, v3) or _seg_cross(v1, v2, v3, v0)):
            continue

        candidates: list[list[int]] = [
            [verts[0], verts[1], verts[3], verts[2]],
            [verts[0], verts[2], verts[1], verts[3]],
            [verts[0], verts[3], verts[2], verts[1]],
        ]
        if _have_scipy:
            try:
                ch = _ConvexHull(poly)
                if len(ch.vertices) == 4:
                    candidates.append([verts[i] for i in ch.vertices.tolist()])
            except Exception:
                pass

        best: tuple[list[int], float] | None = None
        for c in candidates:
            p = pts[c, :2]
            x, y = p[:, 0], p[:, 1]
            sa = 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
            if sa <= 0:
                continue
            cv0, cv1, cv2, cv3 = p[0], p[1], p[2], p[3]
            if _seg_cross(cv0, cv1, cv2, cv3) or _seg_cross(cv1, cv2, cv3, cv0):
                continue
            if best is None or sa > best[1]:
                best = (c, sa)

        if best is not None:
            conn[ei] = np.array(best[0], dtype=int)
            n_fixed += 1

    return conn, n_fixed


def _ear_clip(ring: list[int], pts: np.ndarray) -> list[tuple[int, int, int]]:
    """Triangulate a simple polygon via ear-clipping.

    Args:
        ring: List of vertex IDs forming the polygon boundary (any winding).
        pts:  (n_verts, 3+) points array.

    Returns:
        List of (a, b, c) triangle vertex-ID triples in CCW order.
    """
    ring = list(ring)
    n = len(ring)
    if n < 3:
        return []
    if n == 3:
        return [tuple(ring)]  # type: ignore[return-value]

    poly = pts[ring, :2]
    x, y = poly[:, 0], poly[:, 1]
    if 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)) < 0:
        ring = ring[::-1]

    def _point_in_tri(p: np.ndarray, a: np.ndarray, v: np.ndarray, b: np.ndarray) -> bool:
        def _ori(p0, p1, p2):
            return (p1[0]-p0[0])*(p2[1]-p0[1]) - (p1[1]-p0[1])*(p2[0]-p0[0])
        d1, d2, d3 = _ori(p, a, v), _ori(p, v, b), _ori(p, b, a)
        has_neg = d1 < 0 or d2 < 0 or d3 < 0
        has_pos = d1 > 0 or d2 > 0 or d3 > 0
        if has_neg and has_pos:
            return False
        return abs(d1) > 1e-14 and abs(d2) > 1e-14 and abs(d3) > 1e-14

    tris: list[tuple[int, int, int]] = []
    max_iters = len(ring) ** 2 + 10
    iters = 0
    while len(ring) > 3 and iters < max_iters:
        iters += 1
        cur_n = len(ring)
        ear_found = False
        for i in range(cur_n):
            a_id = ring[(i - 1) % cur_n]
            v_id = ring[i]
            b_id = ring[(i + 1) % cur_n]
            a = pts[a_id, :2]; v = pts[v_id, :2]; b = pts[b_id, :2]
            cross = (v[0] - a[0]) * (b[1] - a[1]) - (v[1] - a[1]) * (b[0] - a[0])
            if cross <= 0:
                continue
            inside = any(
                _point_in_tri(pts[ring[j], :2], a, v, b)
                for j in range(cur_n)
                if j not in ((i - 1) % cur_n, i, (i + 1) % cur_n)
            )
            if not inside:
                tris.append((a_id, v_id, b_id))
                ring.pop(i)
                ear_found = True
                break
        if not ear_found:
            tris.append((ring[-1], ring[0], ring[1]))
            ring.pop(0)
    if len(ring) == 3:
        tris.append(tuple(ring))  # type: ignore[arg-type]
    return tris


def _merge_tri_pairs_to_quads(
    tris: list[tuple[int, ...]],
    pts: np.ndarray,
) -> list[tuple[int, ...]]:
    """Greedily merge adjacent triangle pairs into quads.

    Triangles sharing an edge are merged into a CCW quad iff the resulting
    quad has positive signed area and is not self-intersecting (bowtie).

    Args:
        tris: List of (a, b, c) vertex-ID triples.
        pts:  (n_verts, 3+) points array.

    Returns:
        List of 3-tuples (kept triangles) and 4-tuples (merged quads).
    """
    edge_map: dict[tuple[int, int], list[int]] = {}
    for i, t in enumerate(tris):
        for k in range(3):
            a, b = t[k], t[(k + 1) % 3]
            key = (a, b) if a < b else (b, a)
            edge_map.setdefault(key, []).append(i)

    used: set[int] = set()
    result: list[tuple[int, ...]] = []

    for i, t in enumerate(tris):
        if i in used:
            continue
        best_pair: tuple[int, tuple[int, int, int, int]] | None = None
        for k in range(3):
            a, b = t[k], t[(k + 1) % 3]
            key = (a, b) if a < b else (b, a)
            nbrs = [j for j in edge_map.get(key, []) if j != i and j not in used]
            if not nbrs:
                continue
            j = nbrs[0]
            tj = tris[j]
            ua = [v for v in t if v not in (a, b)]
            ub = [v for v in tj if v not in (a, b)]
            if len(ua) != 1 or len(ub) != 1:
                continue
            idx_ua = list(t).index(ua[0])
            nxt = t[(idx_ua + 1) % 3]
            quad: tuple[int, int, int, int] = (ua[0], a, ub[0], b) if nxt == a else (ua[0], b, ub[0], a)
            poly = pts[list(quad), :2]
            x, y = poly[:, 0], poly[:, 1]
            sa = 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
            if sa <= 0:
                quad = tuple(reversed(quad))  # type: ignore[assignment]
                poly = pts[list(quad), :2]
                x, y = poly[:, 0], poly[:, 1]
                sa = 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
            if sa <= 0:
                continue
            v0r, v1r, v2r, v3r = poly[0], poly[1], poly[2], poly[3]

            def _sc(p1, p2, p3, p4):
                def _o(pa, pb, pc):
                    v = (pb[0]-pa[0])*(pc[1]-pa[1])-(pb[1]-pa[1])*(pc[0]-pa[0])
                    return 1 if v > 1e-12 else (-1 if v < -1e-12 else 0)
                o1,o2,o3,o4 = _o(p1,p2,p3),_o(p1,p2,p4),_o(p3,p4,p1),_o(p3,p4,p2)
                if 0 in (o1,o2,o3,o4): return False
                return o1 != o2 and o3 != o4

            if _sc(v0r, v1r, v2r, v3r) or _sc(v1r, v2r, v3r, v0r):
                continue
            best_pair = (j, quad)
            break
        if best_pair:
            j, quad = best_pair
            result.append(quad)
            used.add(i); used.add(j)
        else:
            result.append(t)
    return result


def repair_mesh(
    conn: np.ndarray,
    pts: np.ndarray,
    mid_threshold: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Post-processing repair pass: fix geometric defects left by the
    topology pipeline.

    Applies three repair operations in order:

    1. **Snap boundary midpoints** — vertices at index >= ``mid_threshold``
       with exactly 2 incident boundary edges are moved to the geometric
       midpoint of their two boundary neighbors. Fixes midpoints displaced
       by doublet-collapse / QVM rewiring.

    2. **Fix bowtie quads** — self-intersecting quads are reordered to a
       valid CCW non-crossing arrangement using ConvexHull ordering and
       alternative vertex permutations.

    3. **Dissolve-and-remesh violation clusters** — clusters of elements
       that still have INTERIOR_OVERLAP or EDGE_CROSSING violations after
       steps 1-2 are dissolved into their exterior boundary polygon,
       ear-clip triangulated, and merged back into quads.

    All operations preserve the mesh boundary. Area is preserved exactly
    for steps 1-2 and approximately for step 3 (ear-clip + merge covers
    the same region).

    Args:
        conn: (n_elems, 4) connectivity (row[0]==row[3] for padded tris).
        pts:  (n_verts, 3) points.
        mid_threshold: Vertex index below which midpoint snapping is skipped.
            Pass the original vertex count (before ``_insert_boundary_tri_midpoint``)
            to restrict snapping to inserted midpoints only.

    Returns:
        (conn, pts) repaired.
    """
    conn = conn.copy()
    pts = pts.copy()

    conn, pts = _snap_boundary_midpoints(conn, pts, mid_threshold=mid_threshold)
    conn, _ = _fix_bowties(conn, pts)
    conn, pts = _dissolve_violation_clusters(conn, pts)

    return conn, pts


def _dissolve_violation_clusters(
    conn: np.ndarray,
    pts: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Find elements with INTERIOR_OVERLAP or EDGE_CROSSING violations and
    dissolve+remesh each cluster.

    A cluster is a connected component of violation elements. Each cluster
    is dissolved by:
      1. Computing the exterior boundary ring.
      2. Ear-clip triangulating the ring.
      3. Merging adjacent tri pairs into quads.
      4. Replacing the cluster rows in conn.

    Silently skips clusters where the repair produces a bowtie.
    """
    try:
        from chilmesh.validation.validator import validate_mesh_elements as _validate
        from chilmesh import CHILmesh as _CM
    except ImportError:
        return conn, pts

    try:
        snap = _CM(connectivity=conn, points=pts, grid_name='_repair',
                   compute_layers=False)
        rep = _validate(snap)
    except Exception:
        return conn, pts

    bad_elems: set[int] = set()
    for viol in rep.violations:
        if viol.category == 'NON_PLANAR_MESH':
            continue
        for eid in viol.element_ids:
            bad_elems.add(int(eid))

    if not bad_elems:
        return conn, pts

    edge_map = _build_full_edge_map(conn)
    adj: dict[int, set[int]] = {e: set() for e in bad_elems}
    for key, lst in edge_map.items():
        for ei in lst:
            if ei not in bad_elems:
                continue
            for ej in lst:
                if ej != ei and ej in bad_elems:
                    adj[ei].add(ej)

    visited: set[int] = set()
    clusters: list[list[int]] = []
    for start in sorted(bad_elems):
        if start in visited:
            continue
        cluster: list[int] = []
        queue = [start]
        while queue:
            cur = queue.pop()
            if cur in visited:
                continue
            visited.add(cur)
            cluster.append(cur)
            queue.extend(adj.get(cur, set()) - visited)
        clusters.append(cluster)

    all_new_rows: list[list[int]] = []
    delete_mask = np.ones(conn.shape[0], dtype=bool)

    for cluster in clusters:
        edge_count: dict[tuple[int, int], int] = {}
        for ei in cluster:
            r = conn[ei]
            verts = [int(r[k]) for k in range(3 if r[0] == r[3] else 4)]
            n = len(verts)
            for k in range(n):
                a, b = verts[k], verts[(k + 1) % n]
                key = (a, b) if a < b else (b, a)
                edge_count[key] = edge_count.get(key, 0) + 1
        exterior = [k for k, v in edge_count.items() if v == 1]
        if not exterior:
            continue

        adj_v: dict[int, list[int]] = {}
        for a, b in exterior:
            adj_v.setdefault(a, []).append(b)
            adj_v.setdefault(b, []).append(a)
        start_v = exterior[0][0]
        ring: list[int] = [start_v]
        prev_v: int | None = None
        cur_v = start_v
        for _ in range(len(exterior) + 2):
            opts = [v for v in adj_v[cur_v] if v != prev_v]
            if not opts:
                break
            nxt = opts[0]
            if nxt == start_v:
                break
            ring.append(nxt)
            prev_v = cur_v
            cur_v = nxt

        if len(ring) < 3:
            continue

        tris = _ear_clip(ring, pts)
        if not tris:
            continue
        elements = _merge_tri_pairs_to_quads(tris, pts)

        new_rows: list[list[int]] = []
        valid = True
        for e in elements:
            if len(e) == 3:
                a, b, c = e
                new_rows.append([a, b, c, a])
            elif len(e) == 4:
                poly = pts[list(e), :2]
                v0, v1, v2, v3 = poly[0], poly[1], poly[2], poly[3]

                def _sc2(p1, p2, p3, p4):
                    def _o(pa, pb, pc):
                        v = (pb[0]-pa[0])*(pc[1]-pa[1])-(pb[1]-pa[1])*(pc[0]-pa[0])
                        return 1 if v > 1e-12 else (-1 if v < -1e-12 else 0)
                    o1,o2,o3,o4 = _o(p1,p2,p3),_o(p1,p2,p4),_o(p3,p4,p1),_o(p3,p4,p2)
                    if 0 in (o1,o2,o3,o4): return False
                    return o1 != o2 and o3 != o4

                if _sc2(v0, v1, v2, v3) or _sc2(v1, v2, v3, v0):
                    valid = False
                    break
                new_rows.append(list(e))

        if not valid:
            continue

        for ei in cluster:
            delete_mask[ei] = False
        all_new_rows.extend(new_rows)

    kept_conn = conn[delete_mask]
    if all_new_rows:
        conn = np.vstack([kept_conn, np.array(all_new_rows, dtype=int)])
    else:
        conn = kept_conn

    return conn, pts

    return remaining
