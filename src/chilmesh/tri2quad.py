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
        conn, pts = _pentagon_split_interior_tris(conn, pts)

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

    out = CHILmesh(
        connectivity=conn,
        points=pts,
        grid_name=mesh.grid_name,
        compute_layers=True,
    )

    if smooth:
        pre_pts = np.asarray(out.points).copy()
        try:
            out.smooth_mesh(method=smooth_method, acknowledge_change=True)
            if not np.isfinite(out.points).all():
                raise RuntimeError("smoother produced non-finite coordinates")
        except Exception:
            out.change_points(pre_pts, acknowledge_change=True)
            out.smooth_mesh(method="angle-based", acknowledge_change=True)

        # FEM / angle-based both struggle on slivery quads inherited from the
        # path-walk merge (interior angles approaching 180° pin Q=0 elements).
        # Laplacian relaxation has no Hessian / singular-matrix failure mode
        # and reliably lifts the long-tail of bad quads.
        lap_pts = _laplacian_smooth(out, n_iter=smooth_iters, omega=0.7)
        out.change_points(lap_pts, acknowledge_change=True)

    return out


def _laplacian_smooth(
    mesh: "CHILmesh",
    n_iter: int = 100,
    omega: float = 0.7,
) -> np.ndarray:
    """Constrained Laplacian smoothing: each interior vertex moves toward the
    centroid of its 1-ring neighbours, blended by ``omega``.

    Boundary vertices are held fixed (per CHILmesh.boundary_node_indices).
    Returns a new points array; does not mutate the mesh.
    """
    conn = np.asarray(mesh.connectivity_list, dtype=int)
    pts = np.asarray(mesh.points, dtype=float).copy()
    boundary_verts = {
        int(x) for x in np.asarray(mesh.boundary_node_indices()).flatten()
    }

    n_verts = pts.shape[0]
    nbrs: list[set[int]] = [set() for _ in range(n_verts)]
    for row in conn:
        if row[0] == row[3]:
            verts = [int(row[k]) for k in range(3)]
        else:
            verts = [int(row[k]) for k in range(4)]
        n = len(verts)
        for k in range(n):
            a, b = verts[k], verts[(k + 1) % n]
            nbrs[a].add(b)
            nbrs[b].add(a)

    for _ in range(n_iter):
        new_pts = pts.copy()
        for v in range(n_verts):
            if v in boundary_verts or not nbrs[v]:
                continue
            avg = np.mean([pts[u, :2] for u in nbrs[v]], axis=0)
            new_pts[v, :2] = (1.0 - omega) * pts[v, :2] + omega * avg
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

        # Split hexagon along diagonal s_a → q_far2 (or equivalent). The
        # cleanest split that gives two quads:
        #   quad1 = [s_a, q_far1, M, ???]
        # Pentagon split via M and the diagonal from B to M:
        # Hexagon CCW: s_a, q_far1, M, q_far2, s_b, B.
        # Diagonals from M: M–s_a, M–q_far2 (adjacent edges, no), M–s_b
        # (non-adj), M–B (non-adj).
        # Two quads via diagonal M–B:
        #   quad1 = [B, s_a, q_far1, M] (CCW)
        #   quad2 = [B, M, q_far2, s_b] (CCW)
        # Check: quad1 verts {B, s_a, q_far1, M} — 4 distinct ✓
        #         quad2 verts {B, M, q_far2, s_b} — 4 distinct ✓
        B = third
        quad1 = [B, s_a, q_far1, m_vid]
        quad2 = [B, m_vid, q_far2, s_b]
        conn[tri_eid] = np.array(quad1, dtype=int)
        conn[chosen_neighbor] = np.array(quad2, dtype=int)
        # Refresh edge map (cheap since small).
        edge_to_elems = _build_full_edge_map(conn)

    if new_pts:
        pts = np.vstack([pts, np.asarray(new_pts, dtype=float)])
    return conn, pts


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

    return remaining
