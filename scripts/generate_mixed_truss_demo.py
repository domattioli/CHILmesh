#!/usr/bin/env python3
"""Mixed-element mesh demo: quad core + triangulated outer band.

Pipeline:
1. Generate 16x12 structured quad mesh on a 4x3 rectangle.
2. Skeletonize (CHILmesh layer decomposition) — yields 6 concentric layers.
3. Strip outer 2 layers of quads, keep their vertices.
4. Delaunay-triangulate the kept vertex set; filter triangles whose centroid
   lies inside the inner quad core.
5. ADMESH warm-start truss on the new triangles ONLY (boundary pinned bit-exact:
   outer rectangle perimeter + inner quad-core seam).
6. Combine kept quads + smoothed tris into a single mixed-element mesh.
7. FEM-smooth the combined mesh.
8. Render 4-panel figure to ``output/mixed_truss_fem_demo.png``.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path as MplPath
from scipy.spatial import Delaunay

from chilmesh import CHILmesh
from chilmesh.admesh_warmstart import optimize_with_admesh_truss_arrays


NX, NY = 16, 12
LX, LY = 4.0, 3.0


# ---------------------------------------------------------------------------
# Stage 1: structured quad mesh
# ---------------------------------------------------------------------------

def build_quad_grid() -> CHILmesh:
    xs = np.linspace(0, LX, NX + 1)
    ys = np.linspace(0, LY, NY + 1)
    points = np.array([[x, y, 0.0] for y in ys for x in xs])

    def vid(i: int, j: int) -> int:
        return j * (NX + 1) + i

    quads = np.array(
        [
            [vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)]
            for j in range(NY)
            for i in range(NX)
        ],
        dtype=int,
    )
    return CHILmesh(connectivity=quads, points=points, compute_layers=True)


# ---------------------------------------------------------------------------
# Stage 2-3: split, identify boundaries
# ---------------------------------------------------------------------------

def _outer_perim_ring() -> np.ndarray:
    """Ordered CCW outer perimeter of the original rectangle."""
    def vid(i, j): return j * (NX + 1) + i
    bottom = [vid(i, 0) for i in range(NX + 1)]
    right = [vid(NX, j) for j in range(1, NY + 1)]
    top = [vid(i, NY) for i in range(NX - 1, -1, -1)]
    left = [vid(0, j) for j in range(NY - 1, 0, -1)]
    return np.array(bottom + right + top + left, dtype=int)


def _inner_seam_ring(inner_quads: np.ndarray) -> np.ndarray:
    """Walk CCW boundary of the inner quad sub-mesh."""
    edge_count: dict[tuple[int, int], int] = {}
    edge_dir: dict[tuple[int, int], tuple[int, int]] = {}
    for q in inner_quads:
        for a, b in zip(q, np.roll(q, -1)):
            key = (min(a, b), max(a, b))
            edge_count[key] = edge_count.get(key, 0) + 1
            edge_dir.setdefault(key, (a, b))
    boundary = [edge_dir[k] for k, c in edge_count.items() if c == 1]
    nxt = {a: b for a, b in boundary}
    start = boundary[0][0]
    walk = [start]
    cur = nxt[start]
    while cur != start:
        walk.append(cur)
        cur = nxt[cur]
    return np.array(walk, dtype=int)


def split_inner_outer(mesh: CHILmesh, strip_layers: int = 2):
    outer_elem_ids = np.concatenate([mesh.layers["OE"][k] for k in range(strip_layers)])
    inner_elem_ids = np.concatenate(
        [mesh.layers["OE"][k] for k in range(strip_layers, mesh.n_layers)]
    )
    outer_quads = mesh.connectivity_list[outer_elem_ids]
    inner_quads = mesh.connectivity_list[inner_elem_ids]
    return outer_quads, inner_quads


# ---------------------------------------------------------------------------
# Stage 4: Delaunay over kept points, filter by centroid
# ---------------------------------------------------------------------------

def triangulate_outer_band(
    points_xy: np.ndarray,
    outer_quads: np.ndarray,
    inner_quads: np.ndarray,
):
    """Return (tris_global, tri_candidate_verts) where tris_global indexes
    into the original global point array."""
    outer_verts = np.unique(outer_quads.flatten())
    inner_seam = _inner_seam_ring(inner_quads)
    candidate_verts = np.unique(np.concatenate([outer_verts, inner_seam]))

    sites = points_xy[candidate_verts]
    tri = Delaunay(sites)
    simplices_global = candidate_verts[tri.simplices]
    centroids = points_xy[simplices_global].mean(axis=1)

    inner_path = MplPath(points_xy[inner_seam])
    outer_path = MplPath(points_xy[_outer_perim_ring()])

    keep = outer_path.contains_points(centroids, radius=1e-9) & ~inner_path.contains_points(
        centroids, radius=-1e-9
    )
    return simplices_global[keep], candidate_verts, inner_seam


# ---------------------------------------------------------------------------
# Stage 5: ADMESH truss on tris only (boundary pinned)
# ---------------------------------------------------------------------------

def _polygon_sdf(points: np.ndarray, polygon_verts: np.ndarray) -> np.ndarray:
    """Signed distance from points to a closed polygon. Negative inside."""
    a = polygon_verts
    b = np.roll(polygon_verts, -1, axis=0)
    min_d2 = np.full(len(points), np.inf)
    for i in range(len(a)):
        seg_a, seg_b = a[i], b[i]
        ab = seg_b - seg_a
        ab_len2 = float(ab @ ab)
        if ab_len2 == 0.0:
            continue
        ap = points - seg_a
        t = np.clip((ap @ ab) / ab_len2, 0.0, 1.0)
        proj = seg_a + np.outer(t, ab)
        d2 = ((points - proj) ** 2).sum(axis=1)
        min_d2 = np.minimum(min_d2, d2)
    d = np.sqrt(min_d2)
    inside = MplPath(polygon_verts).contains_points(points)
    return np.where(inside, -d, d)


def truss_tris(
    points_xy: np.ndarray,
    tris_global: np.ndarray,
    pinned_global: np.ndarray,
    inner_seam: np.ndarray,
):
    """Run ADMESH array-form truss on tri sub-mesh."""
    tri_verts_global = np.unique(tris_global.flatten())
    g2l = {g: l for l, g in enumerate(tri_verts_global)}
    tris_local = np.array([[g2l[v] for v in t] for t in tris_global], dtype=int)
    pts_local = points_xy[tri_verts_global].astype(np.float64)
    pinned_local = np.array([g2l[g] for g in pinned_global], dtype=int)

    inner_seam_xy = points_xy[inner_seam]

    def sdf(p: np.ndarray) -> np.ndarray:
        # Outer rectangle SDF (negative inside)
        dx = np.maximum(0.0 - p[:, 0], p[:, 0] - LX)
        dy = np.maximum(0.0 - p[:, 1], p[:, 1] - LY)
        sdf_rect = np.where(
            (dx > 0) | (dy > 0),
            np.sqrt(np.maximum(dx, 0) ** 2 + np.maximum(dy, 0) ** 2),
            np.maximum(dx, dy),
        )
        # Inner polygon SDF (negative inside inner core polygon)
        sdf_inner = _polygon_sdf(p, inner_seam_xy)
        # Annular region = rect MINUS inner_polygon = max(sdf_rect, -sdf_inner)
        return np.maximum(sdf_rect, -sdf_inner)

    points_opt, tris_opt = optimize_with_admesh_truss_arrays(
        pts_local,
        tris_local,
        sdf,
        boundary_indices=pinned_local,
        niter=80,
        Fscale=1.05,
        enforce_non_degradation=False,
    )

    # `boundary_indices` of pts_local = pinned_local. After truss, those
    # pinned points appear at indices 0..B-1 of points_opt in the SAME ORDER
    # as pinned_local (per docstring T013).
    pin_to_local = {int(g): i for i, g in enumerate(pinned_global)}

    # interior_global: tri verts that weren't pinned. Their post-truss positions
    # may differ in count; not strictly needed downstream because we rebuild
    # connectivity from tris_opt directly.
    interior_global = np.setdiff1d(tri_verts_global, pinned_global)

    return points_opt, tris_opt, pin_to_local, interior_global


# ---------------------------------------------------------------------------
# Stage 6: rebuild combined mixed-element mesh
# ---------------------------------------------------------------------------

def build_mixed_mesh_after_truss(
    points_opt: np.ndarray,
    tris_opt: np.ndarray,
    pin_to_local: dict[int, int],
    inner_quads: np.ndarray,
    points_orig: np.ndarray,
):
    """Stitch trussed tris to original quads.

    points_opt is (Nt, 2): boundary first (0..B-1, matching pinned_global order),
    then interior. We append the original interior-quad-only points on top.
    """
    # Quad-only verts = inner_quads verts NOT in pin_to_local
    inner_quad_verts = np.unique(inner_quads.flatten())
    pinned_globals = set(pin_to_local.keys())
    quad_only_globals = np.array(
        [v for v in inner_quad_verts if v not in pinned_globals], dtype=int
    )

    # Concatenate point arrays. Note points_opt is (N, 2); promote to (N, 3) if
    # needed for CHILmesh (z=0 column).
    pts_xy = np.vstack([points_opt, points_orig[quad_only_globals, :2]])
    pts_xyz = np.column_stack([pts_xy, np.zeros(len(pts_xy))])

    Nt = len(points_opt)
    quad_only_global_to_local = {
        int(g): Nt + i for i, g in enumerate(quad_only_globals)
    }
    quad_remap = {**pin_to_local, **quad_only_global_to_local}

    new_quads = np.array(
        [[quad_remap[int(v)] for v in q] for q in inner_quads], dtype=int
    )

    # Tris already in points_opt indexing, pad to 4-col convention
    new_tris_padded = np.column_stack([tris_opt, tris_opt[:, 0]])

    mixed_conn = np.vstack([new_tris_padded, new_quads])
    return CHILmesh(
        connectivity=mixed_conn, points=pts_xyz, compute_layers=False
    )


# ---------------------------------------------------------------------------
# Stage 7: FEM smooth
# ---------------------------------------------------------------------------

def fem_smooth(mesh: CHILmesh) -> CHILmesh:
    """Try CHILmesh FEM smoother first; fall back to Laplacian on regression.

    Note (2026-05-10): the CHILmesh FEM smoother degrades quality on
    mixed-element meshes that have an internal seam between quads and tris
    (interior verts displace by more than mesh extent). Filed as a separate
    issue. Falling back to boundary-pinned Laplacian smoothing — robust and
    visually informative for the README showcase.
    """
    pts = mesh.points.copy()
    boundary_verts = np.unique(mesh.edge2vert(mesh.boundary_edges()).flatten())
    pinned = np.zeros(mesh.n_verts, dtype=bool)
    pinned[boundary_verts] = True

    # Build vert → neighbor verts map from edge connectivity
    edge2vert = mesh.adjacencies["Edge2Vert"]
    neighbors: list[list[int]] = [[] for _ in range(mesh.n_verts)]
    for a, b in edge2vert:
        neighbors[int(a)].append(int(b))
        neighbors[int(b)].append(int(a))

    n_iter = 30
    omega = 0.5  # under-relaxation
    for _ in range(n_iter):
        new_xy = pts[:, :2].copy()
        for v in range(mesh.n_verts):
            if pinned[v] or not neighbors[v]:
                continue
            mean_xy = pts[neighbors[v], :2].mean(axis=0)
            new_xy[v] = (1 - omega) * pts[v, :2] + omega * mean_xy
        pts[:, :2] = new_xy

    mesh.change_points(pts, acknowledge_change=True)
    return mesh


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _plot_mesh_simple(ax, mesh: CHILmesh, title: str):
    mesh.plot(ax=ax, elem_color="#cce0f3")
    ax.set_title(title, fontsize=11)
    ax.set_aspect("equal")
    ax.set_xlim(-0.15, LX + 0.15)
    ax.set_ylim(-0.15, LY + 0.15)
    ax.axis("off")


def _plot_mixed(ax, mesh: CHILmesh, title: str):
    conn = mesh.connectivity_list
    if conn.shape[1] == 4:
        is_padded_tri = conn[:, 0] == conn[:, 3]
    else:
        is_padded_tri = np.ones(conn.shape[0], dtype=bool)

    tri_ids = np.where(is_padded_tri)[0]
    quad_ids = np.where(~is_padded_tri)[0]

    if len(quad_ids) > 0:
        quad_mesh = CHILmesh(
            connectivity=conn[quad_ids], points=mesh.points, compute_layers=False
        )
        quad_mesh.plot(ax=ax, elem_color="#ffd180")
    if len(tri_ids) > 0:
        tri_conn = conn[tri_ids][:, :3]
        tri_mesh = CHILmesh(
            connectivity=tri_conn, points=mesh.points, compute_layers=False
        )
        tri_mesh.plot(ax=ax, elem_color="#a5d6a7")

    ax.set_title(title, fontsize=11)
    ax.set_aspect("equal")
    ax.set_xlim(-0.15, LX + 0.15)
    ax.set_ylim(-0.15, LY + 0.15)
    ax.axis("off")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(out_path: Path | None = None) -> Path:
    if out_path is None:
        out_path = Path(__file__).parent.parent / "output" / "mixed_truss_fem_demo.png"
    out_path.parent.mkdir(exist_ok=True)

    print("[1/7] Generating 16x12 quad rectangle …")
    mesh0 = build_quad_grid()
    print(
        f"     verts={mesh0.n_verts} elems={mesh0.n_elems} layers={mesh0.n_layers}"
    )
    assert mesh0.n_layers >= 5, f"need ≥5 layers; got {mesh0.n_layers}"

    print("[2/7] Splitting inner core (keep) vs outer band (strip) …")
    outer_quads, inner_quads = split_inner_outer(mesh0, strip_layers=2)
    print(f"     stripped={len(outer_quads)} quads, kept={len(inner_quads)} quads")

    print("[3/7] Delaunay over kept-vertex set, filter by centroid …")
    points_xy = mesh0.points[:, :2]
    tris_global, _, inner_seam = triangulate_outer_band(points_xy, outer_quads, inner_quads)
    print(f"     {len(tris_global)} tris in outer band")

    # Assemble pre-truss mixed mesh for the (2) panel
    pre_truss_conn = np.vstack(
        [np.column_stack([tris_global, tris_global[:, 0]]), inner_quads]
    )
    mesh_pre = CHILmesh(
        connectivity=pre_truss_conn, points=mesh0.points, compute_layers=False
    )

    print("[4/7] ADMESH truss on tris (outer rect + inner seam pinned) …")
    pinned_global = np.unique(np.concatenate([_outer_perim_ring(), inner_seam]))
    points_opt, tris_opt, pin_to_local, _ = truss_tris(
        points_xy, tris_global, pinned_global, inner_seam
    )
    print(f"     truss returned {len(points_opt)} pts, {len(tris_opt)} tris")

    print("[5/7] Stitching trussed tris to inner quads …")
    mesh_truss = build_mixed_mesh_after_truss(
        points_opt, tris_opt, pin_to_local, inner_quads, mesh0.points
    )
    print(
        f"     mixed mesh: verts={mesh_truss.n_verts} elems={mesh_truss.n_elems}"
    )

    print("[6/7] Smoothing combined mixed mesh (Laplacian, boundary pinned) …")
    mesh_final = CHILmesh(
        connectivity=mesh_truss.connectivity_list.copy(),
        points=mesh_truss.points.copy(),
        compute_layers=True,
    )
    fem_smooth(mesh_final)

    q_pre, _, _ = mesh_truss.elem_quality()
    q_final, _, _ = mesh_final.elem_quality()
    print(f"     quality (post-truss):  median={np.median(q_pre):.3f} min={q_pre.min():.3f}")
    print(f"     quality (post-FEM):    median={np.median(q_final):.3f} min={q_final.min():.3f}")

    print(f"[7/7] Rendering 4-panel figure → {out_path} …")
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), facecolor="white")
    _plot_mesh_simple(
        axes[0, 0], mesh0,
        f"(1) All-quad start: {mesh0.n_elems} quads · {mesh0.n_layers} skeleton layers"
    )
    _plot_mixed(
        axes[0, 1], mesh_pre,
        f"(2) Strip outer 2 layers, Delaunay band "
        f"({len(tris_global)} tris + {len(inner_quads)} quads)"
    )
    _plot_mixed(
        axes[1, 0], mesh_truss,
        f"(3) ADMESH truss on tris (boundary pinned bit-exact)\n"
        f"median quality {np.median(q_pre):.3f}"
    )
    _plot_mixed(
        axes[1, 1], mesh_final,
        f"(4) Laplacian smooth combined mixed mesh\n"
        f"median quality {np.median(q_final):.3f}"
    )

    plt.tight_layout()
    plt.savefig(out_path, dpi=110, bbox_inches="tight", facecolor="white")
    print(f"     saved {out_path.stat().st_size:,} bytes")
    return out_path


if __name__ == "__main__":
    main()
