#!/usr/bin/env python3
"""Mixed-element mesh demo: ADMESH ring + Delaunay gap + quad core.

Pipeline:
1. Generate 16x12 structured quad mesh (192 quads, 6 skeleton layers).
2. Skeletonize → split: ADMESH ring (layers 0-1), dropped band (layer 2), quad core (layers 3+).
3. ADMESH full distmesh on the outer ring: grid-sample initial interior, run truss loop.
4. Delaunay-triangulate the gap band (layer-2 region) from ring boundary nodes only.
5. Stitch ADMESH tris + gap tris + quad core into combined mixed-element mesh.
6. FEM smooth combined mesh (symmetric quad stiffness, boundary pinned).
7. Render 4-panel figure to output/mixed_truss_fem_demo.png.
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
from chilmesh._vendor_admesh_truss import distmesh2d_warmstart
from chilmesh.admesh_warmstart import distmesh1d


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
        [[vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)]
         for j in range(NY) for i in range(NX)],
        dtype=int,
    )
    return CHILmesh(connectivity=quads, points=points, compute_layers=True)


# ---------------------------------------------------------------------------
# Boundary helpers
# ---------------------------------------------------------------------------

def _outer_perim_ring() -> np.ndarray:
    """CCW outer perimeter vertex IDs of the original rectangle."""
    def vid(i, j): return j * (NX + 1) + i
    bottom = [vid(i, 0) for i in range(NX + 1)]
    right  = [vid(NX, j) for j in range(1, NY + 1)]
    top    = [vid(i, NY) for i in range(NX - 1, -1, -1)]
    left   = [vid(0, j) for j in range(NY - 1, 0, -1)]
    return np.array(bottom + right + top + left, dtype=int)


def _seam_ring(quads: np.ndarray) -> np.ndarray:
    """Walk CCW outer boundary of a rectangular quad sub-mesh."""
    edge_count: dict[tuple[int, int], int] = {}
    edge_dir: dict[tuple[int, int], tuple[int, int]] = {}
    for q in quads:
        for a, b in zip(q, np.roll(q, -1)):
            key = (min(a, b), max(a, b))
            edge_count[key] = edge_count.get(key, 0) + 1
            edge_dir.setdefault(key, (a, b))
    boundary = [edge_dir[k] for k, c in edge_count.items() if c == 1]
    nxt = {a: b for a, b in boundary}
    # pick start node NOT on outer perimeter if possible, else any
    outer_set = set(_outer_perim_ring().tolist())
    candidates = [a for a, _ in boundary if a not in outer_set]
    start = candidates[0] if candidates else boundary[0][0]
    walk = [start]
    cur = nxt[start]
    while cur != start:
        walk.append(cur)
        cur = nxt[cur]
    return np.array(walk, dtype=int)


def redistribute_outer_perim(points_xy: np.ndarray) -> np.ndarray:
    """Redistribute outer rectangle perim nodes with distmesh1d, corner-dense.

    h(p) = h_min + (h_max - h_min) * (1 - exp(-(d/sigma)^2))
    where d = distance from nearest corner. Smaller h near corners → tighter spacing.
    Returns a new points_xy array (copy) with outer perim positions updated.
    Inner verts unchanged.
    """
    perim_idx = _outer_perim_ring()
    perim_xy = points_xy[perim_idx]
    corners = np.array([[0.0, 0.0], [LX, 0.0], [LX, LY], [0.0, LY]])

    def h_fn(xy: np.ndarray) -> np.ndarray:
        d = np.min(np.linalg.norm(xy[:, None, :] - corners[None, :, :], axis=2), axis=1)
        return 0.05 + 0.45 * (1.0 - np.exp(-((d / 0.5) ** 2)))

    new_perim = distmesh1d(perim_xy, h_fn, closed=True, niter=2000, deltat=0.4, dptol=1e-6)
    out = points_xy.copy()
    out[perim_idx] = new_perim
    return out


def _polygon_sdf(points: np.ndarray, poly_verts: np.ndarray) -> np.ndarray:
    """Signed distance to closed polygon: negative inside."""
    a = poly_verts
    b = np.roll(poly_verts, -1, axis=0)
    min_d2 = np.full(len(points), np.inf)
    for i in range(len(a)):
        seg_a, seg_b = a[i], b[i]
        ab = seg_b - seg_a
        ab_len2 = float(ab @ ab)
        if ab_len2 == 0.0:
            continue
        t = np.clip(((points - seg_a) @ ab) / ab_len2, 0.0, 1.0)
        proj = seg_a + np.outer(t, ab)
        d2 = ((points - proj) ** 2).sum(axis=1)
        min_d2 = np.minimum(min_d2, d2)
    d = np.sqrt(min_d2)
    inside = MplPath(poly_verts).contains_points(points)
    return np.where(inside, -d, d)


# ---------------------------------------------------------------------------
# Stage 2: split into three regions
# ---------------------------------------------------------------------------

def split_mesh(mesh: CHILmesh, n_admesh: int = 2, n_drop: int = 1):
    """Return (kept_quads, seam_1_2, seam_2_3).

    seam_1_2: outer boundary of (layers n_admesh+) = inner boundary of ADMESH ring
    seam_2_3: outer boundary of (layers n_admesh+n_drop+) = inner boundary of gap band
    """
    n = mesh.n_layers
    assert n >= n_admesh + n_drop + 1, f"need ≥{n_admesh+n_drop+1} layers; got {n}"

    drop_plus_kept_ids = np.concatenate(
        [mesh.layers["OE"][k] for k in range(n_admesh, n)]
    )
    kept_ids = np.concatenate(
        [mesh.layers["OE"][k] for k in range(n_admesh + n_drop, n)]
    )
    drop_plus_kept_quads = mesh.connectivity_list[drop_plus_kept_ids]
    kept_quads = mesh.connectivity_list[kept_ids]

    seam_1_2 = _seam_ring(drop_plus_kept_quads)  # outer boundary of layers 2+
    seam_2_3 = _seam_ring(kept_quads)              # outer boundary of layers 3+
    return kept_quads, seam_1_2, seam_2_3


# ---------------------------------------------------------------------------
# Stage 3: ADMESH full distmesh on outer ring
# ---------------------------------------------------------------------------

def run_admesh_ring(
    points_xy: np.ndarray,
    seam_1_2: np.ndarray,
    *,
    niter: int = 400,
    Fscale: float = 1.2,
    dptol: float = 1e-3,
):
    """Grid-sample interior points within ring SDF and run ADMESH truss.

    Returns (admesh_pts, admesh_tris, pinned_global) where admesh_pts[0:B]
    are bit-exact copies of points_xy[pinned_global[0:B]].
    """
    outer_perim = _outer_perim_ring()
    pinned_global = np.unique(np.concatenate([outer_perim, seam_1_2]))
    pfix_xy = points_xy[pinned_global]

    seam_xy = points_xy[seam_1_2]

    def sdf(p: np.ndarray) -> np.ndarray:
        p = np.asarray(p, dtype=float)
        dx = np.maximum(-p[:, 0], p[:, 0] - LX)
        dy = np.maximum(-p[:, 1], p[:, 1] - LY)
        sdf_rect = np.where(
            (dx > 0) | (dy > 0),
            np.sqrt(np.maximum(dx, 0) ** 2 + np.maximum(dy, 0) ** 2),
            np.maximum(dx, dy),
        )
        sdf_inner = _polygon_sdf(p, seam_xy)
        return np.maximum(sdf_rect, -sdf_inner)

    h0 = LX / NX  # 0.25 — matches the quad grid spacing
    bbox = (-0.05, -0.05, LX + 0.05, LY + 0.05)
    geps = 1e-3 * h0

    # Grid-sample initial interior points
    xs = np.arange(h0 * 0.5, LX, h0 * 0.75)
    ys = np.arange(h0 * 0.5, LY, h0 * 0.75)
    xx, yy = np.meshgrid(xs, ys)
    cands = np.column_stack([xx.ravel(), yy.ravel()])
    interior_mask = sdf(cands) < -geps
    interior_init = cands[interior_mask]
    print(f"     ring SDF: {interior_mask.sum()} initial interior pts from {len(cands)}-pt grid")

    if len(interior_init) < 3:
        raise RuntimeError("No interior points sampled in ring domain — check SDF")

    admesh_pts, admesh_tris = distmesh2d_warmstart(
        pfix_xy, interior_init, sdf, None, h0, bbox,
        niter=niter, Fscale=Fscale, dptol=dptol, deltat=0.2,
    )
    return admesh_pts, admesh_tris, pinned_global


# ---------------------------------------------------------------------------
# Stage 4: Delaunay gap fill (layer-2 band)
# ---------------------------------------------------------------------------

def triangulate_gap(
    points_xy: np.ndarray,
    seam_1_2: np.ndarray,
    seam_2_3: np.ndarray,
) -> np.ndarray:
    """Delaunay-triangulate the gap band between seam_1_2 and seam_2_3.

    Returns triangle connectivity in original global vertex IDs.
    """
    gap_verts_global = np.unique(np.concatenate([seam_1_2, seam_2_3]))
    sites = points_xy[gap_verts_global]

    tri = Delaunay(sites)
    simplices_global = gap_verts_global[tri.simplices]
    centroids = points_xy[simplices_global].mean(axis=1)

    outer_path = MplPath(points_xy[seam_1_2])
    inner_path = MplPath(points_xy[seam_2_3])
    keep = (
        outer_path.contains_points(centroids, radius=1e-9)
        & ~inner_path.contains_points(centroids, radius=-1e-9)
    )
    return simplices_global[keep]


# ---------------------------------------------------------------------------
# Stage 5: stitch into combined mesh
# ---------------------------------------------------------------------------

def build_combined_mesh(
    admesh_pts: np.ndarray,
    admesh_tris: np.ndarray,
    pinned_global: np.ndarray,
    gap_tris_global: np.ndarray,
    kept_quads_global: np.ndarray,
    points_orig_xy: np.ndarray,
) -> CHILmesh:
    """Combine ADMESH tris, gap Delaunay tris, and kept quads.

    Point layout:
      0 .. B-1        → pinned boundary nodes (admesh_pts[0:B], bit-exact)
      B .. B+I-1      → ADMESH interior nodes (admesh_pts[B:])
      B+I ..          → original mesh nodes used by gap/kept but not pinned
    """
    B = len(pinned_global)

    all_orig_used = np.unique(np.concatenate([
        gap_tris_global.flatten(),
        kept_quads_global.flatten(),
    ]))
    orig_extra = np.setdiff1d(all_orig_used, pinned_global)

    pts_combined = np.vstack([admesh_pts, points_orig_xy[orig_extra]])
    N_admesh = len(admesh_pts)

    # Global original index → new combined index
    g2new: dict[int, int] = {}
    for i, g in enumerate(pinned_global):
        g2new[int(g)] = i
    for i, g in enumerate(orig_extra):
        g2new[int(g)] = N_admesh + i

    new_gap  = np.array([[g2new[v] for v in t] for t in gap_tris_global], dtype=int)
    new_kept = np.array([[g2new[v] for v in q] for q in kept_quads_global], dtype=int)

    # Pad tris to 4-col degenerate convention (first vert repeated)
    admesh_pad = np.column_stack([admesh_tris, admesh_tris[:, 0]])
    gap_pad    = np.column_stack([new_gap, new_gap[:, 0]])

    mixed_conn = np.vstack([admesh_pad, gap_pad, new_kept])
    pts_xyz = np.column_stack([pts_combined, np.zeros(len(pts_combined))])
    return CHILmesh(connectivity=mixed_conn, points=pts_xyz, compute_layers=False)


# ---------------------------------------------------------------------------
# Stage 6: Angle-based smooth (outer boundary pinned)
# ---------------------------------------------------------------------------

def fem_smooth(mesh: CHILmesh) -> CHILmesh:
    """FEM smoother with symmetric quad stiffness (Zhou & Shimada, boundary pinned)."""
    mesh.smooth_mesh(method='fem', acknowledge_change=True)
    return mesh


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _plot_mesh_simple(ax, mesh: CHILmesh, title: str):
    mesh.plot(ax=ax, elem_color="#cce0f3")
    ax.set_title(title, fontsize=10)
    ax.set_aspect("equal")
    ax.set_xlim(-0.15, LX + 0.15)
    ax.set_ylim(-0.15, LY + 0.15)
    ax.axis("off")


def _plot_mixed(ax, mesh: CHILmesh, title: str):
    conn = mesh.connectivity_list
    is_tri = (conn[:, 0] == conn[:, 3]) if conn.shape[1] == 4 else np.ones(len(conn), bool)
    tri_ids  = np.where(is_tri)[0]
    quad_ids = np.where(~is_tri)[0]

    if len(quad_ids):
        CHILmesh(connectivity=conn[quad_ids], points=mesh.points,
                 compute_layers=False).plot(ax=ax, elem_color="#ffd180")
    if len(tri_ids):
        CHILmesh(connectivity=conn[tri_ids, :3], points=mesh.points,
                 compute_layers=False).plot(ax=ax, elem_color="#a5d6a7")

    ax.set_title(title, fontsize=10)
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

    print("[1/6] Building 16x12 quad grid …")
    mesh0 = build_quad_grid()
    print(f"     verts={mesh0.n_verts} elems={mesh0.n_elems} layers={mesh0.n_layers}")

    print("[2/6] Splitting: ADMESH ring (L0+L1), drop (L2), quad core (L3+) …")
    kept_quads, seam_1_2, seam_2_3 = split_mesh(mesh0, n_admesh=2, n_drop=1)
    points_xy_uniform = mesh0.points[:, :2]
    print(f"     kept quads={len(kept_quads)}  seam_1_2={len(seam_1_2)} nodes  seam_2_3={len(seam_2_3)} nodes")

    print("[2b/6] distmesh1d on outer perim (corner-dense h(p)) …")
    points_xy = redistribute_outer_perim(points_xy_uniform)
    perim_idx = _outer_perim_ring()
    edges_old = np.linalg.norm(np.diff(np.vstack([points_xy_uniform[perim_idx], points_xy_uniform[perim_idx][:1]]), axis=0), axis=1)
    edges_new = np.linalg.norm(np.diff(np.vstack([points_xy[perim_idx], points_xy[perim_idx][:1]]), axis=0), axis=1)
    print(f"     perim edges before: min={edges_old.min():.3f} max={edges_old.max():.3f}")
    print(f"     perim edges after:  min={edges_new.min():.3f} max={edges_new.max():.3f}  (corner-dense)")

    print("[3/6] ADMESH full distmesh on outer ring …")
    admesh_pts, admesh_tris, pinned_global = run_admesh_ring(points_xy, seam_1_2)
    print(f"     ADMESH → {len(admesh_pts)} pts, {len(admesh_tris)} tris")

    print("[4/6] Delaunay gap fill (layer-2 band) …")
    gap_tris = triangulate_gap(points_xy, seam_1_2, seam_2_3)
    print(f"     gap → {len(gap_tris)} tris")

    print("[5/6] Stitching combined mesh …")
    mesh_pre = build_combined_mesh(
        admesh_pts, admesh_tris, pinned_global,
        gap_tris, kept_quads, points_xy,
    )
    n_tris_total = len(admesh_tris) + len(gap_tris)
    print(f"     combined: verts={mesh_pre.n_verts} tris={n_tris_total} quads={len(kept_quads)}")

    # Rebuild with same points/conn for a fresh mesh object for smoothing
    mesh_smooth = CHILmesh(
        connectivity=mesh_pre.connectivity_list.copy(),
        points=mesh_pre.points.copy(),
        compute_layers=True,
    )
    print("[6/6] FEM smooth (symmetric quad stiffness, boundary pinned) …")
    fem_smooth(mesh_smooth)

    q_pre, _, _   = mesh_pre.elem_quality()
    q_post, _, _  = mesh_smooth.elem_quality()
    print(f"     pre-smooth:  median={np.median(q_pre):.3f}  min={q_pre.min():.3f}")
    print(f"     post-smooth: median={np.median(q_post):.3f}  min={q_post.min():.3f}")

    # --- 4-panel figure ---
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), facecolor="white")

    _plot_mesh_simple(
        axes[0, 0], mesh0,
        f"(1) All-quad start: {mesh0.n_elems} quads · {mesh0.n_layers} skeleton layers",
    )

    # Panel (2): ADMESH ring only (just show it on blank background to highlight grading)
    admesh_mesh = CHILmesh(
        connectivity=admesh_tris,
        points=np.column_stack([admesh_pts, np.zeros(len(admesh_pts))]),
        compute_layers=False,
    )
    admesh_mesh.plot(ax=axes[0, 1], elem_color="#a5d6a7")
    axes[0, 1].set_title(
        f"(2) ADMESH ring (L0+L1): {len(admesh_tris)} tris\n"
        f"graded from outer rect → seam 1/2 (both pinned)",
        fontsize=10,
    )
    axes[0, 1].set_aspect("equal")
    axes[0, 1].set_xlim(-0.15, LX + 0.15)
    axes[0, 1].set_ylim(-0.15, LY + 0.15)
    axes[0, 1].axis("off")

    _plot_mixed(
        axes[1, 0], mesh_pre,
        f"(3) Assembled pre-smooth: {n_tris_total} tris + {len(kept_quads)} quads\n"
        f"(ADMESH + Delaunay gap + quad core)  median q={np.median(q_pre):.3f}",
    )
    _plot_mixed(
        axes[1, 1], mesh_smooth,
        f"(4) FEM smooth (symmetric quad stiffness, boundary pinned)\n"
        f"median quality {np.median(q_post):.3f}",
    )

    plt.tight_layout()
    plt.savefig(out_path, dpi=110, bbox_inches="tight", facecolor="white")
    print(f"     saved {out_path.stat().st_size:,} bytes → {out_path}")

    # Single-panel version for README inline use (final angle-based smoothed mesh only)
    single_path = out_path.parent / "mixed_mesh_fem_final.png"
    fig2, ax2 = plt.subplots(1, 1, figsize=(8, 6), facecolor="white")
    _plot_mixed(
        ax2, mesh_smooth,
        f"Mixed-element mesh after FEM smoothing (symmetric quad stiffness)\n"
        f"{n_tris_total} triangles + {len(kept_quads)} quads · median quality {np.median(q_post):.3f}",
    )
    plt.tight_layout()
    plt.savefig(single_path, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"     saved {single_path.stat().st_size:,} bytes → {single_path}")

    # 3-panel showcase: wireframe | layers | quality
    showcase_path = out_path.parent / "mixed_mesh_showcase.png"
    fig3, axes3 = plt.subplots(1, 3, figsize=(18, 6), facecolor="white")
    _plot_mixed(
        axes3[0], mesh_smooth,
        f"Mixed mesh · {n_tris_total} tris + {len(kept_quads)} quads",
    )
    mesh_smooth.plot_layer(ax=axes3[1])
    axes3[1].set_title("Skeletonization layers", fontsize=10)
    axes3[1].set_aspect("equal")
    axes3[1].axis("off")
    mesh_smooth.plot_quality(ax=axes3[2])
    axes3[2].set_title(f"Element quality · median {np.median(q_post):.3f}", fontsize=10)
    axes3[2].set_aspect("equal")
    axes3[2].axis("off")
    plt.tight_layout()
    plt.savefig(showcase_path, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"     saved {showcase_path.stat().st_size:,} bytes → {showcase_path}")

    return out_path


if __name__ == "__main__":
    main()
