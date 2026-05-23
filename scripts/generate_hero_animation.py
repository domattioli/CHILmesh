"""Generate CHILmesh README hero animation as GIF.

Pipeline (4 stages, GIF loop):
  1. Random Delaunay (initial triangulation from random points)
  2. ADMESH truss (spring relaxation solved to convergence)
  3. FEM smoother (Balendran direct method)
  4. Skeletonization layers (medial-axis viz of stage-3 mesh)

Each frame: left = mesh viz with vertex tracking dots, right = colormapped quality histogram.

Output: output/readme_pipeline_annulus.gif (replaces existing).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Circle
import numpy as np
from scipy.spatial import Delaunay

from chilmesh import Mesh, examples
from chilmesh import optimize_with_admesh_truss_arrays

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_PATH = REPO_ROOT / "output" / "readme_pipeline_annulus.gif"

BG = "#0e0e10"
EDGE = "#3a7fbf"
FILL = "#5fb0ff"
TEXT = "#e8e8ee"
DIM = "#888888"
ACCENT = "#ff9f43"
GOOD = "#2ecc71"
VERTEX_DOT = "#ffff00"  # Yellow for vertex tracking

QCMAP = "cool_r"


def _annulus_sdf(points: np.ndarray) -> np.ndarray:
    """Signed distance function for annulus domain (r_inner=0.4, r_outer=1.0)."""
    r = np.linalg.norm(points, axis=1)
    r_inner = 0.4
    r_outer = 1.0
    # Negative inside annulus, positive outside
    return np.maximum(r - r_outer, r_inner - r)


def _annulus_h0(size: float = 0.15) -> float:
    """Target edge length for annulus domain."""
    return size


def _project_to_annulus_boundary(points: np.ndarray, sdf_func) -> np.ndarray:
    """Project points to annulus boundary via binary search on the radial direction."""
    pts_proj = points.copy()
    for i, pt in enumerate(points):
        r = np.linalg.norm(pt)
        if r < 1e-6:
            # Point at origin; move to outer boundary
            pts_proj[i] = np.array([1.0, 0.0])
        else:
            # Binary search to find radius where sdf = 0 along the radial direction
            direction = pt / r
            r_min, r_max = 0.3, 1.1
            for _ in range(20):  # 20 iterations ≈ 1e-6 precision
                r_mid = (r_min + r_max) / 2.0
                sdf_val = sdf_func(np.array([[direction[0] * r_mid, direction[1] * r_mid]]))[0]
                if sdf_val < 0:
                    r_min = r_mid
                else:
                    r_max = r_mid
            r_final = (r_min + r_max) / 2.0
            pts_proj[i] = direction * r_final
    return pts_proj


def _infer_boundary_from_delaunay(points: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Identify boundary vertices from Delaunay triangulation (convex hull)."""
    edges = np.vstack([
        triangles[:, [0, 1]],
        triangles[:, [1, 2]],
        triangles[:, [2, 0]],
    ])
    edges_sorted = np.sort(edges, axis=1)
    edge_count = {}
    for e in edges_sorted:
        key = tuple(e)
        edge_count[key] = edge_count.get(key, 0) + 1
    boundary_vertices = set()
    for (i, j), count in edge_count.items():
        if count == 1:
            boundary_vertices.add(i)
            boundary_vertices.add(j)
    return np.array(sorted(boundary_vertices), dtype=np.int32)


def _stage_data():
    """Build the 4 pipeline stages from random Delaunay + ADMESH truss + FEM + layers."""

    # Stage 1: Random Delaunay in annulus domain
    print("Stage 1: Random Delaunay initialization...")

    # Use annulus fixture as base to ensure proper boundary constraints
    annulus_base = examples.annulus()
    raw_pts_base = annulus_base.points[:, :2].copy()
    raw_elems_base = np.array([list(c[:3]) for c in annulus_base.connectivity_list], dtype=np.int32)

    # Add small perturbations to interior points (not boundary) to make it "random"
    np.random.seed(42)
    boundary_mask = np.zeros(len(raw_pts_base), dtype=bool)
    boundary_indices_base = _infer_boundary_from_delaunay(raw_pts_base, raw_elems_base)
    boundary_mask[boundary_indices_base] = True

    # Perturb interior points slightly
    perturbation = np.random.randn(len(raw_pts_base), 2) * 0.02
    perturbation[boundary_mask] = 0  # Keep boundary fixed
    raw_pts = raw_pts_base + perturbation
    raw_elems = raw_elems_base  # Keep same connectivity

    # Ensure points are valid (sometimes perturbation can push interior points outside domain)
    raw_pts = _project_to_annulus_boundary(raw_pts, _annulus_sdf)

    # Recompute Delaunay to get proper triangulation for perturbed points
    delaunay = Delaunay(raw_pts)
    raw_elems = np.asarray(delaunay.simplices, dtype=np.int32)

    # Identify boundary vertices
    boundary_indices = _infer_boundary_from_delaunay(raw_pts, raw_elems)

    # Stage 2: ADMESH truss optimization to convergence
    print("Stage 2: ADMESH truss warm-start (spring relaxation)...")
    try:
        truss_pts, truss_elems = optimize_with_admesh_truss_arrays(
            raw_pts,
            raw_elems,
            sdf=_annulus_sdf,
            size_fn=None,  # Uniform size
            h0=_annulus_h0(),
            boundary_indices=boundary_indices,
            niter=500,  # Solve to convergence
            enforce_non_degradation=True,
        )
    except Exception as e:
        print(f"  ADMESH truss warn: {e}")
        truss_pts, truss_elems = raw_pts, raw_elems

    truss_elems = np.asarray(truss_elems, dtype=np.int32)

    # Stage 3: FEM smoother
    print("Stage 3: FEM smoother...")
    fem_mesh = None
    try:
        fem_mesh = _arrays_to_mesh(truss_pts, truss_elems)
        fem_mesh.smooth_mesh(method="fem", acknowledge_change=True)
        fem_pts = fem_mesh.points[:, :2].copy()
        fem_elems = np.array([list(c[:3]) for c in fem_mesh.connectivity_list], dtype=np.int32)
        # Check for NaN/Inf
        if not np.isfinite(fem_pts).all():
            raise ValueError("FEM smoother produced non-finite coordinates")
    except Exception as e:
        print(f"  FEM smoother warn: {e}; using truss result as fallback")
        fem_pts = truss_pts.copy()
        fem_elems = truss_elems.copy()
        fem_mesh = None

    # Stage 4: Skeletonization layers
    print("Stage 4: Skeletonization layers...")
    try:
        if fem_mesh is None:
            layer_mesh = _arrays_to_mesh(fem_pts, fem_elems)
        else:
            layer_mesh = fem_mesh
    except Exception as e:
        print(f"  Layer mesh construction warn: {e}")
        layer_mesh = _arrays_to_mesh(truss_pts, truss_elems)
    layers = layer_mesh.Layers
    elem_layer = np.zeros(layer_mesh.n_elems, dtype=np.int32)
    n_layers = layer_mesh.n_layers
    for li in range(n_layers):
        oe = layers["OE"][li]
        ie = layers["IE"][li]
        for e in oe:
            elem_layer[e] = li
        for e in ie:
            elem_layer[e] = li

    return {
        "stages": [
            {
                "name": "Random Delaunay",
                "algo": "Perturbed annulus mesh + Delaunay",
                "pts": raw_pts, "elems": raw_elems,
                "viz": "quality",
            },
            {
                "name": "ADMESH Truss",
                "algo": "Spring relaxation (convergence)",
                "pts": truss_pts, "elems": truss_elems,
                "viz": "quality",
            },
            {
                "name": "FEM Smoother",
                "algo": "Balendran direct method",
                "pts": fem_pts, "elems": fem_elems,
                "viz": "quality",
            },
            {
                "name": "Skeletonization Layers",
                "algo": f"Medial-axis: {n_layers} layers (OE+IE per ring)",
                "pts": fem_pts, "elems": fem_elems,
                "viz": "layers",
                "elem_layer": elem_layer,
                "n_layers": n_layers,
            },
        ],
        "raw_quality": _quality_for(raw_pts, raw_elems),
        "truss_quality": _quality_for(truss_pts, truss_elems),
        "fem_quality": _quality_for(fem_pts, fem_elems),
    }


def _arrays_to_mesh(pts, elems):
    """Build a CHILmesh from points + triangle connectivity."""
    points3d = np.column_stack([pts, np.zeros(len(pts))]).astype(np.float64)
    elems_i = np.asarray(elems, dtype=np.int64)
    conn = np.column_stack([elems_i[:, 0], elems_i[:, 1], elems_i[:, 2], elems_i[:, 0]])
    return Mesh(connectivity=conn, points=points3d,
                compute_layers=True, compute_adjacencies=True)


def _quality_for(pts, elems):
    """Element quality array via temporary mesh."""
    m = _arrays_to_mesh(pts, elems)
    q, _, _ = m.elem_quality()
    return np.asarray(q)


def render_frame(ax_mesh, ax_hist, stage, quality_arr, stage_idx, n_stages,
                 prev_pts=None, current_pts=None):
    """Render a single keyframe with optional vertex tracking dots."""
    ax_mesh.clear()
    ax_hist.clear()

    pts = stage["pts"]
    elems = stage["elems"]
    polys = [pts[elem] for elem in elems]

    if stage["viz"] == "quality":
        # Color by element quality (cool_r colormap).
        q = _quality_for(pts, elems)
        norm = Normalize(vmin=0.0, vmax=1.0)
        colors = matplotlib.colormaps[QCMAP](norm(q))
        pc = PolyCollection(polys, facecolors=colors, edgecolors="#1a1a1f", linewidths=0.4)
        ax_mesh.add_collection(pc)
        cbar_title = "Element quality (cool→good)"
    else:
        # Color by layer index (viridis).
        elem_layer = stage["elem_layer"]
        n_layers = stage["n_layers"]
        norm = Normalize(vmin=0, vmax=max(1, n_layers - 1))
        colors = matplotlib.colormaps["viridis"](norm(elem_layer))
        pc = PolyCollection(polys, facecolors=colors, edgecolors="#1a1a1f", linewidths=0.4)
        ax_mesh.add_collection(pc)
        cbar_title = f"Layer index 0–{n_layers - 1}"

    # Mesh axes setup.
    x_min, x_max = pts[:, 0].min(), pts[:, 0].max()
    y_min, y_max = pts[:, 1].min(), pts[:, 1].max()
    pad = 0.05 * max(x_max - x_min, y_max - y_min)
    ax_mesh.set_xlim(x_min - pad, x_max + pad)
    ax_mesh.set_ylim(y_min - pad, y_max + pad)
    ax_mesh.set_aspect("equal")
    ax_mesh.set_facecolor(BG)
    ax_mesh.set_xticks([])
    ax_mesh.set_yticks([])
    for spine in ax_mesh.spines.values():
        spine.set_color(DIM)

    # Add vertex tracking dots (yellow) if we have current points
    if current_pts is not None:
        for v_idx, (x, y) in enumerate(current_pts):
            circle = Circle((x, y), radius=0.008, color=VERTEX_DOT, alpha=0.6, zorder=10)
            ax_mesh.add_patch(circle)

    # Histogram panel (always colormapped quality of CURRENT stage mesh).
    bins = 40
    counts, edges = np.histogram(quality_arr, bins=bins, range=(0.0, 1.0))
    widths = np.diff(edges)
    midpoints = edges[:-1] + widths / 2.0
    bar_colors = matplotlib.colormaps[QCMAP](Normalize(vmin=0.0, vmax=1.0)(midpoints))
    ax_hist.bar(edges[:-1], counts, width=widths, align="edge",
                color=bar_colors, edgecolor="#1a1a1f", linewidth=0.3)

    med_q = float(np.median(quality_arr))
    min_q = float(np.min(quality_arr))
    mean_q = float(np.mean(quality_arr))

    ax_hist.axvline(med_q, color=GOOD, linestyle="--", linewidth=1.5, alpha=0.85)
    ax_hist.set_xlim(0.0, 1.0)
    ax_hist.set_facecolor(BG)
    ax_hist.set_xlabel("Element quality", color=TEXT, fontsize=11)
    ax_hist.set_ylabel("Count", color=TEXT, fontsize=11)
    ax_hist.tick_params(colors=DIM, labelsize=9)
    for spine in ax_hist.spines.values():
        spine.set_color(DIM)
    ax_hist.set_title(
        f"Median: {med_q:.3f}    Min: {min_q:.3f}    Mean: {mean_q:.3f}",
        color=TEXT, fontsize=11, pad=8,
    )


def main():
    data = _stage_data()
    stages = data["stages"]
    n_stages = len(stages)
    qualities = [data["raw_quality"], data["truss_quality"],
                 data["fem_quality"], data["fem_quality"]]

    fig = plt.figure(figsize=(12, 5.5), facecolor=BG)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.18,
                          left=0.04, right=0.97, top=0.88, bottom=0.12)
    ax_mesh = fig.add_subplot(gs[0, 0])
    ax_hist = fig.add_subplot(gs[0, 1])

    # Place stage title/algo text above figure (outside axes) for readability
    fig.suptitle("CHILmesh pipeline · annulus", color=TEXT, fontsize=14,
                 fontweight="bold", y=0.98)

    # Frames: hold each stage for ~20 frames at 10fps = 2s per stage.
    HOLD_FRAMES = 18
    total_frames = HOLD_FRAMES * n_stages

    def animate(frame_idx):
        stage_idx = frame_idx // HOLD_FRAMES

        # Add stage annotation on the plot area itself (positioned at top-left, outside mesh area)
        ax_mesh.text(
            0.02, 1.08, f"Stage {stage_idx + 1}/{n_stages}: {stages[stage_idx]['name']}",
            transform=ax_mesh.transAxes, color=ACCENT, fontsize=13,
            fontweight="bold", ha="left", va="bottom",
        )
        ax_mesh.text(
            0.02, 1.03, stages[stage_idx]['algo'],
            transform=ax_mesh.transAxes, color=TEXT, fontsize=10,
            ha="left", va="bottom",
        )

        render_frame(ax_mesh, ax_hist, stages[stage_idx],
                     qualities[stage_idx], stage_idx, n_stages,
                     current_pts=stages[stage_idx]["pts"])
        return []

    print(f"Rendering {total_frames} frames → {OUT_PATH}...")
    anim = animation.FuncAnimation(
        fig, animate, frames=total_frames, interval=100, blit=False,
    )
    writer = animation.PillowWriter(fps=10)
    anim.save(str(OUT_PATH), writer=writer, dpi=110)
    plt.close(fig)
    print(f"Done: {OUT_PATH} ({OUT_PATH.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
