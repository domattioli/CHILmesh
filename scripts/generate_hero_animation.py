"""Generate CHILmesh README hero animation as GIF.

Pipeline (4 stages, GIF loop):
  1. Raw annulus (Delaunay input)
  2. ADMESH truss warm-start (spring relaxation)
  3. FEM smoother (Balendran direct)
  4. Skeletonization layers (medial-axis viz of stage-3 mesh)

Each frame: left = mesh viz, right = colormapped quality histogram.

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
import numpy as np

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

QCMAP = "cool_r"


def _stage_data():
    """Build the 4 pipeline stages from the annulus fixture."""
    print("Loading annulus fixture...")
    raw = examples.annulus()
    raw_pts = raw.points[:, :2].copy()
    raw_elems = np.array([list(c[:3]) for c in raw.connectivity_list], dtype=np.int32)

    print("Stage 2: angle-based smoother (1 pass)...")
    ab_mesh = _arrays_to_mesh(raw_pts, raw_elems)
    try:
        ab_mesh.smooth_mesh(method="angle-based", acknowledge_change=True)
    except Exception as e:
        print(f"  Angle-based smoother warn: {e}")
    ab_pts = ab_mesh.points[:, :2].copy()
    ab_elems = np.array([list(c[:3]) for c in ab_mesh.connectivity_list], dtype=np.int32)

    print("Stage 3: FEM smoother...")
    fem_mesh = _arrays_to_mesh(ab_pts, ab_elems)
    try:
        fem_mesh.smooth_mesh(method="fem", acknowledge_change=True)
    except Exception as e:
        print(f"  FEM smoother warn: {e}")

    fem_pts = fem_mesh.points[:, :2].copy()
    fem_elems = np.array([list(c[:3]) for c in fem_mesh.connectivity_list], dtype=np.int32)

    print("Stage 4: skeletonization layers...")
    layer_mesh = _arrays_to_mesh(fem_pts, fem_elems)
    layers = layer_mesh.Layers
    # Build per-element layer-id mapping.
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
                "name": "Raw annulus",
                "algo": "Delaunay input",
                "pts": raw_pts, "elems": raw_elems,
                "viz": "quality",
            },
            {
                "name": "Angle-based smoother",
                "algo": "Zhou-Shimada angle maximisation",
                "pts": ab_pts, "elems": ab_elems,
                "viz": "quality",
            },
            {
                "name": "FEM smoother",
                "algo": "Balendran direct method",
                "pts": fem_pts, "elems": fem_elems,
                "viz": "quality",
            },
            {
                "name": "Skeletonization layers",
                "algo": f"Medial-axis: {n_layers} layers (OE+IE per ring)",
                "pts": fem_pts, "elems": fem_elems,
                "viz": "layers",
                "elem_layer": elem_layer,
                "n_layers": n_layers,
            },
        ],
        "fem_quality": _quality_for(fem_pts, fem_elems),
        "raw_quality": _quality_for(raw_pts, raw_elems),
        "ab_quality": _quality_for(ab_pts, ab_elems),
    }


def _arrays_to_mesh(pts, elems):
    """Build a CHILmesh from points + triangle connectivity.

    Mesh signature is (connectivity, points, ...) — be careful of arg order.
    """
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


def render_frame(ax_mesh, ax_hist, stage, quality_arr, stage_idx, n_stages):
    """Render a single keyframe."""
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

    # Stage name/algo overlay on mesh panel.
    ax_mesh.text(
        0.02, 0.98, f"Stage {stage_idx + 1}/{n_stages}",
        transform=ax_mesh.transAxes, color=DIM, fontsize=11,
        ha="left", va="top",
    )
    ax_mesh.text(
        0.02, 0.93, stage["name"],
        transform=ax_mesh.transAxes, color=ACCENT, fontsize=18,
        fontweight="bold", ha="left", va="top",
    )
    ax_mesh.text(
        0.02, 0.88, stage["algo"],
        transform=ax_mesh.transAxes, color=TEXT, fontsize=10,
        ha="left", va="top",
    )

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
    qualities = [data["raw_quality"], data["ab_quality"],
                 data["fem_quality"], data["fem_quality"]]

    fig = plt.figure(figsize=(12, 5.5), facecolor=BG)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.18,
                          left=0.04, right=0.97, top=0.95, bottom=0.12)
    ax_mesh = fig.add_subplot(gs[0, 0])
    ax_hist = fig.add_subplot(gs[0, 1])

    fig.suptitle("CHILmesh pipeline · annulus", color=TEXT, fontsize=14,
                 fontweight="bold", y=0.985)

    # Frames: hold each stage for ~20 frames at 10fps = 2s per stage.
    HOLD_FRAMES = 18
    total_frames = HOLD_FRAMES * n_stages

    def animate(frame_idx):
        stage_idx = frame_idx // HOLD_FRAMES
        render_frame(ax_mesh, ax_hist, stages[stage_idx],
                     qualities[stage_idx], stage_idx, n_stages)
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
