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
DEBUG_LOG = REPO_ROOT / "output" / "animation_debug.log"

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
    """Signed distance function matching CHILmesh annulus fixture (r_inner=0.3, r_outer=1.0)."""
    r = np.linalg.norm(points, axis=1)
    r_inner = 0.3
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
    """Build the 4 pipeline stages: raw Delaunay → ADMESH truss → FEM → layers.

    Uses examples.annulus() as Stage 1 — it is itself a Delaunay triangulation.
    scipy.spatial.Delaunay cannot represent non-convex (donut) domains without
    a constrained triangulation library, so we use the fixture directly and label
    it accurately: "Delaunay triangulation · raw input".
    """

    # Stage 1: Raw annulus Delaunay triangulation
    print("Stage 1: Annulus Delaunay triangulation (raw input)...")
    annulus_base = examples.annulus()
    raw_pts = annulus_base.points[:, :2].copy()
    raw_elems = np.array([list(c[:3]) for c in annulus_base.connectivity_list], dtype=np.int32)
    boundary_indices = _infer_boundary_from_delaunay(raw_pts, raw_elems)

    # Stage 2: ADMESH truss optimization to convergence
    print("Stage 2: ADMESH truss warm-start (spring relaxation)...")
    try:
        truss_pts, truss_elems = optimize_with_admesh_truss_arrays(
            raw_pts,
            raw_elems,
            sdf=_annulus_sdf,
            size_fn=None,
            boundary_indices=boundary_indices,
            niter=300,
            enforce_non_degradation=False,  # Always return truss output for visualization
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
                "name": "Delaunay triangulation",
                "algo": "Raw input mesh (unsmoothed)",
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
                 prev_pts=None, current_pts=None, debug=False, axis_bounds=None):
    """Render a single keyframe with interpolated mesh morphing and vertex dots.

    prev_pts: interpolated points for mesh rendering (morphing between stages).
    current_pts: positions for yellow vertex tracking dots.
    debug: if True, use fixed color instead of quality colormapping (for debugging morphing).
    axis_bounds: (x_min, x_max, y_min, y_max) tuple to use fixed axes (for smooth transitions).
    """
    ax_mesh.clear()
    ax_hist.clear()

    # Use interpolated mesh points if provided, otherwise use stage points
    pts = prev_pts if prev_pts is not None else stage["pts"]
    elems = stage["elems"]
    polys = [pts[elem] for elem in elems]

    # Debug: verify pts are actually different during transitions
    if debug and prev_pts is not None:
        pt_diff = np.linalg.norm(prev_pts - stage["pts"])
        with open(DEBUG_LOG, 'a') as f:
            f.write(f"DEBUG: point displacement magnitude = {pt_diff:.6f}\n")

    if stage["viz"] == "quality":
        # Color by element quality (cool_r colormap).
        # During debugging, use fixed color to isolate morphing from recoloring
        if debug:
            q = np.ones(len(elems))  # Constant color
            with open(DEBUG_LOG, 'a') as f:
                f.write(f"DEBUG: using constant color for morphing isolation\n")
        else:
            q = _quality_for(pts, elems)
        norm = Normalize(vmin=0.0, vmax=1.0)
        colors = matplotlib.colormaps[QCMAP](norm(q))
        pc = PolyCollection(polys, facecolors=colors, edgecolors="#1a1a1f", linewidths=0.5)
        ax_mesh.add_collection(pc)
        # Draw mesh edges explicitly for clarity
        edges = np.vstack([
            elems[:, [0, 1]],
            elems[:, [1, 2]],
            elems[:, [2, 0]],
        ])
        edge_segments = pts[edges]
        edge_coll = LineCollection(edge_segments, colors=EDGE, linewidths=0.3, alpha=0.4)
        ax_mesh.add_collection(edge_coll)
        cbar_title = "Element quality (cool→good)"
    else:
        # Color by layer index (viridis).
        elem_layer = stage["elem_layer"]
        n_layers = stage["n_layers"]
        norm = Normalize(vmin=0, vmax=max(1, n_layers - 1))
        colors = matplotlib.colormaps["viridis"](norm(elem_layer))
        pc = PolyCollection(polys, facecolors=colors, edgecolors="#1a1a1f", linewidths=0.5)
        ax_mesh.add_collection(pc)
        cbar_title = f"Layer index 0–{n_layers - 1}"

    # Draw mesh edges explicitly for clarity
    edges = np.vstack([
        elems[:, [0, 1]],
        elems[:, [1, 2]],
        elems[:, [2, 0]],
    ])
    edge_segments = pts[edges]
    edge_coll = LineCollection(edge_segments, colors=EDGE, linewidths=0.3, alpha=0.4)
    ax_mesh.add_collection(edge_coll)

    # Mesh axes setup: use fixed bounds if provided (smooth transitions), else auto
    if axis_bounds is not None:
        x_min, x_max, y_min, y_max = axis_bounds
    else:
        x_min, x_max = pts[:, 0].min(), pts[:, 0].max()
        y_min, y_max = pts[:, 1].min(), pts[:, 1].max()
        pad = 0.05 * max(x_max - x_min, y_max - y_min)
        x_min -= pad
        x_max += pad
        y_min -= pad
        y_max += pad
    ax_mesh.set_xlim(x_min, x_max)
    ax_mesh.set_ylim(y_min, y_max)
    ax_mesh.set_aspect("equal")
    ax_mesh.set_facecolor(BG)
    ax_mesh.set_xticks([])
    ax_mesh.set_yticks([])
    for spine in ax_mesh.spines.values():
        spine.set_color(DIM)

    # Add vertex tracking dots (yellow) if provided
    if current_pts is not None:
        for v_idx, (x, y) in enumerate(current_pts):
            circle = Circle((x, y), radius=0.012, color=VERTEX_DOT, alpha=0.7, zorder=10)
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


def _ease(t: float) -> float:
    """Smooth step easing: t in [0,1] → smoothed t."""
    return t * t * (3 - 2 * t)


def main():
    data = _stage_data()
    stages = data["stages"]
    n_stages = len(stages)
    qualities = [data["raw_quality"], data["truss_quality"],
                 data["fem_quality"], data["fem_quality"]]

    # Frame schedule: HOLD at each stage, TRANSITION between stages
    HOLD = 25     # frames to hold each stage (2.5s @ 10fps)
    TRANS = 12    # transition frames between stages (1.2s)

    # Build a flat list of (stage_idx, interp_t, pts_for_dots, quality_for_hist)
    # where interp_t=None means static, or 0.0..1.0 during transition
    frame_schedule = []
    for si in range(n_stages):
        # Extra hold time for final stage (skeletonization view)
        hold_frames = HOLD + (15 if si == n_stages - 1 else 0)
        for _ in range(hold_frames):
            frame_schedule.append({"type": "hold", "stage": si})
        if si < n_stages - 1:
            pts_a = stages[si]["pts"]
            pts_b = stages[si + 1]["pts"]
            can_interp = (len(pts_a) == len(pts_b))
            for fi in range(TRANS):
                t = _ease(fi / (TRANS - 1))
                frame_schedule.append({
                    "type": "trans",
                    "stage": si,          # base stage for mesh/label
                    "stage_to": si + 1,   # target stage
                    "t": t,
                    "can_interp": can_interp,
                })

    total_frames = len(frame_schedule)

    fig = plt.figure(figsize=(12, 5.5), facecolor=BG)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.18,
                          left=0.04, right=0.97, top=0.88, bottom=0.12)
    ax_mesh = fig.add_subplot(gs[0, 0])
    ax_hist = fig.add_subplot(gs[0, 1])

    fig.suptitle("CHILmesh pipeline · annulus", color=TEXT, fontsize=14,
                 fontweight="bold", y=0.98)

    def animate(frame_idx):
        info = frame_schedule[frame_idx]
        si = info["stage"]
        stage = stages[si]

        # Debug: track frame types
        if frame_idx < 5 or frame_idx % 20 == 0:
            print(f"Frame {frame_idx}: {info['type']} stage={si}", flush=True)

        if info["type"] == "hold":
            render_frame(ax_mesh, ax_hist, stage, qualities[si], si, n_stages,
                         current_pts=stage["pts"], debug=False)
        else:
            # Transition: morph mesh + dots between stages
            t = info["t"]
            si_to = info["stage_to"]
            if info["can_interp"]:
                # Interpolate both mesh points and quality
                interp_pts = (1 - t) * stages[si]["pts"] + t * stages[si_to]["pts"]
                interp_qual = (1 - t) * qualities[si] + t * qualities[si_to]
            else:
                # Discrete switch at midpoint
                if t > 0.5:
                    interp_pts = stages[si_to]["pts"]
                    interp_qual = qualities[si_to]
                else:
                    interp_pts = stages[si]["pts"]
                    interp_qual = qualities[si]

            # Compute fixed axis bounds that encompass both start and end stages
            # This prevents axes from rescaling during transitions, making morphing visible
            pts_a = stages[si]["pts"]
            pts_b = stages[si_to]["pts"]
            x_min = min(pts_a[:, 0].min(), pts_b[:, 0].min())
            x_max = max(pts_a[:, 0].max(), pts_b[:, 0].max())
            y_min = min(pts_a[:, 1].min(), pts_b[:, 1].min())
            y_max = max(pts_a[:, 1].max(), pts_b[:, 1].max())
            pad = 0.05 * max(x_max - x_min, y_max - y_min)
            axis_bounds = (x_min - pad, x_max + pad, y_min - pad, y_max + pad)

            render_frame(ax_mesh, ax_hist, stage, interp_qual, si, n_stages,
                         prev_pts=interp_pts, current_pts=interp_pts, debug=False,
                         axis_bounds=axis_bounds)

        # Stage label overlaid above mesh axes (in figure coords above ax_mesh)
        ax_mesh.text(
            0.02, 1.08, f"Stage {si + 1}/{n_stages}: {stage['name']}",
            transform=ax_mesh.transAxes, color=ACCENT, fontsize=13,
            fontweight="bold", ha="left", va="bottom",
        )
        ax_mesh.text(
            0.02, 1.02, stage["algo"],
            transform=ax_mesh.transAxes, color=TEXT, fontsize=10,
            ha="left", va="bottom",
        )
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
