"""Generate CHILmesh README hero animation as GIF.

Sequence (5 stages, loops layers -> layers; #187 lexicon + #198 truss interp):
  1. Peel layers of the RAW mesh -- animated inward reveal (peel_layers())
  2. Crossfade to element quality of the raw Delaunay triangulation
  3. ADMESH truss -- node positions + quality INTERPOLATED to convergence (#198)
  4. FEM smoother (Balendran) -- positions + quality interpolated
  5. Peel layers of the SMOOTHED mesh -- animated reveal, long hold, loop

Each frame: left = mesh viz with vertex tracking dots, right = colormapped
quality histogram of the current geometry.

Output: docs/gallery/readme_pipeline_annulus.gif (replaces existing).
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
OUT_PATH = REPO_ROOT / "docs" / "gallery" / "readme_pipeline_annulus.gif"

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
    print("Stage 4: Peel layers (final mesh)...")
    try:
        if fem_mesh is None:
            layer_mesh = _arrays_to_mesh(fem_pts, fem_elems)
        else:
            layer_mesh = fem_mesh
    except Exception as e:
        print(f"  Layer mesh construction warn: {e}")
        layer_mesh = _arrays_to_mesh(truss_pts, truss_elems)
    elem_layer, n_layers = _elem_layer_for(layer_mesh)

    # Stage 0 (opening view): peel layers of the RAW mesh, so the loop opens
    # and closes on the layer decomposition (peel -> quality -> truss -> FEM -> peel).
    print("Stage 0: Peel layers (raw mesh)...")
    raw_layer_mesh = _arrays_to_mesh(raw_pts, raw_elems)
    raw_elem_layer, raw_n_layers = _elem_layer_for(raw_layer_mesh)

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
                "name": "Peel Layers",
                "algo": f"peel_layers(): {n_layers} concentric layers (OE+IE per ring)",
                "pts": fem_pts, "elems": fem_elems,
                "viz": "layers",
                "elem_layer": elem_layer,
                "n_layers": n_layers,
            },
        ],
        "raw_elem_layer": raw_elem_layer,
        "raw_n_layers": raw_n_layers,
        "raw_quality": _quality_for(raw_pts, raw_elems),
        "truss_quality": _quality_for(truss_pts, truss_elems),
        "fem_quality": _quality_for(fem_pts, fem_elems),
    }


def _elem_layer_for(mesh):
    """Per-element layer index array + layer count from a peeled mesh."""
    layers = mesh.Layers
    elem_layer = np.zeros(mesh.n_elems, dtype=np.int32)
    n_layers = mesh.n_layers
    for li in range(n_layers):
        for e in layers["OE"][li]:
            elem_layer[e] = li
        for e in layers["IE"][li]:
            elem_layer[e] = li
    return elem_layer, n_layers


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
                 prev_pts=None, current_pts=None, quality_array_override=None):
    """Render a single keyframe with interpolated mesh morphing and vertex dots.

    prev_pts: interpolated points for mesh rendering (morphing between stages).
    current_pts: positions for yellow vertex tracking dots.
    quality_array_override: if provided, use this quality array for mesh coloring instead of recomputing.
    """
    ax_mesh.clear()
    ax_hist.clear()

    # Use interpolated mesh points if provided, otherwise use stage points
    pts = prev_pts if prev_pts is not None else stage["pts"]
    elems = stage["elems"]
    polys = [pts[elem] for elem in elems]

    if stage["viz"] == "quality":
        # Color by element quality (cool_r colormap).
        if quality_array_override is not None:
            q = quality_array_override
        else:
            try:
                q = _quality_for(pts, elems)
            except Exception:
                # Fall back to uniform quality=1.0 if quality computation fails
                # (can happen with interpolated degenerate meshes)
                q = np.ones(len(elems))
        norm = Normalize(vmin=0.0, vmax=1.0)
        colors = matplotlib.colormaps[QCMAP](norm(q))
        pc = PolyCollection(polys, facecolors=colors, edgecolors="#1a1a1f", linewidths=0.5)
        ax_mesh.add_collection(pc)
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


DIM_FACE = np.array([0.16, 0.16, 0.18, 1.0])  # un-revealed layer elements


def _quality_colors(q: np.ndarray) -> np.ndarray:
    """RGBA facecolors for a quality array (cool_r, 0..1)."""
    return matplotlib.colormaps[QCMAP](Normalize(vmin=0.0, vmax=1.0)(q))


def _layer_colors(elem_layer: np.ndarray, n_layers: int, reveal: float = 1.0) -> np.ndarray:
    """RGBA facecolors for a peel-layer view with progressive inward reveal.

    reveal in [0, 1]: layers with index <= reveal * (n_layers - 1) get their
    viridis color; deeper (not-yet-peeled) layers render dimmed. reveal=1.0
    shows the full decomposition. This is the "peel_layers interpolation" —
    the animation peels inward one ring at a time.
    """
    norm = Normalize(vmin=0, vmax=max(1, n_layers - 1))
    colors = matplotlib.colormaps["viridis"](norm(elem_layer))
    if reveal < 1.0:
        cutoff = reveal * max(1, n_layers - 1)
        hidden = elem_layer > cutoff
        colors[hidden] = DIM_FACE
    return colors


def _blend(colors_a: np.ndarray, colors_b: np.ndarray, t: float) -> np.ndarray:
    """Linear RGBA crossfade between two facecolor arrays (same length)."""
    return (1.0 - t) * colors_a + t * colors_b


def render_view(ax_mesh, ax_hist, pts, elems, facecolors, quality_arr,
                title, subtitle, cbar_label, dots=True):
    """Render one frame: mesh panel (explicit facecolors) + quality histogram."""
    ax_mesh.clear()
    ax_hist.clear()

    polys = [pts[elem] for elem in elems]
    pc = PolyCollection(polys, facecolors=facecolors,
                        edgecolors="#1a1a1f", linewidths=0.5)
    ax_mesh.add_collection(pc)

    edges = np.vstack([elems[:, [0, 1]], elems[:, [1, 2]], elems[:, [2, 0]]])
    edge_coll = LineCollection(pts[edges], colors=EDGE, linewidths=0.3, alpha=0.4)
    ax_mesh.add_collection(edge_coll)

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

    if dots:
        for (x, y) in pts:
            ax_mesh.add_patch(Circle((x, y), radius=0.012, color=VERTEX_DOT,
                                     alpha=0.7, zorder=10))

    # Histogram panel: colormapped quality of the current geometry.
    bins = 40
    counts, hedges = np.histogram(quality_arr, bins=bins, range=(0.0, 1.0))
    widths = np.diff(hedges)
    midpoints = hedges[:-1] + widths / 2.0
    bar_colors = matplotlib.colormaps[QCMAP](Normalize(0.0, 1.0)(midpoints))
    ax_hist.bar(hedges[:-1], counts, width=widths, align="edge",
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

    ax_mesh.text(0.02, 1.08, title, transform=ax_mesh.transAxes, color=ACCENT,
                 fontsize=13, fontweight="bold", ha="left", va="bottom")
    ax_mesh.text(0.02, 1.02, subtitle, transform=ax_mesh.transAxes, color=TEXT,
                 fontsize=10, ha="left", va="bottom")


def main():
    data = _stage_data()
    stages = data["stages"]
    raw, truss, fem, final_layers = stages
    q_raw = data["raw_quality"]
    q_truss = data["truss_quality"]
    q_fem = data["fem_quality"]
    raw_elem_layer = data["raw_elem_layer"]
    raw_n_layers = data["raw_n_layers"]

    # Pre-computed static facecolors
    fc_raw_quality = _quality_colors(q_raw)
    fc_truss_quality = _quality_colors(q_truss)
    fc_fem_quality = _quality_colors(q_fem)
    fc_raw_layers_full = _layer_colors(raw_elem_layer, raw_n_layers, 1.0)
    fc_final_layers_full = _layer_colors(final_layers["elem_layer"],
                                         final_layers["n_layers"], 1.0)

    can_interp = len(raw["pts"]) == len(truss["pts"]) == len(fem["pts"])
    print(f"Positional interpolation available: {can_interp}")

    # ------------------------------------------------------------------
    # Frame schedule (10 fps). Sequence per #187/#198 + operator direction:
    #   A  peel_layers reveal on the RAW mesh        (the opening)
    #   B  crossfade layers -> quality (raw mesh)
    #   C  hold quality (raw)
    #   D  truss morph  raw -> truss   (positions + quality interpolated)
    #   E  hold quality (truss)
    #   F  FEM morph    truss -> fem   (positions + quality interpolated)
    #   G  hold quality (fem)
    #   H  crossfade quality -> layers (final mesh)
    #   I  peel_layers reveal on the FINAL mesh, long hold -> loops back to A
    # ------------------------------------------------------------------
    REVEAL_A, HOLD_A = 20, 14
    FADE_B = 10
    HOLD_C = 14
    MORPH_D = 24
    HOLD_E = 12
    MORPH_F = 18
    HOLD_G = 14
    FADE_H = 10
    REVEAL_I, HOLD_I = 20, 26

    schedule = []

    def add(n, fn):
        for i in range(n):
            schedule.append((fn, i, n))

    # A: peel reveal on raw mesh
    def frame_A(i, n):
        t = _ease(i / max(1, n - 1))
        fc = _layer_colors(raw_elem_layer, raw_n_layers, t)
        k = int(round(t * max(1, raw_n_layers - 1)))
        return (raw["pts"], raw["elems"], fc, q_raw,
                "Stage 1/5: Peel Layers",
                f"peel_layers() on the input mesh: revealing layer {k + 1}/{raw_n_layers}")

    def frame_A_hold(i, n):
        return (raw["pts"], raw["elems"], fc_raw_layers_full, q_raw,
                "Stage 1/5: Peel Layers",
                f"peel_layers() on the input mesh: {raw_n_layers} concentric layers")

    # B: crossfade layers -> quality on raw mesh
    def frame_B(i, n):
        t = _ease(i / max(1, n - 1))
        fc = _blend(fc_raw_layers_full, fc_raw_quality, t)
        return (raw["pts"], raw["elems"], fc, q_raw,
                "Stage 2/5: Element Quality",
                "elem_quality(): raw Delaunay triangulation")

    def frame_C(i, n):
        return (raw["pts"], raw["elems"], fc_raw_quality, q_raw,
                "Stage 2/5: Element Quality (input)",
                "elem_quality(): raw Delaunay triangulation")

    # D: truss morph raw -> truss
    def frame_D(i, n):
        t = _ease(i / max(1, n - 1))
        if can_interp:
            pts = (1 - t) * raw["pts"] + t * truss["pts"]
            q = (1 - t) * q_raw + t * q_truss
        else:
            pts = truss["pts"] if t > 0.5 else raw["pts"]
            q = q_truss if t > 0.5 else q_raw
        return (pts, truss["elems"], _quality_colors(q), q,
                "Stage 3/5: ADMESH Truss",
                "Spring relaxation — node positions interpolating to convergence")

    def frame_E(i, n):
        return (truss["pts"], truss["elems"], fc_truss_quality, q_truss,
                "Stage 3/5: ADMESH Truss",
                "Spring relaxation (converged)")

    # F: FEM morph truss -> fem
    def frame_F(i, n):
        t = _ease(i / max(1, n - 1))
        if can_interp:
            pts = (1 - t) * truss["pts"] + t * fem["pts"]
            q = (1 - t) * q_truss + t * q_fem
        else:
            pts = fem["pts"] if t > 0.5 else truss["pts"]
            q = q_fem if t > 0.5 else q_truss
        return (pts, fem["elems"], _quality_colors(q), q,
                "Stage 4/5: FEM Smoother",
                "Balendran direct method — interpolating smoothed positions")

    def frame_G(i, n):
        return (fem["pts"], fem["elems"], fc_fem_quality, q_fem,
                "Stage 4/5: FEM Smoother",
                "Balendran direct method (smoothed)")

    # H: crossfade quality -> layers on final mesh
    def frame_H(i, n):
        t = _ease(i / max(1, n - 1))
        fc = _blend(fc_fem_quality, fc_final_layers_full, t)
        return (fem["pts"], fem["elems"], fc, q_fem,
                "Stage 5/5: Peel Layers",
                "peel_layers() on the smoothed mesh")

    # I: peel reveal on final mesh (then long hold; GIF loops back to A)
    def frame_I(i, n):
        t = _ease(i / max(1, n - 1))
        nl = final_layers["n_layers"]
        fc = _layer_colors(final_layers["elem_layer"], nl, t)
        k = int(round(t * max(1, nl - 1)))
        return (fem["pts"], fem["elems"], fc, q_fem,
                "Stage 5/5: Peel Layers",
                f"peel_layers() on the smoothed mesh: revealing layer {k + 1}/{nl}")

    def frame_I_hold(i, n):
        nl = final_layers["n_layers"]
        return (fem["pts"], fem["elems"], fc_final_layers_full, q_fem,
                "Stage 5/5: Peel Layers",
                f"peel_layers() on the smoothed mesh: {nl} layers — loop restarts")

    add(REVEAL_A, frame_A)
    add(HOLD_A, frame_A_hold)
    add(FADE_B, frame_B)
    add(HOLD_C, frame_C)
    add(MORPH_D, frame_D)
    add(HOLD_E, frame_E)
    add(MORPH_F, frame_F)
    add(HOLD_G, frame_G)
    add(FADE_H, frame_H)
    add(REVEAL_I, frame_I)
    add(HOLD_I, frame_I_hold)

    total_frames = len(schedule)

    fig = plt.figure(figsize=(12, 5.5), facecolor=BG)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.18,
                          left=0.04, right=0.97, top=0.88, bottom=0.12)
    ax_mesh = fig.add_subplot(gs[0, 0])
    ax_hist = fig.add_subplot(gs[0, 1])
    fig.suptitle("CHILmesh pipeline · annulus", color=TEXT, fontsize=14,
                 fontweight="bold", y=0.98)

    def animate(frame_idx):
        fn, i, n = schedule[frame_idx]
        pts, elems, fc, q, title, subtitle = fn(i, n)
        render_view(ax_mesh, ax_hist, pts, np.asarray(elems), fc, q,
                    title, subtitle, "")
        return []

    print(f"Rendering {total_frames} frames → {OUT_PATH}...")
    anim = animation.FuncAnimation(fig, animate, frames=total_frames,
                                   interval=100, blit=False)
    writer = animation.PillowWriter(fps=10)
    anim.save(str(OUT_PATH), writer=writer, dpi=110)
    plt.close(fig)
    print(f"Done: {OUT_PATH} ({OUT_PATH.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
