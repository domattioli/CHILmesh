"""Generate CHILmesh README hero animation as GIF.

Sequence (single loop; #179 + #187 + operator direction 2026-07-10, refined
per spec .specify/specs/002-hero-gif-refinement):
  0. Seed      — blank paper, boundary nodes appear, interior nodes FALL onto
                 the paper (scattered), wireframe draws. No colors. A held pause.
  1. Truss     — the scattered seed relaxes toward a quality mesh via
                 equal-motion quantile-resampled playback of real solver
                 iterations (uniform perceived pacing, not fixed-index
                 sampling); nodes AND faces move together, quality coloring +
                 histogram fade in.
  2. FEM       — Balendran direct smoother; the mesh settles further, faces
                 continue to deform, quality climbs.
  3. Peel      — the final peel_layers() concentric layers, revealed ONE
                 LAYER AT A TIME, boundary-inward, with the quality
                 histogram converting to matching stacked-bar segments in
                 lock-step (per-bin totals never change — only the color
                 breakdown does). Loop ends on a held full-layer view.

The truss stage plays back ACTUAL solver iterations: distmesh2d_warmstart is
run once with history_out= capturing (points, triangulation) every iteration
(deltat=0.02 — the converging regime for a cold scattered seed), with
snapshot_retriangulate=True and snapshot_strict_interior=True so every
captured triangulation is a fresh Delaunay of THAT iteration's own points.
This fixes a 1-Euler-step staleness in the raw (p, t) capture that otherwise
renders inverted "overlapping sail" triangles — the connectivity change is
diagnostic/render-only, physics and node positions are untouched (see
_vendor_admesh_truss.py docstring for the kwarg contract).

IMPORTANT — render vs. physics connectivity: distmesh2d_warmstart's RETURN
value (p_out, t_out) still carries the pre-move, one-step-stale connectivity
from its final Euler step, because snapshot_retriangulate only rewrites the
CAPTURED history entries, not the return. FEM smoothing and the converged-
truss frame therefore source BOTH points and connectivity from the solver's
LAST CAPTURED snapshot (hist[-1]), never from the raw return. Do not "fix"
this by reverting to the solver's return value in a future edit — that
reintroduces the staleness this feature removed.

Frames render at 10 fps via matplotlib's PillowWriter; the saved GIF is then
re-encoded through _optimize_gif() with a single global adaptive palette
(quantized, dither=None, optimize=True) to stay under the size budget with
no visual loss — repeated/held frames are free (PillowWriter already merges
runs of pixel-identical consecutive frames into one stored frame with
extended duration), so the budget is driven by unique-frame count only.

Left panel = mesh + vertex dots; right panel = live element-quality
histogram, with median/mean/min called out via _draw_metrics() (green
dotted line + text for Median, red dotted line + text for Mean, neutral
text for Min — see _draw_metrics()).
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
from PIL import Image, ImageSequence

from chilmesh import Mesh, examples

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_PATH = REPO_ROOT / "docs" / "gallery" / "readme_pipeline_annulus.gif"

BG = "#0e0e10"
EDGE = "#3a7fbf"
TEXT = "#e8e8ee"
DIM = "#888888"
ACCENT = "#ff9f43"
GOOD = "#2ecc71"       # Median (green line + text)
RED = "#ff5555"        # Mean (red line + text) — the one new color constant (F-08)
VERTEX_DOT = "#ffff00"      # interior vertex tracking dots
BOUNDARY_DOT = "#ff7a1a"    # boundary nodes (distinct)
WIRE = "#5a5a66"            # colorless wireframe
QCMAP = "cool_r"
R_IN, R_OUT = 0.3, 1.0
HBINS = 40
N_SNAP = 32             # equal-motion quantile snapshot count (US2)
TRUSS_HOLD = 2          # uniform per-snapshot hold, frames (US2)


def _quality(pts, elems):
    """Signed-area triangle quality q = 4√3·A / (l0²+l1²+l2²) in [0,1]."""
    p0, p1, p2 = pts[elems[:, 0]], pts[elems[:, 1]], pts[elems[:, 2]]
    area = 0.5 * np.abs((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
                        - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1]))
    l = (np.sum((p1 - p0) ** 2, 1) + np.sum((p2 - p1) ** 2, 1)
         + np.sum((p0 - p2) ** 2, 1))
    with np.errstate(divide="ignore", invalid="ignore"):
        q = 4.0 * np.sqrt(3.0) * area / l
    return np.clip(np.nan_to_num(q), 0.0, 1.0)


def _boundary_indices(elems, n):
    """Vertices on edges that belong to exactly one triangle."""
    E = np.sort(np.vstack([elems[:, [0, 1]], elems[:, [1, 2]], elems[:, [2, 0]]]), 1)
    uniq, cnt = np.unique(E, axis=0, return_counts=True)
    return np.unique(uniq[cnt == 1])


def _stage_data():
    """Every animation frame comes from REAL pipeline output (operator 2026-07-10):
    the truss stage plays back actual per-iteration solver snapshots captured via
    distmesh2d_warmstart(history_out=...); FEM + peel run on the solver's last
    captured snapshot (hist[-1], never the raw return — see module docstring).
    The annulus fixture supplies only the boundary node rings."""
    print("Boundary rings from the annulus fixture...")
    a = examples.annulus()
    P_fix = a.points[:, :2].copy()
    T_fix = np.array([list(c[:3]) for c in a.connectivity_list], dtype=np.int64)
    bidx_fix = _boundary_indices(T_fix, len(P_fix))
    bnd = P_fix[bidx_fix]
    h0 = float(np.median(np.linalg.norm(P_fix[T_fix[:, 0]] - P_fix[T_fix[:, 1]], axis=1)))

    print("Scattering random interior nodes (rejection sample in the band)...")
    rng = np.random.default_rng(11)
    n_int = len(P_fix) - len(bnd)
    pts = []
    while len(pts) < n_int:
        xy = rng.uniform(-1, 1, size=(500, 2))
        r = np.linalg.norm(xy, axis=1)
        for q in xy[(r < R_OUT - 0.04) & (r > R_IN + 0.04)]:
            pts.append(q)
            if len(pts) >= n_int:
                break
    interior_pts = np.array(pts)

    def fd(p):
        r = np.linalg.norm(p, axis=1)
        return np.maximum(r - R_OUT, R_IN - r)

    print("Running ADMESH truss with per-iteration history capture...")
    from chilmesh._vendor_admesh_truss import distmesh2d_warmstart
    hist = []
    distmesh2d_warmstart(
        bnd, interior_pts, fd, None, h0, (-1, -1, 1, 1),
        deltat=0.02, Fscale=1.2, niter=200,
        track_best_quality=False, history_out=hist, history_every=1,
        snapshot_retriangulate=True, snapshot_strict_interior=True,
    )
    # F-03: the solver's RETURN value is intentionally unused for rendering —
    # see module docstring. Everything downstream sources from hist[-1].
    p_fin, t_fin = hist[-1]
    p_fin = np.asarray(p_fin)[:, :2]
    t_fin = np.asarray(t_fin, dtype=np.int64)

    # Equal-motion quantile resampling (US2): select N_SNAP snapshots at
    # equal quantiles of cumulative node displacement so truss playback
    # speed reads as constant, instead of fixed-index sampling that
    # over-represents the fast-converging tail.
    Pn = [np.asarray(p)[:, :2] for (p, _t) in hist]
    d_iter = np.zeros(len(Pn))
    for i in range(1, len(Pn)):
        d_iter[i] = np.mean(np.linalg.norm(Pn[i] - Pn[i - 1], axis=1))
    Dcum = np.cumsum(d_iter)
    if Dcum[-1] <= 0:
        sel = np.unique(np.linspace(0, len(hist) - 1, N_SNAP).astype(int))
    else:
        targets = np.linspace(0.0, Dcum[-1], N_SNAP)
        sel = np.clip(np.searchsorted(Dcum, targets), 0, len(hist) - 1)
        sel[0] = 0
        sel[-1] = len(hist) - 1
        sel = np.unique(sel)
        if len(sel) < N_SNAP:
            sel = np.unique(np.concatenate(
                [sel, np.linspace(0, len(hist) - 1, N_SNAP).astype(int)]))
    snap_idx = sel.tolist()
    # sel[-1] == len(hist) - 1 already selects the final (re-triangulated)
    # iterate as the converged-truss frame — do not append (p_out, t_out).
    snapshots = [hist[i] for i in snap_idx]
    print(f"  {len(hist)} solver iterations captured; {len(snapshots)} animation snapshots")

    # Node bookkeeping: solver order is [boundary; interior].
    nb = len(bnd)
    P_seed = np.vstack([bnd, interior_pts])
    bidx = np.arange(nb)
    interior = np.arange(nb, len(P_seed))

    # F-11: the seed stage must show the TRUE pre-solver scattered state — a
    # fresh Delaunay of the unrelaxed seed points, filtered by the same
    # centroid-in-domain keep test the solver capture uses — NOT
    # snapshots[0][1], which under snapshot_retriangulate=True is already
    # the post-iteration-0-move triangulation.
    _tri_seed = Delaunay(P_seed)
    _ta_seed = _tri_seed.simplices
    _cen_seed = P_seed[_ta_seed].mean(axis=1)
    T_seed = _ta_seed[fd(_cen_seed) < 0].astype(np.int64)

    print("FEM smoothing the truss output...")
    m = Mesh(connectivity=np.c_[t_fin, t_fin[:, 0]].astype(np.int64),
             points=np.c_[p_fin, np.zeros(len(p_fin))],
             compute_layers=True, compute_adjacencies=True)
    m.smooth_mesh(method="fem", acknowledge_change=True)
    P_fem = m.points[:, :2].copy()

    print("peel_layers() on the smoothed mesh...")
    layers = m.Layers
    n_layers = m.n_layers
    elem_layer = np.zeros(m.n_elems, dtype=np.int32)
    for li in range(n_layers):
        for e in layers["OE"][li]:
            elem_layer[e] = li
        for e in layers["IE"][li]:
            elem_layer[e] = li

    # Boundary-inward reveal-order guard: the peel stage assumes elem_layer
    # increases moving inward from the domain boundary. Verify by mean
    # distance to the nearest annulus rim, NOT raw radius — radius is
    # non-monotonic on this domain (layer 0 hugs both the inner and outer
    # rim). A future chilmesh layer-reordering that breaks this assumption
    # must fail loudly here instead of silently rendering an
    # inward-to-outward reveal.
    _cen_fem = P_fem[t_fin].mean(axis=1)
    _r_fem = np.linalg.norm(_cen_fem, axis=1)
    _rim_dist = np.minimum(np.abs(_r_fem - R_IN), np.abs(_r_fem - R_OUT))
    _mean_rim_by_layer = [float(np.mean(_rim_dist[elem_layer == li])) for li in range(n_layers)]
    assert all(_mean_rim_by_layer[i] < _mean_rim_by_layer[i + 1]
               for i in range(n_layers - 1)), (
        f"peel layers are not boundary-inward by mean rim distance: {_mean_rim_by_layer}"
    )

    q_seed = _quality(P_seed, T_seed)
    q_truss = _quality(p_fin, t_fin)
    q_fem = _quality(P_fem, t_fin)

    ymax = 0
    for q in [q_seed, q_truss, q_fem] + [_quality(p, t) for p, t in snapshots[::6]]:
        c, _ = np.histogram(q, bins=HBINS, range=(0.0, 1.0))
        ymax = max(ymax, int(c.max()))
    ymax = int(ymax * 1.08) + 1

    return dict(T=T_seed, bidx=bidx, interior=interior,
                P_seed=P_seed, P_truss=p_fin, T_truss=t_fin, P_fem=P_fem,
                q_seed=q_seed, q_truss=q_truss, q_fem=q_fem,
                snapshots=snapshots,
                elem_layer=elem_layer, n_layers=n_layers, ymax=ymax,
                hist=hist, fd=fd, snap_idx=snap_idx,
                drop=0.9, fall_delay=rng.uniform(0.0, 0.55, size=len(interior)))


def _ease(t):
    t = min(1.0, max(0.0, t))
    return t * t * (3 - 2 * t)


def _q_colors(q):
    return matplotlib.colormaps[QCMAP](Normalize(0.0, 1.0)(q))


def _layer_color(li, n_layers):
    """Single-layer viridis color for the peel reveal (T013)."""
    return matplotlib.colormaps["viridis"](Normalize(0, max(1, n_layers - 1))(li))


def _peel_facecolors(q, elem_layer, n_layers, k):
    """Mesh facecolors for peel state k: quality-colored everywhere, with
    revealed layers (elem_layer < k) overwritten to their layer color."""
    fc = matplotlib.colormaps[QCMAP](Normalize(0.0, 1.0)(q))
    for li in range(k):
        fc[elem_layer == li] = _layer_color(li, n_layers)
    return fc


def _draw_peel_hist(ax_hist, q, elem_layer, n_layers, k, D):
    """Stacked-bar histogram for peel state k.

    Revealed layers (0..k-1) each render as one stacked segment in their
    layer color; the remainder (elem_layer >= k) is colored PER BIN,
    pinned exactly equal to _draw_hist's bar_colors so that at k=0 this
    renders pixel-for-pixel identical to the plain quality histogram
    (FR-006 convert-in-place). Per-bin totals are invariant across k —
    every element is assigned to exactly one stacked segment per bin.
    """
    counts_total, hedges = np.histogram(q, bins=HBINS, range=(0.0, 1.0))
    widths = np.diff(hedges)
    mids = hedges[:-1] + widths / 2.0
    remainder_colors = matplotlib.colormaps[QCMAP](Normalize(0.0, 1.0)(mids))

    # Per-element bin index (F-10): a q=1.0 element must land in the last
    # bin, never overflow past HBINS.
    b = np.clip(np.floor(q * HBINS).astype(int), 0, HBINS - 1)

    bottom = np.zeros(HBINS)
    for li in range(k):
        mask = elem_layer == li
        counts_li = np.zeros(HBINS, dtype=int)
        np.add.at(counts_li, b[mask], 1)
        ax_hist.bar(hedges[:-1], counts_li, width=widths, align="edge",
                    bottom=bottom, color=_layer_color(li, n_layers),
                    edgecolor="#1a1a1f", linewidth=0.3)
        bottom = bottom + counts_li

    mask_rem = elem_layer >= k
    counts_rem = np.zeros(HBINS, dtype=int)
    np.add.at(counts_rem, b[mask_rem], 1)
    ax_hist.bar(hedges[:-1], counts_rem, width=widths, align="edge",
                bottom=bottom, color=remainder_colors,
                edgecolor="#1a1a1f", linewidth=0.3)


def _setup_axes(ax_mesh, ax_hist, D):
    ax_mesh.set_xlim(-1.12, 1.12); ax_mesh.set_ylim(-1.12, 1.12)
    ax_mesh.set_aspect("equal"); ax_mesh.set_facecolor(BG)
    ax_mesh.set_xticks([]); ax_mesh.set_yticks([])
    for s in ax_mesh.spines.values():
        s.set_color(DIM)
    ax_hist.set_xlim(0.0, 1.0); ax_hist.set_ylim(0, D["ymax"])
    ax_hist.set_facecolor(BG)
    ax_hist.set_xlabel("Element quality", color=TEXT, fontsize=11)
    ax_hist.set_ylabel("Count", color=TEXT, fontsize=11)
    ax_hist.tick_params(colors=DIM, labelsize=9)
    for s in ax_hist.spines.values():
        s.set_color(DIM)


def _draw_metrics(ax, q, alpha=1.0):
    """Color-keyed Median/Mean/Min annotations (US5): green dotted line +
    green text for Median, red dotted line + red text for Mean, neutral
    text for Min. zorder=20 keeps both lines individually distinguishable
    even when median and mean fall close together."""
    med, mn, mean = float(np.median(q)), float(np.min(q)), float(np.mean(q))
    ax.axvline(med, color=GOOD, linestyle=":", linewidth=1.8, alpha=0.9 * alpha, zorder=20)
    ax.axvline(mean, color=RED, linestyle=":", linewidth=1.8, alpha=0.9 * alpha, zorder=20)
    ax.text(0.02, 1.06, f"Median: {med:.3f}", transform=ax.transAxes, color=GOOD,
            fontsize=11, fontweight="bold", ha="left", va="bottom", alpha=alpha)
    ax.text(0.40, 1.06, f"Mean: {mean:.3f}", transform=ax.transAxes, color=RED,
            fontsize=11, fontweight="bold", ha="left", va="bottom", alpha=alpha)
    ax.text(0.76, 1.06, f"Min: {mn:.3f}", transform=ax.transAxes, color=TEXT,
            fontsize=11, ha="left", va="bottom", alpha=alpha)


def _draw_hist(ax_hist, q, D, alpha=1.0):
    if q is None:
        return
    counts, hedges = np.histogram(q, bins=HBINS, range=(0.0, 1.0))
    widths = np.diff(hedges)
    mids = hedges[:-1] + widths / 2.0
    bar_colors = matplotlib.colormaps[QCMAP](Normalize(0.0, 1.0)(mids))
    ax_hist.bar(hedges[:-1], counts, width=widths, align="edge",
                color=bar_colors, edgecolor="#1a1a1f", linewidth=0.3, alpha=alpha)
    _draw_metrics(ax_hist, q, alpha=alpha)


def _draw_mesh(ax, pts, elems, facecolors=None, wire=True,
               interior=None, bidx=None, dot_alpha=1.0, fall=None):
    if facecolors is not None:
        ax.add_collection(PolyCollection([pts[e] for e in elems], facecolors=facecolors,
                                         edgecolors="#1a1a1f", linewidths=0.5))
    if wire:
        E = np.vstack([elems[:, [0, 1]], elems[:, [1, 2]], elems[:, [2, 0]]])
        ax.add_collection(LineCollection(pts[E],
                                         colors=WIRE if facecolors is None else EDGE,
                                         linewidths=0.3, alpha=0.5))
    if bidx is not None:
        for (x, y) in pts[bidx]:
            ax.add_patch(Circle((x, y), 0.012, color=BOUNDARY_DOT, alpha=dot_alpha, zorder=11))
    if interior is not None:
        if fall is None:
            for (x, y) in pts[interior]:
                ax.add_patch(Circle((x, y), 0.012, color=VERTEX_DOT, alpha=dot_alpha, zorder=10))
        else:
            prog, drop, delay = fall
            for k, (x, y) in enumerate(pts[interior]):
                lt = _ease((prog - delay[k]) / 0.45)
                if lt <= 0.0:
                    continue
                ax.add_patch(Circle((x, y + (1.0 - lt) * drop), 0.012,
                                    color=VERTEX_DOT, alpha=lt, zorder=10))


def _optimize_gif(path, colors=256):
    """Re-encode the saved GIF through a single global adaptive palette so
    the output stays under the size budget with no dimension, frame-count,
    or duration loss (plan §D — measured -31.6% on the shipped GIF)."""
    im = Image.open(path)
    frames = [f.convert("RGB").copy() for f in ImageSequence.Iterator(im)]
    durs = [f.info.get("duration", 100) for f in ImageSequence.Iterator(Image.open(path))]
    stack = np.concatenate(
        [np.asarray(frames[i]) for i in range(0, len(frames), max(1, len(frames) // 8))],
        axis=0,
    )
    pal = Image.fromarray(stack).convert("P", palette=Image.ADAPTIVE, colors=colors)
    q = [f.quantize(palette=pal, dither=Image.NONE) for f in frames]
    q[0].save(path, save_all=True, append_images=q[1:],
              duration=durs, loop=0, optimize=True, disposal=1)


def main():
    D = _stage_data()
    T = D["T"]; interior = D["interior"]; bidx = D["bidx"]
    P_seed, P_truss, P_fem = D["P_seed"], D["P_truss"], D["P_fem"]
    T_truss = D["T_truss"]
    q_seed, q_truss, q_fem = D["q_seed"], D["q_truss"], D["q_fem"]
    snapshots = D["snapshots"]
    n_layers = D["n_layers"]
    elem_layer = D["elem_layer"]

    # Uniform-hold truss play list (US2): 3 preroll seed frames, then every
    # equal-motion-resampled snapshot held TRUSS_HOLD frames — replaces the
    # old tiered hold (3/2/1) that over-weighted the fast-converging tail.
    play = [(P_seed, D["T"])] * 3
    for snap in snapshots:
        play += [snap] * TRUSS_HOLD

    fc_fem_q = _q_colors(q_fem)

    fig = plt.figure(figsize=(12, 5.5), facecolor=BG)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.18,
                          left=0.04, right=0.97, top=0.86, bottom=0.12)
    ax_mesh = fig.add_subplot(gs[0, 0])
    ax_hist = fig.add_subplot(gs[0, 1])
    fig.suptitle("CHILmesh pipeline · annulus", color=TEXT, fontsize=14,
                 fontweight="bold", y=0.965)

    def label(title, subtitle):
        ax_mesh.text(0.02, 1.10, title, transform=ax_mesh.transAxes, color=ACCENT,
                     fontsize=13, fontweight="bold", ha="left", va="bottom")
        ax_mesh.text(0.02, 1.035, subtitle, transform=ax_mesh.transAxes, color=TEXT,
                     fontsize=10, ha="left", va="bottom")

    schedule = []

    def add(n, fn):
        for i in range(n):
            schedule.append((fn, i, n))

    def f_bfade(i, n):
        _setup_axes(ax_mesh, ax_hist, D)
        a = _ease(i / max(1, n - 1))
        for (x, y) in P_seed[bidx]:
            ax_mesh.add_patch(Circle((x, y), 0.012, color=BOUNDARY_DOT, alpha=a, zorder=11))
        label("Stage 0/4: Seed", "boundary nodes define the annulus domain")

    def f_fall(i, n):
        _setup_axes(ax_mesh, ax_hist, D)
        prog = i / max(1, n - 1)
        _draw_mesh(ax_mesh, P_seed, T, facecolors=None, wire=False,
                   interior=interior, bidx=bidx, fall=(prog, D["drop"], D["fall_delay"]))
        label("Stage 0/4: Seed", "random interior nodes fall onto the paper")

    def f_settle(i, n):
        _setup_axes(ax_mesh, ax_hist, D)
        _draw_mesh(ax_mesh, P_seed, T, facecolors=None, wire=False,
                   interior=interior, bidx=bidx)
        label("Stage 0/4: Seed", "scattered seed — boundary + interior nodes")

    def f_wire(i, n):
        _setup_axes(ax_mesh, ax_hist, D)
        a = _ease(i / max(1, n - 1))
        E = np.vstack([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]])
        ax_mesh.add_collection(LineCollection(P_seed[E], colors=WIRE,
                                              linewidths=0.3, alpha=0.5 * a))
        _draw_mesh(ax_mesh, P_seed, T, facecolors=None, wire=False,
                   interior=interior, bidx=bidx)
        _draw_hist(ax_hist, q_seed, D, alpha=a)
        label("Stage 0/4: Seed", "unrelaxed triangulation of the scattered seed")

    def f_truss(i, n):
        # REAL solver playback: frame i renders captured iteration snapshot i
        # (points AND that iteration's own triangulation — retriangulations
        # appear exactly when the solver performed them).
        _setup_axes(ax_mesh, ax_hist, D)
        p_i, t_i = play[min(i, len(play) - 1)]
        q = _quality(p_i, t_i)
        _draw_mesh(ax_mesh, p_i, t_i, facecolors=_q_colors(q), wire=True,
                   interior=interior, bidx=bidx)
        _draw_hist(ax_hist, q, D)
        label("Stage 1/4: ADMESH Truss",
              "spring relaxation — actual solver iterations")

    def f_truss_hold(i, n):
        _setup_axes(ax_mesh, ax_hist, D)
        _draw_mesh(ax_mesh, P_truss, T_truss, facecolors=_q_colors(q_truss), wire=True,
                   interior=interior, bidx=bidx)
        _draw_hist(ax_hist, q_truss, D)
        label("Stage 1/4: ADMESH Truss", "converged — relaxed truss equilibrium")

    def f_fem(i, n):
        # Direct method: interpolate between the REAL smoother input (truss
        # output) and its REAL output, on the truss connectivity.
        _setup_axes(ax_mesh, ax_hist, D)
        t = _ease(i / max(1, n - 1))
        pts = (1 - t) * P_truss + t * P_fem
        q = _quality(pts, T_truss)
        _draw_mesh(ax_mesh, pts, T_truss, facecolors=_q_colors(q), wire=True,
                   interior=interior, bidx=bidx)
        _draw_hist(ax_hist, q, D)
        label("Stage 2/4: FEM Smoother",
              "Balendran direct method — faces deform as nodes settle")

    def f_fem_hold(i, n):
        _setup_axes(ax_mesh, ax_hist, D)
        _draw_mesh(ax_mesh, P_fem, T_truss, facecolors=fc_fem_q, wire=True,
                   interior=interior, bidx=bidx)
        _draw_hist(ax_hist, q_fem, D)
        label("Stage 2/4: FEM Smoother", "smoothed — element quality improved")

    def f_peel(k):
        # Boundary-inward layer reveal (US4): mesh recolors instantly at
        # each reveal boundary via _peel_facecolors; the histogram converts
        # to matching stacked-bar segments via _draw_peel_hist, with the
        # median/mean/min overlay drawn explicitly afterward (peel frames
        # do NOT call _draw_hist — F-01).
        def _f(i, n):
            _setup_axes(ax_mesh, ax_hist, D)
            fc = _peel_facecolors(q_fem, elem_layer, n_layers, k)
            _draw_mesh(ax_mesh, P_fem, T_truss, facecolors=fc, wire=True,
                       interior=interior, bidx=bidx)
            _draw_peel_hist(ax_hist, q_fem, elem_layer, n_layers, k, D)
            _draw_metrics(ax_hist, q_fem)
            if 0 < k < n_layers:
                subtitle = f"peel_layers(): revealing layer {k}/{n_layers}, boundary-inward"
            else:
                subtitle = f"peel_layers(): {n_layers} concentric layers (OE+IE per ring)"
            label("Stage 3/4: Peel Layers", subtitle)
        return _f

    add(8, f_bfade)
    add(26, f_fall)
    add(16, f_settle)
    add(10, f_wire)
    add(18, lambda i, n: f_wire(n - 1, n))   # Seed→Truss pause (F-13: 18 f = 1.8 s)
    add(len(play), f_truss)
    add(18, f_truss_hold)                    # Truss→FEM pause (F-13: 18 f = 1.8 s)
    add(24, f_fem)
    add(16, f_fem_hold)
    # Peel reveal schedule (US4): k=0 hold doubles as the FEM→Peel pause
    # (F-13: 18 f = 1.8 s, byte-identical to the FEM-hold histogram), then
    # each layer reveals instantly and holds 8 f (0.8 s), ending on a 50 f
    # (5 s) held full-layer view. No re-peel.
    add(18, f_peel(0))
    for k in range(1, n_layers):
        add(8, f_peel(k))
    add(50, f_peel(n_layers))

    total = len(schedule)

    def animate(idx):
        ax_mesh.clear(); ax_hist.clear()
        fn, i, n = schedule[idx]
        fn(i, n)
        return []

    print(f"Rendering {total} frames → {OUT_PATH}...")
    anim = animation.FuncAnimation(fig, animate, frames=total, interval=100, blit=False)
    anim.save(str(OUT_PATH), writer=animation.PillowWriter(fps=10), dpi=72)
    plt.close(fig)
    _optimize_gif(str(OUT_PATH), colors=256)
    print(f"Done: {OUT_PATH} ({OUT_PATH.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
