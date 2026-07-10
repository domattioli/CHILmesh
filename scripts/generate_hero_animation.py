"""Generate CHILmesh README hero animation as GIF.

Sequence (single loop; #179 + #187 + operator direction 2026-07-10):
  0. Seed      — blank paper, boundary nodes appear, interior nodes FALL onto
                 the paper (scattered), wireframe draws. No colors. A held pause.
  1. Truss     — the scattered seed relaxes toward a quality mesh; nodes AND
                 faces move together (the #179 fix — faces deform with the nodes,
                 not a hard jump), quality coloring + histogram fade in.
  2. FEM       — Balendran direct smoother; the mesh settles further, faces
                 continue to deform, quality climbs.
  3. Peel      — the final peel_layers() concentric layers, shown as a single
                 held frame (no reveal animation, per operator). Loop ends here.

The truss stage plays back ACTUAL solver iterations: distmesh2d_warmstart is
run once with history_out= capturing (points, triangulation) every iteration
(deltat=0.02 — the converging regime for a cold scattered seed), and each
animation frame renders one captured snapshot, retriangulations included.
FEM + peel then run on the solver's real output mesh. The histogram y-axis is
pinned to the global max bin count across all stages.

Left panel = mesh + vertex dots; right panel = live element-quality histogram.
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

from chilmesh import Mesh, examples

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_PATH = REPO_ROOT / "docs" / "gallery" / "readme_pipeline_annulus.gif"

BG = "#0e0e10"
EDGE = "#3a7fbf"
TEXT = "#e8e8ee"
DIM = "#888888"
ACCENT = "#ff9f43"
GOOD = "#2ecc71"
VERTEX_DOT = "#ffff00"      # interior vertex tracking dots
BOUNDARY_DOT = "#ff7a1a"    # boundary nodes (distinct)
WIRE = "#5a5a66"            # colorless wireframe
QCMAP = "cool_r"
R_IN, R_OUT = 0.3, 1.0
HBINS = 40


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
    distmesh2d_warmstart(history_out=...); FEM + peel run on the solver's final
    mesh. The annulus fixture supplies only the boundary node rings."""
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
    p_out, t_out = distmesh2d_warmstart(
        bnd, interior_pts, fd, None, h0, (-1, -1, 1, 1),
        deltat=0.02, Fscale=1.2, niter=200,
        track_best_quality=False, history_out=hist, history_every=1,
    )
    p_out = np.asarray(p_out)[:, :2]
    t_out = np.asarray(t_out, dtype=np.int64)
    # Animation snapshots: the solver converges within ~15 iterations at this
    # regime, so the playback is a slow-motion replay of the real early frames
    # — each early snapshot is HELD for several animation frames (all states
    # are genuine solver iterations; only the pacing is stretched).
    snap_idx = list(range(0, 16)) + list(range(16, 40, 2)) + list(range(40, len(hist), 15))
    snapshots = [hist[i] for i in snap_idx if i < len(hist)] + [(p_out, t_out)]
    print(f"  {len(hist)} solver iterations captured; {len(snapshots)} animation snapshots")

    # Node bookkeeping: solver order is [boundary; interior].
    nb = len(bnd)
    P_seed = np.vstack([bnd, interior_pts])
    T_seed = snapshots[0][1]           # first captured triangulation
    bidx = np.arange(nb)
    interior = np.arange(nb, len(P_seed))

    print("FEM smoothing the truss output...")
    m = Mesh(connectivity=np.c_[t_out, t_out[:, 0]].astype(np.int64),
             points=np.c_[p_out, np.zeros(len(p_out))],
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

    q_seed = _quality(P_seed, T_seed)
    q_truss = _quality(p_out, t_out)
    q_fem = _quality(P_fem, t_out)

    ymax = 0
    for q in [q_seed, q_truss, q_fem] + [_quality(p, t) for p, t in snapshots[::6]]:
        c, _ = np.histogram(q, bins=HBINS, range=(0.0, 1.0))
        ymax = max(ymax, int(c.max()))
    ymax = int(ymax * 1.08) + 1

    return dict(T=T_seed, bidx=bidx, interior=interior,
                P_seed=P_seed, P_truss=p_out, T_truss=t_out, P_fem=P_fem,
                q_seed=q_seed, q_truss=q_truss, q_fem=q_fem,
                snapshots=snapshots,
                elem_layer=elem_layer, n_layers=n_layers, ymax=ymax,
                drop=0.9, fall_delay=rng.uniform(0.0, 0.55, size=len(interior)))


def _ease(t):
    t = min(1.0, max(0.0, t))
    return t * t * (3 - 2 * t)


def _q_colors(q):
    return matplotlib.colormaps[QCMAP](Normalize(0.0, 1.0)(q))


def _layer_colors(elem_layer, n_layers, reveal=1.0):
    colors = matplotlib.colormaps["viridis"](Normalize(0, max(1, n_layers - 1))(elem_layer))
    if reveal < 1.0:
        colors[elem_layer > reveal * max(1, n_layers - 1)] = np.array([0.16, 0.16, 0.18, 1.0])
    return colors


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


def _draw_hist(ax_hist, q, D, alpha=1.0):
    if q is None:
        return
    counts, hedges = np.histogram(q, bins=HBINS, range=(0.0, 1.0))
    widths = np.diff(hedges)
    mids = hedges[:-1] + widths / 2.0
    bar_colors = matplotlib.colormaps[QCMAP](Normalize(0.0, 1.0)(mids))
    ax_hist.bar(hedges[:-1], counts, width=widths, align="edge",
                color=bar_colors, edgecolor="#1a1a1f", linewidth=0.3, alpha=alpha)
    med, mn, mean = float(np.median(q)), float(np.min(q)), float(np.mean(q))
    ax_hist.axvline(med, color=GOOD, linestyle="--", linewidth=1.5, alpha=0.85 * alpha)
    ax_hist.set_title(f"Median: {med:.3f}    Min: {mn:.3f}    Mean: {mean:.3f}",
                      color=TEXT, fontsize=11, pad=8)


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


def main():
    D = _stage_data()
    T = D["T"]; interior = D["interior"]; bidx = D["bidx"]
    P_seed, P_truss, P_fem = D["P_seed"], D["P_truss"], D["P_fem"]
    T_truss = D["T_truss"]
    q_seed, q_truss, q_fem = D["q_seed"], D["q_truss"], D["q_fem"]
    snapshots = D["snapshots"]
    # Slow-motion play list: hold each early snapshot 3 frames, mid 2, late 1.
    play = [(P_seed, D["T"])] * 3
    for k, snap in enumerate(snapshots):
        play += [snap] * (3 if k < 12 else 2 if k < 24 else 1)
    fc_fem_q = _q_colors(q_fem)
    fc_layers = _layer_colors(D["elem_layer"], D["n_layers"], 1.0)

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

    def f_peel_hold(i, n):
        _setup_axes(ax_mesh, ax_hist, D)
        _draw_mesh(ax_mesh, P_fem, T_truss, facecolors=fc_layers, wire=True,
                   interior=interior, bidx=bidx)
        _draw_hist(ax_hist, q_fem, D)
        label("Stage 3/4: Peel Layers",
              f"peel_layers(): {D['n_layers']} concentric layers (OE+IE per ring)")

    add(8, f_bfade)
    add(26, f_fall)
    add(16, f_settle)
    add(10, f_wire)
    add(8, lambda i, n: f_wire(n - 1, n))
    add(len(play), f_truss)
    add(16, f_truss_hold)
    add(24, f_fem)
    add(16, f_fem_hold)
    add(40, f_peel_hold)   # final layers only — no reveal frames (operator)

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
    print(f"Done: {OUT_PATH} ({OUT_PATH.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
