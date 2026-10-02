#!/usr/bin/env python
"""Benchmark the GPU plot backend against matplotlib (CHILmesh #167, Phase A).

Run from a venv that has the ``[gpu]`` extra::

    python scripts/benchmark_gpu_render.py

Measures, as the median of ``--reps`` runs (default 5), on every mesh found:

* first frame: build + upload + render + read back RGBA (fresh view each rep)
* steady re-render: camera pan/zoom, then render + read back (geometry resident)
* scalar update: new per-element values, then render + read back
  (GPU ``set_scalar``; matplotlib ``PolyCollection.set_array`` + ``canvas.draw``)

Both backends draw a scalar-coloured filled mesh with edges into an 800x600
canvas and return pixels to the CPU, so the timings cover the same work.
Mesh loading and connectivity parsing are not timed. A matplotlib mesh whose
first frame exceeds ``--mpl-limit`` seconds is skipped for the rest of its
measurements and reported as such.

Writes ``docs/benchmarks/gpu_render.md``.
"""
from __future__ import annotations

import argparse
import platform
import statistics
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "benchmarks" / "gpu_render.md"
MESHES = [
    ("Block_O", ROOT / "src" / "chilmesh" / "data" / "Block_O.14"),
    ("WNAT_Hagen", Path("/Users/domattioli/Projects/Valence/registry_data/meshes/WNAT_Hagen.14")),
    ("Chesapeake_Bay", Path("/Users/domattioli/Projects/Valence/registry_data/meshes/Chesapeake_Bay.14")),
    ("EasternPacific_ENPAC2003", Path("/Users/domattioli/Projects/Valence/registry_data/meshes/EasternPacific_ENPAC2003.14")),
]
SIZE = (800, 600)


def load(path: Path):
    from chilmesh.fort14_io import read_fort14_raw

    raw = read_fort14_raw(str(path), parse_boundaries=False)
    index = {nid: i for i, nid in enumerate(raw.node_ids)}
    pts = np.array([raw.coords[n][:2] for n in raw.node_ids], dtype=float)
    arity = max(len(raw.elements[e]) for e in raw.elem_ids)
    conn = np.empty((len(raw.elem_ids), arity), dtype=np.int64)
    for i, eid in enumerate(raw.elem_ids):
        row = [index[n] for n in raw.elements[eid]]
        row += [row[-1]] * (arity - len(row))  # pad triangles with a repeated vertex
        conn[i] = row
    return pts, conn


def timed(fn) -> float:
    t0 = time.perf_counter()
    fn()
    return time.perf_counter() - t0


def med(xs):
    return statistics.median(xs) * 1e3  # ms


def bench_gpu(pts, conn, vals, reps):
    from chilmesh.chilplotting_gpu import GPUMeshView

    first = []
    for _ in range(reps):
        first.append(timed(lambda: GPUMeshView(
            pts, conn, values=vals, size=SIZE).render()))
    view = GPUMeshView(pts, conn, values=vals, size=SIZE)
    view.render()  # warm: pipeline compile is a one-time cost, reported via 'first'
    x0, x1, y0, y1 = view._bounds
    dx = (x1 - x0) * 0.01
    steady = []
    for i in range(reps):
        def frame(i=i):
            view.pan(dx * (1 if i % 2 else -1), 0.0)
            view.zoom(1.02 if i % 2 else 1 / 1.02)
            view.render()
        steady.append(timed(frame))
    rng = np.random.default_rng(0)
    scalar = []
    for _ in range(reps):
        new = rng.random(len(conn))
        scalar.append(timed(lambda: (view.set_scalar(new), view.render())))
    return med(first), med(steady), med(scalar)


def bench_mpl(pts, conn, vals, reps, limit):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from chilmesh import chilplotting as cp

    def build():
        fig = Figure(figsize=(SIZE[0] / 100, SIZE[1] / 100), dpi=100)
        FigureCanvasAgg(fig)
        ax = fig.add_axes([0, 0, 1, 1])
        cp.plot_filled(pts, conn, values=vals, ax=ax)
        fig.canvas.draw()
        np.asarray(fig.canvas.buffer_rgba())
        return fig, ax

    t_first = timed(build)
    if t_first > limit:
        return None, t_first
    first = [t_first] + [timed(build) for _ in range(reps - 1)]
    fig, ax = build()
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    steady = []
    for i in range(reps):
        def frame(i=i):
            s = 1.02 if i % 2 else 1 / 1.02
            dx = (x1 - x0) * 0.01 * (1 if i % 2 else -1)
            cx, cy = (x0 + x1) / 2 + dx, (y0 + y1) / 2
            hw, hh = (x1 - x0) / 2 / s, (y1 - y0) / 2 / s
            ax.set_xlim(cx - hw, cx + hw)
            ax.set_ylim(cy - hh, cy + hh)
            fig.canvas.draw()
            np.asarray(fig.canvas.buffer_rgba())
        steady.append(timed(frame))
    pc = ax.collections[0]
    rng = np.random.default_rng(0)
    scalar = []
    for _ in range(reps):
        new = rng.random(len(conn))

        def upd(new=new):
            pc.set_array(new)
            fig.canvas.draw()
            np.asarray(fig.canvas.buffer_rgba())
        scalar.append(timed(upd))
    return (med(first), med(steady), med(scalar)), t_first


def fmt(v):
    return "skipped" if v is None else f"{v:,.1f}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--mpl-limit", type=float, default=60.0,
                    help="skip matplotlib on a mesh whose first frame exceeds this many seconds")
    args = ap.parse_args()
    if args.reps < 5:
        ap.error("--reps must be >= 5")

    sys.path.insert(0, str(ROOT / "src"))
    from chilmesh import chilplotting_gpu as gpu

    info = gpu.gpu_adapter_info()
    if info is None:
        print("No wgpu adapter available; cannot benchmark.", file=sys.stderr)
        return 1
    print(f"adapter: {info['summary']}")

    rows, notes = [], []
    for name, path in MESHES:
        if not path.exists():
            print(f"skip {name}: {path} missing")
            notes.append(f"{name}: file not found, skipped.")
            continue
        pts, conn = load(path)
        vals = np.random.default_rng(1).random(len(conn))
        print(f"{name}: {len(conn):,} elems, {len(pts):,} nodes", flush=True)
        g = bench_gpu(pts, conn, vals, args.reps)
        m, t_first = bench_mpl(pts, conn, vals, args.reps, args.mpl_limit)
        if m is None:
            notes.append(f"{name}: matplotlib first frame took {t_first:.0f} s "
                         f"(> {args.mpl_limit:.0f} s limit); matplotlib skipped.")
        rows.append((name, len(conn), len(pts), g, m))
        print(f"  gpu first/steady/scalar ms: {g}", flush=True)
        print(f"  mpl: {m}", flush=True)

    lines = [
        "# GPU plot backend benchmark (CHILmesh #167, Phase A)",
        "",
        "Generated by `scripts/benchmark_gpu_render.py`. Do not edit by hand.",
        "",
        "## Environment",
        "",
        f"- Adapter: `{info['summary']}`",
        f"- Adapter info: `{ {k: v for k, v in info.items() if k not in ('summary',)} }`",
        f"- pygfx {info['pygfx_version']}, wgpu {info['wgpu_version']}",
        f"- Python {platform.python_version()}, numpy {np.__version__}, "
        f"{platform.platform()}, machine {platform.machine()}",
        "",
        "## Method",
        "",
        f"- Canvas {SIZE[0]}x{SIZE[1]} px. Both backends draw a per-element scalar "
        "(viridis) fill plus edges and copy the RGBA pixels back to a numpy array, "
        f"so the work compared is the same. Median of {args.reps} runs, milliseconds.",
        "- **First frame**: build the view, upload buffers, render, read back. "
        "GPU: new `GPUMeshView` each run (pipelines are "
        "cached in-process, so only the first run pays shader compilation; the median excludes it). "
        "matplotlib: new `Figure` + `plot_filled` + `canvas.draw`.",
        "- **Steady re-render**: pan and zoom the camera, then render and read back. "
        "GPU geometry stays resident. matplotlib: change axes limits + `canvas.draw`.",
        "- **Scalar update**: new random per-element values, then render and read back. "
        "GPU: `GPUMeshView.set_scalar`. matplotlib: `PolyCollection.set_array` + `canvas.draw`.",
        "- Mesh parsing is not timed. The matplotlib path draws every polygon outline; "
        "the GPU path draws each unique edge once as a line object. Elements are "
        "triangulated for the GPU (quads become two triangles).",
        "- GPU timings include the readback stall that offscreen capture needs, so "
        "they are an upper bound on an interactive on-screen frame.",
        "",
        "## Results (ms, median)",
        "",
        "| Mesh | Elements | Nodes | GPU first | GPU steady | GPU set_scalar "
        "| mpl first | mpl steady | mpl scalar | First speedup | Steady speedup |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name, ne, nn, g, m in rows:
        mf, ms, mc = (m if m else (None, None, None))
        sp1 = f"{mf / g[0]:.1f}x" if m else "n/a"
        sp2 = f"{ms / g[1]:.1f}x" if m else "n/a"
        lines.append(
            f"| {name} | {ne:,} | {nn:,} | {fmt(g[0])} | {fmt(g[1])} | {fmt(g[2])} "
            f"| {fmt(mf)} | {fmt(ms)} | {fmt(mc)} | {sp1} | {sp2} |")
    if notes:
        lines += ["", "## Notes", ""] + [f"- {n}" for n in notes]
    lines += [
        "",
        "Speedup is matplotlib time divided by GPU time; a value below 1x means the "
        "GPU path is slower.",
        "",
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines))
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
