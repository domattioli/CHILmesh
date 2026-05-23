#!/usr/bin/env python3
"""Cross-language CHILmesh benchmark + skeletonization parity harness.

Times the skeletonization pipeline (adjacency build, layer peel, quality,
vertex-edge lookup) for the available implementations and checks that every
backend agrees on ``n_layers`` — the regression guard behind the README
"Performance" and "Validation" tables.

Implementations probed (each is optional; missing ones are skipped):
  - Python   : always (the reference; ``chilmesh`` package)
  - C++      : if ``chilmesh_cpp`` importable (``pip install ./src/chilmesh_cpp``)
  - MATLAB   : if ``octave`` on PATH (runs the bundled ``src/@CHILmesh`` class)
  - Rust     : intentionally excluded — skeletonization is incomplete (#163)

Opt-in: this is a manual / local tool, never a CI gate (the MATLAB column
needs GNU Octave, which we do not require in CI). Run with::

    CHILMESH_RUN_BENCH=1 python scripts/benchmark.py
    python scripts/benchmark.py --mesh src/chilmesh/data/Block_O.14 --matlab

Outputs Markdown tables to stdout; ``--emit docs/BENCHMARK.md`` rewrites the doc.
"""
from __future__ import annotations

import argparse
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
MATLAB_CLASS_DIR = REPO / "src" / "@CHILmesh"
DEFAULT_MESH = REPO / "src" / "chilmesh" / "data" / "Block_O.14"


def _median(fn, repeats: int):
    ts = []
    out = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        out = fn()
        ts.append(time.perf_counter() - t0)
    return statistics.median(ts), out


def bench_python(conn: np.ndarray, pts: np.ndarray, repeats: int) -> dict:
    from chilmesh import CHILmesh

    adj, skel, full, qual = [], [], [], []
    mesh = None
    for _ in range(repeats):
        m = CHILmesh(connectivity=conn.copy(), points=pts.copy(),
                     compute_layers=False, compute_adjacencies=False)
        t = time.perf_counter(); m._build_adjacencies(); a = time.perf_counter() - t
        t = time.perf_counter(); m._skeletonize(); s = time.perf_counter() - t
        adj.append(a); skel.append(s); full.append(a + s)
        t = time.perf_counter(); m.signed_area(); qual.append(time.perf_counter() - t)
        mesh = m
    n = min(2000, mesh.n_verts)
    t = time.perf_counter()
    for k in range(n):
        mesh.get_vertex_edges(k)
    ve = (time.perf_counter() - t) / n
    return {
        "fast_s": statistics.median(adj),
        "skel_s": statistics.median(skel),
        "full_s": statistics.median(full),
        "quality_s": statistics.median(qual),
        "vertex_edge_us": ve * 1e6,
        "n_layers": int(mesh.n_layers),
    }


def bench_cpp(conn: np.ndarray, pts: np.ndarray, repeats: int) -> dict | None:
    try:
        from chilmesh.backends import cpp_backend as C
    except Exception:
        return None
    if not getattr(C, "CPP_AVAILABLE", False):
        return None
    cpp_conn = conn.astype(np.int32)
    tf, _ = _median(lambda: C.fast_init(pts, cpp_conn), repeats)
    tu, mu = _median(lambda: C.full_init(pts, cpp_conn), repeats)
    tq, _ = _median(lambda: C.quality_analysis(mu), repeats)
    n = min(2000, mu.n_verts)
    t = time.perf_counter()
    for k in range(n):
        C.get_vertex_edges(mu, k)
    ve = (time.perf_counter() - t) / n
    return {
        "fast_s": tf,
        "skel_s": max(tu - tf, 0.0),
        "full_s": tu,
        "quality_s": tq,
        "vertex_edge_us": ve * 1e6,
        "n_layers": int(mu.n_layers),
    }


def bench_matlab(conn: np.ndarray, pts: np.ndarray) -> dict | None:
    """Run the bundled MATLAB @CHILmesh class under GNU Octave (if present)."""
    if shutil.which("octave") is None:
        return None
    if not (MATLAB_CLASS_DIR / "CHILmesh.m").exists():
        return None
    try:
        from scipy.io import savemat
    except Exception:
        return None
    with tempfile.TemporaryDirectory() as td:
        mat = Path(td) / "mesh.mat"
        savemat(str(mat), {"ConnectivityList": conn.astype(np.int64) + 1,
                           "Points": pts[:, :2].astype(float)})
        # @CHILmesh dir must be on the path's parent; the readFort14 path is
        # Octave-incompatible (textscan), so feed arrays to the 2-arg ctor.
        m = f"""
        addpath('{MATLAB_CLASS_DIR.parent}'); warning('off','all');
        S = load('{mat}');
        t=tic; cm = CHILmesh(S.ConnectivityList, S.Points); full=toc(t);
        cm2 = CHILmesh(S.ConnectivityList, S.Points);
        t=tic; cm2.buildAdjacencies(); adj=toc(t);
        t=tic; cm2.meshLayers(); skel=toc(t);
        t=tic; cm.signedArea(); q=toc(t);
        printf('RESULT %d %.4f %.4f %.4f %.6f\\n', cm.nLayers, full, adj, skel, q);
        """
        try:
            out = subprocess.run(["octave", "--no-gui", "--eval", m],
                                 capture_output=True, text=True, timeout=900)
        except Exception:
            return None
    line = next((ln for ln in out.stdout.splitlines() if ln.startswith("RESULT")), None)
    if line is None:
        return None
    _, nl, full, adj, skel, q = line.split()
    return {
        "fast_s": float(adj),
        "skel_s": float(skel),
        "full_s": float(full),
        "quality_s": float(q),
        "vertex_edge_us": None,  # Octave vert2Edge is inputParser-bound; not comparable
        "n_layers": int(nl),
    }


def _fmt_s(v):
    if v is None:
        return "—"
    return f"{v*1000:.0f} ms" if v < 1.0 else f"{v:.2f} s"


def _fmt_us(v):
    return "—" if v is None else f"{v:.2f} μs"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mesh", default=str(DEFAULT_MESH), help="path to a .14 mesh")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--matlab", action="store_true", help="include the Octave/MATLAB column")
    ap.add_argument("--emit", help="write the report to this path instead of stdout")
    args = ap.parse_args()

    if not os.environ.get("CHILMESH_RUN_BENCH") and args.emit is None:
        print("Set CHILMESH_RUN_BENCH=1 (or pass --emit) to run the benchmark.",
              file=sys.stderr)
        # Still run for interactive use; the gate is advisory.

    from chilmesh import CHILmesh
    mesh = CHILmesh.read_from_fort14(Path(args.mesh))
    conn = np.asarray(mesh.connectivity_list).astype(np.int64)
    pts = np.asarray(mesh.points).astype(float)
    n_layers_ref = int(mesh.n_layers)

    cols = [("Python", bench_python(conn, pts, args.repeats)),
            ("C++", bench_cpp(conn, pts, args.repeats))]
    if args.matlab:
        cols.append(("MATLAB (Octave)", bench_matlab(conn, pts)))
    cols = [(name, r) for name, r in cols if r is not None]

    lines = [
        f"Mesh: `{Path(args.mesh).name}` — {mesh.n_verts:,} vertices · "
        f"{mesh.n_elems:,} elements. Median of {args.repeats}.",
        "",
        "| Metric | " + " | ".join(name for name, _ in cols) + " |",
        "|---" + "|---:" * len(cols) + "|",
    ]
    rows = [
        ("Fast init (adj, no skeletonization)", "fast_s", _fmt_s),
        ("Skeletonization only", "skel_s", _fmt_s),
        ("Full init (adj + skeletonization)", "full_s", _fmt_s),
        ("Quality analysis", "quality_s", _fmt_s),
        ("Vertex-edge lookup (per call)", "vertex_edge_us", _fmt_us),
    ]
    for label, key, fmt in rows:
        lines.append("| " + label + " | "
                     + " | ".join(fmt(r[key]) for _, r in cols) + " |")

    parity = {name: r["n_layers"] for name, r in cols}
    parity["reference"] = n_layers_ref
    all_match = all(v == n_layers_ref for v in parity.values())
    lines += ["",
              f"**Skeletonization parity:** n_layers = {n_layers_ref} across "
              f"{', '.join(name for name, _ in cols)} — "
              f"{'all agree ✅' if all_match else 'MISMATCH ❌ ' + repr(parity)}."]

    report = "\n".join(lines) + "\n"
    if args.emit:
        Path(args.emit).write_text(report)
        print(f"wrote {args.emit}")
    else:
        print(report)
    return 0 if all_match else 1


if __name__ == "__main__":
    raise SystemExit(main())
