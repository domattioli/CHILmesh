#!/usr/bin/env python3
"""CHILmesh benchmark — cross-language parity + full mesh lifecycle.

Two reports are produced:

1. **Cross-language skeletonization parity** (Python / C++ / MATLAB). Times the
   skeletonization pipeline (adjacency build, layer peel, quality, vertex-edge
   lookup) for the available implementations and checks every backend agrees on
   ``n_layers`` — the regression guard behind the README "Performance" /
   "Validation" tables.

2. **Python mesh lifecycle** (#155). Times every post-generation stage a mesh
   passes through — *generation itself is out of scope*:

     adjacency build → skeletonization → quality → FEM direct smoother →
     angle-based smoother → ADMESH warm-start truss optimizer

   The truss stage needs a signed-distance function (SDF) for the domain. SDFs
   for the bundled ``annulus`` / ``donut`` fixtures are built in and auto-matched
   by mesh filename; for any other mesh pass ``--sdf {annulus,donut,square}`` or
   the stage is reported as skipped. The smoothing/skeletonization stages run on
   any mesh and scale to large inputs — point ``--mesh`` at a WNAT-scale ``.14``
   to profile the lifecycle at production size.

Implementations probed for report 1 (each optional; missing ones skipped):
  - Python   : always (the reference; ``chilmesh`` package)
  - C++      : if ``chilmesh_cpp`` importable (``pip install ./src/chilmesh_cpp``)
  - MATLAB   : if ``octave`` on PATH (runs the bundled ``src/@CHILmesh`` class)
  - Rust     : [result shown below]

Opt-in: this is a manual / local tool, never a CI gate (the MATLAB column
needs GNU Octave, which we do not require in CI). Run with::

    CHILMESH_RUN_BENCH=1 python scripts/benchmark.py
    python scripts/benchmark.py --mesh src/chilmesh/data/Block_O.14 --matlab
    # full lifecycle incl. truss on a domain with a known SDF:
    python scripts/benchmark.py --mesh src/chilmesh/data/donut_domain.fort.14
    # keep heavy meshes manual-only — auto-skip the OOM-prone stages in CI:
    python scripts/benchmark.py --mesh wnat.14 --max-elements 500000

The memory/time-heavy lifecycle stages (FEM direct smoother, ADMESH truss) are
gated by ``--max-elements N`` (#155). On a mesh whose element count exceeds
``N`` those two stages are skipped with a reason instead of being run — the FEM
direct sparse solver OOM-kills at WNAT scale (~4M DOF, see #168), so a heavy
mesh would otherwise abort the whole report. ``N=0`` (the default) disables the
gate and runs every stage; CI passes a finite ``N`` so large meshes stay
manual-only profiling targets while the scalable stages still report.

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


def _load_fort14_arrays(path: Path):
    """Fast fort.14 reader — returns (conn, pts, n_elems, n_verts, grid_name).

    Bypasses CHILmesh.__init__ (which builds adjacencies even with
    compute_layers=False) so loading a 531k-element mesh takes ~seconds
    instead of hanging on _build_adjacencies().
    """
    with open(path) as f:
        lines = f.readlines()
    grid_name = lines[0].strip()
    n_elems, n_verts = map(int, lines[1].split()[:2])
    # Node lines: index 2..2+n_verts  format: id lon lat depth
    node_lines = lines[2:2 + n_verts]
    pts = np.array([ln.split()[1:4] for ln in node_lines], dtype=float)
    # Element lines: index 2+n_verts..  format: id nv v1 v2 v3
    elem_lines = lines[2 + n_verts:2 + n_verts + n_elems]
    conn = np.array([ln.split()[2:5] for ln in elem_lines], dtype=np.int64) - 1
    return conn, pts, n_elems, n_verts, grid_name


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


def bench_rust(conn: np.ndarray, pts: np.ndarray, repeats: int) -> dict | None:
    try:
        from chilmesh.backends.rust_backend import RUST_AVAILABLE, full_init
    except Exception:
        return None
    if not RUST_AVAILABLE:
        return None
    rust_conn = conn.astype(np.int32)
    # fast_init not exposed from Rust backend; measure full_init directly
    tu, mu = _median(lambda: full_init(pts, rust_conn), repeats)
    tq, _ = _median(lambda: mu.quality_analysis(), repeats)
    n = min(2000, mu.n_verts)
    t = time.perf_counter()
    for k in range(n):
        mu.get_vertex_edges(k)
    ve = (time.perf_counter() - t) / n
    return {
        "fast_s": None,  # Rust backend doesn't expose fast_init separately
        "skel_s": None,  # Skeletonization timing not available (#163)
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


# --- Analytic signed-distance functions for the bundled fixture domains. ---
# Mirror the values used by tests/test_admesh_warmstart.py so the truss-stage
# numbers are comparable with the warm-start regression suite. Keyed by the mesh
# filename stem; used by the lifecycle report's truss stage.
def _annulus_sdf(p: np.ndarray) -> np.ndarray:
    r = np.linalg.norm(p, axis=1)
    return np.maximum(r - 1.0, 0.3 - r)


def _square_sdf(p: np.ndarray) -> np.ndarray:
    return np.max(np.abs(p) - 1.0, axis=1)


_SDF_REGISTRY = {
    "annulus": _annulus_sdf,
    "donut": _annulus_sdf,   # donut fixture shares the annulus geometry (R=1, r=0.3)
    "square": _square_sdf,
}


def _resolve_sdf(mesh_path: str, override: str | None):
    """Return (name, sdf) for the truss stage, or (reason, None) if unavailable."""
    if override:
        if override == "none":
            return "disabled (--sdf none)", None
        if override in _SDF_REGISTRY:
            return override, _SDF_REGISTRY[override]
        return f"unknown sdf '{override}'", None
    stem = Path(mesh_path).name.lower()
    for key, fn in _SDF_REGISTRY.items():
        if stem.startswith(key):
            return key, fn
    return "no built-in SDF for this domain (pass --sdf)", None


def _gate_heavy(n_elems: int, max_elements: int) -> tuple[bool, str]:
    """Decide whether the OOM-prone lifecycle stages should be skipped (#155).

    ``max_elements <= 0`` disables the gate (every stage runs — the manual /
    local default). Otherwise a mesh with more than ``max_elements`` elements
    has its FEM direct smoother + ADMESH truss stages skipped, because the FEM
    direct sparse solver OOM-kills at WNAT scale (#168). Returns
    ``(skip, reason)``; ``reason`` is empty when not skipping.
    """
    if max_elements <= 0 or n_elems <= max_elements:
        return False, ""
    return True, (f"mesh has {n_elems:,} elements > --max-elements "
                  f"{max_elements:,} gate (FEM solver OOMs at scale, #168)")


def bench_lifecycle(conn: np.ndarray, pts: np.ndarray, repeats: int,
                    smooth_iters: int, skip_fem: bool = False,
                    skip_skel: bool = False) -> dict:
    """Time the post-generation mesh lifecycle stages (Python; #155).

    Generation is intentionally excluded. Each repeat rebuilds the mesh so the
    adjacency/skeletonization timings are cold-start; the smoothers run on the
    fully-initialised mesh and do not mutate it (they return new point arrays).

    ``skip_fem`` omits the FEM direct smoother stage (``fem_smooth_s`` is then
    ``None``) — used by the ``--max-elements`` gate to keep the OOM-prone solver
    off heavy meshes while the scalable stages still report (#155, #168).

    ``skip_skel`` omits skeletonization (``skel_s`` is then ``None``) — used by
    the ``--max-elements`` gate when C++ backend is unavailable and Python
    skeletonization would hang on large meshes.
    """
    from chilmesh import CHILmesh

    adj, skel, qual, fem, ang = [], [], [], [], []
    mesh = None
    for _ in range(repeats):
        m = CHILmesh(connectivity=conn.copy(), points=pts.copy(),
                     compute_layers=False, compute_adjacencies=False)
        t = time.perf_counter(); m._build_adjacencies(); adj.append(time.perf_counter() - t)
        if not skip_skel:
            t = time.perf_counter(); m._skeletonize(); skel.append(time.perf_counter() - t)
        t = time.perf_counter(); m.signed_area(); qual.append(time.perf_counter() - t)
        if not skip_fem:
            t = time.perf_counter(); m.direct_smoother(); fem.append(time.perf_counter() - t)
        t = time.perf_counter(); m.angle_based_smoother(n_iter=smooth_iters)
        ang.append(time.perf_counter() - t)
        mesh = m
    return {
        "adj_s": statistics.median(adj),
        "skel_s": None if skip_skel else statistics.median(skel),
        "skel_skipped": skip_skel,
        "quality_s": statistics.median(qual),
        "fem_smooth_s": None if skip_fem else statistics.median(fem),
        "fem_skipped": skip_fem,
        "angle_smooth_s": statistics.median(ang),
        "smooth_iters": smooth_iters,
        "n_verts": int(mesh.n_verts),
    }


def bench_truss(conn: np.ndarray, pts: np.ndarray, sdf, repeats: int,
                niter: int) -> dict | None:
    """Time the ADMESH warm-start truss optimizer (#155).

    Returns ``None`` when the mesh is not triangle-only (the optimizer requires
    a pure-triangle input). ``sdf`` is the domain signed-distance function.
    """
    from chilmesh import CHILmesh, optimize_with_admesh_truss_arrays

    m = CHILmesh(connectivity=conn.copy(), points=pts.copy(),
                 compute_layers=False, compute_adjacencies=True)
    tri_idx, quad_idx = m._detect_element_types()
    if len(quad_idx) > 0:
        return None  # truss optimizer is triangle-only
    tris = np.asarray(m.connectivity_list)[:, :3].astype(np.int64)
    edge2vert = m.adjacencies["Edge2Vert"]
    boundary_indices = np.unique(edge2vert[m.boundary_edges()].flatten())
    xy = pts[:, :2].astype(float)

    t, _ = _median(
        lambda: optimize_with_admesh_truss_arrays(
            xy.copy(), tris, sdf, boundary_indices=boundary_indices, niter=niter),
        repeats)
    return {"truss_s": t, "niter": niter}


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
    ap.add_argument("--smooth-iters", type=int, default=20,
                    help="iterations for the angle-based smoother lifecycle stage")
    ap.add_argument("--truss-iters", type=int, default=50,
                    help="iterations for the ADMESH truss lifecycle stage")
    ap.add_argument("--sdf", choices=["annulus", "donut", "square", "none"],
                    help="signed-distance function for the truss stage "
                         "(auto-detected from the mesh name if omitted)")
    ap.add_argument("--max-elements", type=int, default=0,
                    help="skip the OOM-prone heavy lifecycle stages (FEM direct "
                         "smoother + ADMESH truss) when the mesh exceeds this "
                         "element count; 0 (default) runs every stage. CI sets a "
                         "finite value so WNAT-scale meshes stay manual-only (#155, #168)")
    args = ap.parse_args()

    if not os.environ.get("CHILMESH_RUN_BENCH") and args.emit is None:
        print("Set CHILMESH_RUN_BENCH=1 (or pass --emit) to run the benchmark.",
              file=sys.stderr)
        # Still run for interactive use; the gate is advisory.

    from chilmesh import CHILmesh
    conn, pts, n_elems_raw, n_verts_raw, grid_name = _load_fort14_arrays(Path(args.mesh))
    n_layers_ref = None

    skip_heavy_bench, bench_gate_reason = _gate_heavy(n_elems_raw,
                                                       args.max_elements)
    if skip_heavy_bench:
        cols = []
    else:
        cols = []
        # Run C++ first (fast), then Rust, then Python (slow), then MATLAB
        print("  [cross-lang] Running C++ bench...", file=sys.stderr, flush=True)
        cpp_r = bench_cpp(conn, pts, args.repeats)
        if cpp_r is not None:
            cols.append(("C++", cpp_r))
            print(f"  [cross-lang] C++ done: full_s={cpp_r['full_s']:.3f}s, n_layers={cpp_r['n_layers']}", file=sys.stderr, flush=True)
        else:
            print("  [cross-lang] C++ not available, skipping.", file=sys.stderr, flush=True)

        print("  [cross-lang] Running Rust bench...", file=sys.stderr, flush=True)
        rust_r = bench_rust(conn, pts, args.repeats)
        if rust_r is not None:
            cols.append(("Rust", rust_r))
            print(f"  [cross-lang] Rust done: full_s={rust_r['full_s']:.3f}s, n_layers={rust_r['n_layers']}", file=sys.stderr, flush=True)
        else:
            print("  [cross-lang] Rust not available, skipping.", file=sys.stderr, flush=True)

        print("  [cross-lang] Running Python bench (may take minutes on large meshes)...", file=sys.stderr, flush=True)
        py_r = bench_python(conn, pts, args.repeats)
        cols.append(("Python", py_r))
        print(f"  [cross-lang] Python done: full_s={py_r['full_s']:.3f}s, n_layers={py_r['n_layers']}", file=sys.stderr, flush=True)

        if args.matlab:
            print("  [cross-lang] Running MATLAB/Octave bench...", file=sys.stderr, flush=True)
            mat_r = bench_matlab(conn, pts)
            if mat_r is not None:
                cols.append(("MATLAB (Octave)", mat_r))
                print(f"  [cross-lang] MATLAB done: full_s={mat_r['full_s']:.3f}s", file=sys.stderr, flush=True)

        # n_layers_ref from C++ if available, else Rust if available, else Python
        n_layers_ref = next((r["n_layers"] for _, r in cols if r is not None), None)

    if skip_heavy_bench:
        lines = [
            f"Mesh: `{Path(args.mesh).name}` — {n_verts_raw:,} vertices · "
            f"{n_elems_raw:,} elements.",
            "",
            f"> **Cross-language bench skipped** — {bench_gate_reason}.",
        ]
        all_match = True
    else:
        lines = [
            f"Mesh: `{Path(args.mesh).name}` — {n_verts_raw:,} vertices · "
            f"{n_elems_raw:,} elements. Median of {args.repeats}.",
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

    # --- Report 2: Python mesh lifecycle (#155) — generation excluded. ---
    skip_heavy, gate_reason = _gate_heavy(n_elems_raw, args.max_elements)
    print("  [lifecycle] Running lifecycle bench...", file=sys.stderr, flush=True)
    life = bench_lifecycle(conn, pts, args.repeats, args.smooth_iters,
                           skip_fem=skip_heavy, skip_skel=skip_heavy)
    sdf_name, sdf = _resolve_sdf(args.mesh, args.sdf)
    truss = (bench_truss(conn, pts, sdf, args.repeats, args.truss_iters)
             if (sdf and not skip_heavy) else None)

    lines += [
        "",
        "### Mesh lifecycle (Python, post-generation)",
        "",
        f"Median of {args.repeats}. Smoother iterations: {life['smooth_iters']}.",
    ]
    if skip_heavy:
        lines.append(f"\n> **Heavy stages skipped** — {gate_reason}.")
    lines += [
        "",
        "| Stage | Time |",
        "|---|---:|",
        f"| Adjacency build | {_fmt_s(life['adj_s'])} |",
        f"| Skeletonization | {_fmt_s(life['skel_s'])} {'_(skipped: --max-elements gate)_' if life.get('skel_skipped') else ''} |",
        f"| Quality (signed area) | {_fmt_s(life['quality_s'])} |",
    ]
    if life["fem_skipped"]:
        lines.append("| FEM direct smoother | — _(skipped: --max-elements gate)_ |")
    else:
        lines.append(f"| FEM direct smoother | {_fmt_s(life['fem_smooth_s'])} |")
    lines.append(f"| Angle-based smoother | {_fmt_s(life['angle_smooth_s'])} |")
    if truss is not None:
        lines.append(
            f"| ADMESH truss optimizer (sdf=`{sdf_name}`, {truss['niter']} iters) "
            f"| {_fmt_s(truss['truss_s'])} |")
    elif skip_heavy:
        lines.append("| ADMESH truss optimizer | — _(skipped: --max-elements gate)_ |")
    else:
        lines.append(f"| ADMESH truss optimizer | — _(skipped: {sdf_name})_ |")

    report = "\n".join(lines) + "\n"
    if args.emit:
        Path(args.emit).write_text(report)
        print(f"wrote {args.emit}")
    else:
        print(report)
    return 0 if all_match else 1


if __name__ == "__main__":
    raise SystemExit(main())
