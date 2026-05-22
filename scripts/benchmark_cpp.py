#!/usr/bin/env python3
"""Benchmark CHILmesh backends: Python EdgeMap, Rust, and C++ half-edge.

Usage:
    python scripts/benchmark_cpp.py [--trials N] [--json OUTPUT.json]

Runs on annulus, donut, structured fixtures (and WNAT_Hagen if available).
Reports mean ± std dev for fast_init, full_init, quality_analysis, query_latency.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

import numpy as np
import chilmesh


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fmt_ms(s: float) -> str:
    if s >= 1.0:
        return f"{s:.3f} s"
    if s >= 1e-3:
        return f"{s*1000:.1f} ms"
    return f"{s*1e6:.1f} μs"


def _mean_std(times: list[float]) -> tuple[float, float]:
    if len(times) == 1:
        return times[0], 0.0
    return statistics.mean(times), statistics.stdev(times)


def _ratio(val: float, ref: float) -> str:
    if ref <= 0:
        return "—"
    return f"{val/ref:.2f}×"


# ---------------------------------------------------------------------------
# Python (EdgeMap) backend
# ---------------------------------------------------------------------------

def bench_python(pts, conn, fixture_path, n_trials):
    results = {}

    # fast_init
    times = []
    for _ in range(n_trials):
        t0 = time.perf_counter()
        m = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=False)
        times.append(time.perf_counter() - t0)
    results["fast_init"] = _mean_std(times)

    # full_init
    times = []
    for _ in range(n_trials):
        t0 = time.perf_counter()
        m = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=True)
        times.append(time.perf_counter() - t0)
    results["full_init"] = _mean_std(times)
    m_full = m

    # quality_analysis
    times = []
    for _ in range(n_trials):
        t0 = time.perf_counter()
        m_full.elem_quality()
        times.append(time.perf_counter() - t0)
    results["quality"] = _mean_std(times)

    # query latency: vert2edge lookups
    v2e = m_full.adjacencies.get("Vert2Edge", {})
    n_q = min(5000, m_full.n_verts)
    times = []
    for _ in range(n_trials):
        t0 = time.perf_counter()
        for v in range(n_q):
            _ = v2e[v]
        times.append((time.perf_counter() - t0) / n_q * 1e6)
    results["query_us"] = _mean_std(times)

    return results


# ---------------------------------------------------------------------------
# Rust backend
# ---------------------------------------------------------------------------

def bench_rust(pts, conn, fixture_path, n_trials):
    try:
        import chilmesh_core as _rust
    except ImportError:
        return None

    results = {}
    conn32 = conn.astype(np.int32)
    pts64 = pts.astype(np.float64)

    # fast_init: just read + no adjacency
    times = []
    for _ in range(n_trials):
        t0 = time.perf_counter()
        m = _rust.RustMesh()
        m.read_from_fort14(str(fixture_path))
        times.append(time.perf_counter() - t0)
    results["fast_init"] = _mean_std(times)

    # full_init: read + build adjacency
    times = []
    for _ in range(n_trials):
        t0 = time.perf_counter()
        m = _rust.RustMesh()
        m.read_from_fort14(str(fixture_path))
        m.build_adjacencies()
        times.append(time.perf_counter() - t0)
    results["full_init"] = _mean_std(times)
    m_full_rust = m

    # quality_analysis (Rust doesn't expose this directly; time Python fallback on Rust data)
    e2v = m_full_rust.get_edge2vert()
    elem2v = m_full_rust.get_elem2vert()
    times = []
    for _ in range(n_trials):
        t0 = time.perf_counter()
        n_el = m_full_rust.n_elems
        # Compute signed area from elem2vert
        v0 = pts64[elem2v[:, 0]]
        v1 = pts64[elem2v[:, 1]]
        v2 = pts64[elem2v[:, 2]]
        _areas = 0.5 * ((v1[:,0]-v0[:,0])*(v2[:,1]-v0[:,1]) - (v2[:,0]-v0[:,0])*(v1[:,1]-v0[:,1]))
        times.append(time.perf_counter() - t0)
    results["quality"] = _mean_std(times)

    # query latency
    v2e = m_full_rust.get_vert2edge()
    n_q = min(5000, m_full_rust.n_verts)
    times = []
    for _ in range(n_trials):
        t0 = time.perf_counter()
        for v in range(n_q):
            _ = v2e[v]
        times.append((time.perf_counter() - t0) / n_q * 1e6)
    results["query_us"] = _mean_std(times)

    return results


# ---------------------------------------------------------------------------
# C++ backend
# ---------------------------------------------------------------------------

def bench_cpp(pts, conn, fixture_path, n_trials):
    try:
        import chilmesh_cpp as _cpp
    except ImportError:
        return None

    results = {}
    pts64 = np.asarray(pts[:, :2], dtype=np.float64, order='C')
    conn32 = np.asarray(conn, dtype=np.int32, order='C')

    # fast_init
    times = []
    for _ in range(n_trials):
        t0 = time.perf_counter()
        m = _cpp.fast_init(pts64, conn32)
        times.append(time.perf_counter() - t0)
    results["fast_init"] = _mean_std(times)

    # full_init
    times = []
    for _ in range(n_trials):
        t0 = time.perf_counter()
        m = _cpp.full_init(pts64, conn32)
        times.append(time.perf_counter() - t0)
    results["full_init"] = _mean_std(times)
    m_full = m

    # quality_analysis
    times = []
    for _ in range(n_trials):
        t0 = time.perf_counter()
        _areas = _cpp.quality_analysis(m_full)
        times.append(time.perf_counter() - t0)
    results["quality"] = _mean_std(times)

    # query latency
    n_q = min(5000, m_full.n_verts)
    times = []
    for _ in range(n_trials):
        t0 = time.perf_counter()
        for v in range(n_q):
            _ = _cpp.get_vertex_edges(m_full, v)
        times.append((time.perf_counter() - t0) / n_q * 1e6)
    results["query_us"] = _mean_std(times)

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_fixture(name: str, fixture_fn, n_trials: int) -> dict:
    print(f"\n{'='*60}")
    print(f"Fixture: {name}")
    print(f"{'='*60}")

    # Load reference mesh to get points/connectivity + fort14 path
    fpath = chilmesh.examples.fixture_path(
        {"annulus": "annulus_200pts.fort.14",
         "donut": "donut_domain.fort.14",
         "structured": "structuredMesh1.14"}[name]
    )
    ref = fixture_fn()
    pts = ref.points[:, :2].astype(np.float64)
    conn = ref.connectivity_list.astype(np.int32)
    print(f"  n_verts={ref.n_verts}  n_elems={ref.n_elems}")

    print(f"\n  Running Python (EdgeMap) [{n_trials} trials]...")
    py_res = bench_python(pts, conn, fpath, n_trials)

    print(f"  Running Rust [{n_trials} trials]...")
    rust_res = bench_rust(pts, conn, fpath, n_trials)

    print(f"  Running C++ [{n_trials} trials]...")
    cpp_res = bench_cpp(pts, conn, fpath, n_trials)

    return {"name": name, "n_verts": ref.n_verts, "n_elems": ref.n_elems,
            "python": py_res, "rust": rust_res, "cpp": cpp_res}


def print_table(results: list[dict]) -> None:
    print("\n" + "="*80)
    print("BENCHMARK RESULTS")
    print("="*80)

    ops = [("fast_init", "Fast init"), ("full_init", "Full init"),
           ("quality", "Quality"), ("query_us", "Query (μs/call)")]

    for res in results:
        name = res["name"]
        print(f"\n### {name} (n_verts={res['n_verts']}, n_elems={res['n_elems']})")
        print(f"{'Operation':<25} {'Python':>15} {'Rust':>15} {'C++':>15} {'C++/Py ratio':>14}")
        print("-"*85)
        py = res["python"]
        rust = res.get("rust") or {}
        cpp = res.get("cpp") or {}

        for key, label in ops:
            def fmt(d, k):
                if not d or k not in d:
                    return "N/A"
                m, s = d[k]
                if k == "query_us":
                    return f"{m:.3f}±{s:.3f}μs"
                return f"{_fmt_ms(m)}"

            py_val = py[key][0] if key in py else None
            cpp_val = cpp[key][0] if key in cpp else None
            ratio = _ratio(cpp_val, py_val) if py_val and cpp_val else "—"

            print(f"  {label:<23} {fmt(py,key):>15} {fmt(rust,key):>15} {fmt(cpp,key):>15} {ratio:>14}")


def main():
    parser = argparse.ArgumentParser(description="Benchmark CHILmesh backends")
    parser.add_argument("--trials", type=int, default=3, help="Number of trials per operation")
    parser.add_argument("--json", type=str, default=None, help="Save JSON results")
    args = parser.parse_args()

    fixtures = [
        ("annulus", chilmesh.examples.annulus),
        ("donut", chilmesh.examples.donut),
        ("structured", chilmesh.examples.structured),
    ]

    all_results = []
    for name, fn in fixtures:
        try:
            r = run_fixture(name, fn, args.trials)
            all_results.append(r)
        except Exception as e:
            print(f"  ERROR on {name}: {e}")

    print_table(all_results)

    if args.json:
        # Convert tuples to lists for JSON serialization
        def _ser(obj):
            if isinstance(obj, tuple):
                return list(obj)
            raise TypeError(f"Not serializable: {type(obj)}")
        with open(args.json, "w") as f:
            json.dump(all_results, f, default=_ser, indent=2)
        print(f"\nJSON saved to {args.json}")

    return all_results


if __name__ == "__main__":
    main()
