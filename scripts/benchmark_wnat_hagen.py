#!/usr/bin/env python3
"""Benchmark CHILmesh on the WNAT_Hagen reference mesh.

Usage:
    python scripts/benchmark_wnat_hagen.py [MESH_PATH] [--json OUTPUT.json]

If MESH_PATH is omitted the script looks for the mesh in the default
ADMESH-Domains clone location (/tmp/admesh-domains/registry_data/meshes/).

Outputs a markdown table suitable for pasting into docs/BENCHMARK.md and,
optionally, a JSON file for CI archival / diff-over-time comparisons.
"""
from __future__ import annotations

import argparse
import json
import platform
import statistics
import sys
import time
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fmt(seconds: float, unit: str = "s") -> str:
    if unit == "s":
        return f"{seconds:.3f}s"
    if unit == "ms":
        return f"{seconds * 1000:.1f}ms"
    if unit == "us":
        return f"{seconds * 1e6:.1f}μs"
    return str(seconds)


def _measure(label: str, fn, *, n: int = 1) -> tuple[str, float]:
    """Run *fn* *n* times and return (label, mean_seconds)."""
    t0 = time.perf_counter()
    for _ in range(n):
        fn()
    elapsed = (time.perf_counter() - t0) / n
    return label, elapsed


def _locate_mesh(explicit: str | None) -> Path:
    candidates = [
        explicit,
        "/tmp/admesh-domains/registry_data/meshes/WNAT_Hagen.14",
        "/tmp/admesh-domains/WNAT_Hagen.14",
        str(Path.home() / "WNAT_Hagen.14"),
    ]
    for c in candidates:
        if c and Path(c).exists():
            return Path(c)
    raise FileNotFoundError(
        "WNAT_Hagen mesh not found. Clone ADMESH-Domains:\n"
        "  git clone https://github.com/domattioli/ADMESH-Domains "
        "/tmp/admesh-domains\n"
        "then re-run this script."
    )


# ---------------------------------------------------------------------------
# Benchmark sections
# ---------------------------------------------------------------------------

def bench_initialization(mesh_path: Path) -> dict:
    from chilmesh import CHILmesh

    results: dict[str, float] = {}

    print("[1/4] Fast init (compute_layers=False) ...", end=" ", flush=True)
    t0 = time.perf_counter()
    mesh_fast = CHILmesh.read_from_fort14(mesh_path, compute_layers=False)
    results["fast_init_s"] = time.perf_counter() - t0
    print(_fmt(results["fast_init_s"]))

    print("[2/4] Full init (compute_layers=True)  ...", end=" ", flush=True)
    t0 = time.perf_counter()
    mesh_full = CHILmesh.read_from_fort14(mesh_path, compute_layers=True)
    results["full_init_s"] = time.perf_counter() - t0
    print(_fmt(results["full_init_s"]))

    print("[3/4] Quality analysis                  ...", end=" ", flush=True)
    t0 = time.perf_counter()
    _q, _a, _s = mesh_full.elem_quality()
    results["quality_s"] = time.perf_counter() - t0
    print(_fmt(results["quality_s"]))

    results["total_s"] = results["full_init_s"] + results["quality_s"]
    results["n_verts"] = int(mesh_full.n_verts)
    results["n_elems"] = int(mesh_full.n_elems)
    results["n_edges"] = int(mesh_full.n_edges)
    results["n_layers"] = int(mesh_full.n_layers)

    return results, mesh_fast, mesh_full


SAMPLE_COUNT = 5_000


def bench_queries(mesh) -> dict:
    from chilmesh import CHILmesh  # noqa: F401 — ensures import path works

    results: dict[str, float] = {}
    n_elems = mesh.n_elems
    n_verts = mesh.n_verts
    sample_elems = min(SAMPLE_COUNT, n_elems)
    sample_verts = min(SAMPLE_COUNT, n_verts)

    # elem2edge
    print(f"[4/4] elem2edge queries ({sample_elems:,} samples)  ...", end=" ", flush=True)
    elem_ids = list(range(sample_elems))
    t0 = time.perf_counter()
    for e in elem_ids:
        mesh.elem2edge(e)
    elapsed = time.perf_counter() - t0
    results["elem2edge_total_s"] = elapsed
    results["elem2edge_us"] = elapsed / sample_elems * 1e6
    print(_fmt(elapsed) + f"  ({results['elem2edge_us']:.1f}μs/query)")

    # Vert2Edge
    print(f"     Vert2Edge lookups  ({sample_verts:,} samples)   ...", end=" ", flush=True)
    vert2edge = mesh.adjacencies["Vert2Edge"]
    t0 = time.perf_counter()
    for v in range(sample_verts):
        _ = vert2edge[v]
    elapsed = time.perf_counter() - t0
    results["vert2edge_total_s"] = elapsed
    results["vert2edge_us"] = elapsed / sample_verts * 1e6
    print(_fmt(elapsed) + f"  ({results['vert2edge_us']:.2f}μs/query)")

    # Elem2Edge bulk (element-neighbor proxy via Elem2Edge ndarray index)
    sample_nn = min(1_000, n_elems)
    print(f"     Elem2Edge bulk     ({sample_nn:,} samples)      ...", end=" ", flush=True)
    elem2edge_arr = mesh.adjacencies["Elem2Edge"]
    t0 = time.perf_counter()
    for e in range(sample_nn):
        _ = elem2edge_arr[e]
    elapsed = time.perf_counter() - t0
    results["elem2edge_bulk_total_s"] = elapsed
    results["elem2edge_bulk_us"] = elapsed / sample_nn * 1e6
    print(_fmt(elapsed) + f"  ({results['elem2edge_bulk_us']:.2f}μs/query)")

    return results


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

BENCHMARK_V011 = {
    "fast_init_s": 3200,
    "full_init_s": 5400,
    "quality_s": 4800,
    "total_s": 13400,
    "elem2edge_us": 2000,
    "vert2edge_us": 3500,
    "elem2edge_bulk_us": 4500,
}


def _speedup(current: float, baseline: float) -> str:
    if current <= 0:
        return "N/A"
    return f"{baseline / current:.0f}x"


def print_report(init: dict, query: dict) -> None:
    print()
    print("=" * 70)
    print("CHILmesh WNAT_Hagen Benchmark Results")
    print(f"  Mesh:    {init['n_verts']:,} vertices, {init['n_elems']:,} elements, {init['n_edges']:,} edges")
    print(f"  Layers:  {init['n_layers']}")
    print(f"  Python:  {platform.python_version()}")
    print(f"  OS:      {platform.system()} {platform.machine()}")
    print("=" * 70)
    print()

    print("### Initialization Performance")
    print()
    rows = [
        ("Fast init (no layers)",  "fast_init_s",  "s"),
        ("Full init (with layers)", "full_init_s",  "s"),
        ("Quality analysis",        "quality_s",    "s"),
        ("Total (full + quality)",  "total_s",      "s"),
    ]
    print(f"{'Operation':<30} {'v0.1.1 (est)':>15} {'Current':>12} {'Speedup':>10}")
    print("-" * 70)
    for label, key, unit in rows:
        cur = init[key]
        est = BENCHMARK_V011.get(key)
        est_str = _fmt(est, unit) if est else "N/A"
        print(f"{label:<30} {est_str:>15} {_fmt(cur, unit):>12} {_speedup(cur, est):>10}")
    print()

    print("### Query Performance")
    print()
    q_rows = [
        ("elem2edge (5k queries)",  "elem2edge_us",       "us"),
        ("Vert2Edge (5k lookups)",  "vert2edge_us",       "us"),
        ("Elem2Edge bulk (1k)",     "elem2edge_bulk_us",  "us"),
    ]
    print(f"{'Operation':<30} {'v0.1.1 (est)':>15} {'Current':>12} {'Speedup':>10}")
    print("-" * 70)
    for label, key, unit in q_rows:
        cur = query[key]
        est = BENCHMARK_V011.get(key)
        est_str = _fmt(est, unit) if est else "N/A"
        print(f"{label:<30} {est_str:>15} {_fmt(cur, unit):>12} {_speedup(cur, est):>10}")
    print()
    print("=" * 70)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark CHILmesh on WNAT_Hagen mesh.")
    parser.add_argument("mesh_path", nargs="?", help="Path to WNAT_Hagen.14 fort.14 file")
    parser.add_argument("--json", dest="json_out", metavar="FILE",
                        help="Write results to JSON file for archival")
    args = parser.parse_args()

    try:
        mesh_path = _locate_mesh(args.mesh_path)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Mesh: {mesh_path}")
    print()

    init_results, mesh_fast, mesh_full = bench_initialization(mesh_path)
    query_results = bench_queries(mesh_fast)

    print_report(init_results, query_results)

    if args.json_out:
        payload = {
            "chilmesh_version": _get_version(),
            "python": platform.python_version(),
            "os": f"{platform.system()} {platform.machine()}",
            "mesh": str(mesh_path),
            "init": {k: v for k, v in init_results.items() if isinstance(v, (int, float))},
            "query": query_results,
        }
        Path(args.json_out).write_text(json.dumps(payload, indent=2))
        print(f"Results written to {args.json_out}")

    return 0


def _get_version() -> str:
    try:
        from importlib.metadata import version
        return version("chilmesh")
    except Exception:
        return "unknown"


if __name__ == "__main__":
    sys.exit(main())
