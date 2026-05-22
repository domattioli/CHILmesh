#!/usr/bin/env python3
"""Benchmark Rust backend vs EdgeMap Python backend on WNAT_Hagen and built-in fixtures.

Measures:
1. fast_init (parse + adjacency, no layers)
2. full_init (parse + adjacency + skeletonization)
3. quality_analysis (signed_area computation)
4. query_latency (get_vertex_edges for sampled vertices)

Protocol: 3 trials per operation; reports mean ± std dev.

Usage:
    python scripts/benchmark_rust.py
    python scripts/benchmark_rust.py --fixture annulus
    python scripts/benchmark_rust.py --no-wnat  # skip WNAT_Hagen (uses built-in fixtures only)
"""

import argparse
import sys
import time
import tracemalloc
from pathlib import Path
from typing import Dict, Tuple, Optional

import numpy as np

WNAT_HAGEN_PATH = "/tmp/admesh-domains/registry_data/meshes/WNAT_Hagen.14"
N_TRIALS = 3


def measure_operation(op_fn, n_trials: int = N_TRIALS) -> Tuple[float, float, float]:
    """Run op_fn n_trials times; return (mean_sec, std_sec, peak_memory_bytes)."""
    times = []
    peak_memory = 0
    for _ in range(n_trials):
        tracemalloc.start()
        t0 = time.perf_counter()
        try:
            op_fn()
        finally:
            elapsed = time.perf_counter() - t0
            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
        times.append(elapsed)
        peak_memory = max(peak_memory, peak)
    arr = np.array(times)
    return float(arr.mean()), float(arr.std()), peak_memory


def benchmark_edgemap_on_path(mesh_path: str) -> Dict[str, Tuple[float, float, float]]:
    """Benchmark EdgeMap (Python) backend using fort14 file path."""
    from chilmesh import CHILmesh

    results: Dict[str, Tuple[float, float, float]] = {}

    # fast_init
    def fast_init():
        CHILmesh.read_from_fort14(mesh_path, compute_layers=False, compute_adjacencies=True,
                                  topology_backend='edgemap')

    results['fast_init'] = measure_operation(fast_init)

    # full_init
    def full_init():
        CHILmesh.read_from_fort14(mesh_path, compute_layers=True, topology_backend='edgemap')

    results['full_init'] = measure_operation(full_init)

    # Load once for quality + query
    mesh = CHILmesh.read_from_fort14(mesh_path, compute_layers=True, topology_backend='edgemap')

    # quality_analysis (use already-loaded mesh)
    def quality_analysis():
        _ = mesh.signed_area()

    results['quality_analysis'] = measure_operation(quality_analysis)

    # query_latency
    sample_verts = list(range(0, mesh.n_verts, max(1, mesh.n_verts // 100)))

    def query_latency():
        for v in sample_verts:
            _ = mesh.get_vertex_edges(v)

    results['query_latency'] = measure_operation(query_latency)

    return results


def benchmark_rust_on_path(mesh_path: str) -> Dict[str, Tuple[float, float, float]]:
    """Benchmark Rust backend using fort14 file path."""
    from chilmesh import CHILmesh

    results: Dict[str, Tuple[float, float, float]] = {}

    # fast_init: parse + Rust adjacency
    def fast_init():
        CHILmesh(
            **_load_fort14_raw(mesh_path),
            compute_layers=False,
            compute_adjacencies=True,
            use_rust_backend=True,
        )

    results['fast_init'] = measure_operation(fast_init)

    # full_init
    def full_init():
        CHILmesh(
            **_load_fort14_raw(mesh_path),
            compute_layers=True,
            use_rust_backend=True,
        )

    results['full_init'] = measure_operation(full_init)

    # Load once for quality + query
    raw = _load_fort14_raw(mesh_path)
    mesh_rust = CHILmesh(**raw, compute_layers=True, use_rust_backend=True)

    # quality_analysis (use already-loaded mesh)
    def quality_analysis():
        _ = mesh_rust.signed_area()

    results['quality_analysis'] = measure_operation(quality_analysis)

    # query_latency
    sample_verts = list(range(0, mesh_rust.n_verts, max(1, mesh_rust.n_verts // 100)))

    def query_latency():
        for v in sample_verts:
            _ = mesh_rust.get_vertex_edges(v)

    results['query_latency'] = measure_operation(query_latency)

    return results


def _load_fort14_raw(mesh_path: str) -> dict:
    """Parse fort14 file and return dict of {connectivity, points, grid_name}."""
    from chilmesh import CHILmesh
    # Use edgemap (fast, no Rust needed) to parse
    # Extract raw arrays from parsed mesh (avoids duplicating parser logic)
    m = CHILmesh.read_from_fort14(mesh_path, compute_layers=False,
                                   compute_adjacencies=False, topology_backend='edgemap')
    return {
        'connectivity': m.connectivity_list.copy(),
        'points': m.points.copy(),
        'grid_name': m.grid_name,
    }


def benchmark_edgemap_fixture(fixture_name: str) -> Dict[str, Tuple[float, float, float]]:
    """Benchmark EdgeMap on a built-in example fixture."""
    from chilmesh import examples, CHILmesh

    results: Dict[str, Tuple[float, float, float]] = {}
    loader = getattr(examples, fixture_name)
    raw_mesh = loader(compute_layers=False, compute_adjacencies=False)

    def fast_init():
        CHILmesh(connectivity=raw_mesh.connectivity_list.copy(),
                 points=raw_mesh.points.copy(),
                 grid_name=raw_mesh.grid_name,
                 compute_layers=False, compute_adjacencies=True,
                 topology_backend='edgemap')

    results['fast_init'] = measure_operation(fast_init)

    def full_init():
        CHILmesh(connectivity=raw_mesh.connectivity_list.copy(),
                 points=raw_mesh.points.copy(),
                 grid_name=raw_mesh.grid_name,
                 compute_layers=True, topology_backend='edgemap')

    results['full_init'] = measure_operation(full_init)

    mesh = loader(compute_layers=True)

    def quality_analysis():
        _ = mesh.signed_area()

    results['quality_analysis'] = measure_operation(quality_analysis)

    sample_verts = list(range(0, mesh.n_verts, max(1, mesh.n_verts // 100)))

    def query_latency():
        for v in sample_verts:
            _ = mesh.get_vertex_edges(v)

    results['query_latency'] = measure_operation(query_latency)

    return results


def benchmark_rust_fixture(fixture_name: str) -> Dict[str, Tuple[float, float, float]]:
    """Benchmark Rust backend on a built-in example fixture."""
    from chilmesh import examples, CHILmesh

    results: Dict[str, Tuple[float, float, float]] = {}
    loader = getattr(examples, fixture_name)
    raw_mesh = loader(compute_layers=False, compute_adjacencies=False)

    def fast_init():
        CHILmesh(connectivity=raw_mesh.connectivity_list.copy(),
                 points=raw_mesh.points.copy(),
                 grid_name=raw_mesh.grid_name,
                 compute_layers=False, compute_adjacencies=True,
                 use_rust_backend=True)

    results['fast_init'] = measure_operation(fast_init)

    def full_init():
        CHILmesh(connectivity=raw_mesh.connectivity_list.copy(),
                 points=raw_mesh.points.copy(),
                 grid_name=raw_mesh.grid_name,
                 compute_layers=True, use_rust_backend=True)

    results['full_init'] = measure_operation(full_init)

    mesh = CHILmesh(connectivity=raw_mesh.connectivity_list.copy(),
                    points=raw_mesh.points.copy(),
                    grid_name=raw_mesh.grid_name,
                    compute_layers=True, use_rust_backend=True)

    def quality_analysis():
        _ = mesh.signed_area()

    results['quality_analysis'] = measure_operation(quality_analysis)

    sample_verts = list(range(0, mesh.n_verts, max(1, mesh.n_verts // 100)))

    def query_latency():
        for v in sample_verts:
            _ = mesh.get_vertex_edges(v)

    results['query_latency'] = measure_operation(query_latency)

    return results


OPERATIONS = ['fast_init', 'full_init', 'quality_analysis', 'query_latency']


def _fmt(median: float, std: float) -> str:
    if median >= 1.0:
        return f"{median:.3f}s ± {std:.3f}s"
    elif median >= 0.001:
        return f"{median*1000:.1f}ms ± {std*1000:.1f}ms"
    else:
        return f"{median*1e6:.1f}µs ± {std*1e6:.1f}µs"


def _ratio(rust_t: float, base_t: float) -> str:
    if base_t == 0:
        return "N/A"
    r = rust_t / base_t
    return f"{r:.2f}×"


def print_table(title: str, edgemap_r: dict, rust_r: dict) -> None:
    print(f"\n### {title}\n")
    print(f"| Operation       | EdgeMap (baseline)     | Rust backend           | Ratio (Rust/EdgeMap) |")
    print(f"|-----------------|------------------------|------------------------|----------------------|")
    for op in OPERATIONS:
        if op in edgemap_r and op in rust_r:
            em_t, em_s, em_m = edgemap_r[op]
            rs_t, rs_s, rs_m = rust_r[op]
            ratio = _ratio(rs_t, em_t)
            print(f"| {op:<15} | {_fmt(em_t, em_s):<22} | {_fmt(rs_t, rs_s):<22} | {ratio:<20} |")
    print()


def main():
    parser = argparse.ArgumentParser(description="Benchmark Rust vs EdgeMap backend")
    parser.add_argument('--fixture', default=None, help="Built-in fixture name (annulus/donut/structured)")
    parser.add_argument('--no-wnat', action='store_true', help="Skip WNAT_Hagen benchmark")
    args = parser.parse_args()

    results_all = {}

    # WNAT_Hagen benchmark
    if not args.no_wnat and not args.fixture:
        wnat_path = WNAT_HAGEN_PATH
        if not Path(wnat_path).exists():
            print(f"WARNING: WNAT_Hagen not found at {wnat_path}, skipping.", file=sys.stderr)
        else:
            print(f"Benchmarking on WNAT_Hagen ({N_TRIALS} trials each)...", file=sys.stderr)
            print("  EdgeMap...", file=sys.stderr, end=" ", flush=True)
            em_r = benchmark_edgemap_on_path(wnat_path)
            print("done", file=sys.stderr)
            print("  Rust...", file=sys.stderr, end=" ", flush=True)
            rust_r = benchmark_rust_on_path(wnat_path)
            print("done", file=sys.stderr)
            results_all['WNAT_Hagen'] = {'edgemap': em_r, 'rust': rust_r}
            print_table("WNAT_Hagen", em_r, rust_r)

            # Performance gate check
            full_init_rust = rust_r['full_init'][0]
            TARGET = 3.5
            print(f"## Performance Gate: full_init WNAT_Hagen")
            if full_init_rust <= TARGET:
                print(f"  PASS: {full_init_rust:.3f}s ≤ {TARGET}s target")
            elif full_init_rust <= TARGET * 1.2:
                print(f"  WARN: {full_init_rust:.3f}s > {TARGET}s target but ≤ {TARGET*1.2:.1f}s (1.2× threshold)")
            else:
                print(f"  FAIL: {full_init_rust:.3f}s > {TARGET*1.2:.1f}s (1.2× threshold). Profiling recommended.")

    # Built-in fixture benchmarks
    fixture_names = [args.fixture] if args.fixture else ['annulus', 'donut', 'structured']
    for fname in fixture_names:
        print(f"\nBenchmarking on {fname} ({N_TRIALS} trials each)...", file=sys.stderr)
        print("  EdgeMap...", file=sys.stderr, end=" ", flush=True)
        try:
            em_r = benchmark_edgemap_fixture(fname)
            print("done", file=sys.stderr)
            print("  Rust...", file=sys.stderr, end=" ", flush=True)
            rust_r = benchmark_rust_fixture(fname)
            print("done", file=sys.stderr)
            results_all[fname] = {'edgemap': em_r, 'rust': rust_r}
            print_table(fname, em_r, rust_r)
        except Exception as e:
            print(f"ERROR: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()

    # Summary
    if 'WNAT_Hagen' in results_all:
        em_fi = results_all['WNAT_Hagen']['edgemap']['full_init'][0]
        ru_fi = results_all['WNAT_Hagen']['rust']['full_init'][0]
        em_mi = results_all['WNAT_Hagen']['edgemap']['full_init'][2] / 1e6
        ru_mi = results_all['WNAT_Hagen']['rust']['full_init'][2] / 1e6
        print("\n## Summary (WNAT_Hagen)\n")
        print(f"  EdgeMap full_init: {em_fi:.3f}s, peak memory: {em_mi:.0f} MB")
        print(f"  Rust full_init:    {ru_fi:.3f}s, peak memory: {ru_mi:.0f} MB")
        print(f"  Rust/EdgeMap ratio: {ru_fi/em_fi:.2f}×")
        mem_ratio = ru_mi / em_mi if em_mi > 0 else 0
        print(f"  Memory ratio: {mem_ratio:.2f}× ({'OK' if mem_ratio <= 1.25 else 'OVER BUDGET'})")

    return results_all


if __name__ == '__main__':
    main()
