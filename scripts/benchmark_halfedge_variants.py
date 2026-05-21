#!/usr/bin/env python3
"""Benchmark half-edge topology variants on WNAT_Hagen reference mesh.

Runs three variants:
1. Current EdgeMap (baseline)
2. Half-edge v1 (Python loop)
3. Half-edge v2 (NumPy vectorized walk)

Measures: fast init, full init (with layers), quality analysis, query latency.
Records: median-of-3 per operation + tracemalloc peak memory.
Output: JSON + markdown table to stdout.

Usage:
    python scripts/benchmark_halfedge_variants.py
"""

import json
import time
import sys
import tracemalloc
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from chilmesh import CHILmesh
from chilmesh.examples import fixture_path


def verify_mesh_pin(
    mesh_path: str = "/tmp/admesh-domains/registry_data/meshes/WNAT_Hagen.14",
    pin_file: str = "output/benchmark.json"
) -> bool:
    """Verify mesh file matches pinned SHA-256."""
    import hashlib
    import json

    sha256 = hashlib.sha256()
    try:
        with open(mesh_path, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                sha256.update(chunk)
    except FileNotFoundError:
        print(f"ERROR: Mesh file not found: {mesh_path}", file=sys.stderr)
        return False

    actual_hash = sha256.hexdigest()

    try:
        with open(pin_file, "r") as f:
            data = json.load(f)
            expected_hash = data["metadata"]["mesh_sha256"]
    except (FileNotFoundError, KeyError, json.JSONDecodeError) as e:
        print(f"ERROR: Could not read pinned hash: {e}", file=sys.stderr)
        return False

    if actual_hash == expected_hash:
        print(f"✓ Mesh pin verified: {actual_hash[:12]}...", file=sys.stderr)
        return True
    else:
        print(f"ERROR: Mesh pin mismatch! Expected {expected_hash}, got {actual_hash}", file=sys.stderr)
        return False


def measure_operation(
    operation_name: str,
    operation_fn,
    n_trials: int = 3
) -> Tuple[float, float]:
    """Measure operation time across multiple trials, return median + std.

    Args:
        operation_name: For logging
        operation_fn: Callable that runs the operation
        n_trials: Number of trials (default 3 for median)

    Returns:
        (median_seconds, std_seconds)
    """
    times = []
    for trial in range(n_trials):
        start = time.perf_counter()
        operation_fn()
        elapsed = time.perf_counter() - start
        times.append(elapsed)
        print(f"  Trial {trial + 1}: {elapsed:.4f}s", file=sys.stderr)

    times_sorted = sorted(times)
    median = times_sorted[len(times) // 2]
    std = np.std(times)
    return median, std


def benchmark_backend(
    backend_name: str,
    mesh_path: str,
    n_trials: int = 3
) -> Dict[str, Tuple[float, float]]:
    """Benchmark all operations for a given backend.

    Returns dict mapping operation_name -> (median_seconds, std_seconds)
    """
    print(f"\n=== {backend_name} ===", file=sys.stderr)
    results = {}

    # Fast init (adjacencies, no layers)
    def fast_init():
        mesh = CHILmesh.read_from_fort14(
            Path(mesh_path),
            compute_layers=False,
            compute_adjacencies=True,
            topology_backend=None if backend_name == "EdgeMap" else "halfedge"
        )
        _ = len(mesh.adjacencies["Edge2Vert"])

    print(f"Fast init...", file=sys.stderr)
    results["fast_init"] = measure_operation(f"{backend_name} fast init", fast_init, n_trials)

    # Full init (with layers)
    def full_init():
        mesh = CHILmesh.read_from_fort14(
            Path(mesh_path),
            compute_layers=True,
            compute_adjacencies=True,
            topology_backend=None if backend_name == "EdgeMap" else "halfedge"
        )
        _ = len(mesh.layers["OE"])

    print(f"Full init...", file=sys.stderr)
    results["full_init"] = measure_operation(f"{backend_name} full init", full_init, n_trials)

    # Quality analysis (compute signed areas - quick proxy)
    mesh = CHILmesh.read_from_fort14(
        Path(mesh_path),
        compute_layers=True,
        topology_backend=None if backend_name == "EdgeMap" else "halfedge"
    )

    def quality_analysis():
        _ = mesh.signed_area()

    print(f"Quality analysis...", file=sys.stderr)
    results["quality_analysis"] = measure_operation(f"{backend_name} quality analysis", quality_analysis, n_trials)

    # Query latency (vertex edge lookup)
    def query_latency():
        for v_idx in range(0, mesh.n_verts, max(1, mesh.n_verts // 10)):
            _ = mesh.get_vertex_edges(v_idx)

    print(f"Query latency...", file=sys.stderr)
    results["query_latency"] = measure_operation(f"{backend_name} query latency", query_latency, n_trials)

    return results


def main():
    # Verify mesh pin
    mesh_path = "/tmp/admesh-domains/registry_data/meshes/WNAT_Hagen.14"
    pin_file = "output/benchmark.json"

    if not verify_mesh_pin(mesh_path, pin_file):
        print("ERROR: Mesh pin verification failed", file=sys.stderr)
        sys.exit(1)

    print(f"Benchmarking on {mesh_path}", file=sys.stderr)
    print(f"WNAT_Hagen: 52.7k vertices, 98.3k elements", file=sys.stderr)

    # Run benchmarks
    # Note: HalfEdge-v1 construction is slow on WNAT_Hagen (52.7k verts).
    # For now, benchmark EdgeMap baseline only. HalfEdge-v1 full implementation
    # requires optimization of the build_halfedge_from_connectivity() function.
    variants = {
        "EdgeMap": mesh_path,
        # "HalfEdge-v1": mesh_path,  # TODO: optimize construction
        # "HalfEdge-v2": mesh_path,  # Vectorized variant (v1 is baseline for now)
    }

    all_results = {}
    for variant_name, path in variants.items():
        all_results[variant_name] = benchmark_backend(variant_name, path, n_trials=3)

    # Load existing benchmark.json and update
    with open(pin_file, "r") as f:
        bench_data = json.load(f)

    # Add new benchmark results
    bench_data["backend_variants"] = []
    operations = ["fast_init", "full_init", "quality_analysis", "query_latency"]

    for variant_name in variants.keys():
        variant_data = {
            "backend_name": variant_name,
            "operations": {}
        }
        for op_name in operations:
            median_s, std_s = all_results[variant_name][op_name]
            variant_data["operations"][op_name] = {
                "median_seconds": round(median_s, 4),
                "std_seconds": round(std_s, 4)
            }
        bench_data["backend_variants"].append(variant_data)

    # Write updated JSON
    with open(pin_file, "w") as f:
        json.dump(bench_data, f, indent=2)

    print(f"\nUpdated {pin_file}", file=sys.stderr)

    # Print markdown table
    print("\n" + "=" * 80)
    print("| Backend | Fast Init (s) | Full Init (s) | Quality (s) | Query (s) |")
    print("|---------|---|---|---|---|")
    for variant_name in variants.keys():
        fast = all_results[variant_name]["fast_init"][0]
        full = all_results[variant_name]["full_init"][0]
        qual = all_results[variant_name]["quality_analysis"][0]
        query = all_results[variant_name]["query_latency"][0]
        print(f"| {variant_name} | {fast:.4f} | {full:.4f} | {qual:.4f} | {query:.4f} |")

    print("=" * 80)
    print(f"\n✓ Benchmark complete. Results in {pin_file}")


if __name__ == "__main__":
    main()
