#!/usr/bin/env python3
"""Benchmark quad-edge (v1) alongside EdgeMap, HE-v1, HE-v2 on WNAT_Hagen.

Measures:
1. Fast init (adjacencies only, no layers)
2. Full init (adjacencies + skeletonization layers)

Records: median-of-3 per operation + peak memory (tracemalloc)
Output: JSON (docs/gallery/benchmark.json) + markdown table (stdout)

Usage:
    python scripts/benchmark_quadegg_variants.py
"""

import json
import time
import sys
import tracemalloc
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
from chilmesh import CHILmesh

BACKENDS = ['edgemap', 'halfedge', 'quadegg']
OPERATIONS = ['fast_init', 'full_init', 'quality_analysis', 'query_latency']
WNAT_HAGEN_PATH = "/tmp/admesh-domains/registry_data/meshes/WNAT_Hagen.14"


def measure_operation(op_name: str, op_fn, n_trials: int = 2) -> Tuple[float, float, float]:
    """Run operation n_trials times, return (median, std, peak_memory).
    
    Uses n_trials=2 (median of 2 is average) to keep runtime reasonable.
    """
    times = []
    peak_memory = 0
    
    for trial in range(n_trials):
        tracemalloc.start()
        start = time.perf_counter()
        try:
            op_fn()
        finally:
            elapsed = time.perf_counter() - start
            current, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
        
        times.append(elapsed)
        peak_memory = max(peak_memory, peak)
    
    times_sorted = sorted(times)
    median = times_sorted[len(times) // 2]
    std = np.std(times) if len(times) > 1 else 0.0
    return median, std, peak_memory


def benchmark_backend(backend_name: str, mesh_path: str) -> Dict[str, Tuple[float, float, float]]:
    """Benchmark all operations for one backend.
    
    Returns: {operation: (median_sec, std_sec, peak_memory_bytes)}
    """
    results = {}
    
    print(f"  {backend_name}...", file=sys.stderr, end=" ", flush=True)
    
    # Operation 1: Fast init (adjacencies only)
    def fast_init():
        mesh = CHILmesh.read_from_fort14(
            mesh_path,
            compute_layers=False,
            compute_adjacencies=True,
            topology_backend=backend_name
        )
    
    median, std, peak_mem = measure_operation('fast_init', fast_init, n_trials=2)
    results['fast_init'] = (median, std, peak_mem)
    
    # Operation 2: Full init (adjacencies + layers)
    def full_init():
        mesh = CHILmesh.read_from_fort14(
            mesh_path,
            compute_layers=True,
            topology_backend=backend_name
        )
    
    median, std, peak_mem = measure_operation('full_init', full_init, n_trials=2)
    results['full_init'] = (median, std, peak_mem)

    # Load mesh for quality analysis and query latency
    mesh = CHILmesh.read_from_fort14(
        mesh_path,
        compute_layers=True,
        topology_backend=backend_name
    )

    # Operation 3: Quality analysis (compute signed areas)
    def quality_analysis():
        _ = mesh.signed_area()

    median, std, peak_mem = measure_operation('quality_analysis', quality_analysis, n_trials=2)
    results['quality_analysis'] = (median, std, peak_mem)

    # Operation 4: Query latency (vertex edge lookup)
    def query_latency():
        for v_idx in range(0, mesh.n_verts, max(1, mesh.n_verts // 10)):
            _ = mesh.get_vertex_edges(v_idx)

    median, std, peak_mem = measure_operation('query_latency', query_latency, n_trials=2)
    results['query_latency'] = (median, std, peak_mem)

    print(f"✓", file=sys.stderr)

    return results


def main():
    # Verify mesh file exists
    if not Path(WNAT_HAGEN_PATH).exists():
        print(f"ERROR: Mesh not found: {WNAT_HAGEN_PATH}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Benchmarking {len(BACKENDS)} backends on WNAT_Hagen", file=sys.stderr)
    print(f"Mesh: {WNAT_HAGEN_PATH}", file=sys.stderr)
    
    # Run benchmarks for all backends
    all_results = {}
    for backend in BACKENDS:
        try:
            all_results[backend] = benchmark_backend(backend, WNAT_HAGEN_PATH)
        except Exception as e:
            print(f"\nERROR benchmarking {backend}: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
            continue
    
    print(f"", file=sys.stderr)  # newline
    
    # Write JSON
    output_data = {
        'metadata': {
            'mesh': 'WNAT_Hagen',
            'mesh_path': WNAT_HAGEN_PATH,
            'backends': BACKENDS,
            'operations': OPERATIONS,
        },
        'results': {
            backend: {
                op: {
                    'median_sec': results[op][0],
                    'std_sec': results[op][1],
                    'peak_memory_bytes': results[op][2]
                }
                for op, results in all_results[backend].items()
            }
            for backend in all_results
        }
    }
    
    output_path = Path('docs/gallery/benchmark.json')
    output_path.parent.mkdir(exist_ok=True, parents=True)
    
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"Benchmark results written to {output_path}")
    
    # Print markdown table
    print("\n## Benchmark Results: WNAT_Hagen (3 Backends)\n")
    print("| Operation | EdgeMap | Half-Edge v1 | Quad-Edge |")
    print("|-----------|---------|--------------|-----------|")
    
    for op in OPERATIONS:
        row = [op]
        for backend in BACKENDS:
            if backend in all_results and op in all_results[backend]:
                median, std, peak_mem = all_results[backend][op]
                row.append(f"{median:.2f}s")
            else:
                row.append("N/A")
        print("| " + " | ".join(row) + " |")
    
    print("\n## Peak Memory Usage (MB)\n")
    print("| Operation | EdgeMap | Half-Edge v1 | Quad-Edge |")
    print("|-----------|---------|--------------|-----------|")
    
    for op in OPERATIONS:
        row = [op]
        for backend in BACKENDS:
            if backend in all_results and op in all_results[backend]:
                median, std, peak_mem = all_results[backend][op]
                row.append(f"{peak_mem/1e6:.0f}MB")
            else:
                row.append("N/A")
        print("| " + " | ".join(row) + " |")


if __name__ == '__main__':
    main()
