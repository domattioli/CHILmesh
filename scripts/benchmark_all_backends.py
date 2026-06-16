#!/usr/bin/env python3
"""Benchmark all backends (Python, Rust, C++) on WNAT_Hagen.

.. deprecated::
    Superseded by ``scripts/benchmark.py --matlab`` for cross-language
    comparison. This script is NOT like-for-like: its Python path calls
    ``read_from_fort14`` (includes fort.14 parsing) and ``elem_quality``
    (skew metric), while its C++ path uses in-memory arrays and signed
    area — so Python rows are inflated by file I/O and a heavier quality
    op. Use ``benchmark.py``, which times identical operations on shared
    ``(connectivity, points)`` arrays across all backends. See #216.

Usage:
    python scripts/benchmark_all_backends.py [MESH_PATH]
"""
import json
import platform
import statistics
import sys
import time
from pathlib import Path

def _fmt(seconds: float, unit: str = "s") -> str:
    if unit == "s":
        return f"{seconds:.3f}s"
    if unit == "ms":
        return f"{seconds * 1000:.1f}ms"
    if unit == "us":
        return f"{seconds * 1e6:.1f}μs"
    return str(seconds)

def _locate_mesh(explicit: str | None) -> Path:
    candidates = [
        explicit,
        "/tmp/valence-domains/registry_data/meshes/WNAT_Hagen.14",
        "/tmp/WNAT_Hagen.14",
        str(Path.home() / "WNAT_Hagen.14"),
    ]
    for c in candidates:
        if c and Path(c).exists():
            return Path(c)
    raise FileNotFoundError("WNAT_Hagen mesh not found")

# ============================================================================
# PYTHON Backend
# ============================================================================

def bench_python(mesh_path: Path) -> dict:
    from chilmesh import CHILmesh
    results = {}

    # Fast init
    t0 = time.perf_counter()
    mesh_fast = CHILmesh.read_from_fort14(mesh_path, compute_layers=False)
    results["fast_init_s"] = time.perf_counter() - t0

    # Skeletonization delta (full - fast)
    t0 = time.perf_counter()
    mesh_full = CHILmesh.read_from_fort14(mesh_path, compute_layers=True)
    full_time = time.perf_counter() - t0
    results["skel_only_s"] = full_time - results["fast_init_s"]
    results["full_init_s"] = full_time

    # Quality
    t0 = time.perf_counter()
    mesh_full.elem_quality()
    results["quality_s"] = time.perf_counter() - t0

    # Vertex-edge query (5k samples)
    vert2edge = mesh_full.adjacencies["Vert2Edge"]
    t0 = time.perf_counter()
    for v in range(5000):
        _ = vert2edge[v]
    ve_elapsed = time.perf_counter() - t0
    results["vert2edge_us"] = ve_elapsed / 5000 * 1e6

    results["n_verts"] = int(mesh_full.n_verts)
    results["n_elems"] = int(mesh_full.n_elems)
    results["n_layers"] = int(mesh_full.n_layers)
    return results

# ============================================================================
# RUST Backend
# ============================================================================

def bench_rust(mesh_path: Path) -> dict:
    try:
        from chilmesh_core import RustMesh
    except ImportError:
        return {"error": "RustMesh not available"}

    results = {}

    # Fast init (read only)
    t0 = time.perf_counter()
    mesh = RustMesh()
    mesh.read_from_fort14(str(mesh_path))
    results["fast_init_s"] = time.perf_counter() - t0

    # Full init (read + skeletonize)
    t0 = time.perf_counter()
    mesh2 = RustMesh()
    mesh2.read_from_fort14(str(mesh_path))
    mesh2.skeletonize()
    full_time = time.perf_counter() - t0
    results["skel_only_s"] = full_time - results["fast_init_s"]
    results["full_init_s"] = full_time

    # Quality
    t0 = time.perf_counter()
    mesh2.compute_quality()
    results["quality_s"] = time.perf_counter() - t0

    # Vertex-edge query (5k samples)
    vert2edge = mesh2.get_vert2edge()
    t0 = time.perf_counter()
    for v in range(5000):
        _ = vert2edge[v]
    ve_elapsed = time.perf_counter() - t0
    results["vert2edge_us"] = ve_elapsed / 5000 * 1e6

    # Metadata
    results["n_verts"] = mesh.n_verts
    results["n_elems"] = mesh.n_elems
    results["n_layers"] = mesh2.get_num_layers()

    return results

# ============================================================================
# C++ Backend
# ============================================================================

def bench_cpp(mesh_path: Path) -> dict:
    try:
        import chilmesh_cpp
        import numpy as np
    except ImportError:
        return {"error": "chilmesh_cpp not available"}

    # Load mesh with Python first
    from chilmesh import CHILmesh
    py_mesh = CHILmesh.read_from_fort14(mesh_path, compute_layers=False)
    points = py_mesh.points[:, :2].astype(np.float64)  # Extract x, y only
    connectivity = np.array([list(c[:3]) for c in py_mesh.connectivity_list], dtype=np.int32)

    results = {}

    # Fast init (benchmark C++ fast_init on arrays)
    t0 = time.perf_counter()
    cpp_mesh_fast = chilmesh_cpp.fast_init(points, connectivity)
    results["fast_init_s"] = time.perf_counter() - t0

    # Full init
    t0 = time.perf_counter()
    cpp_mesh_full = chilmesh_cpp.full_init(points, connectivity)
    full_time = time.perf_counter() - t0
    results["skel_only_s"] = full_time - results["fast_init_s"]
    results["full_init_s"] = full_time

    # Quality
    t0 = time.perf_counter()
    chilmesh_cpp.quality_analysis(cpp_mesh_full)
    results["quality_s"] = time.perf_counter() - t0

    # Vertex-edge query (5k samples)
    vert2edge = cpp_mesh_full.vert2edge
    t0 = time.perf_counter()
    for v in range(5000):
        _ = vert2edge[v]
    ve_elapsed = time.perf_counter() - t0
    results["vert2edge_us"] = ve_elapsed / 5000 * 1e6

    results["n_verts"] = cpp_mesh_full.n_verts
    results["n_elems"] = cpp_mesh_full.n_elems
    results["n_layers"] = cpp_mesh_full.n_layers

    return results

# ============================================================================
# Report
# ============================================================================

def print_report(python_results: dict, rust_results: dict, cpp_results: dict) -> None:
    print()
    print("=" * 90)
    print("CHILmesh WNAT_Hagen Multi-Backend Benchmark")
    print("=" * 90)
    print()

    # Mesh info
    if "error" not in python_results:
        print(f"Mesh:    {python_results['n_verts']:,} vertices, {python_results['n_elems']:,} elements")
        print(f"Layers:  {python_results['n_layers']}")
    print(f"Python:  {platform.python_version()}")
    print(f"OS:      {platform.system()} {platform.machine()}")
    print()

    # Initialization
    print("### Initialization Performance (seconds)")
    print()
    print(f"{'Operation':<25} {'Python':>15} {'Rust':>15} {'C++':>15}")
    print("-" * 90)

    metrics = [
        ("Fast init (no layers)", "fast_init_s"),
        ("Skeletonization only", "skel_only_s"),
        ("Full init (with layers)", "full_init_s"),
        ("Quality analysis", "quality_s"),
    ]

    for label, key in metrics:
        py_val = python_results.get(key)
        rs_val = rust_results.get(key)
        cpp_val = cpp_results.get(key)

        py_str = _fmt(py_val) if py_val else "N/A"
        rs_str = _fmt(rs_val) if rs_val else "N/A"
        cpp_str = _fmt(cpp_val) if cpp_val else "N/A"

        print(f"{label:<25} {py_str:>15} {rs_str:>15} {cpp_str:>15}")

    print()
    print("### Query Performance (microseconds per call)")
    print()
    print(f"{'Operation':<25} {'Python':>15} {'Rust':>15} {'C++':>15}")
    print("-" * 90)

    py_ve = python_results.get("vert2edge_us")
    rs_ve = rust_results.get("vert2edge_us")
    cpp_ve = cpp_results.get("vert2edge_us")

    py_str = f"{py_ve:.2f}μs" if py_ve else "N/A"
    rs_str = f"{rs_ve:.2f}μs" if rs_ve else "N/A"
    cpp_str = f"{cpp_ve:.2f}μs" if cpp_ve else "N/A"

    print(f"{'Vert2Edge lookup':<25} {py_str:>15} {rs_str:>15} {cpp_str:>15}")

    print()
    print("### Speedup vs Python")
    print()
    print(f"{'Operation':<25} {'Rust':>15} {'C++':>15}")
    print("-" * 55)

    for label, key in metrics:
        py_val = python_results.get(key)
        rs_val = rust_results.get(key)
        cpp_val = cpp_results.get(key)

        if py_val and py_val > 0:
            rs_speedup = f"{py_val / rs_val:.1f}x" if rs_val and rs_val > 0 else "N/A"
            cpp_speedup = f"{py_val / cpp_val:.1f}x" if cpp_val and cpp_val > 0 else "N/A"
        else:
            rs_speedup = "N/A"
            cpp_speedup = "N/A"

        print(f"{label:<25} {rs_speedup:>15} {cpp_speedup:>15}")

    print()
    print("=" * 90)

def main() -> int:
    import argparse
    parser = argparse.ArgumentParser(description="Benchmark all CHILmesh backends.")
    parser.add_argument("mesh_path", nargs="?", help="Path to WNAT_Hagen.14")
    args = parser.parse_args()

    try:
        mesh_path = _locate_mesh(args.mesh_path)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Mesh: {mesh_path}")
    print()

    print("Benchmarking Python backend...", flush=True)
    py_results = bench_python(mesh_path)
    print("Benchmarking Rust backend...", flush=True)
    rs_results = bench_rust(mesh_path)
    print("Benchmarking C++ backend...", flush=True)
    cpp_results = bench_cpp(mesh_path)

    print_report(py_results, rs_results, cpp_results)

    # JSON output
    output = {
        "mesh": str(mesh_path),
        "python": py_results,
        "rust": rs_results,
        "cpp": cpp_results,
    }
    json_path = Path("/tmp/wnat_all_backends.json")
    json_path.write_text(json.dumps(output, indent=2))
    print(f"Results saved: {json_path}")

    return 0

if __name__ == "__main__":
    sys.exit(main())
