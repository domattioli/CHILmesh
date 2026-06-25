"""Tests for the post-generation mesh-lifecycle benchmark stages (#155).

The benchmark script lives under ``scripts/`` (not an installed package), so it
is loaded by file path. These tests guard that every lifecycle stage runs and
reports a positive wall-clock time, and that the truss stage resolves a built-in
SDF for the bundled domain fixtures.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
DONUT = REPO / "src" / "chilmesh" / "data" / "donut_domain.fort.14"


def _load_benchmark():
    path = REPO / "scripts" / "benchmark.py"
    spec = importlib.util.spec_from_file_location("chilmesh_benchmark", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def bench():
    return _load_benchmark()


@pytest.fixture(scope="module")
def donut_arrays():
    from chilmesh import CHILmesh

    m = CHILmesh.read_from_fort14(DONUT, compute_layers=False, compute_adjacencies=True)
    conn = np.asarray(m.connectivity_list).astype(np.int64)
    pts = np.asarray(m.points).astype(float)
    return conn, pts


@pytest.mark.slow
def test_lifecycle_reports_all_stages(bench, donut_arrays):
    conn, pts = donut_arrays
    out = bench.bench_lifecycle(conn, pts, repeats=1, smooth_iters=3)
    for key in ("adj_s", "skel_s", "quality_s", "fem_smooth_s", "angle_smooth_s"):
        assert key in out, f"missing lifecycle stage: {key}"
        assert out[key] > 0.0, f"stage {key} reported non-positive time"
    assert out["smooth_iters"] == 3
    assert out["n_verts"] == pts.shape[0]
    assert out["fem_skipped"] is False


def test_gate_heavy_decision(bench):
    # gate disabled (0) or mesh within budget → run everything
    assert bench._gate_heavy(1_000, 0) == (False, "")
    assert bench._gate_heavy(40, 50)[0] is False
    assert bench._gate_heavy(50, 50)[0] is False  # boundary: equal is allowed
    # mesh over the gate → skip, with both counts in the reason
    skip, reason = bench._gate_heavy(4_000_000, 500_000)
    assert skip is True
    assert "4,000,000" in reason and "500,000" in reason


def test_fem_gated_is_solver_aware(bench):
    # heavy mesh + direct solver → FEM gated out (OOM risk)
    assert bench._fem_gated(True, "direct") is True
    # heavy mesh + iterative MINRES (#168) → FEM NOT gated (memory-safe)
    assert bench._fem_gated(True, "iterative") is False
    # within budget (not heavy) → never gated regardless of solver
    assert bench._fem_gated(False, "direct") is False
    assert bench._fem_gated(False, "iterative") is False


@pytest.mark.slow
def test_max_elements_gate_skips_fem(bench, donut_arrays):
    conn, pts = donut_arrays
    out = bench.bench_lifecycle(conn, pts, repeats=1, smooth_iters=3, skip_fem=True)
    assert out["fem_skipped"] is True
    assert out["fem_smooth_s"] is None
    # scalable stages still report positive times
    for key in ("adj_s", "skel_s", "quality_s", "angle_smooth_s"):
        assert out[key] > 0.0, f"stage {key} reported non-positive time"


def test_sdf_registry_resolves_bundled_domains(bench):
    name, sdf = bench._resolve_sdf("src/chilmesh/data/donut_domain.fort.14", None)
    assert name == "donut" and callable(sdf)

    name, sdf = bench._resolve_sdf("whatever/Block_O.14", None)
    assert sdf is None and "no built-in SDF" in name

    name, sdf = bench._resolve_sdf("Block_O.14", "square")
    assert name == "square" and callable(sdf)

    name, sdf = bench._resolve_sdf("donut_domain.fort.14", "none")
    assert sdf is None


def test_truss_stage_times_triangle_mesh(bench, donut_arrays):
    conn, pts = donut_arrays
    _, sdf = bench._resolve_sdf("donut_domain.fort.14", None)
    out = bench.bench_truss(conn, pts, sdf, repeats=1, niter=10)
    assert out is not None, "truss stage skipped on a triangle-only donut mesh"
    assert out["truss_s"] > 0.0
    assert out["niter"] == 10
