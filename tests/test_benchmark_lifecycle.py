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

    m = CHILmesh.read_from_fort14(DONUT)
    conn = np.asarray(m.connectivity_list).astype(np.int64)
    pts = np.asarray(m.points).astype(float)
    return conn, pts


def test_lifecycle_reports_all_stages(bench, donut_arrays):
    conn, pts = donut_arrays
    out = bench.bench_lifecycle(conn, pts, repeats=1, smooth_iters=3)
    for key in ("adj_s", "skel_s", "quality_s", "fem_smooth_s", "angle_smooth_s"):
        assert key in out, f"missing lifecycle stage: {key}"
        assert out[key] > 0.0, f"stage {key} reported non-positive time"
    assert out["smooth_iters"] == 3
    assert out["n_verts"] == pts.shape[0]


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
