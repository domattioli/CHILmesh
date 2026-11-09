"""W1: CHILmesh smoothers under Track 1 (frozen) and Track 2, external tools degrade to not-run."""
from __future__ import annotations

import numpy as np
import pytest

from benchmarks.competitors.w1 import tools, track1, track2

MESH = dict(name="synthetic", version="v0", sha256="0" * 64, provenance="synthetic", crs="cartesian",
            units="m", n_nodes=0, n_elems=0)


def _perturbed(pts, conn, seed=0):
    from benchmarks.competitors.metrics.quality import boundary_nodes

    rng = np.random.default_rng(seed)
    out = pts.copy()
    inner = np.setdiff1d(np.arange(len(pts)), boundary_nodes(conn))
    out[inner] += rng.uniform(-0.15, 0.15, size=(len(inner), 2))
    return out


def test_chilmesh_has_all_three_smoothers():
    assert set(tools.CHILMESH_SMOOTHERS) == {"smooth", "angle_based_smoother", "direct_smoother"}


@pytest.mark.parametrize("smoother", sorted(tools.CHILMESH_SMOOTHERS))
def test_track1_chilmesh_record_is_ok_with_zero_boundary_displacement(grid, smoother):
    kind, pts, conn = grid
    pts = _perturbed(pts, conn)
    depths = np.full(len(pts), 10.0)
    rec = track1.run_chilmesh(smoother, pts, conn, depths, mesh=MESH, crs="cartesian", dt=1.0)
    mv = rec["metric_values"]
    if smoother == "direct_smoother" and rec["status"] == "invalid":
        # direct_smoother pins the boundary with a kinf=1e12 penalty, so a mixed mesh can leave a
        # ~1e-11 residual. The harness check is exact on purpose and reports it instead of hiding it.
        assert "boundary" in rec["reason"]
        assert 0.0 < mv["boundary_displacement"]["value"] < 1e-8
    else:
        assert rec["status"] == "ok", rec["reason"]
        assert mv["boundary_displacement"]["value"] == 0.0
    assert rec["track"] == 1 and rec["workflow"] == "w1"
    for key in ("q_min_before", "q_min_after", "q_mean_minus_3sigma_after", "pct_angles_lt_10_after",
                "pct_irregular_valence_after", "max_courant_after", "dt_max_after",
                "boundary_displacement", "wall_s"):
        assert key in mv, key
    assert rec["tool"] == f"chilmesh:{smoother}"


def test_track1_boundary_violation_marks_invalid(monkeypatch):
    pts, conn = __import__("_grids").tri_grid()

    def bad(points, connectivity):
        out = points.copy()
        out[0] += 1e-9  # corner node is on the boundary
        return out

    monkeypatch.setitem(tools.CHILMESH_SMOOTHERS, "smooth", bad)
    rec = track1.run_chilmesh("smooth", pts, conn, np.full(len(pts), 10.0), mesh=MESH, crs="cartesian", dt=1.0)
    assert rec["status"] == "invalid"
    assert "boundary" in rec["reason"]


def test_external_tool_without_env_is_not_run(monkeypatch):
    monkeypatch.delenv("CHILMESH_BENCH_PY_OCEANMESH", raising=False)
    pts, conn = __import__("_grids").tri_grid()
    rec = track1.run_external("oceanmesh", pts, conn, np.full(len(pts), 10.0), mesh=MESH, crs="cartesian",
                              dt=1.0, workdir=None)
    assert rec["status"] == "not-run"
    assert "CHILMESH_BENCH_PY_OCEANMESH" in rec["reason"]


def test_octave_missing_is_not_run(monkeypatch):
    monkeypatch.setattr(tools.shutil, "which", lambda name: None)
    pts, conn = __import__("_grids").tri_grid()
    rec = track1.run_external("oceanmesh2d", pts, conn, np.full(len(pts), 10.0), mesh=MESH, crs="cartesian",
                              dt=1.0, workdir=None)
    assert rec["status"] == "not-run"
    assert "octave" in rec["reason"].lower()


def test_track2_notes_chilmesh_deletes_no_vertices(grid):
    kind, pts, conn = grid
    pts = _perturbed(pts, conn)
    rec = track2.run_chilmesh("smooth", pts, conn, np.full(len(pts), 10.0), mesh=MESH, crs="cartesian", dt=1.0)
    assert rec["track"] == 2
    assert rec["metric_values"]["vertices_deleted"]["value"] == 0
    assert "deletes no vertices" in rec["reason"] or "deletes none" in rec["reason"]
    assert rec["status"] == "ok"
