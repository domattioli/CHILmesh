"""Regression tests for W1 timeout, cache labels, and exact boundary pinning."""
import time

import numpy as np
import pytest

from _grids import tri_grid
from benchmarks.competitors import cli
from benchmarks.competitors.w1 import run, tools, track1


MESH = dict(name="synthetic", version="v0", sha256="0" * 64, provenance="synthetic",
            crs="cartesian", units="m", n_nodes=0, n_elems=0)


def test_smoother_timeout(monkeypatch):
    def slow(points, connectivity):
        time.sleep(5)
        return points

    monkeypatch.setitem(tools.CHILMESH_SMOOTHERS, "slow", slow)
    pts, conn = tri_grid()
    rec = track1.run_chilmesh("slow", pts, conn, np.ones(len(pts)), mesh=MESH,
                              crs="cartesian", dt=1, timeout_s=0.1, cache="cold")
    assert rec["status"] == "failed"
    assert rec["reason"] == "timeout after 0.1 s"
    assert rec["cache"] == "cold"


def test_w1_cli_cache_reaches_records(monkeypatch, tmp_path):
    pts, conn = tri_grid()
    monkeypatch.setenv("CHILMESH_RUN_BENCH", "1")
    monkeypatch.setattr(run.valence, "fetch_mesh_file", lambda *args: {"status": "ok", "path": "unused"})
    monkeypatch.setattr(run, "_load", lambda path: (pts, conn, np.ones(len(pts))))
    monkeypatch.setattr(tools, "CHILMESH_SMOOTHERS", {"smooth": lambda p, c: p})
    monkeypatch.setattr(tools, "EXTERNAL", {})
    assert cli.main(["w1", "--mesh", "synthetic", "--out", str(tmp_path), "--cache", "warm"]) == 0
    recs = run.records.read_records(tmp_path)
    assert len(recs) == 2
    assert all(rec["cache"] == "warm" for rec in recs)


@pytest.mark.parametrize("solver", ["direct", "iterative"])
def test_direct_smoother_keeps_boundary_bit_identical(solver):
    from chilmesh import CHILmesh
    from benchmarks.competitors.metrics.quality import boundary_nodes
    pts, conn = tri_grid()
    mesh = CHILmesh(connectivity=conn, points=np.column_stack([pts[:, :2], np.zeros(len(pts))]))
    out = mesh.direct_smoother(solver=solver)
    b = boundary_nodes(conn)
    assert np.array_equal(out[b, :2], mesh.points[b, :2])
