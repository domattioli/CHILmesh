"""Track 1 (headline): smoothing with frozen topology and an exactly pinned boundary.

``build_record`` scores any tool's output with the shared metrics, so every tool is
judged by the same code. ``run_chilmesh`` and ``run_external`` produce one run record.
"""
from __future__ import annotations

import os
import shutil
import tempfile
import time
from pathlib import Path
from typing import Optional

import numpy as np

import chilmesh

from .. import machine, records
from ..measure import run_measured
from ..metrics import courant, displacement, quality
from . import tools

TRACK = 1
FROZEN = True


def _cell(value, na_reason: str):
    return records.metric(None, na=True, reason=na_reason) if value is None else records.metric(value)


def build_record(*, tool, tool_version, license_, run_mode, mesh, crs, dt, courant_max, points_before,
                 conn_before, depths, points_after, conn_after, wall_s, peak_rss_bytes, command,
                 track, frozen, cache="unknown", repeat_index=0, extra_reason="", extra_metrics=None):
    disp = displacement.check_boundary(points_before, conn_before, points_after, conn_after,
                                       frozen_topology=frozen)
    before = quality.quality_summary(points_before, conn_before)
    after = quality.quality_summary(points_after, conn_after)
    # Depths follow nodes; if a tool deleted vertices the harness cannot align them, so
    # Courant is computed only when the node count is unchanged.
    if len(points_after) == len(depths):
        c_max = courant.max_courant(points_after, conn_after, depths, dt, crs=crs)
        dt_max = courant.dt_max(points_after, conn_after, depths, courant_max=courant_max, crs=crs)
        c_reason = "every edge is dry (depth <= 0)"
    else:
        c_max = dt_max = None
        c_reason = "vertices were deleted, nodal depths cannot be aligned to the output mesh"
    mv = {
        "q_min_before": records.metric(before["q_min"]),
        "q_min_after": records.metric(after["q_min"]),
        "q_mean_minus_3sigma_after": records.metric(after["q_mean_minus_3sigma"]),
        "meets_0_75_bar_after": records.metric(after["meets_0_75_bar"]),
        "pct_angles_lt_10_after": records.metric(after["pct_angles_lt_10"]),
        "pct_irregular_valence_after": _cell(after["pct_irregular_valence"], "mesh has no interior nodes"),
        "max_courant_after": _cell(c_max, c_reason),
        "dt_max_after": _cell(dt_max, c_reason),
        "boundary_displacement": records.metric(disp["max_boundary_displacement"]),
        "wall_s": records.metric(wall_s),
    }
    mv.update(extra_metrics or {})
    reason = "; ".join(x for x in (disp["reason"], extra_reason) if x)
    mesh = {**mesh, "n_nodes": int(len(points_before)), "n_elems": int(len(np.asarray(conn_before)))}
    return records.make_record(
        tool=tool, tool_version=tool_version, license=license_, run_mode=run_mode, mesh=mesh,
        workflow="w1", track=track, metric_values=mv, status=disp["status"], reason=reason, cache=cache,
        command=command, wall_s=wall_s, peak_rss_bytes=peak_rss_bytes, repeat_index=repeat_index,
        machine=machine.describe_machine(), chilmesh_git_sha=machine.chilmesh_git_sha())


def _failed(tool, license_, run_mode, mesh, track, reason, command, *, status="failed", cache="unknown",
            repeat_index=0, tool_version="unknown", wall_s=None, peak_rss_bytes=None):
    return records.make_record(
        tool=tool, tool_version=tool_version, license=license_, run_mode=run_mode, mesh=mesh, workflow="w1",
        track=track, metric_values={}, status=status, reason=reason, cache=cache, command=command,
        wall_s=wall_s, peak_rss_bytes=peak_rss_bytes, repeat_index=repeat_index,
        machine=machine.describe_machine(), chilmesh_git_sha=machine.chilmesh_git_sha())


def run_chilmesh(smoother, points, conn, depths, *, mesh, crs, dt, courant_max=1.0, track=TRACK, frozen=FROZEN,
                 cache="unknown", repeat_index=0, extra_metrics_fn=None):
    tool = f"chilmesh:{smoother}"
    mesh = {**mesh, "crs": crs}
    t0 = time.perf_counter()
    try:
        new_points = tools.CHILMESH_SMOOTHERS[smoother](points, conn)
    except Exception as exc:  # recorded, never silently dropped
        return _failed(tool, "PolyForm-NC", "in-process", mesh, track, f"{type(exc).__name__}: {exc}",
                       [smoother], cache=cache, repeat_index=repeat_index, tool_version=chilmesh.__version__)
    wall = time.perf_counter() - t0
    extra = extra_metrics_fn(points, new_points) if extra_metrics_fn else None
    return build_record(
        tool=tool, tool_version=chilmesh.__version__, license_="PolyForm-NC", run_mode="in-process", mesh=mesh,
        crs=crs, dt=dt, courant_max=courant_max, points_before=points, conn_before=conn, depths=depths,
        points_after=new_points, conn_after=conn, wall_s=wall, peak_rss_bytes=None, command=[smoother],
        track=track, frozen=frozen, cache=cache, repeat_index=repeat_index, extra_metrics=extra)


def run_external(tool, points, conn, depths, *, mesh, crs, dt, workdir, courant_max=1.0, track=TRACK,
                 frozen=FROZEN, fort14_path: Optional[str] = None, timeout_s=3600, cache="unknown",
                 repeat_index=0):
    """Run oceanmesh or OceanMesh2D as a subprocess; missing environment gives ``not-run``."""
    cfg = tools.EXTERNAL[tool]
    mesh = {**mesh, "crs": crs}
    lic, mode = cfg["license"], cfg["run_mode"]
    if tool == "oceanmesh":
        py = tools.external_python(tool)
        if not py:
            return _failed(tool, lic, mode, mesh, track,
                           "no competitor environment configured: set CHILMESH_BENCH_PY_OCEANMESH to the "
                           "oceanmesh environment's python", [], status="not-run", cache=cache,
                           repeat_index=repeat_index)
        work = Path(workdir or tempfile.mkdtemp(prefix="bench_oceanmesh_"))
        work.mkdir(parents=True, exist_ok=True)
        tri = np.asarray(conn)[:, :3]
        np.savez(work / "in.npz", points=np.asarray(points)[:, :2], cells=tri)
        out_path = work / "out.npz"
        cmd = [py, str(cfg["script"]), "--in", str(work / "in.npz"), "--out", str(out_path),
               "--mode", "frozen" if frozen else "full"]
        res = run_measured(cmd, timeout_s=timeout_s, env=dict(os.environ))
        if res["status"] != "ok":
            return _failed(tool, lic, mode, mesh, track, res["reason"] + " | " + res["stderr_tail"][-300:], cmd,
                           cache=cache, repeat_index=repeat_index, wall_s=res["wall_s"],
                           peak_rss_bytes=res["peak_rss_bytes"])
        out = np.load(out_path)
        pts_out, cells_out = out["points"], out["cells"]
    else:  # oceanmesh2d via GNU Octave
        octave = shutil.which("octave")
        if not octave:
            return _failed(tool, lic, mode, mesh, track,
                           "GNU Octave not found on PATH; OceanMesh2D not run", [], status="not-run",
                           cache=cache, repeat_index=repeat_index)
        if not fort14_path:
            return _failed(tool, lic, mode, mesh, track, "no fort.14 path given for OceanMesh2D", [],
                           status="not-run", cache=cache, repeat_index=repeat_index)
        work = Path(workdir or tempfile.mkdtemp(prefix="bench_om2d_"))
        work.mkdir(parents=True, exist_ok=True)
        out_base = work / "om2d_out"
        expr = (f"bench_mesh='{fort14_path}'; bench_out='{out_base}'; run('{cfg['script']}')")
        cmd = [octave, "--no-gui", "--eval", expr]
        res = run_measured(cmd, timeout_s=timeout_s, env=dict(os.environ))
        if res["status"] != "ok":
            return _failed(tool, lic, mode, mesh, track, "Octave run failed: " + res["stderr_tail"][-300:], cmd,
                           status="not-run", cache=cache, repeat_index=repeat_index, wall_s=res["wall_s"],
                           peak_rss_bytes=res["peak_rss_bytes"])
        pts_out = np.loadtxt(str(out_base) + ".points", delimiter=",")
        cells_out = np.loadtxt(str(out_base) + ".cells", delimiter=",", dtype=int) - 1  # Octave is 1-based
    extra = None
    if not frozen:
        extra = {"vertices_deleted": records.metric(int(len(points) - len(pts_out)))}
    return build_record(
        tool=tool, tool_version="see adapter output", license_=lic, run_mode=mode, mesh=mesh, crs=crs, dt=dt,
        courant_max=courant_max, points_before=points, conn_before=conn, depths=depths, points_after=pts_out,
        conn_after=cells_out, wall_s=res["wall_s"], peak_rss_bytes=res["peak_rss_bytes"], command=cmd,
        track=track, frozen=frozen, cache=cache, repeat_index=repeat_index, extra_metrics=extra)
