"""W2 runner: parse, adjacency, point location and render, one fresh subprocess per measurement."""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Optional

import chilmesh

from .. import machine, noaa, records, valence
from ..measure import run_measured
from . import tools
from .scaling import fit_scaling

REPO_ROOT = Path(__file__).resolve().parents[3]


def _json_tail(stdout_tail: str) -> Optional[dict]:
    for line in reversed((stdout_tail or "").strip().splitlines()):
        line = line.strip()
        if line.startswith("{"):
            try:
                return json.loads(line)
            except ValueError:
                return None
    return None


def _worker(op, mesh_path, py, env_extra, timeout_s, extra_args=()):
    cmd = [py, "-m", "benchmarks.competitors.w2.measures", op, "--mesh", str(mesh_path), *extra_args]
    env = {**os.environ, **env_extra, "PYTHONPATH": str(REPO_ROOT) + os.pathsep + os.environ.get("PYTHONPATH", "")}
    res = run_measured(cmd, timeout_s=timeout_s, env=env, cwd=str(REPO_ROOT))
    res["json"] = _json_tail(res["stdout_tail"]) if res["status"] == "ok" else None
    if res["status"] == "ok" and res["json"] is None:
        res["status"], res["reason"] = "failed", "worker printed no JSON result"
    return res


def _mesh_meta(name, fetch, crs, units, n_nodes=0, n_elems=0, provenance="Valence registry"):
    return dict(name=name, version=fetch.get("sha256", "unknown")[:12], sha256=fetch.get("sha256", "unknown"),
                provenance=provenance, crs=crs, units=units, n_nodes=n_nodes, n_elems=n_elems)


def _record(tool, version, license_, mode, mesh, metric_values, status, reason, cache, command, res, idx):
    return records.make_record(
        tool=tool, tool_version=version, license=license_, run_mode=mode, mesh=mesh, workflow="w2",
        metric_values=metric_values, status=status, reason=reason, cache=cache, command=command,
        wall_s=res.get("wall_s") if res else None, peak_rss_bytes=res.get("peak_rss_bytes") if res else None,
        stdout_tail=(res or {}).get("stdout_tail"), stderr_tail=(res or {}).get("stderr_tail"),
        repeat_index=idx, machine=machine.describe_machine(), chilmesh_git_sha=machine.chilmesh_git_sha())


def _na_secondary(tool):
    why = tools.COMPETITORS[tool]["na"]
    return {m: records.metric(None, na=True, reason=why) for m in tools.SECONDARY}


def not_run_records(name, reason, cache, crs="unknown", units="unknown"):
    """One not-run record per tool, so every cell carries a visible reason (SC-001)."""
    mesh = _mesh_meta(name, {}, crs, units)
    out = [_record("chilmesh", chilmesh.__version__, "PolyForm-NC", "subprocess", mesh, {}, "not-run", reason,
                   cache, [], None, 0)]
    for t, cfg in tools.COMPETITORS.items():
        out.append(_record(t, "unknown", cfg["license"], "subprocess", mesh, {}, "not-run", reason, cache, [],
                           None, 0))
    return out


def run_chilmesh_once(mesh_path, mesh, cache, idx, timeout_s, aux_dir=None):
    mv, res_all, reason, status = {}, [], "", "ok"
    locate_args = ("--aux-dir", str(aux_dir)) if aux_dir else ()
    parse = _worker("parse", mesh_path, sys.executable, {}, timeout_s)
    res_all.append(parse)
    if parse["status"] != "ok":
        return _record("chilmesh", chilmesh.__version__, "PolyForm-NC", "subprocess", mesh, {}, "failed",
                       parse["reason"], cache, parse["command"], parse, idx)
    mv["parse_s"] = records.metric(parse["json"]["parse_s"])
    mesh = {**mesh, "n_nodes": parse["json"]["n_nodes"], "n_elems": parse["json"]["n_elems"]}
    plan = [("adjacency", {"CHILMESH_BACKEND": "cpp"}, "adj_cpp_s", "adjacency_s"),
            ("adjacency", {"CHILMESH_BACKEND": "python"}, "adj_python_s", "adjacency_s"),
            ("locate", {}, "tree_build_s", "tree_build_s"),
            ("locate", {}, "query_per_s", "query_per_s"),
            ("render", {}, "render_s", "render_s")]
    cache_locate = None
    for op, env_extra, metric_name, key in plan:
        if op == "locate" and cache_locate is not None:
            res = cache_locate
        else:
            res = _worker(op, mesh_path, sys.executable, env_extra, timeout_s,
                          locate_args if op in ("locate", "render") else ())
            if op == "locate":
                cache_locate = res
        res_all.append(res)
        js = res["json"] if res["status"] == "ok" else {}
        if res["status"] == "ok" and key in js and js[key] is not None:
            mv[metric_name] = records.metric(js[key])
        elif res["status"] == "ok" and js.get("status") == "not-run":
            mv[metric_name] = records.metric(None, na=True, reason=js.get("reason") or f"{op} not run")
        else:
            why = res["reason"] or js.get("reason") or f"{op} produced no {key}"
            if metric_name == "adj_cpp_s" and res["status"] == "ok":
                mv[metric_name] = records.metric(None, na=True, reason="C++ backend not available in this build")
            else:
                mv[metric_name] = records.metric(None, na=True, reason=f"{op} failed: {why}")
                status, reason = "failed", (reason + "; " if reason else "") + f"{op} failed: {why}"
    last = res_all[-1]
    return _record("chilmesh", chilmesh.__version__, "PolyForm-NC", "subprocess", mesh, mv, status, reason,
                   cache, last["command"], parse, idx)


def run_competitor_once(tool, mesh_path, mesh, cache, idx, timeout_s):
    cfg = tools.COMPETITORS[tool]
    if cfg["input"] != "fort.14":
        return _record(tool, "unknown", cfg["license"], "subprocess", mesh, {}, "not-run",
                       tools.XARRAY_FORT14_REASON, cache, [], None, idx)
    py = tools.competitor_python(tool)
    if not py:
        return _record(tool, "unknown", cfg["license"], "subprocess", mesh, {}, "not-run",
                       f"no competitor environment configured: set {cfg['env']} to its python", cache, [],
                       None, idx)
    cmd = [py, str(cfg["script"]), "--mesh", str(mesh_path)]
    res = run_measured(cmd, timeout_s=timeout_s, env=dict(os.environ))
    out = _json_tail(res["stdout_tail"]) if res["status"] == "ok" else None
    if res["status"] != "ok" or not out or "parse_s" not in out:
        return _record(tool, "unknown", cfg["license"], "subprocess", mesh, {}, "failed",
                       res["reason"] or "adapter printed no parse_s", cache, cmd, res, idx)
    mv = {"parse_s": records.metric(out["parse_s"]), **_na_secondary(tool)}
    return _record(tool, out.get("version", "unknown"), cfg["license"], "subprocess", mesh, mv, "ok", "", cache,
                   cmd, res, idx)


def station_record(mesh_path, mesh, cache, aux_dir, timeout_s, idx=0, worker=_worker):
    """CHILmesh station sampling versus ``points.cwl.nc`` (T023): one record, not-run without NOAA inputs."""
    def rec(mv, status, reason, res=None):
        return _record("chilmesh", chilmesh.__version__, "PolyForm-NC", "subprocess", mesh, mv, status, reason,
                       cache, (res or {}).get("command", []), res, idx)

    stations_path, points_path = noaa.aux_paths(aux_dir)
    if stations_path is None or points_path is None:
        return rec({}, "not-run", noaa.NOAA_INPUTS_REASON)
    if noaa.maxele_path(aux_dir) is None:
        return rec({}, "not-run", noaa.MAXELE_REASON)
    res = worker("stations", mesh_path, sys.executable, {}, timeout_s, ("--aux-dir", str(aux_dir)))
    out = res.get("json") or {}
    if res["status"] != "ok":
        return rec({}, "failed", res["reason"] or "stations worker failed", res)
    if out.get("status") == "failed":
        return rec({}, "failed", out.get("reason") or "station sampling failed", res)
    if out.get("status") != "ok":
        return rec({}, "not-run", out.get("reason") or noaa.NOAA_INPUTS_REASON, res)
    if out.get("max_abs_diff") is None:
        return rec({}, "failed", f"all {out.get('n_outside', 0)} stations fell outside the mesh", res)
    return rec({"station_max_abs_diff": records.metric(out["max_abs_diff"])}, "ok", "", res)


def run_mesh(name, out_dir, *, repeats=3, timeout_s=3600, cache="unknown", crs="unknown", units="unknown",
             aux_dir=None):
    """Fetch ``name`` from Valence and write one record per tool and repeat. Returns the records.

    ``aux_dir`` holds the fetched CO-OPS station list, ``points.cwl.nc`` and ``maxele.63.nc``. Without the
    first two, station sampling is not-run ("NOAA inputs not fetched") and point location uses centroids.
    Without ``maxele.63.nc``, station sampling is not-run ("NOAA inputs not fetched (maxele)") and the
    render metric is not applicable with that reason.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fetch = valence.fetch_mesh_file(name, out_dir / "meshes")
    if fetch["status"] != "ok":
        recs = not_run_records(name, fetch["reason"], cache, crs, units)
    else:
        mesh = _mesh_meta(name, fetch, crs, units)
        recs = []
        for idx in range(repeats):
            recs.append(run_chilmesh_once(fetch["path"], mesh, cache, idx, timeout_s, aux_dir))
            for t in sorted(tools.COMPETITORS):
                recs.append(run_competitor_once(t, fetch["path"], mesh, cache, idx, timeout_s))
        recs.append(station_record(fetch["path"], mesh, cache, aux_dir, timeout_s))
    for i, rec in enumerate(recs):
        records.write_record(rec, out_dir / f"w2_{name}_{rec['tool']}_{rec['repeat_index']}_{i}.json")
    return recs


def scaling_summary(sizes, times):
    """Fit only with at least 3 sizes; otherwise say so instead of fitting two points."""
    if len(sizes) < 3:
        return {"best": None, "note": f"only {len(sizes)} mesh size(s) served; need at least 3 (SC-004)"}
    return fit_scaling(sizes, times)


def _query_time(rec):
    cell = rec["metric_values"].get("query_per_s")
    return cell["value"] if rec["status"] == "ok" and cell and not cell["na"] else None


def collect_scaling(recs, out_dir):
    """Collect CHILmesh query time per mesh size from W2 records and fit it (T025, FR-013, SC-004).

    Size is the element count; time is the mean ``query_per_s`` over the ok repeats. A mesh with no
    usable query time (could not load, failed, not run) goes into ``not_run`` with its reason and,
    unless it already has a not-run record, gets one written to ``out_dir``. The fit goes to
    ``out_dir/scaling/w2_scaling.json``, outside the run-record directory listing.
    """
    out_dir = Path(out_dir)
    by_mesh = {}
    for r in recs:
        if r["tool"] == "chilmesh" and r["workflow"] == "w2":
            by_mesh.setdefault(r["mesh"]["name"], []).append(r)
    points, not_run = [], []
    for name, rs in by_mesh.items():
        usable = [(r["mesh"]["n_elems"], t) for r in rs if r["mesh"]["n_elems"] > 0
                  for t in [_query_time(r)] if t is not None]
        if usable:
            points.append((usable[0][0], sum(t for _, t in usable) / len(usable)))
            continue
        why = next((r["reason"] for r in rs if r["reason"]), "no query time collected")
        not_run.append({"mesh": name, "reason": f"scaling: {why}"})
        if not any(r["status"] == "not-run" for r in rs):
            rec = _record("chilmesh", chilmesh.__version__, "PolyForm-NC", "subprocess", rs[0]["mesh"], {},
                          "not-run", f"scaling: {why}", rs[0]["cache"], [], None, 0)
            records.write_record(rec, out_dir / f"w2_{name}_chilmesh_scaling_not-run.json")
    points.sort()
    sizes, times = [p[0] for p in points], [p[1] for p in points]
    result = {"sizes": sizes, "times": times, "summary": scaling_summary(sizes, times), "not_run": not_run}
    path = out_dir / "scaling" / "w2_scaling.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    return result
