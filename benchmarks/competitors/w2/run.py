"""W2 runner: parse, adjacency, point location and render, one fresh subprocess per measurement."""
from __future__ import annotations

import json
import hashlib
import importlib.metadata
import os
import sys
from pathlib import Path
from typing import Optional

import chilmesh

from .. import machine, noaa, records, valence
from ..measure import run_measured
from . import tools
from .scaling import fit_scaling
from .measures import install_mb

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


_consumer_worker = _worker


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
    steps = {"parse": parse}
    for op, names in (("write", ("write_fort14_s", "roundtrip_ok")),
                      ("full-init", ("full_init_s",))):
        result = _consumer_worker(op, mesh_path, sys.executable, {}, timeout_s)
        steps[op] = result
        out = result.get("json") or {}
        for name in names:
            if result["status"] == "ok" and name in out:
                mv[name] = records.metric(out[name], reason=out.get("reason", "") if name == "roundtrip_ok" else "")
                if name == "full_init_s":
                    mv[name]["backend_info"] = out.get("backend_info")
            else:
                mv[name] = records.metric(None, na=True, reason=result.get("reason") or f"{op} produced no {name}")
    imported = run_measured([sys.executable, "-c", "import chilmesh"], timeout_s=timeout_s,
                            env=dict(os.environ), cwd=str(REPO_ROOT))
    steps["import"] = imported
    mv["import_s"] = (records.metric(imported["wall_s"]) if imported["status"] == "ok" else
                      records.metric(None, na=True, reason=imported["reason"] or "import failed"))
    try:
        mv["install_mb"] = records.metric(install_mb())
    except (OSError, ValueError, importlib.metadata.PackageNotFoundError) as exc:
        mv["install_mb"] = records.metric(None, na=True, reason=str(exc))
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
        steps[op + ("-" + env_extra["CHILMESH_BACKEND"] if "CHILMESH_BACKEND" in env_extra else "")] = res
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
    rss = {name: result["peak_rss_bytes"] / (mesh["n_elems"] / 1e6)
           for name, result in steps.items() if result.get("peak_rss_bytes") is not None and mesh["n_elems"] > 0}
    mv["peak_rss_per_melem"] = (records.metric(rss) if rss else
                                records.metric(None, na=True, reason="peak RSS or element count unavailable"))
    last = res_all[-1]
    return _record("chilmesh", chilmesh.__version__, "PolyForm-NC", "subprocess", mesh, mv, status, reason,
                   cache, last["command"], parse, idx)


def run_competitor_once(tool, mesh_path, mesh, cache, idx, timeout_s):
    cfg = tools.COMPETITORS[tool]
    if cfg["input"] != "fort.14":
        return _record(tool, "unknown", cfg["license"], "subprocess", mesh, {}, "not-run",
                       tools.NETCDF_ONLY_REASON.format(tool=tool), cache, [], None, idx)
    py = tools.competitor_python(tool)
    if not py:
        return _record(tool, "unknown", cfg["license"], "subprocess", mesh, {}, "not-run",
                       f"no competitor environment configured: set {cfg['env']} to its python", cache, [],
                       None, idx)
    cmd = [py, str(cfg["script"]), "--mesh", str(mesh_path)]
    res = run_measured(cmd, timeout_s=timeout_s, env={**os.environ, "PYTHONPATH": str(REPO_ROOT)})
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
        if "interp_max_abs_diff" not in mv:
            mv["interp_max_abs_diff"] = (mv["station_max_abs_diff"].copy() if "station_max_abs_diff" in mv else
                                         records.metric(None, na=True, reason=reason or "station sampling unavailable"))
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


def fair_evidence(tool, operation, reference, result, *, options=None, conversion=None):
    """Classify one operation from injected full evidence, never stdout tails.

    US2 will supply measured sections. Missing controller evidence cannot become
    an eligible attempt. Native/conversion contracts have distinct identities.
    """
    from ..verification import verify_output, verify_written_output
    from ..campaign import validate_conversion
    contract = tools.operation_contract(tool, operation, options=options)
    summary = dict(contract=contract, equivalence=contract['equivalence'], outcome='invalid',
                   verification='incomplete', reason='', output_verification=None, result=result)
    if contract['equivalence'] == 'not-comparable':
        return dict(summary, outcome='not-comparable', verification='not-applicable', reason=contract['reason'])
    if reference is None or result is None:
        return dict(summary, outcome='not-run', reason='missing independent reference or configured supported tool/measurement result')
    if result.get('status') == 'not-run':
        return dict(summary, outcome='not-run', reason=result.get('reason') or 'operation prerequisite unavailable')
    if result.get('status') == 'failed' or result.get('returncode', 0) != 0:
        return dict(summary, outcome='failed', reason=result.get('reason') or 'operation child failed')
    reasons = []
    if not result.get('version') or result.get('version') == 'unknown': reasons.append('missing metadata distribution version')
    if result.get('input_format') is not None and result['input_format'] != contract['input_format']: reasons.append('operation input format mismatch')
    proof = verify_output(result.get('output') or {}, reference)
    summary.update(verification=proof['verification_state'], output_verification=proof)
    reasons.extend(proof['reasons'])
    if result.get('new_modules') != []: reasons.append('deferred imports or missing import evidence')
    if result.get('materialization') is not True: reasons.append('missing eager materialization proof')
    if operation == 'parse-netcdf':
        conversion = validate_conversion(conversion or {}, reference)
        summary['conversion'] = conversion
        reasons.extend(conversion['reasons'])
    if not operation.startswith('parse'):
        operation_proof = result.get('operation_output') or {}
        if operation == 'write':
            operation_proof = verify_written_output(operation_proof.get('path'), reference)
        if operation_proof.get('verification_state') != 'verified': reasons.append('operation-specific output/prerequisite proof incomplete')
        summary['operation_output'] = operation_proof
        actual = result.get('actual_backend')
        summary['actual_backend'] = actual
        if not actual: reasons.append('actual backend/work evidence missing')
    validity = result.get('measurement_validity') or {}
    if not all(validity.get(k, {}).get('valid') is True for k in ('time', 'memory', 'cache', 'load', 'runtime')):
        reasons.append('measurement evidence incomplete; US2 controller required')
    summary.update(outcome='invalid' if reasons else 'ok', reason='; '.join(reasons))
    return summary


def fair_record(identity, tool, operation, reference, result, *, options=None, conversion=None):
    """Produce a single validated v2 operation record from frozen identity/evidence.

    identity contains the existing records.py fields. Controller sections are
    injected without a bundled-operation or last-command substitution.
    """
    import copy
    summary = fair_evidence(tool, operation, reference, result, options=options, conversion=conversion)
    fields = copy.deepcopy(identity)
    reason = summary['reason'] or 'evidence unavailable'
    if reference is not None:
        fields['reference_id'] = reference['reference_id']
        fields['input_sha256'] = (conversion or {}).get('converted_sha256') if operation=='parse-netcdf' else reference.get('raw_file_sha256')
    fields['input_format'] = summary['contract']['input_format']
    fields.update(tool=tool, operation=operation, contract_id=summary['contract']['contract_id'],
                  backend=summary['contract']['backend_selected'], options=dict(options or {}),
                  outcome=summary['outcome'], equivalence=summary['equivalence'], verification=summary['verification'],
                  reason=summary['reason'], publication_eligible=False)
    for name, names in records.EVIDENCE_FIELDS.items():
        section = copy.deepcopy((result or {}).get(name) or {})
        if name == 'output_verification' and summary['output_verification'] is not None:
            section = copy.deepcopy(summary['output_verification'])
        section.setdefault('missing', {})
        for key in names.split():
            section.setdefault(key, None)
        for key, value in list(section.items()):
            if value is None: section['missing'][key] = reason
            elif key != 'missing': section['missing'].pop(key, None)
        fields[name] = section
    fields['measurement_validity'] = copy.deepcopy((result or {}).get('measurement_validity') or {
        key: {'valid': False, 'reason': reason} for key in ('time', 'memory', 'cache', 'load', 'runtime')})
    if summary.get('conversion') is not None: fields['conversion'] = copy.deepcopy(summary['conversion'])
    if result is not None:
        # A structured child result is retained whole, with a digest; bounded
        # stdout/stderr strings never serve as output transport.
        fields['adapter_result'] = copy.deepcopy(result)
        fields['adapter_result_sha256'] = hashlib.sha256(json.dumps(result, sort_keys=True, allow_nan=False).encode()).hexdigest()
        if result.get('argv'): fields['argv'] = list(result['argv'])
        if result.get('version'):
            fields['tool_version'] = result['version']
            if result['version'] != identity['tool_version'] and fields['outcome'] != 'failed':
                fields.update(outcome='invalid', reason=(fields['reason'] + '; version differs from frozen environment').lstrip('; '))
                fields['measurement_validity']['runtime'] = {'valid': False, 'reason': 'metadata version differs from frozen environment'}
    return records.make_fair_record(**fields)


def prepare_fair_mesh(name, version, asset, fetched, output_dir):
    """Generate independent reference evidence before any repeat preparation.

    Acquisition must be a Valence fetch result. This helper never searches for a
    local mesh or invokes another source when acquisition is unavailable.
    """
    return valence.reference_from_fetch(name, version, asset, fetched, output_dir)


def fair_operation_records(identities, reference, results, *, conversions=None):
    """US1 integration seam: one supplied result per predeclared operation.

    No child launches, cache operations or scheduling occur here. Missing keys
    remain explicit not-run records. US2 will provide the common controller.
    """
    conversions = conversions or {}
    return [fair_record(identity, identity['tool'], identity['operation'], reference,
                        results.get(identity['repeat_id']), options=identity.get('options'),
                        conversion=conversions.get(identity['tool'])) for identity in identities]


def measure_fair_attempt(identity, reference, command, *, mesh_class, expected_runtime,
                         artifact_dir, cache_reset=None, load_evidence=None, launcher=None,
                         observer=None, conversion=None):
    """Common fresh-child path for CHILmesh workers and competitor scripts.

    The caller prepares references/conversions/commands before reset. Nothing in
    this path reads the measured input between reset and child START. A callable
    load_evidence supplies synthetic evidence in offline integration tests; live
    observation covers baseline, imports, operation and verification.
    """
    import copy
    from ..cache import reset_cache
    from ..campaign import THREAD_ENV, LIMITS
    from ..machine import LoadObserver, evaluate_load
    limits=LIMITS[mesh_class]
    cache=(cache_reset or reset_cache)()
    if load_evidence is None:
        # Baseline after the reset: the blocking reset would otherwise leave a gap in the
        # 0.1 s observation cadence. The baseline reads /proc only, never the measured input.
        observer=observer or LoadObserver()
        observer.baseline()
    env={**os.environ,**THREAD_ENV,'PYTHONPATH':str(REPO_ROOT)}
    if identity['tool']=='chilmesh' and identity['backend'] in ('cpp','python'): env['CHILMESH_BACKEND']=identity['backend']
    try:
        measured=(launcher or run_measured)(command,timeout_s=limits['wall_s'],rss_limit_bytes=limits['rss_bytes'],
                                           env=env,cwd=str(REPO_ROOT),repeat_id=identity['repeat_id'],
                                           artifact_dir=artifact_dir,load_observer=observer)
    except KeyboardInterrupt:
        # The sampler kills/reaps in finally and saves partial controller evidence.
        # Consume it so the interrupted active attempt survives finalization.
        path=Path(artifact_dir)/'controller.json'
        if not path.exists(): raise
        measured=json.loads(path.read_text())
        measured.update(status='failed',reason='controller abort')
    child=copy.deepcopy(measured.get('adapter_result') or {})
    child.update({key:copy.deepcopy(measured.get(key) or {}) for key in ('measurement','protocol','termination')})
    child.update(status=measured['status'],returncode=measured.get('returncode'),reason=measured.get('reason',''),argv=list(command))
    launch=child['termination'].get('launch_s')
    cache=copy.deepcopy(cache);cache['launch_s']=launch
    cache.setdefault('missing',{})
    if launch is not None: cache['missing'].pop('launch_s',None)
    child['cache']=cache
    if load_evidence is not None: load=load_evidence()
    else:
        observer.poll()
        load=evaluate_load(observer.observations,launch_s=launch if launch is not None else observer.clock(),
                           host_id=identity['host_id'],allocation=identity['host_allocation'],
                           allowed_pids=(os.getpid(),child['protocol'].get('pid')))
    child['load']=copy.deepcopy(load)
    actual=child.get('runtime') or {}
    keys=('version','implementation','implementation_fields','platform','architecture','realpath','sha256')
    runtime_ok=(str(actual.get('version','')).startswith('3.11.') and all(actual.get(k)==expected_runtime.get(k) and actual.get(k) is not None for k in keys)
                and actual.get('thread_environment')==THREAD_ENV and expected_runtime.get('thread_environment')==THREAD_ENV)
    p=child['protocol'];m=child['measurement'];t=child['termination']
    time_ok=(all(type(p.get(k)) in (int,float) for k in ('start_s','end_s','ack_s'))
             and p['start_s']<=p['end_s']<=p['ack_s'] and m.get('operation_elapsed_s')==p['end_s']-p['start_s']
             and child.get('new_modules')==[] and child.get('materialization') is True)
    operations=cache.get('operations') or []
    cache_ok=(cache.get('success') is True and type(cache.get('cached_before_kib')) is int
              and type(cache.get('cached_after_kib')) is int and cache['cached_after_kib']<cache['cached_before_kib']
              and launch is not None and type(cache.get('reset_s')) in (int,float) and cache['reset_s']<launch
              and [o.get('operation') for o in operations]==['sync','drop_caches']
              and all(o.get('exit_code')==0 for o in operations) and operations[1].get('value')=='3')
    load_ok=(load.get('verdict')=='clean' and load.get('host_id')==identity['host_id'] and load.get('allocation')==identity['host_allocation']
             and load.get('pre_run_s',0)>=1 and load.get('interval_s')==.100
             and load.get('detector_policy')=={'cpu_percentage_points':5,'io_bytes_per_s':1048576})
    memory_ok=(m.get('coverage') is True and t.get('memory_events_before') is not None
               and t.get('memory_events_after') is not None and t.get('memory_events_delta') is not None)
    child['measurement_validity']={k:{'valid':bool(ok),'reason':'' if ok else why} for k,ok,why in (
        ('time',time_ok,'missing protocol/materialization or deferred imports'),('memory',memory_ok,'missing sampler/cgroup evidence'),
        ('cache',cache_ok,'reset evidence failed or wrong ordering'),('load',load_ok,'unknown/contaminated load or host/allocation drift'),
        ('runtime',runtime_ok,'child build/digest/thread policy differs from frozen Python 3.11'))}
    if identity['operation']=='parse-netcdf' and conversion is not None:
        conversion=copy.deepcopy(conversion)
        conversion['read_s']=child['measurement'].get('operation_elapsed_s')
    return fair_record(identity,identity['tool'],identity['operation'],reference,child,
                       options=identity.get('options'),conversion=conversion)


def run_fair_campaign(manifest, identities, references, commands, out_dir, *, conversions=None,
                      host_probe=None, runtime_probe=None, cache_reset=None, load_evidence=None,
                      launcher=None, lock_check=None, cgroup_check=None, lock_path=None):
    """Join a prepared frozen schedule to the shared measurement path.

    All maps are keyed by repeat ID except references (mesh alias) and conversions
    (repeat ID, with legacy tool fallback). Acquisition/conversion is deliberately completed before entering this
    controller. Missing supported commands retain reasoned not-run attempts.
    """
    from ..campaign import execute_schedule, runtime_identity, validate_lock_files, CampaignError
    conversions=conversions or {};out_dir=Path(out_dir)
    host=(host_probe or machine.host_identity)()
    runtime=(runtime_probe or runtime_identity)()
    environments={e['environment_id']:e for e in manifest['environment_manifests']}
    # Validate all caller-supplied identities and commands before the first reset.
    for entry in manifest['schedule']:
        ident=identities.get(entry['repeat_id'])
        if ident is None: raise CampaignError('missing prepared repeat identity')
        expected={k:manifest[k] for k in ('campaign_id','manifest_sha256')}
        expected.update({k:entry[k] for k in ('cell_id','series_id','repeat_id','repeat_index','execution_order','tool','operation','backend')})
        if any(ident.get(k)!=v for k,v in expected.items()): raise CampaignError('prepared identity differs from frozen schedule')
        contract=tools.operation_contract(entry['tool'],entry['operation'],options=ident.get('options'))
        if ident['contract_id']!=contract['contract_id']: raise CampaignError('operation contract differs from prepared identity')
        cmd=commands.get(entry['repeat_id'])
        if cmd and contract['equivalence']!='not-comparable':
            environment=environments.get(ident['environment_id'])
            if not environment or cmd[0]!=environment['executable'] or ident['tool_version']!=environment['version']:
                raise CampaignError('command/version differs from frozen environment')
            if ident['host_id']!=host['host_id'] or ident['host_allocation']!={k:host[k] for k in ('affinity','cpu_quota','cpuset','memory_limits')}:
                raise CampaignError('prepared host/allocation differs from frozen host')
    def conversion_for(entry):
        return conversions.get(entry['repeat_id'],conversions.get(entry['tool']))
    def execute(entry):
        ident=identities[entry['repeat_id']]
        reference=references.get(entry['mesh']);cmd=commands.get(entry['repeat_id'])
        contract=tools.operation_contract(entry['tool'],entry['operation'],options=ident.get('options'))
        if contract['equivalence']=='not-comparable' or not cmd or reference is None:
            # No child launches for this cell; the record still carries the (empty or prepared) argv.
            return fair_record(dict(ident,argv=list(cmd or [])),entry['tool'],entry['operation'],reference,None,options=ident.get('options'))
        return measure_fair_attempt(ident,reference,cmd,mesh_class='stofs' if 'stofs' in entry['mesh'].lower() else 'small',
            expected_runtime=runtime,artifact_dir=out_dir/'raw'/entry['repeat_id'],cache_reset=cache_reset,
            load_evidence=load_evidence,launcher=launcher,conversion=conversion_for(entry))
    def finalize(entry,reason,active=False):
        # Preserve the same frozen cell identity as any completed repeats. An
        # interrupted series must remain aggregatable as descriptive evidence.
        ident=dict(identities[entry['repeat_id']])
        if commands.get(entry['repeat_id']): ident['argv']=list(commands[entry['repeat_id']])
        row=fair_record(ident,entry['tool'],entry['operation'],references.get(entry['mesh']),None,
                        options=ident.get('options'),conversion=conversion_for(entry))
        if entry['operation']=='parse-netcdf' and conversion_for(entry) is not None:
            from ..campaign import validate_conversion
            lineage=dict(conversion_for(entry),read_s=None)
            row['conversion']=validate_conversion(lineage,references.get(entry['mesh']))
        if row['outcome']!='not-comparable':
            row['reason']=reason
            if active:
                row['outcome']='failed'
                row['termination']['trigger']='launch-error'
                row['termination']['missing'].pop('trigger',None)
        return records.validate_fair_record(row)
    rows=execute_schedule(manifest,out_dir,execute=execute,host_manifest=host,runtime_manifest=runtime,
        lock_check=lock_check or validate_lock_files,cgroup_check=cgroup_check,lock_path=lock_path,finalize=finalize)
    from ..records import _write_new
    refs={r['reference_id']:r for r in references.values() if r is not None}
    for entry,row in zip(manifest['schedule'],rows):
        refs.setdefault(row['reference_id'],dict(reference_id=row['reference_id'],name=entry['mesh'],
            exact_version='unavailable',raw_file_sha256=None))
    from ..campaign import write_references
    write_references(refs,out_dir)
    return rows
