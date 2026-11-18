"""Run-record schema for the competitor benchmark harness.

One JSON record per tool x mesh x workflow x track x repeat. Machine statuses
are ``ok | failed | not-run | invalid``. A metric that does not apply carries
``value: null`` and ``na: true`` with a justification; the report labels it
"not applicable". The status ``not-run`` is displayed as "not run".
"""
from __future__ import annotations

import json
import hashlib
import os
import tempfile
import math
from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional

SCHEMA_VERSION = 1
STATUSES = ("ok", "failed", "not-run", "invalid")
CACHE_STATES = ("cold", "warm", "unknown")
MESH_FIELDS = ("name", "version", "sha256", "provenance", "crs", "units", "n_nodes", "n_elems")
REQUIRED_FIELDS = (
    "tool", "tool_version", "license", "run_mode", "mesh", "workflow", "track", "metric_values",
    "status", "reason", "cache", "command", "wall_s", "peak_rss_bytes", "repeat_index", "machine",
    "chilmesh_git_sha",
)


class RecordError(ValueError):
    """A run record or metric cell breaks the schema."""


class MissingCellError(RuntimeError):
    """A tool x metric cell has no value, failure reason or justified na."""


def metric(value, *, na: bool = False, reason: str = "") -> dict:
    """Build one metric cell. ``na=True`` needs ``value=None`` and a reason."""
    if na:
        if value is not None:
            raise RecordError("a not-applicable metric must have value None")
        if not str(reason).strip():
            raise RecordError("a not-applicable metric needs a justification")
        return {"value": None, "na": True, "reason": str(reason)}
    return {"value": value, "na": False, "reason": str(reason)}


def _check_cell(name: str, cell) -> None:
    if not isinstance(cell, dict) or "value" not in cell or "na" not in cell:
        raise RecordError(f"metric {name!r} is not a metric() cell")
    if cell["na"]:
        metric(cell["value"], na=True, reason=cell.get("reason", ""))
    elif cell["value"] is None:
        raise RecordError(f"metric {name!r} is null but not marked na")


def validate_record(rec: dict) -> dict:
    if rec.get("schema_version") == 2:
        return validate_fair_record(rec)
    if rec.get("schema_version") != 1:
        raise RecordError(f"unsupported schema_version {rec.get('schema_version')!r}")
    missing = [f for f in REQUIRED_FIELDS if f not in rec]
    if missing:
        raise RecordError(f"record is missing fields: {missing}")
    if rec["schema_version"] != SCHEMA_VERSION:
        raise RecordError(f"unsupported schema_version {rec['schema_version']!r}")
    if rec["status"] not in STATUSES:
        raise RecordError(f"status {rec['status']!r} not in {STATUSES}")
    if rec["cache"] not in CACHE_STATES:
        raise RecordError(f"cache {rec['cache']!r} not in {CACHE_STATES}")
    mesh = rec["mesh"]
    absent = [f for f in MESH_FIELDS if f not in mesh]
    if absent:
        raise RecordError(f"mesh is missing fields: {absent}")
    if rec["status"] in ("failed", "not-run", "invalid") and not str(rec["reason"]).strip():
        raise RecordError(f"status {rec['status']!r} needs a reason")
    for name, cell in rec["metric_values"].items():
        _check_cell(name, cell)
    return rec


def make_record(**fields) -> dict:
    """Create and validate a record; optional fields default to None or empty."""
    rec = {"schema_version": SCHEMA_VERSION}
    rec.update(fields)
    for opt in ("track", "wall_s", "peak_rss_bytes", "stdout_tail", "stderr_tail"):
        rec.setdefault(opt, None)
    rec.setdefault("reason", "")
    rec.setdefault("metric_values", {})
    return validate_record(rec)


def write_record(rec: dict, path) -> Path:
    if rec.get("schema_version") == 2:
        write_fair_record(rec, path)
        return Path(path)
    validate_record(rec)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")
    return path


def read_record(path, *, expected_sha256=None) -> dict:
    raw = Path(path).read_bytes()
    if expected_sha256 is not None and hashlib.sha256(raw).hexdigest() != expected_sha256:
        raise RecordError("record hash mismatch")
    try:
        return validate_record(json.loads(raw, parse_constant=lambda x: (_ for _ in ()).throw(RecordError(f"nonfinite {x}"))))
    except (ValueError, UnicodeDecodeError) as exc:
        raise RecordError(f"invalid record JSON: {exc}") from exc


def read_records(directory) -> list[dict]:
    return [read_record(p) for p in sorted(Path(directory).glob("*.json"))]


def cell_for(records: Iterable[dict], tool: str, metric_name: str) -> Optional[dict]:
    """Describe the cell for one tool x metric, or None when nothing covers it.

    Returns ``{"kind": "value"|"na"|"failed"|"not-run"|"invalid", ...}``.
    """
    candidates = [r for r in records if r["tool"] == tool]
    for rec in candidates:
        cell = rec["metric_values"].get(metric_name)
        if cell is not None and rec["status"] in ("ok", "invalid"):
            if cell["na"]:
                return {"kind": "na", "reason": cell["reason"], "record": rec}
            kind = "value" if rec["status"] == "ok" else "invalid"
            return {"kind": kind, "value": cell["value"], "reason": rec["reason"], "record": rec}
    for rec in candidates:
        if rec["status"] in ("failed", "not-run", "invalid"):
            return {"kind": rec["status"], "reason": rec["reason"], "record": rec}
    return None


def check_coverage(records: list[dict], tools: Iterable[str], metrics: Iterable[str]) -> None:
    """Fail loudly when any tool x metric cell has no value, reason or justified na (SC-001)."""
    metrics = list(metrics)
    missing = [(t, m) for t in tools for m in metrics if cell_for(records, t, m) is None]
    if missing:
        listing = ", ".join(f"{t} x {m}" for t, m in missing)
        raise MissingCellError(f"no value, failure reason or justified na for: {listing}")


FAIR_REQUIRED = ("campaign_id", "manifest_sha256", "cell_id", "series_id", "repeat_id", "repeat_index", "execution_order", "reference_id", "contract_id", "tool", "environment_id", "tool_version", "chilmesh_git_sha", "argv", "thread_environment", "host_id", "host_allocation", "backend", "outcome", "reason", "equivalence", "verification", "measurement_validity", "publication_eligible", "protocol", "measurement", "cache", "load", "output_verification", "termination")


# All evidence slots exist even when no child launches. Null slots require reasons.
FAIR_REQUIRED += ("operation", "options", "created_at")
EVIDENCE_FIELDS = {
    "protocol": "pid process_group_id ready_s start_s end_s ack_s clock receipt_timestamps new_modules materialization",
    "measurement": "operation_elapsed_s process_wall_s rss_baseline_bytes rss_operation_peak_bytes rss_whole_process_peak_bytes rss_growth_bytes last_rss_bytes wait4_lifetime_maxrss_bytes policy_version interval_s scope coverage gaps sample_count samples samples_path samples_sha256",
    "cache": "operations cached_before_kib cached_after_kib cached_before_bytes cached_after_bytes reset_utc reset_s launch_s success",
    "load": "pre_run_s interval_s observations allowlist detector_policy verdict host_id allocation",
    "output_verification": "n_nodes n_elems node_sha256 element_sha256 expected mapping_provenance normalization_revision reasons",
    "termination": "launch_utc exit_utc launch_s exit_s returncode signal limits trigger cgroup memory_events_before memory_events_after memory_events_delta memory_peak_bytes memory_peak_scope kill_diagnostics reap_diagnostics stdout_path stdout_sha256 stderr_path stderr_sha256 stdout_tail stderr_tail",
}
MEASUREMENT_VALUES = ("operation_elapsed_s", "process_wall_s", "rss_baseline_bytes", "rss_operation_peak_bytes",
                      "rss_whole_process_peak_bytes", "rss_growth_bytes", "last_rss_bytes", "wait4_lifetime_maxrss_bytes")
THREAD_ENV = {k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")}


def _require(condition, message):
    if not condition:
        raise RecordError(message)


def _hash(value):
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def _number(value):
    return type(value) in (int, float) and math.isfinite(value) and value >= 0


def _utc(value):
    try:
        return isinstance(value, str) and datetime.fromisoformat(value.replace("Z", "+00:00")).utcoffset().total_seconds() == 0
    except (ValueError, AttributeError):
        return False


def _section(rec, name, fields):
    section = rec[name]
    _require(isinstance(section, dict), f"{name} must be an evidence object")
    _require(all(k in section for k in fields.split()), f"{name} missing evidence fields")
    missing = section.get("missing", {})
    _require(isinstance(missing, dict), f"{name} missing reasons must be an object")
    for key, value in section.items():
        if value is None:
            _require(isinstance(missing.get(key), str) and bool(missing[key].strip()), f"{name}.{key} needs missing reason")
        elif key in missing:
            raise RecordError(f"{name}.{key} has both value and missing flag")
        if value is not None and key.endswith("sha256"):
            _require(_hash(value), f"invalid {name}.{key} digest")
    return section


def _validate_fair_record(rec):
    _require(isinstance(rec, dict) and rec.get("schema_version") == 2, "expected schema version 2")
    _require(all(k in rec for k in FAIR_REQUIRED), "fair record missing required fields")
    try:
        json.dumps(rec, allow_nan=False)
    except (ValueError, TypeError) as exc:
        raise RecordError(f"invalid JSON value: {exc}") from exc
    for key in ("campaign_id", "cell_id", "series_id", "repeat_id", "reference_id", "contract_id", "operation",
                "tool", "environment_id", "tool_version", "chilmesh_git_sha", "host_id", "backend"):
        _require(isinstance(rec[key], str) and bool(rec[key].strip()), f"missing identity: {key}")
    _require(_hash(rec["manifest_sha256"]), "invalid manifest hash")
    _require(_utc(rec["created_at"]), "created_at must be UTC")
    _require(isinstance(rec["options"], dict) and isinstance(rec["host_allocation"], dict), "invalid options/allocation")
    # A cell that launched no child (not-run / not-comparable) has no command; every other record needs one.
    _require(isinstance(rec["argv"], list) and all(isinstance(x, str) and x for x in rec["argv"])
             and (bool(rec["argv"]) or rec.get("outcome") in ("not-run", "not-comparable")), "invalid argv")
    _require(rec["thread_environment"] == THREAD_ENV, "invalid thread environment")
    _require(rec["outcome"] in ("ok", "invalid", "failed", "not-run", "not-comparable"), "invalid outcome")
    _require(rec["equivalence"] in ("like-for-like", "CHILmesh-only", "not-comparable"), "invalid equivalence")
    _require(rec["verification"] in ("verified", "failed", "incomplete", "not-applicable"), "invalid verification")
    _require(isinstance(rec["reason"], str) and (rec["outcome"] == "ok" or bool(rec["reason"].strip())), "non-ok outcome needs a reason")
    _require(rec["outcome"] != "not-comparable" or rec["equivalence"] == "not-comparable", "inconsistent not-comparable state")
    _require(rec["publication_eligible"] is False, "repeat cannot establish series publication eligibility")
    _require("stable" not in rec and "stability" not in rec, "stability belongs to series")
    for key in ("repeat_index", "execution_order"):
        _require(type(rec[key]) is int and rec[key] >= 0, f"invalid {key}")
    validity = rec["measurement_validity"]
    _require(isinstance(validity, dict), "invalid validity decisions")
    for key in ("time", "memory", "cache", "load", "runtime"):
        decision = validity.get(key)
        _require(isinstance(decision, dict) and type(decision.get("valid")) is bool and isinstance(decision.get("reason"), str), f"missing {key} validity")
        _require(decision["valid"] or bool(decision["reason"].strip()), f"invalid {key} needs a reason")
    sections = {key: _section(rec, key, names) for key, names in EVIDENCE_FIELDS.items()}
    m, p, c, load, output, t = (sections[k] for k in ("measurement", "protocol", "cache", "load", "output_verification", "termination"))
    for key in MEASUREMENT_VALUES:
        if m[key] is not None:
            _require(_number(m[key]) and (not key.endswith("bytes") or type(m[key]) is int), f"invalid measurement {key}")
    for key in ("n_nodes", "n_elems"):
        if output[key] is not None:
            _require(type(output[key]) is int and output[key] >= 0, f"invalid {key}")
    for key in ("pid", "process_group_id"):
        if p[key] is not None:
            _require(type(p[key]) is int and p[key] > 0, f"invalid protocol {key}")
    for key in ("ready_s", "start_s", "end_s", "ack_s"):
        if p[key] is not None:
            _require(_number(p[key]), f"invalid protocol {key}")
    for key in ("returncode", "signal"):
        if t[key] is not None:
            _require(type(t[key]) is int, f"invalid termination {key}")
    for key in ("cached_before_kib", "cached_after_kib", "cached_before_bytes", "cached_after_bytes"):
        if c[key] is not None:
            _require(type(c[key]) is int and c[key] >= 0, f"invalid cache {key}")
    if output["mapping_provenance"] is not None:
        _require(isinstance(output["mapping_provenance"], dict), "invalid output mapping provenance")
    if load["verdict"] is not None:
        _require(load["verdict"] in ("clean", "contaminated", "unknown"), "invalid load verdict")
    if m["samples"] is not None:
        _require(isinstance(m["samples"], list), "samples must be a list")
        for sample in m["samples"]:
            _require(isinstance(sample, dict) and _number(sample.get("monotonic_s")) and type(sample.get("pid")) is int and sample["pid"] > 0 and isinstance(sample.get("descendants"), list), "malformed RSS sample")
            if sample.get("rss_bytes") is None:
                _require(isinstance(sample.get("missing_reason"), str) and bool(sample["missing_reason"].strip()), "missing sample RSS needs reason")
            else:
                _require(type(sample["rss_bytes"]) is int and sample["rss_bytes"] >= 0, "invalid sample RSS")
    if t["memory_events_delta"] is not None:
        _require(all(isinstance(t[k], dict) for k in ("memory_events_before", "memory_events_after", "memory_events_delta")), "invalid cgroup counters")
        for key, value in t["memory_events_delta"].items():
            _require(type(value) is int and value >= 0 and value == t["memory_events_after"].get(key, 0) - t["memory_events_before"].get(key, 0), "inconsistent cgroup counter delta")
    if t["trigger"] is not None:
        _require(t["trigger"] in ("timeout", "rss-limit", "kernel-oom", "signal-unknown", "nonzero-exit", "launch-error", "none"), "invalid termination trigger")
    if rec["outcome"] in ("not-run", "not-comparable"):
        _require(p["pid"] is None and t["launch_s"] is None and all(m[k] is None for k in MEASUREMENT_VALUES), "unstarted record has execution evidence")
        _require(not any(v["valid"] for v in validity.values()), "unstarted record has valid measurements")
    else:
        _require(t["launch_s"] is not None or t["trigger"] == "launch-error", "executed outcome needs launch evidence")
        _require(t["trigger"] is not None, "executed outcome needs termination trigger")
    if rec["outcome"] == "ok":
        _require(rec["verification"] == "verified" and all(v["valid"] for v in validity.values()), "ok requires verified valid evidence")
        _require(rec["equivalence"] != "not-comparable", "ok needs a measured operation contract")
        _require(bool(rec["host_allocation"]), "missing allocation")
        _require(all(m[k] is not None for k in MEASUREMENT_VALUES if k != "wait4_lifetime_maxrss_bytes"), "ok missing measurements")
        _require(t["returncode"] == 0 and t["trigger"] == "none", "ok process did not succeed")
        _require(type(p["pid"]) is int and p["pid"] > 0 and type(p["process_group_id"]) is int and p["process_group_id"] > 0, "invalid process identity")
        stamps = [p[k] for k in ("ready_s", "start_s", "end_s", "ack_s")]
        _require(all(_number(v) for v in stamps) and stamps == sorted(stamps), "invalid protocol timestamps")
        _require(isinstance(p["clock"], dict) and p["clock"].get("monotonic") is True and p["clock"].get("implementation") and _number(p["clock"].get("resolution_s")), "missing clock metadata")
        _require(isinstance(p["receipt_timestamps"], dict) and all(_number(p["receipt_timestamps"].get(k)) for k in ("READY", "START", "END", "ACK")), "missing protocol receipts")
        _require(p["new_modules"] == [] and p["materialization"] is True, "deferred imports or lazy output")
        _require(m["operation_elapsed_s"] == p["end_s"] - p["start_s"], "elapsed does not bracket region")
        _require(m["rss_growth_bytes"] == max(0, m["rss_operation_peak_bytes"] - m["rss_baseline_bytes"]), "invalid RSS growth")
        _require(m["rss_whole_process_peak_bytes"] >= max(m["rss_operation_peak_bytes"], m["rss_baseline_bytes"]), "invalid whole RSS peak")
        _require(m["interval_s"] == .010 and m["scope"] == "child-only" and m["coverage"] is True and bool(m["policy_version"]), "invalid sampler policy")
        _require(isinstance(m["gaps"], list) and all(_number(gap) and gap <= .020 for gap in m["gaps"]), "incomplete sample gaps")
        _require(isinstance(m["samples"], list) and m["sample_count"] == len(m["samples"]) and bool(m["samples"]), "missing samples")
        for sample in m["samples"]:
            _require(isinstance(sample, dict) and _number(sample.get("monotonic_s")) and sample.get("pid") == p["pid"] and type(sample.get("rss_bytes")) is int and sample["rss_bytes"] >= 0 and sample.get("descendants") == [], "invalid/descendant sample")
        before = [s for s in m["samples"] if s["monotonic_s"] < p["start_s"]]
        region = [s for s in m["samples"] if p["start_s"] <= s["monotonic_s"] <= p["end_s"]]
        whole = [s for s in m["samples"] if t["launch_s"] <= s["monotonic_s"] <= p["end_s"]]
        _require(bool(before) and bool(region) and bool(whole), "missing baseline/region samples")
        _require(m["rss_baseline_bytes"] == before[-1]["rss_bytes"] and m["rss_operation_peak_bytes"] == max(s["rss_bytes"] for s in region) and m["rss_whole_process_peak_bytes"] == max(s["rss_bytes"] for s in whole), "samples do not support peaks")
        timestamps = [sample["monotonic_s"] for sample in m["samples"]]
        _require(timestamps == sorted(timestamps), "samples out of order")
        gaps = [b - a for a, b in zip(timestamps, timestamps[1:])]
        _require(len(gaps) == len(m["gaps"]) and all(math.isclose(a, b, abs_tol=1e-12) for a, b in zip(gaps, m["gaps"])), "reported gaps do not match samples")
        _require(all(p["receipt_timestamps"][event] >= p[key] for event, key in (("READY", "ready_s"), ("START", "start_s"), ("END", "end_s"), ("ACK", "ack_s"))), "receipt precedes child marker")
        _require(math.isclose(m["process_wall_s"], t["exit_s"] - t["launch_s"], abs_tol=1e-12), "process wall must span launch to exit")
        _require(c["success"] is True and c["cached_after_kib"] < c["cached_before_kib"] and c["cached_before_bytes"] == c["cached_before_kib"] * 1024 and c["cached_after_bytes"] == c["cached_after_kib"] * 1024, "invalid cache evidence")
        _require(c["reset_s"] < c["launch_s"] == t["launch_s"] and _utc(c["reset_utc"]), "cache reset must precede launch")
        _require(isinstance(c["operations"], list) and [x.get("operation") for x in c["operations"]] == ["sync", "drop_caches"] and all(x.get("exit_code") == 0 and _number(x.get("monotonic_s")) and _utc(x.get("utc")) for x in c["operations"]) and c["operations"][1].get("value") == "3", "invalid reset operations")
        _require(c["operations"][0]["monotonic_s"] <= c["operations"][1]["monotonic_s"] <= c["reset_s"], "reset operation timestamps out of order")
        _require(load["detector_policy"] == {"cpu_percentage_points": 5, "io_bytes_per_s": 1048576}, "load detector policy changed")
        _require(t["limits"] in ({"wall_s": 600, "rss_bytes": 8589934592}, {"wall_s": 1800, "rss_bytes": 51539607552}), "invalid watchdog limits")
        _require(load["pre_run_s"] >= 1 and load["interval_s"] == .100 and load["verdict"] == "clean" and bool(load["observations"]) and isinstance(load["allowlist"], list) and bool(load["detector_policy"]), "invalid load evidence")
        _require(load["host_id"] == rec["host_id"] and load["allocation"] == rec["host_allocation"], "host/allocation drift")
        _require(all(output[k] is not None for k in ("n_nodes", "n_elems", "node_sha256", "element_sha256", "mapping_provenance", "normalization_revision")) and output["reasons"] == [], "incomplete output verification")
        _require(isinstance(output["expected"], dict) and all(output["expected"].get(k) == output[k] for k in ("n_nodes", "n_elems", "node_sha256", "element_sha256")), "output mismatch")
        _require(_utc(t["launch_utc"]) and _utc(t["exit_utc"]) and _number(t["exit_s"]) and t["exit_s"] >= p["ack_s"], "invalid termination timing")
        _require(t["memory_events_before"] is not None and t["memory_events_after"] is not None and t["memory_events_delta"] is not None and bool(t["cgroup"]), "missing cgroup evidence")
        _require(all(t[k] is not None for k in ("limits", "kill_diagnostics", "reap_diagnostics", "stdout_path", "stdout_sha256", "stderr_path", "stderr_sha256", "stdout_tail", "stderr_tail")), "missing termination diagnostics")
    if rec["verification"] == "verified":
        _require(output["expected"] is not None and bool(output["mapping_provenance"]), "verified needs actual output evidence")
    if t["trigger"] == "kernel-oom":
        _require(isinstance(t["memory_events_delta"], dict) and t["memory_events_delta"].get("oom_kill", 0) > 0 and t["signal"] == 9, "kernel OOM needs attributable counter evidence")
    return rec


def validate_fair_record(rec):
    """Validate v2 without filling absent historical evidence or changing state."""
    try:
        return _validate_fair_record(rec)
    except (KeyError, TypeError, AttributeError, OverflowError) as exc:
        raise RecordError(f"malformed fair evidence: {exc}") from exc


def make_fair_record(**fields):
    return validate_fair_record({"schema_version": 2, **fields})


def write_fair_record(rec, path):
    validate_fair_record(rec)
    raw = (json.dumps(rec, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    return _write_new(raw, path)


def _write_new(raw, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=path.parent, prefix=".record-")
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, path)
    finally:
        os.unlink(temporary)
    return hashlib.sha256(raw).hexdigest()


def write_record_hashes(paths, index_path):
    """Finalize an append-only hash inventory outside the record-directory glob."""
    index_path = Path(index_path)
    rows, identities = [], set()
    for path in paths:
        path = Path(path)
        _require(path.parent.resolve() != index_path.parent.resolve(), "hash inventory must be in an evidence subdirectory")
        raw = path.read_bytes()
        rec = read_record(path, expected_sha256=hashlib.sha256(raw).hexdigest())
        validate_fair_record(rec)
        identity = (rec["campaign_id"], rec["repeat_id"])
        _require(identity not in identities, "duplicate record identity")
        identities.add(identity)
        rows.append({"path": os.path.relpath(path.resolve(), index_path.parent.resolve()),
                     "sha256": hashlib.sha256(raw).hexdigest(), "campaign_id": identity[0], "repeat_id": identity[1]})
    return _write_new((json.dumps({"artifact": "record-hashes-v1", "records": rows}, sort_keys=True, allow_nan=False) + "\n").encode(), index_path)


def read_record_hashes(index_path):
    index_path = Path(index_path)
    try:
        index = json.loads(index_path.read_text(encoding="utf-8"))
        _require(index["artifact"] == "record-hashes-v1", "unsupported hash inventory")
        records, identities = [], set()
        for row in index["records"]:
            _require(_hash(row["sha256"]), "invalid record hash")
            rec = read_record(index_path.parent / row["path"], expected_sha256=row["sha256"])
            validate_fair_record(rec)
            identity = (rec["campaign_id"], rec["repeat_id"])
            _require(identity == (row["campaign_id"], row["repeat_id"]) and identity not in identities, "record hash identity mismatch")
            identities.add(identity)
            records.append(rec)
        return records
    except (KeyError, TypeError, ValueError) as exc:
        raise RecordError(f"invalid record hash inventory: {exc}") from exc


def _normalize_argv(argv):
    """Copy of argv with the value after each "--output" replaced by a placeholder."""
    out = list(argv)
    for i in range(len(out) - 1):
        if argv[i] == "--output":
            out[i + 1] = "<output>"
    return out


def check_fair_coverage(records, expected):
    """Check matrix identity; mapping values additionally declare exact attempt IDs."""
    actual, ids, orders = {}, set(), set()
    bindings, indices = {}, set()
    for rec in records:
        validate_fair_record(rec)
        identity = (rec["campaign_id"], rec["repeat_id"])
        order = (rec["campaign_id"], rec["execution_order"])
        if identity in ids or order in orders:
            raise RecordError("duplicate attempt ID or execution order")
        ids.add(identity)
        orders.add(order)
        key = tuple(rec[k] for k in ("campaign_id", "cell_id", "series_id", "contract_id", "backend"))
        view = dict(rec, argv=_normalize_argv(rec["argv"]))
        binding = tuple(json.dumps(view[k], sort_keys=True) for k in ("manifest_sha256", "reference_id", "tool", "environment_id", "tool_version", "operation", "host_id", "host_allocation", "argv", "options"))
        _require(key not in bindings or bindings[key] == binding, "cell identity collision; cannot pool different evidence")
        bindings[key] = binding
        index = (key, rec["repeat_index"])
        _require(index not in indices, "duplicate repeat index")
        indices.add(index)
        actual.setdefault(key, []).append(rec["repeat_id"])
    missing = set(expected) - actual.keys()
    if missing:
        raise MissingCellError(f"missing fair coverage: {sorted(missing)}")
    if isinstance(expected, dict):
        for key, attempts in expected.items():
            if actual[key] != list(attempts):
                raise MissingCellError(f"missing or out-of-order attempts for {key}")


def fair_eligible(record):
    try:
        validate_record(record)
    except (RecordError, TypeError, KeyError):
        return False
    return record.get("schema_version") == 2 and record.get("publication_eligible") is True
