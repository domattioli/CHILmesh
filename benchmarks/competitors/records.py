"""Run-record schema for the competitor benchmark harness.

One JSON record per tool x mesh x workflow x track x repeat. Machine statuses
are ``ok | failed | not-run | invalid``. A metric that does not apply carries
``value: null`` and ``na: true`` with a justification; the report labels it
"not applicable". The status ``not-run`` is displayed as "not run".
"""
from __future__ import annotations

import json
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
    validate_record(rec)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")
    return path


def read_record(path) -> dict:
    return validate_record(json.loads(Path(path).read_text(encoding="utf-8")))


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
