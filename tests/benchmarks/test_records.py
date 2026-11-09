"""Run-record schema, status enum, na handling and coverage check (SC-001)."""
from __future__ import annotations

import pytest

from benchmarks.competitors import records as R


def _mesh(**kw):
    base = dict(name="synthetic", version="v0", sha256="0" * 64,
                provenance="synthetic test mesh", crs="cartesian", units="m",
                n_nodes=25, n_elems=32)
    base.update(kw)
    return base


def _rec(**kw):
    base = dict(tool="chilmesh", tool_version="1.0", license="PolyForm-NC",
                run_mode="subprocess", mesh=_mesh(), workflow="w2", track=None,
                metric_values={"parse_s": R.metric(1.5)}, status="ok", reason="",
                cache="cold", command=["python", "x"], wall_s=1.6,
                peak_rss_bytes=1000, repeat_index=0,
                machine={"os": "test"}, chilmesh_git_sha="abc123")
    base.update(kw)
    return R.make_record(**base)


def test_valid_record_has_schema_version_and_mesh_fields():
    rec = _rec()
    assert rec["schema_version"] == R.SCHEMA_VERSION
    for key in ("version", "provenance", "crs", "units", "sha256", "n_nodes", "n_elems"):
        assert key in rec["mesh"]
    R.validate_record(rec)


def test_status_enum():
    assert set(R.STATUSES) == {"ok", "failed", "not-run", "invalid"}
    with pytest.raises(R.RecordError):
        _rec(status="done")


def test_cache_enum():
    assert set(R.CACHE_STATES) == {"cold", "warm", "unknown"}
    with pytest.raises(R.RecordError):
        _rec(cache="lukewarm")


@pytest.mark.parametrize("missing", ["version", "provenance", "crs", "units"])
def test_missing_mesh_field_rejected(missing):
    mesh = _mesh()
    del mesh[missing]
    with pytest.raises(R.RecordError):
        _rec(mesh=mesh)


def test_failed_and_not_run_need_reason():
    with pytest.raises(R.RecordError):
        _rec(status="failed", reason="", metric_values={})
    with pytest.raises(R.RecordError):
        _rec(status="not-run", reason="  ", metric_values={})
    _rec(status="not-run", reason="mesh not yet served by Valence", metric_values={})


def test_na_metric_is_null_with_flag_and_justification():
    cell = R.metric(None, na=True, reason="adcircpy has no point-location API")
    assert cell == {"value": None, "na": True, "reason": "adcircpy has no point-location API"}
    with pytest.raises(R.RecordError):
        R.metric(None, na=True, reason="")
    with pytest.raises(R.RecordError):
        R.metric(3.0, na=True, reason="has a value")


def test_null_value_without_na_rejected():
    with pytest.raises(R.RecordError):
        _rec(metric_values={"parse_s": {"value": None, "na": False, "reason": ""}})


def test_coverage_passes_with_value_failure_or_justified_na():
    ok = _rec(tool="a", metric_values={"parse_s": R.metric(1.0), "adj_s": R.metric(2.0)})
    na = _rec(tool="b", metric_values={"parse_s": R.metric(1.0),
                                       "adj_s": R.metric(None, na=True, reason="no adjacency API")})
    failed = _rec(tool="c", status="failed", reason="OOM", metric_values={})
    R.check_coverage([ok, na, failed], tools=["a", "b", "c"], metrics=["parse_s", "adj_s"])


def test_coverage_fails_loudly_on_missing_cell():
    partial = _rec(tool="a", metric_values={"parse_s": R.metric(1.0)})
    with pytest.raises(R.MissingCellError) as exc:
        R.check_coverage([partial], tools=["a"], metrics=["parse_s", "adj_s"])
    assert "a" in str(exc.value) and "adj_s" in str(exc.value)


def test_coverage_fails_for_tool_with_no_record():
    rec = _rec(tool="a")
    with pytest.raises(R.MissingCellError):
        R.check_coverage([rec], tools=["a", "ghost"], metrics=["parse_s"])


def test_record_roundtrips_through_json(tmp_path):
    rec = _rec()
    path = tmp_path / "r.json"
    R.write_record(rec, path)
    assert R.read_record(path) == rec
