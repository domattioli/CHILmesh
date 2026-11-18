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


def test_fair_record_rejects_unsupported_schema():
    from benchmarks.competitors import records as records_module
    with pytest.raises(records_module.RecordError):
        records_module.validate_fair_record({"schema_version": 1})


def fair_data():
    return dict(campaign_id="c", manifest_sha256="a" * 64, cell_id="cell", series_id="s",
                repeat_id="r", repeat_index=0, execution_order=0, reference_id="mesh", contract_id="parse",
                tool="chilmesh", environment_id="env", tool_version="1", chilmesh_git_sha="abc",
                argv=["python"], thread_environment={k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
                host_id="host", host_allocation={}, backend="python", outcome="not-comparable", reason="unsupported",
                equivalence="not-comparable", verification="not-applicable", measurement_validity={},
                publication_eligible=False, protocol={}, measurement={"operation_elapsed_s": None, "missing": {"operation_elapsed_s": "not run"}},
                cache={}, load={}, output_verification={}, termination={})


def test_v2_roundtrip_hash_and_collision(tmp_path):
    rec = R.make_fair_record(**fair_data())
    path = tmp_path / "r.json"
    digest = R.write_fair_record(rec, path)
    assert R.read_record(path, expected_sha256=digest) == rec
    with pytest.raises(FileExistsError): R.write_fair_record(rec, path)
    with pytest.raises(R.RecordError): R.read_record(path, expected_sha256="0" * 64)


@pytest.mark.parametrize("field", ["campaign_id", "cell_id", "series_id", "repeat_id", "reference_id", "contract_id", "environment_id", "protocol", "cache", "load", "termination"])
def test_v2_requires_identity_and_evidence(field):
    data = fair_data(); del data[field]
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


def test_v2_rejects_null_without_missing_reason():
    data = fair_data(); data["measurement"] = {"operation_elapsed_s": None}
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


def test_v2_rejects_nonfinite_json_but_accepts_explicit_infinite_spread():
    data = fair_data(); data["measurement"]["spread"] = float("inf")
    with pytest.raises(R.RecordError): R.make_fair_record(**data)
    data["measurement"]["spread"] = None
    data["measurement"]["missing"]["spread"] = "infinite"
    data["measurement"]["spread_state"] = "infinite"
    assert R.make_fair_record(**data)["measurement"]["spread_state"] == "infinite"


def test_v2_states_are_separate_and_not_comparable_is_explicit():
    data = fair_data(); data["equivalence"] = "like-for-like"
    with pytest.raises(R.RecordError): R.make_fair_record(**data)
    data = fair_data(); data["verification"] = "nonsense"
    with pytest.raises(R.RecordError): R.make_fair_record(**data)
    data = fair_data(); data["publication_eligible"] = True
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


def test_v2_coverage_includes_operation_backend_and_series():
    record = R.make_fair_record(**fair_data())
    R.check_fair_coverage([record], [("c", "cell", "s", "parse", "python")])
    for expected in [("c", "cell", "s", "write", "python"), ("c", "cell", "s", "parse", "cpp"), ("c", "cell", "other", "parse", "python")]:
        with pytest.raises(R.MissingCellError): R.check_fair_coverage([record], [expected])


def test_legacy_v1_roundtrip_never_fair_eligible(tmp_path):
    path = tmp_path / "legacy.json"
    R.write_record(_rec(), path)
    assert not R.fair_eligible(R.read_record(path))


def test_v2_cannot_be_created_with_another_schema():
    data = fair_data(); data["schema_version"] = 1
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


def test_invalid_record_cannot_grant_fair_eligibility():
    data = fair_data(); data["schema_version"] = 2; data["publication_eligible"] = True
    assert not R.fair_eligible(data)

# Independent inventory from RepeatRecord v2; absence is evidence, not a default.
EVIDENCE_FIELDS = {
    "protocol": "pid process_group_id ready_s start_s end_s ack_s clock receipt_timestamps new_modules materialization",
    "measurement": "operation_elapsed_s process_wall_s rss_baseline_bytes rss_operation_peak_bytes rss_whole_process_peak_bytes rss_growth_bytes last_rss_bytes wait4_lifetime_maxrss_bytes policy_version interval_s scope coverage gaps sample_count samples samples_path samples_sha256",
    "cache": "operations cached_before_kib cached_after_kib cached_before_bytes cached_after_bytes reset_utc reset_s launch_s success",
    "load": "pre_run_s interval_s observations allowlist detector_policy verdict host_id allocation",
    "output_verification": "n_nodes n_elems node_sha256 element_sha256 expected mapping_provenance normalization_revision reasons",
    "termination": "launch_utc exit_utc launch_s exit_s returncode signal limits trigger cgroup memory_events_before memory_events_after memory_events_delta memory_peak_bytes memory_peak_scope kill_diagnostics reap_diagnostics stdout_path stdout_sha256 stderr_path stderr_sha256 stdout_tail stderr_tail",
}
_original_fair_data = fair_data


def fair_data():
    data = _original_fair_data()
    data.update(operation="parse", options={}, created_at="2026-11-10T00:00:00Z")
    for section, names in EVIDENCE_FIELDS.items():
        data[section] = {name: None for name in names.split()}
        data[section]["missing"] = {name: "unsupported operation" for name in names.split()}
    data["measurement_validity"] = {key: {"valid": False, "reason": "unsupported operation"}
                                    for key in ("time", "memory", "cache", "load", "runtime")}
    return data


@pytest.mark.parametrize("section,field", [(s, f) for s, fields in EVIDENCE_FIELDS.items() for f in fields.split()])
def test_v2_every_evidence_field_is_explicit(section, field):
    data = fair_data()
    del data[section][field]
    with pytest.raises(R.RecordError):
        R.make_fair_record(**data)


@pytest.mark.parametrize("section", list(EVIDENCE_FIELDS))
def test_v2_nested_null_needs_missing_reason(section):
    data = fair_data(); data[section]["missing"] = {}
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


@pytest.mark.parametrize("decision", ["time", "memory", "cache", "load", "runtime"])
def test_v2_validity_requires_each_decision(decision):
    data = fair_data(); del data["measurement_validity"][decision]
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


@pytest.mark.parametrize("field", ["operation", "options", "created_at"])
def test_v2_requires_operation_options_and_utc(field):
    data = fair_data(); del data[field]
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


@pytest.mark.parametrize("change", [dict(outcome="ok", verification="verified"),
                                    dict(outcome="not-run", equivalence="like-for-like", publication_eligible=True),
                                    dict(outcome="failed", equivalence="like-for-like", verification="failed")])
def test_v2_states_cannot_bypass_execution_evidence(change):
    data = fair_data(); data.update(change)
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


def test_v2_record_id_collision_and_duplicate_coverage():
    rec = R.make_fair_record(**fair_data())
    with pytest.raises(R.RecordError):
        R.check_fair_coverage([rec, rec], [("c", "cell", "s", "parse", "python")])


def test_v2_coverage_validates_records_and_attempt_counts():
    rec = R.make_fair_record(**fair_data())
    key = ("c", "cell", "s", "parse", "python")
    with pytest.raises(R.MissingCellError): R.check_fair_coverage([rec], {key: ["r", "r2"]})
    bad = dict(rec, manifest_sha256="broken")
    with pytest.raises(R.RecordError): R.check_fair_coverage([bad], [key])


def test_record_artifact_subdirectories_are_not_legacy_records(tmp_path):
    R.write_fair_record(R.make_fair_record(**fair_data()), tmp_path / "repeat.json")
    (tmp_path / "samples").mkdir()
    (tmp_path / "samples" / "samples.json").write_text('{"samples": []}')
    assert len(R.read_records(tmp_path)) == 1


def successful_fair_data():
    data = fair_data()
    data.update(outcome="ok", reason="", equivalence="like-for-like", verification="verified", host_allocation={"cpuset": [0], "memory_bytes": 8589934592})
    data["measurement_validity"] = {k: {"valid": True, "reason": ""} for k in ("time", "memory", "cache", "load", "runtime")}
    data["protocol"] = dict(pid=42, process_group_id=42, ready_s=1.99, start_s=2., end_s=2.01, ack_s=2.02,
                            clock={"monotonic": True, "implementation": "clock_gettime(CLOCK_MONOTONIC)", "resolution_s": 1e-9},
                            receipt_timestamps={"READY": 1.99, "START": 2., "END": 2.01, "ACK": 2.02}, new_modules=[], materialization=True)
    samples = [dict(monotonic_s=t, pid=42, rss_bytes=rss, descendants=[]) for t, rss in [(1.99, 100), (2., 150), (2.01, 140)]]
    data["measurement"] = dict(operation_elapsed_s=2.01-2., process_wall_s=2.03-1.98, rss_baseline_bytes=100,
                               rss_operation_peak_bytes=150, rss_whole_process_peak_bytes=150, rss_growth_bytes=50,
                               last_rss_bytes=140, wait4_lifetime_maxrss_bytes=None, missing={"wait4_lifetime_maxrss_bytes": "unavailable on this host"},
                               policy_version="procfs-v1", interval_s=.010, scope="child-only", coverage=True, gaps=[.010, .010],
                               sample_count=3, samples=samples, samples_path="samples/r.json", samples_sha256="b" * 64)
    utc = "2026-11-10T00:00:00Z"
    data["cache"] = dict(operations=[dict(operation="sync", exit_code=0, monotonic_s=1., utc=utc),
                                     dict(operation="drop_caches", value="3", exit_code=0, monotonic_s=1.1, utc=utc)],
                         cached_before_kib=100, cached_after_kib=50, cached_before_bytes=102400, cached_after_bytes=51200,
                         reset_utc=utc, reset_s=1.2, launch_s=1.98, success=True)
    data["load"] = dict(pre_run_s=1, interval_s=.100, observations=[{"proc_stat": "cpu 0", "loadavg": "0", "diskstats": "0", "processes": []}],
                        allowlist=[], detector_policy={"cpu_percentage_points": 5, "io_bytes_per_s": 1048576}, verdict="clean",
                        host_id="host", allocation=data["host_allocation"])
    expected = dict(n_nodes=4, n_elems=2, node_sha256="c" * 64, element_sha256="d" * 64)
    data["output_verification"] = dict(expected, expected=expected.copy(), mapping_provenance={"node_ids": [10, 20, 40, 80], "element_ids": [7, 9]}, normalization_revision="binary64-repr-v1", reasons=[])
    data["termination"] = dict(launch_utc=utc, exit_utc=utc, launch_s=1.98, exit_s=2.03, returncode=0, signal=None,
                               limits={"wall_s": 600, "rss_bytes": 8589934592}, trigger="none", cgroup="/fake",
                               memory_events_before={"oom_kill": 0}, memory_events_after={"oom_kill": 0}, memory_events_delta={"oom_kill": 0},
                               memory_peak_bytes=None, memory_peak_scope=None, kill_diagnostics=[], reap_diagnostics={"reaped": True},
                               stdout_path="diagnostics/r.stdout", stdout_sha256="e" * 64, stderr_path="diagnostics/r.stderr", stderr_sha256="f" * 64,
                               stdout_tail="", stderr_tail="", missing={"signal": "normal exit", "memory_peak_bytes": "unavailable", "memory_peak_scope": "unavailable"})
    return data


def test_success_requires_actual_verified_output_and_independent_decisions(tmp_path):
    rec = R.make_fair_record(**successful_fair_data())
    digest = R.write_fair_record(rec, tmp_path / "r.json")
    assert R.read_record(tmp_path / "r.json", expected_sha256=digest) == rec
    assert not R.fair_eligible(rec)  # Stability and rights still belong to the series.


@pytest.mark.parametrize("section,field,value", [
    ("protocol", "new_modules", ["lazy_import"]), ("protocol", "materialization", False),
    ("protocol", "clock", {"monotonic": False}), ("protocol", "ack_s", 1.),
    ("measurement", "rss_growth_bytes", 999), ("measurement", "rss_whole_process_peak_bytes", 100),
    ("measurement", "rss_baseline_bytes", 90), ("measurement", "sample_count", 0),
    ("measurement", "samples", []), ("measurement", "coverage", False),
    ("measurement", "gaps", [.021]), ("measurement", "interval_s", .020),
    ("measurement", "operation_elapsed_s", 99.), ("measurement", "rss_growth_bytes", -1),
    ("cache", "success", False), ("cache", "cached_after_kib", 100),
    ("cache", "cached_after_bytes", 0), ("cache", "reset_s", 2.), ("cache", "operations", []),
    ("load", "verdict", "unknown"), ("load", "pre_run_s", .9), ("load", "allocation", {}),
    ("output_verification", "n_nodes", 0), ("output_verification", "node_sha256", "0"*64),
    ("output_verification", "mapping_provenance", {}), ("output_verification", "reasons", ["corruption"]),
    ("termination", "returncode", 1), ("termination", "trigger", "signal-unknown"),
])
def test_success_rejects_contradictory_evidence(section, field, value):
    data = successful_fair_data(); data[section][field] = value
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


@pytest.mark.parametrize("decision", ["time", "memory", "cache", "load", "runtime"])
def test_one_invalid_decision_prevents_ok(decision):
    data = successful_fair_data(); data["measurement_validity"][decision] = {"valid": False, "reason": "unavailable"}
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


def test_descendant_sample_invalidates_child_only_memory():
    data = successful_fair_data(); data["measurement"]["samples"][0]["descendants"] = [43]
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


def test_sample_missing_rss_is_explicit_and_invalid_for_ok():
    data = successful_fair_data(); data["measurement"]["samples"][0].update(rss_bytes=None, missing_reason="gone")
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


def test_kernel_oom_requires_counter_delta():
    data = successful_fair_data(); data.update(outcome="failed", reason="killed")
    data["termination"].update(trigger="kernel-oom", returncode=-9, signal=9)
    del data["termination"]["missing"]["signal"]
    with pytest.raises(R.RecordError): R.make_fair_record(**data)
    data["termination"]["memory_events_delta"]["oom_kill"] = 1
    data["termination"]["memory_events_after"]["oom_kill"] = 1
    assert R.make_fair_record(**data)["outcome"] == "failed"


def test_invalid_attempt_retains_returncode_and_partial_evidence():
    data = successful_fair_data(); data.update(outcome="invalid", reason="deferred import")
    data["protocol"]["new_modules"] = ["late"]
    data["measurement_validity"]["time"] = {"valid": False, "reason": "deferred import"}
    rec = R.make_fair_record(**data)
    assert rec["termination"]["returncode"] == 0
    assert rec["measurement"]["rss_operation_peak_bytes"] == 150


def test_failed_launch_has_no_fabricated_measurements():
    data = fair_data(); data.update(outcome="failed", equivalence="like-for-like", reason="missing executable")
    data["termination"]["trigger"] = "launch-error"
    del data["termination"]["missing"]["trigger"]
    assert R.make_fair_record(**data)["measurement"]["operation_elapsed_s"] is None


def test_atomic_collision_does_not_modify_original_or_leave_temporary(tmp_path):
    rec = R.make_fair_record(**fair_data()); path = tmp_path / "r.json"
    R.write_fair_record(rec, path); original = path.read_bytes()
    rec["reason"] = "second attempt"
    with pytest.raises(FileExistsError): R.write_fair_record(rec, path)
    assert path.read_bytes() == original
    assert list(tmp_path.iterdir()) == [path]


def test_record_hash_manifest_roundtrip_and_tamper_detection(tmp_path):
    rec = R.make_fair_record(**fair_data()); path = tmp_path / "r.json"
    R.write_fair_record(rec, path)
    index = tmp_path / "evidence" / "record-hashes.json"
    R.write_record_hashes([path], index)
    assert R.read_record_hashes(index) == [rec]
    path.write_text(path.read_text() + " ")
    with pytest.raises(R.RecordError): R.read_record_hashes(index)


@pytest.mark.parametrize("section,field,value", [
    ("measurement", "samples", [{"monotonic_s": 1.99, "pid": 42, "rss_bytes": None, "descendants": []}]),
    ("protocol", "pid", -1), ("termination", "returncode", "zero"),
    ("cache", "cached_before_kib", -1), ("load", "verdict", "made-up"),
    ("output_verification", "mapping_provenance", "source reparse"),
])
def test_invalid_attempt_still_validates_available_evidence(section, field, value):
    data = successful_fair_data(); data.update(outcome="invalid", reason="incomplete")
    data[section][field] = value
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


@pytest.mark.parametrize("section,field,value", [
    ("measurement", "process_wall_s", 100), ("measurement", "gaps", [0, 0]),
    ("protocol", "receipt_timestamps", {"READY": 0, "START": 0, "END": 0, "ACK": 0}),
    ("cache", "operations", [dict(operation="sync", exit_code=0, monotonic_s=1.1, utc="2026-11-10T00:00:00Z"), dict(operation="drop_caches", value="3", exit_code=0, monotonic_s=1., utc="2026-11-10T00:00:00Z")]),
    ("load", "detector_policy", {"cpu_percentage_points": 99, "io_bytes_per_s": 999999999}),
    ("termination", "limits", {"wall_s": 1, "rss_bytes": 1}),
    ("termination", "memory_events_delta", {"oom_kill": 3}),
])
def test_success_evidence_must_match_frozen_definitions(section, field, value):
    data = successful_fair_data(); data[section][field] = value
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


@pytest.mark.parametrize("state", ["failed", "invalid", "not-run", "not-comparable"])
def test_non_ok_states_require_reason(state):
    data = fair_data(); data.update(outcome=state, reason=" ")
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


def test_series_stability_cannot_be_copied_onto_repeat():
    data = fair_data(); data["stable"] = True
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


def test_strict_json_nan_is_rejected_on_read(tmp_path):
    path = tmp_path / "r.json"; path.write_text('{"schema_version": 2, "value": NaN}')
    with pytest.raises(R.RecordError): R.read_record(path)


def test_cell_cannot_pool_different_environment_input_or_host():
    rec = R.make_fair_record(**fair_data())
    other = dict(rec, repeat_id="r2", repeat_index=1, execution_order=1, environment_id="different")
    with pytest.raises(R.RecordError): R.check_fair_coverage([rec, other], [("c", "cell", "s", "parse", "python")])


def test_cell_cannot_repeat_same_index_with_another_attempt_id():
    rec = R.make_fair_record(**fair_data())
    other = dict(rec, repeat_id="r2", execution_order=1)
    with pytest.raises(R.RecordError): R.check_fair_coverage([rec, other], [("c", "cell", "s", "parse", "python")])


def test_hash_inventory_is_append_only_outside_legacy_glob(tmp_path):
    rec = R.make_fair_record(**fair_data()); path = tmp_path / "r.json"
    R.write_fair_record(rec, path)
    with pytest.raises(R.RecordError): R.write_record_hashes([path], tmp_path / "index.json")
    index = tmp_path / "evidence" / "index.json"
    R.write_record_hashes([path], index)
    with pytest.raises(FileExistsError): R.write_record_hashes([path], index)
    assert R.read_records(tmp_path) == [rec]


def test_cell_pools_repeats_differing_only_in_output_path():
    base = fair_data()
    a = R.make_fair_record(**dict(base, argv=["x", "--output", "raw/r/output.json"]))
    b = dict(a, repeat_id="r2", repeat_index=1, execution_order=1, argv=["x", "--output", "raw/r2/output.json"])
    R.check_fair_coverage([a, b], [("c", "cell", "s", "parse", "python")])
    assert a["argv"][2] == "raw/r/output.json"
    c = dict(b, repeat_id="r3", repeat_index=2, execution_order=2, argv=["y", "--output", "raw/r3/output.json"])
    with pytest.raises(R.RecordError): R.check_fair_coverage([a, c], [("c", "cell", "s", "parse", "python")])
