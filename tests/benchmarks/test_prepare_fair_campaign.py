"""Offline tests for the fair-campaign preparer: synthetic hosts, locks and Valence stubs only."""
import hashlib
import json
from argparse import Namespace
from pathlib import Path

import pytest

from benchmarks.competitors import campaign
from benchmarks.competitors.scripts import prepare_fair_campaign as prep

TOOLS = prep.TOOLS


def _locks(tmp_path):
    d = tmp_path / "locks"
    d.mkdir()
    envs = []
    for tool in TOOLS:
        lock = d / f"requirements-{tool}-fair.lock"
        lock.write_text(f"{tool}==1.0\nnumpy==2.0.0\nnetCDF4==1.7.0\n")
        envs.append(dict(tool=tool, lock_path=str(lock), lock_sha256=hashlib.sha256(lock.read_bytes()).hexdigest(),
                         compatible=True))
    (d / "environments.json").write_text(json.dumps({"environments": envs}))
    return d, {"environments": envs}


def _host():
    return {**{k: "x" for k in campaign.HOST_FIELDS}, "host_id": "h" * 64, "logical_count": 16, "physical_ram_bytes": 1,
            "affinity": [0, 1], "available_cores": 2, "raw_evidence": {"a": 1}, "cpu_quota": {"raw": "max 100000"},
            "cpuset": "0-1", "memory_limits": {"max": "max", "high": "max"}}


def _runtime():
    return dict(version="3.11.18", implementation="cpython", implementation_fields={"name": "cpython", "cache_tag": "cpython-311", "version": [3, 11, 18], "hexversion": 1},
                thread_environment=dict(campaign.THREAD_ENV), platform="linux-x86_64", architecture="x86_64",
                executable="/usr/local/bin/python3.11", realpath="/usr/local/bin/python3.11", sha256="a" * 64)


def _args(tmp_path, **kw):
    d = dict(out_dir=str(tmp_path / "out"), campaign_id="c1", locks_dir=str(tmp_path / "locks"), mesh=["Small"],
             valence_meta=None, extra_mesh=None, extra_sha256=None, extra_note="")
    d.update(kw)
    return Namespace(**d)


def _fake_reference(monkeypatch, status="ok"):
    monkeypatch.setattr(prep.valence, "reference_from_fetch",
                        lambda name, v, a, f, out: {"status": status, "reference_id": "r-" + name, "raw_file_sha256": "b" * 64,
                                                    "verification_state": "verified", "reason": "nope"})
    monkeypatch.setattr(prep.campaign, "prepare_conversion", lambda ref, reader, out, **k: {"status": "failed", "reason": "stub"})


def _fetch_ok(name, work):
    return {"status": "ok", "path": f"/tmp/{name}.14", "sha256": "b" * 64, "bytes": 1}


def _build(tmp_path, monkeypatch, **kw):
    d, envs = _locks(tmp_path)
    _fake_reference(monkeypatch, kw.pop("status", "ok"))
    return prep.build(_args(tmp_path, **kw), host=_host(), runtime=_runtime(), fetch=kw.pop("fetch", _fetch_ok) if False else _fetch_ok,
                      environments=envs, runner=lambda *a, **k: type("R", (), {"returncode": 1, "stdout": ""})())


def test_parse_lock_reads_exact_pins(tmp_path):
    p = tmp_path / "l"
    p.write_text("# c\na==1\nb==2\n")
    assert prep.parse_lock(p) == {"a": "1", "b": "2"}


def test_environment_manifest_fields(tmp_path):
    d, envs = _locks(tmp_path)
    env = prep.environment_manifest("xarray", d, envs, runner=lambda *a, **k: type("R", (), {"returncode": 1, "stdout": ""})())
    assert env["executable"] == "/opt/fair/xarray/bin/python" and env["version"] == "1.0"
    assert set(env["direct_dependencies"]) | set(env["transitive_dependencies"]) == set(env["distributions"])
    assert env["converter"]["required"] is True and env["native"] == {}


def test_nonreader_env_has_no_converter(tmp_path):
    d, envs = _locks(tmp_path)
    env = prep.environment_manifest("adcircpy", d, envs, runner=lambda *a, **k: type("R", (), {"returncode": 1, "stdout": ""})())
    assert env["converter"]["required"] is False and env["converter"]["reason"]


def test_manifest_valid_and_missing_native_blocks_ready(tmp_path, monkeypatch):
    manifest, prepared = _build(tmp_path, monkeypatch)
    campaign.validate_manifest(manifest)
    assert manifest["meshes"] == ["Small"] and not manifest["missing"]
    with pytest.raises(campaign.CampaignError):
        campaign.validate_manifest(manifest, require_ready=True)  # native pins unestablished: left absent


def test_commands_only_for_comparable_cells(tmp_path, monkeypatch):
    manifest, prepared = _build(tmp_path, monkeypatch)
    by_id = {e["repeat_id"]: e for e in manifest["schedule"]}
    assert set(prepared["identities"]) == set(by_id)
    for rid, cmd in prepared["commands"].items():
        e = by_id[rid]
        assert e["operation"] != "parse-netcdf"  # conversions failed in this stub
        assert cmd[0] == f"/opt/fair/{e['tool']}/bin/python"
    assert not any(by_id[r]["tool"] in ("thalassa", "xarray") for r in prepared["commands"])


def test_repeat_counts_small(tmp_path, monkeypatch):
    manifest, _ = _build(tmp_path, monkeypatch)
    assert sum(1 for e in manifest["schedule"] if e["tool"] == "chilmesh" and e["operation"] == "parse") == 5


def test_unserved_mesh_excluded_with_reason(tmp_path, monkeypatch):
    with pytest.raises(SystemExit, match="nope"):
        _build(tmp_path, monkeypatch, status="not-run")


def test_extra_mesh_sha_mismatch_refused(tmp_path, monkeypatch):
    f = tmp_path / "m.14"
    f.write_text("x")
    with pytest.raises(SystemExit):
        _build(tmp_path, monkeypatch, extra_mesh=f"EC2015={f}", extra_sha256="0" * 64)


def test_extra_mesh_registered_internal_only(tmp_path, monkeypatch):
    f = tmp_path / "m.14"
    f.write_text("x")
    sha = hashlib.sha256(b"x").hexdigest()
    manifest, prepared = _build(tmp_path, monkeypatch, extra_mesh=f"EC2015={f}", extra_sha256=sha, extra_note="EC2015")
    assert "EC2015" in manifest["meshes"] and prepared["notes"]["EC2015"].startswith(prep.RIGHTS_MARKER)


def test_token_never_written(tmp_path, monkeypatch):
    monkeypatch.setenv("VALENCE_TOKEN", "SECRET-TOKEN-VALUE")
    manifest, prepared = _build(tmp_path, monkeypatch)
    assert "SECRET-TOKEN-VALUE" not in json.dumps([manifest, prepared])
