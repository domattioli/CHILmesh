import copy
import pytest
from benchmarks.competitors import campaign as C


def runtime():
    return dict(version="3.11.9", implementation="cpython", platform="linux", architecture="x86_64", realpath="/python", sha256="a"*64, executable="/venv/bin/python",
                implementation_fields={"name": "cpython", "cache_tag": "cpython-311", "version": [3, 11, 9, "final", 0], "hexversion": 51054832},
                thread_environment=dict(C.THREAD_ENV))


def test_schedule_budgets_and_frozen_manifest(fair_fixture):
    manifest = C.make_manifest("c", fair_fixture["meshes"], ["parse"], fair_fixture["tools"])
    assert len(manifest["schedule"]) == (5*3+3)*6
    assert manifest["limits_by_mesh_class"]["small"] == {"wall_s": 600, "rss_bytes": 8*1024**3}
    assert manifest["limits_by_mesh_class"]["stofs"] == {"wall_s": 1800, "rss_bytes": 48*1024**3}
    C.validate_manifest(manifest)
    manifest["schedule"].reverse()
    with pytest.raises(C.CampaignError): C.validate_manifest(manifest)


def test_preflight_runtime_pins_threads_and_lineage():
    env = locked_environment()
    C.validate_preflight(runtime(), [runtime()], [env], C.THREAD_ENV)
    bad = runtime(); bad["sha256"] = "c"*64
    with pytest.raises(C.CampaignError): C.validate_preflight(runtime(), [bad], [env], C.THREAD_ENV)
    badenv = copy.deepcopy(env); badenv["native"] = {}
    with pytest.raises(C.CampaignError): C.validate_preflight(runtime(), [runtime()], [badenv], C.THREAD_ENV)
    with pytest.raises(C.CampaignError): C.validate_preflight(runtime(), [runtime()], [env], {})
    with pytest.raises(C.CampaignError): C.validate_replacement("old", "")


def locked_environment():
    return dict(tool_distribution="tool", environment_id="env", lock_path="locks/tool.lock",
                converter={"required": False, "reason": "native input"}, executable="/venv/bin/python", source_revision="abc",
                version="1", distributions={"tool": "1", "dependency": "2"},
                direct_dependencies=["tool"], transitive_dependencies=["dependency"],
                native={"lib": "1"}, lock_sha256="b" * 64, compatible=True)


@pytest.mark.parametrize("key", ["version", "implementation", "platform", "architecture", "realpath", "sha256"])
def test_parent_and_child_build_identity_must_match(key):
    child = runtime()
    child[key] = "different"
    with pytest.raises(C.CampaignError):
        C.validate_preflight(runtime(), [child], [locked_environment()], C.THREAD_ENV)


def test_parent_must_be_python_311():
    parent = runtime(); parent["version"] = "3.12.0"
    with pytest.raises(C.CampaignError):
        C.validate_preflight(parent, [parent], [locked_environment()], C.THREAD_ENV)


def test_distribution_version_uses_metadata(monkeypatch):
    monkeypatch.setattr(C.importlib.metadata, "version", lambda name: {"tool": "4.2"}[name])
    assert C.distribution_version("tool") == "4.2"


@pytest.mark.parametrize("field", ["direct_dependencies", "transitive_dependencies", "native"])
def test_declared_pins_must_be_present(field):
    env = locked_environment()
    if field == "native":
        env[field] = {}
    else:
        env[field].append("missing")
    with pytest.raises(C.CampaignError):
        C.validate_preflight(runtime(), [runtime()], [env], C.THREAD_ENV)


def test_manifest_schedule_is_deterministic_serial_and_repeat_counts(fair_fixture):
    a = C.make_manifest("c", fair_fixture["meshes"], ["parse", "write"], fair_fixture["tools"])
    assert a == C.make_manifest("c", fair_fixture["meshes"], ["parse", "write"], fair_fixture["tools"])
    for mesh in fair_fixture["meshes"]:
        entries = [x for x in a["schedule"] if x["mesh"] == mesh and x["operation"] == "parse"]
        count = 3 if mesh == "STOFS" else 5
        assert [(x["repeat_index"], x["tool"]) for x in entries] == [(i, t) for i in range(count) for t in fair_fixture["tools"]]


def test_replacement_requires_link_and_rationale():
    for previous, reason in [("", "reason"), ("old", " ")]:
        with pytest.raises(C.CampaignError):
            C.validate_replacement(previous, reason)
    assert C.validate_replacement("old", "imports deferred")


@pytest.mark.parametrize("field", ["direct_dependencies", "transitive_dependencies"])
def test_missing_dependency_inventory_is_not_a_complete_lock(field):
    env = locked_environment(); del env[field]
    with pytest.raises(C.CampaignError):
        C.validate_preflight(runtime(), [runtime()], [env], C.THREAD_ENV)


def test_invalid_digest_is_not_runtime_evidence():
    parent = runtime(); parent["sha256"] = "not-a-digest"
    with pytest.raises(C.CampaignError):
        C.validate_preflight(parent, [parent], [locked_environment()], C.THREAD_ENV)


@pytest.mark.parametrize("digest", ["z" * 64, "a" * 63, ""])
def test_invalid_lock_digest_is_rejected(digest):
    env = locked_environment(); env["lock_sha256"] = digest
    with pytest.raises(C.CampaignError):
        C.validate_preflight(runtime(), [runtime()], [env], C.THREAD_ENV)


def test_empty_direct_dependency_inventory_is_rejected():
    env = locked_environment(); env["direct_dependencies"] = []
    with pytest.raises(C.CampaignError):
        C.validate_preflight(runtime(), [runtime()], [env], C.THREAD_ENV)


def test_noncanonical_interpreter_path_is_rejected():
    parent = runtime(); parent["realpath"] = "/path/../python"
    with pytest.raises(C.CampaignError):
        C.validate_preflight(parent, [parent], [locked_environment()], C.THREAD_ENV)


@pytest.mark.parametrize("field", ["created_at", "chilmesh_git_sha", "host_id", "host_manifest", "runtime_manifest", "environment_manifests", "sampler_policy", "load_policy", "cache_policy", "canonicalization_revision", "ordered_cells"])
def test_manifest_requires_every_frozen_evidence_identity(fair_fixture, field):
    manifest = C.make_manifest("c", fair_fixture["meshes"], ["parse"], fair_fixture["tools"])
    assert field in manifest
    del manifest[field]
    with pytest.raises(C.CampaignError): C.validate_manifest(manifest)


@pytest.mark.parametrize("field", ["tool_distribution", "environment_id", "lock_path", "converter", "executable", "source_revision"])
def test_environment_requires_inventory_identity(field):
    env = locked_environment(); del env[field]
    with pytest.raises(C.CampaignError): C.validate_preflight(runtime(), [runtime()], [env], C.THREAD_ENV)


@pytest.mark.parametrize("change", [dict(distributions={"tool": "1", "dependency": "2", "undeclared": "3"}),
                                    dict(distributions={"tool": "1", "dependency": ""}),
                                    dict(version="wrong"), dict(compatible="yes"),
                                    dict(native={"lib": "unknown"})])
def test_incomplete_or_inconsistent_lock_is_rejected(change):
    env = locked_environment(); env.update(change)
    with pytest.raises(C.CampaignError): C.validate_preflight(runtime(), [runtime()], [env], C.THREAD_ENV)


def test_converter_lock_is_required_and_pinned_when_used():
    env = locked_environment(); env["converter"] = {"required": True}
    with pytest.raises(C.CampaignError): C.validate_preflight(runtime(), [runtime()], [env], C.THREAD_ENV)


def test_child_thread_policy_is_checked():
    child = runtime(); child["thread_environment"] = dict(C.THREAD_ENV, OMP_NUM_THREADS="2")
    with pytest.raises(C.CampaignError): C.validate_preflight(runtime(), [child], [locked_environment()], C.THREAD_ENV)


def test_missing_child_thread_policy_is_not_evidence():
    child = runtime(); del child["thread_environment"]
    with pytest.raises(C.CampaignError): C.validate_preflight(runtime(), [child], [locked_environment()], C.THREAD_ENV)


def test_backend_cells_and_replacement_lineage_are_frozen(fair_fixture):
    manifest = C.make_manifest("c", ["STOFS"], [("adjacency", "cpp"), ("adjacency", "python")], ["chilmesh"],
                               previous_series_id="old", replacement_rationale="deferred import")
    assert [x["backend"] for x in manifest["ordered_cells"]] == ["cpp", "python"]
    assert manifest["previous_series_id"] == "old"
    C.validate_manifest(manifest)


def test_manifest_source_collections_cannot_mutate_frozen_content(fair_fixture):
    meshes = fair_fixture["meshes"]
    manifest = C.make_manifest("c", meshes, ["parse"], ["chilmesh"])
    meshes.append("another")
    assert "another" not in manifest["meshes"]
    C.validate_manifest(manifest)


def test_manifest_self_consistent_hash_cannot_change_policy(fair_fixture):
    manifest = C.make_manifest("c", fair_fixture["meshes"], ["parse"], ["chilmesh"])
    manifest["repeats_by_class"]["small"] = 1
    manifest["manifest_sha256"] = C._digest({k: v for k, v in manifest.items() if k != "manifest_sha256"})
    with pytest.raises(C.CampaignError): C.validate_manifest(manifest)


def test_incomplete_frozen_proposal_cannot_pass_execution_preflight(fair_manifest):
    C.validate_manifest(fair_manifest)
    with pytest.raises(C.CampaignError): C.validate_manifest(fair_manifest, require_ready=True)


@pytest.mark.parametrize("field", ["implementation_fields", "executable", "thread_environment"])
def test_runtime_requires_build_executable_and_thread_evidence(field):
    parent = runtime(); del parent[field]
    with pytest.raises(C.CampaignError): C.validate_preflight(parent, [runtime()], [locked_environment()], C.THREAD_ENV)


def test_runtime_collector_hashes_real_binary_and_collects_full_build(monkeypatch, tmp_path):
    import hashlib
    binary = tmp_path / "python"; binary.write_bytes(b"synthetic interpreter")
    link = tmp_path / "venv-python"; link.symlink_to(binary)
    monkeypatch.setattr(C.sys, "executable", str(link))
    for name, value in C.THREAD_ENV.items(): monkeypatch.setenv(name, value)
    identity = C.runtime_identity()
    assert identity["realpath"] == str(binary)
    assert identity["sha256"] == hashlib.sha256(binary.read_bytes()).hexdigest()
    assert identity["version"] == C.sys.version
    assert identity["implementation_fields"]["cache_tag"] == C.sys.implementation.cache_tag
    assert identity["thread_environment"] == C.THREAD_ENV


def test_fixture_is_synthetic_isolated_and_offline(request, monkeypatch):
    import os
    import socket
    def forbidden(*args, **kwargs): raise AssertionError("network forbidden")
    import subprocess
    monkeypatch.setattr(socket, "create_connection", forbidden)
    monkeypatch.setattr(subprocess, "run", forbidden)
    monkeypatch.setattr(os, "system", forbidden)
    fair_fixture = request.getfixturevalue("fair_fixture")
    lines = fair_fixture["mesh"].read_text().splitlines()
    assert lines[1] == "2 4"
    assert [int(line.split()[0]) for line in lines[2:6]] == [10, 20, 40, 80]
    assert [int(line.split()[1]) for line in lines[6:]] == [3, 4]
    assert len(fair_fixture["meshes"]) == 4 and len(fair_fixture["tools"]) == 6
    assert os.environ["HOME"].startswith(str(fair_fixture["output"].parent))
    assert os.environ["XDG_CACHE_HOME"].startswith(os.environ["HOME"].rsplit('/', 1)[0])
    assert fair_fixture["clock"]() < fair_fixture["clock"]()
    assert (fair_fixture["proc"] / "meminfo").read_text() == "Cached: 100 kB\n"
    assert (fair_fixture["cgroup"] / "memory.events").read_text() == "oom_kill 0\n"
    assert fair_fixture["adapter"]["node_ids"] == [10, 20, 40, 80]


@pytest.mark.parametrize("value", ["unknown", ">=1", "1.*", ""])
def test_exact_dependency_pin_rejects_unknown_or_ranges(value):
    env = locked_environment(); env["distributions"]["dependency"] = value
    with pytest.raises(C.CampaignError): C.validate_preflight(runtime(), [runtime()], [env], C.THREAD_ENV)


def test_dependency_version_collection_uses_distribution_metadata_for_every_pin(monkeypatch):
    called = []
    def version(name):
        called.append(name)
        return {"tool": "1", "dependency": "2"}[name]
    monkeypatch.setattr(C.importlib.metadata, "version", version)
    assert C.distribution_versions(["tool", "dependency"]) == {"tool": "1", "dependency": "2"}
    assert called == ["tool", "dependency"]


def test_missing_installed_distribution_is_explicit_preflight_failure(monkeypatch):
    def absent(name): raise C.importlib.metadata.PackageNotFoundError(name)
    monkeypatch.setattr(C.importlib.metadata, "version", absent)
    with pytest.raises(C.CampaignError): C.distribution_versions(["missing"])


def ready_manifest():
    host = dict(host_id="host", cpu_model="synthetic", logical_count=1, affinity=[0], available_cores=1,
                physical_ram_bytes=16*1024**3, kernel="Linux synthetic", os="Linux", cgroup_path="/fake", cgroup_mount="/fake",
                cpu_quota={"raw": "max 100000"}, cpuset=[0], memory_limits={"raw": "8589934592"}, raw_evidence={"boot": "hashed", "container": "hashed"})
    return C.make_manifest("c", ["STOFS"], [("parse", "python")], ["tool"], created_at="2026-11-10T00:00:00Z",
                           chilmesh_git_sha="abc", host_id="host", host_manifest=host, runtime_manifest=runtime(), environment_manifests=[locked_environment()])


def rehash(manifest):
    manifest["manifest_sha256"] = C._digest({k: v for k, v in manifest.items() if k != "manifest_sha256"})
    return manifest


def test_complete_synthetic_manifest_is_ready_but_not_real_compatibility_evidence():
    manifest = ready_manifest()
    assert C.validate_manifest(manifest, require_ready=True) == manifest


@pytest.mark.parametrize("field", ["host_id", "cpu_model", "logical_count", "affinity", "available_cores", "physical_ram_bytes", "kernel", "os", "cgroup_path", "cgroup_mount", "cpu_quota", "cpuset", "memory_limits", "raw_evidence"])
def test_ready_manifest_requires_every_host_allocation_field(field):
    manifest = ready_manifest(); del manifest["host_manifest"][field]
    with pytest.raises(C.CampaignError): C.validate_manifest(rehash(manifest), require_ready=True)


@pytest.mark.parametrize("field,value", [("logical_count", 0), ("physical_ram_bytes", -1), ("affinity", []), ("available_cores", 0), ("raw_evidence", {})])
def test_ready_host_evidence_has_valid_values(field, value):
    manifest = ready_manifest(); manifest["host_manifest"][field] = value
    with pytest.raises(C.CampaignError): C.validate_manifest(rehash(manifest), require_ready=True)


def test_host_allocation_drift_requires_new_campaign():
    manifest = ready_manifest()
    changed = copy.deepcopy(manifest["host_manifest"]); changed["memory_limits"] = {"raw": "4294967296"}
    with pytest.raises(C.CampaignError): C.validate_continuation(manifest, changed, runtime())
    assert C.validate_continuation(manifest, manifest["host_manifest"], runtime())


def test_new_campaign_can_retain_previous_campaign_identity():
    manifest = C.make_manifest("new", ["STOFS"], ["parse"], ["tool"], previous_campaign_id="old")
    assert C.validate_manifest(manifest)["previous_campaign_id"] == "old"


@pytest.mark.parametrize("change", [{"name": "pypy"}, {"version": [3, 12, 0, "final", 0]}, {"cache_tag": ""}])
def test_parent_build_fields_must_agree_with_python_311(change):
    parent = runtime(); parent["implementation_fields"].update(change)
    with pytest.raises(C.CampaignError): C.validate_preflight(parent, [parent], [locked_environment()], C.THREAD_ENV)


def test_environment_executable_must_bind_corresponding_child():
    env = locked_environment(); env["executable"] = "/other/python"
    with pytest.raises(C.CampaignError): C.validate_preflight(runtime(), [runtime()], [env], C.THREAD_ENV)


def test_physical_lock_validation_rejects_missing_placeholder_and_hash_mismatch(tmp_path):
    path=tmp_path/'tool.lock';env=locked_environment();env['lock_path']=str(path)
    with pytest.raises(C.CampaignError): C.validate_lock_files([env])
    path.write_text('# compatibility not established\n')
    env['lock_sha256']=__import__('hashlib').sha256(path.read_bytes()).hexdigest()
    with pytest.raises(C.CampaignError): C.validate_lock_files([env])
    path.write_text('tool==1\ndependency==2\n')
    with pytest.raises(C.CampaignError): C.validate_lock_files([env])
    env['lock_sha256']=__import__('hashlib').sha256(path.read_bytes()).hexdigest()
    assert C.validate_lock_files([env])


def test_serial_executor_keeps_order_append_only_and_replacement(tmp_path):
    from test_records import successful_fair_data
    manifest=ready_manifest();seen=[]
    def execute(entry):
        seen.append(entry['repeat_id']);r=successful_fair_data()
        r.update(campaign_id=manifest['campaign_id'],manifest_sha256=manifest['manifest_sha256'],
                 cell_id=entry['cell_id'],series_id=entry['series_id'],repeat_id=entry['repeat_id'],
                 repeat_index=entry['repeat_index'],execution_order=entry['execution_order'],tool=entry['tool'],backend=entry['backend'])
        return dict(schema_version=2,**r)
    rows=C.execute_schedule(manifest,tmp_path,execute=execute,host_manifest=manifest['host_manifest'],runtime_manifest=runtime(),lock_check=lambda env:True,lock_path=tmp_path/'host.lock',cgroup_check=lambda p:{'memory_events':{'oom_kill':0}})
    assert seen==[e['repeat_id'] for e in manifest['schedule']]
    assert len(rows)==3 and len(list((tmp_path/'attempts').glob('*.json')))==3
    with pytest.raises(FileExistsError):
        C.execute_schedule(manifest,tmp_path,execute=execute,host_manifest=manifest['host_manifest'],runtime_manifest=runtime(),lock_check=lambda env:True,lock_path=tmp_path/'host.lock',cgroup_check=lambda p:{'memory_events':{'oom_kill':0}})
    assert len(seen)==3


def test_ready_manifest_supports_multiple_isolated_environments():
    manifest=ready_manifest()
    environments=[locked_environment(),locked_environment()]
    environments[1].update(environment_id='env-two',executable='/other-venv/bin/python')
    m=C.make_manifest('multi',['STOFS'],['parse'],['one','two'],
        **{k:manifest[k] for k in C.MANIFEST_EVIDENCE if k!='environment_manifests'},environment_manifests=environments)
    assert C.validate_manifest(m,require_ready=True)


def test_cgroup_unknown_blocks_execution_before_callback(tmp_path):
    manifest=ready_manifest()
    with pytest.raises(C.CampaignError,match='memory.events'):
        C.execute_schedule(manifest,tmp_path,execute=lambda e:pytest.fail('must not launch'),
          host_manifest=manifest['host_manifest'],runtime_manifest=runtime(),lock_check=lambda e:True,
          cgroup_check=lambda p:{'memory_events':None},lock_path=tmp_path/'host.lock')
