"""Frozen, offline campaign contracts; execution is a separate phase.

Dependency inventory shape (data-model.md Environment manifest): distributions
maps distribution names to exact metadata versions; direct_dependencies and
transitive_dependencies are disjoint name lists whose union is that mapping.
native maps native library/runtime names to exact versions. tool_distribution,
version, environment_id, lock_path/SHA-256, executable and source_revision bind
the inventory. converter is {required: False, reason: ...} for native input or
{required: True, revision, sha256, lock_sha256, native, options} for conversion.
compatible is a boolean recording a later install/import check; offline helpers
never infer it. Empty inventories/unknown versions block execution preflight.

make_manifest can freeze an incomplete proposal with explicit missing reasons;
validate_manifest(require_ready=True) rejects it before execution. Policies and
ordering are immutable by digest, rather than by Python object mutability.
"""
import copy
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import sys
import sysconfig
from datetime import datetime

THREAD_ENV = {key: "1" for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")}
LIMITS = {"small": {"wall_s": 600, "rss_bytes": 8 * 1024**3}, "stofs": {"wall_s": 1800, "rss_bytes": 48 * 1024**3}}
REPEATS = {"small": 5, "stofs": 3}


class CampaignError(ValueError):
    """A frozen campaign or its prerequisites are invalid."""


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


SAMPLER_POLICY = {"version": "procfs-v1", "interval_s": .010, "max_gap_s": .020, "scope": "child-only", "unit": "bytes"}
LOAD_POLICY = {"interval_s": .100, "pre_run_s": 1, "cpu_percentage_points": 5, "io_bytes_per_s": 1048576, "allowlist": []}
CACHE_POLICY = {"reset": "sync/drop_caches=3", "before_launch": True, "require_cached_decrease": True}
MANIFEST_EVIDENCE = ("created_at", "chilmesh_git_sha", "host_id", "host_manifest", "runtime_manifest", "environment_manifests")
HOST_FIELDS = ("host_id", "cpu_model", "logical_count", "affinity", "available_cores", "physical_ram_bytes", "kernel", "os", "cgroup_path", "cgroup_mount", "cpu_quota", "cpuset", "memory_limits", "raw_evidence")


def _require(condition, message):
    if not condition:
        raise CampaignError(message)


def _hash(value):
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def _text(value):
    return isinstance(value, str) and bool(value.strip()) and value.lower() not in ("unknown", "missing", "unavailable")


def _pin(value):
    return _text(value) and not any(c in value for c in "*<>=!~,;")


def _utc(value):
    try:
        return isinstance(value, str) and datetime.fromisoformat(value.replace("Z", "+00:00")).utcoffset().total_seconds() == 0
    except (ValueError, AttributeError):
        return False


def _operation(value):
    if isinstance(value, str):
        _require(_text(value), "empty operation")
        return value, "default"
    _require(isinstance(value, (list, tuple)) and len(value) == 2 and all(_text(v) for v in value), "invalid operation/backend")
    return tuple(value)


def _layout(meshes, operations, tools, series_id):
    cells, schedule = [], []
    for mesh in meshes:
        mesh_class = "stofs" if "stofs" in mesh.lower() else "small"
        for declared in operations:
            operation, backend = _operation(declared)
            for tool in tools:
                cell = dict(mesh=mesh, operation=operation, backend=backend, tool=tool, series_id=series_id)
                cell["cell_id"] = _digest(cell)
                cells.append(cell)
            for repeat in range(REPEATS[mesh_class]):
                for cell in cells[-len(tools):]:
                    entry = dict(cell, repeat_index=repeat, execution_order=len(schedule))
                    entry["repeat_id"] = _digest(entry)
                    schedule.append(entry)
    return cells, schedule


def make_manifest(campaign_id, meshes, operations, tools, *, series_id="original",
                  previous_series_id=None, replacement_rationale=None, **evidence):
    _require(_text(campaign_id) and _text(series_id), "missing campaign/series identity")
    meshes, operations, tools = (copy.deepcopy(list(v)) for v in (meshes, operations, tools))
    for values in (meshes, tools):
        _require(bool(values) and all(_text(v) for v in values) and len(set(values)) == len(values), "missing/duplicate inventory")
    normalized = [_operation(op) for op in operations]
    _require(bool(normalized) and len(set(normalized)) == len(normalized), "missing/duplicate operation/backend")
    _require(set(evidence) <= set(MANIFEST_EVIDENCE) | {"previous_campaign_id", "canonicalization_revision"}, "unknown manifest evidence")
    if previous_series_id is not None or replacement_rationale is not None:
        validate_replacement(previous_series_id, replacement_rationale)
        _require(previous_series_id != series_id, "replacement cannot supersede itself")
    cells, schedule = _layout(meshes, operations, tools, series_id)
    manifest = dict(schema_version=2, campaign_id=campaign_id, meshes=meshes, operations=operations,
                    tools=tools, ordered_cells=cells, schedule=schedule, series_id=series_id,
                    previous_series_id=previous_series_id, replacement_rationale=replacement_rationale,
                    previous_campaign_id=evidence.get("previous_campaign_id"),
                    limits_by_mesh_class=copy.deepcopy(LIMITS), repeats_by_class=dict(REPEATS),
                    thread_environment=dict(THREAD_ENV), sampler_policy=copy.deepcopy(SAMPLER_POLICY),
                    load_policy=copy.deepcopy(LOAD_POLICY), cache_policy=copy.deepcopy(CACHE_POLICY),
                    canonicalization_revision=evidence.get("canonicalization_revision", "binary64-repr-v1"))
    manifest["missing"] = {}
    for key in MANIFEST_EVIDENCE:
        manifest[key] = copy.deepcopy(evidence.get(key))
        if manifest[key] is None:
            manifest["missing"][key] = "not yet collected; blocks execution preflight"
    manifest["manifest_sha256"] = _digest(manifest)
    return manifest


def validate_manifest(manifest, *, require_ready=False):
    try:
        body = {k: v for k, v in manifest.items() if k != "manifest_sha256"}
        _require(_digest(body) == manifest["manifest_sha256"], "manifest changed")
        evidence = {k: manifest[k] for k in MANIFEST_EVIDENCE}
        evidence.update(previous_campaign_id=manifest["previous_campaign_id"], canonicalization_revision=manifest["canonicalization_revision"])
        expected = make_manifest(manifest["campaign_id"], manifest["meshes"], manifest["operations"], manifest["tools"],
                                 series_id=manifest["series_id"], previous_series_id=manifest["previous_series_id"],
                                 replacement_rationale=manifest["replacement_rationale"], **evidence)
        _require(manifest == expected, "manifest schedule or policy invalid")
        _require(_text(manifest["canonicalization_revision"]), "missing canonicalization revision")
        if manifest["created_at"] is not None:
            _require(_utc(manifest["created_at"]), "created_at must be UTC")
        if require_ready:
            _require(not manifest["missing"], "manifest missing execution evidence")
            host = manifest["host_manifest"]
            _require(isinstance(host, dict) and all(host.get(k) is not None for k in HOST_FIELDS), "incomplete host manifest")
            _require(host["host_id"] == manifest["host_id"] and _text(manifest["chilmesh_git_sha"]), "host/revision identity mismatch")
            for key in ("logical_count", "physical_ram_bytes"):
                _require(type(host[key]) is int and host[key] > 0, f"invalid host {key}")
            _require(type(host["available_cores"]) in (int, float) and host["available_cores"] > 0, "invalid available cores")
            _require(isinstance(host["affinity"], list) and bool(host["affinity"]) and all(type(cpu) is int and cpu >= 0 for cpu in host["affinity"]), "invalid affinity")
            for key in ("cpu_model", "kernel", "os", "cgroup_path", "cgroup_mount"):
                _require(_text(host[key]), f"invalid host {key}")
            _require(isinstance(host["raw_evidence"], dict) and bool(host["raw_evidence"]), "missing raw host evidence")
            runtime = manifest["runtime_manifest"]
            # This validates the frozen expected builds. Actual child evidence is
            # compared again after every launch, never inferred from these paths.
            environments = manifest["environment_manifests"]
            expected_children = [dict(runtime, executable=e.get("executable")) for e in environments]
            validate_preflight(runtime, expected_children, environments, manifest["thread_environment"])
            python_tools = [t for t in manifest["tools"] if t not in ("oceanmesh", "OceanMesh2D")]
            _require(len(environments) == len(python_tools), "missing tool environments")
    except (KeyError, TypeError, ValueError) as exc:
        raise CampaignError(f"invalid manifest: {exc}") from exc
    return manifest


def validate_continuation(manifest, host_manifest, runtime_manifest):
    """Reject host/allocation or runtime drift under a frozen campaign identity."""
    validate_manifest(manifest, require_ready=True)
    _require(manifest["host_manifest"] == host_manifest, "host/allocation drift; create a new campaign")
    _require(manifest["runtime_manifest"] == runtime_manifest, "runtime drift; create a new campaign")
    return True


def runtime_identity():
    path = os.path.realpath(sys.executable)
    with open(path, "rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return dict(version=sys.version, implementation=sys.implementation.name,
                implementation_fields={"name": sys.implementation.name, "cache_tag": sys.implementation.cache_tag,
                                       "version": list(sys.implementation.version), "hexversion": sys.hexversion},
                thread_environment={key: os.environ.get(key) for key in THREAD_ENV}, platform=sysconfig.get_platform(), architecture=platform.machine(), executable=sys.executable, realpath=path, sha256=digest)


def distribution_version(name):
    return importlib.metadata.version(name)


def distribution_versions(names):
    """Collect every declared direct/transitive distribution without importing it."""
    _require(isinstance(names, (list, tuple)) and bool(names) and all(_text(n) for n in names) and len(set(names)) == len(names), "invalid distribution inventory")
    try:
        versions = {name: distribution_version(name) for name in names}
    except importlib.metadata.PackageNotFoundError as exc:
        raise CampaignError(f"missing installed distribution: {exc}") from exc
    _require(all(_pin(v) for v in versions.values()), "missing exact metadata versions")
    return versions


def validate_preflight(parent, children, environments, thread_environment):
    keys = ("version", "implementation", "implementation_fields", "platform", "architecture", "realpath", "sha256")
    _require(isinstance(parent, dict) and str(parent.get("version", "")).startswith("3.11.") and all(parent.get(k) for k in keys), "missing Python 3.11 runtime evidence")
    _require(parent["implementation"] == "cpython" and isinstance(parent["implementation_fields"], dict) and all(parent["implementation_fields"].get(k) is not None for k in ("name", "cache_tag", "version", "hexversion")), "missing implementation build identity")
    implementation = parent["implementation_fields"]
    _require(implementation["name"] == parent["implementation"] and _text(implementation["cache_tag"]) and isinstance(implementation["version"], (list, tuple)) and list(implementation["version"][:2]) == [3, 11], "inconsistent Python build identity")
    _require(_text(parent.get("executable")) and os.path.isabs(parent["executable"]), "missing executable path")
    _require(os.path.isabs(parent["realpath"]) and os.path.normpath(parent["realpath"]) == parent["realpath"], "invalid real interpreter path")
    _require(_hash(parent["sha256"]), "invalid runtime digest")
    _require(bool(children) and bool(environments) and len(children) == len(environments), "missing child/environment evidence")
    _require(thread_environment == THREAD_ENV and parent.get("thread_environment") == THREAD_ENV, "parent thread policy mismatch")
    for child in children:
        _require(isinstance(child, dict) and all(child.get(k) == parent[k] for k in keys), "runtime mismatch")
        _require(_text(child.get("executable")) and os.path.isabs(child["executable"]), "missing child executable")
        _require(child.get("thread_environment") == THREAD_ENV, "child thread policy mismatch")
    identities = []
    for child, env in zip(children, environments):
        _require(isinstance(env, dict), "invalid environment")
        for key in ("tool_distribution", "environment_id", "lock_path", "executable", "source_revision", "version"):
            _require(_text(env.get(key)), f"missing environment {key}")
        _require(env["executable"] == child["executable"], "environment executable does not bind child")
        identities.append(env["environment_id"])
        _require(_hash(env.get("lock_sha256")), "invalid lock digest")
        distributions = env.get("distributions")
        _require(isinstance(distributions, dict) and bool(distributions) and all(_text(k) and _pin(v) for k, v in distributions.items()), "missing distribution pins")
        direct, transitive = env.get("direct_dependencies"), env.get("transitive_dependencies")
        _require(isinstance(direct, list) and bool(direct) and isinstance(transitive, list), "missing dependency inventories")
        _require(all(_text(n) for n in direct + transitive) and len(set(direct + transitive)) == len(direct + transitive), "duplicate dependency inventory")
        _require(set(direct + transitive) == set(distributions), "unlocked direct/transitive dependency")
        _require(env["tool_distribution"] in direct and distributions.get(env["tool_distribution"]) == env["version"], "tool metadata pin mismatch")
        native = env.get("native")
        _require(isinstance(native, dict) and bool(native) and all(_text(k) and _pin(v) for k, v in native.items()), "missing native pins")
        _require(env.get("compatible") is True, "install/import compatibility not established")
        converter = env.get("converter")
        _require(isinstance(converter, dict) and type(converter.get("required")) is bool, "missing converter inventory")
        if converter["required"]:
            _require(_text(converter.get("revision")) and _hash(converter.get("sha256")) and _hash(converter.get("lock_sha256")) and isinstance(converter.get("options"), dict), "incomplete converter identity/options")
            _require(isinstance(converter.get("native"), dict) and bool(converter["native"]) and all(_pin(v) for v in converter["native"].values()), "missing converter native pins")
        else:
            _require(_text(converter.get("reason")), "converter absence needs reason")
    _require(len(set(identities)) == len(identities), "duplicate environment identity")
    return True


def validate_replacement(previous_series_id, rationale):
    if not _text(previous_series_id) or not _text(rationale):
        raise CampaignError("replacement needs superseded series and rationale")
    return True


def validate_conversion(conversion, reference=None):
    """Missing lineage never certifies a conversion, even if dimensions match."""
    required = ('source_reference_id', 'source_sha256', 'converted_sha256', 'converter_revision',
                'converter_sha256', 'netcdf4_version', 'native_versions', 'environment_lock_sha256',
                'argv', 'options', 'contract', 'source_to_row_maps', 'source_to_converted_id_map',
                'source_content_verification', 'conversion_s')
    reasons = [f'missing {key}' for key in required if conversion.get(key) is None]
    for key in ('source_sha256', 'converted_sha256', 'converter_sha256', 'environment_lock_sha256'):
        if not _hash(conversion.get(key)): reasons.append(f'invalid {key}')
    if conversion.get('converter_revision') != 'fair-netcdf-v1': reasons.append('unknown converter revision')
    if not _pin(conversion.get('netcdf4_version')): reasons.append('missing exact netCDF4 version')
    native = conversion.get('native_versions') or {}
    if not all(_pin(native.get(k)) for k in ('netcdf', 'hdf5')): reasons.append('missing native converter pins')
    if not isinstance(conversion.get('argv'), list) or not conversion.get('argv'): reasons.append('missing converter invocation')
    if not isinstance(conversion.get('options'), dict) or not conversion.get('options'): reasons.append('missing converter options')
    if not conversion.get('source_to_row_maps'): reasons.append('missing source-to-row map evidence')
    proof = conversion.get('source_content_verification') or {}
    if proof.get('verification_state') != 'verified': reasons.append('missing independently reopened content proof')
    if reference is not None:
        bindings = {'source_reference_id': reference['reference_id'], 'source_sha256': reference['raw_file_sha256']}
        for key, value in bindings.items():
            if conversion.get(key) != value: reasons.append(f'{key} mismatch')
        expected = dict(n_nodes=reference['header_node_count'], n_elems=reference['header_element_count'],
                        node_sha256=reference['canonical_node_sha256'], element_sha256=reference['canonical_element_sha256'])
        if any(proof.get(k) != v for k, v in expected.items()): reasons.append('converted content mismatch')
        if reference['verification_state'] != 'verified': reasons.append('incomplete source reference')
    if not conversion.get('source_to_converted_id_map', {}).get('independently_reopened'): reasons.append('unproven ID map')
    contract = conversion.get('contract') or {}
    reader = conversion.get('reader')
    base = {'node_id': 'i8', 'x': 'f8', 'y': 'f8', 'depth': 'f8', 'element_id': 'i8', 'arity': 'i8'}
    layouts = {'xarray': ({'node', 'nele', 'max_face_nodes'}, {**base, 'connectivity': 'i8'}, 'original-file-ids'),
               'thalassa': ({'node', 'nele', 'nvertex'}, {**base, 'element': 'i8'}, 'one-based-row-index')}
    layout = layouts.get(reader)
    if layout is None or set(contract.get('dimensions') or {}) != layout[0]:
        reasons.append('missing reader-specific dimension contract')
    if layout is None or contract.get('variables') != layout[1]:
        reasons.append('missing binary64/int64 variable contract')
    if contract.get('fill_value') != 0 or layout is None or contract.get('indexing') != layout[2]: reasons.append('missing fill/ID convention')
    for key in ('conversion_s', 'read_s'):
        value = conversion.get(key)
        if value is not None and (type(value) not in (float, int) or not math.isfinite(value) or value < 0): reasons.append(f'invalid {key}')
    result = copy.deepcopy(conversion)
    result.update(verification_state='incomplete' if reasons else 'verified', reasons=reasons,
                  read_s=conversion.get('read_s'), end_to_end_s=None,
                  end_to_end_definition='conversion_s + operation read_s; conversion outside repeat')
    if not reasons and result['read_s'] is not None:
        result['end_to_end_s'] = result['conversion_s'] + result['read_s']
    return result


def prepare_conversion(reference, reader, output_dir, *, converter_evidence=None, runner=None,
                       python=None, clock=None):
    """Prepare once, before repeat flushes; injected runner supports offline tests."""
    from pathlib import Path
    import subprocess
    import time
    if reference is None or reference.get('verification_state') != 'verified':
        return {'status': 'not-run', 'reason': 'missing verified independent source reference'}
    evidence = converter_evidence or {}
    if reader not in ('xarray', 'thalassa') or not _hash(evidence.get('lock_sha256')) or not evidence.get('compatible'):
        return {'status': 'not-run', 'reason': 'missing compatible pinned converter environment'}
    directory = Path(output_dir); directory.mkdir(parents=True, exist_ok=True)
    ref_path = directory / f"{reference['reference_id']}-{reader}-reference.json"
    with ref_path.open('x') as stream: json.dump(reference, stream, allow_nan=False)
    script = Path(__file__).parent / 'scripts' / 'convert_fair_netcdf.py'
    output = directory / f"{reference['reference_id']}-{reader}.nc"
    argv = [python or sys.executable, str(script), '--reference', str(ref_path), '--reader', reader, '--output', str(output)]
    clock = clock or time.perf_counter
    start = clock()
    try:
        result = (runner or subprocess.run)(argv, capture_output=True, text=True, check=False)
        elapsed = clock() - start
        if result.returncode != 0:
            return {'status': 'failed', 'reason': 'independent converter failed', 'stderr': result.stderr, 'argv': argv}
        lineage = json.loads(result.stdout)
        converted_hash = None
        if output.is_file():
            with output.open('rb') as stream: converted_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
        if converted_hash != lineage.get('converted_sha256') or converted_hash is None:
            return {'status': 'failed', 'reason': 'converted artifact digest mismatch', 'argv': argv}
        lineage.update(conversion_s=elapsed, environment_lock_sha256=evidence['lock_sha256'], argv=argv,
                       options={'reader': reader, 'contract_revision': 'fair-netcdf-v1'})
        with script.open('rb') as stream: expected_script = hashlib.file_digest(stream, 'sha256').hexdigest()
        if lineage.get('converter_sha256') != expected_script:
            return {'status': 'failed', 'reason': 'converter script digest mismatch'}
        checked = validate_conversion(lineage, reference)
        return dict(checked, status='ok' if checked['verification_state'] == 'verified' else 'invalid')
    except (OSError, ValueError, TypeError) as exc:
        return {'status': 'failed', 'reason': f'conversion failed: {exc}', 'argv': argv}


def validate_lock_files(environments):
    """Execution gate: inventory hashes alone cannot certify a physical lock."""
    from pathlib import Path
    for environment in environments:
        try: raw=Path(environment['lock_path']).read_bytes()
        except OSError as exc: raise CampaignError(f'lock unavailable: {exc}') from exc
        _require(hashlib.sha256(raw).hexdigest()==environment['lock_sha256'], 'lock file digest mismatch')
        pins={}
        for line in raw.decode().splitlines():
            line=line.strip()
            if not line or line.startswith('#'): continue
            _require(line.count('==')==1, 'lock must contain exact distribution==version entries')
            name,version=line.split('==')
            _require(_text(name) and _pin(version) and name not in pins, 'invalid or duplicate lock pin')
            pins[name]=version
        _require(bool(pins) and pins==environment['distributions'], 'missing real compatible lock or distribution mismatch')
        _require(environment.get('compatible') is True, 'install/import compatibility not established')
    return True


def execute_schedule(manifest, out_dir, *, execute, host_manifest, runtime_manifest,
                     lock_check=validate_lock_files, lock_path=None, cgroup_check=None, finalize=None):
    """Serial append-only controller seam; execute receives one frozen entry.

    Inputs, references and commands must be prepared before this function. A new
    schedule with explicit replacement lineage is required to replace a series.
    """
    from pathlib import Path
    from .machine import host_lock
    from .records import write_fair_record, write_record_hashes, _write_new
    validate_continuation(manifest,host_manifest,runtime_manifest)
    lock_check(manifest['environment_manifests'])
    from .machine import cgroup_snapshot
    counters=(cgroup_check or cgroup_snapshot)(host_manifest['cgroup_mount'])
    _require(counters.get('memory_events') is not None and 'oom_kill' in counters['memory_events'], 'memory.events unavailable; execution blocked')
    out_dir=Path(out_dir);attempts=out_dir/'attempts';attempts.mkdir(parents=True,exist_ok=True)
    # A common host path is mandatory across output directories by default.
    path=lock_path or '/tmp/chilmesh-fair-host.lock'
    with host_lock(path):
        paths=[attempts/(entry['repeat_id']+'.json') for entry in manifest['schedule']]
        if any(p.exists() for p in paths): raise FileExistsError('immutable attempts already exist; declare a replacement series')
        _write_new((json.dumps(manifest,sort_keys=True,allow_nan=False)+'\n').encode(),out_dir/'evidence'/'campaign.json')
        finish=finalize or (lambda entry,reason,active=False: unstarted_record(manifest,entry,reason,active=active))
        rows=[]
        aborted=False
        for entry,path in zip(manifest['schedule'],paths):
            if aborted:
                row=finish(entry,'controller abort; scheduled attempt never started')
            else:
                try:
                    row=execute(copy.deepcopy(entry))
                    aborted='controller abort' in row.get('reason','')
                except (KeyboardInterrupt, Exception) as exc:
                    aborted=True
                    import sys, traceback
                    traceback.print_exc(file=sys.stderr)  # the record keeps only the exception text
                    row=finish(entry,f'controller abort: {type(exc).__name__}: {exc}',active=True)
            for key in ('campaign_id','manifest_sha256'):
                _require(row[key]==manifest[key],f'attempt {key} differs from manifest')
            for key in ('cell_id','series_id','repeat_id','repeat_index','execution_order','tool','operation','backend'):
                _require(row[key]==entry[key],f'attempt {key} differs from schedule')
            write_fair_record(row,path);rows.append(row)
        write_record_hashes(paths,out_dir/'evidence'/'record-hashes.json')
        return rows


def unstarted_record(manifest, entry, reason, *, active=False):
    """Keep a scheduled cell when preparation/launch could not supply evidence.

    'unavailable' identities explicitly denote absent evidence; they cannot
    acquire a valid measurement or rights ruling. An interrupted callback with
    no retained child evidence is a failed launch, never an invented execution.
    """
    from .records import EVIDENCE_FIELDS, make_fair_record
    from .w2.tools import operation_contract
    try:
        contract = operation_contract(entry['tool'], entry['operation'])
    except ValueError:
        contract = dict(contract_id='unavailable', equivalence='like-for-like', reason=reason)
    unsupported = contract['equivalence'] == 'not-comparable'
    reason = contract['reason'] if unsupported else reason
    envs = manifest.get('environment_manifests') or []
    env = next((e for e in envs if e.get('tool_distribution', '').lower() == entry['tool'].lower()), {})
    row = {k: entry[k] for k in ('cell_id','series_id','repeat_id','repeat_index','execution_order','tool','operation','backend')}
    row.update(campaign_id=manifest['campaign_id'], manifest_sha256=manifest['manifest_sha256'],
        reference_id=entry['mesh']+'-unavailable', contract_id=contract['contract_id'], environment_id=env.get('environment_id','unavailable'),
        tool_version=env.get('version','unavailable'), chilmesh_git_sha=manifest.get('chilmesh_git_sha') or 'unavailable',
        argv=[env.get('executable') or 'unavailable'], options={}, thread_environment=dict(THREAD_ENV),
        host_id=manifest.get('host_id') or 'unavailable', host_allocation={}, created_at=manifest.get('created_at') or '1970-01-01T00:00:00Z',
        outcome='not-comparable' if unsupported else ('failed' if active else 'not-run'), reason=reason,
        equivalence=contract['equivalence'], verification='not-applicable' if unsupported else 'incomplete', publication_eligible=False,
        measurement_validity={k: dict(valid=False, reason=reason) for k in ('time','memory','cache','load','runtime')})
    for section, fields in EVIDENCE_FIELDS.items():
        row[section] = {k: None for k in fields.split()}
        row[section]['missing'] = {k: reason for k in fields.split()}
    if active and not unsupported:
        row['termination']['trigger']='launch-error'
        row['termination']['missing'].pop('trigger')
    return make_fair_record(**row)


def finalize_unstarted_campaign(manifest, out_dir, reason):
    """Append full-matrix reasoned outcomes after failed execution preflight."""
    from pathlib import Path
    from .records import _write_new, write_fair_record, write_record_hashes
    validate_manifest(manifest)
    directory=Path(out_dir)
    if (directory/'evidence'/'campaign.json').exists() or any((directory/'attempts').glob('*.json')):
        raise FileExistsError('campaign already finalized; declare new campaign/series')
    _write_new((json.dumps(manifest,sort_keys=True,allow_nan=False)+'\n').encode(),directory/'evidence'/'campaign.json')
    refs={entry['mesh']+'-unavailable': dict(reference_id=entry['mesh']+'-unavailable', name=entry['mesh'],
        exact_version='unavailable', raw_file_sha256=None) for entry in manifest['ordered_cells']}
    write_references(refs,directory)
    paths=[];rows=[]
    for entry in manifest['schedule']:
        row=unstarted_record(manifest,entry,reason)
        path=directory/'attempts'/(entry['repeat_id']+'.json')
        write_fair_record(row,path);paths.append(path);rows.append(row)
    write_record_hashes(paths,directory/'evidence'/'record-hashes.json')
    return rows


def execution_preflight(manifest, *, host_probe=None, runtime_probe=None, capability_probe=None):
    """Read-only host/cgroup/reset/load checks; never resets caches or launches tools."""
    from pathlib import Path
    from .machine import host_identity, cgroup_snapshot
    reasons=[]
    try:
        validate_manifest(manifest,require_ready=True)
        validate_lock_files(manifest['environment_manifests'])
        validate_continuation(manifest,(host_probe or host_identity)(),(runtime_probe or runtime_identity)())
        if capability_probe:
            capabilities=capability_probe()
        else:
            mount=manifest['host_manifest']['cgroup_mount']
            counters=cgroup_snapshot(mount)
            capabilities=dict(cgroup=counters.get('memory_events') is not None and 'oom_kill' in counters['memory_events'],
                cache=os.access('/proc/sys/vm/drop_caches',os.W_OK),
                load=all(os.access('/proc/'+p,os.R_OK) for p in ('stat','loadavg','diskstats','self/stat','self/io')))
        _require(all(capabilities.get(k) is True for k in ('cgroup','cache','load')), 'missing cgroup/cache/load capabilities')
        # Verify declared binaries bind the frozen base; no alternative interpreter.
        for path in [manifest['runtime_manifest']['executable']] + [e['executable'] for e in manifest['environment_manifests']]:
            real=Path(path).resolve(strict=True)
            _require(str(real)==manifest['runtime_manifest']['realpath'] and os.access(path,os.X_OK), 'base Python realpath mismatch/unavailable')
            with real.open('rb') as stream:
                _require(hashlib.file_digest(stream,'sha256').hexdigest()==manifest['runtime_manifest']['sha256'], 'base Python digest mismatch')
    except (OSError, ValueError, TypeError, KeyError) as exc:
        reasons.append(str(exc))
    return dict(ready=not reasons,reasons=reasons,manifest_sha256=manifest.get('manifest_sha256'))


def write_references(references, out_dir):
    """Freeze exact-version reference metadata with an offline integrity hash."""
    from pathlib import Path
    from .records import _write_new
    raw=(json.dumps(references,sort_keys=True,allow_nan=False)+'\n').encode()
    directory=Path(out_dir)/'evidence'
    digest=_write_new(raw,directory/'references.json')
    _write_new((json.dumps({'references.json':digest},sort_keys=True)+'\n').encode(),directory/'reference-hashes.json')
