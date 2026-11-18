"""Prepare a frozen fair campaign (campaign.json + evidence/prepared.json) on the Linux host.

Run with the base interpreter the campaign will use, thread env set to 1:
    CHILMESH_FAIR_PYTHON=... python -m benchmarks.competitors.scripts.prepare_fair_campaign --out-dir D

Meshes come only through ``valence.fetch_mesh_file`` (token in env VALENCE_TOKEN,
never written anywhere). An optional local fort.14 is registered internal-only.
Values the host cannot supply stay absent, so preflight reports them.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from benchmarks.competitors import campaign, machine, valence
from benchmarks.competitors.w2 import tools

# (operation, backend): the schedule backend must equal the contract's selected backend for each cell.
OPERATIONS = ("parse", "parse-netcdf", ("adjacency-cpp", "cpp"), ("adjacency-python", "python"))
TOOLS = ("chilmesh", "adcircpy", "thalassa", "xarray")
TOOL_DISTRIBUTION = {"chilmesh": "chilmesh", "adcircpy": "adcircpy", "thalassa": "thalassa", "xarray": "xarray"}
READERS = ("xarray", "thalassa")
DEFAULT_MESHES = ("ShinnecockInlet_coarse.14", "EasternPacific_ENPAC2003.14", "STOFS2DGlobal_v2.1.14")
RIGHTS_MARKER = "no ruling, internal only"
ADAPTERS = Path(__file__).resolve().parents[1] / "adapters"
CONVERTER = Path(__file__).resolve().parent / "convert_fair_netcdf.py"


def _sha(path):
    with open(path, "rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def parse_lock(path):
    pins = {}
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            name, _, version = line.partition("==")
            pins[name] = version
    return pins


def native_pins(python, runner=subprocess.run):
    """Ask an environment's interpreter for libc/netCDF/HDF5 pins; absent when unknown."""
    code = ("import json,platform;d={'glibc':'.'.join(platform.libc_ver()[1].split('.')[:2]) or None}\n"
            "try:\n import netCDF4;d['netcdf']=netCDF4.__netcdf4libversion__;d['hdf5']=netCDF4.__hdf5libversion__\n"
            "except Exception: pass\nprint(json.dumps({k:v for k,v in d.items() if v}))")
    try:
        out = runner([python, "-I", "-c", code], capture_output=True, text=True, check=False)
        return json.loads(out.stdout) if out.returncode == 0 else {}
    except (OSError, ValueError):
        return {}


def environment_manifest(tool, locks_dir, environments, *, runner=subprocess.run, root="/opt/fair"):
    entry = next(e for e in environments["environments"] if e["tool"] == tool) if "environments" in environments \
        else environments[tool]
    lock = Path(entry["lock_path"]) if Path(entry["lock_path"]).is_absolute() else Path(locks_dir) / Path(entry["lock_path"]).name
    pins = parse_lock(lock)
    dist = TOOL_DISTRIBUTION[tool]
    direct = [dist] + (["netCDF4"] if tool in READERS and "netCDF4" in pins else [])
    executable = f"{root}/{tool}/bin/python"
    native = native_pins(executable, runner)
    env = dict(tool_distribution=dist, environment_id=f"{tool}-{entry['lock_sha256'][:12]}", lock_path=str(lock),
               lock_sha256=_sha(lock), executable=executable, source_revision=f"{dist}=={pins.get(dist)}",
               version=pins.get(dist), distributions=pins, direct_dependencies=direct,
               transitive_dependencies=sorted(set(pins) - set(direct)), native=native,
               compatible=bool(entry.get("compatible")))
    if tool in READERS:
        env["converter"] = dict(required=True, revision="fair-netcdf-v1", sha256=_sha(CONVERTER), lock_sha256=env["lock_sha256"],
                                options={"reader": tool, "contract_revision": "fair-netcdf-v1"}, native=native)
    else:
        env["converter"] = dict(required=False, reason="reads fort.14 natively; no conversion")
    return env


def command_for(tool, operation, python, mesh_path, output_path):
    if tool == "chilmesh":
        return [python, "-m", "benchmarks.competitors.w2.measures", operation, "--mesh", mesh_path, "--fair", "--output", output_path]
    script = ADAPTERS / {"adcircpy": "adcircpy_read.py", "thalassa": "thalassa_read.py", "xarray": "xarray_read.py"}[tool]
    return [python, str(script), "--mesh", mesh_path, "--fair"]


def build(args, *, host=None, runtime=None, fetch=None, environments=None, runner=subprocess.run, now=None, environ=None):
    environ = os.environ if environ is None else environ
    out = Path(args.out_dir)
    host = host or machine.host_identity()
    runtime = runtime or campaign.runtime_identity()
    environments = environments if environments is not None else json.loads((Path(args.locks_dir) / "environments.json").read_text())
    env_manifests = [environment_manifest(t, args.locks_dir, environments, runner=runner) for t in TOOLS]
    by_tool = {e["tool_distribution"]: e for e in env_manifests}
    meta = json.loads(Path(args.valence_meta).read_text()) if args.valence_meta else {}
    fetch = fetch or valence.fetch_mesh_file
    work = out / "meshes"
    references, paths, notes = {}, {}, {}
    for name in args.mesh:
        fetched = fetch(name, work)
        info = meta.get(name, {})
        # The registry API exposes no version field: default to the served SHA-256 as the
        # exact version and the API path as the asset, unless --valence-meta supplies both.
        version = info.get("version") or fetched.get("sha256")
        asset = info.get("asset") or f"meshes/{name}"
        ref = valence.reference_from_fetch(name, version, asset, fetched, out / "references")
        if ref.get("status") != "ok":
            notes[name] = ref.get("reason", "not run")
            continue
        references[name], paths[name] = ref, fetched["path"]
    if args.extra_mesh:
        name, _, path = args.extra_mesh.partition("=")
        if _sha(path) != args.extra_sha256:
            raise SystemExit("extra mesh SHA-256 mismatch")
        fetched = {"status": "ok", "path": path, "sha256": args.extra_sha256, "bytes": Path(path).stat().st_size}
        ref = valence.reference_from_fetch(name, "local", "local-file", fetched, out / "references")
        if ref.get("status") == "ok":
            references[name], paths[name] = ref, path
            notes[name] = f"{RIGHTS_MARKER}; {args.extra_note or ''}".strip("; ")
    meshes = [m for m in list(args.mesh) + ([args.extra_mesh.partition("=")[0]] if args.extra_mesh else []) if m in references]
    if not meshes:
        raise SystemExit("no mesh acquired: " + json.dumps(notes, sort_keys=True))
    created = (now or datetime.now(timezone.utc)).strftime("%Y-%m-%dT%H:%M:%SZ")
    manifest = campaign.make_manifest(
        args.campaign_id, meshes, list(OPERATIONS), list(TOOLS), created_at=created,
        chilmesh_git_sha=machine.chilmesh_git_sha(), host_id=host["host_id"],
        host_manifest=dict(host), runtime_manifest=runtime,  # full record: preflight compares it whole
        environment_manifests=env_manifests)
    conversions, converted = {}, {}
    for mesh in meshes:
        for reader in READERS:
            conv = campaign.prepare_conversion(references[mesh], reader, out / "converted" / mesh,
                                               converter_evidence=by_tool[reader], python=by_tool[reader]["executable"])
            conversions[f"{mesh}:{reader}"] = conv
            if conv.get("status") == "ok":
                converted[(mesh, reader)] = conv
    identities, commands, repeat_conversions = {}, {}, {}
    allocation = {k: host[k] for k in ("affinity", "cpu_quota", "cpuset", "memory_limits")}
    for entry in manifest["schedule"]:
        env = by_tool[entry["tool"]]
        contract = tools.operation_contract(entry["tool"], entry["operation"])
        ident = {k: entry[k] for k in ("cell_id", "series_id", "repeat_id", "repeat_index", "execution_order", "tool", "operation", "backend")}
        ident.update(campaign_id=manifest["campaign_id"], manifest_sha256=manifest["manifest_sha256"], mesh=entry["mesh"],
                     environment_id=env["environment_id"], tool_version=env["version"], chilmesh_git_sha=manifest["chilmesh_git_sha"],
                     thread_environment=dict(campaign.THREAD_ENV), host_id=host["host_id"], host_allocation=allocation,
                     contract_id=contract["contract_id"], options={}, created_at=created)
        identities[entry["repeat_id"]] = ident
        if contract["equivalence"] == "not-comparable" or entry["mesh"] not in paths:
            continue
        op = entry["operation"]
        if op == "parse-netcdf":
            conv = converted.get((entry["mesh"], entry["tool"]))
            if conv is None:
                continue
            repeat_conversions[entry["repeat_id"]] = conv
            target = conv["converted_path"]
        else:
            target = paths[entry["mesh"]]
        commands[entry["repeat_id"]] = command_for(entry["tool"], op, env["executable"], target,
                                                   str(out / "raw" / entry["repeat_id"] / "output.json"))
    return manifest, dict(identities=identities, references=references, commands=commands,
                          conversions=repeat_conversions, notes=notes)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--campaign-id", default="fair-" + datetime.now(timezone.utc).strftime("%Y%m%d"))
    ap.add_argument("--locks-dir", default="/opt/fair/locks")
    ap.add_argument("--mesh", action="append", default=None)
    ap.add_argument("--valence-meta", help="JSON {name: {version, asset}} for exact Valence versions")
    ap.add_argument("--extra-mesh", help="NAME=PATH local fort.14, internal-only")
    ap.add_argument("--extra-sha256")
    ap.add_argument("--extra-note", default="")
    args = ap.parse_args(argv)
    args.mesh = args.mesh or list(DEFAULT_MESHES)
    if args.extra_mesh and not (args.extra_sha256 and re.fullmatch("[0-9a-f]{64}", args.extra_sha256)):
        ap.error("--extra-mesh needs --extra-sha256 (64 hex)")
    manifest, prepared = build(args)
    out = Path(args.out_dir)
    (out / "evidence").mkdir(parents=True, exist_ok=True)
    (out / "campaign.json").write_text(json.dumps(manifest, sort_keys=True, allow_nan=False))
    (out / "evidence" / "prepared.json").write_text(json.dumps(prepared, sort_keys=True, allow_nan=False))
    missing = sorted(manifest["missing"])
    print(f"manifest: {out / 'campaign.json'}  missing evidence: {missing or 'none'}")
    print(f"next: CHILMESH_FAIR_PYTHON={sys.executable} bash benchmarks/competitors/scripts/run_fair_codespace.sh "
          f"--preflight --manifest {out / 'campaign.json'} --out {out / 'run'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
