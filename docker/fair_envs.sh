#!/usr/bin/env bash
# Build the four fair-campaign environments (chilmesh, adcircpy, thalassa,
# xarray) from ONE CPython 3.11 base, import-test each tool, and freeze exact
# locks. Runs inside the container or Codespace only, never on the Mac.
#
# Usage: docker/fair_envs.sh <repo-root> <out-dir>
# Writes <out-dir>/requirements-<tool>-fair.lock, <out-dir>/<tool>.import.txt
# and <out-dir>/environments.json (base interpreter identity, lock SHA-256,
# importlib.metadata versions, compatibility result per tool).
# Copy the locks into benchmarks/competitors/envs/ only after review.
set -euo pipefail
repo=${1:?repo root}
out=${2:?output dir}
base=${CHILMESH_FAIR_PYTHON:-$(command -v python3.11)}
"$base" -c 'import sys; assert sys.version_info[:2] == (3, 11), sys.version' \
    || { echo "base interpreter must be CPython 3.11: $base" >&2; exit 2; }
root=${FAIR_VENV_ROOT:-/opt/fair}
mkdir -p "$out" "$root"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PIP_DISABLE_PIP_VERSION_CHECK=1

declare -A IMPORTS=(
    [chilmesh]="chilmesh chilmesh_cpp"
    [adcircpy]="adcircpy"
    [thalassa]="thalassa"
    [xarray]="xarray netCDF4"
)

for tool in chilmesh adcircpy thalassa xarray; do
    venv="$root/$tool"
    rm -rf "$venv"
    "$base" -m venv "$venv"
    py="$venv/bin/python"
    if [[ "$tool" == chilmesh ]]; then
        "$py" -m pip install "$repo" "$repo/src/chilmesh_cpp"
    else
        "$py" -m pip install -r "$repo/benchmarks/competitors/envs/requirements-$tool.txt"
    fi
    status=ok
    for mod in ${IMPORTS[$tool]}; do
        "$py" -c "import $mod" >>"$out/$tool.import.txt" 2>&1 || status=failed
    done
    echo "$tool import: $status" >>"$out/$tool.import.txt"
    # "pip list" gives name==version for every distribution, including chilmesh installed from the
    # working tree; the campaign gate needs the lock to equal the installed set exactly.
    "$py" -m pip list --format=freeze \
        | grep -v -E '^(pip|setuptools|wheel)==' \
        | LC_ALL=C sort >"$out/requirements-$tool-fair.lock"
    echo "$tool $status" >>"$out/status.txt"
done

"$base" - "$out" "$root" "$repo" <<'PY'
import hashlib, json, platform, subprocess, sys
from pathlib import Path
out, root, repo = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
base = Path(sys.executable).resolve()
status = dict(line.split() for line in (out / "status.txt").read_text().splitlines())
envs = []
for tool, result in status.items():
    lock = out / f"requirements-{tool}-fair.lock"
    raw = lock.read_bytes()
    pins = dict(l.split("==", 1) for l in raw.decode().splitlines() if "==" in l)
    envs.append({"tool": tool, "lock_path": str(lock), "lock_sha256": hashlib.sha256(raw).hexdigest(),
                 "distributions": pins, "compatible": result == "ok",
                 "venv": str(root / tool)})
rev = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
dirty = bool(subprocess.run(["git", "-C", str(repo), "status", "--porcelain"], capture_output=True, text=True).stdout.strip())
doc = {"base_python": str(base), "base_sha256": hashlib.sha256(base.read_bytes()).hexdigest(),
       "base_version": sys.version, "machine": platform.machine(), "platform": platform.platform(),
       "chilmesh_revision": rev, "chilmesh_tree_dirty": dirty, "environments": envs}
(out / "environments.json").write_text(json.dumps(doc, indent=2))
print(json.dumps({e["tool"]: e["compatible"] for e in envs}))
PY
