#!/usr/bin/env bash
# No installation, download or publication. Python owns frozen schedule checks.
set -euo pipefail
mode=preflight
manifest=
out=
while (($#)); do
    case "$1" in
        --preflight) mode=preflight; shift ;;
        --run) mode=run; shift ;;
        --manifest|--out)
            (($# >= 2)) || { echo "missing argument for $1" >&2; exit 2; }
            if [[ "$1" == --manifest ]]; then manifest=$2; else out=$2; fi
            shift 2 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done
[[ -n "$manifest" && -n "$out" ]] || { echo '--manifest and --out required' >&2; exit 2; }
[[ "$mode" != run || "${CHILMESH_RUN_BENCH:-}" == 1 ]] || { echo 'CHILMESH_RUN_BENCH=1 required with --run' >&2; exit 2; }
base_python=${CHILMESH_FAIR_PYTHON:-}
[[ -n "$base_python" && -x "$base_python" ]] || { echo 'declare executable CHILMESH_FAIR_PYTHON from frozen Python 3.11 base' >&2; exit 2; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
# Host-wide execution lock lives in campaign.py; this lock prevents two scripts
# from writing one output directory, including during preflight/archiving.
mkdir -p "$out"
script_lock="$out/.controller-lock"
mkdir "$script_lock" 2>/dev/null || { echo 'output controller already active' >&2; exit 1; }
trap 'rmdir "$script_lock"' EXIT
preflight_status=0
"$base_python" -m benchmarks.competitors fair-preflight --manifest "$manifest" --out "$out" || preflight_status=$?
if [[ "$mode" == run ]]; then
    status=0
    "$base_python" -m benchmarks.competitors fair-run --manifest "$manifest" --out "$out" || status=$?
    if [[ -d "$out/attempts" && -d "$out/evidence" ]]; then
        [[ ! -e "$out/evidence-archive.tar.gz" ]] || { echo 'archive exists; use a replacement series/output' >&2; exit 1; }
        archive_paths=(evidence attempts)
        [[ ! -d "$out/raw" ]] || archive_paths+=(raw)
        tar -czf "$out/evidence-archive.tar.gz" -C "$out" "${archive_paths[@]}"
    fi
    exit "$status"
fi

exit "$preflight_status"
