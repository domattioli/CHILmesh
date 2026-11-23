#!/usr/bin/env bash
# Install the W1 GPL competitors in their own container image only:
#   oceanmesh (GPLv3, Python) into its own venv, and
#   OceanMesh2D (GPLv3, MATLAB toolbox) as a source checkout run by Octave.
# The benchmark harness calls both as separate processes. Never import either
# into chilmesh or copy their code into this repository.
#
# Usage: docker/oceanmesh_env.sh <out-dir>
# Writes <out-dir>/requirements-oceanmesh-fair.lock, <out-dir>/oceanmesh.import.txt,
# <out-dir>/oceanmesh2d.txt (git revision + Octave version + smoke result).
set -euo pipefail
out=${1:?output dir}
base=${CHILMESH_FAIR_PYTHON:-$(command -v python3.11)}
root=${FAIR_VENV_ROOT:-/opt/fair}
mkdir -p "$out" "$root"
export PIP_DISABLE_PIP_VERSION_CHECK=1

venv="$root/oceanmesh"
rm -rf "$venv"
"$base" -m venv "$venv"
# oceanmesh builds a pybind11 extension against CGAL; the image provides
# libcgal-dev, libgdal-dev and a C++ toolchain.
"$venv/bin/python" -m pip install "oceanmesh @ git+https://github.com/CHLNDDEV/oceanmesh.git@${OCEANMESH_REF:-master}"
status=ok
"$venv/bin/python" -c 'import oceanmesh, importlib.metadata as m; print(m.version("oceanmesh"))' \
    >"$out/oceanmesh.import.txt" 2>&1 || status=failed
echo "oceanmesh import: $status" >>"$out/oceanmesh.import.txt"
"$venv/bin/python" -m pip freeze --all \
    | grep -v -E '^(pip|setuptools|wheel)==' | LC_ALL=C sort >"$out/requirements-oceanmesh-fair.lock"

om2d=${OCEANMESH2D_DIR:-/opt/fair/OceanMesh2D}
rm -rf "$om2d"
git clone --depth 1 ${OCEANMESH2D_REF:+--branch "$OCEANMESH2D_REF"} \
    https://github.com/CHLNDDEV/OceanMesh2D.git "$om2d"
{
    echo "revision: $(git -C "$om2d" rev-parse HEAD)"
    echo "octave: $(octave --version | head -n 1)"
    if octave --no-gui --quiet --eval "addpath(genpath('$om2d')); exist('msh','class') || exist('msh','file')" >/dev/null 2>&1; then
        echo "smoke: ok"
    else
        echo "smoke: failed"
    fi
} >"$out/oceanmesh2d.txt"
cat "$out/oceanmesh.import.txt" "$out/oceanmesh2d.txt"
