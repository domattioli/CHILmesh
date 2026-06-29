#!/usr/bin/env bash
# Build + install the chilmesh_cpp C++ backend into the active Python
# environment, then verify the cpp backend actually activates. Reproduces
# the CI `cpp-equivalence` gate locally. Refs CHILmesh #163 / #229.
#
# Usage (from an activated venv that already has chilmesh installed):
#   bash scripts/build_cpp.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python -m pip install "${HERE}/src/chilmesh_cpp"
python - <<'PY'
import chilmesh
info = chilmesh.backend_info()
print("backend_info:", info)
assert "cpp" in info["available"], "C++ backend not active after install"
print("C++ backend OK ->", info["selected"])
PY
