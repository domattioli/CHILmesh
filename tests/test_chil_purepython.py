"""Pure-Python path and dependency tests for .chil (T009b; SC-006, FR-070).

The flow test runs in a subprocess with the C++ backend blocked. It fails
until the implementation phases land; T054 only runs it.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

_FLOW = r"""
import sys
sys.modules['chilmesh_cpp'] = None
import importlib, os, tempfile
from pathlib import Path
chil_io = importlib.import_module("chilmesh.chil_io")
from chilmesh import examples
mesh = examples.annulus()
d = Path(tempfile.mkdtemp())
f14 = d / "fort.14"
from chilmesh.CHILmesh import write_fort14
write_fort14(mesh, str(f14))
packed = d / "packed.chil"
chil_io.pack(packed, f14)
chil_io.validate_chil(packed)
assert chil_io.summary_chil(packed)["format"] == "chil"
loaded = chil_io.load_chil(packed)
assert loaded.n_elems == mesh.n_elems
saved = d / "saved.chil"
mesh.save(saved)
chil_io.save_chil(mesh, d / "saved2.chil")
assert chil_io.read_member(saved, "fort14")
chil_io.rewrite_chil(saved, d / "rw.chil", ext={"valence": {"a": 1}})
assert chil_io.validate_chil(d / "rw.chil").ext_namespaces == ("valence",)
import chilmesh
assert chilmesh.CHILmesh.load(packed).n_elems == mesh.n_elems
print("OK")
"""


def _env():
    env = dict(os.environ)
    env["CHILMESH_BACKEND"] = "python"
    env["PYTHONPATH"] = os.pathsep.join(
        [str(REPO / "src"), str(REPO / "tests"), env.get("PYTHONPATH", "")])
    env.setdefault("MPLCONFIGDIR", env.get("TMPDIR", "/tmp"))
    return env


def test_full_flow_without_cpp_backend():
    r = subprocess.run([sys.executable, "-c", _FLOW], env=_env(),
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr[-2000:]
    assert r.stdout.strip().endswith("OK")


def test_import_chil_io_does_not_import_numpy():
    code = ("import sys; import chilmesh.chil_io; "
            "sys.exit(1 if 'numpy' in sys.modules else 0)")
    r = subprocess.run([sys.executable, "-c", code], env=_env(),
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr[-2000:]


def test_runtime_dependencies_are_exactly_numpy_scipy_matplotlib():
    import re
    text = (REPO / "pyproject.toml").read_text(encoding="utf-8")
    block = re.search(r"^dependencies\s*=\s*\[(.*?)\]", text, re.S | re.M).group(1)
    names = sorted(re.match(r"[A-Za-z0-9_.-]+", x.strip().strip('",')).group(0).lower()
                   for x in block.splitlines() if x.strip().strip(",") and not x.strip().startswith("#"))
    assert names == ["matplotlib", "numpy", "scipy"]
