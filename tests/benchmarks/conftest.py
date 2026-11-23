"""Fixtures for the competitor-benchmark harness tests.

All tests here run without competitor tools, network or Valence. Meshes are
small synthetic triangle, quad and mixed grids (see ``_grids.py``).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from _grids import GRIDS  # noqa: E402  (sibling module; pytest puts this dir on sys.path)


@pytest.fixture(autouse=True)
def isolated_benchmark_environment(tmp_path_factory, monkeypatch):
    """Give subprocess tests source imports and disposable HOME/cache directories."""
    home = tmp_path_factory.mktemp("benchmark-home")
    cache = tmp_path_factory.mktemp("benchmark-cache")
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("XDG_CACHE_HOME", str(cache))
    paths = [str(REPO_ROOT / "src"), str(REPO_ROOT)]
    if os.environ.get("PYTHONPATH"):
        paths.append(os.environ["PYTHONPATH"])
    monkeypatch.setenv("PYTHONPATH", os.pathsep.join(paths))


@pytest.fixture(params=sorted(GRIDS))
def grid(request):
    """(kind, points, connectivity) for each element combination."""
    pts, conn = GRIDS[request.param]()
    return request.param, pts, conn

@pytest.fixture
def fair_fixture(tmp_path, monkeypatch):
    from _grids import FAIR_FORT14
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    mesh = tmp_path / "mesh.14"
    mesh.write_text(FAIR_FORT14)
    cache = tmp_path / "cache"
    cache.mkdir()
    output = tmp_path / "out"
    output.mkdir()
    proc = tmp_path / "proc"
    proc.mkdir()
    (proc / "meminfo").write_text("Cached: 100 kB\n")
    (proc / "42").mkdir()
    (proc / "42" / "status").write_text("Pid: 42\nVmRSS: 100 kB\n")
    (proc / "42" / "cgroup").write_text("0::/synthetic\n")
    cgroup = tmp_path / "cgroup"
    cgroup.mkdir()
    (cgroup / "memory.events").write_text("oom_kill 0\n")
    (cgroup / "memory.max").write_text("8589934592\n")
    return {"mesh": mesh, "proc": proc, "cgroup": cgroup, "output": output,
            "meshes": ["Shinnecock", "WNAT_Hagen", "ENPAC2003", "STOFS"],
            "tools": ["chilmesh", "adcircpy", "oceanmesh", "OceanMesh2D", "thalassa", "xarray"],
            "adapter": {"n_nodes": 4, "n_elems": 2, "node_ids": [10, 20, 40, 80],
                        "element_ids": [7, 9], "connectivity": [[10, 20, 40], [10, 20, 40, 80]],
                        "coordinates": [[0., -0., 1.], [1., 0., 2.], [1., 1., 3.], [0., 1., 4.]]}, "clock": iter(range(100)).__next__}


@pytest.fixture
def fair_manifest(fair_fixture):
    from benchmarks.competitors.campaign import make_manifest
    return make_manifest("synthetic", fair_fixture["meshes"], ["parse"], fair_fixture["tools"])


# The measurement harness uses POSIX process APIs (fork, killpg, wait4, fcntl,
# pass_fds, cgroups). These modules cannot run on Windows.
_POSIX_ONLY_MODULES = {
    "test_campaign.py", "test_fair_adapters.py", "test_fair_campaign.py",
    "test_fair_codespace.py", "test_machine_evidence.py",
    "test_measure_lifecycle.py", "test_w1_fixes.py", "test_w1_tracks.py",
    "test_w2_consumer_metrics.py", "test_w2_wiring.py",
}


def pytest_collection_modifyitems(config, items):
    if sys.platform != "win32":
        return
    skip = pytest.mark.skip(reason="benchmark harness needs POSIX process APIs")
    for item in items:
        if item.path.name in _POSIX_ONLY_MODULES:
            item.add_marker(skip)
