"""W1 tool entries: the three CHILmesh smoothers (plan D2) and the external-tool launch config.

Every CHILmesh smoother takes ``(points, connectivity)`` and returns new ``(n, 2)``
coordinates. The boundary is pinned by the smoother itself; the harness verifies it
afterwards instead of trusting it.
"""
from __future__ import annotations

import os
import shutil  # noqa: F401  (tests patch tools.shutil.which)
from pathlib import Path

import numpy as np

SMOOTH_ITERS = 5
ANGLE_ITERS = 100

ADAPTER_DIR = Path(__file__).resolve().parents[1] / "adapters"

# External tools run in their own environment: an interpreter path from an env var.
EXTERNAL = {
    "oceanmesh": {"env": "CHILMESH_BENCH_PY_OCEANMESH", "script": ADAPTER_DIR / "oceanmesh_clean.py",
                  "license": "GPL (per plan; confirm in research.md)", "run_mode": "subprocess"},
    "oceanmesh2d": {"env": None, "script": ADAPTER_DIR / "oceanmesh2d.m",
                    "license": "unverified (research.md T002)", "run_mode": "subprocess (GNU Octave)"},
}


def _chil(points, connectivity):
    from chilmesh import CHILmesh

    return CHILmesh(connectivity=np.array(connectivity, copy=True), points=np.array(points, dtype=float, copy=True),
                    grid_name="bench")


def _smooth(points, connectivity):
    mesh = _chil(points, connectivity)
    mesh.smooth(n_iter=SMOOTH_ITERS, fixed_boundary=True)
    return np.asarray(mesh.points)[:, :2]


def _angle_based(points, connectivity):
    mesh = _chil(points, connectivity)
    return np.asarray(mesh.angle_based_smoother(n_iter=ANGLE_ITERS))[:, :2]


def _direct(points, connectivity):
    mesh = _chil(points, connectivity)
    return np.asarray(mesh.direct_smoother())[:, :2]


CHILMESH_SMOOTHERS = {
    "smooth": _smooth,
    "angle_based_smoother": _angle_based,
    "direct_smoother": _direct,
}


def external_python(tool: str):
    """Interpreter path for an external tool, or None when its env var is unset."""
    var = EXTERNAL[tool]["env"]
    return os.environ.get(var) if var else None
