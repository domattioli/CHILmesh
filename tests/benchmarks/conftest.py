"""Fixtures for the competitor-benchmark harness tests.

All tests here run without competitor tools, network or Valence. Meshes are
small synthetic triangle, quad and mixed grids (see ``_grids.py``).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from _grids import GRIDS  # noqa: E402  (sibling module; pytest puts this dir on sys.path)


@pytest.fixture(params=sorted(GRIDS))
def grid(request):
    """(kind, points, connectivity) for each element combination."""
    pts, conn = GRIDS[request.param]()
    return request.param, pts, conn
