"""Fixture names and loader shared by test modules.

Test modules import these names from here, not from ``conftest``. Two
``conftest.py`` files exist (``tests/`` and ``tests/benchmarks/``), and both
register as the top-level module ``conftest``. A test module that runs
``from conftest import ...`` can get the wrong one, depending on collection
order.
"""
from __future__ import annotations

from chilmesh import examples as _examples

FIXTURE_NAMES = ["annulus", "donut", "block_o", "structured", "quad_2x2"]
TRI_FIXTURE_NAMES = [n for n in FIXTURE_NAMES if n != "quad_2x2"]


def _load(name: str):
    """Load a fixture mesh (memoized when ``tests/conftest.py`` has run)."""
    return getattr(_examples, name)()
