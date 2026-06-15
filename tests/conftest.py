"""Shared pytest fixtures for the chilmesh test suite.

The four ``.fort.14`` mesh fixtures shipped in :mod:`chilmesh.data` are
parametrised here so individual tests can declare ``mesh`` and receive each
fixture in turn. The fixture names match the ``chilmesh.examples`` factory
function names so failures point straight at the helper that loads them.
"""
from __future__ import annotations

import pytest

import chilmesh
from chilmesh import examples as _examples

FIXTURE_NAMES = ["annulus", "donut", "block_o", "structured", "quad_2x2"]
TRI_FIXTURE_NAMES = [n for n in FIXTURE_NAMES if n != "quad_2x2"]


# Memoize the example loaders for the duration of the test session so repeated
# fixture loads don't dominate the test runtime. The pure-Python full init is
# linear and fast on these fixtures (Block_O ~5k elems is ~0.2-0.3s; see #202
# for the perf-regression guard), but several parametrized tests load each
# fixture many times, so caching still saves wall-clock. Tests that mutate a
# mesh in place must call ``.copy()`` first.
#
# xdist safety (#122): pytest-xdist uses process-per-worker, so each worker
# owns an independent ``_MESH_CACHE`` dict — no cross-worker race on
# initialization. Cost: each worker rebuilds the cache once. Acceptable
# given the alternative (shared-memory cache) needs file locks.
_MESH_CACHE: dict = {}


def _make_cached(name):
    original = getattr(_examples, name)

    def cached():
        if name not in _MESH_CACHE:
            _MESH_CACHE[name] = original()
        return _MESH_CACHE[name]

    cached.__wrapped__ = original
    return cached


for _n in FIXTURE_NAMES:
    if hasattr(_examples, _n):
        setattr(_examples, _n, _make_cached(_n))


def _load(name: str):
    return getattr(_examples, name)()


@pytest.fixture(params=FIXTURE_NAMES)
def mesh(request):
    """Yield a cached :class:`chilmesh.CHILmesh` for each fixture."""
    return _load(request.param)


@pytest.fixture(params=FIXTURE_NAMES)
def fresh_mesh(request):
    """Yield a fresh (uncached) mesh — for tests that need to mutate."""
    return getattr(_examples, request.param).__wrapped__()


@pytest.fixture(params=FIXTURE_NAMES)
def fixture_name(request):
    """Yield the bare fixture name for tests that want to load it lazily."""
    return request.param
