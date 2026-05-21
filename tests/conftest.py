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


# Memoize the example loaders for the duration of the test session so the
# O(n^2) adjacency build (deferred to the 0.1.2 perf release) doesn't
# dominate the test runtime. Block_O alone takes ~30s on first load and
# is touched by many parametrized tests. Tests that mutate a mesh in place
# must call ``.copy()`` first.
#
# xdist safety (#122): pytest-xdist uses process-per-worker, so each worker
# owns an independent ``_MESH_CACHE`` dict — no cross-worker race on
# initialization. Cost: each worker rebuilds the cache once. Acceptable
# given the alternative (shared-memory cache) needs file locks.
_MESH_CACHE: dict = {}


def _make_cached(name):
    original = getattr(_examples, name)

    def cached(**kwargs):
        if kwargs:
            # If kwargs provided, create a fresh mesh with those kwargs
            return original(**kwargs)
        # Otherwise use cached version
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


@pytest.fixture(params=[None, "halfedge"])
def topology_backend_env(request, monkeypatch):
    """Monkeypatch CHILMESH_TOPOLOGY_BACKEND for backend switching tests.

    Parametrized over [None, 'halfedge']. None = unset env var (uses default
    EdgeMap). 'halfedge' = switch to half-edge backend. Monkeypatch ensures
    cleanup on teardown to prevent env-var leakage between tests.

    Used by test_halfedge_basic.py and test_halfedge_equivalence.py to verify
    both backends produce identical results without modifying test code.
    """
    backend = request.param
    if backend is None:
        # Ensure env var is not set
        monkeypatch.delenv("CHILMESH_TOPOLOGY_BACKEND", raising=False)
    else:
        monkeypatch.setenv("CHILMESH_TOPOLOGY_BACKEND", backend)
    return backend
