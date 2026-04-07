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

FIXTURE_NAMES = ["annulus", "donut", "block_o", "structured"]


# Memoize the example loaders for the duration of the test session so the
# O(n^2) adjacency build (deferred to the 0.1.2 perf release) doesn't
# dominate the test runtime. Block_O alone takes ~30s on first load and
# is touched by many parametrized tests. Tests that mutate a mesh in place
# must call ``.copy()`` first.
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
