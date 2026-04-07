"""Shared pytest fixtures for the chilmesh test suite.

The four ``.fort.14`` mesh fixtures shipped in :mod:`chilmesh.data` are
parametrised here so individual tests can declare ``mesh`` and receive each
fixture in turn. The fixture names match the ``chilmesh.examples`` factory
function names so failures point straight at the helper that loads them.
"""
from __future__ import annotations

import pytest

import chilmesh

FIXTURE_NAMES = ["annulus", "donut", "block_o", "structured"]


# Cache built meshes across tests so the O(n^2) adjacency build (deferred to
# the 0.1.2 perf release) doesn't dominate the test runtime. Tests that mutate
# the mesh in place must call ``.copy()`` first.
_MESH_CACHE: dict = {}


def _load(name: str):
    if name not in _MESH_CACHE:
        _MESH_CACHE[name] = getattr(chilmesh.examples, name)()
    return _MESH_CACHE[name]


@pytest.fixture(params=FIXTURE_NAMES)
def mesh(request):
    """Yield a cached :class:`chilmesh.CHILmesh` for each fixture."""
    return _load(request.param)


@pytest.fixture(params=FIXTURE_NAMES)
def fresh_mesh(request):
    """Yield a fresh (uncached) mesh — for tests that need to mutate."""
    return getattr(chilmesh.examples, request.param)()


@pytest.fixture(params=FIXTURE_NAMES)
def fixture_name(request):
    """Yield the bare fixture name for tests that want to load it lazily."""
    return request.param
