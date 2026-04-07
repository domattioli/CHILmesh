"""Bundled example meshes for the chilmesh quickstart and tests.

Each helper returns a fully-initialised :class:`chilmesh.CHILmesh` instance
loaded from a ``.fort.14`` fixture that ships inside the wheel.

Example::

    >>> import chilmesh
    >>> mesh = chilmesh.examples.annulus()
    >>> mesh.n_elems > 0
    True
"""
from __future__ import annotations

from importlib import resources
from pathlib import Path

from .CHILmesh import CHILmesh

__all__ = ["annulus", "donut", "block_o", "structured", "fixture_path"]


def fixture_path(name: str) -> Path:
    """Return the on-disk path of a fixture shipped in ``chilmesh.data``.

    Uses :mod:`importlib.resources` so it works for both editable installs and
    installed wheels.
    """
    ref = resources.files("chilmesh.data").joinpath(name)
    # ``as_file`` would be ideal for zip-safe access, but every file in the
    # ``chilmesh.data`` package is unpacked on disk in practice. Resolve to a
    # ``Path`` so callers can pass it straight to :func:`open`.
    return Path(str(ref))


def _load(name: str, grid_name: str | None = None) -> CHILmesh:
    mesh = CHILmesh.read_from_fort14(fixture_path(name))
    if grid_name is not None:
        mesh.grid_name = grid_name
    return mesh


def annulus() -> CHILmesh:
    """Annular triangular mesh (~580 elements). Good for layer demos."""
    return _load("annulus_200pts.fort.14", grid_name="annulus")


def donut() -> CHILmesh:
    """Multi-ring donut triangular mesh (~276 elements)."""
    return _load("donut_domain.fort.14", grid_name="donut")


def block_o() -> CHILmesh:
    """Block-O triangular mesh (~5,200 elements). Larger stress test."""
    return _load("Block_O.14", grid_name="block_o")


def structured() -> CHILmesh:
    """Small structured triangular mesh (~660 elements)."""
    return _load("structuredMesh1.14", grid_name="structured")
