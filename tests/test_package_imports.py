"""Smoke test that core import paths and the bundled examples module work.

Catches B6 (missing ``chilmesh/utils/__init__.py``) and verifies that the
``chilmesh.examples`` API and the data fixtures actually ship in the wheel.
"""
from __future__ import annotations


def test_top_level_import():
    import chilmesh

    assert hasattr(chilmesh, "CHILmesh")
    assert hasattr(chilmesh, "write_fort14")
    assert hasattr(chilmesh, "examples")
    assert isinstance(chilmesh.__version__, str)


def test_examples_annulus_loads():
    import chilmesh

    mesh = chilmesh.examples.annulus()
    assert mesh.n_elems > 0
    assert mesh.n_verts > 0


def test_all_example_factories_load():
    import chilmesh

    for name in ("annulus", "donut", "block_o", "structured"):
        mesh = getattr(chilmesh.examples, name)()
        assert mesh.n_elems > 0, f"{name} produced empty mesh"


def test_utils_subpackage_importable():
    """B6: ``chilmesh.utils`` must be a real package and ``plot_utils`` must
    be reachable from it."""
    from chilmesh.utils import plot_utils  # noqa: F401

    assert hasattr(plot_utils, "CHILmeshPlotMixin")
