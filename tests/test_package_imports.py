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


def test_mesh_alias_is_chilmesh():
    """``Mesh`` alias must be the same class as ``CHILmesh``."""
    from chilmesh import Mesh, CHILmesh

    assert Mesh is CHILmesh
    assert hasattr(Mesh, "read_from_fort14")
    assert hasattr(Mesh, "from_admesh_domain")


def test_version_string():
    import chilmesh

    assert chilmesh.__version__ != "0.0.0", "package not installed properly"
    assert chilmesh.__version__ == "1.0.0", f"expected 1.0.0, got {chilmesh.__version__}"


def test_backend_info_structure():
    from chilmesh import backend_info

    info = backend_info()
    assert "available" in info
    assert "selected" in info
    assert "versions" in info
    assert isinstance(info["available"], list)
    assert len(info["available"]) >= 1
    assert info["selected"] in info["available"]
    assert "python" in info["available"]
    assert "python" in info["versions"]


def test_bridge_adapters_importable():
    from chilmesh import MeshAdapterForMADMESHR, MeshAdapterForADMESH, MeshAdapterForADMESHDomains  # noqa: F401

    assert MeshAdapterForMADMESHR is not None
    assert MeshAdapterForADMESH is not None
    assert MeshAdapterForADMESHDomains is not None


def test_all_exports_importable():
    """Every name in __all__ must be importable from chilmesh."""
    import chilmesh

    for name in chilmesh.__all__:
        assert hasattr(chilmesh, name), f"__all__ lists '{name}' but it's not importable"
