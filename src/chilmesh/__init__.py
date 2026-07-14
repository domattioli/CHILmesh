"""CHILmesh — 2D mesh library for hydrodynamic domains.

Primary import:
    from chilmesh import Mesh
    mesh = Mesh.read_from_fort14('mesh.14')

Backend selection (C++ / Rust / Python) is automatic; introspect with
``chilmesh.backend_info()`` and force a specific backend via the
``CHILMESH_BACKEND`` environment variable.

The legacy name ``CHILmesh`` is kept as an alias of ``Mesh`` for backward
compatibility with code that imported the class directly.

Lazy loading (#255): the heavy mesh/plotting/geometry surface (which pulls in
numpy + matplotlib) is loaded on first attribute access via PEP 562
``__getattr__``. Only the stdlib-pure ``fort14_io`` names are imported eagerly,
so lightweight consumers can ``from chilmesh import read_fort14_raw`` without
dragging in numpy or matplotlib.
"""
from __future__ import annotations

import importlib
from importlib import metadata

# Eagerly export the stdlib-pure fort.14 raw I/O surface only. fort14_io imports
# nothing heavier than the standard library, so this path stays numpy/matplotlib
# free (#255 — keeps Valence's pure-stdlib fort.14 delegation contract intact).
from .fort14_io import (
    Fort14Raw,
    OpenBoundary,
    FlowBoundary,
    read_fort14_raw,
    write_fort14_raw,
    Fort14ParseError,
)

try:
    __version__ = metadata.version("chilmesh")
except metadata.PackageNotFoundError:  # pragma: no cover - package not installed
    __version__ = "0.0.0"


# --- Lazy attribute loading (PEP 562) --------------------------------------
# name -> relative module that defines it. Accessing any of these triggers the
# heavy import chain (numpy / matplotlib) only when actually needed.
_LAZY_ATTRS = {
    # .CHILmesh  (the CHILmesh/Mesh class names are handled by the
    # _ChilmeshModule __getattribute__ guard below, not here — see #255.)
    "write_fort14": ".CHILmesh",
    # .gmsh_io
    "read_msh": ".gmsh_io",
    "write_msh": ".gmsh_io",
    "GmshParseError": ".gmsh_io",
    # .fort13_io
    "Fort13": ".fort13_io",
    "NodalAttribute": ".fort13_io",
    "read_fort13": ".fort13_io",
    "write_fort13": ".fort13_io",
    "Fort13ParseError": ".fort13_io",
    # .fort15_io
    "Fort15": ".fort15_io",
    "read_fort15": ".fort15_io",
    "write_fort15": ".fort15_io",
    "Fort15ParseError": ".fort15_io",
    # .summary_io
    "summary": ".summary_io",
    "SummaryError": ".summary_io",
    # .mesh_topology
    "EdgeMap": ".mesh_topology",
    "quad_from_tri_pair": ".mesh_topology",
    "quads_from_tri_pairs": ".mesh_topology",
    # .mutations
    "MutableMesh": ".mutations",
    # .quality
    "element_quality": ".quality",
    "courant_number": ".quality",
    "cfl_gate": ".quality",
    # .geometry
    "haversine_m": ".geometry",
    "edge_lengths": ".geometry",
    "EARTH_RADIUS_M": ".geometry",
    "convex_hull": ".geometry",
    "is_antimeridian_wrapping": ".geometry",
    "split_antimeridian_bbox": ".geometry",
    "bbox_iou": ".geometry",
    "hausdorff_distance": ".geometry",
    # .node_match
    "NodeMatch": ".node_match",
    "match_nodes": ".node_match",
    "derive_tolerance": ".node_match",
    "nodal_field_delta": ".node_match",
    # .layer_paths
    "paths_on_outer_vertices": ".layer_paths",
    # .admesh_warmstart
    "optimize_with_admesh_truss": ".admesh_warmstart",
    "optimize_with_admesh_truss_arrays": ".admesh_warmstart",
    # .bridge
    "MeshAdapterForMADMESHR": ".bridge",
    "MeshAdapterForADMESH": ".bridge",
    "MeshAdapterForADMESHDomains": ".bridge",
}

# Submodules re-exported lazily as attributes of the package.
_LAZY_SUBMODULES = ("examples", "bridge", "chilplotting", "layer_paths")


def __getattr__(name):  # PEP 562 lazy attribute access
    # ``Mesh`` / ``CHILmesh`` are resolved by the _ChilmeshModule
    # __getattribute__ guard below (submodule name-collision, #255) and never
    # reach here.
    if name in _LAZY_SUBMODULES:
        mod = importlib.import_module(f".{name}", __name__)
        globals()[name] = mod
        return mod
    target = _LAZY_ATTRS.get(name)
    if target is not None:
        mod = importlib.import_module(target, __name__)
        value = getattr(mod, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))


def backend_info() -> dict:
    """Return information about the available and selected mesh backends.

    Returns:
        dict with keys:
            available: list of backend names available in this environment
            selected: the backend currently used by default
            versions: dict of backend name to version string

    Example:
        >>> import chilmesh
        >>> info = chilmesh.backend_info()
        >>> info['selected']
        'cpp'
    """
    from .backends.cpp_backend import CPP_AVAILABLE, _cpp as _cpp_module
    from .backends.rust_backend import RUST_AVAILABLE, _rust as _rust_module

    available = ["python"]
    versions = {"python": __version__}

    if CPP_AVAILABLE:
        available.insert(0, "cpp")
        versions["cpp"] = getattr(_cpp_module, "__version__", "0.6.0.dev0")

    if RUST_AVAILABLE:
        available.insert(-1 if "cpp" in available else 0, "rust")
        versions["rust"] = "0.5.0.dev0"

    import os
    env_override = os.environ.get("CHILMESH_BACKEND", "").lower()
    if env_override in available:
        selected = env_override
    else:
        selected = available[0]

    return {
        "available": available,
        "selected": selected,
        "versions": versions,
    }


__all__ = [
    # Primary class
    "Mesh",
    "CHILmesh",  # legacy alias, kept for backward compat
    # Topology / I/O
    "EdgeMap",
    "quad_from_tri_pair",
    "quads_from_tri_pairs",
    "MutableMesh",
    "write_fort14",
    "read_msh",
    "write_msh",
    "GmshParseError",
    "Fort13",
    "NodalAttribute",
    "read_fort13",
    "write_fort13",
    "Fort13ParseError",
    "Fort15",
    "read_fort15",
    "write_fort15",
    "Fort15ParseError",
    "Fort14Raw",
    "OpenBoundary",
    "FlowBoundary",
    "read_fort14_raw",
    "write_fort14_raw",
    "Fort14ParseError",
    "summary",
    "SummaryError",
    # Standalone quality computation
    "element_quality",
    "courant_number",
    "cfl_gate",  # CFL / Courant gate
    # Geodesic / planar geometry helpers
    "haversine_m",
    "edge_lengths",
    "EARTH_RADIUS_M",
    "convex_hull",
    "is_antimeridian_wrapping",
    "split_antimeridian_bbox",
    "bbox_iou",
    "hausdorff_distance",
    # Cross-mesh node matching + nodal-field deltas (#239)
    "NodeMatch",
    "match_nodes",
    "derive_tolerance",
    "nodal_field_delta",
    # Backend introspection
    "backend_info",
    # Submodules
    "examples",
    "bridge",
    "chilplotting",
    "layer_paths",
    # Convenience re-exports
    "paths_on_outer_vertices",
    "optimize_with_admesh_truss",
    "optimize_with_admesh_truss_arrays",
    # Downstream adapters
    "MeshAdapterForMADMESHR",
    "MeshAdapterForADMESH",
    "MeshAdapterForADMESHDomains",
]


# --- Class/submodule name-collision guard (#255) ---------------------------
# ``CHILmesh`` names both the primary class and the ``chilmesh.CHILmesh``
# submodule that defines it. Importing that submodule (directly, or via any
# ``from .CHILmesh import ...`` / ``from . import CHILmesh``) makes the import
# system bind the *package attribute* ``chilmesh.CHILmesh`` to the module,
# shadowing the class for ``chilmesh.CHILmesh`` and ``from chilmesh import
# CHILmesh``. A PEP 562 ``__getattr__`` can't override this (it only fires when
# the attribute is missing, and the submodule binding makes it present). A
# module ``__getattribute__`` override resolves both public names to the class
# unconditionally, independent of import order — while still importing the
# (numpy-heavy) class module only on first access, so the lightweight fort.14
# path stays stdlib-only.
import sys as _sys
from types import ModuleType as _ModuleType


class _ChilmeshModule(_ModuleType):
    def __getattribute__(self, name):
        if name in ("CHILmesh", "Mesh"):
            cache = super().__getattribute__("__dict__")
            cls = cache.get("_mesh_cls")
            if cls is None:
                from chilmesh.CHILmesh import CHILmesh as cls
                cache["_mesh_cls"] = cls
            return cls
        return super().__getattribute__(name)


_sys.modules[__name__].__class__ = _ChilmeshModule
