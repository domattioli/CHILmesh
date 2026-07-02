"""CHILmesh — 2D mesh library for hydrodynamic domains.

Primary import:
    from chilmesh import Mesh
    mesh = Mesh.read_from_fort14('mesh.14')

Backend selection (C++ / Rust / Python) is automatic; introspect with
``chilmesh.backend_info()`` and force a specific backend via the
``CHILMESH_BACKEND`` environment variable.

The legacy name ``CHILmesh`` is kept as an alias of ``Mesh`` for backward
compatibility with code that imported the class directly.
"""
from __future__ import annotations

from importlib import metadata

from .CHILmesh import CHILmesh, write_fort14
from .gmsh_io import read_msh, write_msh, GmshParseError
from .fort13_io import Fort13, NodalAttribute, read_fort13, write_fort13, Fort13ParseError
from .fort15_io import Fort15, read_fort15, write_fort15, Fort15ParseError
from .summary_io import summary, SummaryError
from .mesh_topology import EdgeMap, quad_from_tri_pair, quads_from_tri_pairs
from .mutations import MutableMesh
from .quality import element_quality
from . import examples
from . import bridge
from . import chilplotting
from . import layer_paths
from .layer_paths import paths_on_outer_vertices
from .admesh_warmstart import optimize_with_admesh_truss, optimize_with_admesh_truss_arrays
from .bridge import (
    MeshAdapterForMADMESHR,
    MeshAdapterForADMESH,
    MeshAdapterForADMESHDomains,
)

Mesh = CHILmesh

try:
    __version__ = metadata.version("chilmesh")
except metadata.PackageNotFoundError:  # pragma: no cover - package not installed
    __version__ = "0.0.0"


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
    "summary",
    "SummaryError",
    # Standalone quality computation
    "element_quality",
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
