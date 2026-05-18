from __future__ import annotations

from importlib import metadata

from .CHILmesh import CHILmesh, write_fort14
from .mesh_topology import EdgeMap
from .mutations import MutableMesh
from . import examples
from . import bridge
from . import layer_paths
from .layer_paths import paths_on_outer_vertices
from .admesh_warmstart import optimize_with_admesh_truss, optimize_with_admesh_truss_arrays

try:
    __version__ = metadata.version("chilmesh")
except metadata.PackageNotFoundError:  # pragma: no cover - package not installed
    __version__ = "0.0.0"

__all__ = [
    "CHILmesh",
    "write_fort14",
    "EdgeMap",
    "MutableMesh",
    "examples",
    "bridge",
    "layer_paths",
    "paths_on_outer_vertices",
    "optimize_with_admesh_truss",
    "optimize_with_admesh_truss_arrays",
]
