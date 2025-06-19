from importlib import metadata

from .CHILmesh import CHILmesh, write_fort14

try:
    __version__ = metadata.version("chilmesh")
except metadata.PackageNotFoundError:  # pragma: no cover - package not installed
    __version__ = "0.0.0"

__all__ = ["CHILmesh", "write_fort14"]
