"""``.chil`` container (C1): zip with manifest, exact float64 points, connectivity, boundary segments."""
from __future__ import annotations

import io
import json
import zipfile
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

C1_VERSION = 1


def _npy_bytes(arr: np.ndarray) -> bytes:
    buf = io.BytesIO()
    np.save(buf, arr, allow_pickle=False)
    return buf.getvalue()


def write_chil(mesh, path, metadata: Optional[Dict[str, Any]] = None) -> None:
    """Write ``mesh`` to a ``.chil`` zip container at ``path``."""
    points = np.ascontiguousarray(mesh.points, dtype=np.float64)
    conn = np.asarray(mesh.connectivity_list)
    manifest = {
        "format": "chil",
        "c1_version": C1_VERSION,
        "kind": "mesh",
        "n_points": int(points.shape[0]),
        "n_elements": int(conn.shape[0]),
        "crs": {"native": "unknown"},
        "metadata": dict(metadata) if metadata else {},
        "grid_name": getattr(mesh, "grid_name", None),
    }
    segs = json.dumps(getattr(mesh, "boundary_segments", []) or [], default=_json_default)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("manifest.json", json.dumps(manifest))
        z.writestr("points.npy", _npy_bytes(points))
        z.writestr("connectivity.npy", _npy_bytes(conn))
        z.writestr("boundary_segments.json", segs)


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.generic):
        return o.item()
    raise TypeError(f"not JSON serializable: {type(o)}")


def read_chil(path, **load_kwargs):
    """Read a ``.chil`` container into a CHILmesh (kwargs go to the constructor)."""
    from .CHILmesh import CHILmesh

    with zipfile.ZipFile(path) as z:
        manifest = json.loads(z.read("manifest.json"))
        if manifest.get("format") != "chil":
            raise ValueError("Not a .chil container")
        ver = manifest.get("c1_version")
        if ver != C1_VERSION:
            raise ValueError(f"Unsupported .chil c1_version {ver!r}; supported: {C1_VERSION}")
        points = np.load(io.BytesIO(z.read("points.npy")), allow_pickle=False)
        conn = np.load(io.BytesIO(z.read("connectivity.npy")), allow_pickle=False)
        segs = json.loads(z.read("boundary_segments.json"))
    mesh = CHILmesh(
        connectivity=conn, points=points, grid_name=manifest.get("grid_name"), **load_kwargs
    )
    if segs:
        mesh.boundary_segments = segs
    mesh.metadata = manifest.get("metadata", {})
    mesh.chil_c1_version = ver
    return mesh
