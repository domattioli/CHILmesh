"""Valence mesh access for the harness: streaming download, not-served handling, rights.

Meshes come only through the Valence API (``chilmesh.registry``). There is no
local-clone and no direct-source fallback: a mesh that is not served is a
``not-run`` result. ``chilmesh.registry.fetch_mesh_bytes`` holds the whole
payload in RAM, so this module streams the download in chunks instead and
verifies SHA-256 on the fly. ``registry.py`` itself is not changed.
"""
from __future__ import annotations

import hashlib
import json
import os
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Optional

NOT_SERVED_REASON = "mesh not yet served by Valence"
DEFAULT_CHUNK = 1 << 20
_TIMEOUT_S = 60.0


def _default_registry():
    from chilmesh import registry

    return registry


def fetch_mesh_file(name: str, dest_dir, *, registry=None, opener=None,
                    chunk_bytes: int = DEFAULT_CHUNK) -> dict:
    """Stream mesh ``name`` into ``dest_dir/<name>.14``.

    Returns ``{"status": "ok", "path", "sha256", "bytes"}`` or a ``not-run`` /
    ``failed`` dict carrying ``reason``. Never raises for registry errors.
    """
    reg = registry or _default_registry()
    opener = opener or urllib.request.urlopen
    try:
        reg._check_name(name)
        expected = reg.fetch_hash(name)
    except reg.RegistryAuthError:
        return {"status": "not-run", "reason": "Valence token missing or invalid; set VALENCE_TOKEN"}
    except (reg.RegistryNotFoundError, reg.RegistryUnavailableError):
        return {"status": "not-run", "reason": NOT_SERVED_REASON}
    except reg.RegistryError as exc:
        return {"status": "failed", "reason": f"registry error: {exc}"}

    base, token = reg._config()
    req = urllib.request.Request(
        f"{base}/meshes/{urllib.parse.quote(name, safe='')}",
        headers={"Authorization": f"Bearer {token}", "User-Agent": "chilmesh-bench/1"},
    )
    dest = Path(dest_dir) / f"{name}.14"
    dest.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    total = 0
    try:
        with opener(req, timeout=_TIMEOUT_S) as resp, open(dest, "wb") as fh:
            while True:
                chunk = resp.read(chunk_bytes)
                if not chunk:
                    break
                digest.update(chunk)
                fh.write(chunk)
                total += len(chunk)
    except OSError as exc:
        _unlink(dest)
        return {"status": "failed", "reason": f"download failed: {exc}"}
    if digest.hexdigest() != expected:
        _unlink(dest)
        return {"status": "failed", "reason": "payload failed SHA-256 hash check; mesh rejected"}
    return {"status": "ok", "path": str(dest), "sha256": expected, "bytes": total}


def _unlink(path: Path) -> None:
    try:
        os.unlink(path)
    except OSError:
        pass


def load_rights(path) -> dict:
    """Read the manual ``rights.json``. A missing file means no mesh is cleared."""
    path = Path(path)
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def apply_rights(names: list, rights: dict) -> tuple[list, list]:
    """Split mesh names into those cleared to publish and an exclusion list (FR-008)."""
    allowed, excluded = [], []
    for name in names:
        entry: Optional[dict] = rights.get(name)
        if entry and entry.get("publish_results") is True:
            allowed.append(name)
        elif entry is None:
            excluded.append({"mesh": name, "reason": "no rights.json entry; results not publishable"})
        else:
            excluded.append({"mesh": name, "reason": entry.get("note") or "rights do not allow publishing results"})
    return allowed, excluded
