"""Thin client for a token-gated mesh registry API.

Uses only the standard library. The token comes from ``VALENCE_TOKEN`` and
the base URL from ``VALENCE_API_URL``. Every download is checked against the
SHA-256 hash the API reports. The token never appears in messages or reprs.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import socket
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:  # pragma: no cover
    from .CHILmesh import CHILmesh

__all__ = [
    "RegistryError",
    "RegistryAuthError",
    "RegistryNotFoundError",
    "RegistryUnavailableError",
    "RegistryIntegrityError",
    "list_meshes",
    "fetch_hash",
    "fetch_mesh_bytes",
    "fetch_mesh",
    "download_mesh",
]

# Placeholder: the operator sets the production URL after deployment. While
# empty, the VALENCE_API_URL environment variable is required.
_DEFAULT_BASE_URL = ""
_TIMEOUT_S = 30.0
_TOKEN_VAR = "VALENCE_TOKEN"
_URL_VAR = "VALENCE_API_URL"
_NAME_RE = re.compile(r"^[A-Za-z0-9._-]+$")
_HEX64_RE = re.compile(r"^[0-9a-fA-F]{64}$")


class RegistryError(RuntimeError):
    """Base class for registry failures."""


class RegistryAuthError(RegistryError):
    """Token missing, invalid or expired."""


class RegistryNotFoundError(RegistryError, ValueError):
    """The mesh name is not in the registry or is not a valid name."""


class RegistryUnavailableError(RegistryError):
    """Registry API unreachable, rate limited or failing."""


class RegistryIntegrityError(RegistryError):
    """The payload failed its hash check or has an unexpected shape."""


def _config() -> tuple[str, str]:
    token = os.environ.get(_TOKEN_VAR, "").strip()
    if not token:
        raise RegistryAuthError(
            f"Registry access needs a token. Set the {_TOKEN_VAR} "
            "environment variable."
        )
    base = os.environ.get(_URL_VAR, "").strip() or _DEFAULT_BASE_URL
    if not base:
        raise RegistryUnavailableError(
            "Registry API unavailable: no API URL configured. "
            f"Set the {_URL_VAR} environment variable."
        )
    return base.rstrip("/"), token


def _check_name(name: str) -> str:
    if not isinstance(name, str) or not _NAME_RE.match(name) or name in (".", ".."):
        raise RegistryNotFoundError("Invalid mesh name.")
    return name


def _get(path: str) -> bytes:
    base, token = _config()
    req = urllib.request.Request(
        base + path, headers={"Authorization": f"Bearer {token}"}
    )
    try:
        with urllib.request.urlopen(req, timeout=_TIMEOUT_S) as resp:
            return resp.read()
    except urllib.error.HTTPError as e:
        code = e.code
        if code == 401:
            raise RegistryAuthError(
                f"Registry token is invalid or expired. Check {_TOKEN_VAR}."
            ) from None
        if code == 404:
            raise RegistryNotFoundError("Mesh not found in the registry.") from None
        if code == 429:
            retry = None
            try:
                retry = e.headers.get("Retry-After") if e.headers else None
            except Exception:
                retry = None
            msg = "Registry API unavailable: rate limit reached."
            if retry:
                msg += f" Retry after {retry} seconds."
            else:
                msg += " Retry later."
            raise RegistryUnavailableError(msg) from None
        raise RegistryUnavailableError(
            f"Registry API unavailable (HTTP {code})."
        ) from None
    except (urllib.error.URLError, socket.timeout, TimeoutError, OSError):
        raise RegistryUnavailableError(
            "Registry API unavailable: connection failed or timed out."
        ) from None


def _json(raw: bytes, what: str):
    try:
        return json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, ValueError):
        raise RegistryIntegrityError(f"Registry returned malformed {what}.") from None


def list_meshes() -> list[str]:
    """Return mesh names from ``GET /meshes``."""
    data = _json(_get("/meshes"), "mesh list")
    if not isinstance(data, list) or not all(isinstance(n, str) for n in data):
        raise RegistryIntegrityError("Registry returned a malformed mesh list.")
    return list(data)


def fetch_hash(name: str) -> str:
    """Return the SHA-256 hex digest from ``GET /meshes/{name}/hash``."""
    _check_name(name)
    data = _json(_get(f"/meshes/{urllib.parse.quote(name, safe='')}/hash"), "hash")
    if (
        not isinstance(data, dict)
        or set(data) != {"sha256"}
        or not isinstance(data["sha256"], str)
        or not _HEX64_RE.match(data["sha256"])
    ):
        raise RegistryIntegrityError("Registry returned a malformed hash.")
    return data["sha256"].lower()


def fetch_mesh_bytes(name: str) -> bytes:
    """Return the file bytes of a mesh after checking their SHA-256 hash."""
    _check_name(name)
    expected = fetch_hash(name)
    payload = _get(f"/meshes/{urllib.parse.quote(name, safe='')}")
    if hashlib.sha256(payload).hexdigest() != expected:
        raise RegistryIntegrityError(
            "Registry payload failed hash check; mesh rejected."
        )
    return payload


def download_mesh(name: str, dest) -> Path:
    """Write the verified mesh bytes to ``dest`` and return its path."""
    payload = fetch_mesh_bytes(name)
    path = Path(dest)
    path.write_bytes(payload)
    return path


def fetch_mesh(
    name: str,
    *,
    compute_layers: bool = True,
    compute_adjacencies: Optional[bool] = None,
) -> "CHILmesh":
    """Download a verified fort.14 mesh and return it as a ``CHILmesh``."""
    from .CHILmesh import CHILmesh

    payload = fetch_mesh_bytes(name)
    fd, tmp = tempfile.mkstemp(suffix=".14")
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(payload)
        return CHILmesh.read_from_fort14(
            Path(tmp),
            compute_layers=compute_layers,
            compute_adjacencies=compute_adjacencies,
        )
    finally:
        try:
            os.unlink(tmp)
        except OSError:
            pass
