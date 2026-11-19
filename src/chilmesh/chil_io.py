"""``.chil`` container: a zip archive holding a manifest plus ADCIRC parts.

EXPERIMENTAL. ``format_version = 0``; the layout may change in later minor
releases. This module imports only the standard library.

See ``docs/CHIL_FORMAT.md`` for the format and its rules.
"""
from __future__ import annotations

__all__ = [
    "pack",
    "rewrite_chil",
    "validate_chil",
    "read_member",
    "summary_chil",
    "save_chil",
    "load_chil",
    "ChilFormatError",
    "ChilLimits",
    "ChilReport",
    "PartInfo",
    "MemberSpec",
]

import collections.abc
import contextlib
import dataclasses
import hashlib
import io
import json
import math
import os
import re
import stat
import tempfile
import unicodedata
import zipfile
import zlib
from typing import Any

_UNSET_UMASK: int | None = None
_MAX_JSON_DEPTH = 64

_CORE_KEYS = (
    "format_version", "crs", "parts", "domain_boundary", "summary",
    "source_archive_sha256", "ext",
)
_PART_KEYS = ("name", "role", "kind", "size", "sha256")
_SUMMARY_KEYS = ("n_nodes", "n_elements", "n_triangles", "n_quads", "grid_name", "bbox")
_KINDS = ("fort14", "fort13", "fort15", "domain_boundary", "ext_member")
_ROLES = ("source", "generated", "extension")
_NS_RE = re.compile(r"^[a-z][a-z0-9_]{0,63}$")
_SHA_RE = re.compile(r"^[0-9a-f]{64}$")
_DRIVE_RE = re.compile(r"^[A-Za-z]:")


class ChilFormatError(ValueError):
    """Refusal raised for a malformed archive or an invalid write request.

    ``check`` is the stable slug (see ``data-model.md``); ``str(exc)`` is
    ``"<check>: <detail>"``.

    Example:
        >>> try:
        ...     raise ChilFormatError("request", "bad input")
        ... except ValueError as exc:
        ...     print(exc.check, "|", exc)
        request | request: bad input
    """

    def __init__(self, check: str, detail: str) -> None:
        super().__init__(f"{check}: {detail}")
        self.check = check
        self.detail = detail


@dataclasses.dataclass(frozen=True)
class ChilLimits:
    """Read and write limits (FR-058, FR-069).

    Example:
        >>> ChilLimits().max_members
        64
        >>> ChilLimits(max_members=8).max_members
        8
    """

    max_members: int = 64
    max_member_bytes: int = 2 * 1024**3
    max_total_bytes: int = 4 * 1024**3
    max_manifest_bytes: int = 16 * 1024**2


@dataclasses.dataclass(frozen=True)
class PartInfo:
    """One manifest part as reported by ``validate_chil``.

    Example:
        >>> PartInfo("fort.14", "source", "fort14", 10, "0" * 64).kind
        'fort14'
    """

    name: str
    role: str
    kind: str
    size: int
    sha256: str


@dataclasses.dataclass(frozen=True)
class ChilReport:
    """Result of validating or writing a ``.chil`` archive.

    Example:
        >>> import tempfile, os
        >>> from chilmesh import examples
        >>> with tempfile.TemporaryDirectory() as d:
        ...     p = os.path.join(d, "a.chil")
        ...     examples.annulus().save(p)
        ...     rep = validate_chil(p)
        >>> (rep.format_version, rep.parts[0].name, len(rep.archive_sha256))
        (0, 'fort.14', 64)
    """

    path: str
    format_version: int
    crs: str
    parts: tuple
    domain_boundary: str | None
    summary: dict | None
    ext_namespaces: tuple
    source_archive_sha256: str | None
    archive_sha256: str


@dataclasses.dataclass(frozen=True)
class MemberSpec:
    """One member input for ``rewrite_chil``: a regular-file path or bytes.

    Example:
        >>> spec = MemberSpec(b"polygon", name="shapes.shp")
        >>> (spec.source, spec.name)
        (b'polygon', 'shapes.shp')
    """

    source: Any
    name: str | None = None


# ----------------------------------------------------------- member names

def _check_member_name(name: Any) -> str:
    """Return ``name`` if it passes the FR-051 rules, else raise ``member-name``."""
    def bad(why: str) -> ChilFormatError:
        return ChilFormatError("member-name", f"unsafe member name {name!r}: {why}")

    if not isinstance(name, str):
        raise bad("not a string")
    if not name:
        raise bad("empty")
    if name.startswith("/"):
        raise bad("absolute")
    if "\\" in name:
        raise bad("backslash")
    if any(ord(c) < 0x20 or ord(c) == 0x7F for c in name):
        raise bad("control character")
    if _DRIVE_RE.match(name):
        raise bad("drive letter")
    if name.endswith("/"):
        raise bad("trailing slash")
    if any(seg in ("", ".", "..") for seg in name.split("/")):
        raise bad("empty, '.' or '..' segment")
    try:
        encoded = name.encode("utf-8")
    except UnicodeEncodeError:
        raise bad("not encodable as UTF-8") from None
    if len(encoded) > 255:
        raise bad("longer than 255 UTF-8 bytes")
    return name


_MAX_DB_BASENAME_BYTES = 200


def _check_db_basename(base: str) -> None:
    """Domain Boundary basename: at most 200 UTF-8 bytes (FR-067); else ``member-name``."""
    if len(base.encode("utf-8")) > _MAX_DB_BASENAME_BYTES:
        raise ChilFormatError(
            "member-name",
            f"unsafe member name {base!r}: basename longer than {_MAX_DB_BASENAME_BYTES} UTF-8 bytes")


def _fspath(path, what: str) -> str:
    """``os.fspath`` that refuses a ``bytes`` path with ``request`` (L4)."""
    p = os.fspath(path)
    if isinstance(p, bytes):
        raise ChilFormatError("request", f"{what}: a bytes path is not accepted; pass str or PathLike")
    return p


def _check_crs_request(crs) -> None:
    """Caller-supplied ``crs``: 1 to 256 characters, no control characters (L6)."""
    if (not isinstance(crs, str) or not 1 <= len(crs) <= 256
            or any(ord(c) < 0x20 or ord(c) == 0x7F for c in crs)):
        raise ChilFormatError(
            "request", "crs must be a 1 to 256 character string without control characters")


def _collision_key(name: str) -> str:
    """Key under which two member names collide (FR-052)."""
    return unicodedata.normalize("NFC", unicodedata.normalize("NFC", name).casefold())


# ------------------------------------------------------------ strict JSON

def _no_dup_pairs(pairs):
    out = {}
    for k, v in pairs:
        if k in out:
            raise ChilFormatError("duplicate-key", f"duplicate key {k!r}")
        out[k] = v
    return out


def _no_constant(name):
    raise ChilFormatError("non-finite", f"non-finite constant {name}")


def _finite_float(s):
    v = float(s)
    if not math.isfinite(v):
        raise ChilFormatError("non-finite", f"non-finite number {s}")
    return v


def _check_strings(obj) -> None:
    """Refuse lone surrogates and excessive nesting in the parsed document.

    Iterative, so a deeply nested document cannot exhaust the stack here;
    interpreters whose JSON parser accepts any depth are held to
    ``_MAX_JSON_DEPTH`` explicitly.
    """
    stack = [(obj, 1)]
    while stack:
        cur, depth = stack.pop()
        if depth > _MAX_JSON_DEPTH:
            raise ChilFormatError("manifest-json", f"nesting deeper than {_MAX_JSON_DEPTH}")
        if isinstance(cur, str):
            try:
                cur.encode("utf-8")
            except UnicodeEncodeError:
                raise ChilFormatError("manifest-json", "string is not valid Unicode") from None
        elif isinstance(cur, dict):
            for k, v in cur.items():
                stack.append((k, depth + 1))
                stack.append((v, depth + 1))
        elif isinstance(cur, list):
            stack.extend((v, depth + 1) for v in cur)


def _loads_strict(data: bytes) -> dict:
    """Parse manifest bytes with the strict rules of R4 (FR-004, FR-054)."""
    try:
        text = bytes(data).decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ChilFormatError("manifest-json", f"not valid UTF-8: {exc}") from None
    if text.startswith("﻿"):
        raise ChilFormatError("manifest-json", "byte order mark is not allowed")
    try:
        obj = json.loads(
            text,
            object_pairs_hook=_no_dup_pairs,
            parse_constant=_no_constant,
            parse_float=_finite_float,
        )
        _check_strings(obj)
    except ChilFormatError:
        raise
    except (ValueError, RecursionError) as exc:
        raise ChilFormatError("manifest-json", f"invalid JSON: {exc}") from None
    if not isinstance(obj, dict):
        raise ChilFormatError("manifest-json", "top level is not an object")
    return obj


def _check_request_json(obj, what: str, base_depth: int = 1) -> None:
    """Strict caller-input check for a value headed into the manifest (review M4).

    ``base_depth`` is the manifest depth of ``obj`` itself, so a value that passes
    here also passes the reader's nesting limit when re-parsed.
    """
    stack = [(obj, base_depth)]
    while stack:
        cur, depth = stack.pop()
        if depth > _MAX_JSON_DEPTH:
            raise ChilFormatError("request", f"{what}: nesting deeper than {_MAX_JSON_DEPTH}")
        if cur is None or isinstance(cur, bool):
            continue
        if isinstance(cur, int):
            try:
                str(cur)
            except ValueError:
                raise ChilFormatError("request", f"{what}: integer is too large") from None
            continue
        if isinstance(cur, float):
            if not math.isfinite(cur):
                raise ChilFormatError("request", f"{what}: non-finite number")
        elif isinstance(cur, str):
            try:
                cur.encode("utf-8")
            except UnicodeEncodeError:
                raise ChilFormatError("request", f"{what}: string is not valid Unicode") from None
        elif isinstance(cur, dict):
            for k, v in cur.items():
                if not isinstance(k, str):
                    raise ChilFormatError("request", f"{what}: key {k!r} is not a string")
                stack.append((k, depth + 1))
                stack.append((v, depth + 1))
        elif isinstance(cur, (list, tuple)):
            stack.extend((v, depth + 1) for v in cur)
        else:
            raise ChilFormatError(
                "request", f"{what}: {type(cur).__name__} value is not JSON serializable")


def _dumps_manifest(obj: dict) -> bytes:
    """Deterministic manifest bytes: sorted keys, 2-space indent, trailing newline."""
    try:
        text = json.dumps(obj, sort_keys=True, ensure_ascii=True, allow_nan=False,
                          indent=2, separators=(",", ": "))
    except ValueError as exc:
        raise ChilFormatError("non-finite", str(exc)) from None
    return (text + "\n").encode("utf-8")


# --------------------------------------------------------- manifest schema

def _is_int(v) -> bool:
    return type(v) is int


def _schema(detail: str) -> ChilFormatError:
    return ChilFormatError("manifest-schema", detail)


def _check_summary(s) -> None:
    if not isinstance(s, dict):
        raise _schema("summary is not an object")
    for k in s:
        if k not in _SUMMARY_KEYS:
            raise ChilFormatError("unknown-key", f"unknown summary key {k!r}")
    for k in _SUMMARY_KEYS:
        if k != "bbox" and k not in s:
            raise _schema(f"summary lacks {k!r}")
    for k in ("n_nodes", "n_elements", "n_triangles", "n_quads"):
        if not _is_int(s[k]) or s[k] < 0:
            raise _schema(f"summary.{k} is not a non-negative integer")
    if s["n_triangles"] + s["n_quads"] != s["n_elements"]:
        raise _schema("summary n_triangles + n_quads != n_elements")
    if not isinstance(s["grid_name"], str):
        raise _schema("summary.grid_name is not a string")
    if "bbox" in s:
        b = s["bbox"]
        if (not isinstance(b, list) or len(b) != 4
                or any(type(x) not in (int, float) or not math.isfinite(x) for x in b)):
            raise _schema("summary.bbox is not four finite numbers")
        if b[0] > b[2] or b[1] > b[3]:
            raise _schema("summary.bbox is not [xmin, ymin, xmax, ymax]")


def _check_ext_member_name(name: str) -> str:
    """Return the owner of ``ext/<owner>/<rest>`` or raise ``ext-member-name``."""
    segs = name.split("/")
    if (len(segs) < 3 or segs[0] != "ext" or not _NS_RE.match(segs[1])
            or any(seg in ("", ".", "..") for seg in segs[2:])):
        raise ChilFormatError("ext-member-name",
                              f"ext member {name!r} is not named ext/<owner>/<rest>")
    return segs[1]


_FIXED_NAMES = {"fort14": "fort.14", "fort13": "fort.13", "fort15": "fort.15"}


def _check_name_binding(kind: str, name: str) -> None:
    """Bind a non-ext part name to its kind (review H1, L3); ``member-name`` on a mismatch."""
    segs = name.split("/")
    if kind in _FIXED_NAMES:
        if name != _FIXED_NAMES[kind]:
            raise ChilFormatError(
                "member-name", f"a {kind} part must be named {_FIXED_NAMES[kind]!r}, got {name!r}")
    elif kind == "domain_boundary":
        if len(segs) != 2 or segs[0] != "domain_boundary":
            raise ChilFormatError(
                "member-name",
                f"a domain_boundary part must be named domain_boundary/<basename>, got {name!r}")
    if segs[0] == "ext":
        raise ChilFormatError("member-name", f"a {kind} part must not be named under ext/: {name!r}")


def _validate_manifest(m: dict, limits: ChilLimits | None = None) -> None:
    """Validate a parsed manifest (check order steps 5 to 7, minus entry matching).

    Raises ``ChilFormatError``; returns ``None`` when the manifest is valid.
    """
    limits = limits or ChilLimits()
    if "c1_version" in m:
        raise ChilFormatError(
            "prototype-container",
            "this is a pre-release prototype container (c1_version); it is not "
            "format_version 0 and cannot be read")
    if "format_version" not in m or not _is_int(m["format_version"]) or m["format_version"] != 0:
        raise ChilFormatError("format-version",
                              f"format_version must be the integer 0, got {m.get('format_version')!r}")
    for k in m:
        if k not in _CORE_KEYS:
            raise ChilFormatError("unknown-key", f"unknown top-level key {k!r}")

    crs = m.get("crs")
    if (not isinstance(crs, str) or not 1 <= len(crs) <= 256
            or any(ord(c) < 0x20 or ord(c) == 0x7F for c in crs)):
        raise _schema("crs must be a 1 to 256 character string without control characters")

    parts = m.get("parts")
    if not isinstance(parts, list) or not 1 <= len(parts) <= max(limits.max_members - 1, 1):
        raise _schema("parts must be a list of 1 to max_members - 1 entries")

    if "source_archive_sha256" in m:
        v = m["source_archive_sha256"]
        if not isinstance(v, str) or not _SHA_RE.match(v):
            raise _schema("source_archive_sha256 is not 64 lowercase hex characters")
    if "summary" in m:
        _check_summary(m["summary"])

    ext = m.get("ext", {})
    if not isinstance(ext, dict):
        raise ChilFormatError("ext-table", "ext is not an object")
    for ns, table in ext.items():
        if not _NS_RE.match(ns) or ns in _CORE_KEYS:
            raise ChilFormatError("ext-namespace", f"bad extension namespace {ns!r}")
        if not isinstance(table, dict):
            raise ChilFormatError("ext-table", f"ext.{ns} is not an object")

    seen = {_collision_key("manifest.json")}
    by_kind: dict[str, list[dict]] = {}
    for p in parts:
        if not isinstance(p, dict):
            raise _schema("a part is not an object")
        for k in p:
            if k not in _PART_KEYS:
                raise ChilFormatError("unknown-key", f"unknown part key {k!r}")
        for k in _PART_KEYS:
            if k not in p:
                raise _schema(f"part lacks {k!r}")
        if not all(isinstance(p[k], str) for k in ("name", "role", "kind", "sha256")):
            raise _schema("part name, role, kind and sha256 must be strings")
        if not _is_int(p["size"]) or p["size"] < 0:
            raise _schema("part size is not a non-negative integer")
        if not _SHA_RE.match(p["sha256"]):
            raise _schema("part sha256 is not 64 lowercase hex characters")
        if p["kind"] not in _KINDS:
            raise ChilFormatError("part-kind", f"unknown part kind {p['kind']!r}")
        if p["role"] not in _ROLES:
            raise ChilFormatError("part-role", f"unknown part role {p['role']!r}")
        kind, role = p["kind"], p["role"]
        if (role == "extension") != (kind == "ext_member"):
            raise ChilFormatError("part-role", f"role {role!r} does not pair with kind {kind!r}")
        if kind != "fort14" and kind != "ext_member" and role != "source":
            raise ChilFormatError("part-role", f"kind {kind!r} must have role 'source'")
        if kind == "ext_member":
            _check_member_name(p["name"])  # safety first (R5)
            owner = _check_ext_member_name(p["name"])
            if owner not in ext:
                raise ChilFormatError("ext-member-owner",
                                      f"ext member {p['name']!r} has no ext.{owner} table")
        else:
            _check_member_name(p["name"])
        key = _collision_key(p["name"])
        if key in seen:
            raise ChilFormatError("duplicate-name", f"member name {p['name']!r} collides")
        seen.add(key)
        by_kind.setdefault(kind, []).append(p)

    if len(by_kind.get("fort14", [])) != 1:
        raise ChilFormatError("part-count", "exactly one fort14 part is required")
    for kind in ("fort13", "fort15", "domain_boundary"):
        if len(by_kind.get(kind, [])) > 1:
            raise ChilFormatError("part-count", f"more than one {kind} part")
    if by_kind["fort14"][0]["role"] == "generated":
        for p in parts:
            if p["kind"] not in ("fort14", "ext_member"):
                raise ChilFormatError(
                    "part-count", "a generated archive may only add ext members")

    db_parts = by_kind.get("domain_boundary", [])
    if "domain_boundary" in m:
        if not isinstance(m["domain_boundary"], str):
            raise ChilFormatError("domain-boundary", "domain_boundary is not a string")
        if not db_parts or db_parts[0]["name"] != m["domain_boundary"]:
            raise ChilFormatError(
                "domain-boundary", "domain_boundary does not name a domain_boundary part")
    elif db_parts:
        raise ChilFormatError("domain-boundary", "domain_boundary part without a domain_boundary key")

    # Names are bound to kinds last, so structural slugs keep their priority.
    for p in parts:
        if p["kind"] != "ext_member":
            _check_name_binding(p["kind"], p["name"])


# ---------------------------------------------------------------- writing

def _write_archive_bytes(entries) -> bytes:
    """Deterministic zip of ``(name, bytes)`` pairs in the order given (R3)."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for name, data in entries:
            zi = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_DEFLATED
            zi.create_system = 3
            zi.external_attr = 0o644 << 16
            zi.comment = b""
            zi.extra = b""
            zf.writestr(zi, data, compress_type=zipfile.ZIP_DEFLATED, compresslevel=6)
    return buf.getvalue()


def _umask() -> int:
    global _UNSET_UMASK
    if _UNSET_UMASK is None:
        _UNSET_UMASK = os.umask(0)
        os.umask(_UNSET_UMASK)
    return _UNSET_UMASK


def _atomic_write(dst, data: bytes) -> None:
    """Write ``data`` to ``dst`` through a temp file and ``os.replace`` (R7)."""
    dst = os.fspath(dst)
    d = os.path.dirname(os.path.abspath(dst))
    fd, tmp = tempfile.mkstemp(prefix="." + os.path.basename(dst) + ".", suffix=".tmp", dir=d)
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(data)
            f.flush()
            os.fsync(f.fileno())
        os.chmod(tmp, 0o666 & ~_umask())
        os.replace(tmp, dst)
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(tmp)
        raise


# ----------------------------------------------------------------- reading

_CHUNK = 1 << 20


class _Fort14Scanner:
    """Line scanner over a packed ``fort.14`` fed in chunks (R6, FR-057, FR-063).

    Checks element rows for node ids outside ``1..n_nodes`` (id 0 has its own
    slug). The boundary block after the elements is scanned line by line, as it
    streams in and in the reader's order, for the same ids up to the first row
    that does not parse; rows after that row are not scanned. Nothing is
    buffered, so memory does not depend on the size of the block.
    """

    def __init__(self) -> None:
        self._pending: list[bytes] = []  # pieces of the unfinished last line
        self._state = 0  # 0 name, 1 header, 2 node rows, 3 element rows, 4 tail
        self._n_elems = 0
        self._n_nodes = 0
        self._left = 0
        # boundary state machine: 0 block count, 1 total, 2 segment header,
        # 3 node row, 4 done (reader stops or block complete)
        self._b = 0
        self._label = 0  # 0 open, 1 land
        self._blocks_left = 0
        self._seg = 0
        self._seg_left = 0

    def feed(self, chunk: bytes) -> None:
        """Take the next chunk; only new data is searched for line ends (M3).

        A line may end in ``\\n``, ``\\r\\n`` or ``\\r`` (the reader reads text
        with universal newlines). A trailing ``\\r`` is held back, since the next
        chunk may begin with its ``\\n``. The unfinished tail is kept as a list of
        pieces and joined once, when its line ends, so a very long line costs
        time linear in its length.
        """
        if not chunk:
            return
        end = len(chunk) - 1 if chunk.endswith(b"\r") else len(chunk)
        cut = max(chunk.rfind(b"\n", 0, end), chunk.rfind(b"\r", 0, end))
        if cut < 0:
            self._pending.append(chunk)
            return
        head = chunk[:cut + 1]
        if self._pending:
            self._pending.append(head)
            head = b"".join(self._pending)
        rest = chunk[cut + 1:]
        self._pending = [rest] if rest else []
        self._lines(head.splitlines())

    def _lines(self, lines: list) -> None:
        i = 0
        n = len(lines)
        while i < n:
            st = self._state
            if st == 2:
                take = min(self._left, n - i)
                self._left -= take
                i += take
                if self._left == 0:
                    self._enter_elements()
                continue
            ln = lines[i]
            i += 1
            if st == 0:
                self._state = 1
            elif st == 1:
                try:
                    p = ln.decode("utf-8").split()
                    self._n_elems, self._n_nodes = int(p[0]), int(p[1])
                except (ValueError, IndexError):
                    raise ChilFormatError("fort14-parse", "bad fort.14 header line") from None
                self._left = self._n_nodes
                if self._left > 0:
                    self._state = 2
                else:
                    self._enter_elements()
            elif st == 3:
                self._element(ln)
            elif self._b != 4:
                self._boundary_line(ln)

    def _enter_elements(self) -> None:
        self._left = self._n_elems
        self._state = 3 if self._left > 0 else 4

    def _element(self, ln: bytes) -> None:
        try:
            p = ln.decode("utf-8").split()
            elem_id = int(p[0])
            k = int(p[1])
            ids = [int(float(x)) for x in p[2:2 + k]]
        except (ValueError, IndexError, OverflowError):
            raise ChilFormatError("fort14-parse", "bad fort.14 element row") from None
        if k not in (3, 4):
            raise ChilFormatError(
                "fort14-parse",
                f"element {elem_id} has {k} nodes; .chil carries triangles and quads only")
        if len(p) - 2 != k:
            raise ChilFormatError(
                "fort14-parse",
                f"element {elem_id} declares {k} nodes but its row has {len(p) - 2} ids")
        if any(v == 0 for v in ids):
            raise ChilFormatError("node-id-zero", f"element {elem_id} references node id 0")
        bad = sorted({v for v in ids if v < 1 or v > self._n_nodes})
        if bad:
            raise ChilFormatError(
                "node-id-range",
                f"element {elem_id} references node id(s) {bad}; valid ids are 1..{self._n_nodes}")
        self._left -= 1
        if self._left == 0:
            self._state = 4

    def finish(self) -> None:
        if self._pending:
            data = b"".join(self._pending)
            self._pending = []
            self._lines(data.splitlines())
        if self._state != 4:
            raise ChilFormatError("fort14-parse", "fort.14 ends before the element block does")


    def _boundary_line(self, ln: bytes) -> None:
        """Feed one boundary row (M1; same rules as the reader).

        An id read before the first row that fails to parse (or is missing) is
        checked; once a row fails to parse the scan stops, matching the reader
        keeping the segments it parsed before it stops.
        """
        try:
            tok = ln.decode("utf-8").split()
            st = self._b
            if st == 0:
                self._blocks_left = int(tok[0])
                self._b = 1
            elif st == 1:
                int(tok[0])
                self._b = 2
                self._next_segment()
            elif st == 2:
                n_seg = int(tok[0])
                if len(tok) > 1:
                    int(tok[1])
                self._seg += 1
                self._blocks_left -= 1
                self._seg_left = n_seg
                self._b = 3
                self._after_nodes()
            else:
                nid = int(tok[0])
        except (IndexError, ValueError):
            self._b = 4  # rest of the block does not parse: reader behavior (R14.4)
            return
        if st == 3:
            label = ("open", "land")[self._label]
            if nid == 0:
                raise ChilFormatError(
                    "node-id-zero", f"{label} boundary segment {self._seg} references node id 0")
            if nid < 1 or nid > self._n_nodes:
                raise ChilFormatError(
                    "node-id-range",
                    f"{label} boundary segment {self._seg} references node id {nid}; "
                    f"valid ids are 1..{self._n_nodes}")
            self._seg_left -= 1
            self._after_nodes()

    def _next_segment(self) -> None:
        """Move to the next segment header, the next block, or the end."""
        if self._blocks_left > 0:
            self._b = 2
        elif self._label == 0:
            self._label = 1
            self._seg = 0
            self._b = 0
        else:
            self._b = 4

    def _after_nodes(self) -> None:
        if self._b == 3 and self._seg_left <= 0:
            self._next_segment()


class _Budget:
    __slots__ = ("used",)

    def __init__(self) -> None:
        self.used = 0


def _stream_member(zf, info, *, limit, limit_slug, budget, max_total, sink=None):
    """Read one member in chunks, counting real bytes (R5, FR-058).

    Returns ``(byte_count, sha256_hex)``. Declared header sizes are never used.
    """
    sha = hashlib.sha256()
    n = 0
    try:
        with zf.open(info) as f:
            while True:
                chunk = f.read(_CHUNK)
                if not chunk:
                    break
                n += len(chunk)
                budget.used += len(chunk)
                if n > limit:
                    raise ChilFormatError(
                        limit_slug, f"member {info.filename!r} exceeds {limit} bytes")
                if budget.used > max_total:
                    raise ChilFormatError(
                        "limit-total-bytes", f"archive exceeds {max_total} bytes in total")
                sha.update(chunk)
                if sink is not None:
                    sink(chunk)
    except (zipfile.BadZipFile, zlib.error, EOFError, NotImplementedError,
            UnicodeDecodeError) as exc:
        raise ChilFormatError("not-zip", f"member {info.filename!r}: {exc}") from None
    return n, sha.hexdigest()


def _validate_stream(path, limits, want=None):
    """Run the R5 check order. Returns ``(report, captured_bytes_or_None)``.

    ``want`` is a member name or part kind whose verified bytes are captured.
    """
    limits = limits or ChilLimits()
    path = _fspath(path, "path")
    st = os.stat(path)
    if not stat.S_ISREG(st.st_mode):
        raise ChilFormatError("request", f"{path!r} is not a regular file")
    with open(path, "rb", buffering=0) as fh:
        # Unbuffered, so the second hash pass reads the file, not a cache.
        # One open handle for both the hash and the zip reader; the hash is
        # streamed before and after validation and must match, so the
        # reported archive_sha256 always describes the bytes that were checked.
        archive_sha = _hash_handle(fh)
        zf = _open_zip(fh)
        with zf:
            result = _validate_open(zf, path, archive_sha, limits, want)
        if _hash_handle(fh) != archive_sha:
            raise ChilFormatError(
                "archive-changed", f"{path!r} changed while it was being validated")
        return result


_REWRITE_KINDS_WORDS = ("fort14", "fort13", "fort15", "domain_boundary")


def _open_zip(fh):
    """Open ``fh`` as a zip; hostile headers surface as ``ChilFormatError`` (review M2)."""
    try:
        return zipfile.ZipFile(fh)
    except (zipfile.BadZipFile, EOFError, OSError, UnicodeDecodeError,
            NotImplementedError, ValueError) as exc:
        raise ChilFormatError("not-zip", f"not a zip archive: {exc}") from None


def _hash_handle(fh) -> str:
    """Stream the sha256 of an open binary file from offset 0, then rewind."""
    fh.seek(0)
    h = hashlib.sha256()
    for blk in iter(lambda: fh.read(_CHUNK), b""):
        h.update(blk)
    fh.seek(0)
    return h.hexdigest()


def _validate_open(zf, path, archive_sha, limits, want, entries_only=False):
    infos = zf.infolist()
    if len(infos) > limits.max_members:
        raise ChilFormatError(
            "limit-members", f"{len(infos)} members exceed the limit of {limits.max_members}")
    seen: set[str] = set()
    for info in infos:
        if info.flag_bits & 0x41:
            raise ChilFormatError("encrypted", f"member {info.filename!r} is encrypted")
        if info.compress_type not in (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED):
            raise ChilFormatError(
                "compression", f"member {info.filename!r} uses compression method {info.compress_type}")
        name = _check_member_name(info.orig_filename)
        key = _collision_key(name)
        if key in seen:
            raise ChilFormatError("duplicate-name", f"member name {name!r} collides")
        seen.add(key)
    if not infos or infos[0].filename != "manifest.json":
        if any(i.filename == "manifest.json" for i in infos):
            raise ChilFormatError("manifest-not-first", "manifest.json is not the first member")
        raise ChilFormatError("manifest-missing", "manifest.json is missing")
    if infos[0].header_offset != 0 or any(i.header_offset < infos[0].header_offset
                                          for i in infos[1:]):
        raise ChilFormatError(
            "manifest-not-first", "manifest.json is not the first entry in the archive bytes")

    budget = _Budget()
    mchunks: list[bytes] = []
    _stream_member(zf, infos[0], limit=limits.max_manifest_bytes,
                   limit_slug="limit-manifest-bytes", budget=budget,
                   max_total=limits.max_total_bytes, sink=mchunks.append)
    m = _loads_strict(b"".join(mchunks))
    _validate_manifest(m, limits)

    parts = m["parts"]
    listed = {p["name"] for p in parts}
    by_name = {i.filename: i for i in infos[1:]}
    for name in by_name:
        if name not in listed:
            raise ChilFormatError("unlisted-member", f"member {name!r} is not listed in parts")
    for p in parts:
        if p["name"] not in by_name:
            raise ChilFormatError("missing-member", f"part {p['name']!r} has no zip member")
    if entries_only:
        return m

    target = None
    if want is not None:
        if want in _REWRITE_KINDS_WORDS:
            for p in parts:
                if p["kind"] == want:
                    target = p["name"]
                    break
        else:
            for p in parts:
                if p["name"] == want:
                    target = p["name"]
                    break

    captured = None
    for p in parts:
        scanner = _Fort14Scanner() if p["kind"] == "fort14" else None
        deferred: list[ChilFormatError] = []
        keep: list[bytes] = []

        def sink(chunk, scanner=scanner, deferred=deferred, keep=keep, grab=p["name"] == target):
            if grab:
                keep.append(chunk)
            if scanner is not None and not deferred:
                try:
                    scanner.feed(chunk)
                except ChilFormatError as exc:
                    deferred.append(exc)

        n, digest = _stream_member(
            zf, by_name[p["name"]], limit=limits.max_member_bytes,
            limit_slug="limit-member-bytes", budget=budget,
            max_total=limits.max_total_bytes, sink=sink)
        if n != p["size"]:
            raise ChilFormatError(
                "size-mismatch", f"part {p['name']!r} has {n} bytes, manifest says {p['size']}")
        if digest != p["sha256"]:
            raise ChilFormatError("digest-mismatch", f"part {p['name']!r} sha256 differs")
        if scanner is not None:
            if deferred:
                raise deferred[0]
            scanner.finish()
        if p["name"] == target:
            captured = b"".join(keep)

    report = ChilReport(
        path=path,
        format_version=0,
        crs=m["crs"],
        parts=tuple(PartInfo(p["name"], p["role"], p["kind"], p["size"], p["sha256"])
                    for p in parts),
        domain_boundary=m.get("domain_boundary"),
        summary=m.get("summary"),
        ext_namespaces=tuple(sorted(m.get("ext", {}))),
        source_archive_sha256=m.get("source_archive_sha256"),
        archive_sha256=archive_sha,
    )
    return report, (captured, target)


def validate_chil(path, *, limits: ChilLimits | None = None) -> ChilReport:
    """Fully validate a ``.chil`` archive and return its report (EXPERIMENTAL).

    Raises ``ChilFormatError`` on any refusal; writes nothing.

    Example:
        >>> import tempfile, os
        >>> from chilmesh import examples
        >>> with tempfile.TemporaryDirectory() as d:
        ...     p = os.path.join(d, "a.chil")
        ...     examples.annulus().save(p)
        ...     report = validate_chil(p)
        >>> [part.kind for part in report.parts]
        ['fort14']
    """
    return _validate_stream(path, limits)[0]


def read_member(path, name_or_kind: str, *, limits: ChilLimits | None = None) -> bytes:
    """Validate the archive and return the verified bytes of one part.

    ``name_or_kind`` is a part name or one of the kinds ``fort14``, ``fort13``,
    ``fort15``, ``domain_boundary``. Anything else raises ``request``. A kind
    word always means the part of that kind; a part name never shadows it.

    Example:
        >>> import tempfile, os
        >>> from chilmesh import examples
        >>> with tempfile.TemporaryDirectory() as d:
        ...     p = os.path.join(d, "a.chil")
        ...     examples.annulus().save(p)
        ...     data = read_member(p, "fort14")
        >>> isinstance(data, bytes) and len(data) > 0
        True
    """
    _report, (captured, target) = _validate_stream(path, limits, name_or_kind)
    if target is None:
        raise ChilFormatError("request", f"no part matches {name_or_kind!r}")
    return captured


# ------------------------------------------------------------ summary

def _summary_manifest_only(path, limits) -> tuple[dict, int]:
    """Steps 1 to 7 only: opens just ``manifest.json``; returns (manifest, file bytes)."""
    limits = limits or ChilLimits()
    path = _fspath(path, "path")
    st = os.stat(path)
    if not stat.S_ISREG(st.st_mode):
        raise ChilFormatError("request", f"{path!r} is not a regular file")
    with open(path, "rb") as fh:
        zf = _open_zip(fh)
        with zf:
            m = _validate_open(zf, path, None, limits, None, entries_only=True)
    return m, st.st_size


def _bbox_dict(b) -> dict:
    return {"min_x": b[0], "max_x": b[2], "min_y": b[1], "max_y": b[3]}


def _fill_counts(result: dict, s: dict) -> None:
    result["n_nodes"] = s["n_nodes"]
    result["n_elems"] = s["n_elements"]
    result["n_triangles"] = s["n_triangles"]
    result["n_quads"] = s["n_quads"]
    result["grid_name"] = s["grid_name"]
    if s["n_quads"] == 0:
        result["element_type"] = "Triangular"
    elif s["n_triangles"] == 0:
        result["element_type"] = "Quadrilateral"
    else:
        result["element_type"] = "Mixed-Element"
    if "bbox" in s:
        result["bbox"] = _bbox_dict(s["bbox"])


def summary_chil(path, *, deep: bool = False, limits: ChilLimits | None = None) -> dict:
    """Summarize a ``.chil`` archive (EXPERIMENTAL).

    Non-deep mode reads only ``manifest.json`` (checks steps 1 to 7) and never
    parses or checksums members. ``deep=True`` validates fully, parses the
    ``fort.14`` and recomputes the summary. Raises ``ChilFormatError``.

    Example:
        >>> import tempfile, os
        >>> from chilmesh import examples
        >>> with tempfile.TemporaryDirectory() as d:
        ...     p = os.path.join(d, "a.chil")
        ...     examples.annulus().save(p)
        ...     info = summary_chil(p)
        >>> (info["verified"], info["n_nodes"] > 0)
        (False, True)
    """
    if deep:
        from .CHILmesh import CHILmesh

        report, (captured, target) = _validate_stream(path, limits, "fort14")
        if target is None or captured is None:
            raise ChilFormatError("request", "archive has no fort14 part")
        try:
            with io.TextIOWrapper(io.BytesIO(captured), encoding="utf-8") as t:
                lines = t.readlines()
            mesh = CHILmesh._from_fort14_lines(
                lines, compute_layers=False, compute_adjacencies=False)
            recomputed = _summary_fields(mesh, captured)
        except ChilFormatError:
            raise
        except (ValueError, IndexError, UnicodeDecodeError) as exc:
            raise ChilFormatError("fort14-parse", f"cannot parse fort.14: {exc}") from None
        cached = report.summary
        result = {
            "path": os.fspath(path), "format": "chil",
            "file_bytes": os.stat(path).st_size,
            "format_version": report.format_version, "crs": report.crs,
            "parts": [dataclasses.asdict(p) for p in report.parts],
            "domain_boundary": report.domain_boundary,
            "ext_namespaces": sorted(report.ext_namespaces),
            "source_archive_sha256": report.source_archive_sha256,
            "cached_summary": cached,
            "recomputed_summary": recomputed,
            "verified": True,
        }
        _fill_counts(result, recomputed)
        result["mismatches"] = sorted(
            k for k in cached if cached[k] != recomputed.get(k)) if cached else []
        return result

    m, size = _summary_manifest_only(path, limits)
    cached = m.get("summary")
    result = {
        "path": os.fspath(path), "format": "chil", "file_bytes": size,
        "format_version": m["format_version"], "crs": m["crs"],
        "parts": [{k: p[k] for k in _PART_KEYS} for p in m["parts"]],
        "domain_boundary": m.get("domain_boundary"),
        "ext_namespaces": sorted(m.get("ext", {})),
        "source_archive_sha256": m.get("source_archive_sha256"),
        "cached_summary": cached,
        "verified": False,
    }
    if cached is not None:
        _fill_counts(result, cached)
    return result


# ------------------------------------------------------------ save and load

_GRID_NAME_MAX_BYTES = 256
_LINE_END_RE = re.compile(rb"\r\n|\r|\n")


def _grid_name(f14: bytes) -> str:
    """First line of a fort.14 (ends: ``\\n``, ``\\r\\n``, ``\\r``), stripped, cut to
    256 UTF-8 bytes at a character boundary (L1, FR-030)."""
    m = _LINE_END_RE.search(f14)
    line = f14 if m is None else f14[:m.start()]
    text = line.decode("utf-8", "replace").strip()
    raw = text.encode("utf-8")
    if len(raw) > _GRID_NAME_MAX_BYTES:
        text = raw[:_GRID_NAME_MAX_BYTES].decode("utf-8", "ignore")
    return text


def _summary_fields(mesh, f14: bytes) -> dict:
    """Cached summary fields (FR-030) for ``mesh`` and its generated fort.14."""
    import numpy as np

    conn = np.asarray(mesh.connectivity_list)
    if conn.ndim == 2 and conn.shape[0] and conn.shape[1] not in (3, 4):
        raise ChilFormatError(
            "fort14-parse",
            f"element rows have {conn.shape[1]} nodes; .chil carries triangles and quads only")
    if conn.ndim == 2 and conn.shape[1] == 4:
        n_tri = int(np.count_nonzero(conn[:, 3] == conn[:, 0]))
    elif conn.ndim == 2 and conn.shape[0] and conn.shape[1] == 3:
        n_tri = int(conn.shape[0])
    else:
        n_tri = 0
    n_elems = int(conn.shape[0]) if conn.ndim == 2 else 0
    first = _grid_name(f14)
    summary = {
        "n_nodes": int(mesh.n_verts),
        "n_elements": n_elems,
        "n_triangles": n_tri,
        "n_quads": n_elems - n_tri,
        "grid_name": first,
    }
    pts = np.asarray(mesh.points)
    if pts.ndim == 2 and pts.shape[0] and pts.shape[1] >= 2:
        xy = pts[:, :2].astype(float)
        if np.isfinite(xy).all():
            summary["bbox"] = [float(xy[:, 0].min()), float(xy[:, 1].min()),
                               float(xy[:, 0].max()), float(xy[:, 1].max())]
    return summary


def _check_write_limits(sizes, limits: ChilLimits) -> None:
    """Raise the ``limit-*`` check the validator would raise (FR-069).

    ``sizes`` is the byte length of each member in archive order (manifest first).
    """
    if len(sizes) > limits.max_members:
        raise ChilFormatError(
            "limit-members", f"{len(sizes)} members exceed the limit of {limits.max_members}")
    used = 0
    for i, n in enumerate(sizes):
        if i == 0:
            if n > limits.max_manifest_bytes:
                raise ChilFormatError(
                    "limit-manifest-bytes", f"manifest exceeds {limits.max_manifest_bytes} bytes")
        elif n > limits.max_member_bytes:
            raise ChilFormatError(
                "limit-member-bytes", f"member exceeds {limits.max_member_bytes} bytes")
        used += n
        if used > limits.max_total_bytes:
            raise ChilFormatError(
                "limit-total-bytes", f"archive exceeds {limits.max_total_bytes} bytes in total")


def save_chil(mesh, dst, *, crs: str | None = None,
              limits: ChilLimits | None = None) -> ChilReport:
    """Write ``mesh`` to a ``.chil`` archive (EXPERIMENTAL).

    One generated ``fort.14`` part plus the manifest. Nothing is written when a
    limit is exceeded. Returns the report of the written archive.

    Example:
        >>> import tempfile, os
        >>> from chilmesh import examples
        >>> with tempfile.TemporaryDirectory() as d:
        ...     report = save_chil(examples.annulus(), os.path.join(d, "a.chil"), crs="EPSG:4326")
        >>> (report.crs, report.parts[0].role)
        ('EPSG:4326', 'generated')
    """
    from .CHILmesh import _write_fort14_stream

    limits = limits or ChilLimits()
    if crs is not None:
        _check_crs_request(crs)
    dst = _fspath(dst, "dst")
    buf = io.StringIO(newline="\n")
    _write_fort14_stream(mesh, buf)
    f14 = buf.getvalue().encode("utf-8")
    sha = hashlib.sha256(f14).hexdigest()
    part = {"name": "fort.14", "role": "generated", "kind": "fort14",
            "size": len(f14), "sha256": sha}
    manifest = {
        "format_version": 0,
        "crs": "unknown" if crs is None else crs,
        "parts": [part],
        "summary": _summary_fields(mesh, f14),
    }
    _validate_manifest(manifest, limits)
    mbytes = _dumps_manifest(manifest)
    _check_write_limits([len(mbytes), len(f14)], limits)
    data = _write_archive_bytes([("manifest.json", mbytes), ("fort.14", f14)])
    _atomic_write(dst, data)
    return ChilReport(
        path=os.fspath(dst),
        format_version=0,
        crs=manifest["crs"],
        parts=(PartInfo("fort.14", "generated", "fort14", len(f14), sha),),
        domain_boundary=None,
        summary=manifest["summary"],
        ext_namespaces=(),
        source_archive_sha256=None,
        archive_sha256=hashlib.sha256(data).hexdigest(),
    )


def load_chil(path, *, compute_layers: bool = True, compute_adjacencies=None,
              limits: ChilLimits | None = None):
    """Validate a ``.chil`` archive, then parse its fort.14 part (EXPERIMENTAL).

    Writes nothing to disk.

    Example:
        >>> import tempfile, os
        >>> from chilmesh import examples
        >>> with tempfile.TemporaryDirectory() as d:
        ...     p = os.path.join(d, "a.chil")
        ...     mesh = examples.annulus()
        ...     _ = save_chil(mesh, p)
        ...     back = load_chil(p, compute_layers=False)
        >>> back.n_verts == mesh.n_verts
        True
    """
    from .CHILmesh import CHILmesh

    _report, (captured, target) = _validate_stream(path, limits, "fort14")
    if target is None or captured is None:
        raise ChilFormatError("request", "archive has no fort14 part")
    try:
        with io.TextIOWrapper(io.BytesIO(captured), encoding="utf-8") as t:
            lines = t.readlines()
        return CHILmesh._from_fort14_lines(lines, compute_layers, compute_adjacencies)
    except ChilFormatError:
        raise
    except (ValueError, IndexError, UnicodeDecodeError) as exc:
        raise ChilFormatError("fort14-parse", f"cannot parse fort.14: {exc}") from None


# ------------------------------------------------------------------- pack

def pack(dst, fort14, *, fort13=None, fort15=None, domain_boundary=None,
         crs: str = "unknown", limits: ChilLimits | None = None) -> ChilReport:
    """Pack ADCIRC source files into a ``.chil`` archive, byte for byte (EXPERIMENTAL).

    Inputs are read in binary mode only and must be regular files. Each is
    hashed in a first pass and read again in a second; a change between the
    passes raises ``request``. Nothing is written when a limit is exceeded.

    Example:
        >>> import tempfile, os
        >>> from chilmesh import examples
        >>> with tempfile.TemporaryDirectory() as d:
        ...     f14 = os.path.join(d, "fort.14")
        ...     examples.annulus().save(f14 + ".14")
        ...     report = pack(os.path.join(d, "a.chil"), f14 + ".14", crs="EPSG:4326")
        >>> [(part.name, part.role) for part in report.parts]
        [('fort.14', 'source')]
    """
    limits = limits or ChilLimits()
    _check_crs_request(crs)
    dst = _fspath(dst, "dst")
    fort14 = _fspath(fort14, "fort14")
    fort13 = None if fort13 is None else _fspath(fort13, "fort13")
    fort15 = None if fort15 is None else _fspath(fort15, "fort15")
    if domain_boundary is not None:
        domain_boundary = _fspath(domain_boundary, "domain_boundary")
    inputs = [("fort.14", "fort14", fort14)]
    if fort13 is not None:
        inputs.append(("fort.13", "fort13", fort13))
    if fort15 is not None:
        inputs.append(("fort.15", "fort15", fort15))
    if domain_boundary is not None:
        base = os.path.basename(domain_boundary)
        _check_member_name(base)
        _check_db_basename(base)
        inputs.append((_check_member_name(f"domain_boundary/{base}"),
                       "domain_boundary", domain_boundary))

    paths = []
    for _name, _kind, src in inputs:
        src = os.fspath(src)
        if not stat.S_ISREG(os.stat(src).st_mode):
            raise ChilFormatError("request", f"{src!r} is not a regular file")
        paths.append(src)

    # Check every input size against the limits before any byte is read.
    _check_write_limits([0] + [os.stat(p).st_size for p in paths], limits)

    # Pass 1: hash every input; scan the fort.14 for bad node ids.
    first = []
    for (_name, kind, _src), path in zip(inputs, paths):
        sha = hashlib.sha256()
        n = 0
        scanner = _Fort14Scanner() if kind == "fort14" else None
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(_CHUNK), b""):
                n += len(chunk)
                if n > limits.max_member_bytes:
                    raise ChilFormatError(
                        "limit-member-bytes", f"member exceeds {limits.max_member_bytes} bytes")
                sha.update(chunk)
                if scanner is not None:
                    scanner.feed(chunk)
        if scanner is not None:
            scanner.finish()
        first.append((n, sha.hexdigest()))
    _check_write_limits([0] + [n for n, _ in first], limits)

    # Pass 2: read the bytes for the archive; they must match pass 1.
    blobs = []
    for (name, _kind, _src), path, (n, digest) in zip(inputs, paths, first):
        buf = bytearray()
        sha = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(_CHUNK), b""):
                buf += chunk
                if len(buf) > n:
                    break
                sha.update(chunk)
        if len(buf) != n or sha.hexdigest() != digest:
            raise ChilFormatError("request", f"input {path!r} changed while it was being packed")
        blobs.append(bytes(buf))

    from .CHILmesh import CHILmesh

    try:
        with io.TextIOWrapper(io.BytesIO(blobs[0]), encoding="utf-8") as t:
            lines = t.readlines()
        mesh = CHILmesh._from_fort14_lines(
            lines, compute_layers=False, compute_adjacencies=False)
        summary = _summary_fields(mesh, blobs[0])
    except ChilFormatError:
        raise
    except (ValueError, IndexError, UnicodeDecodeError) as exc:
        raise ChilFormatError("fort14-parse", f"cannot parse fort.14: {exc}") from None

    parts = [{"name": name, "role": "source", "kind": kind, "size": n, "sha256": digest}
             for (name, kind, _src), (n, digest) in zip(inputs, first)]
    manifest = {"format_version": 0, "crs": crs, "parts": parts, "summary": summary}
    if domain_boundary is not None:
        manifest["domain_boundary"] = parts[-1]["name"]
    _validate_manifest(manifest, limits)
    mbytes = _dumps_manifest(manifest)
    _check_write_limits([len(mbytes)] + [len(b) for b in blobs], limits)
    data = _write_archive_bytes(
        [("manifest.json", mbytes)] + [(p["name"], b) for p, b in zip(parts, blobs)])
    _atomic_write(dst, data)
    return ChilReport(
        path=os.fspath(dst),
        format_version=0,
        crs=crs,
        parts=tuple(PartInfo(p["name"], p["role"], p["kind"], p["size"], p["sha256"])
                    for p in parts),
        domain_boundary=manifest.get("domain_boundary"),
        summary=summary,
        ext_namespaces=(),
        source_archive_sha256=None,
        archive_sha256=hashlib.sha256(data).hexdigest(),
    )


# ---------------------------------------------------------------- rewrite

_REWRITE_KINDS = ("fort14", "fort13", "fort15", "domain_boundary")


def _as_spec(value) -> MemberSpec:
    return value if isinstance(value, MemberSpec) else MemberSpec(value)


def _plan_member(key: str, value):
    """Resolve one ``set_parts`` entry to ``(name, kind, role, source)``."""
    spec = _as_spec(value)
    src = spec.source
    is_bytes = isinstance(src, (bytes, bytearray, memoryview))
    if not is_bytes and not isinstance(src, (str, os.PathLike)):
        raise ChilFormatError("request", f"member source for {key!r} is not a path or bytes")
    if key in _FIXED_NAMES:
        if spec.name is not None:
            raise ChilFormatError("request", f"a {key} part has a fixed name; MemberSpec.name must be None")
        return _FIXED_NAMES[key], key, "source", src
    if key == "domain_boundary":
        if spec.name is None:
            if is_bytes:
                raise ChilFormatError(
                    "request", "domain_boundary bytes need MemberSpec(name=...); no basename is known")
            base = os.path.basename(_fspath(src, "domain_boundary"))
        else:
            base = spec.name
        _check_member_name(base)
        _check_db_basename(base)
        if "/" in base:
            raise ChilFormatError("member-name", f"unsafe member name {base!r}: contains '/'")
        return _check_member_name(f"domain_boundary/{base}"), "domain_boundary", "source", src
    _check_ext_member_name(key)
    _check_member_name(key)
    if spec.name is not None and spec.name != key:
        raise ChilFormatError("request", f"MemberSpec.name {spec.name!r} differs from the member name {key!r}")
    return key, "ext_member", "extension", src


def _input_size(src) -> int:
    """Size of a member input from ``os.stat`` (or ``len`` for bytes); reads nothing (M5)."""
    if isinstance(src, (bytes, bytearray, memoryview)):
        return len(src)
    path = _fspath(src, "member source")
    st = os.stat(path)
    if not stat.S_ISREG(st.st_mode):
        raise ChilFormatError("request", f"{path!r} is not a regular file")
    return st.st_size


def _read_input(src, limit: int):
    """Bytes of a member input: stat check, binary open, two passes (FR-066).

    Reading stops with ``limit-member-bytes`` once more than ``limit`` bytes
    have been seen, so a file that grows after the stat check stays bounded.
    """
    if isinstance(src, (bytes, bytearray, memoryview)):
        return bytes(src)
    path = _fspath(src, "member source")
    if not stat.S_ISREG(os.stat(path).st_mode):
        raise ChilFormatError("request", f"{path!r} is not a regular file")
    sha = hashlib.sha256()
    n = 0
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            n += len(chunk)
            if n > limit:
                raise ChilFormatError("limit-member-bytes", f"member exceeds {limit} bytes")
            sha.update(chunk)
    digest = sha.hexdigest()
    buf = bytearray()
    sha = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            buf += chunk
            sha.update(chunk)
    if len(buf) != n or sha.hexdigest() != digest:
        raise ChilFormatError("request", f"input {path!r} changed while it was being read")
    return buf


def rewrite_chil(src, dst, *, ext=None, set_parts=None, remove_parts=(),
                 crs: str | None = None, limits: ChilLimits | None = None) -> ChilReport:
    """Rewrite a ``.chil`` archive, keeping what is not changed (EXPERIMENTAL).

    ``src`` is fully validated first. ``ext`` maps a namespace to a table
    (``None`` removes it). ``set_parts`` adds or replaces ``fort13``,
    ``fort15``, ``domain_boundary`` or ``ext/<owner>/<rest>`` members;
    ``remove_parts`` removes them. The ``fort14`` part cannot be added,
    replaced or removed (``request``): replacing the mesh would leave extension
    tables, signatures and companions describing a mesh that is no longer in
    the archive. ``dst == src`` is allowed. Nothing is written on refusal.

    Caller input is checked as strictly as the reader checks a manifest: ``ext``
    tables must be JSON (string keys, finite numbers, valid Unicode, nesting at
    most 64) and the new manifest is re-parsed and re-validated before any
    byte is written; bad input raises ``request``.

    Example:
        >>> import tempfile, os
        >>> from chilmesh import examples
        >>> with tempfile.TemporaryDirectory() as d:
        ...     a, b = os.path.join(d, "a.chil"), os.path.join(d, "b.chil")
        ...     examples.annulus().save(a)
        ...     report = rewrite_chil(a, b, ext={"demo": {"note": "hi"}})
        >>> report.ext_namespaces
        ('demo',)
    """
    limits = limits or ChilLimits()
    if crs is not None:
        _check_crs_request(crs)
    dst = _fspath(dst, "dst")
    src_path = _fspath(src, "src")
    if not stat.S_ISREG(os.stat(src_path).st_mode):
        raise ChilFormatError("request", f"{src_path!r} is not a regular file")

    # Validate src fully; keep its manifest and the verified bytes of every part.
    with open(src_path, "rb", buffering=0) as fh:
        archive_sha = _hash_handle(fh)
        zf = _open_zip(fh)
        with zf:
            _validate_open(zf, src_path, archive_sha, limits, None)
            infos = zf.infolist()
            mchunks: list[bytes] = []
            _stream_member(zf, infos[0], limit=limits.max_manifest_bytes,
                           limit_slug="limit-manifest-bytes", budget=_Budget(),
                           max_total=limits.max_total_bytes, sink=mchunks.append)
            old = _loads_strict(b"".join(mchunks))
            by_name = {i.filename: i for i in infos[1:]}
            old_bytes: dict[str, bytes] = {}
            for p in old["parts"]:
                chunks: list[bytes] = []
                _n, digest = _stream_member(
                    zf, by_name[p["name"]], limit=limits.max_member_bytes,
                    limit_slug="limit-member-bytes", budget=_Budget(),
                    max_total=limits.max_total_bytes, sink=chunks.append)
                if digest != p["sha256"]:
                    raise ChilFormatError("digest-mismatch", f"part {p['name']!r} sha256 differs")
                old_bytes[p["name"]] = b"".join(chunks)
        if _hash_handle(fh) != archive_sha:
            raise ChilFormatError(
                "archive-changed", f"{src_path!r} changed while it was being validated")

    # Request checks.
    if ext is not None and not isinstance(ext, dict):
        raise ChilFormatError("request", "ext must be a dict of namespace to table")
    if set_parts is not None and not isinstance(set_parts, dict):
        raise ChilFormatError("request", "set_parts must be a dict")
    if not isinstance(remove_parts, (list, tuple)):
        raise ChilFormatError("request", "remove_parts must be a list or tuple of names")
    new_ext = dict(old.get("ext", {}))
    for ns, table in (ext or {}).items():
        if not isinstance(ns, str) or not _NS_RE.match(ns) or ns in _CORE_KEYS:
            raise ChilFormatError("ext-namespace", f"bad extension namespace {ns!r}")
        if table is None:
            new_ext.pop(ns, None)
        elif isinstance(table, dict):
            _check_request_json(table, f"ext.{ns}", base_depth=3)
            new_ext[ns] = dict(table)
        else:
            raise ChilFormatError("ext-table", f"ext.{ns} is not an object")

    sets = dict(set_parts or {})
    removes = list(remove_parts)
    for key in sets:
        if key == "fort14":
            raise ChilFormatError(
                "request", "the fort14 part cannot be added or replaced; use pack or save for a new mesh")
        if not isinstance(key, str):
            raise ChilFormatError("request", f"set_parts key {key!r} is not a string")
    for key in removes:
        if not isinstance(key, str):
            raise ChilFormatError("request", f"remove_parts entry {key!r} is not a string")
        if key == "fort14":
            raise ChilFormatError("request", "the fort14 part cannot be removed")
        if key in sets:
            raise ChilFormatError("request", f"{key!r} is both set and removed")
    plan = {key: _plan_member(key, value) for key, value in sets.items()}

    by_kind = {p["kind"]: p for p in old["parts"] if p["kind"] != "ext_member"}
    ext_names = {p["name"] for p in old["parts"] if p["kind"] == "ext_member"}
    drop: set[str] = set()
    for key in removes:
        if key in _REWRITE_KINDS:
            if key not in by_kind:
                raise ChilFormatError("request", f"no {key} part to remove")
            drop.add(by_kind[key]["name"])
        else:
            if not isinstance(key, str):
                raise ChilFormatError("request", f"remove_parts entry {key!r} is not a string")
            _check_ext_member_name(key)
            if key not in ext_names:
                raise ChilFormatError("request", f"no part named {key!r} to remove")
            drop.add(key)

    replace_at: dict[str, tuple] = {}
    appended = []
    for key, entry in plan.items():
        name, kind = entry[0], entry[1]
        existing = by_kind.get(kind)["name"] if kind in by_kind else (
            name if name in ext_names else None)
        (replace_at.__setitem__(existing, entry) if existing is not None
         else appended.append(entry))

    # Check member input sizes from os.stat before any input is read (M5).
    sizes = [0]
    for p in old["parts"]:
        if p["name"] in drop:
            continue
        sizes.append(_input_size(replace_at[p["name"]][3]) if p["name"] in replace_at
                     else p["size"])
    sizes.extend(_input_size(entry[3]) for entry in appended)
    _check_write_limits(sizes, limits)

    # Assemble the new parts; read member inputs only after every request check.
    parts: list[dict] = []
    blobs: list[bytes] = []

    def add(name, kind, role, data):
        parts.append({"name": name, "role": role, "kind": kind, "size": len(data),
                      "sha256": hashlib.sha256(data).hexdigest()})
        blobs.append(data)

    for p in old["parts"]:
        if p["name"] in drop:
            continue
        if p["name"] in replace_at:
            name, kind, role, source = replace_at[p["name"]]
            add(name, kind, role, _read_input(source, limits.max_member_bytes))
        else:
            parts.append(dict(p))
            blobs.append(old_bytes[p["name"]])
    for name, kind, role, source in appended:
        add(name, kind, role, _read_input(source, limits.max_member_bytes))

    _check_write_limits([0] + [len(b) for b in blobs], limits)

    manifest: dict = {
        "format_version": 0,
        "crs": old["crs"] if crs is None else crs,
        "parts": parts,
        "source_archive_sha256": archive_sha,
    }
    if "summary" in old:
        manifest["summary"] = old["summary"]
    if new_ext:
        manifest["ext"] = new_ext
    db = [p for p in parts if p["kind"] == "domain_boundary"]
    if db:
        manifest["domain_boundary"] = db[0]["name"]
    _validate_manifest(manifest, limits)
    mbytes = _dumps_manifest(manifest)
    # Round trip through the reader's own parser and checks (review M4).
    manifest = _loads_strict(mbytes)
    _validate_manifest(manifest, limits)
    _check_write_limits([len(mbytes)] + [len(b) for b in blobs], limits)
    data = _write_archive_bytes(
        [("manifest.json", mbytes)] + [(p["name"], b) for p, b in zip(parts, blobs)])
    _atomic_write(dst, data)
    return ChilReport(
        path=os.fspath(dst),
        format_version=0,
        crs=manifest["crs"],
        parts=tuple(PartInfo(p["name"], p["role"], p["kind"], p["size"], p["sha256"])
                    for p in parts),
        domain_boundary=manifest.get("domain_boundary"),
        summary=manifest.get("summary"),
        ext_namespaces=tuple(sorted(new_ext)),
        source_archive_sha256=archive_sha,
        archive_sha256=hashlib.sha256(data).hexdigest(),
    )
