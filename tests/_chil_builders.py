"""Builders shared by the .chil container tests (specs/010-chil-container-mvp)."""
from __future__ import annotations

import functools
import hashlib
import io
import json
import struct
import zipfile
from pathlib import Path

import numpy as np

GOLDEN_DIR = Path(__file__).parent / "_golden"


def fresh(name: str):
    """Fresh fixture mesh; bypasses the session cache that other tests may mutate."""
    from chilmesh import examples

    fn = getattr(examples, name)
    return getattr(fn, "__wrapped__", fn)()


def good_fort14_bytes() -> bytes:
    """Bytes of the annulus fixture as written by ``write_fort14``."""
    import tempfile

    from chilmesh.CHILmesh import write_fort14

    mesh = fresh("annulus")
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "fort.14"
        write_fort14(mesh, str(p))
        return p.read_bytes()


def mixed_mesh():
    """quad_2x2 plus one triangle, built through the constructor."""
    from chilmesh import CHILmesh

    base = fresh("quad_2x2")
    pts = np.vstack([np.asarray(base.points), [[3.0, 0.5, 0.0]]])
    conn = [list(map(int, row)) for row in base.connectivity_list]
    conn.append([2, 9, 5])
    return CHILmesh(connectivity=np.array([r + [r[0]] if len(r) == 3 else r for r in conn]),
                    points=pts, grid_name="mixed", compute_layers=False,
                    compute_adjacencies=False)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def default_manifest(parts: dict, **extra) -> dict:
    """Manifest dict for ``parts`` mapping member name to bytes (kind inferred)."""
    entries = []
    for name, data in parts.items():
        if name == "fort.14":
            kind, role = "fort14", "source"
        elif name == "fort.13":
            kind, role = "fort13", "source"
        elif name == "fort.15":
            kind, role = "fort15", "source"
        elif name.startswith("domain_boundary/"):
            kind, role = "domain_boundary", "source"
        elif name.startswith("ext/"):
            kind, role = "ext_member", "extension"
        else:
            kind, role = "fort13", "source"
        entries.append({"name": name, "role": role, "kind": kind,
                        "size": len(data), "sha256": _sha(data)})
    m = {"format_version": 0, "crs": "unknown", "parts": entries}
    m.update(extra)
    return m


def build_archive(parts, *, manifest=None, order=None, extra_entries=(), raw_names=None,
                  ext_member=None, ext=None, compress=zipfile.ZIP_DEFLATED) -> bytes:
    """Build a .chil-shaped zip in memory.

    parts: dict name -> bytes.
    manifest: dict, bytes (written verbatim) or None (derived from parts).
    order: list of member names giving write order (manifest.json may be placed
        anywhere); default manifest first, then parts in dict order.
    extra_entries: iterable of (name, bytes) written after, not in the manifest.
    raw_names: dict name -> raw zip filename to store (for unsafe names).
    ext_member: dict name -> bytes of ext members, added to parts with an
        ``ext.<owner>`` table (empty object) in the derived manifest.
    ext: dict namespace -> table placed in the derived manifest.
    """
    parts = dict(parts)
    if ext_member:
        parts.update(ext_member)
    if manifest is None:
        manifest = default_manifest(parts)
        tables = dict(ext or {})
        for name in (ext_member or {}):
            tables.setdefault(name.split("/")[1], {})
        if tables:
            manifest["ext"] = tables
    if isinstance(manifest, (bytes, bytearray)):
        mbytes = bytes(manifest)
    else:
        mbytes = (json.dumps(manifest, sort_keys=True, indent=2) + "\n").encode("utf-8")
    entries = {"manifest.json": mbytes, **parts}
    names = list(order) if order is not None else ["manifest.json", *parts]
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for name in names:
            zi = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            # Assign after construction: ZipInfo() truncates at NUL.
            zi.filename = (raw_names or {}).get(name, name)
            zi.compress_type = compress
            zi.external_attr = 0o644 << 16
            zf.writestr(zi, entries[name])
        for name, data in extra_entries:
            zi = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            zi.compress_type = compress
            zf.writestr(zi, data)
    return buf.getvalue()


def patch_zip_sizes(archive: bytes, name: str, *, local=None, central=None) -> bytes:
    """Overwrite the uncompressed size fields of entry ``name``.

    ``local`` and ``central`` give the size to write into the local file header
    and the central directory record (``None`` leaves that header alone). Used
    to build archives whose headers lie about decompressed size.
    """
    data = bytearray(archive)
    nb = name.encode("utf-8")
    pos = 0
    while True:
        pos = data.find(b"PK\x03\x04", pos)
        if pos < 0:
            break
        fn_len, ex_len = struct.unpack_from("<HH", data, pos + 26)
        if bytes(data[pos + 30:pos + 30 + fn_len]) == nb and local is not None:
            struct.pack_into("<I", data, pos + 24, local)
        pos += 4
    pos = 0
    while True:
        pos = data.find(b"PK\x01\x02", pos)
        if pos < 0:
            break
        fn_len = struct.unpack_from("<H", data, pos + 28)[0]
        if bytes(data[pos + 46:pos + 46 + fn_len]) == nb and central is not None:
            struct.pack_into("<I", data, pos + 24, central)
        pos += 4
    return bytes(data)


@functools.lru_cache(maxsize=8)
def _summary_fields_cached(f14: bytes) -> tuple:
    return tuple(_summary_fields_raw(f14).items())


def summary_fields(f14: bytes) -> dict:
    d = dict(_summary_fields_cached(f14))
    d["bbox"] = list(d["bbox"]) if "bbox" in d else None
    return d


def _summary_fields_raw(f14: bytes) -> dict:
    """FR-030 cached summary computed from fort.14 bytes (independent of chil_io)."""
    from chilmesh import CHILmesh

    lines = f14.decode("utf-8").splitlines()
    mesh = CHILmesh._from_fort14_lines(lines, compute_layers=False,
                                       compute_adjacencies=False)
    conn = np.asarray(mesh.connectivity_list)
    if conn.shape[1] == 3:
        n_tri = conn.shape[0]
    else:
        n_tri = int(np.sum(conn[:, 3] == conn[:, 0]))
    pts = np.asarray(mesh.points)[:, :2]
    return {
        "n_nodes": int(mesh.n_verts), "n_elements": int(mesh.n_elems),
        "n_triangles": n_tri, "n_quads": int(mesh.n_elems) - n_tri,
        "grid_name": lines[0].strip(),
        "bbox": (float(pts[:, 0].min()), float(pts[:, 1].min()),
                 float(pts[:, 0].max()), float(pts[:, 1].max())),
    }


def archive_with_summary(f14: bytes, *, summary="auto", overrides=None) -> bytes:
    """Archive of one fort.14 whose manifest carries a cached summary.

    summary: "auto" computes it, None omits it, or a dict used as given.
    overrides: fields replaced in the cached summary (to make it wrong).
    """
    parts = {"fort.14": f14}
    m = default_manifest(parts)
    if summary == "auto":
        summary = summary_fields(f14)
    if summary is not None:
        summary = {**summary, **(overrides or {})}
        m["summary"] = summary
    return build_archive(parts, manifest=m)
