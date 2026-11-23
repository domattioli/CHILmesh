"""Archive limit and unsupported-entry tests for .chil (T022, FR-058, FR-065).

Limits count bytes actually decompressed, never header-declared sizes.
"""
from __future__ import annotations

import functools
import importlib
import struct
import warnings
import zipfile

import pytest

from _chil_builders import build_archive, default_manifest, good_fort14_bytes, patch_zip_sizes

chil_io = importlib.import_module("chilmesh.chil_io")

MIB = 1 << 20


@functools.lru_cache(maxsize=1)
def _f14() -> bytes:
    return good_fort14_bytes()


def _limits(**kw):
    return getattr(chil_io, "ChilLimits")(**kw)


def _write(tmp_path, parts=None, name="x.chil", **kw) -> str:
    parts = {"fort.14": _f14()} if parts is None else parts
    p = tmp_path / name
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p.write_bytes(build_archive(parts, **kw))
    return str(p)


def _refuses(path, slug, **limit_kw):
    err = getattr(chil_io, "ChilFormatError")
    validate = getattr(chil_io, "validate_chil")
    with pytest.raises(err) as ei:
        validate(path, limits=_limits(**limit_kw))
    assert ei.value.check == slug
    return ei.value


def _set_flag_bits(archive: bytes, name: str, flag: int) -> bytes:
    """OR ``flag`` into the general-purpose bit flag of entry ``name`` (both headers)."""
    data = bytearray(archive)
    nb = name.encode()
    pos = 0
    while (pos := data.find(b"PK\x03\x04", pos)) >= 0:
        fn_len = struct.unpack_from("<H", data, pos + 26)[0]
        if bytes(data[pos + 30:pos + 30 + fn_len]) == nb:
            f = struct.unpack_from("<H", data, pos + 6)[0]
            struct.pack_into("<H", data, pos + 6, f | flag)
        pos += 4
    pos = 0
    while (pos := data.find(b"PK\x01\x02", pos)) >= 0:
        fn_len = struct.unpack_from("<H", data, pos + 28)[0]
        if bytes(data[pos + 46:pos + 46 + fn_len]) == nb:
            f = struct.unpack_from("<H", data, pos + 8)[0]
            struct.pack_into("<H", data, pos + 8, f | flag)
        pos += 4
    return bytes(data)


def _count_decompressed(monkeypatch):
    """Count bytes returned by ZipExtFile.read and read1."""
    total = [0]
    for attr in ("read", "read1"):
        orig = getattr(zipfile.ZipExtFile, attr)

        def wrap(self, *a, _orig=orig, **k):
            out = _orig(self, *a, **k)
            total[0] += len(out)
            return out

        monkeypatch.setattr(zipfile.ZipExtFile, attr, wrap)
    return total


def test_member_count_limit(tmp_path):
    parts = {"fort.14": _f14(), "fort.13": b"a\n"}
    _refuses(_write(tmp_path, parts), "limit-members", max_members=2)


def test_member_count_at_limit_is_accepted(tmp_path):
    parts = {"fort.14": _f14(), "fort.13": b"a\n"}
    getattr(chil_io, "validate_chil")(_write(tmp_path, parts), limits=_limits(max_members=3))


def test_per_member_byte_limit(tmp_path):
    _refuses(_write(tmp_path), "limit-member-bytes", max_member_bytes=100)


def test_total_byte_limit(tmp_path):
    parts = {"fort.14": _f14(), "fort.13": b"x" * 4000}
    n = len(_f14())
    _refuses(_write(tmp_path, parts), "limit-total-bytes", max_total_bytes=n + 2000)


def test_manifest_byte_limit(tmp_path):
    _refuses(_write(tmp_path), "limit-manifest-bytes", max_manifest_bytes=10)


def test_manifest_is_excluded_from_per_member_limit(tmp_path):
    """A manifest longer than max_member_bytes is fine when fort.14 fits."""
    path = _write(tmp_path)
    getattr(chil_io, "validate_chil")(path, limits=_limits(max_member_bytes=len(_f14())))


def test_lying_small_headers_do_not_widen_a_limit(tmp_path, monkeypatch):
    """Headers that declare 10 bytes must not let a bigger member through.

    zipfile stops at the declared size, so the refusal may surface as the
    limit itself, a size mismatch or a zip error; it must never validate.
    """
    parts = {"fort.14": _f14(), "fort.13": b"z" * (3 * MIB)}
    raw = build_archive(parts)
    raw = patch_zip_sizes(raw, "fort.13", local=10, central=10)
    p = tmp_path / "lie.chil"
    p.write_bytes(raw)
    err = getattr(chil_io, "ChilFormatError")
    total = _count_decompressed(monkeypatch)
    with pytest.raises(err) as ei:
        getattr(chil_io, "validate_chil")(str(p), limits=_limits(max_member_bytes=MIB))
    assert ei.value.check in {"limit-member-bytes", "size-mismatch", "not-zip"}
    assert total[0] <= len(_f14()) + MIB + MIB + 64 * 1024, total[0]


def test_compressible_bomb_stops_within_one_chunk(tmp_path, monkeypatch):
    """50 MiB of zeros against a 1 MiB limit: stop within limit plus one chunk."""
    assert len(_f14()) < MIB
    parts = {"fort.14": _f14(), "fort.13": b"\x00" * (50 * MIB)}
    path = _write(tmp_path, parts)
    total = _count_decompressed(monkeypatch)
    _refuses(path, "limit-member-bytes", max_member_bytes=MIB)
    assert total[0] <= len(_f14()) + MIB + MIB + 64 * 1024, total[0]


def test_encrypted_entry_refused(tmp_path):
    raw = _set_flag_bits(build_archive({"fort.14": _f14()}), "fort.14", 0x1)
    p = tmp_path / "enc.chil"
    p.write_bytes(raw)
    _refuses(str(p), "encrypted")


@pytest.mark.parametrize("method", [zipfile.ZIP_BZIP2, zipfile.ZIP_LZMA],
                         ids=["bzip2", "lzma"])
def test_unsupported_compression_refused(tmp_path, method):
    path = _write(tmp_path, compress=method)
    _refuses(path, "compression")


def test_stored_entries_are_accepted(tmp_path):
    path = _write(tmp_path, compress=zipfile.ZIP_STORED)
    report = getattr(chil_io, "validate_chil")(path)
    assert report.parts[0].size == len(_f14())


def test_default_limits_accept_a_normal_archive(tmp_path):
    parts = {"fort.14": _f14(), "fort.13": b"a\n"}
    m = default_manifest(parts)
    getattr(chil_io, "validate_chil")(_write(tmp_path, parts, manifest=m))
