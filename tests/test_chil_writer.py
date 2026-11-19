"""Deterministic zip writer and atomic write tests (T009d; R3, R7, FR-007).

Fixed internal names (see T014): ``_write_archive_bytes(entries)`` takes a
list of ``(name, bytes)`` pairs, manifest first, and returns the archive
bytes; ``_atomic_write(dst, data)`` writes ``data`` to ``dst`` atomically.
"""
from __future__ import annotations

import importlib
import io
import os
import time
import zipfile

import pytest

chil_io = importlib.import_module("chilmesh.chil_io")

_ENTRIES = [
    ("manifest.json", b'{"format_version": 0}\n'),
    ("fort.14", b"grid\n1 3\n" * 50),
    ("fort.13", b"companion\n" * 20),
    ("domain_boundary/c.txt", b"0,0\n1,1\n"),
]


def _write(entries=None) -> bytes:
    return getattr(chil_io, "_write_archive_bytes")(list(entries or _ENTRIES))


def test_member_order_is_order_given_manifest_first():
    names = zipfile.ZipFile(io.BytesIO(_write())).namelist()
    assert names == [n for n, _ in _ENTRIES]
    assert names[0] == "manifest.json"


def test_member_order_follows_input_for_other_orders():
    entries = [_ENTRIES[0], _ENTRIES[3], _ENTRIES[1]]
    assert zipfile.ZipFile(io.BytesIO(_write(entries))).namelist() == [n for n, _ in entries]


def test_members_roundtrip_through_stdlib_reader():
    zf = zipfile.ZipFile(io.BytesIO(_write()))
    assert {n: zf.read(n) for n in zf.namelist()} == dict(_ENTRIES)
    assert zf.testzip() is None


def test_same_entries_same_bytes_twice():
    assert _write() == _write()


def test_bytes_independent_of_clock_and_utime(monkeypatch, tmp_path):
    first = _write()
    real = time.time
    monkeypatch.setattr(time, "time", lambda: real() + 10_000_000)
    src = tmp_path / "f"
    src.write_bytes(b"x")
    os.utime(src, (1, 1))
    assert _write() == first


def test_bytes_independent_of_destination_directory(tmp_path):
    atomic = getattr(chil_io, "_atomic_write")
    data = _write()
    d1, d2 = tmp_path / "one", tmp_path / "two"
    d1.mkdir()
    d2.mkdir()
    atomic(str(d1 / "a.chil"), data)
    atomic(str(d2 / "b.chil"), data)
    assert (d1 / "a.chil").read_bytes() == (d2 / "b.chil").read_bytes() == data


def test_every_entry_has_fixed_metadata():
    for zi in zipfile.ZipFile(io.BytesIO(_write())).infolist():
        assert zi.date_time == (1980, 1, 1, 0, 0, 0), zi.filename
        assert zi.create_system == 3, zi.filename
        assert (zi.external_attr >> 16) & 0o777 == 0o644, zi.filename
        assert zi.compress_type == zipfile.ZIP_DEFLATED, zi.filename
        assert zi.comment == b""


def test_atomic_write_creates_file(tmp_path):
    dst = tmp_path / "out.chil"
    getattr(chil_io, "_atomic_write")(str(dst), b"payload")
    assert dst.read_bytes() == b"payload"
    assert sorted(os.listdir(tmp_path)) == ["out.chil"]


def test_atomic_write_replaces_existing(tmp_path):
    dst = tmp_path / "out.chil"
    dst.write_bytes(b"old")
    getattr(chil_io, "_atomic_write")(str(dst), b"new")
    assert dst.read_bytes() == b"new"


def test_failed_replace_keeps_old_file_and_leaves_no_temp(tmp_path, monkeypatch):
    dst = tmp_path / "out.chil"
    dst.write_bytes(b"old bytes")
    atomic = getattr(chil_io, "_atomic_write")

    def boom(*a, **k):
        raise OSError("replace failed")

    monkeypatch.setattr(os, "replace", boom)
    with pytest.raises(OSError, match="replace failed"):
        atomic(str(dst), b"new bytes")
    monkeypatch.undo()
    assert dst.read_bytes() == b"old bytes"
    assert os.listdir(tmp_path) == ["out.chil"]


def test_failed_replace_with_no_old_file_leaves_nothing(tmp_path, monkeypatch):
    dst = tmp_path / "out.chil"
    atomic = getattr(chil_io, "_atomic_write")
    monkeypatch.setattr(os, "replace", lambda *a, **k: (_ for _ in ()).throw(OSError("x")))
    with pytest.raises(OSError):
        atomic(str(dst), b"data")
    monkeypatch.undo()
    assert os.listdir(tmp_path) == []
