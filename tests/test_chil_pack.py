"""Failing tests for ``chil_io.pack`` (specs/010-chil-container-mvp T044 to T046).

``pack`` is fetched per test with ``getattr`` so a missing name fails that test
only. Refusals assert ``ChilFormatError`` and the exact ``.check`` slug.
"""
from __future__ import annotations

import builtins
import functools
import hashlib
import importlib
import json
import os
import pathlib
import zipfile

import pytest

from _chil_builders import good_fort14_bytes

chil_io = importlib.import_module("chilmesh.chil_io")


def _pack(*a, **k):
    return getattr(chil_io, "pack")(*a, **k)


@functools.lru_cache(maxsize=1)
def _f14() -> bytes:
    return good_fort14_bytes()


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _tiny(elem_row="1 3 1 2 3", n_nodes=3) -> bytes:
    lines = ["tiny", f"1 {n_nodes}"]
    lines += [f"{i + 1} {float(i)} 0.0 0.0" for i in range(n_nodes)]
    lines += [elem_row, "0", "0", "0", "0"]
    return ("\n".join(lines) + "\n").encode()


def _bnd(open_rows=("2", "1"), land_rows=None, nodes=3) -> bytes:
    lines = ["tiny", f"1 {nodes}"]
    lines += [f"{i + 1} {float(i)} 0.0 0.0" for i in range(nodes)]
    lines += ["1 3 1 2 3", "1", str(len(open_rows)), f"{len(open_rows)} 0", *open_rows]
    if land_rows is None:
        lines += ["0", "0"]
    else:
        lines += ["1", str(len(land_rows)), f"{len(land_rows)} 0", *land_rows]
    return ("\n".join(lines) + "\n").encode()


def _src(tmp_path, name, data: bytes):
    p = tmp_path / name
    p.write_bytes(data)
    return p


def _refuses(check, fn, *a, **k):
    with pytest.raises(chil_io.ChilFormatError) as ei:
        fn(*a, **k)
    assert ei.value.check == check, ei.value
    return ei.value


def _no_output(tmp_path, dst):
    assert not dst.exists()
    assert sorted(p.name for p in tmp_path.iterdir() if p != dst and p.name.startswith(".")
                  or p.suffix == ".tmp") == []


def _members(path):
    with zipfile.ZipFile(path) as zf:
        return {n: zf.read(n) for n in zf.namelist()}


F13_CRLF = b"grid\r\n3\r\n1\r\nattr\r\n"
F15_MIXED = b"line one\r\nline two\nline three\r\nlast no newline"
F14_CRLF = _tiny().replace(b"\n", b"\r\n")
BND = b"boundary\r\ndata\n\xff\xfe raw"


# ------------------------------------------------------------ T044 FR-011

def test_pack_members_equal_input_bytes_and_fixed_names(tmp_path):
    f14 = _src(tmp_path, "any_name.14", _f14())
    f13 = _src(tmp_path, "a.13", F13_CRLF)
    f15 = _src(tmp_path, "b.15", F15_MIXED)
    db = _src(tmp_path, "Coast.bin", BND)
    dst = tmp_path / "o.chil"
    _pack(dst, f14, fort13=f13, fort15=f15, domain_boundary=db)
    m = _members(dst)
    assert set(m) == {"manifest.json", "fort.14", "fort.13", "fort.15",
                      "domain_boundary/Coast.bin"}
    assert m["fort.14"] == _f14()
    assert m["fort.13"] == F13_CRLF
    assert m["fort.15"] == F15_MIXED
    assert m["domain_boundary/Coast.bin"] == BND


def test_pack_crlf_fort14_survives_stdlib_zip(tmp_path):
    dst = tmp_path / "o.chil"
    _pack(dst, _src(tmp_path, "x", F14_CRLF))
    assert _members(dst)["fort.14"] == F14_CRLF
    assert b"\r\n" in _members(dst)["fort.14"]


def test_pack_no_trailing_newline_and_mixed_endings_fort14(tmp_path):
    data = _tiny().rstrip(b"\n")
    assert not data.endswith(b"\n")
    dst = tmp_path / "o.chil"
    _pack(dst, _src(tmp_path, "x", data))
    assert _members(dst)["fort.14"] == data


def test_pack_read_member_returns_same_bytes(tmp_path):
    dst = tmp_path / "o.chil"
    _pack(dst, _src(tmp_path, "x", _f14()), fort13=_src(tmp_path, "a", F13_CRLF),
          fort15=_src(tmp_path, "b", F15_MIXED), domain_boundary=_src(tmp_path, "c.bin", BND))
    rm = getattr(chil_io, "read_member")
    assert rm(dst, "fort14") == _f14()
    assert rm(dst, "fort13") == F13_CRLF
    assert rm(dst, "fort15") == F15_MIXED
    assert rm(dst, "domain_boundary") == BND


def test_pack_domain_boundary_pointer_and_name(tmp_path):
    dst = tmp_path / "o.chil"
    rep = _pack(dst, _src(tmp_path, "x", _f14()), domain_boundary=_src(tmp_path, "c.bin", BND))
    assert rep.domain_boundary == "domain_boundary/c.bin"
    man = json.loads(_members(dst)["manifest.json"])
    assert man["domain_boundary"] == "domain_boundary/c.bin"


def test_pack_without_domain_boundary_has_no_pointer(tmp_path):
    dst = tmp_path / "o.chil"
    _pack(dst, _src(tmp_path, "x", _f14()))
    assert "domain_boundary" not in json.loads(_members(dst)["manifest.json"])


def test_pack_never_uses_text_readers(tmp_path, monkeypatch):
    def boom(*a, **k):
        raise AssertionError("text reader called")

    for modname in ("chilmesh.fort13_io", "chilmesh.fort15_io"):
        mod = importlib.import_module(modname)
        for n in dir(mod):
            if n.startswith("read_"):
                monkeypatch.setattr(mod, n, boom)
    real_open = builtins.open

    def guarded(file, mode="r", *a, **k):
        if "b" not in mode and str(file).endswith((".13", ".15")):
            raise AssertionError(f"text-mode open of {file}")
        return real_open(file, mode, *a, **k)

    monkeypatch.setattr(builtins, "open", guarded)
    dst = tmp_path / "o.chil"
    _pack(dst, _src(tmp_path, "x.14", _f14()), fort13=_src(tmp_path, "a.13", F13_CRLF),
          fort15=_src(tmp_path, "b.15", F15_MIXED))
    assert _members(dst)["fort.13"] == F13_CRLF


# ------------------------------------------------------- T045 FR-010 FR-013

def test_pack_parts_role_source_size_sha(tmp_path):
    dst = tmp_path / "o.chil"
    rep = _pack(dst, _src(tmp_path, "x", _f14()), fort13=_src(tmp_path, "a", F13_CRLF),
                fort15=_src(tmp_path, "b", F15_MIXED), domain_boundary=_src(tmp_path, "c.bin", BND))
    want = {"fort.14": _f14(), "fort.13": F13_CRLF, "fort.15": F15_MIXED,
            "domain_boundary/c.bin": BND}
    man = json.loads(_members(dst)["manifest.json"])
    assert {p["name"] for p in man["parts"]} == set(want)
    for p in man["parts"]:
        assert p["role"] == "source"
        assert p["size"] == len(want[p["name"]])
        assert p["sha256"] == _sha(want[p["name"]])
    assert {p.name: p.role for p in rep.parts} == {n: "source" for n in want}
    assert rep.archive_sha256 == _sha(dst.read_bytes())


def test_pack_default_crs_unknown_and_explicit(tmp_path):
    a, b = tmp_path / "a.chil", tmp_path / "b.chil"
    s = _src(tmp_path, "x", _f14())
    assert _pack(a, s).crs == "unknown"
    assert _pack(b, s, crs="EPSG:4326").crs == "EPSG:4326"
    assert json.loads(_members(a)["manifest.json"])["crs"] == "unknown"
    assert json.loads(_members(b)["manifest.json"])["crs"] == "EPSG:4326"


def test_pack_cached_summary_present(tmp_path):
    dst = tmp_path / "o.chil"
    rep = _pack(dst, _src(tmp_path, "x", _tiny()))
    s = json.loads(_members(dst)["manifest.json"])["summary"]
    for k in ("n_nodes", "n_elements", "n_triangles", "n_quads", "grid_name"):
        assert k in s
    assert (s["n_nodes"], s["n_elements"], s["n_triangles"], s["n_quads"]) == (3, 1, 1, 0)
    assert rep.summary == s
    assert getattr(chil_io, "validate_chil")(dst).summary == s


def test_pack_then_load_honors_flags(tmp_path):
    from chilmesh import CHILmesh

    dst = tmp_path / "o.chil"
    _pack(dst, _src(tmp_path, "x", _f14()))
    full = CHILmesh.load(dst)
    assert full.n_layers > 0 and full.adjacencies
    bare = CHILmesh.load(dst, compute_layers=False, compute_adjacencies=False)
    assert bare.n_layers == 0 and not bare.adjacencies
    assert bare.n_verts == full.n_verts


# ------------------------------------------------- T046 refusals FR-012 etc.

def test_pack_node_zero_in_element_row(tmp_path):
    dst = tmp_path / "o.chil"
    exc = _refuses("node-id-zero", _pack, dst, _src(tmp_path, "x", _tiny("1 3 0 2 3")))
    assert "0" in str(exc)
    _no_output(tmp_path, dst)


def test_pack_node_range_in_element_row(tmp_path):
    dst = tmp_path / "o.chil"
    _refuses("node-id-range", _pack, dst, _src(tmp_path, "x", _tiny("1 3 1 2 9")))
    _no_output(tmp_path, dst)


def test_pack_boundary_node_zero_open_and_land(tmp_path):
    for i, data in enumerate((_bnd(("1", "0")), _bnd(("1", "2"), ("3", "0")))):
        dst = tmp_path / f"o{i}.chil"
        _refuses("node-id-zero", _pack, dst, _src(tmp_path, f"x{i}", data))
        assert not dst.exists()


def test_pack_boundary_node_range(tmp_path):
    dst = tmp_path / "o.chil"
    _refuses("node-id-range", _pack, dst, _src(tmp_path, "x", _bnd(("1", "9"))))
    assert not dst.exists()


def test_pack_unparseable_boundary_block_not_refused(tmp_path):
    base = _tiny().decode().rstrip("\n").split("\n")[:-4]
    # Executive correction 2026-11-19: the bad id comes AFTER the unparseable row.
    data = ("\n".join(base + ["1", "2", "2 0", "garbage", "0"]) + "\n").encode()
    dst = tmp_path / "o.chil"
    _pack(dst, _src(tmp_path, "x", data))
    assert _members(dst)["fort.14"] == data


def test_pack_unparseable_fort14(tmp_path):
    dst = tmp_path / "o.chil"
    _refuses("fort14-parse", _pack, dst, _src(tmp_path, "x", b"not a mesh\nat all\n"))
    assert not dst.exists()


@pytest.mark.parametrize("name", ["a\\b", "x" * 250, "a\x01b"], ids=["backslash", "too_long", "control_char"])
def test_pack_unsafe_domain_boundary_basename(tmp_path, name):
    """FR-067: an unsafe basename gives member-name. "x" * 250 makes the member
    name domain_boundary/<250 bytes>, over the 255-byte limit."""
    dst = tmp_path / "o.chil"
    p = tmp_path / name
    try:
        p.write_bytes(BND)
    except OSError:
        pytest.skip("filesystem rejects this name")
    _refuses("member-name", _pack, dst, _src(tmp_path, "x", _f14()), domain_boundary=p)
    assert not dst.exists()


def test_pack_domain_boundary_directory_input_refused(tmp_path):
    """A directory path as the boundary input is not a regular file (request)."""
    d = tmp_path / "dir"
    d.mkdir()
    dst = tmp_path / "o.chil"
    _refuses("request", _pack, dst, _src(tmp_path, "x", _f14()), domain_boundary=d)
    assert not dst.exists()


def test_pack_refuses_ext_member_keyword(tmp_path):
    dst = tmp_path / "o.chil"
    with pytest.raises(TypeError):
        _pack(dst, _src(tmp_path, "x", _f14()), ext_members={"ext/a/b": b"1"})
    assert not dst.exists()


def test_pack_never_writes_ext_member(tmp_path):
    dst = tmp_path / "o.chil"
    _pack(dst, _src(tmp_path, "x", _f14()), fort13=_src(tmp_path, "a", F13_CRLF))
    assert not [n for n in _members(dst) if n.startswith("ext/")]


def test_pack_refuses_directory_input(tmp_path):
    d = tmp_path / "d"
    d.mkdir()
    dst = tmp_path / "o.chil"
    _refuses("request", _pack, dst, d)
    assert not dst.exists()


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="no mkfifo")
def test_pack_refuses_fifo_input(tmp_path):
    fifo = tmp_path / "f"
    os.mkfifo(fifo)
    dst = tmp_path / "o.chil"
    _refuses("request", _pack, dst, fifo)  # must not block on open
    assert not dst.exists()


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="no mkfifo")
def test_pack_refuses_fifo_companion(tmp_path):
    fifo = tmp_path / "f"
    os.mkfifo(fifo)
    dst = tmp_path / "o.chil"
    _refuses("request", _pack, dst, _src(tmp_path, "x", _f14()), fort13=fifo)
    assert not dst.exists()


def test_pack_missing_input_is_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        _pack(tmp_path / "o.chil", tmp_path / "nope")


def test_pack_input_changed_between_passes(tmp_path, monkeypatch):
    """Implementation-sensitive: appends to the input after its first read-open closes.

    Hooks builtins.open and pathlib.Path.open; an implementation using os.open
    directly will not trigger the change (and the test then fails: pack succeeds).
    """
    src = _src(tmp_path, "x.14", _f14())
    dst = tmp_path / "o.chil"
    state = {"closes": 0}
    real_open = builtins.open
    target = str(src)

    class Wrap:
        def __init__(self, fh):
            self._fh = fh

        def __getattr__(self, n):
            return getattr(self._fh, n)

        def __enter__(self):
            self._fh.__enter__()
            return self

        def __exit__(self, *a):
            r = self._fh.__exit__(*a)
            self._after()
            return r

        def close(self):
            self._fh.close()
            self._after()

        def _after(self):
            if state["closes"] == 0:
                state["closes"] = 1
                with real_open(target, "ab") as w:
                    w.write(b"\n")

        def __iter__(self):
            return iter(self._fh)

    def hooked(file, mode="r", *a, **k):
        fh = real_open(file, mode, *a, **k)
        if str(file) == target and "b" in mode and "r" in mode and "+" not in mode:
            return Wrap(fh)
        return fh

    monkeypatch.setattr(builtins, "open", hooked)
    monkeypatch.setattr(pathlib.Path, "open",
                        lambda self, mode="r", *a, **k: hooked(str(self), mode, *a, **k))
    _refuses("request", _pack, dst, src)
    assert not dst.exists()


def test_pack_input_grown_between_passes_reads_bounded(tmp_path, monkeypatch):
    """An input that grows 64 MiB after pass 1 is refused; pass 2 reads at most size + one chunk."""
    src = _src(tmp_path, "x.14", _f14())
    size1 = src.stat().st_size
    dst = tmp_path / "o.chil"
    state = {"opens": 0, "pass2_read": 0}
    real_open = builtins.open
    target = str(src)

    class Wrap:
        def __init__(self, fh, idx):
            self._fh = fh
            self._idx = idx

        def __getattr__(self, n):
            return getattr(self._fh, n)

        def read(self, *a):
            data = self._fh.read(*a)
            if self._idx == 2:
                state["pass2_read"] += len(data)
            return data

        def __enter__(self):
            self._fh.__enter__()
            return self

        def __exit__(self, *a):
            r = self._fh.__exit__(*a)
            if self._idx == 1:
                with real_open(target, "r+b") as w:
                    w.truncate(size1 + 64 * 1024 * 1024)
            return r

    def hooked(file, mode="r", *a, **k):
        fh = real_open(file, mode, *a, **k)
        if str(file) == target and "b" in mode and "r" in mode and "+" not in mode:
            state["opens"] += 1
            return Wrap(fh, state["opens"])
        return fh

    monkeypatch.setattr(builtins, "open", hooked)
    monkeypatch.setattr(pathlib.Path, "open",
                        lambda self, mode="r", *a, **k: hooked(str(self), mode, *a, **k))
    _refuses("request", _pack, dst, src)
    assert state["opens"] == 2
    assert state["pass2_read"] <= size1 + chil_io._CHUNK
    assert not dst.exists()
    assert [p.name for p in tmp_path.iterdir()] == ["x.14"]


def test_pack_failed_replace_leaves_no_output(tmp_path, monkeypatch):
    dst = tmp_path / "o.chil"

    def boom(*a, **k):
        raise OSError("replace failed")

    monkeypatch.setattr(os, "replace", boom)
    with pytest.raises(OSError):
        _pack(dst, _src(tmp_path, "x", _f14()))
    assert not dst.exists()
    assert not [p for p in tmp_path.iterdir() if p.name not in ("x",)]


def test_pack_failed_replace_keeps_existing_destination(tmp_path, monkeypatch):
    dst = tmp_path / "o.chil"
    dst.write_bytes(b"old")
    monkeypatch.setattr(os, "replace", lambda *a, **k: (_ for _ in ()).throw(OSError("x")))
    with pytest.raises(OSError):
        _pack(dst, _src(tmp_path, "x", _f14()))
    assert dst.read_bytes() == b"old"


@pytest.mark.parametrize("limits,slug", [
    ({"max_members": 1}, "limit-members"),
    ({"max_member_bytes": 10}, "limit-member-bytes"),
    ({"max_total_bytes": 10}, "limit-total-bytes"),
    ({"max_manifest_bytes": 10}, "limit-manifest-bytes"),
])
def test_pack_over_limit_refused_before_output(tmp_path, limits, slug):
    dst = tmp_path / "o.chil"
    src = _src(tmp_path, "x", _f14())
    before = sorted(p.name for p in tmp_path.iterdir())
    _refuses(slug, _pack, dst, src, limits=chil_io.ChilLimits(**limits))
    assert sorted(p.name for p in tmp_path.iterdir()) == before
