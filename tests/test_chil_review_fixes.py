"""Independent review fixes for .chil (specs/010-chil-container-mvp T051b to T051d, T056a)."""
from __future__ import annotations

import doctest
import functools
import struct
import warnings

import pytest

import chilmesh.chil_io as chil_io
from _chil_builders import build_archive, default_manifest, good_fort14_bytes

Err = chil_io.ChilFormatError


@functools.lru_cache(maxsize=1)
def _f14() -> bytes:
    return good_fort14_bytes()


def _write(tmp_path, parts, manifest=None, name="x.chil", **kw) -> str:
    p = tmp_path / name
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p.write_bytes(build_archive(parts, manifest=manifest, **kw))
    return str(p)


def _custom(tmp_path, entries, **extra):
    """entries: list of (name, kind, role, bytes)."""
    from _chil_builders import _sha

    parts = {n: d for n, _k, _r, d in entries}
    m = {"format_version": 0, "crs": "unknown",
         "parts": [{"name": n, "role": r, "kind": k, "size": len(d), "sha256": _sha(d)}
                   for n, k, r, d in entries]}
    m.update(extra)
    return _write(tmp_path, parts, manifest=m)


# ------------------------------------------------------------ H1 / L3 (T051b)

def _shadow_archive(tmp_path):
    return _custom(tmp_path, [
        ("fort.14", "fort14", "source", _f14()),
        ("fort14", "fort13", "source", b"fort13 bytes\n1 0\n"),
    ])


def test_shadow_name_refused_on_validate(tmp_path):
    with pytest.raises(Err) as ei:
        chil_io.validate_chil(_shadow_archive(tmp_path))
    assert ei.value.check == "member-name"


def test_shadow_name_refused_on_load(tmp_path):
    with pytest.raises(Err) as ei:
        chil_io.load_chil(_shadow_archive(tmp_path))
    assert ei.value.check == "member-name"


def test_shadow_name_refused_on_read_member(tmp_path):
    with pytest.raises(Err) as ei:
        chil_io.read_member(_shadow_archive(tmp_path), "fort14")
    assert ei.value.check == "member-name"


@pytest.mark.parametrize("name,kind,role", [
    ("fort.13", "fort14", "source"),
    ("fort.15", "fort13", "source"),
    ("fort.13", "domain_boundary", "source"),
    ("domain_boundary/a/b", "domain_boundary", "source"),
    ("domain_boundary", "domain_boundary", "source"),
    ("other/x", "domain_boundary", "source"),
    ("ext/zz/sig", "fort13", "source"),
    ("ext/zz/sig", "domain_boundary", "source"),
])
def test_mis_bound_name_refused(tmp_path, name, kind, role):
    entries = [("fort.14", "fort14", "source", _f14())]
    if kind == "fort14":
        entries = []
    entries.append((name, kind, role, b"x\n"))
    if kind == "fort14":
        entries.append(("fort.14x", "fort13", "source", b"y"))  # keep other part
    extra = {"domain_boundary": name} if kind == "domain_boundary" else {}
    path = _custom(tmp_path, entries, **extra)
    with pytest.raises(Err) as ei:
        chil_io.validate_chil(path)
    assert ei.value.check in ("member-name", "part-count")
    if kind != "fort14":
        assert ei.value.check == "member-name"


def test_fort14_wrong_name_is_member_name(tmp_path):
    path = _custom(tmp_path, [("fort.13", "fort14", "source", _f14())])
    with pytest.raises(Err) as ei:
        chil_io.validate_chil(path)
    assert ei.value.check == "member-name"


def test_ext_member_outside_ext_still_refused(tmp_path):
    from _chil_builders import _sha

    d = b"x"
    path = _custom(tmp_path, [
        ("fort.14", "fort14", "source", _f14()),
        ("notext/zz/a", "ext_member", "extension", d),
    ], ext={"zz": {}})
    with pytest.raises(Err) as ei:
        chil_io.validate_chil(path)
    assert ei.value.check in ("ext-member-name", "member-name")


def test_good_bound_names_still_validate(tmp_path):
    path = _custom(tmp_path, [
        ("fort.14", "fort14", "source", _f14()),
        ("fort.13", "fort13", "source", b"a"),
        ("fort.15", "fort15", "source", b"b"),
        ("domain_boundary/x.shp", "domain_boundary", "source", b"c"),
    ], domain_boundary="domain_boundary/x.shp")
    assert len(chil_io.validate_chil(path).parts) == 4


# ------------------------------------------------------------------ M2 (T051c)

def _patch_central(data: bytes, index: int, flag_or: int = 0, name_byte=None) -> bytes:
    b = bytearray(data)
    pos = -1
    for _ in range(index + 1):
        pos = b.index(b"PK\x01\x02", pos + 1)
    flags = struct.unpack_from("<H", b, pos + 8)[0] | flag_or
    struct.pack_into("<H", b, pos + 8, flags)
    if name_byte is not None:
        b[pos + 46] = name_byte
    return bytes(b)


def _hostile(tmp_path, **kw) -> str:
    raw = build_archive({"fort.14": _f14()})
    p = tmp_path / "h.chil"
    p.write_bytes(_patch_central(raw, **kw))
    return str(p)


@pytest.mark.parametrize("fn", ["validate", "summary"])
def test_utf8_flag_bad_name_byte_is_not_zip(tmp_path, fn):
    path = _hostile(tmp_path, index=1, flag_or=0x800, name_byte=0xFF)
    call = chil_io.validate_chil if fn == "validate" else chil_io.summary_chil
    with pytest.raises(Err) as ei:
        call(path)
    assert ei.value.check == "not-zip"


@pytest.mark.parametrize("idx", [0, 1])
@pytest.mark.parametrize("fn", ["validate", "summary"])
def test_strong_encryption_flag_is_encrypted(tmp_path, fn, idx):
    path = _hostile(tmp_path, index=idx, flag_or=0x40)
    call = chil_io.validate_chil if fn == "validate" else chil_io.summary_chil
    with pytest.raises(Err) as ei:
        call(path)
    assert ei.value.check == "encrypted"


def test_hostile_headers_in_rewrite(tmp_path):
    path = _hostile(tmp_path, index=0, flag_or=0x800, name_byte=0xFF)
    with pytest.raises(Err) as ei:
        chil_io.rewrite_chil(path, tmp_path / "o.chil", ext={"a": {}})
    assert ei.value.check == "not-zip"
    assert not (tmp_path / "o.chil").exists()


# ------------------------------------------------------------------ M4 (T051d)

def _src(tmp_path):
    p = tmp_path / "src.chil"
    p.write_bytes(build_archive({"fort.14": _f14()}))
    return str(p)


def _nested(depth):
    d: dict = {}
    cur = d
    for _ in range(depth):
        cur["a"] = {}
        cur = cur["a"]
    return d


@pytest.mark.parametrize("label,kw", [
    ("deep", {"ext": {"zz": _nested(70)}}),
    ("surrogate", {"ext": {"zz": {"k": "\ud800"}}}),
    ("set_value", {"ext": {"zz": {"k": {1, 2}}}}),
    ("int_key", {"ext": {"zz": {1: "v"}}}),
    ("nested_int_key", {"ext": {"zz": {"a": [{2: "v"}]}}}),
    ("bytes_value", {"ext": {"zz": {"k": b"x"}}}),
    ("nan", {"ext": {"zz": {"k": float("nan")}}}),
    ("remove_list", {"remove_parts": [["a"]]}),
])
def test_rewrite_bad_caller_input_refused_no_output(tmp_path, label, kw):
    out = tmp_path / "out.chil"
    with pytest.raises(Err) as ei:
        chil_io.rewrite_chil(_src(tmp_path), out, **kw)
    assert ei.value.check in ("request", "non-finite", "manifest-json")
    if label not in ("nan",):
        assert ei.value.check in ("request", "manifest-json")
    assert not out.exists()
    assert not [p for p in tmp_path.iterdir() if p.name not in ("src.chil",)]


def test_rewrite_deep_but_legal_ext_round_trips(tmp_path):
    out = tmp_path / "out.chil"
    chil_io.rewrite_chil(_src(tmp_path), out, ext={"zz": _nested(10)})
    assert chil_io.validate_chil(str(out)).ext_namespaces == ("zz",)


# ------------------------------------------------------------ T056a doctests

def test_chil_io_doctests():
    res = doctest.testmod(chil_io)
    assert res.failed == 0
    assert res.attempted >= 12
