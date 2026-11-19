"""Failing tests for rewrite_chil (specs/010-chil-container-mvp T049, T050, T050a, T051)."""
from __future__ import annotations

import hashlib
import importlib
import json
import os
import stat
import zipfile

import pytest

import _chil_builders as B

F13 = b"fort13 line one\r\nline two\r\n"
F15 = b"fort15 control\n"
DB = b"boundary polygon bytes\x00\xff\n"
EXT_A = b"\x00\x01signature-a\xff"
EXT_B = b"signature-b"


def _m():
    return importlib.import_module("chilmesh.chil_io")


def _g(name):
    return getattr(_m(), name)


def _f14():
    return B.good_fort14_bytes()


def _sha(b):
    return hashlib.sha256(b).hexdigest()


def _write(tmp_path, data, name="in.chil"):
    p = tmp_path / name
    p.write_bytes(data)
    return str(p)


def _members(path):
    with zipfile.ZipFile(path) as zf:
        return {i.filename: zf.read(i) for i in zf.infolist()}


def _manifest(path):
    return json.loads(_members(path)["manifest.json"])


def _full_parts():
    return {"fort.14": _f14(), "fort.13": F13, "fort.15": F15, "domain_boundary/bnd.poly": DB}


def _full(tmp_path, ext=None, ext_member=None, name="in.chil"):
    parts = {**_full_parts(), **(ext_member or {})}
    m = B.default_manifest(parts, domain_boundary="domain_boundary/bnd.poly")
    tables = dict(ext or {})
    for n in (ext_member or {}):
        tables.setdefault(n.split("/")[1], {})
    if tables:
        m["ext"] = tables
    return _write(tmp_path, B.build_archive(parts, manifest=m), name)


def _generated(tmp_path, ext=None, ext_member=None):
    parts = {"fort.14": _f14()}
    m = B.default_manifest({**parts, **(ext_member or {})})
    m["parts"][0]["role"] = "generated"
    tables = dict(ext or {})
    for n in (ext_member or {}):
        tables.setdefault(n.split("/")[1], {})
    if tables:
        m["ext"] = tables
    return _write(tmp_path, B.build_archive(parts, manifest=m, ext_member=ext_member))


def _dst(tmp_path, name="out.chil"):
    return str(tmp_path / name)


def _refused(check, fn):
    with pytest.raises(_g("ChilFormatError")) as ei:
        fn()
    assert ei.value.check == check
    return ei.value


def _only(tmp_path, *names):
    assert sorted(os.listdir(tmp_path)) == sorted(names)


NESTED = {"obj": {"a": [1, 2.5, True, None, "héllo 世界"], "b": {"c": False}},
          "n": None, "i": 7, "f": 0.25, "s": "café", "l": [[1], {"z": 0}, []]}


# ---------------------------------------------------------------- T049

def test_ext_tables_identical_after_rewrite(tmp_path):
    """FR-040, SC-005: ext tables keep keys, types, nested values, unknown owners."""
    ext = {"valence": NESTED, "somebody_else": {"k": [1, {"x": None}], "u": "ü"}}
    src = _full(tmp_path, ext=ext)
    out = _dst(tmp_path)
    _g("rewrite_chil")(src, out)
    got = _manifest(out)["ext"]
    assert got == ext
    assert json.dumps(got["valence"], sort_keys=True) == json.dumps(NESTED, sort_keys=True)
    assert type(got["valence"]["obj"]["a"][0]) is int
    assert type(got["valence"]["f"]) is float
    assert got["valence"]["obj"]["a"][2] is True


def test_ext_key_order_irrelevant(tmp_path):
    """FR-040: a manifest with unsorted ext keys rewrites to the same parsed tables."""
    parts = {"fort.14": _f14()}
    m = B.default_manifest(parts)
    m["ext"] = {"zeta": {"b": 1, "a": 2}, "alpha": {"y": {"q": 1, "p": 2}}}
    raw = json.dumps(m, indent=1).encode()  # not sorted
    src = _write(tmp_path, B.build_archive(parts, manifest=raw))
    out = _dst(tmp_path)
    _g("rewrite_chil")(src, out)
    assert _manifest(out)["ext"] == m["ext"]


def test_untouched_members_byte_identical(tmp_path):
    """FR-040, FR-042, SC-005: fort14, fort13 CRLF, fort15, boundary, ext members byte for byte."""
    src = _full(tmp_path, ext_member={"ext/valence/sig.bin": EXT_A})
    out = _dst(tmp_path)
    _g("rewrite_chil")(src, out, crs="EPSG:4326")
    a, b = _members(src), _members(out)
    for n in ("fort.14", "fort.13", "fort.15", "domain_boundary/bnd.poly", "ext/valence/sig.bin"):
        assert a[n] == b[n]
    assert b["fort.13"] == F13 and b"\r\n" in b["fort.13"]
    assert _manifest(out)["crs"] == "EPSG:4326"


def test_source_archive_sha256_is_input_file_hash(tmp_path):
    """FR-040: source_archive_sha256 equals sha256 of the src file."""
    src = _full(tmp_path)
    out = _dst(tmp_path)
    rep = _g("rewrite_chil")(src, out)
    raw = open(src, "rb").read()
    assert _manifest(out)["source_archive_sha256"] == _sha(raw)
    assert rep.source_archive_sha256 == _sha(raw)
    assert rep.archive_sha256 == _sha(open(out, "rb").read())


def test_source_sha_names_immediate_input(tmp_path):
    """FR-040: rewrite of a rewrite names the first output, not the original."""
    src = _full(tmp_path)
    mid, out = _dst(tmp_path, "mid.chil"), _dst(tmp_path, "out.chil")
    _g("rewrite_chil")(src, mid)
    _g("rewrite_chil")(mid, out)
    assert _manifest(out)["source_archive_sha256"] == _sha(open(mid, "rb").read())
    assert _manifest(out)["source_archive_sha256"] != _sha(open(src, "rb").read())


def test_output_validates_and_reports(tmp_path):
    """FR-040: result passes validate_chil; report parts match."""
    src = _full(tmp_path, ext_member={"ext/valence/sig.bin": EXT_A})
    out = _dst(tmp_path)
    rep = _g("rewrite_chil")(src, out)
    again = _g("validate_chil")(out)
    assert [(p.name, p.role, p.kind, p.size, p.sha256) for p in rep.parts] == \
        [(p.name, p.role, p.kind, p.size, p.sha256) for p in again.parts]


def test_role_and_kind_of_untouched_parts_unchanged(tmp_path):
    """Contract: role, kind, size, sha256 of untouched parts unchanged."""
    src = _full(tmp_path, ext_member={"ext/valence/sig.bin": EXT_A})
    out = _dst(tmp_path)
    _g("rewrite_chil")(src, out, set_parts={"fort13": b"new13\n"})
    a = {p["name"]: p for p in _manifest(src)["parts"]}
    b = {p["name"]: p for p in _manifest(out)["parts"]}
    for n in ("fort.14", "fort.15", "domain_boundary/bnd.poly", "ext/valence/sig.bin"):
        assert a[n] == b[n]
    assert (b["fort.13"]["role"], b["fort.13"]["kind"]) == ("source", "fort13")
    assert b["ext/valence/sig.bin"]["role"] == "extension"
    assert b["ext/valence/sig.bin"]["kind"] == "ext_member"


def test_add_member(tmp_path):
    """US4: add fort13 and fort15 to a fort14-only archive."""
    src = _write(tmp_path, B.build_archive({"fort.14": _f14()}))
    out = _dst(tmp_path)
    _g("rewrite_chil")(src, out, set_parts={"fort13": F13, "fort15": F15})
    mem = _members(out)
    assert mem["fort.13"] == F13 and mem["fort.15"] == F15
    assert {p["kind"] for p in _manifest(out)["parts"]} == {"fort14", "fort13", "fort15"}
    _g("validate_chil")(out)


def test_add_member_from_path(tmp_path):
    """US4: a bare path source is accepted as shorthand for MemberSpec."""
    src = _write(tmp_path, B.build_archive({"fort.14": _f14()}))
    extra = tmp_path / "x13.bin"
    extra.write_bytes(F13)
    out = _dst(tmp_path)
    _g("rewrite_chil")(src, out, set_parts={"fort13": str(extra)})
    assert _members(out)["fort.13"] == F13


def test_replace_member(tmp_path):
    """US4: replace fort13; others untouched, size and digest updated."""
    src = _full(tmp_path)
    out = _dst(tmp_path)
    new = b"replaced\r\n"
    _g("rewrite_chil")(src, out, set_parts={"fort13": new})
    assert _members(out)["fort.13"] == new
    p = {x["name"]: x for x in _manifest(out)["parts"]}["fort.13"]
    assert p["size"] == len(new) and p["sha256"] == _sha(new)
    assert _members(out)["fort.15"] == F15


def test_remove_domain_boundary_clears_key(tmp_path):
    """FR-064: removing the Domain Boundary part clears the manifest key."""
    src = _full(tmp_path)
    assert "domain_boundary" in _manifest(src)
    out = _dst(tmp_path)
    _g("rewrite_chil")(src, out, remove_parts=["domain_boundary"])
    m = _manifest(out)
    assert "domain_boundary" not in m
    assert "domain_boundary/bnd.poly" not in _members(out)
    assert all(p["kind"] != "domain_boundary" for p in m["parts"])
    _g("validate_chil")(out)


def test_add_domain_boundary_sets_key(tmp_path):
    """FR-064: adding a boundary with MemberSpec name sets the key."""
    src = _write(tmp_path, B.build_archive({"fort.14": _f14()}))
    out = _dst(tmp_path)
    _g("rewrite_chil")(src, out, set_parts={"domain_boundary": _g("MemberSpec")(DB, name="b.poly")})
    assert _manifest(out)["domain_boundary"] == "domain_boundary/b.poly"
    assert _members(out)["domain_boundary/b.poly"] == DB


def test_fort14_replace_or_remove_refused(tmp_path):
    """FR-041: add, replace or remove of the fort14 part -> request; no output, no temp."""
    src = _write(tmp_path, B.archive_with_summary(_f14()))
    for i, kw in enumerate(({"set_parts": {"fort14": _f14()}},
                            {"set_parts": {"fort14": _g("MemberSpec")(_f14())}},
                            {"remove_parts": ["fort14"]})):
        out = _dst(tmp_path, f"o{i}.chil")
        _refused("request", lambda: _g("rewrite_chil")(src, out, **kw))
        _only(tmp_path, "in.chil")


def test_summary_kept_when_fort14_untouched(tmp_path):
    """Contract: summary kept unless fort14 replaced."""
    src = _write(tmp_path, B.archive_with_summary(_f14()))
    out = _dst(tmp_path)
    _g("rewrite_chil")(src, out, set_parts={"fort13": F13})
    assert _manifest(out)["summary"] == _manifest(src)["summary"]


def test_ext_add_replace_remove_namespace_leaves_others(tmp_path):
    """FR-041, SC-005: add, replace, remove one namespace; others identical."""
    base = {"valence": NESTED, "other": {"keep": [1, "x"]}}
    src = _full(tmp_path, ext=base)
    r = _g("rewrite_chil")
    o1, o2, o3 = (_dst(tmp_path, n) for n in ("o1.chil", "o2.chil", "o3.chil"))
    r(src, o1, ext={"fresh": {"a": 1}})
    assert _manifest(o1)["ext"] == {**base, "fresh": {"a": 1}}
    r(o1, o2, ext={"valence": {"replaced": True}})
    assert _manifest(o2)["ext"] == {**base, "fresh": {"a": 1}, "valence": {"replaced": True}}
    r(o2, o3, ext={"fresh": None})
    assert _manifest(o3)["ext"] == {"valence": {"replaced": True}, "other": base["other"]}
    for o in (o1, o2, o3):
        _g("validate_chil")(o)


# ---------------------------------------------------------------- T050

@pytest.mark.parametrize("ns", ["crs", "parts", "summary", "format_version", "ext",
                                "domain_boundary", "source_archive_sha256"])
def test_core_field_through_ext_refused(tmp_path, ns):
    """FR-041: a core name as ext namespace is refused (ext-namespace); no output."""
    src = _full(tmp_path)
    out = _dst(tmp_path)
    _refused("ext-namespace", lambda: _g("rewrite_chil")(src, out, ext={ns: {"x": 1}}))
    _only(tmp_path, "in.chil")


@pytest.mark.parametrize("bad", ["Bad", "1abc", "has-dash", "", "a" * 65, "with space"])
def test_bad_namespace_shape_refused(tmp_path, bad):
    """FR-041: namespace failing the shape rule is refused (ext-namespace)."""
    src = _full(tmp_path)
    out = _dst(tmp_path)
    _refused("ext-namespace", lambda: _g("rewrite_chil")(src, out, ext={bad: {"x": 1}}))
    _only(tmp_path, "in.chil")


@pytest.mark.parametrize("table", [[1], "s", 3, 1.5, True, [{"a": 1}]])
def test_non_object_table_refused(tmp_path, table):
    """FR-041: non-object table refused (ext-table); no output."""
    src = _full(tmp_path)
    out = _dst(tmp_path)
    _refused("ext-table", lambda: _g("rewrite_chil")(src, out, ext={"valence": table}))
    _only(tmp_path, "in.chil")


# ---------------------------------------------------------------- T050a

def test_untouched_ext_members_byte_identical(tmp_path):
    """FR-043: untouched ext members survive byte for byte, with role and kind."""
    src = _full(tmp_path, ext_member={"ext/valence/a.bin": EXT_A, "ext/other/b.bin": EXT_B})
    out = _dst(tmp_path)
    _g("rewrite_chil")(src, out, set_parts={"fort13": b"x\n"})
    mem = _members(out)
    assert mem["ext/valence/a.bin"] == EXT_A and mem["ext/other/b.bin"] == EXT_B
    assert set(_manifest(out)["ext"]) == {"valence", "other"}


def test_ext_member_add_replace_remove(tmp_path):
    """FR-043: add, replace, remove an ext member on request."""
    r = _g("rewrite_chil")
    src = _full(tmp_path, ext={"valence": {}})
    o1, o2, o3 = (_dst(tmp_path, n) for n in ("o1.chil", "o2.chil", "o3.chil"))
    r(src, o1, set_parts={"ext/valence/sig.bin": EXT_A})
    assert _members(o1)["ext/valence/sig.bin"] == EXT_A
    p = {x["name"]: x for x in _manifest(o1)["parts"]}["ext/valence/sig.bin"]
    assert (p["kind"], p["role"], p["size"], p["sha256"]) == ("ext_member", "extension",
                                                              len(EXT_A), _sha(EXT_A))
    r(o1, o2, set_parts={"ext/valence/sig.bin": EXT_B})
    assert _members(o2)["ext/valence/sig.bin"] == EXT_B
    r(o2, o3, remove_parts=["ext/valence/sig.bin"])
    assert "ext/valence/sig.bin" not in _members(o3)
    assert "valence" in _manifest(o3)["ext"]
    for o in (o1, o2, o3):
        _g("validate_chil")(o)


def test_ext_member_ops_on_generated_archive(tmp_path):
    """FR-043, FR-006: add, replace, remove ext members on a generated archive."""
    r = _g("rewrite_chil")
    src = _generated(tmp_path, ext={"valence": {"v": 1}})
    o1, o2, o3 = (_dst(tmp_path, n) for n in ("o1.chil", "o2.chil", "o3.chil"))
    r(src, o1, set_parts={"ext/valence/sig.bin": EXT_A})
    assert _members(o1)["ext/valence/sig.bin"] == EXT_A
    assert {p["name"]: p["role"] for p in _manifest(o1)["parts"]}["fort.14"] == "generated"
    r(o1, o2, set_parts={"ext/valence/sig.bin": EXT_B})
    assert _members(o2)["ext/valence/sig.bin"] == EXT_B
    r(o2, o3, remove_parts=["ext/valence/sig.bin"])
    assert "ext/valence/sig.bin" not in _members(o3)
    _g("validate_chil")(o3)


def test_ext_member_without_owner_table_refused(tmp_path):
    """FR-068: adding a member whose ext.<owner> table is absent -> ext-member-owner."""
    src = _full(tmp_path)
    out = _dst(tmp_path)
    _refused("ext-member-owner",
             lambda: _g("rewrite_chil")(src, out, set_parts={"ext/valence/s.bin": EXT_A}))
    _only(tmp_path, "in.chil")


def test_remove_table_with_members_refused(tmp_path):
    """FR-068: removing ext.<owner> while members remain -> ext-member-owner."""
    src = _full(tmp_path, ext_member={"ext/valence/s.bin": EXT_A})
    out = _dst(tmp_path)
    _refused("ext-member-owner", lambda: _g("rewrite_chil")(src, out, ext={"valence": None}))
    _only(tmp_path, "in.chil")


def test_remove_table_and_members_together_ok(tmp_path):
    """FR-068: removing the table and its members in one request succeeds."""
    src = _full(tmp_path, ext_member={"ext/valence/s.bin": EXT_A})
    out = _dst(tmp_path)
    _g("rewrite_chil")(src, out, ext={"valence": None}, remove_parts=["ext/valence/s.bin"])
    assert "valence" not in _manifest(out).get("ext", {})
    assert "ext/valence/s.bin" not in _members(out)
    _g("validate_chil")(out)


def test_add_table_and_member_together_ok(tmp_path):
    """FR-068: adding the table and the member in one request succeeds."""
    src = _full(tmp_path)
    out = _dst(tmp_path)
    _g("rewrite_chil")(src, out, ext={"valence": {"v": 1}},
                       set_parts={"ext/valence/s.bin": EXT_A})
    assert _manifest(out)["ext"]["valence"] == {"v": 1}
    assert _members(out)["ext/valence/s.bin"] == EXT_A
    _g("validate_chil")(out)


@pytest.mark.parametrize("name", ["ext/valence", "ext/valence/", "valence/s.bin", "ext/s.bin",
                                  "other/valence/s.bin", "ext/Bad-Owner/s.bin"])
def test_ext_member_name_outside_owner_dir_refused(tmp_path, name):
    """FR-068: a name outside ext/<owner>/<rest> -> ext-member-name."""
    src = _full(tmp_path, ext={"valence": {}})
    out = _dst(tmp_path)
    _refused("ext-member-name", lambda: _g("rewrite_chil")(src, out, set_parts={name: EXT_A}))
    _only(tmp_path, "in.chil")


def test_ext_members_count_toward_member_limit(tmp_path):
    """FR-043, FR-069: ext members count toward max_members; no output."""
    src = _full(tmp_path, ext={"valence": {}})
    out = _dst(tmp_path)
    add = {f"ext/valence/m{i}.bin": b"x" * 3 for i in range(8)}
    lim = _g("ChilLimits")(max_members=6)
    _g("validate_chil")(src, limits=lim)  # source fits
    _refused("limit-members", lambda: _g("rewrite_chil")(src, out, set_parts=add, limits=lim))
    _only(tmp_path, "in.chil")


def test_ext_member_counts_toward_member_bytes(tmp_path):
    """FR-043, FR-069: oversized ext member -> limit-member-bytes; no output."""
    f14 = _f14()
    src = _write(tmp_path, B.build_archive({"fort.14": f14}, ext={"valence": {}}))
    out = _dst(tmp_path)
    lim = _g("ChilLimits")(max_member_bytes=len(f14) + 100)
    _refused("limit-member-bytes", lambda: _g("rewrite_chil")(
        src, out, set_parts={"ext/valence/big.bin": b"z" * (len(f14) + 1000)}, limits=lim))
    _only(tmp_path, "in.chil")


# ---------------------------------------------------------------- T051

def _bad_digest_archive():
    parts = {"fort.14": _f14()}
    m = B.default_manifest(parts)
    m["parts"][0]["sha256"] = "0" * 64
    return B.build_archive(parts, manifest=m)


def test_invalid_digest_source_refused_no_output(tmp_path):
    """FR-040: digest-mismatch input refused before any output."""
    src = _write(tmp_path, _bad_digest_archive())
    _refused("digest-mismatch", lambda: _g("rewrite_chil")(src, _dst(tmp_path)))
    _only(tmp_path, "in.chil")


def test_prototype_source_refused_no_output(tmp_path):
    """FR-060: c1_version input refused (prototype-container) before any output."""
    parts = {"fort.14": _f14()}
    m = B.default_manifest(parts, c1_version=1)
    src = _write(tmp_path, B.build_archive(parts, manifest=m))
    _refused("prototype-container", lambda: _g("rewrite_chil")(src, _dst(tmp_path)))
    _only(tmp_path, "in.chil")


def test_not_zip_source_refused_no_output(tmp_path):
    """FR-040: non-zip input refused (not-zip) before any output."""
    src = _write(tmp_path, b"this is not a zip archive at all")
    _refused("not-zip", lambda: _g("rewrite_chil")(src, _dst(tmp_path)))
    _only(tmp_path, "in.chil")


def test_existing_dst_untouched_when_source_invalid(tmp_path):
    """FR-066: an invalid source leaves an existing dst byte-identical."""
    src = _write(tmp_path, _bad_digest_archive())
    dst = tmp_path / "out.chil"
    dst.write_bytes(b"precious")
    with pytest.raises(_g("ChilFormatError")):
        _g("rewrite_chil")(src, str(dst))
    assert dst.read_bytes() == b"precious"


def test_in_place_rewrite_succeeds(tmp_path):
    """Contract: dst == src allowed; result valid, source sha names the old bytes."""
    src = _full(tmp_path, ext_member={"ext/valence/s.bin": EXT_A})
    old = open(src, "rb").read()
    _g("rewrite_chil")(src, src, crs="EPSG:4326")
    assert _manifest(src)["source_archive_sha256"] == _sha(old)
    assert _manifest(src)["crs"] == "EPSG:4326"
    assert _members(src)["ext/valence/s.bin"] == EXT_A
    _g("validate_chil")(src)
    _only(tmp_path, "in.chil")


def test_in_place_failed_replace_keeps_old_file(tmp_path, monkeypatch):
    """FR-066: dst == src with a failing os.replace keeps old bytes, no temp file."""
    src = _full(tmp_path)
    old = open(src, "rb").read()

    def boom(*a, **k):
        raise OSError("replace failed")

    monkeypatch.setattr(os, "replace", boom)
    with pytest.raises(OSError):
        _g("rewrite_chil")(src, src, crs="EPSG:4326")
    assert open(src, "rb").read() == old
    _only(tmp_path, "in.chil")


def test_failed_replace_new_dst_no_output_no_temp(tmp_path, monkeypatch):
    """FR-066: failing os.replace to a new dst leaves neither dst nor temp."""
    src = _full(tmp_path)
    monkeypatch.setattr(os, "replace", lambda *a, **k: (_ for _ in ()).throw(OSError("x")))
    with pytest.raises(OSError):
        _g("rewrite_chil")(src, _dst(tmp_path))
    _only(tmp_path, "in.chil")


def test_add_companion_to_generated_refused(tmp_path):
    """FR-006, FR-064: companion on a generated archive -> part-count; no output."""
    src = _generated(tmp_path)
    for kind in ("fort13", "fort15"):
        out = _dst(tmp_path, f"{kind}.chil")
        _refused("part-count", lambda: _g("rewrite_chil")(src, out, set_parts={kind: F13}))
    out = _dst(tmp_path, "db.chil")
    _refused("part-count", lambda: _g("rewrite_chil")(
        src, out, set_parts={"domain_boundary": _g("MemberSpec")(DB, name="b.poly")}))
    _only(tmp_path, "in.chil")


def test_remove_fort14_refused(tmp_path):
    """FR-041: removing fort14 -> request; no output."""
    src = _full(tmp_path)
    _refused("request", lambda: _g("rewrite_chil")(src, _dst(tmp_path), remove_parts=["fort14"]))
    _only(tmp_path, "in.chil")


def test_fifo_input_refused(tmp_path):
    """FR-066: a FIFO source is refused (request) without blocking."""
    fifo = tmp_path / "fifo.chil"
    os.mkfifo(fifo)
    _refused("request", lambda: _g("rewrite_chil")(str(fifo), _dst(tmp_path)))
    assert not os.path.exists(_dst(tmp_path))


def test_directory_input_refused(tmp_path):
    """FR-066: a directory source is refused (request)."""
    d = tmp_path / "adir"
    d.mkdir()
    _refused("request", lambda: _g("rewrite_chil")(str(d), _dst(tmp_path)))
    assert not os.path.exists(_dst(tmp_path))


def test_fifo_and_directory_member_inputs_refused(tmp_path):
    """FR-066: a FIFO or directory given as a set_parts path is refused (request)."""
    src = _full(tmp_path)
    fifo = tmp_path / "f.fifo"
    os.mkfifo(fifo)
    d = tmp_path / "d"
    d.mkdir()
    for i, bad in enumerate((fifo, d)):
        out = _dst(tmp_path, f"o{i}.chil")
        _refused("request", lambda: _g("rewrite_chil")(src, out, set_parts={"fort13": str(bad)}))
        assert not os.path.exists(out)


def test_member_input_changed_between_passes_refused(tmp_path, monkeypatch):
    """FR-066: an input changing between hash pass and copy pass -> request; no output.

    Same `open` hook as the pack changed-input test: appends to the input after
    its first binary read-open closes.
    """
    import builtins
    import pathlib

    src = _full(tmp_path)
    extra = tmp_path / "x13.bin"
    extra.write_bytes(b"A" * 5000)
    state = {"closes": 0}
    real_open = builtins.open
    target = str(extra)

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
                    w.write(b"changed")

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
    out = _dst(tmp_path)
    _refused("request", lambda: _g("rewrite_chil")(src, out, set_parts={"fort13": target}))
    assert not os.path.exists(out)
    assert sorted(os.listdir(tmp_path)) == ["in.chil", "x13.bin"]


def test_input_over_limits_refused_before_output(tmp_path):
    """FR-058, FR-069: the same limits apply to the input archive; an input over them
    is refused with its limit-* check before any output exists."""
    f14 = _f14()
    src = _write(tmp_path, B.build_archive({"fort.14": f14}))
    out = _dst(tmp_path)
    for lim, check in ((_g("ChilLimits")(max_members=1), "limit-members"),
                       (_g("ChilLimits")(max_member_bytes=len(f14) - 1), "limit-member-bytes"),
                       (_g("ChilLimits")(max_total_bytes=len(f14) - 1), "limit-total-bytes")):
        _refused(check, lambda: _g("rewrite_chil")(src, out, limits=lim))
        _only(tmp_path, "in.chil")


def test_limit_members_refused_no_output(tmp_path):
    """FR-069: result over max_members -> limit-members; no output."""
    src = _write(tmp_path, B.build_archive({"fort.14": _f14()}))
    out = _dst(tmp_path)
    lim = _g("ChilLimits")(max_members=3)
    _g("validate_chil")(src, limits=lim)
    _refused("limit-members", lambda: _g("rewrite_chil")(
        src, out, set_parts={"fort13": F13, "fort15": F15,
                             "domain_boundary": _g("MemberSpec")(DB, name="b.poly")}, limits=lim))
    _only(tmp_path, "in.chil")


def test_limit_member_bytes_refused_no_output(tmp_path):
    """FR-069: result member over max_member_bytes -> limit-member-bytes; no output."""
    f14 = _f14()
    src = _write(tmp_path, B.build_archive({"fort.14": f14}))
    out = _dst(tmp_path)
    lim = _g("ChilLimits")(max_member_bytes=len(f14) + 10)
    _refused("limit-member-bytes", lambda: _g("rewrite_chil")(
        src, out, set_parts={"fort13": b"q" * (len(f14) + 500)}, limits=lim))
    _only(tmp_path, "in.chil")


def test_limit_total_bytes_refused_no_output(tmp_path):
    """FR-069: result over max_total_bytes -> limit-total-bytes; no output."""
    f14 = _f14()
    src = _write(tmp_path, B.build_archive({"fort.14": f14}))
    out = _dst(tmp_path)
    lim = _g("ChilLimits")(max_total_bytes=len(f14) + 4000)
    _g("validate_chil")(src, limits=lim)
    _refused("limit-total-bytes", lambda: _g("rewrite_chil")(
        src, out, set_parts={"fort13": b"q" * 3000, "fort15": b"r" * 3000}, limits=lim))
    _only(tmp_path, "in.chil")


def test_limit_manifest_bytes_refused_no_output(tmp_path):
    """FR-069: result manifest over max_manifest_bytes -> limit-manifest-bytes; no output."""
    src = _write(tmp_path, B.build_archive({"fort.14": _f14()}))
    msize = len(_members(src)["manifest.json"])
    lim = _g("ChilLimits")(max_manifest_bytes=msize + 200)
    _g("validate_chil")(src, limits=lim)
    out = _dst(tmp_path)
    _refused("limit-manifest-bytes", lambda: _g("rewrite_chil")(
        src, out, ext={"valence": {"blob": "x" * 2000}}, limits=lim))
    _only(tmp_path, "in.chil")
