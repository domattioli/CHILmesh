"""Manifest, name-rule, strict-JSON, schema, public-type tests for .chil.

specs/010-chil-container-mvp tasks T007, T008, T009, T009c. Written before the
implementation: ``chilmesh.chil_io`` is an empty module, so every test here
fails with AttributeError (missing name) until Phase 2 lands.
"""
from __future__ import annotations

import dataclasses
import functools
import hashlib
import importlib
import json
import unicodedata

import pytest

from _chil_builders import build_archive, default_manifest, good_fort14_bytes

chil_io = importlib.import_module("chilmesh.chil_io")


@functools.lru_cache(maxsize=1)
def _f14() -> bytes:
    return good_fort14_bytes()


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _base_manifest(**extra) -> dict:
    return default_manifest({"fort.14": _f14()}, **extra)


def _archive(tmp_path, *, parts=None, manifest=None, **kw) -> str:
    parts = {"fort.14": _f14()} if parts is None else parts
    p = tmp_path / "x.chil"
    p.write_bytes(build_archive(parts, manifest=manifest, **kw))
    return str(p)


def _refuses(path, slug):
    err = getattr(chil_io, "ChilFormatError")
    with pytest.raises(err) as ei:
        getattr(chil_io, "validate_chil")(path)
    assert ei.value.check == slug
    return ei.value


def _manifest_refused(tmp_path, manifest, slug, parts=None):
    _refuses(_archive(tmp_path, parts=parts, manifest=manifest), slug)


# --------------------------------------------------------------------- T007

_BAD_NAMES = [
    ("dotdot", "../evil"),
    ("dotdot_mid", "a/../b"),
    ("absolute", "/abs"),
    ("backslash", "a\\b"),
    ("drive", "C:x"),
    ("nul", "a\x00b"),
    ("empty_segment", "a//b"),
    ("control_char", "a\x1fb"),
    ("del_char", "a\x7fb"),
    ("trailing_slash", "dir/"),
    ("dot_segment", "./a"),
    ("too_long", "a" * 256),
]


@pytest.mark.parametrize("label,name", _BAD_NAMES, ids=[b[0] for b in _BAD_NAMES])
def test_unsafe_member_name_refused(tmp_path, label, name):
    """Each unsafe name in the FR-051 table raises member-name."""
    parts = {"fort.14": _f14(), name: b"x"}
    path = _archive(tmp_path, parts=parts)
    _refuses(path, "member-name")


def test_name_of_255_bytes_is_not_a_name_rule_violation(tmp_path):
    """255 UTF-8 bytes is the longest allowed; the rule is 'over 255'."""
    name = "domain_boundary/" + "a" * (255 - len("domain_boundary/"))
    assert len(name.encode()) == 255
    parts = {"fort.14": _f14(), name: b"x"}
    manifest = default_manifest(parts, domain_boundary=name)
    report = getattr(chil_io, "validate_chil")(_archive(tmp_path, parts=parts, manifest=manifest))
    assert name in [p.name for p in report.parts]


def test_multibyte_name_over_255_bytes_refused(tmp_path):
    """The limit counts UTF-8 bytes, not characters."""
    name = "é" * 128  # 256 bytes, 128 characters
    parts = {"fort.14": _f14(), name: b"x"}
    _refuses(_archive(tmp_path, parts=parts), "member-name")


_PAIRS = [
    ("case_fold", "domain_boundary/Coast.txt", "domain_boundary/COAST.txt"),
    ("case_fold_companion", "fort.13", "FORT.13"),
    ("nfc", "domain_boundary/é.txt", "domain_boundary/é.txt"),
    ("nfc_plus_case", "domain_boundary/É.txt", "domain_boundary/é.txt"),
    ("sharp_s", "domain_boundary/straße", "domain_boundary/STRASSE"),
]


@pytest.mark.parametrize("label,a,b", _PAIRS, ids=[p[0] for p in _PAIRS])
def test_collision_pair_refused(tmp_path, label, a, b):
    """Names equal under NFC(casefold(NFC(x))) collide (duplicate-name)."""
    assert unicodedata.normalize("NFC", unicodedata.normalize("NFC", a).casefold()) == \
        unicodedata.normalize("NFC", unicodedata.normalize("NFC", b).casefold())
    parts = {"fort.14": _f14(), a: b"1", b: b"2"}
    _refuses(_archive(tmp_path, parts=parts), "duplicate-name")


def test_manifest_json_is_in_the_collision_set(tmp_path):
    """A part named MANIFEST.JSON collides with manifest.json."""
    parts = {"fort.14": _f14(), "MANIFEST.JSON": b"{}"}
    _refuses(_archive(tmp_path, parts=parts), "duplicate-name")


# --------------------------------------------------------------------- T008

def _json_bytes(obj_text: str) -> bytes:
    return obj_text.encode("utf-8")


def _good_text() -> str:
    return json.dumps(_base_manifest(), sort_keys=True, indent=2)


def test_duplicate_key_top_level(tmp_path):
    text = _good_text().rstrip("}") + ',\n  "crs": "other"\n}'
    _manifest_refused(tmp_path, _json_bytes(text), "duplicate-key")


def test_duplicate_key_nested(tmp_path):
    text = _good_text().replace('"kind": "fort14"', '"kind": "fort14", "kind": "fort14"', 1)
    assert text != _good_text()
    _manifest_refused(tmp_path, _json_bytes(text), "duplicate-key")


@pytest.mark.parametrize("literal", ["NaN", "Infinity", "-Infinity"])
def test_bare_nonfinite_constant(tmp_path, literal):
    text = _good_text().rstrip("}") + f',\n  "ext": {{"x": {{"v": {literal}}}}}\n}}'
    _manifest_refused(tmp_path, _json_bytes(text), "non-finite")


def test_overflowing_float_is_nonfinite(tmp_path):
    text = _good_text().rstrip("}") + ',\n  "ext": {"x": {"v": 1e999}}\n}'
    _manifest_refused(tmp_path, _json_bytes(text), "non-finite")


def test_deep_nesting_is_manifest_json(tmp_path):
    depth = 100_000
    text = _good_text().rstrip("}") + ',\n  "ext": {"x": {"v": ' + "[" * depth + "]" * depth + "}}\n}"
    _manifest_refused(tmp_path, _json_bytes(text), "manifest-json")


def test_utf8_bom_refused(tmp_path):
    _manifest_refused(tmp_path, b"\xef\xbb\xbf" + _json_bytes(_good_text()), "manifest-json")


def test_invalid_utf8_refused(tmp_path):
    _manifest_refused(tmp_path, b'{"crs": "\xff"}', "manifest-json")


def test_top_level_not_an_object(tmp_path):
    _manifest_refused(tmp_path, b"[]", "manifest-json")


@pytest.mark.parametrize("value", [True, False, 0.0, "0", None, 1, -1])
def test_format_version_must_be_integer_zero(tmp_path, value):
    """Bool, float, string, null and unknown versions are refused."""
    m = _base_manifest()
    m["format_version"] = value
    _manifest_refused(tmp_path, m, "format-version")


def test_format_version_missing(tmp_path):
    m = _base_manifest()
    del m["format_version"]
    _manifest_refused(tmp_path, m, "format-version")


def test_strict_json_accepts_valid_manifest(tmp_path):
    report = getattr(chil_io, "validate_chil")(_archive(tmp_path, manifest=_base_manifest()))
    assert report.format_version == 0


# --------------------------------------------------------------------- T009

def test_unknown_top_level_key(tmp_path):
    m = _base_manifest(surprise=1)
    _manifest_refused(tmp_path, m, "unknown-key")


def test_unknown_part_key(tmp_path):
    m = _base_manifest()
    m["parts"][0]["extra"] = 1
    _manifest_refused(tmp_path, m, "unknown-key")


def test_unknown_summary_key(tmp_path):
    m = _base_manifest(summary=_good_summary(extra=1))
    _manifest_refused(tmp_path, m, "unknown-key")


def _good_summary(**over) -> dict:
    s = {"n_nodes": 1, "n_elements": 1, "n_triangles": 1, "n_quads": 0,
         "grid_name": "g", "bbox": [0.0, 0.0, 1.0, 1.0]}
    s.update(over)
    return s


@pytest.mark.parametrize("sha", ["abc", "A" * 64, "g" * 64, "a" * 63, "a" * 65, 5, None])
def test_bad_sha256_shape(tmp_path, sha):
    m = _base_manifest()
    m["parts"][0]["sha256"] = sha
    _manifest_refused(tmp_path, m, "manifest-schema")


@pytest.mark.parametrize("size", [-1, True, False, 1.5, "10", None])
def test_bad_size(tmp_path, size):
    m = _base_manifest()
    m["parts"][0]["size"] = size
    _manifest_refused(tmp_path, m, "manifest-schema")


@pytest.mark.parametrize("field,value", [
    ("n_nodes", "3"), ("n_nodes", -1), ("n_nodes", True), ("n_elements", 1.0),
    ("grid_name", 5), ("bbox", [0, 0, 1]), ("bbox", [1, 0, 0, 1]),
    ("bbox", [0, 0, "a", 1]), ("bbox", "box"),
])
def test_bad_summary_field(tmp_path, field, value):
    m = _base_manifest(summary=_good_summary(**{field: value}))
    _manifest_refused(tmp_path, m, "manifest-schema")


def test_summary_missing_required_key(tmp_path):
    s = _good_summary()
    del s["n_quads"]
    _manifest_refused(tmp_path, _base_manifest(summary=s), "manifest-schema")


def test_summary_counts_inconsistent(tmp_path):
    s = _good_summary(n_elements=3, n_triangles=1, n_quads=1)
    _manifest_refused(tmp_path, _base_manifest(summary=s), "manifest-schema")


def test_summary_not_object(tmp_path):
    _manifest_refused(tmp_path, _base_manifest(summary=[1]), "manifest-schema")


@pytest.mark.parametrize("crs", ["", "a" * 257, "bad\x00", "tab\there", 5, None])
def test_bad_crs(tmp_path, crs):
    m = _base_manifest()
    m["crs"] = crs
    _manifest_refused(tmp_path, m, "manifest-schema")


def test_bad_source_archive_sha256(tmp_path):
    _manifest_refused(tmp_path, _base_manifest(source_archive_sha256="XYZ"), "manifest-schema")


@pytest.mark.parametrize("ns", ["Valence", "1abc", "a-b", "", "a" * 65,
                                "crs", "parts", "summary", "format_version", "ext",
                                "domain_boundary", "source_archive_sha256"])
def test_bad_ext_namespace(tmp_path, ns):
    _manifest_refused(tmp_path, _base_manifest(ext={ns: {}}), "ext-namespace")


@pytest.mark.parametrize("table", [[], "s", 1, None, True])
def test_ext_table_must_be_object(tmp_path, table):
    _manifest_refused(tmp_path, _base_manifest(ext={"valence": table}), "ext-table")


def test_ext_not_object(tmp_path):
    _manifest_refused(tmp_path, _base_manifest(ext=[]), "ext-table")


def test_valid_ext_namespace_accepted(tmp_path):
    m = _base_manifest(ext={"valence": {"a": [1, 2]}, "x_1": {}})
    report = getattr(chil_io, "validate_chil")(_archive(tmp_path, manifest=m))
    assert report.ext_namespaces == ("valence", "x_1")


def test_unknown_part_kind(tmp_path):
    m = _base_manifest()
    m["parts"][0]["kind"] = "fort99"
    _manifest_refused(tmp_path, m, "part-kind")


def test_unknown_part_role(tmp_path):
    m = _base_manifest()
    m["parts"][0]["role"] = "owner"
    _manifest_refused(tmp_path, m, "part-role")


def test_five_kinds_and_three_roles_accepted(tmp_path):
    """All kinds and roles in their legal pairings validate."""
    parts = {
        "fort.14": _f14(), "fort.13": b"a\n", "fort.15": b"b\n",
        "domain_boundary/c.txt": b"c", "ext/valence/sig.bin": b"d",
    }
    m = default_manifest(parts, domain_boundary="domain_boundary/c.txt", ext={"valence": {}})
    report = getattr(chil_io, "validate_chil")(_archive(tmp_path, parts=parts, manifest=m))
    assert {p.kind for p in report.parts} == {
        "fort14", "fort13", "fort15", "domain_boundary", "ext_member"}
    assert {p.role for p in report.parts} == {"source", "extension"}
    assert report.domain_boundary == "domain_boundary/c.txt"


@pytest.mark.parametrize("kind,role", [
    ("fort13", "generated"), ("fort13", "extension"), ("fort15", "generated"),
    ("domain_boundary", "generated"), ("domain_boundary", "extension"),
    ("fort14", "extension"),
])
def test_role_kind_pairing_core_kinds(tmp_path, kind, role):
    parts = {"fort.14": _f14(), "fort.13": b"a\n"}
    m = default_manifest(parts)
    m["parts"][1]["kind"] = kind
    m["parts"][1]["role"] = role
    if kind == "fort14":
        m["parts"][0]["kind"] = "fort13"
    _manifest_refused(tmp_path, m, "part-role", parts=parts)


def test_ext_member_with_source_role(tmp_path):
    parts = {"fort.14": _f14(), "ext/valence/a": b"1"}
    m = default_manifest(parts, ext={"valence": {}})
    m["parts"][1]["role"] = "source"
    _manifest_refused(tmp_path, m, "part-role", parts=parts)


def test_extension_role_on_non_ext_kind(tmp_path):
    parts = {"fort.14": _f14(), "fort.13": b"a"}
    m = default_manifest(parts)
    m["parts"][1]["role"] = "extension"
    _manifest_refused(tmp_path, m, "part-role", parts=parts)


def test_two_fort14_parts(tmp_path):
    parts = {"fort.14": _f14(), "other.14": _f14()}
    m = default_manifest(parts)
    m["parts"][1]["kind"] = "fort14"
    _manifest_refused(tmp_path, m, "part-count", parts=parts)


def test_zero_fort14_parts(tmp_path):
    parts = {"fort.13": b"a"}
    _manifest_refused(tmp_path, default_manifest(parts), "part-count", parts=parts)


def test_two_parts_of_one_companion_kind(tmp_path):
    parts = {"fort.14": _f14(), "fort.13": b"a", "second.13": b"b"}
    m = default_manifest(parts)
    m["parts"][2]["kind"] = "fort13"
    _manifest_refused(tmp_path, m, "part-count", parts=parts)


def test_generated_archive_with_companion(tmp_path):
    parts = {"fort.14": _f14(), "fort.13": b"a"}
    m = default_manifest(parts)
    m["parts"][0]["role"] = "generated"
    _manifest_refused(tmp_path, m, "part-count", parts=parts)


def test_generated_archive_with_ext_member_is_valid(tmp_path):
    parts = {"fort.14": _f14(), "ext/valence/a": b"1"}
    m = default_manifest(parts, ext={"valence": {}})
    m["parts"][0]["role"] = "generated"
    report = getattr(chil_io, "validate_chil")(_archive(tmp_path, parts=parts, manifest=m))
    assert report.parts[0].role == "generated"


def test_domain_boundary_part_without_key(tmp_path):
    parts = {"fort.14": _f14(), "domain_boundary/c.txt": b"c"}
    _manifest_refused(tmp_path, default_manifest(parts), "domain-boundary", parts=parts)


def test_domain_boundary_key_without_part(tmp_path):
    m = _base_manifest(domain_boundary="domain_boundary/c.txt")
    _manifest_refused(tmp_path, m, "domain-boundary")


def test_domain_boundary_key_names_wrong_kind(tmp_path):
    parts = {"fort.14": _f14(), "domain_boundary/c.txt": b"c", "fort.13": b"a"}
    m = default_manifest(parts, domain_boundary="fort.13")
    _manifest_refused(tmp_path, m, "domain-boundary", parts=parts)


def test_domain_boundary_key_names_unknown_part(tmp_path):
    parts = {"fort.14": _f14(), "domain_boundary/c.txt": b"c"}
    m = default_manifest(parts, domain_boundary="domain_boundary/other.txt")
    _manifest_refused(tmp_path, m, "domain-boundary", parts=parts)


def test_ext_member_name_outside_owner_dir(tmp_path):
    # Executive correction 2026-11-19: "ext//a" has an empty segment, so the
    # member-name safety check (run first, R5) refuses it as member-name.
    for name, slug in (("ext/x", "ext-member-name"), ("other/valence/a", "ext-member-name"),
                       ("ext//a", "member-name")):
        parts = {"fort.14": _f14(), name: b"1"}
        m = default_manifest(parts, ext={"valence": {}})
        m["parts"][1]["kind"] = "ext_member"
        m["parts"][1]["role"] = "extension"
        _manifest_refused(tmp_path, m, slug, parts=parts)


def test_ext_member_owner_table_missing(tmp_path):
    parts = {"fort.14": _f14(), "ext/valence/a": b"1"}
    _manifest_refused(tmp_path, default_manifest(parts), "ext-member-owner", parts=parts)


def test_c1_version_refused_first(tmp_path):
    """c1_version wins even when format_version is absent (prototype manifest)."""
    _manifest_refused(tmp_path, {"c1_version": 1, "parts": []}, "prototype-container")


def test_unlisted_and_missing_members(tmp_path):
    parts = {"fort.14": _f14()}
    _refuses(_archive(tmp_path, parts=parts, extra_entries=[("stray.txt", b"x")]),
             "unlisted-member")


# -------------------------------------------------------------------- T009c

def test_chil_format_error_contract():
    err = getattr(chil_io, "ChilFormatError")
    assert issubclass(err, ValueError)
    exc = err("digest-mismatch", "member fort.14 changed")
    assert exc.check == "digest-mismatch"
    assert str(exc) == "digest-mismatch: member fort.14 changed"
    with pytest.raises(ValueError):
        raise exc


def test_chil_limits_defaults_and_overrides():
    limits = getattr(chil_io, "ChilLimits")
    d = limits()
    assert (d.max_members, d.max_member_bytes, d.max_total_bytes, d.max_manifest_bytes) == (
        64, 2**31, 2**32, 2**24)
    o = limits(max_members=3, max_member_bytes=4, max_total_bytes=5, max_manifest_bytes=6)
    assert (o.max_members, o.max_member_bytes, o.max_total_bytes, o.max_manifest_bytes) == (
        3, 4, 5, 6)


def test_chil_limits_is_frozen():
    limits = getattr(chil_io, "ChilLimits")
    assert dataclasses.is_dataclass(limits)
    with pytest.raises(dataclasses.FrozenInstanceError):
        limits().max_members = 1


def test_part_info_fields():
    part = getattr(chil_io, "PartInfo")
    assert [f.name for f in dataclasses.fields(part)] == ["name", "role", "kind", "size", "sha256"]


def test_chil_report_fields():
    rep = getattr(chil_io, "ChilReport")
    assert [f.name for f in dataclasses.fields(rep)] == [
        "path", "format_version", "crs", "parts", "domain_boundary", "summary",
        "ext_namespaces", "source_archive_sha256", "archive_sha256"]


def test_report_parts_is_tuple_of_part_info(tmp_path):
    report = getattr(chil_io, "validate_chil")(_archive(tmp_path))
    assert isinstance(report.parts, tuple)
    assert all(isinstance(p, chil_io.PartInfo) for p in report.parts)
    assert report.archive_sha256 == _sha((tmp_path / "x.chil").read_bytes())


def _rewrite(tmp_path, parts, manifest, spec_map):
    src = _archive(tmp_path, parts=parts, manifest=manifest)
    return getattr(chil_io, "rewrite_chil")(src, str(tmp_path / "out.chil"), set_parts=spec_map)


def test_memberspec_domain_boundary_bytes_need_name(tmp_path):
    spec = getattr(chil_io, "MemberSpec")
    err = getattr(chil_io, "ChilFormatError")
    with pytest.raises(err) as ei:
        _rewrite(tmp_path, {"fort.14": _f14()}, None, {"domain_boundary": spec(b"data")})
    assert ei.value.check == "request"


def test_memberspec_domain_boundary_bytes_with_name(tmp_path):
    spec = getattr(chil_io, "MemberSpec")
    rep = _rewrite(tmp_path, {"fort.14": _f14()}, None,
                   {"domain_boundary": spec(b"data", name="coast.txt")})
    assert rep.domain_boundary == "domain_boundary/coast.txt"


def test_memberspec_domain_boundary_path_defaults_to_basename(tmp_path):
    spec = getattr(chil_io, "MemberSpec")
    f = tmp_path / "shore.csv"
    f.write_bytes(b"1,2\n")
    rep = _rewrite(tmp_path, {"fort.14": _f14()}, None, {"domain_boundary": spec(f)})
    assert rep.domain_boundary == "domain_boundary/shore.csv"


@pytest.mark.parametrize("kind", ["fort14", "fort13", "fort15"])
def test_memberspec_fixed_kinds_reject_a_name(tmp_path, kind):
    spec = getattr(chil_io, "MemberSpec")
    err = getattr(chil_io, "ChilFormatError")
    with pytest.raises(err) as ei:
        _rewrite(tmp_path, {"fort.14": _f14()}, None, {kind: spec(_f14(), name="custom")})
    assert ei.value.check == "request"


def test_memberspec_ext_member_name_none_or_equal(tmp_path):
    spec = getattr(chil_io, "MemberSpec")
    parts = {"fort.14": _f14(), "ext/valence/a": b"1"}
    m = default_manifest(parts, ext={"valence": {}})
    key = "ext/valence/b"
    for name in (None, key):
        rep = _rewrite(tmp_path, parts, m, {key: spec(b"2", name=name)})
        assert key in [p.name for p in rep.parts]
    err = getattr(chil_io, "ChilFormatError")
    with pytest.raises(err) as ei:
        _rewrite(tmp_path, parts, m, {key: spec(b"2", name="ext/valence/other")})
    assert ei.value.check == "request"


@pytest.mark.parametrize("bad", ["../x", "a\\b", "a/b"])
def test_memberspec_unsafe_domain_boundary_name(tmp_path, bad):
    spec = getattr(chil_io, "MemberSpec")
    err = getattr(chil_io, "ChilFormatError")
    with pytest.raises(err) as ei:
        _rewrite(tmp_path, {"fort.14": _f14()}, None,
                 {"domain_boundary": spec(b"d", name=bad)})
    assert ei.value.check == "member-name"


def test_memberspec_unsafe_ext_member_name(tmp_path):
    spec = getattr(chil_io, "MemberSpec")
    err = getattr(chil_io, "ChilFormatError")
    parts = {"fort.14": _f14(), "ext/valence/a": b"1"}
    m = default_manifest(parts, ext={"valence": {}})
    with pytest.raises(err) as ei:
        _rewrite(tmp_path, parts, m, {"ext/x": spec(b"2")})
    assert ei.value.check == "ext-member-name"
