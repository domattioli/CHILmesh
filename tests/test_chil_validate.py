"""Validator refusal tests for .chil (specs/010-chil-container-mvp T015 to T024b).

Every refusal asserts ``ChilFormatError`` and the exact ``.check`` slug of
data-model.md. Written before ``validate_chil`` exists: names are fetched with
``getattr`` on the module so a missing name fails that test only.
"""
from __future__ import annotations

import functools
import hashlib
import importlib
import warnings

import pytest

from _chil_builders import build_archive, default_manifest, good_fort14_bytes

chil_io = importlib.import_module("chilmesh.chil_io")


@functools.lru_cache(maxsize=1)
def _f14() -> bytes:
    return good_fort14_bytes()


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _tiny_fort14(elem_row: str = "1 3 1 2 3", n_nodes: int = 3) -> bytes:
    """Smallest fort.14 whose single element row is ``elem_row``."""
    lines = ["tiny", f"1 {n_nodes}"]
    lines += [f"{i + 1} {float(i)} 0.0 0.0" for i in range(n_nodes)]
    lines += [elem_row, "0", "0", "0", "0"]
    return ("\n".join(lines) + "\n").encode()


def _write(tmp_path, *, parts=None, manifest=None, name="x.chil", **kw) -> str:
    parts = {"fort.14": _f14()} if parts is None else parts
    p = tmp_path / name
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # duplicate zip names warn
        p.write_bytes(build_archive(parts, manifest=manifest, **kw))
    return str(p)


def _refuses(path, slug, limits=None):
    err = getattr(chil_io, "ChilFormatError")
    validate = getattr(chil_io, "validate_chil")
    with pytest.raises(err) as ei:
        validate(path) if limits is None else validate(path, limits=limits)
    assert ei.value.check == slug
    return ei.value


# --------------------------------------------------------------- T015 FR-051

@pytest.mark.parametrize("label,name", [
    ("absolute", "/abs/fort.13"),
    ("dotdot", "../fort.13"),
    ("dotdot_inner", "a/../fort.13"),
    ("backslash", "a\\fort.13"),
    ("drive_letter", "C:fort.13"),
    ("nul_byte", "for\x00t.13"),
    ("control_char", "for\x1ft.13"),
])
def test_unsafe_name_is_member_name(tmp_path, label, name):
    path = _write(tmp_path, parts={"fort.14": _f14(), name: b"x"})
    _refuses(path, "member-name")


def test_directory_entry_is_member_name(tmp_path):
    path = _write(tmp_path, extra_entries=[("sub/", b"")])
    _refuses(path, "member-name")


# --------------------------------------------------------------- T016 FR-052

def test_exact_duplicate_name(tmp_path):
    path = _write(tmp_path, extra_entries=[("fort.14", _f14())])
    _refuses(path, "duplicate-name")


def test_case_fold_collision(tmp_path):
    parts = {"fort.14": _f14(), "domain_boundary/Coast.txt": b"1",
             "domain_boundary/COAST.txt": b"2"}
    _refuses(_write(tmp_path, parts=parts), "duplicate-name")


def test_nfc_collision(tmp_path):
    parts = {"fort.14": _f14(), "domain_boundary/é.txt": b"1",
             "domain_boundary/é.txt": b"2"}
    _refuses(_write(tmp_path, parts=parts), "duplicate-name")


# --------------------------------------------------------------- T017 FR-053

def test_member_not_in_parts_is_unlisted(tmp_path):
    path = _write(tmp_path, extra_entries=[("stray.txt", b"x")])
    _refuses(path, "unlisted-member")


def test_part_without_member_is_missing(tmp_path):
    m = default_manifest({"fort.14": _f14(), "fort.13": b"a\n"})
    path = _write(tmp_path, parts={"fort.14": _f14()}, manifest=m)
    _refuses(path, "missing-member")


# --------------------------------------------------------------- T018 FR-054

def test_missing_manifest(tmp_path):
    path = _write(tmp_path, order=["fort.14"])
    _refuses(path, "manifest-missing")


def test_manifest_not_first(tmp_path):
    path = _write(tmp_path, order=["fort.14", "manifest.json"])
    _refuses(path, "manifest-not-first")


def test_manifest_invalid_json(tmp_path):
    _refuses(_write(tmp_path, manifest=b"{not json"), "manifest-json")


def test_manifest_duplicate_key(tmp_path):
    text = b'{"format_version": 0, "format_version": 0, "crs": "unknown", "parts": []}'
    _refuses(_write(tmp_path, manifest=text), "duplicate-key")


def test_manifest_non_finite_number(tmp_path):
    text = b'{"format_version": 0, "ext": {"x": {"v": NaN}}}'
    _refuses(_write(tmp_path, manifest=text), "non-finite")


def test_not_a_zip(tmp_path):
    p = tmp_path / "x.chil"
    p.write_bytes(b"this is not a zip archive at all" * 10)
    _refuses(str(p), "not-zip")


# --------------------------------------------------------------- T019 FR-055

def _version_manifest(value=..., drop=False):
    m = default_manifest({"fort.14": _f14()})
    if drop:
        del m["format_version"]
    else:
        m["format_version"] = value
    return m


def test_format_version_missing(tmp_path):
    _refuses(_write(tmp_path, manifest=_version_manifest(drop=True)), "format-version")


@pytest.mark.parametrize("value", ["0", 0.0, True, 1, 7, None])
def test_format_version_wrong(tmp_path, value):
    """String, float, bool, null and unknown versions are all format-version."""
    _refuses(_write(tmp_path, manifest=_version_manifest(value)), "format-version")


# --------------------------------------------------------------- T020 FR-056

def test_size_mismatch(tmp_path):
    m = default_manifest({"fort.14": _f14()})
    m["parts"][0]["size"] += 1
    _refuses(_write(tmp_path, manifest=m), "size-mismatch")


def test_digest_mismatch(tmp_path):
    m = default_manifest({"fort.14": _f14()})
    m["parts"][0]["sha256"] = "0" * 64
    _refuses(_write(tmp_path, manifest=m), "digest-mismatch")


# --------------------------------------------------------------- T021 FR-057

def test_node_id_zero_names_the_element(tmp_path):
    data = _tiny_fort14("7 3 0 2 3")
    path = _write(tmp_path, parts={"fort.14": data})
    exc = _refuses(path, "node-id-zero")
    assert "7" in str(exc)


@pytest.mark.parametrize("row", ["1 3 1 2 4", "1 3 1 2 -2"])
def test_node_id_out_of_range(tmp_path, row):
    path = _write(tmp_path, parts={"fort.14": _tiny_fort14(row)})
    _refuses(path, "node-id-range")


# ------------------------------------------------------------- T021a FR-057

def _fort14_with_boundary(open_rows=("2", "1"), land_rows=None, *, nodes=3) -> bytes:
    """Tiny fort.14 with one open segment and an optional land segment.

    ``open_rows`` is the node id list of the open segment; ``land_rows`` the
    node id list of one land segment (``None`` for NBOU 0).
    """
    lines = ["tiny", f"1 {nodes}"]
    lines += [f"{i + 1} {float(i)} 0.0 0.0" for i in range(nodes)]
    lines += ["1 3 1 2 3"]
    lines += ["1", str(len(open_rows)), f"{len(open_rows)} 0", *open_rows]
    if land_rows is None:
        lines += ["0", "0"]
    else:
        lines += ["1", str(len(land_rows)), f"{len(land_rows)} 0", *land_rows]
    return ("\n".join(lines) + "\n").encode()


def test_boundary_scan_accepts_valid_ids(tmp_path):
    data = _fort14_with_boundary(("1", "2"), ("2", "3"))
    getattr(chil_io, "validate_chil")(_write(tmp_path, parts={"fort.14": data}))


def test_open_boundary_node_zero(tmp_path):
    data = _fort14_with_boundary(("1", "0"))
    exc = _refuses(_write(tmp_path, parts={"fort.14": data}), "node-id-zero")
    assert "open" in str(exc).lower() and "boundary" in str(exc).lower()


def test_land_boundary_node_zero(tmp_path):
    data = _fort14_with_boundary(("1", "2"), ("3", "0"))
    exc = _refuses(_write(tmp_path, parts={"fort.14": data}), "node-id-zero")
    assert "land" in str(exc).lower() or "flow" in str(exc).lower()


def test_boundary_node_id_above_n_nodes(tmp_path):
    data = _fort14_with_boundary(("1", "9"))
    exc = _refuses(_write(tmp_path, parts={"fort.14": data}), "node-id-range")
    assert "boundary" in str(exc).lower()


def test_unparseable_boundary_block_passes_the_scan(tmp_path):
    base = _tiny_fort14().decode().rstrip("\n").split("\n")[:-4]
    # Executive correction 2026-11-19: the bad id comes AFTER the unparseable row.
    data = ("\n".join(base + ["1", "2", "2 0", "garbage", "0"]) + "\n").encode()
    report = getattr(chil_io, "validate_chil")(_write(tmp_path, parts={"fort.14": data}))
    assert report.parts[0].kind == "fort14"


# ----------------------------------------------------------- T023 FR-059..061

def _m(parts, **extra):
    return default_manifest(parts, **extra)


def test_prototype_container_c1_version(tmp_path):
    m = {"c1_version": 1, "parts": []}
    exc = _refuses(_write(tmp_path, manifest=m), "prototype-container")
    assert "prototype" in str(exc).lower()


def test_unknown_role(tmp_path):
    m = _m({"fort.14": _f14()})
    m["parts"][0]["role"] = "owner"
    _refuses(_write(tmp_path, manifest=m), "part-role")


def test_unknown_kind(tmp_path):
    m = _m({"fort.14": _f14()})
    m["parts"][0]["kind"] = "fort99"
    _refuses(_write(tmp_path, manifest=m), "part-kind")


def test_two_fort14_parts(tmp_path):
    parts = {"fort.14": _f14(), "other.14": _f14()}
    m = _m(parts)
    m["parts"][1]["kind"] = "fort14"
    _refuses(_write(tmp_path, parts=parts, manifest=m), "part-count")


def test_zero_fort14_parts(tmp_path):
    parts = {"fort.13": b"a\n"}
    _refuses(_write(tmp_path, parts=parts, manifest=_m(parts)), "part-count")


def test_generated_archive_with_companion(tmp_path):
    parts = {"fort.14": _f14(), "fort.13": b"a\n"}
    m = _m(parts)
    m["parts"][0]["role"] = "generated"
    _refuses(_write(tmp_path, parts=parts, manifest=m), "part-count")


def test_two_parts_of_same_companion_kind(tmp_path):
    parts = {"fort.14": _f14(), "fort.13": b"a\n", "second.13": b"b\n"}
    m = _m(parts)
    m["parts"][2]["kind"] = "fort13"
    _refuses(_write(tmp_path, parts=parts, manifest=m), "part-count")


def test_companion_with_generated_role(tmp_path):
    parts = {"fort.14": _f14(), "fort.13": b"a\n"}
    m = _m(parts)
    m["parts"][1]["role"] = "generated"
    _refuses(_write(tmp_path, parts=parts, manifest=m), "part-role")


def test_unknown_top_level_key(tmp_path):
    m = _m({"fort.14": _f14()}, surprise=1)
    _refuses(_write(tmp_path, manifest=m), "unknown-key")


def test_domain_boundary_part_without_key(tmp_path):
    parts = {"fort.14": _f14(), "domain_boundary/c.txt": b"c"}
    _refuses(_write(tmp_path, parts=parts, manifest=_m(parts)), "domain-boundary")


def test_domain_boundary_key_without_part(tmp_path):
    m = _m({"fort.14": _f14()}, domain_boundary="domain_boundary/c.txt")
    _refuses(_write(tmp_path, manifest=m), "domain-boundary")


def test_domain_boundary_key_names_other_kind(tmp_path):
    parts = {"fort.14": _f14(), "domain_boundary/c.txt": b"c", "fort.13": b"a\n"}
    m = _m(parts, domain_boundary="fort.13")
    _refuses(_write(tmp_path, parts=parts, manifest=m), "domain-boundary")


# ---------------------------------------------------------- T023a FR-068, 043

def _ext_parts(**members):
    parts = {"fort.14": _f14()}
    parts.update(members)
    return parts


def test_valid_ext_member_passes_and_is_reported(tmp_path):
    parts = _ext_parts(**{"ext/valence/signature.sig": b"sig"})
    path = _write(tmp_path, parts=parts, manifest=_m(parts, ext={"valence": {}}))
    report = getattr(chil_io, "validate_chil")(path)
    part = {p.name: p for p in report.parts}["ext/valence/signature.sig"]
    assert (part.kind, part.role, part.size) == ("ext_member", "extension", 3)
    assert "valence" in report.ext_namespaces


def test_ext_member_without_owner_table(tmp_path):
    parts = _ext_parts(**{"ext/valence/a": b"1"})
    _refuses(_write(tmp_path, parts=parts, manifest=_m(parts)), "ext-member-owner")


@pytest.mark.parametrize("name", ["ext/x", "other/valence/a"])
def test_ext_member_name_outside_owner_dir(tmp_path, name):
    parts = _ext_parts(**{name: b"1"})
    m = _m(parts, ext={"valence": {}})
    m["parts"][1]["kind"] = "ext_member"
    m["parts"][1]["role"] = "extension"
    _refuses(_write(tmp_path, parts=parts, manifest=m), "ext-member-name")


def test_extension_role_with_fort13_kind(tmp_path):
    parts = {"fort.14": _f14(), "fort.13": b"a\n"}
    m = _m(parts)
    m["parts"][1]["role"] = "extension"
    _refuses(_write(tmp_path, parts=parts, manifest=m), "part-role")


def test_ext_member_kind_with_source_role(tmp_path):
    parts = _ext_parts(**{"ext/valence/a": b"1"})
    m = _m(parts, ext={"valence": {}})
    m["parts"][1]["role"] = "source"
    _refuses(_write(tmp_path, parts=parts, manifest=m), "part-role")


def test_generated_archive_with_ext_member_passes(tmp_path):
    parts = _ext_parts(**{"ext/valence/a": b"1"})
    m = _m(parts, ext={"valence": {}})
    m["parts"][0]["role"] = "generated"
    report = getattr(chil_io, "validate_chil")(_write(tmp_path, parts=parts, manifest=m))
    assert report.parts[0].role == "generated"


def test_generated_archive_with_companion_still_refused(tmp_path):
    parts = {"fort.14": _f14(), "fort.13": b"a\n", "ext/valence/a": b"1"}
    m = _m(parts, ext={"valence": {}})
    m["parts"][0]["role"] = "generated"
    _refuses(_write(tmp_path, parts=parts, manifest=m), "part-count")


def test_unlisted_entry_under_ext(tmp_path):
    parts = _ext_parts(**{"ext/valence/a": b"1"})
    path = _write(tmp_path, parts=parts, manifest=_m(parts, ext={"valence": {}}),
                  extra_entries=[("ext/valence/stray", b"x")])
    _refuses(path, "unlisted-member")


def _limits(**kw):
    return getattr(chil_io, "ChilLimits")(**kw)


def test_ext_members_count_toward_member_limit(tmp_path):
    parts = _ext_parts(**{"ext/valence/a": b"1", "ext/valence/b": b"2"})
    path = _write(tmp_path, parts=parts, manifest=_m(parts, ext={"valence": {}}))
    _refuses(path, "limit-members", limits=_limits(max_members=3))


def test_ext_member_counts_toward_member_bytes(tmp_path):
    # Executive correction 2026-11-19: the ext member must exceed the limit
    # (len(f14) + 100) so that only it can trip limit-member-bytes.
    parts = _ext_parts(**{"ext/valence/a": b"x" * (len(_f14()) + 500)})
    path = _write(tmp_path, parts=parts, manifest=_m(parts, ext={"valence": {}}))
    _refuses(path, "limit-member-bytes", limits=_limits(max_member_bytes=len(_f14()) + 100))


def test_ext_members_count_toward_total_bytes(tmp_path):
    parts = _ext_parts(**{"ext/valence/a": b"x" * 500})
    path = _write(tmp_path, parts=parts, manifest=_m(parts, ext={"valence": {}}))
    _refuses(path, "limit-total-bytes", limits=_limits(max_total_bytes=len(_f14()) + 400))


# --------------------------------------------------------------- T024 FR-062

def test_valid_archive_reports_version_and_parts(tmp_path):
    parts = {"fort.14": _f14(), "fort.13": b"a\n"}
    report = getattr(chil_io, "validate_chil")(_write(tmp_path, parts=parts))
    assert report.format_version == 0
    assert [(p.name, p.kind, p.role, p.size, p.sha256) for p in report.parts] == [
        ("fort.14", "fort14", "source", len(_f14()), _sha(_f14())),
        ("fort.13", "fort13", "source", 2, _sha(b"a\n")),
    ]


def _bad_archives(tmp_path):
    p = tmp_path / "notzip.chil"
    p.write_bytes(b"garbage" * 20)
    m = default_manifest({"fort.14": _f14()})
    m["parts"][0]["sha256"] = "0" * 64
    return [
        str(p),
        _write(tmp_path, name="a.chil", order=["fort.14"]),
        _write(tmp_path, name="b.chil", manifest=b"{"),
        _write(tmp_path, name="c.chil", manifest=m),
        _write(tmp_path, name="d.chil", manifest={"c1_version": 1}),
        _write(tmp_path, name="e.chil", parts={"fort.14": _tiny_fort14("1 3 0 2 3")}),
    ]


def test_every_refusal_is_a_value_error(tmp_path):
    validate = getattr(chil_io, "validate_chil")
    err = getattr(chil_io, "ChilFormatError")
    for path in _bad_archives(tmp_path):
        with pytest.raises(ValueError) as ei:
            validate(path)
        assert isinstance(ei.value, err)
        assert str(ei.value).startswith(f"{ei.value.check}: ")


# -------------------------------------------------------------- T024b read_member

def test_read_member_returns_verified_bytes(tmp_path):
    parts = {"fort.14": _f14(), "fort.13": b"companion\n"}
    path = _write(tmp_path, parts=parts)
    read_member = getattr(chil_io, "read_member")
    assert read_member(path, "fort.13") == b"companion\n"
    assert read_member(path, "fort13") == b"companion\n"
    assert read_member(path, "fort14") == _f14()


def test_read_member_unlisted_name_is_request(tmp_path):
    path = _write(tmp_path)
    _read_refuses(path, "nothing.txt", "request")


def test_read_member_absent_kind_is_request(tmp_path):
    path = _write(tmp_path)
    _read_refuses(path, "fort13", "request")


def _read_refuses(path, which, slug):
    err = getattr(chil_io, "ChilFormatError")
    with pytest.raises(err) as ei:
        getattr(chil_io, "read_member")(path, which)
    assert ei.value.check == slug


def test_read_member_digest_mismatch(tmp_path):
    m = default_manifest({"fort.14": _f14()})
    m["parts"][0]["sha256"] = "0" * 64
    _read_refuses(_write(tmp_path, manifest=m), "fort.14", "digest-mismatch")


def test_read_member_prototype_container(tmp_path):
    _read_refuses(_write(tmp_path, manifest={"c1_version": 1}), "fort.14",
                  "prototype-container")


# ------------------------------------------------- executive ruling 2026-11-19

def test_archive_changed_during_validation_refused(tmp_path, monkeypatch):
    """archive_sha256 must describe the bytes checked: an in-place change made
    during validation (simulated deterministically, no race) is refused."""
    chil_io = importlib.import_module("chilmesh.chil_io")
    parts = {"fort.14": _f14()}
    path = _write(tmp_path, parts=parts, manifest=_m(parts))
    real = chil_io._validate_open

    def tamper(*args, **kwargs):
        result = real(*args, **kwargs)
        with open(path, "r+b") as f:
            f.seek(-1, 2)
            last = f.read(1)
            f.seek(-1, 2)
            f.write(bytes([last[0] ^ 0xFF]))
        return result

    monkeypatch.setattr(chil_io, "_validate_open", tamper)
    with pytest.raises(chil_io.ChilFormatError) as ei:
        chil_io.validate_chil(path)
    assert ei.value.check == "archive-changed"


def test_unparseable_boundary_block_with_bad_id_before_it_is_refused(tmp_path):
    """Executive correction 2026-11-19: ids read before the first bad row are scanned."""
    base = _tiny_fort14().decode().rstrip("\n").split("\n")[:-4]
    data = ("\n".join(base + ["1", "2", "2 0", "0", "garbage"]) + "\n").encode()
    _refuses(_write(tmp_path, parts={"fort.14": data}), "node-id-zero")
