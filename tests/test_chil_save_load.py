"""Save/load tests for the .chil container (specs/010-chil-container-mvp).

T004 guard: the fort.14 reader/writer refactors (R1, R2) are pure.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from _chil_builders import GOLDEN_DIR, fresh, mixed_mesh
from chilmesh.CHILmesh import CHILmesh, write_fort14

_MESHES = {
    "annulus": lambda: fresh("annulus"),
    "quad": lambda: fresh("quad_2x2"),
    "mixed": mixed_mesh,
}


@pytest.mark.parametrize("name", sorted(_MESHES))
def test_write_fort14_bytes_equal_golden(name, tmp_path):
    out = tmp_path / f"{name}.14"
    assert write_fort14(_MESHES[name](), str(out)) is True
    assert out.read_bytes() == (GOLDEN_DIR / f"{name}.fort14").read_bytes()


@pytest.mark.parametrize("name", sorted(_MESHES))
def test_read_from_fort14_parses_golden(name):
    golden = GOLDEN_DIR / f"{name}.fort14"
    mesh = CHILmesh.read_from_fort14(golden)
    ref = _MESHES[name]()
    assert mesh.n_verts == ref.n_verts
    assert mesh.n_elems == ref.n_elems


# ---------------------------------------------------------------------------
# US1: save/load (T030 to T034b)
# ---------------------------------------------------------------------------
import hashlib
import importlib
import json
import os
import time
import zipfile
from unittest import mock

import numpy as np

import chilmesh
from _chil_builders import build_archive, default_manifest, good_fort14_bytes


def _chil_io():
    return importlib.import_module("chilmesh.chil_io")


def _err():
    return _chil_io().ChilFormatError


def _golden(name: str) -> bytes:
    return (GOLDEN_DIR / f"{name}.fort14").read_bytes()


def _manifest_of(path) -> dict:
    with zipfile.ZipFile(path) as zf:
        return json.loads(zf.read("manifest.json"))


def _segs_equal(a, b) -> bool:
    if len(a) != len(b):
        return False
    for sa, sb in zip(a, b):
        if set(sa) != set(sb):
            return False
        for k in sa:
            if not np.array_equal(np.asarray(sa[k]), np.asarray(sb[k])):
                return False
    return True


def _source_archive_bytes() -> bytes:
    f14 = good_fort14_bytes()
    return build_archive(
        {"fort.14": f14, "fort.13": b"companion\n"},
        ext_member={"ext/example/note.bin": b"\x00\x01note"},
        ext={"example": {"k": 1}},
        manifest=None,
    )


# ---- T030 (FR-015, SC-002) ------------------------------------------------

@pytest.mark.parametrize("name", sorted(_MESHES))
def test_save_chil_layout_manifest_first_one_generated_fort14(name, tmp_path):
    """FR-015: manifest first, exactly one fort14 part with role generated."""
    out = tmp_path / "m.chil"
    _MESHES[name]().save(str(out))
    with zipfile.ZipFile(out) as zf:
        names = zf.namelist()
    assert names[0] == "manifest.json"
    parts = _manifest_of(out)["parts"]
    assert len(parts) == 1
    assert parts[0]["kind"] == "fort14" and parts[0]["role"] == "generated"
    assert names == ["manifest.json", parts[0]["name"]]


@pytest.mark.parametrize("name", sorted(_MESHES))
def test_save_chil_part_bytes_equal_write_fort14(name, tmp_path):
    """SC-002: part bytes, size and sha256 equal the write_fort14 file bytes."""
    out = tmp_path / "m.chil"
    mesh = _MESHES[name]()
    mesh.save(str(out))
    ref = tmp_path / "ref.14"
    assert write_fort14(_MESHES[name](), str(ref)) is True
    expected = ref.read_bytes()
    part = _manifest_of(out)["parts"][0]
    with zipfile.ZipFile(out) as zf:
        assert zf.read(part["name"]) == expected
    assert part["size"] == len(expected)
    assert part["sha256"] == hashlib.sha256(expected).hexdigest()


# ---- T031 (FR-018..FR-021) ------------------------------------------------

@pytest.mark.parametrize("name", sorted(_MESHES))
def test_load_restores_mesh(name, tmp_path):
    """FR-018: coordinates, connectivity, boundaries, element types and counts survive."""
    src = _MESHES[name]()
    out = tmp_path / "m.chil"
    src.save(str(out))
    ref_path = tmp_path / "ref.14"
    write_fort14(src, str(ref_path))
    ref = CHILmesh.read_from_fort14(ref_path)
    got = CHILmesh.load(str(out))
    assert got.n_verts == ref.n_verts and got.n_elems == ref.n_elems
    assert np.array_equal(np.asarray(got.points), np.asarray(ref.points))
    assert np.array_equal(np.asarray(got.connectivity_list), np.asarray(ref.connectivity_list))
    assert _segs_equal(got.boundary_segments, ref.boundary_segments)
    gt, gq = got._elem_type()
    rt, rq = ref._elem_type()
    assert (len(gt), len(gq)) == (len(rt), len(rq))


@pytest.mark.parametrize("layers", [True, False])
@pytest.mark.parametrize("adj", [True, False])
def test_load_flag_combinations_match_read_from_fort14(layers, adj, tmp_path):
    """FR-021: all four compute flag combinations behave as read_from_fort14."""
    out = tmp_path / "m.chil"
    fresh("annulus").save(str(out))
    ref = CHILmesh.read_from_fort14(
        GOLDEN_DIR / "annulus.fort14", compute_layers=layers, compute_adjacencies=adj)
    got = CHILmesh.load(str(out), compute_layers=layers, compute_adjacencies=adj)
    assert np.array_equal(np.asarray(got.points), np.asarray(ref.points))
    assert np.array_equal(np.asarray(got.connectivity_list), np.asarray(ref.connectivity_list))
    assert hasattr(got, "adjacencies") == hasattr(ref, "adjacencies")
    assert getattr(got, "n_layers", None) == getattr(ref, "n_layers", None)
    for attr in ("Elem2Vert", "Edge2Vert", "Vert2Edge", "Vert2Elem", "Elem2Edge"):
        assert (attr in getattr(got, "adjacencies", {})) == (attr in getattr(ref, "adjacencies", {}))


def _assert_refused_without_parse(path, slug):
    with mock.patch.object(CHILmesh, "_from_fort14_lines",
                           side_effect=AssertionError("parser must not run")) as m:
        with pytest.raises(_err()) as ei:
            CHILmesh.load(str(path))
    assert ei.value.check == slug
    m.assert_not_called()


def test_load_digest_mismatch_refused_before_parse(tmp_path):
    """FR-020: validator runs before any member is parsed (digest-mismatch)."""
    good = good_fort14_bytes()
    man = default_manifest({"fort.14": good})
    man["parts"][0]["sha256"] = "0" * 64
    p = tmp_path / "bad.chil"
    p.write_bytes(build_archive({"fort.14": good}, manifest=man))
    _assert_refused_without_parse(p, "digest-mismatch")


def test_load_prototype_container_refused_before_parse(tmp_path):
    """FR-020: a c1_version archive is refused with prototype-container."""
    good = good_fort14_bytes()
    man = default_manifest({"fort.14": good})
    man["c1_version"] = 1
    p = tmp_path / "proto.chil"
    p.write_bytes(build_archive({"fort.14": good}, manifest=man))
    _assert_refused_without_parse(p, "prototype-container")


# ---- T032 (FR-016, FR-017, FR-019) ----------------------------------------

def test_save_after_load_of_source_archive_drops_everything(tmp_path):
    """FR-016: re-save of a source archive has one generated part and no carry-over."""
    src = tmp_path / "src.chil"
    src.write_bytes(_source_archive_bytes())
    mesh = CHILmesh.load(str(src))
    out = tmp_path / "out.chil"
    mesh.save(str(out))
    m = _manifest_of(out)
    assert [(p["kind"], p["role"]) for p in m["parts"]] == [("fort14", "generated")]
    assert "source_archive_sha256" not in m
    assert "ext" not in m
    with zipfile.ZipFile(out) as zf:
        names = zf.namelist()
    assert not any(n.startswith("ext/") or n == "fort.13" for n in names)
    assert len(names) == 2


def test_save_default_crs_is_unknown(tmp_path):
    """FR-017: crs=None records 'unknown'."""
    out = tmp_path / "m.chil"
    fresh("annulus").save(str(out))
    assert _manifest_of(out)["crs"] == "unknown"


# ---- T033 (FR-007, SC-007, FR-017) ----------------------------------------

def test_two_saves_byte_identical(tmp_path):
    """SC-007: same mesh saved twice gives identical bytes."""
    mesh = fresh("annulus")
    a, b = tmp_path / "a.chil", tmp_path / "b.chil"
    mesh.save(str(a))
    mesh.save(str(b))
    assert a.read_bytes() == b.read_bytes()


def test_save_bytes_independent_of_time_utime_and_directory(tmp_path):
    """SC-007: bytes do not depend on mtime, time.time or destination directory."""
    d1, d2 = tmp_path / "one", tmp_path / "two"
    d1.mkdir()
    d2.mkdir()
    a, b = d1 / "a.chil", d2 / "other_name.chil"
    fresh("annulus").save(str(a))
    os.utime(a, (1_000_000_000, 1_000_000_000))
    with mock.patch.object(time, "time", return_value=4_102_444_800.0):
        fresh("annulus").save(str(b))
    os.utime(b, (1_700_000_000, 1_700_000_000))
    assert a.read_bytes() == b.read_bytes()


def test_save_with_crs_records_label_and_keeps_coordinates(tmp_path):
    """FR-017: crs label is recorded; coordinates are not transformed."""
    mesh = fresh("annulus")
    before = np.array(mesh.points, copy=True)
    out = tmp_path / "m.chil"
    mesh.save(str(out), crs="EPSG:4326")
    assert _manifest_of(out)["crs"] == "EPSG:4326"
    assert np.array_equal(np.asarray(mesh.points), before)
    assert np.array_equal(np.asarray(CHILmesh.load(str(out)).points), before)


# ---- T034 (FR-017, FR-022, FR-069, R7) ------------------------------------

def test_failed_replace_keeps_old_file_and_leaves_no_temp(tmp_path):
    """R7: if os.replace fails the old file is intact and no temp file remains."""
    out = tmp_path / "m.chil"
    out.write_bytes(b"OLD CONTENT")
    with mock.patch("os.replace", side_effect=OSError("boom")):
        with pytest.raises(OSError):
            fresh("annulus").save(str(out))
    assert out.read_bytes() == b"OLD CONTENT"
    assert [p.name for p in tmp_path.iterdir()] == ["m.chil"]


def test_load_writes_nothing_to_working_directory(tmp_path, monkeypatch):
    """FR-022: load extracts nothing to disk."""
    out = tmp_path / "m.chil"
    fresh("annulus").save(str(out))
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.chdir(work)
    CHILmesh.load(str(out))
    assert list(work.iterdir()) == []
    assert sorted(p.name for p in tmp_path.iterdir()) == ["m.chil", "work"]


@pytest.mark.parametrize("suffix", [".14", ".fort14", ".2dm"])
def test_crs_with_non_chil_target_raises_value_error(suffix, tmp_path):
    """FR-017: a crs label is never silently dropped."""
    out = tmp_path / f"m{suffix}"
    with pytest.raises(ValueError):
        fresh("annulus").save(str(out), crs="EPSG:4326")
    assert not out.exists()


@pytest.mark.parametrize("limits_kw,slug", [
    ({"max_member_bytes": 100}, "limit-member-bytes"),
    ({"max_total_bytes": 100}, "limit-total-bytes"),
    ({"max_members": 1}, "limit-members"),
    ({"max_manifest_bytes": 10}, "limit-manifest-bytes"),
])
def test_save_over_limit_refused_without_output(limits_kw, slug, tmp_path):
    """FR-069: over-limit save raises the limit-* check, no output, old file intact."""
    limits = _chil_io().ChilLimits(**limits_kw)
    new = tmp_path / "new.chil"
    with pytest.raises(_err()) as ei:
        fresh("annulus").save(str(new), limits=limits)
    assert ei.value.check == slug
    assert not new.exists()
    old = tmp_path / "old.chil"
    old.write_bytes(b"OLD CONTENT")
    with pytest.raises(_err()) as ei:
        fresh("annulus").save(str(old), limits=limits)
    assert ei.value.check == slug
    assert old.read_bytes() == b"OLD CONTENT"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["old.chil"]


# ---- T034b (Q8, FR-071) ---------------------------------------------------

def test_top_level_names_unchanged():
    """Q8: no new top-level name in chilmesh."""
    base = json.loads((GOLDEN_DIR / "top_level_names.json").read_text())
    assert list(chilmesh.__all__) == base["all"]
    # dir() depends on which submodules other tests imported, so check it in a
    # fresh interpreter: importing chilmesh.chil_io adds only the chil_io
    # submodule attribute, and every baseline name still resolves.
    import subprocess, sys
    code = (
        "import json, chilmesh\n"
        "before = set(dir(chilmesh))\n"
        "import chilmesh.chil_io\n"
        "print(json.dumps(sorted(set(dir(chilmesh)) - before)))\n"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert set(json.loads(out.stdout)) <= {"chil_io"}
    for n in base["dir"]:
        assert hasattr(chilmesh, n), n
    for n in ("pack", "rewrite_chil", "validate_chil", "read_member", "summary_chil",
              "save_chil", "load_chil", "ChilFormatError", "ChilLimits", "ChilReport",
              "PartInfo", "MemberSpec"):
        assert n not in dir(chilmesh)


def test_chil_io_all_is_exactly_the_twelve_names():
    """contracts/python-api.md: chil_io.__all__ lists the 12 public names."""
    mod = _chil_io()
    assert sorted(mod.__all__) == sorted([
        "pack", "rewrite_chil", "validate_chil", "read_member", "summary_chil",
        "save_chil", "load_chil", "ChilFormatError", "ChilLimits", "ChilReport",
        "PartInfo", "MemberSpec"])
    assert len(mod.__all__) == 12
