"""Summary tests for .chil (T038-T040; FR-030..FR-033, SC-004, specs/010-chil-container-mvp)."""
from __future__ import annotations

import functools
import importlib
import zipfile

import pytest

from _chil_builders import (archive_with_summary, build_archive, fresh,
                            good_fort14_bytes, summary_fields)


def _chil_io():
    return importlib.import_module("chilmesh.chil_io")


def _summary_chil():
    return getattr(_chil_io(), "summary_chil")


def _err():
    return getattr(_chil_io(), "ChilFormatError")


@functools.lru_cache(maxsize=1)
def _f14() -> bytes:
    return good_fort14_bytes()


@functools.lru_cache(maxsize=1)
def _cache() -> dict:
    return summary_fields(_f14())


@functools.lru_cache(maxsize=1)
def _annulus_type() -> str:
    return fresh("annulus").type


def _write(tmp_path, data: bytes, name="a.chil") -> str:
    p = tmp_path / name
    p.write_bytes(data)
    return str(p)


@pytest.fixture
def opened(monkeypatch):
    names: list[str] = []
    orig = zipfile.ZipFile.open

    def spy(self, name, mode="r", *a, **k):
        if mode == "r":
            names.append(getattr(name, "filename", name))
        return orig(self, name, mode, *a, **k)

    monkeypatch.setattr(zipfile.ZipFile, "open", spy)
    return names


@pytest.fixture
def no_parse(monkeypatch):
    from chilmesh import CHILmesh

    _cache()  # prime before patching
    _annulus_type()

    def boom(*a, **k):
        raise AssertionError("fort.14 payload was parsed")

    monkeypatch.setattr(CHILmesh, "_from_fort14_lines", classmethod(boom))


# ---- T038: non-deep returns cached fields, opens only the manifest (FR-031, SC-004)

def test_nondeep_returns_cached_fields(tmp_path, no_parse):
    path = _write(tmp_path, archive_with_summary(_f14()))
    r = _summary_chil()(path)
    c = _cache()
    assert r["format"] == "chil"
    assert r["n_nodes"] == c["n_nodes"]
    assert r["n_elems"] == c["n_elements"]
    assert r["n_triangles"] == c["n_triangles"] and r["n_quads"] == c["n_quads"]
    assert r["grid_name"] == c["grid_name"]
    assert r["element_type"] == _annulus_type()
    assert r["cached_summary"] == c
    xmin, ymin, xmax, ymax = c["bbox"]
    assert r["bbox"] == {"min_x": xmin, "max_x": xmax, "min_y": ymin, "max_y": ymax}
    assert "mismatches" not in r and "recomputed_summary" not in r
    assert r["verified"] is False


def test_nondeep_returns_core_fields(tmp_path):
    path = _write(tmp_path, archive_with_summary(_f14()))
    r = _summary_chil()(path)
    assert r["format_version"] == 0 and r["crs"] == "unknown"
    assert [p["name"] for p in r["parts"]] == ["fort.14"]
    assert set(r["parts"][0]) >= {"name", "role", "kind", "size", "sha256"}
    assert r["domain_boundary"] is None
    assert r["ext_namespaces"] == []
    assert r["source_archive_sha256"] is None
    assert r["file_bytes"] == (tmp_path / "a.chil").stat().st_size
    assert r["path"] == path


def test_nondeep_opens_only_manifest(tmp_path, opened):
    path = _write(tmp_path, archive_with_summary(_f14()))
    _summary_chil()(path)
    assert opened, "manifest.json was never opened through ZipFile.open"
    assert set(opened) == {"manifest.json"}


def test_nondeep_with_companions_opens_only_manifest(tmp_path, opened):
    parts = {"fort.14": _f14(), "fort.13": b"attrs\n"}
    from _chil_builders import default_manifest
    m = default_manifest(parts, summary=_cache())
    path = _write(tmp_path, build_archive(parts, manifest=m))
    r = _summary_chil()(path)
    assert set(opened) == {"manifest.json"}
    assert [p["name"] for p in r["parts"]] == ["fort.14", "fort.13"]


# ---- T039: no cache; payload not parsed; corruption not detected (FR-031)

def test_nondeep_no_cache_core_only(tmp_path, opened, no_parse):
    path = _write(tmp_path, archive_with_summary(_f14(), summary=None))
    r = _summary_chil()(path)
    assert r["cached_summary"] is None
    for k in ("n_nodes", "n_elems", "n_triangles", "n_quads", "element_type",
              "bbox", "grid_name"):
        assert k not in r, k
    assert r["format_version"] == 0 and r["crs"] == "unknown"
    assert [p["name"] for p in r["parts"]] == ["fort.14"]
    assert "fort.14" not in opened
    assert set(opened) == {"manifest.json"}


def test_nondeep_does_not_detect_corrupt_payload(tmp_path):
    bad = b"garbage, not a fort.14\n"
    path = _write(tmp_path, archive_with_summary(bad, summary=_cache()))
    r = _summary_chil()(path)  # must not raise
    assert r["n_nodes"] == _cache()["n_nodes"]
    # deep mode does notice
    with pytest.raises(_err()):
        _summary_chil()(path, deep=True)


@pytest.mark.parametrize("manifest,slug", [
    (b"{not json", "manifest-json"),
    ({"c1_version": 1}, "prototype-container"),
    ({"format_version": 1, "crs": "unknown", "parts": []}, "format-version"),
])
def test_nondeep_malformed_manifest_refused(tmp_path, manifest, slug):
    path = _write(tmp_path, build_archive({"fort.14": _f14()}, manifest=manifest))
    with pytest.raises(_err()) as ei:
        _summary_chil()(path)
    assert ei.value.check == slug


def test_nondeep_not_zip_refused(tmp_path):
    path = _write(tmp_path, b"not a zip" * 20)
    with pytest.raises(_err()) as ei:
        _summary_chil()(path)
    assert ei.value.check == "not-zip"


def test_nondeep_bad_summary_schema_refused(tmp_path):
    bad = {**_cache(), "n_nodes": "many"}
    path = _write(tmp_path, archive_with_summary(_f14(), summary=bad))
    with pytest.raises(_err()) as ei:
        _summary_chil()(path)
    assert ei.value.check == "manifest-schema"


# ---- T040: deep (FR-032), bbox omission (FR-030), dispatch (FR-033)

def test_deep_agreeing_cache_no_mismatch(tmp_path):
    path = _write(tmp_path, archive_with_summary(_f14()))
    r = _summary_chil()(path, deep=True)
    assert r["mismatches"] == []
    assert r["recomputed_summary"] == _cache()
    assert r["n_nodes"] == _cache()["n_nodes"]
    assert r["cached_summary"] == _cache()
    assert r["verified"] is True


def test_deep_reports_wrong_n_nodes(tmp_path):
    path = _write(tmp_path, archive_with_summary(_f14(), overrides={"n_nodes": 7,
                                                                    "n_elements": 7,
                                                                    "n_triangles": 7}))
    r = _summary_chil()(path, deep=True)
    assert r["mismatches"] == ["n_elements", "n_nodes", "n_triangles"]
    assert r["recomputed_summary"]["n_nodes"] == _cache()["n_nodes"]
    assert r["n_nodes"] == _cache()["n_nodes"]  # recomputed value wins
    assert r["cached_summary"]["n_nodes"] == 7


def test_deep_reports_wrong_grid_name_and_bbox(tmp_path):
    path = _write(tmp_path, archive_with_summary(
        _f14(), overrides={"grid_name": "other", "bbox": [0.0, 0.0, 1.0, 1.0]}))
    r = _summary_chil()(path, deep=True)
    assert r["mismatches"] == ["bbox", "grid_name"]


def test_deep_no_cache_computes_and_no_mismatches(tmp_path):
    path = _write(tmp_path, archive_with_summary(_f14(), summary=None))
    r = _summary_chil()(path, deep=True)
    assert r["cached_summary"] is None
    assert r["mismatches"] == []
    assert r["recomputed_summary"] == _cache()
    assert r["n_nodes"] == _cache()["n_nodes"]
    assert r["element_type"] == _annulus_type()
    assert set(r["bbox"]) == {"min_x", "max_x", "min_y", "max_y"}


def test_deep_validates_first(tmp_path):
    from _chil_builders import default_manifest
    m = default_manifest({"fort.14": _f14()})
    m["parts"][0]["sha256"] = "0" * 64
    path = _write(tmp_path, build_archive({"fort.14": _f14()}, manifest=m))
    with pytest.raises(_err()) as ei:
        _summary_chil()(path, deep=True)
    assert ei.value.check == "digest-mismatch"


def test_deep_parses_fort14_without_layers(tmp_path, monkeypatch):
    from chilmesh import CHILmesh
    seen = {}
    orig = CHILmesh._from_fort14_lines.__func__

    def spy(cls, lines, compute_layers=True, compute_adjacencies=None, **kw):
        seen["layers"], seen["adj"] = compute_layers, compute_adjacencies
        return orig(cls, lines, compute_layers, compute_adjacencies, **kw)

    monkeypatch.setattr(CHILmesh, "_from_fort14_lines", classmethod(spy))
    path = _write(tmp_path, archive_with_summary(_f14()))
    _summary_chil()(path, deep=True)
    assert seen == {"layers": False, "adj": False}


def test_deep_nonfinite_coordinate_refused(tmp_path):
    lines = _f14().decode().splitlines()
    parts = lines[2].split()  # first node row: id x y z
    parts[1] = "nan"
    lines[2] = " ".join(parts)
    bad = ("\n".join(lines) + "\n").encode()
    cache = {k: v for k, v in _cache().items() if k != "bbox"}
    path = _write(tmp_path, archive_with_summary(bad, summary=cache))
    with pytest.raises(_err()) as ei:
        _summary_chil()(path, deep=True)
    assert ei.value.check == "fort14-parse"


def test_deep_cache_without_bbox_is_not_a_mismatch(tmp_path):
    cache = {k: v for k, v in _cache().items() if k != "bbox"}
    path = _write(tmp_path, archive_with_summary(_f14(), summary=cache))
    r = _summary_chil()(path, deep=True)
    assert r["mismatches"] == []
    assert "bbox" in r["recomputed_summary"]


def test_nondeep_cache_without_bbox_has_no_bbox(tmp_path):
    cache = {k: v for k, v in _cache().items() if k != "bbox"}
    path = _write(tmp_path, archive_with_summary(_f14(), summary=cache))
    r = _summary_chil()(path)
    assert "bbox" not in r
    assert r["n_nodes"] == cache["n_nodes"]


def test_summary_io_dispatches_chil(tmp_path):
    summary_io = importlib.import_module("chilmesh.summary_io")
    path = _write(tmp_path, archive_with_summary(_f14()))
    r = summary_io.summary(path)
    assert r["format"] == "chil"
    assert r["n_nodes"] == _cache()["n_nodes"]
    rd = summary_io.summary(path, deep=True)
    assert rd["mismatches"] == [] and "recomputed_summary" in rd


def test_chilmesh_summary_dispatches_chil_uppercase_suffix(tmp_path):
    import chilmesh
    path = _write(tmp_path, archive_with_summary(_f14()), name="A.CHIL")
    assert chilmesh.summary(path)["format"] == "chil"


def test_summary_io_missing_chil_is_file_not_found(tmp_path):
    summary_io = importlib.import_module("chilmesh.summary_io")
    with pytest.raises(FileNotFoundError):
        summary_io.summary(tmp_path / "nope.chil")


def test_summary_io_wraps_chil_error_as_summary_error(tmp_path):
    summary_io = importlib.import_module("chilmesh.summary_io")
    path = _write(tmp_path, b"not a zip" * 20)
    with pytest.raises(summary_io.SummaryError) as ei:
        summary_io.summary(path)
    assert isinstance(ei.value.__cause__, _err())
    assert ei.value.__cause__.check == "not-zip"
    assert str(ei.value).startswith("not-zip: ")
