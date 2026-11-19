"""CLI tests for the ``validate`` verb and ``.chil`` handling (T024c, contracts/cli.md)."""
from __future__ import annotations

import functools
import hashlib
import importlib
import re

import pytest

from _chil_builders import build_archive, default_manifest, good_fort14_bytes

cli = importlib.import_module("chilmesh.cli")


@functools.lru_cache(maxsize=1)
def _f14() -> bytes:
    return good_fort14_bytes()


def _good(tmp_path, name="ok.chil") -> str:
    p = tmp_path / name
    p.write_bytes(build_archive({"fort.14": _f14()}))
    return str(p)


def _malformed(tmp_path) -> list[tuple[str, str]]:
    """(path, expected slug) for three malformed archives."""
    out = []
    p = tmp_path / "notzip.chil"
    p.write_bytes(b"not a zip" * 20)
    out.append((str(p), "not-zip"))
    m = default_manifest({"fort.14": _f14()})
    m["parts"][0]["sha256"] = "0" * 64
    p = tmp_path / "digest.chil"
    p.write_bytes(build_archive({"fort.14": _f14()}, manifest=m))
    out.append((str(p), "digest-mismatch"))
    p = tmp_path / "proto.chil"
    p.write_bytes(build_archive({"fort.14": _f14()}, manifest={"c1_version": 1}))
    out.append((str(p), "prototype-container"))
    return out


def _flat(text: str) -> str:
    return " ".join(text.split())


def test_validate_valid_archive(tmp_path, capsys):
    path = _good(tmp_path)
    assert cli.main(["validate", path]) == 0
    out = capsys.readouterr()
    assert out.err == ""
    assert f"OK: {path}" in out.out
    assert "format_version: 0" in out.out
    sha = hashlib.sha256(_f14()).hexdigest()
    assert re.search(rf"fort\.14\s+fort14\s+source\s+{len(_f14())} bytes\s+sha256 {sha[:4]}", out.out)


def test_validate_malformed_exit_1_and_stderr(tmp_path, capsys):
    for path, slug in _malformed(tmp_path):
        assert cli.main(["validate", path]) == 1
        out = capsys.readouterr()
        assert out.out == ""
        assert out.err.startswith(f"error: {slug}: "), (slug, out.err)


def test_validate_missing_file_exit_2(tmp_path, capsys):
    assert cli.main(["validate", str(tmp_path / "nope.chil")]) == 2
    out = capsys.readouterr()
    assert out.out == ""
    assert out.err.startswith("error: ")


def test_top_level_help_lists_validate(capsys):
    with pytest.raises(SystemExit) as ei:
        cli.main(["--help"])
    assert ei.value.code == 0
    assert "validate" in capsys.readouterr().out


def test_validate_help_says_experimental(capsys):
    with pytest.raises(SystemExit) as ei:
        cli.main(["validate", "--help"])
    assert ei.value.code == 0
    out = _flat(capsys.readouterr().out)
    assert "Validate a .chil archive (experimental format, format_version 0)." in out


@pytest.mark.parametrize("argv", [
    ["info", "{p}"],
    ["convert", "{p}", "{o}.14"],
    ["smooth", "{p}", "-o", "{o}.14"],
    ["plot", "{p}", "-o", "{o}.png"],
])
def test_other_verbs_reject_chil_input(tmp_path, capsys, argv):
    path = _good(tmp_path)
    args = [a.format(p=path, o=str(tmp_path / "out")) for a in argv]
    with pytest.raises(SystemExit) as ei:
        rc = cli.main(args)
        pytest.fail(f"expected a non-zero exit, got return {rc}")
    code = ei.value.code
    assert code not in (0, None)
    assert "cannot infer mesh format" in f"{code} {capsys.readouterr().err}"


# ---- T040b: `chilmesh summary x.chil` (contracts/cli.md, FR-033) ----

from _chil_builders import archive_with_summary  # noqa: E402


def _wr(tmp_path, data: bytes, name="s.chil") -> str:
    p = tmp_path / name
    p.write_bytes(data)
    return str(p)


def _lines(text: str) -> list[str]:
    return [" ".join(ln.split()) for ln in text.splitlines()]


def test_summary_plain_prints_chil_lines(tmp_path, capsys):
    path = _wr(tmp_path, archive_with_summary(_f14()))
    assert cli.main(["summary", path]) == 0
    out = capsys.readouterr()
    assert out.err == ""
    ln = _lines(out.out)
    assert "Format: chil" in ln
    assert "Version: 0" in ln
    assert "CRS: unknown" in ln
    assert "Parts: fort.14 (generated)" not in ln  # built archive part role is source
    assert "Parts: fort.14 (source)" in ln
    assert "Cache: yes (declared in manifest, not checked; use --deep)" in ln
    assert not any(x.startswith("Mismatch:") for x in ln)
    assert any(x.startswith("Nodes: ") and x != "Nodes: N/A" for x in ln)
    assert out.out.count("Cache:      yes") == 1  # exact column alignment
    assert "  Cache:      yes (declared in manifest, not checked; use --deep)\n" in out.out


def test_summary_plain_no_cache_hint(tmp_path, capsys):
    path = _wr(tmp_path, archive_with_summary(_f14(), summary=None))
    assert cli.main(["summary", path]) == 0
    out = capsys.readouterr().out
    assert "  Cache:      none (run with --deep to compute from the payload)" in out
    assert "Mismatch:" not in out


def test_summary_deep_wrong_cache_lists_field(tmp_path, capsys):
    path = _wr(tmp_path, archive_with_summary(_f14(), overrides={"n_nodes": 7}))
    assert cli.main(["summary", path, "--deep"]) == 0
    out = capsys.readouterr().out
    assert "  Mismatch:   n_nodes" in out
    assert "Cache:      yes" in out
    assert "Mismatch:   none" not in out


def test_summary_deep_agreeing_cache_prints_none(tmp_path, capsys):
    path = _wr(tmp_path, archive_with_summary(_f14()))
    assert cli.main(["summary", path, "--deep"]) == 0
    out = capsys.readouterr().out
    assert "  Mismatch:   none" in out
    assert "Cache:      yes" in out


def test_summary_deep_no_cache_output(tmp_path, capsys):
    path = _wr(tmp_path, archive_with_summary(_f14(), summary=None))
    assert cli.main(["summary", path, "--deep"]) == 0
    out = capsys.readouterr().out
    assert "  Cache:      none\n" in out
    assert "run with --deep" not in out
    assert "  Mismatch:   n/a (no cached summary)" in out
    ln = _lines(out)
    assert any(x.startswith("Nodes: ") and x != "Nodes: N/A" for x in ln)
    assert any(x.startswith("Elements: ") and x != "Elements: N/A" for x in ln)
    assert any(x.startswith("Type: ") for x in ln)
    assert any(x.startswith("Bbox: ") for x in ln)


@pytest.mark.parametrize("extra", [[], ["--deep"]])
def test_summary_invalid_archive_exit_1(tmp_path, capsys, extra):
    for path, slug in _malformed(tmp_path):
        if extra == [] and slug == "digest-mismatch":
            continue  # non-deep does not check digests (steps 1 to 7 only)
        assert cli.main(["summary", path, *extra]) == 1, slug
        out = capsys.readouterr()
        assert out.out == ""
        assert out.err.startswith(f"error: {slug}: "), (slug, out.err)


@pytest.mark.parametrize("extra", [[], ["--deep"]])
def test_summary_missing_chil_exit_2(tmp_path, capsys, extra):
    assert cli.main(["summary", str(tmp_path / "nope.chil"), *extra]) == 2
    out = capsys.readouterr()
    assert out.out == ""
    assert out.err.startswith("error: ")


def test_summary_help_mentions_chil_experimental(capsys):
    with pytest.raises(SystemExit) as ei:
        cli.main(["summary", "--help"])
    assert ei.value.code == 0
    out = _flat(capsys.readouterr().out)
    assert "Input mesh path (.14, .fort.14, .grd, .2dm, .chil; .chil is experimental)." in out
    assert "For .chil: validate and recompute from the payload." in out
