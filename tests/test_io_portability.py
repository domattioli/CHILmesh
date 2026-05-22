"""Cross-platform I/O byte-invariant tests (#121).

Locks in the behaviour added by the explicit ``encoding`` / ``newline``
arguments on the three ``open()`` sites in ``chilmesh.CHILmesh``: text
is read as UTF-8 with universal-newline translation, and written as
UTF-8 with a forced ``\\n`` terminator so the output is byte-identical
on Linux, macOS, and Windows.
"""
from __future__ import annotations

import numpy as np
import pytest

import chilmesh
from chilmesh import CHILmesh
from conftest import TRI_FIXTURE_NAMES as FIXTURES


@pytest.mark.parametrize("name", FIXTURES)
def test_write_fort14_uses_lf_line_endings(name, tmp_path):
    mesh = getattr(chilmesh.examples, name)()
    out = tmp_path / f"{name}.fort.14"
    assert mesh.write_to_fort14(str(out), grid_name=mesh.grid_name or name)

    raw = out.read_bytes()
    assert b"\r\n" not in raw, f"{name}: CRLF found in output (Windows-only artifact)"
    assert b"\r" not in raw, f"{name}: stray CR byte in output"
    assert raw.endswith(b"\n"), f"{name}: output missing trailing LF"


@pytest.mark.parametrize("name", FIXTURES)
def test_read_fort14_handles_crlf_line_endings(name, tmp_path):
    mesh = getattr(chilmesh.examples, name)()
    lf_path = tmp_path / f"{name}.lf.fort.14"
    crlf_path = tmp_path / f"{name}.crlf.fort.14"
    assert mesh.write_to_fort14(str(lf_path), grid_name=mesh.grid_name or name)

    crlf_path.write_bytes(lf_path.read_bytes().replace(b"\n", b"\r\n"))
    assert b"\r\n" in crlf_path.read_bytes()

    reloaded_lf = CHILmesh.read_from_fort14(lf_path)
    reloaded_crlf = CHILmesh.read_from_fort14(crlf_path)

    np.testing.assert_array_equal(reloaded_lf.points, reloaded_crlf.points)
    np.testing.assert_array_equal(
        reloaded_lf.connectivity_list, reloaded_crlf.connectivity_list
    )


def test_write_fort14_is_utf8_encoded(tmp_path):
    mesh = chilmesh.examples.annulus()
    out = tmp_path / "annulus.utf8.fort.14"
    header = "Annulus — UTF-8 header µ"
    assert mesh.write_to_fort14(str(out), grid_name=header)

    raw = out.read_bytes()
    assert raw.startswith(header.encode("utf-8")), "header not UTF-8 encoded"

    reloaded = CHILmesh.read_from_fort14(out)
    assert reloaded.grid_name.startswith("Annulus")
