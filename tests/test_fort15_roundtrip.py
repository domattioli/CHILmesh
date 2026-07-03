"""Tests for chilmesh.fort15_io (ADCIRC fort.15 byte-preserving I/O)."""
from __future__ import annotations

from pathlib import Path

import pytest

from chilmesh.fort15_io import Fort15, read_fort15, write_fort15, Fort15ParseError

FIXTURE = Path(__file__).parent / "fixtures" / "fort15" / "sample.15"


def test_header_extraction():
    f15 = read_fort15(FIXTURE)
    assert f15.rundes == "Sample WNAT tidal run"
    assert f15.runid == "sample01"
    assert f15.nfover == 1
    assert f15.nabout == 1
    assert f15.nscreen == 1
    assert f15.ihot == 0
    assert f15.ics == 2
    assert f15.im == 0


def test_roundtrip_byte_identical(tmp_path):
    f15 = read_fort15(FIXTURE)
    out = tmp_path / "out.15"
    write_fort15(f15, out)
    assert out.read_bytes() == FIXTURE.read_bytes()


def test_comment_stripping_on_int_field(tmp_path):
    # A control line with an inline comment must parse the leading integer only.
    content = (
        "desc\n"
        "id\n"
        "1 ! NFOVER\n"
        "0 ! NABOUT\n"
        "1 ! NSCREEN\n"
        "0 ! IHOT\n"
        "22 ! ICS spherical\n"
        "0 ! IM\n"
    )
    p = tmp_path / "c.15"
    p.write_text(content)
    f15 = read_fort15(p)
    assert f15.ics == 22
    assert f15.rundes == "desc"


def test_too_short_raises(tmp_path):
    p = tmp_path / "short.15"
    p.write_text("one\ntwo\nthree\n")
    with pytest.raises(Fort15ParseError):
        read_fort15(p)


def test_bad_int_field_raises(tmp_path):
    content = "d\ni\nX\n1\n1\n0\n2\n0\n"  # NFOVER = "X" is non-integer
    p = tmp_path / "bad.15"
    p.write_text(content)
    with pytest.raises(Fort15ParseError):
        read_fort15(p)


def test_empty_raw_text_write_raises(tmp_path):
    f15 = Fort15("d", "i", 1, 1, 1, 0, 2, 0, raw_text="")
    with pytest.raises(Fort15ParseError):
        write_fort15(f15, tmp_path / "x.15")
