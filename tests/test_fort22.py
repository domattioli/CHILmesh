"""fort.22 best-track read/write (spec 005 US7, FR-015)."""
from __future__ import annotations

import datetime as dt
from pathlib import Path

import pytest

from chilmesh import read_fort22_besttrack, write_fort22_besttrack

ATCF = Path(__file__).parent / "fixtures" / "atcf" / "bal092011.dat"

ROW8 = "AL, 09, 2011082100,   , BEST,   0, 150N,  590W,  45, 1006, TS,  34, NEQ,  105,    0,    0,   45, 1010,  175,  60,  55,   0,   L,   0,    ,   0,   0,      IRENE, M, 12, NEQ,  480,    0,    0,    0"
ROW_SHORT = "AL, 09, 2011082106,   , BEST,   0, 160S,  606E,    ,     , TS,  34, NEQ"


def test_typed_fields_from_real_atcf_file():
    tr = read_fort22_besttrack(ATCF)
    assert len(tr.records) == 99
    r = tr.records[0]
    assert r.time == dt.datetime(2011, 8, 21, 0)
    assert r.lat == pytest.approx(15.0) and r.lon == pytest.approx(-59.0)
    assert r.vmax == 45 and r.pmin == 1006 and r.isotach == 34
    assert r.radii == [105.0, 0.0, 0.0, 45.0]
    assert r.rmw == 60
    times = [x.time for x in tr.records]
    assert times == sorted(times)


def test_untyped_columns_kept_verbatim(tmp_path):
    p = tmp_path / "fort.22"
    p.write_text(ROW8 + "\n")
    r = read_fort22_besttrack(p).records[0]
    assert "IRENE" in r.extra and "NEQ" in r.extra
    assert r.extra.count("NEQ") == 2  # windcode and the NWS 20 style repeat, in order
    assert r.extra[:2] == ["AL", "09"]


def test_blank_and_short_fields_are_missing_not_zero(tmp_path):
    p = tmp_path / "fort.22"
    p.write_text(ROW_SHORT + "\n")
    r = read_fort22_besttrack(p).records[0]
    assert r.vmax is None and r.pmin is None
    assert r.rmw is None and r.radii == [None] * 4
    assert r.lat == pytest.approx(-16.0) and r.lon == pytest.approx(60.6)


@pytest.mark.parametrize("nws", [8, 20])
def test_round_trip_real_file(tmp_path, nws):
    tr = read_fort22_besttrack(ATCF, nws=nws)
    out = tmp_path / "fort.22"
    write_fort22_besttrack(tr, out, nws=nws)
    assert read_fort22_besttrack(out, nws=nws) == tr


def test_round_trip_synthetic_rows(tmp_path):
    p = tmp_path / "in"
    p.write_text(ROW8 + "\n" + ROW_SHORT + "\n")
    tr = read_fort22_besttrack(p)
    out = tmp_path / "out"
    write_fort22_besttrack(tr, out, nws=20)
    again = read_fort22_besttrack(out, nws=20)
    assert again == tr
    assert again.records[1].extra == tr.records[1].extra


@pytest.mark.parametrize("bad", [0, 3, 12, 30, "x"])
def test_unsupported_nws_error(tmp_path, bad):
    with pytest.raises(ValueError, match="NWS"):
        read_fort22_besttrack(ATCF, nws=bad)
    tr = read_fort22_besttrack(ATCF)
    with pytest.raises(ValueError, match="NWS"):
        write_fort22_besttrack(tr, tmp_path / "o", nws=bad)


def test_not_a_besttrack_file_error(tmp_path):
    p = tmp_path / "fort.22"
    p.write_text("Oceanweather WIN/PRE Format   2005082500   2005083121\n")
    with pytest.raises(ValueError, match="best-track"):
        read_fort22_besttrack(p)
