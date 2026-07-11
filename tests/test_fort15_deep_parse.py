"""Tests for the opt-in fort.15 deep-parse leading scalar block (#249)."""
from __future__ import annotations

from pathlib import Path

import pytest

from chilmesh.fort15_io import read_fort15, write_fort15

FIXTURE = Path(__file__).parent / "fixtures" / "fort15" / "sample.15"


def test_deep_false_leaves_fields_none():
    f15 = read_fort15(FIXTURE)
    assert f15.im == 0
    assert f15.dtdp is None
    assert f15.tau0 is None
    assert f15.nodal_attribute_names is None


def test_deep_true_parses_leading_scalar_block():
    f15 = read_fort15(FIXTURE, deep=True)
    # header still correct
    assert f15.ics == 2
    assert f15.im == 0
    # leading scalar block through REFTIM
    assert f15.nolibf == 0
    assert f15.nolifa == 1
    assert f15.nolica == 1
    assert f15.nolicat == 1
    assert f15.nwp == 0
    assert f15.nodal_attribute_names == []
    assert f15.ncor == 1
    assert f15.ntip == 0
    assert f15.nws == 0
    assert f15.nramp == 1
    assert f15.g == pytest.approx(9.81)
    assert f15.tau0 == pytest.approx(0.0)
    assert f15.dtdp == pytest.approx(2.0)
    assert f15.statim == pytest.approx(0.0)
    assert f15.reftim == pytest.approx(0.0)


def test_deep_parse_preserves_byte_exact_write(tmp_path):
    f15 = read_fort15(FIXTURE, deep=True)
    out = tmp_path / "out.15"
    write_fort15(f15, out)
    assert out.read_text() == FIXTURE.read_text()


def test_deep_parse_reads_nwp_attribute_name_block(tmp_path):
    """NWP>0 → the following NWP lines are nodal-attribute names, not scalars."""
    src = tmp_path / "attr.15"
    src.write_text(
        "run desc                 ! RUNDES\n"
        "runid                    ! RUNID\n"
        "1                        ! NFOVER\n"
        "1                        ! NABOUT\n"
        "1                        ! NSCREEN\n"
        "0                        ! IHOT\n"
        "2                        ! ICS\n"
        "0                        ! IM\n"
        "0                        ! NOLIBF\n"
        "1                        ! NOLIFA\n"
        "1                        ! NOLICA\n"
        "1                        ! NOLICAT\n"
        "2                        ! NWP\n"
        "mannings_n_at_sea_floor\n"
        "primitive_weighting_in_continuity_equation\n"
        "1                        ! NCOR\n"
        "0                        ! NTIP\n"
        "0                        ! NWS\n"
        "1                        ! NRAMP\n"
        "9.81                     ! G\n"
        "0.0                      ! TAU0\n"
        "1.0                      ! DTDP\n"
        "0.0                      ! STATIM\n"
        "0.0                      ! REFTIM\n"
    )
    f15 = read_fort15(src, deep=True)
    assert f15.nwp == 2
    assert f15.nodal_attribute_names == [
        "mannings_n_at_sea_floor",
        "primitive_weighting_in_continuity_equation",
    ]
    # scalars after the name block must line up (proves index bookkeeping)
    assert f15.ncor == 1
    assert f15.dtdp == pytest.approx(1.0)
    assert f15.reftim == pytest.approx(0.0)
