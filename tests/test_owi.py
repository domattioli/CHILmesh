"""OWI fort.221 / fort.222 read/write (spec 005 US7, FR-016)."""
from __future__ import annotations

import datetime as dt
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent / "fixtures"))
import synthetic_adcirc as syn  # noqa: E402

from chilmesh import read_owi, write_owi  # noqa: E402


def _snaps(wind):
    out = []
    for k, t in enumerate(["200508250000", "200508250600"]):
        a = np.round(np.arange(12, dtype=float).reshape(3, 4) * 1.25 + 1000 + k, 4)
        b = np.round(-a / 7.0, 4)
        out.append((t, 3, 4, 0.25, 0.5, 17.0, -99.0, a, b))
    return out


def test_pressure_grid_header_and_values(tmp_path):
    p = syn.write_owi_synth(tmp_path / "fort.221", _snaps(False))
    snaps = read_owi(p)
    assert len(snaps) == 2
    s = snaps[0]
    assert (s.ilat, s.ilon) == (3, 4)
    assert (s.dx, s.dy, s.swlat, s.swlon) == (0.25, 0.5, 17.0, -99.0)
    assert s.time == dt.datetime(2005, 8, 25, 0, 0)
    assert snaps[1].time == dt.datetime(2005, 8, 25, 6, 0)
    assert s.values.shape == (3, 4) and s.u is None and s.v is None
    assert s.values[0, 0] == 1000.0 and s.values[2, 3] == pytest.approx(1000 + 11 * 1.25)


def test_wind_u_v_blocks(tmp_path):
    src = _snaps(True)
    p = syn.write_owi_synth(tmp_path / "fort.222", src, wind=True)
    s = read_owi(p)[1]
    assert s.values is None
    np.testing.assert_allclose(s.u, src[1][7], atol=1e-4)
    np.testing.assert_allclose(s.v, src[1][8], atol=1e-4)


@pytest.mark.parametrize("wind", [False, True])
def test_round_trip(tmp_path, wind):
    p = syn.write_owi_synth(tmp_path / "in", _snaps(wind), wind=wind)
    first = read_owi(p)
    out = tmp_path / "out"
    write_owi(first, out)
    again = read_owi(out)
    assert len(again) == len(first)
    for x, y in zip(first, again):
        assert (x.time, x.ilat, x.ilon, x.dx, x.dy, x.swlat, x.swlon) == \
               (y.time, y.ilat, y.ilon, y.dx, y.dy, y.swlat, y.swlon)
        for name in ("values", "u", "v"):
            a, b = getattr(x, name), getattr(y, name)
            assert (a is None) == (b is None)
            if a is not None:
                np.testing.assert_array_equal(a, b)


def test_grid_not_multiple_of_eight(tmp_path):
    # 3 x 4 = 12 values: second line holds 4, so partial last line must parse
    p = syn.write_owi_synth(tmp_path / "f", _snaps(False))
    assert read_owi(p)[0].values.size == 12


def test_truncated_snapshot_error(tmp_path):
    p = syn.write_owi_synth(tmp_path / "f", _snaps(False))
    lines = p.read_text().splitlines()
    p.write_text("\n".join(lines[:-1]) + "\n")
    with pytest.raises(ValueError, match="values"):
        read_owi(p)


def test_not_owi_error(tmp_path):
    p = tmp_path / "f"
    p.write_text("hello\nworld\n")
    with pytest.raises(ValueError, match="OWI"):
        read_owi(p)
