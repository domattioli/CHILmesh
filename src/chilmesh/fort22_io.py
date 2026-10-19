"""Best-track fort.22 reader and writer for ADCIRC NWS 8 and NWS 20 (spec 005 US7).

A best-track fort.22 is an ATCF b-deck: comma separated columns, one record
per line, several records per time (one per wind-radius threshold). Only the
columns CHILmesh interprets are typed; every other column is kept as text in
record order and written back unchanged, so the writer needs no per-version
column list (NWS 20 simply has more untyped columns than NWS 8).
"""
from __future__ import annotations

import datetime as _dt
from dataclasses import dataclass, field
from typing import List, Optional

__all__ = ["StormRecord", "StormTrack", "read_fort22_besttrack", "write_fort22_besttrack"]

SUPPORTED_NWS = (8, 20)

# ATCF column positions of the typed fields.
_C_TIME, _C_LAT, _C_LON, _C_VMAX, _C_PMIN, _C_ISO = 2, 6, 7, 8, 9, 11
_C_RADII = (13, 14, 15, 16)
_C_RMW = 19
_TYPED = {_C_TIME, _C_LAT, _C_LON, _C_VMAX, _C_PMIN, _C_ISO, *_C_RADII, _C_RMW}


@dataclass
class StormRecord:
    """One best-track line. ``None`` marks a blank or absent field.

    ``extra`` holds the untyped columns (basin, cyclone number, technique,
    wind code, storm name, NWS 20 extras, ...) as stripped text in column
    order; ``_ncols`` (private) is the original column count so the writer reproduces
    short lines without padding them.
    """

    time: Optional[_dt.datetime]
    lat: Optional[float]
    lon: Optional[float]
    vmax: Optional[float]
    pmin: Optional[float]
    isotach: Optional[float]
    radii: List[Optional[float]]
    rmw: Optional[float]
    extra: List[str] = field(default_factory=list)
    _ncols: int = field(default=0, repr=False)


@dataclass
class StormTrack:
    """Ordered list of :class:`StormRecord` read from a best-track fort.22."""

    records: List[StormRecord]


def _check_nws(nws) -> int:
    try:
        n = abs(int(nws))
    except (TypeError, ValueError):
        n = None
    if n not in SUPPORTED_NWS:
        raise ValueError(
            f"unsupported NWS {nws!r}: best-track fort.22 is supported for NWS 8 and 20 only")
    return n


def _num(s: str) -> Optional[float]:
    return float(s) if s else None


def _parse_latlon(s: str, pos: str, neg: str) -> Optional[float]:
    if not s:
        return None
    hemi = s[-1].upper()
    if hemi not in (pos, neg):
        raise ValueError(f"bad coordinate {s!r}")
    v = int(s[:-1]) / 10.0
    return v if hemi == pos else -v


def _fmt_latlon(v: Optional[float], pos: str, neg: str) -> str:
    if v is None:
        return ""
    return f"{round(abs(v) * 10):d}{pos if v >= 0 else neg}"


def _fmt_num(v: Optional[float]) -> str:
    if v is None:
        return ""
    return str(int(v)) if float(v) == int(v) else repr(float(v))


def _parse_line(line: str) -> StormRecord:
    cols = [c.strip() for c in line.rstrip("\r\n").split(",")]
    if len(cols) < 8 or len(cols[_C_TIME]) != 10 or not cols[_C_TIME].isdigit():
        raise ValueError("not a best-track fort.22 line (expected ATCF columns with a YYYYMMDDHH time)")

    def col(i):
        return cols[i] if i < len(cols) else ""

    t = col(_C_TIME)
    time = _dt.datetime(int(t[:4]), int(t[4:6]), int(t[6:8]), int(t[8:10])) if t else None
    return StormRecord(
        time=time,
        lat=_parse_latlon(col(_C_LAT), "N", "S"),
        lon=_parse_latlon(col(_C_LON), "E", "W"),
        vmax=_num(col(_C_VMAX)),
        pmin=_num(col(_C_PMIN)),
        isotach=_num(col(_C_ISO)),
        radii=[_num(col(i)) for i in _C_RADII],
        rmw=_num(col(_C_RMW)),
        extra=[c for i, c in enumerate(cols) if i not in _TYPED],
        _ncols=len(cols),
    )


def read_fort22_besttrack(path, nws: int = 8) -> StormTrack:
    """Read a best-track fort.22 (NWS 8 or NWS 20).

    Typed fields are time, latitude and longitude (degrees, west negative),
    maximum wind, central pressure, isotach and its four wind radii, and the
    radius of maximum wind. Blank or missing fields read as ``None``. All
    other columns stay in ``StormRecord.extra``.

    Raises ``ValueError`` for an unsupported ``nws`` or a non-ATCF file.

    Example
    -------
    >>> import os, tempfile
    >>> p = os.path.join(tempfile.mkdtemp(), "fort.22")
    >>> _ = open(p, "w").write("AL, 09, 2011082100,   , BEST,   0, 150N,  590W,  45, 1006, TS,"
    ...                        "  34, NEQ,   60,    0,    0,   60, 1012,  150,  40,\\n")
    >>> track = read_fort22_besttrack(p)
    >>> track.records[0].vmax, track.records[0].lon
    (45.0, -59.0)

    Tested by: tests/test_fort22.py::test_typed_fields_from_real_atcf_file.
    """
    _check_nws(nws)
    records = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for lineno, line in enumerate(f, 1):
            if not line.strip():
                continue
            try:
                records.append(_parse_line(line))
            except ValueError as exc:
                raise ValueError(f"{path}:{lineno}: {exc}") from None
    if not records:
        raise ValueError(f"{path}: no best-track records found")
    return StormTrack(records)


def write_fort22_besttrack(track: StormTrack, path, nws: int = 8) -> None:
    """Write ``track`` as a best-track fort.22. Untyped columns are written unchanged.

    Raises ``ValueError`` for an unsupported ``nws``.

    Example
    -------
    >>> import os, tempfile
    >>> d = tempfile.mkdtemp()
    >>> line = "AL, 09, 2011082100,   , BEST,   0, 150N,  590W,  45, 1006, TS,  34\\n"
    >>> _ = open(os.path.join(d, "in.22"), "w").write(line)
    >>> track = read_fort22_besttrack(os.path.join(d, "in.22"))
    >>> write_fort22_besttrack(track, os.path.join(d, "fort.22"), nws=8)
    >>> read_fort22_besttrack(os.path.join(d, "fort.22")).records[0].pmin
    1006.0

    Tested by: tests/test_fort22.py::test_round_trip_real_file.
    """
    _check_nws(nws)
    lines = []
    for r in track.records:
        typed = {
            _C_TIME: r.time.strftime("%Y%m%d%H") if r.time else "",
            _C_LAT: _fmt_latlon(r.lat, "N", "S"),
            _C_LON: _fmt_latlon(r.lon, "E", "W"),
            _C_VMAX: _fmt_num(r.vmax),
            _C_PMIN: _fmt_num(r.pmin),
            _C_ISO: _fmt_num(r.isotach),
            _C_RMW: _fmt_num(r.rmw),
        }
        typed.update({c: _fmt_num(v) for c, v in zip(_C_RADII, r.radii)})
        extra = iter(r.extra)
        cols = [typed[i] if i in _TYPED else next(extra) for i in range(r._ncols)]
        lines.append(", ".join(cols))
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
