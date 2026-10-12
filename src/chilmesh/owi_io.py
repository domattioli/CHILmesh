"""Oceanweather (OWI) WIN/PRE grid files: ADCIRC fort.221 and fort.222 (spec 005 US7).

fort.221 holds pressure and fort.222 holds wind. Each snapshot has one header
line followed by ``ilat * ilon`` values in ``F10.4``, eight per line, with
longitude varying fastest. A wind snapshot carries the U block and then the V
block under a single header.
"""
from __future__ import annotations

import datetime as _dt
import re
from dataclasses import dataclass
from typing import List, Optional

import numpy as np

__all__ = ["OWISnapshot", "read_owi", "write_owi"]

_HEAD = re.compile(
    r"iLat=\s*(\d+)iLong=\s*(\d+)DX=\s*([-\d.]+)DY=\s*([-\d.]+)"
    r"SWLat=\s*([-\d.]+)SWLon=\s*([-\d.]+)DT=\s*(\d{12})")
_FMT_TITLE = "Oceanweather WIN/PRE Format"
_W = 10  # field width


@dataclass
class OWISnapshot:
    """One time level of an OWI grid.

    Pressure files fill ``values`` (shape ``(ilat, ilon)``); wind files fill
    ``u`` and ``v`` of the same shape. Unused arrays are ``None``.
    """

    time: _dt.datetime
    ilat: int
    ilon: int
    dx: float
    dy: float
    swlat: float
    swlon: float
    values: Optional[np.ndarray] = None
    u: Optional[np.ndarray] = None
    v: Optional[np.ndarray] = None


def _read_block(lines, pos: int, n: int, path):
    """Read ``n`` fixed-width values starting at ``pos``; return (array, new pos)."""
    nlines = -(-n // 8)
    if pos + nlines > len(lines):
        raise ValueError(f"{path}: truncated OWI snapshot, expected {n} values")
    text = "".join(ln.rstrip("\r") for ln in lines[pos:pos + nlines])
    try:
        arr = np.array([float(text[i:i + _W]) for i in range(0, n * _W, _W)])
    except ValueError:
        raise ValueError(f"{path}: truncated or malformed OWI snapshot, expected {n} values") from None
    return arr, pos + nlines


def read_owi(path) -> List[OWISnapshot]:
    """Read an OWI pressure (fort.221) or wind (fort.222) file.

    Returns one :class:`OWISnapshot` per time level. A snapshot followed by a
    second value block (before the next header) is a wind snapshot.

    Example
    -------
    >>> snaps = read_owi("fort.221")  # doctest: +SKIP
    >>> snaps[0].values.shape  # doctest: +SKIP
    (100, 100)

    Tested by: tests/test_owi.py::test_wind_u_v_blocks.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.read().split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines or not lines[0].lstrip().startswith("Oceanweather"):
        raise ValueError(f"{path}: not an OWI WIN/PRE file (first line must start with 'Oceanweather')")
    out, pos = [], 1
    while pos < len(lines):
        if not lines[pos].strip():
            pos += 1
            continue
        m = _HEAD.match(lines[pos])
        if not m:
            raise ValueError(f"{path}:{pos + 1}: expected an OWI grid header (iLat=...)")
        ilat, ilon = int(m.group(1)), int(m.group(2))
        dx, dy, swlat, swlon = (float(m.group(i)) for i in (3, 4, 5, 6))
        t = m.group(7)
        time = _dt.datetime(int(t[:4]), int(t[4:6]), int(t[6:8]), int(t[8:10]), int(t[10:12]))
        pos += 1
        n = ilat * ilon
        a, pos = _read_block(lines, pos, n, path)
        a = a.reshape(ilat, ilon)
        if pos < len(lines) and lines[pos].strip() and not lines[pos].startswith("iLat="):
            b, pos = _read_block(lines, pos, n, path)
            out.append(OWISnapshot(time, ilat, ilon, dx, dy, swlat, swlon, u=a, v=b.reshape(ilat, ilon)))
        else:
            out.append(OWISnapshot(time, ilat, ilon, dx, dy, swlat, swlon, values=a))
    return out


def _write_block(lines, arr, shape):
    arr = np.asarray(arr, dtype=float)
    if arr.shape != shape:
        raise ValueError(f"array shape {arr.shape} does not match grid {shape}")
    flat = arr.ravel()
    for i in range(0, flat.size, 8):
        lines.append("".join(f"{x:{_W}.4f}" for x in flat[i:i + 8]))


def write_owi(snapshots: List[OWISnapshot], path) -> None:
    """Write snapshots as an OWI WIN/PRE file (fort.221 or fort.222).

    A snapshot with ``values`` is written as pressure; one with ``u`` and
    ``v`` as wind. All snapshots in one file must be of the same kind.

    Example
    -------
    >>> write_owi(snaps, "fort.221")  # doctest: +SKIP

    Tested by: tests/test_owi.py::test_round_trip.
    """
    if not snapshots:
        raise ValueError("no snapshots to write")
    kinds = {s.values is not None for s in snapshots}
    if len(kinds) != 1:
        raise ValueError("cannot mix pressure and wind snapshots in one file")
    if any(s.values is None and (s.u is None or s.v is None) for s in snapshots):
        raise ValueError("wind snapshots need both u and v")
    fmt = "%Y%m%d%H"
    lines = [f"{_FMT_TITLE:<59}{snapshots[0].time.strftime(fmt)}     {snapshots[-1].time.strftime(fmt)}"]
    for s in snapshots:
        lines.append(f"iLat={s.ilat:4d}iLong={s.ilon:4d}DX={s.dx:6.4f}DY={s.dy:6.4f}"
                     f"SWLat={s.swlat:8.5f}SWLon={s.swlon:8.4f}DT={s.time.strftime('%Y%m%d%H%M')}")
        shape = (s.ilat, s.ilon)
        if s.values is not None:
            _write_block(lines, s.values, shape)
        else:
            _write_block(lines, s.u, shape)
            _write_block(lines, s.v, shape)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
