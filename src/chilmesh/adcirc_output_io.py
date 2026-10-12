"""Readers for ADCIRC ASCII output files (spec 005 US4).

Covers the global node files ``fort.63`` (scalar) and ``fort.64`` (vector),
the station files ``fort.61`` and ``fort.62``, and the two-data-set
``maxele.63`` / ``maxvel.63`` files. Full and sparse record forms are read by
one parser.
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from typing import Optional

import numpy as np

__all__ = ["OutputSeries", "read_adcirc_output"]

FILL = -99999.0
_CHUNK = 20000  # record lines parsed per slab; bounds single-step memory
_STATION_RE = re.compile(r"fort\.6[12](\.nc)?$")


@dataclass
class OutputSeries:
    """Time series read from an ADCIRC output file.

    ``values`` is ``(nt, n)`` for scalar or ``(nt, n, 2)`` for vector data;
    NaN marks dry or fill points. For maxele/maxvel files record 1 is the time
    of the maximum. ``iterations`` is -1 for NetCDF files (not stored there).
    """

    times: np.ndarray
    iterations: np.ndarray
    values: np.ndarray
    kind: str
    location: str
    dropped_records: int = 0


def _location(path) -> str:
    return "stations" if _STATION_RE.search(str(path)) else "nodes"


def _read_header(f):
    f.readline()  # free-text title
    tok = f.readline().split()
    if len(tok) < 5:
        raise ValueError("not an ADCIRC output file: second line needs NDSETS NP DT NSPOOL IRTYPE")
    ndsets, n, irtype = int(tok[0]), int(tok[1]), int(tok[4])
    if irtype not in (1, 2):
        raise ValueError(f"unsupported IRTYPE {irtype} (only 1 = scalar, 2 = vector)")
    return ndsets, n, irtype


def _scan(f, ndsets: int, n: int, ncomp: int):
    """Index data sets: returns offsets, header tokens, n complete sets."""
    offsets, heads = [], []
    for _ in range(ndsets):
        off = f.tell()
        line = f.readline()
        if not line:
            break
        tok = line.split()
        if len(tok) < 2:
            break
        nrec = int(tok[2]) if len(tok) >= 4 else n
        ok = True
        for k in range(nrec):
            rec = f.readline()
            if not rec:
                ok = False
                break
            if k == nrec - 1 and not rec.endswith(b"\n") and len(rec.split()) != 1 + ncomp:
                ok = False  # last line cut mid-record
        if not ok:
            break
        offsets.append(off)
        heads.append(tok)
    return offsets, heads


def _read_set(f, off: int, head, n: int, ncomp: int) -> np.ndarray:
    """Read one data set into an (n, ncomp) array (NaN where fill or unlisted)."""
    sparse = len(head) >= 4
    nrec = int(head[2]) if sparse else n
    defval = float(head[3]) if sparse else FILL
    out = np.full((n, ncomp), np.nan)
    f.seek(off)
    f.readline()
    left = nrec
    while left > 0:
        lines = [f.readline() for _ in range(min(_CHUNK, left))]
        left -= len(lines)
        a = np.array(b"".join(lines).split(), dtype=float).reshape(-1, 1 + ncomp)
        ids = a[:, 0].astype(np.int64) - 1
        vals = a[:, 1:]
        vals[(vals == FILL) | (vals == defval)] = np.nan
        out[ids] = vals
    return out


def read_adcirc_output(path, mesh=None, step: Optional[int] = None) -> OutputSeries:
    """Read an ADCIRC ASCII output file (fort.61/62/63/64, maxele.63, maxvel.63).

    ``step`` selects one record by index (negative counts from the end) and
    reads only that record; the result then has ``nt == 1``. With ``mesh``
    given, node files must match ``mesh.n_verts`` or ``ValueError`` names both
    counts. A file that ends early yields its complete records, a
    ``UserWarning`` and ``dropped_records``. Fill value -99999 becomes NaN, as
    do nodes unlisted in the sparse form.

    Example:
        >>> import tempfile, os
        >>> p = os.path.join(tempfile.mkdtemp(), "fort.63")
        >>> _ = open(p, "w").write("t\\n 1 2 1.0 1 1\\n 5.0 1\\n 1 0.5\\n 2 -99999\\n")
        >>> s = read_adcirc_output(p)
        >>> s.values.tolist()
        [[0.5, nan]]

    Tested by: tests/test_adcirc_output_ascii.py::test_full_scalar_values.
    """
    location = _location(path)
    with open(path, "rb") as f:
        ndsets, n, irtype = _read_header(f)
        if mesh is not None and location == "nodes" and n != mesh.n_verts:
            raise ValueError(
                f"file has {n} nodes but mesh has {mesh.n_verts} nodes")
        ncomp = irtype
        offsets, heads = _scan(f, ndsets, n, ncomp)
        dropped = ndsets - len(offsets)
        if dropped:
            warnings.warn(
                f"{path}: file ends early; {dropped} record(s) dropped, "
                f"{len(offsets)} complete record(s) returned", UserWarning, stacklevel=2)
        idx = range(len(offsets))
        if step is not None:
            if not -len(offsets) <= step < len(offsets):
                raise IndexError(f"step {step} out of range for {len(offsets)} complete record(s)")
            idx = [step % len(offsets)]
        vals = np.stack([_read_set(f, offsets[i], heads[i], n, ncomp) for i in idx]) \
            if len(idx) else np.empty((0, n, ncomp))
    if irtype == 1:
        vals = vals[..., 0]
    times = np.array([float(heads[i][0]) for i in idx])
    its = np.array([int(float(heads[i][1])) for i in idx], dtype=int)
    return OutputSeries(times, its, vals, "scalar" if irtype == 1 else "vector",
                        location, dropped)
