"""Synthetic ADCIRC-format writers for tests (skeleton; extended in Stage B2/B3)."""
from __future__ import annotations

from pathlib import Path


def write_text(path: Path, lines: list[str]) -> Path:
    """Write ``lines`` to ``path`` with LF endings and return the path."""
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# --------------------------------------------------------------------------
# Stage B2: ADCIRC output files (fort.61/62/63/64, maxele.63, maxvel.63, .nc)
# --------------------------------------------------------------------------
import numpy as np

FILL = -99999.0


def _fmt_rec(i: int, row) -> str:
    return f"{i:10d}" + "".join(f" {v:20.10E}" for v in row)


def write_output_ascii(path, times, values, *, sparse: bool = False,
                       iterations=None, dt: float = 1.0, nspool: int = 1,
                       defval: float = FILL, header: str = "synthetic ADCIRC output") -> Path:
    """Write a global or station output file in ADCIRC ASCII layout.

    ``values`` has shape ``(nt, n)`` (scalar, IRTYPE 1) or ``(nt, n, 2)``
    (vector, IRTYPE 2). Station files use the same layout (``n`` = stations).
    NaN or ``defval`` entries are written as ``defval`` in full form and
    omitted in sparse form (header ``TIME IT NNONDEF DEFVAL``).
    """
    values = np.asarray(values, dtype=float)
    nt, n = values.shape[:2]
    irtype = 2 if values.ndim == 3 else 1
    if iterations is None:
        iterations = [(k + 1) * nspool for k in range(nt)]
    out = [header, f"{nt:10d} {n:10d} {dt:16.7E} {nspool:5d} {irtype:5d} FileFmtVersion: 1050624"]
    for k in range(nt):
        v = values[k].reshape(n, -1)
        bad = np.isnan(v).any(axis=1) | (v == defval).any(axis=1)
        if sparse:
            idx = np.nonzero(~bad)[0]
            out.append(f"{float(times[k]):20.10E} {int(iterations[k]):10d} {len(idx):10d} {defval:16.7E}")
        else:
            idx = np.arange(n)
            out.append(f"{float(times[k]):20.10E} {int(iterations[k]):10d}")
        for i in idx:
            row = np.full(v.shape[1], defval) if bad[i] else v[i]
            out.append(_fmt_rec(int(i) + 1, row))
    return write_text(Path(path), out)


def write_max_ascii(path, max_values, time_of_max, *, t_end: float = 100.0) -> Path:
    """Write a maxele.63 / maxvel.63 style file: two data sets (max, time of max)."""
    vals = np.vstack([np.asarray(max_values, float), np.asarray(time_of_max, float)])
    return write_output_ascii(path, [t_end, t_end], vals, iterations=[7, 7], dt=1.0, nspool=1)


def truncate_lines(path, drop: int) -> Path:
    """Remove the last ``drop`` lines of a text file (simulates an interrupted run)."""
    p = Path(path)
    lines = p.read_text().splitlines()
    p.write_text("\n".join(lines[: len(lines) - drop]) + "\n")
    return p


def write_output_nc(path, times, values, *, variable: str = "zeta", stations: bool = False) -> Path:
    """Write a minimal NetCDF output file (needs netCDF4; imported here, not at module load).

    ``values`` shape ``(nt, n)`` for ``variable='zeta'``, or ``(nt, n, 2)`` for
    ``u-vel``/``v-vel``. For max files pass ``variable='zeta_max'`` (or
    ``'vel_max'``) with ``values`` shape ``(2, n)``: row 0 the maximum, row 1 the
    time of maximum. NaN is stored as the ``_FillValue`` -99999.
    """
    import netCDF4
    values = np.asarray(values, dtype=float)
    dim = "station" if stations else "node"
    with netCDF4.Dataset(str(path), "w") as d:
        n = values.shape[1] if values.ndim >= 2 else values.shape[0]
        d.createDimension(dim, n)
        if variable in ("zeta_max", "vel_max"):
            d.createDimension("time", 1)
            t = d.createVariable("time", "f8", ("time",))
            t[:] = [times[0]]
            for name, row in ((variable, values[0]), ("time_of_" + variable, values[1])):
                v = d.createVariable(name, "f8", (dim,), fill_value=FILL)
                v[:] = np.ma.masked_invalid(row)
        else:
            d.createDimension("time", values.shape[0])
            t = d.createVariable("time", "f8", ("time",))
            t[:] = times
            if values.ndim == 3:
                for name, comp in (("u-vel", 0), ("v-vel", 1)):
                    v = d.createVariable(name, "f8", ("time", dim), fill_value=FILL)
                    v[:] = np.ma.masked_invalid(values[:, :, comp])
            else:
                v = d.createVariable(variable, "f8", ("time", dim), fill_value=FILL)
                v[:] = np.ma.masked_invalid(values)
        if stations:
            d.createDimension("namelen", 8)
            d.createVariable("station_name", "S1", (dim, "namelen"))
    return Path(path)
