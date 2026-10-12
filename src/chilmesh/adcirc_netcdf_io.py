"""NetCDF ADCIRC output readers (spec 005 US5). Needs the optional extra ``netcdf``."""
from __future__ import annotations

from typing import Optional

import numpy as np

from .adcirc_output_io import FILL, OutputSeries, _location

__all__ = ["read_adcirc_output_nc"]

_MAX_VARS = {"zeta_max": "time_of_zeta_max", "vel_max": "time_of_vel_max"}


def _import_netcdf4():
    try:
        import netCDF4
    except ImportError as exc:
        raise ImportError(
            "reading NetCDF ADCIRC output needs netCDF4; install the optional "
            "extra: pip install 'chilmesh[netcdf]'") from exc
    return netCDF4


def _clean(a) -> np.ndarray:
    a = np.ma.filled(np.ma.asarray(a, dtype=float), np.nan)
    a[a == FILL] = np.nan
    return a


def read_adcirc_output_nc(path, mesh=None, step: Optional[int] = None,
                          variable: Optional[str] = None) -> OutputSeries:
    """Read fort.61/62/63/64 .nc, maxele.63.nc or maxvel.63.nc output.

    Returns the same :class:`~chilmesh.adcirc_output_io.OutputSeries` as
    :func:`~chilmesh.read_adcirc_output`. ``step`` reads one record by slicing
    the variable, so other steps are not loaded. ``variable`` picks one named
    variable (returned as scalar); by default ``zeta``, ``u-vel`` + ``v-vel``
    (vector), ``zeta_max`` or ``vel_max`` (with the time-of-maximum as record 1)
    is chosen from what the file holds. NetCDF files store no iteration count,
    so ``iterations`` is -1. Needs the ``netcdf`` extra.

    Example:
        >>> from chilmesh import read_adcirc_output_nc  # doctest: +SKIP
        >>> series = read_adcirc_output_nc("fort.63.nc", step=0)  # doctest: +SKIP
    """
    nc4 = _import_netcdf4()
    location = _location(path)
    with nc4.Dataset(str(path)) as d:
        names = set(d.variables)
        if "station_name" in names:
            location = "stations"
        tvar = d.variables["time"]
        nt = len(tvar)

        def pick(name):
            return d.variables[name]

        if variable is not None:
            if variable not in names:
                raise ValueError(f"variable {variable!r} not in file; has {sorted(names)}")
            comps, kind = [variable], "scalar"
        elif "u-vel" in names and "v-vel" in names:
            comps, kind = ["u-vel", "v-vel"], "vector"
        else:
            comps = next(([m] for m in ("zeta", "zeta_max", "vel_max") if m in names), None)
            if comps is None:
                raise ValueError(f"no ADCIRC output variable found in {path}; has {sorted(names)}")
            kind = "scalar"
        first = comps[0]
        is_max = first in _MAX_VARS and _MAX_VARS[first] in names
        if is_max:
            arrs = [_clean(pick(first)[:]), _clean(pick(_MAX_VARS[first])[:])]
            values = np.stack(arrs)
            tval = float(np.asarray(tvar[:]).ravel()[0])
            times = np.array([tval, tval])
        else:
            idx = slice(None)
            if step is not None:
                if not -nt <= step < nt:
                    raise IndexError(f"step {step} out of range for {nt} record(s)")
                idx = step % nt
            cols = [_clean(pick(c)[idx, :]) for c in comps]
            values = np.stack(cols, axis=-1) if len(cols) > 1 else cols[0]
            times = np.atleast_1d(np.asarray(tvar[idx], dtype=float))
            if step is not None:
                values = values[None]
        n = values.shape[1]
        if mesh is not None and location == "nodes" and n != mesh.n_verts:
            raise ValueError(f"file has {n} nodes but mesh has {mesh.n_verts} nodes")
    if is_max and step is not None:
        if not -2 <= step < 2:
            raise IndexError(f"step {step} out of range for 2 record(s)")
        values, times = values[[step % 2]], times[[step % 2]]
    return OutputSeries(times, np.full(len(times), -1, dtype=int), values, kind, location, 0)
