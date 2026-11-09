"""Direct NOAA fetch for auxiliary data only (plan D1): ``points.cwl.nc`` and the CO-OPS station list.

Meshes never come from here; they come from Valence. Every fetch records the
URL, SHA-256 and fetch date.
"""
from __future__ import annotations

import datetime
import hashlib
import urllib.request
from pathlib import Path
from typing import Optional

DEFAULT_CHUNK = 1 << 20

# Names the fetched auxiliary files carry inside the aux directory (``--aux-dir``).
STATIONS_FILE = "coops_stations.json"  # CO-OPS metadata API station list: {"stations": [{"lat", "lng"}, ...]}
POINTS_FILE = "points.cwl.nc"
NOAA_INPUTS_REASON = "NOAA inputs not fetched"
# STOFS max-elevation field: one value per mesh node, the quantity ``points.cwl.nc`` peak zeta is compared with.
MAXELE_FILE = "maxele.63.nc"
MAXELE_VAR = "zeta_max"
MAXELE_REASON = "NOAA inputs not fetched (maxele)"


class MaxeleError(RuntimeError):
    """The max-elevation file is unreadable or lacks the expected variable; the message is the failure reason."""


def aux_paths(aux_dir):
    """``(stations_path, points_path)`` inside ``aux_dir``; each is None when the file is absent."""
    if not aux_dir:
        return None, None
    base = Path(aux_dir)
    found = [base / name for name in (STATIONS_FILE, POINTS_FILE)]
    return tuple(p if p.is_file() else None for p in found)


def load_station_coords(path):
    """``(n, 2)`` array of ``(lng, lat)`` from a CO-OPS station list, or None when the file is absent.

    Rows without numeric ``lat`` and ``lng`` are skipped.
    """
    import json

    import numpy as np

    path = Path(path)
    if not path.is_file():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = data.get("stations", []) if isinstance(data, dict) else data
    xy = []
    for row in rows:
        try:
            xy.append((float(row["lng"]), float(row["lat"])))
        except (KeyError, TypeError, ValueError):
            continue
    return np.asarray(xy, dtype=float).reshape(-1, 2)


def maxele_path(aux_dir):
    """Path of the max-elevation file inside ``aux_dir``, or None when it is absent."""
    if not aux_dir:
        return None
    p = Path(aux_dir) / MAXELE_FILE
    return p if p.is_file() else None


def reduce_time_max(zeta):
    """Max over time (axis 0) of a ``(time, station)`` array; a 1-D array is returned unchanged."""
    import numpy as np

    zeta = np.asarray(zeta, dtype=float)
    return zeta.max(axis=0) if zeta.ndim > 1 else zeta


def _read_cwl_xarray(path):
    """Default ``points.cwl.nc`` reader. Assumes variables ``x``, ``y`` and ``zeta``; value is max over time."""
    import xarray as xr

    with xr.open_dataset(path) as ds:
        xy = list(zip(ds["x"].values.ravel(), ds["y"].values.ravel()))
        vals = reduce_time_max(ds["zeta"].values)
    return xy, vals


def _read_maxele_xarray(path):
    import xarray as xr

    try:
        with xr.open_dataset(path) as ds:
            if MAXELE_VAR not in ds.variables:
                raise MaxeleError(f"variable '{MAXELE_VAR}' not found in {Path(path).name}")
            vals = ds[MAXELE_VAR].values
    except MaxeleError:
        raise
    except Exception as exc:  # unreadable or corrupt file
        raise MaxeleError(f"cannot read {Path(path).name}: {exc}") from exc
    return vals


def read_maxele(path, *, reader=None):
    """Nodal max-elevation array from the STOFS maxele file. ``reader(path) -> values`` is injectable.

    Raises ``MaxeleError`` when the ``zeta_max`` variable is missing or the file cannot be read.
    """
    import numpy as np

    return np.asarray((reader or _read_maxele_xarray)(path), dtype=float).ravel()


def read_points_cwl(path, *, reader=None):
    """``(xy, values)`` from ``points.cwl.nc``. ``reader(path) -> (xy, values)`` is injectable for tests."""
    import numpy as np

    xy, vals = (reader or _read_cwl_xarray)(path)
    return np.asarray(xy, dtype=float).reshape(-1, 2), np.asarray(vals, dtype=float).ravel()


def fetch_aux(url: str, dest, *, opener=None, today: Optional[str] = None,
              chunk_bytes: int = DEFAULT_CHUNK) -> dict:
    opener = opener or urllib.request.urlopen
    today = today or datetime.date.today().isoformat()
    if not url.lower().startswith("https://"):
        return {"status": "failed", "url": url, "fetch_date": today, "reason": "only https URLs are fetched"}
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    try:
        with opener(urllib.request.Request(url, headers={"User-Agent": "chilmesh-bench/1"}), timeout=120) as resp, \
                open(dest, "wb") as fh:
            while True:
                chunk = resp.read(chunk_bytes)
                if not chunk:
                    break
                digest.update(chunk)
                fh.write(chunk)
    except OSError as exc:
        return {"status": "failed", "url": url, "fetch_date": today, "reason": f"fetch failed: {exc}"}
    return {"status": "ok", "url": url, "sha256": digest.hexdigest(), "fetch_date": today, "path": str(dest)}
