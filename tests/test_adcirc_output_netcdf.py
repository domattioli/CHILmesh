"""NetCDF ADCIRC output readers (spec 005 US5)."""
import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent / "fixtures"))
import synthetic_adcirc as syn  # noqa: E402

from chilmesh import examples, read_adcirc_output, read_adcirc_output_nc  # noqa: E402

N = 7
TIMES = [10.0, 20.0, 30.0]


def _scalar(nt=3, n=N):
    v = np.arange(nt * n, dtype=float).reshape(nt, n) / 10.0 - 0.5
    v[0, 2] = np.nan
    return v


def test_import_chilmesh_without_netcdf4():
    code = ("import sys; sys.modules['netCDF4'] = None; import chilmesh; "
            "f = chilmesh.read_adcirc_output_nc; print('ok')")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0 and "ok" in r.stdout, r.stderr


def test_missing_netcdf4_names_extra(monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, "netCDF4", None)
    with pytest.raises(ImportError, match=r"chilmesh\[netcdf\]"):
        read_adcirc_output_nc(tmp_path / "fort.63.nc")


needs_nc = pytest.mark.skipif(importlib.util.find_spec("netCDF4") is None,
                           reason="netCDF4 (extra: netcdf) not installed")


@needs_nc
def test_scalar_same_as_ascii(tmp_path):
    vals = _scalar()
    nc = read_adcirc_output_nc(syn.write_output_nc(tmp_path / "fort.63.nc", TIMES, vals))
    asc = read_adcirc_output(syn.write_output_ascii(tmp_path / "fort.63", TIMES, vals))
    assert nc.kind == asc.kind == "scalar" and nc.location == "nodes"
    np.testing.assert_allclose(nc.times, asc.times)
    np.testing.assert_allclose(nc.values, asc.values, atol=1e-9, equal_nan=True)


@needs_nc
def test_vector_same_as_ascii(tmp_path):
    vals = np.arange(3 * N * 2, dtype=float).reshape(3, N, 2) / 7.0
    vals[1, 4, :] = np.nan
    nc = read_adcirc_output_nc(syn.write_output_nc(tmp_path / "fort.64.nc", TIMES, vals))
    asc = read_adcirc_output(syn.write_output_ascii(tmp_path / "fort.64", TIMES, vals))
    assert nc.kind == "vector" and nc.values.shape == (3, N, 2)
    np.testing.assert_allclose(nc.values, asc.values, atol=1e-9, equal_nan=True)


@needs_nc
def test_station_files(tmp_path):
    vals = _scalar(nt=4, n=3)
    nc = read_adcirc_output_nc(
        syn.write_output_nc(tmp_path / "fort.61.nc", [1, 2, 3, 4], vals, stations=True))
    assert nc.location == "stations" and nc.values.shape == (4, 3)
    v2 = np.ones((2, 3, 2))
    nc2 = read_adcirc_output_nc(
        syn.write_output_nc(tmp_path / "fort.62.nc", [1, 2], v2, stations=True))
    assert nc2.location == "stations" and nc2.kind == "vector"


@needs_nc
@pytest.mark.parametrize("var,name", [("zeta_max", "maxele.63.nc"), ("vel_max", "maxvel.63.nc")])
def test_max_files_match_ascii(tmp_path, var, name):
    mx, tom = np.linspace(0.1, 2.0, N), np.linspace(100.0, 900.0, N)
    mx[3] = np.nan
    tom[3] = np.nan
    nc = read_adcirc_output_nc(
        syn.write_output_nc(tmp_path / name, [500.0], np.vstack([mx, tom]), variable=var))
    asc = read_adcirc_output(syn.write_max_ascii(tmp_path / name[:-3], mx, tom, t_end=500.0))
    assert nc.values.shape == (2, N)
    np.testing.assert_allclose(nc.values, asc.values, atol=1e-6, equal_nan=True)
    np.testing.assert_allclose(nc.times, asc.times)


@needs_nc
def test_single_step_is_lazy(tmp_path, monkeypatch):
    import netCDF4
    vals = _scalar()
    p = syn.write_output_nc(tmp_path / "fort.63.nc", TIMES, vals)
    seen = []
    orig_dataset = netCDF4.Dataset

    # netCDF4.Variable is an immutable extension type in recent releases, so
    # spy on reads through a proxy Dataset instead of patching the class.
    class SpyVar:
        def __init__(self, var):
            self._var = var

        def __getitem__(self, key):
            seen.append((self._var.name, key))
            return self._var[key]

        def __len__(self):
            return len(self._var)

        def __getattr__(self, attr):
            return getattr(self._var, attr)

    class SpyDataset:
        def __init__(self, *args, **kwargs):
            self._ds = orig_dataset(*args, **kwargs)
            self.variables = {k: SpyVar(v) for k, v in self._ds.variables.items()}

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            self._ds.close()

        def __getattr__(self, attr):
            return getattr(self._ds, attr)

    monkeypatch.setattr(netCDF4, "Dataset", SpyDataset)
    s = read_adcirc_output_nc(p, step=1)
    assert s.values.shape == (1, N)
    np.testing.assert_allclose(s.values[0], vals[1], atol=1e-9, equal_nan=True)
    np.testing.assert_allclose(s.times, [20.0])
    zeta_reads = [k for n, k in seen if n == "zeta"]
    assert zeta_reads and all(k != slice(None) and k != (slice(None), slice(None)) for k in zeta_reads)
    with pytest.raises(IndexError):
        read_adcirc_output_nc(p, step=3)


@needs_nc
def test_node_count_mismatch(tmp_path):
    mesh = examples.annulus.__wrapped__()
    p = syn.write_output_nc(tmp_path / "fort.63.nc", TIMES, _scalar())
    with pytest.raises(ValueError) as ei:
        read_adcirc_output_nc(p, mesh=mesh)
    assert str(N) in str(ei.value) and str(mesh.n_verts) in str(ei.value)


@needs_nc
def test_explicit_variable(tmp_path):
    vals = np.arange(3 * N * 2, dtype=float).reshape(3, N, 2)
    p = syn.write_output_nc(tmp_path / "fort.64.nc", TIMES, vals)
    s = read_adcirc_output_nc(p, variable="v-vel")
    assert s.kind == "scalar"
    np.testing.assert_allclose(s.values, vals[:, :, 1])
    with pytest.raises(ValueError, match="nope"):
        read_adcirc_output_nc(p, variable="nope")
