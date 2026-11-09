"""W2 wiring: station sampling (T022, T023) and scaling collection (T025).

No competitor tools, network or netCDF library: station lists are tiny JSON files,
the points.cwl.nc reader is injected, and worker subprocesses are monkeypatched.
"""
from __future__ import annotations

import json
import math

import numpy as np
import pytest

from benchmarks.competitors import cli, noaa, records
from benchmarks.competitors.w2 import measures, run


# ---------------------------------------------------------------- noaa helpers

def _write_stations(path, rows):
    path.write_text(json.dumps({"stations": rows}), encoding="utf-8")
    return path


def test_load_station_coords_reads_lng_lat_and_skips_bad_rows(tmp_path):
    p = _write_stations(tmp_path / noaa.STATIONS_FILE, [
        {"id": "1", "lat": 40.0, "lng": -74.0},
        {"id": "2", "lat": "x", "lng": -73.0},
        {"id": "3"},
        {"id": "4", "lat": 41.5, "lng": -70.25},
    ])
    xy = noaa.load_station_coords(p)
    assert xy.shape == (2, 2)
    assert xy.tolist() == [[-74.0, 40.0], [-70.25, 41.5]]


def test_load_station_coords_missing_file_is_none(tmp_path):
    assert noaa.load_station_coords(tmp_path / "nope.json") is None


def test_aux_paths_needs_both_files(tmp_path):
    assert noaa.aux_paths(None) == (None, None)
    assert noaa.aux_paths(tmp_path) == (None, None)
    _write_stations(tmp_path / noaa.STATIONS_FILE, [])
    st, pt = noaa.aux_paths(tmp_path)
    assert st is not None and pt is None
    (tmp_path / noaa.POINTS_FILE).write_bytes(b"")
    st, pt = noaa.aux_paths(tmp_path)
    assert st is not None and pt is not None


def test_read_points_cwl_uses_injected_reader(tmp_path):
    f = tmp_path / noaa.POINTS_FILE
    f.write_bytes(b"")
    xy, vals = noaa.read_points_cwl(f, reader=lambda p: ([[1, 2], [3, 4]], [0.5, 0.25]))
    assert xy.shape == (2, 2) and vals.tolist() == [0.5, 0.25]


# ---------------------------------------------------------------- station sampling

def _flat_mesh(z=2.0):
    from chilmesh import CHILmesh

    pts = np.array([[0.0, 0.0, z], [3.0, 0.0, z], [0.0, 3.0, z]])
    return CHILmesh(connectivity=np.array([[0, 1, 2]]), points=pts, grid_name="flat")


def test_station_query_points_cycles_to_requested_count():
    xy = np.array([[0.0, 0.0], [1.0, 1.0]])
    q = measures.station_query_points(xy, 5)
    assert q.shape == (5, 2)
    assert q[2].tolist() == [0.0, 0.0] and q[3].tolist() == [1.0, 1.0]
    assert measures.station_query_points(np.arange(20.0).reshape(10, 2), 4).shape == (4, 2)


def _aux_with_maxele(tmp_path, *, maxele=True):
    _write_stations(tmp_path / noaa.STATIONS_FILE, [{"id": "1", "lat": 1.0, "lng": 1.0}])
    (tmp_path / noaa.POINTS_FILE).write_bytes(b"")
    if maxele:
        (tmp_path / noaa.MAXELE_FILE).write_bytes(b"")
    return tmp_path


_CWL_READER = lambda p: (np.array([[1.0, 1.0], [9.0, 9.0]]), np.array([10.5, 0.0]))  # noqa: E731
_MAXELE_READER = lambda p: np.full(3, 10.0)  # noqa: E731  (nodal max elevation; bathymetry z is 2.0)


def test_maxele_names_and_reason():
    assert noaa.MAXELE_FILE == "maxele.63.nc"
    assert noaa.MAXELE_VAR == "zeta_max"
    assert noaa.MAXELE_REASON == "NOAA inputs not fetched (maxele)"


def test_maxele_path_is_none_when_absent(tmp_path):
    assert noaa.maxele_path(None) is None
    assert noaa.maxele_path(tmp_path) is None
    (tmp_path / noaa.MAXELE_FILE).write_bytes(b"")
    assert noaa.maxele_path(tmp_path) == tmp_path / noaa.MAXELE_FILE


def test_reduce_time_max_takes_max_over_time_axis():
    zeta = np.array([[0.1, 0.5], [0.4, 0.2], [0.3, 0.9]])  # (time, station)
    assert noaa.reduce_time_max(zeta).tolist() == [0.4, 0.9]
    assert noaa.reduce_time_max(np.array([0.1, 0.5])).tolist() == [0.1, 0.5]  # already one value per station


def test_read_maxele_uses_injected_reader(tmp_path):
    f = tmp_path / noaa.MAXELE_FILE
    f.write_bytes(b"")
    assert noaa.read_maxele(f, reader=lambda p: [1, 2, 3]).tolist() == [1.0, 2.0, 3.0]


def test_read_maxele_default_reader_names_missing_variable(tmp_path):
    xr = pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    f = tmp_path / noaa.MAXELE_FILE
    xr.Dataset({"other": ("node", np.zeros(3))}).to_netcdf(f)
    with pytest.raises(noaa.MaxeleError, match="zeta_max"):
        noaa.read_maxele(f)


def test_read_maxele_default_reader_reads_zeta_max(tmp_path):
    xr = pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    f = tmp_path / noaa.MAXELE_FILE
    xr.Dataset({"zeta_max": ("node", np.array([1.0, 2.0, 3.0]))}).to_netcdf(f)
    assert noaa.read_maxele(f).tolist() == [1.0, 2.0, 3.0]


def test_read_points_cwl_default_reader_uses_max_over_time(tmp_path):
    xr = pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    f = tmp_path / noaa.POINTS_FILE
    xr.Dataset({"x": ("station", np.array([1.0, 2.0])), "y": ("station", np.array([3.0, 4.0])),
                "zeta": (("time", "station"), np.array([[0.1, 0.5], [0.4, 0.2], [0.3, 0.1]]))}).to_netcdf(f)
    xy, vals = noaa.read_points_cwl(f)
    assert xy.tolist() == [[1.0, 3.0], [2.0, 4.0]]
    assert vals.tolist() == [0.4, 0.5]


def test_run_station_check_samples_maxele_not_bathymetry(tmp_path):
    aux = _aux_with_maxele(tmp_path)
    out = measures.run_station_check(_flat_mesh(2.0), aux, reader=_CWL_READER, maxele_reader=_MAXELE_READER)
    assert out["status"] == "ok"
    assert out["n_sampled"] == 1 and out["n_outside"] == 1
    assert out["max_abs_diff"] == pytest.approx(0.5)  # |10.0 - 10.5|; bathymetry z=2.0 would give 8.5


def test_run_station_check_without_inputs_is_not_run(tmp_path):
    out = measures.run_station_check(_flat_mesh(), tmp_path)
    assert out == {"status": "not-run", "reason": "NOAA inputs not fetched"}


def test_run_station_check_without_maxele_is_not_run(tmp_path):
    aux = _aux_with_maxele(tmp_path, maxele=False)
    out = measures.run_station_check(_flat_mesh(), aux, reader=_CWL_READER)
    assert out == {"status": "not-run", "reason": "NOAA inputs not fetched (maxele)"}


def test_run_station_check_missing_variable_is_failed_with_reason(tmp_path):
    aux = _aux_with_maxele(tmp_path)

    def bad(path):
        raise noaa.MaxeleError("variable 'zeta_max' not found in maxele.63.nc")

    out = measures.run_station_check(_flat_mesh(), aux, reader=_CWL_READER, maxele_reader=bad)
    assert out["status"] == "failed" and "zeta_max" in out["reason"]


def test_run_station_check_node_count_mismatch_is_failed(tmp_path):
    aux = _aux_with_maxele(tmp_path)
    out = measures.run_station_check(_flat_mesh(), aux, reader=_CWL_READER, maxele_reader=lambda p: np.zeros(5))
    assert out["status"] == "failed" and "5" in out["reason"] and "3" in out["reason"]


def test_render_values_none_without_maxele_never_bathymetry(tmp_path):
    assert measures.render_values(_flat_mesh(2.0), tmp_path) is None
    assert measures.render_values(_flat_mesh(2.0), None) is None
    (tmp_path / noaa.MAXELE_FILE).write_bytes(b"")
    vals = measures.render_values(_flat_mesh(2.0), tmp_path, maxele_reader=_MAXELE_READER)
    assert vals.tolist() == [10.0, 10.0, 10.0]


def test_run_render_check_not_run_without_maxele(tmp_path):
    out = measures.run_render_check(_flat_mesh(2.0), tmp_path)
    assert out["status"] == "not-run" and out["reason"] == "NOAA inputs not fetched (maxele)"
    assert "render_s" not in out or out["render_s"] is None


def test_run_render_check_times_render_with_maxele(tmp_path):
    pytest.importorskip("matplotlib")
    (tmp_path / noaa.MAXELE_FILE).write_bytes(b"")
    out = measures.run_render_check(_flat_mesh(2.0), tmp_path, maxele_reader=_MAXELE_READER)
    assert out["status"] == "ok" and out["render_s"] >= 0


def test_station_sample_max_diff_docstring_is_readable():
    doc = measures.station_sample_max_diff.__doc__
    assert "barycentric/inverse-distance-free" not in doc
    assert "inverse-distance" in doc


def test_locate_uses_station_coords_when_present(tmp_path):
    _write_stations(tmp_path / noaa.STATIONS_FILE, [{"lat": 1.0, "lng": 1.0}, {"lat": 0.5, "lng": 0.5}])
    q, src = measures.locate_queries(_flat_mesh(), 6, aux_dir=tmp_path)
    assert src == "coops-stations" and q.shape == (6, 2)
    assert q[0].tolist() == [1.0, 1.0]


def test_locate_falls_back_to_centroids_without_stations(tmp_path):
    q, src = measures.locate_queries(_flat_mesh(), 6, aux_dir=tmp_path)
    assert src == "element-centroids" and q.shape == (6, 2)


def _meta(name="M", n_elems=10):
    return run._mesh_meta(name, {}, "c", "u", n_nodes=n_elems, n_elems=n_elems)


def test_station_record_not_run_when_inputs_absent(tmp_path):
    def boom(*a, **k):
        raise AssertionError("worker must not start without NOAA inputs")

    rec = run.station_record("m.14", _meta(), "unknown", tmp_path, 10, worker=boom)
    assert rec["status"] == "not-run" and rec["reason"] == "NOAA inputs not fetched"
    assert run.station_record("m.14", _meta(), "unknown", None, 10, worker=boom)["status"] == "not-run"
    records.validate_record(rec)


def test_station_record_not_run_when_maxele_absent(tmp_path):
    _aux_with_maxele(tmp_path, maxele=False)

    def boom(*a, **k):
        raise AssertionError("worker must not start without maxele")

    rec = run.station_record("m.14", _meta(), "unknown", tmp_path, 10, worker=boom)
    assert rec["status"] == "not-run" and rec["reason"] == "NOAA inputs not fetched (maxele)"
    records.validate_record(rec)


def test_station_record_failed_when_worker_reports_failed(tmp_path):
    _aux_with_maxele(tmp_path)

    def fake(*a, **k):
        return {"status": "ok", "reason": "", "command": ["x"], "stdout_tail": "", "stderr_tail": "",
                "json": {"status": "failed", "reason": "variable 'zeta_max' not found in maxele.63.nc"}}

    rec = run.station_record("m.14", _meta(), "cold", tmp_path, 10, worker=fake)
    assert rec["status"] == "failed" and "zeta_max" in rec["reason"]
    records.validate_record(rec)


def test_station_record_ok_when_worker_reports_diff(tmp_path):
    _aux_with_maxele(tmp_path)
    seen = {}

    def fake(op, mesh_path, py, env_extra, timeout_s, extra_args=()):
        seen["op"], seen["args"] = op, list(extra_args)
        return {"status": "ok", "reason": "", "command": ["x"], "stdout_tail": "", "stderr_tail": "",
                "json": {"status": "ok", "max_abs_diff": 0.5, "n_sampled": 3, "n_outside": 1}}

    rec = run.station_record("m.14", _meta(), "cold", tmp_path, 10, worker=fake)
    assert seen["op"] == "stations" and "--aux-dir" in seen["args"]
    assert rec["status"] == "ok"
    assert rec["metric_values"]["station_max_abs_diff"]["value"] == 0.5


def test_station_record_failed_when_every_station_outside(tmp_path):
    _aux_with_maxele(tmp_path)

    def fake(*a, **k):
        return {"status": "ok", "reason": "", "command": ["x"], "stdout_tail": "", "stderr_tail": "",
                "json": {"status": "ok", "max_abs_diff": None, "n_sampled": 0, "n_outside": 4}}

    rec = run.station_record("m.14", _meta(), "cold", tmp_path, 10, worker=fake)
    assert rec["status"] == "failed" and "outside" in rec["reason"]


# ---------------------------------------------------------------- render op in the chilmesh run

def _fake_worker(render_json, seen):
    def fake(op, mesh_path, py, env_extra, timeout_s, extra_args=()):
        seen.append((op, list(extra_args)))
        js = {"parse": {"parse_s": 0.1, "n_nodes": 3, "n_elems": 1},
              "adjacency": {"adjacency_s": 0.2},
              "locate": {"tree_build_s": 0.01, "query_per_s": 1e-5},
              "render": render_json}[op]
        return {"status": "ok", "reason": "", "command": ["x"], "stdout_tail": "", "stderr_tail": "", "json": js}
    return fake


def test_render_not_run_without_maxele_is_na_cell_not_failure(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(run, "_worker", _fake_worker(
        {"status": "not-run", "reason": "NOAA inputs not fetched (maxele)", "render_s": None}, seen))
    rec = run.run_chilmesh_once("m.14", _meta(), "cold", 0, 10, aux_dir=tmp_path)
    assert rec["status"] == "ok"
    cell = rec["metric_values"]["render_s"]
    assert cell["na"] and "maxele" in cell["reason"]
    assert ("render", ["--aux-dir", str(tmp_path)]) in seen  # render worker gets the aux dir


def test_render_ok_records_value(monkeypatch):
    monkeypatch.setattr(run, "_worker", _fake_worker({"status": "ok", "render_s": 0.7}, []))
    rec = run.run_chilmesh_once("m.14", _meta(), "cold", 0, 10)
    assert rec["status"] == "ok" and rec["metric_values"]["render_s"]["value"] == 0.7


# ---------------------------------------------------------------- scaling collection

def _chil_rec(name, n_elems, q=None, status="ok", reason=""):
    mv = {} if q is None else {"query_per_s": records.metric(q)}
    return run._record("chilmesh", "0", "PolyForm-NC", "subprocess", _meta(name, n_elems), mv, status, reason,
                       "cold", [], None, 0)


def test_collect_scaling_fits_four_sizes(tmp_path):
    sizes = [1_000, 10_000, 100_000, 1_000_000]
    recs = [_chil_rec(f"m{n}", n, 2e-5 + 1e-5 * math.log(n)) for n in sizes]
    out = run.collect_scaling(recs, tmp_path)
    assert out["sizes"] == sizes
    assert out["summary"]["best"] == "log-n"
    assert out["not_run"] == []
    saved = json.loads((tmp_path / "scaling" / "w2_scaling.json").read_text(encoding="utf-8"))
    assert saved["summary"]["best"] == "log-n"


def test_collect_scaling_records_unloadable_size_as_not_run(tmp_path):
    recs = [_chil_rec("a", 1_000, 1e-5), _chil_rec("b", 10_000, 2e-5),
            _chil_rec("c", 0, None, status="failed", reason="parse failed: boom")]
    out = run.collect_scaling(recs, tmp_path)
    assert out["sizes"] == [1_000, 10_000]
    assert out["summary"]["best"] is None and "only 2" in out["summary"]["note"]
    assert out["not_run"] == [{"mesh": "c", "reason": "scaling: parse failed: boom"}]
    nr = [r for r in records.read_records(tmp_path) if r["status"] == "not-run"]
    assert len(nr) == 1 and nr[0]["mesh"]["name"] == "c" and nr[0]["reason"].startswith("scaling:")


def test_collect_scaling_reuses_existing_not_run_record(tmp_path):
    existing = [_chil_rec("x", 0, None, status="not-run", reason="Valence fetch failed")]
    out = run.collect_scaling(existing, tmp_path)
    assert out["not_run"][0]["mesh"] == "x"
    assert list(tmp_path.glob("*.json")) == []  # nothing duplicated


def test_cli_w2_collects_scaling_across_meshes(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("CHILMESH_RUN_BENCH", "1")
    fake_sizes = {"A": 1_000, "B": 10_000, "C": 100_000}
    calls = []

    def fake_run_mesh(name, out_dir, **kw):
        calls.append((name, kw.get("aux_dir")))
        return [_chil_rec(name, fake_sizes[name], 1e-5 * math.log(fake_sizes[name]))]

    monkeypatch.setattr(run, "run_mesh", fake_run_mesh)
    rc = cli.main(["w2", "--mesh", "A", "--mesh", "B", "--mesh", "C", "--out", str(tmp_path),
                   "--aux-dir", str(tmp_path / "aux")])
    assert rc == 0
    assert [c[0] for c in calls] == ["A", "B", "C"]
    assert all(c[1] == str(tmp_path / "aux") for c in calls)
    saved = json.loads((tmp_path / "scaling" / "w2_scaling.json").read_text(encoding="utf-8"))
    assert saved["sizes"] == [1_000, 10_000, 100_000]
    assert "scaling" in capsys.readouterr().out.lower()
