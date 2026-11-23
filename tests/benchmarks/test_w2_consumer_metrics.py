"""Consumer metrics on synthetic meshes, without external benchmark inputs."""
import pytest

from chilmesh import CHILmesh
from benchmarks.competitors.w2 import measures, run, tools
from tests.benchmarks._grids import tri_grid


@pytest.fixture
def grid_file(tmp_path):
    points, elements = tri_grid(2, 2)
    mesh = CHILmesh(connectivity=elements, points=points, compute_layers=False, compute_adjacencies=False)
    path = tmp_path / 'grid.14'
    mesh.write_to_fort14(str(path))
    return path


def test_write_roundtrip(grid_file):
    result = measures.write_roundtrip(str(grid_file), measures._load(str(grid_file)))
    assert result['write_fort14_s'] >= 0
    assert result['roundtrip_ok'] is True and result['reason'] == ''


def test_parse_timer_starts_after_class_import(monkeypatch, capsys):
    import builtins
    from types import SimpleNamespace

    events = []
    original_import = builtins.__import__

    def record_import(name, *args, **kwargs):
        result = original_import(name, *args, **kwargs)
        if name == 'chilmesh':
            events.append('import')
        return result

    def load(path, *, mesh_class):
        events.append('load')
        return SimpleNamespace(n_verts=3, n_elems=1)

    def clock():
        events.append('clock')
        return 1.0

    monkeypatch.setattr(builtins, '__import__', record_import)
    monkeypatch.setattr(measures, '_load', load)
    monkeypatch.setattr(measures.time, 'perf_counter', clock)
    assert measures.main(['parse', '--mesh', 'unused']) == 0
    assert events == ['import', 'clock', 'load', 'clock']
    assert '"parse_s": 0.0' in capsys.readouterr().out


def test_roundtrip_detects_changed_nodes(grid_file):
    mesh = measures._load(str(grid_file))
    mesh.points[0, 0] += 1
    result = measures.write_roundtrip(str(grid_file), mesh)
    assert result['roundtrip_ok'] is False and 'nodes differ' in result['reason']


def test_full_init_records_backend(grid_file):
    result = measures.full_init(measures._load(str(grid_file)))
    assert result['full_init_s'] >= 0
    assert 'selected' in result['backend_info']


def test_install_size_and_metric_registration(tmp_path, monkeypatch):
    file = tmp_path / 'installed.py'
    file.write_bytes(b'x' * 1024)
    class Distribution:
        files = ['installed.py']
        def locate_file(self, name):
            return tmp_path / name
    monkeypatch.setattr(measures.importlib.metadata, 'distribution', lambda name: Distribution())
    assert measures.install_mb() == 1 / 1024
    assert {'import_s', 'install_mb', 'write_fort14_s', 'roundtrip_ok', 'full_init_s',
            'peak_rss_per_melem', 'interp_max_abs_diff'} <= set(tools.METRICS)


def test_subprocess_metrics_and_rss(grid_file):
    mesh = run._mesh_meta('tiny', {}, 'unknown', 'unknown')
    record = run.run_chilmesh_once(grid_file, mesh, 'unknown', 0, 30)
    cells = record['metric_values']
    assert cells['import_s']['value'] > 0
    assert cells['write_fort14_s']['value'] >= 0
    assert cells['roundtrip_ok']['value'] is True
    assert cells['full_init_s']['backend_info']['selected']
    assert cells['peak_rss_per_melem']['value']['parse'] > 0
    assert cells['peak_rss_per_melem']['value']['full-init'] > 0


def test_interpolation_missing_inputs_is_na(tmp_path):
    mesh = run._mesh_meta('tiny', {}, 'unknown', 'unknown')
    cell = run.station_record('unused', mesh, 'unknown', tmp_path, 10)['metric_values']['interp_max_abs_diff']
    assert cell == {'value': None, 'na': True, 'reason': 'NOAA inputs not fetched'}
