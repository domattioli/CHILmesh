"""Fake tool modules executed by standalone adapter scripts, never imported."""
import json
import os
import subprocess
import sys
from pathlib import Path
import pytest

ADAPTERS = Path(__file__).resolve().parents[2] / 'benchmarks/competitors/adapters'

# Minimal writer for the pre-READY warm-up file (protocol.warm_netcdf_reader); the fake
# readers return a no-op object for it, so the measured load still owns its deferred imports.
NETCDF4_STUB = '''
class _Var:
    def __setitem__(self, key, value): pass
class Dataset:
    def __init__(self, path, mode='r', **kwargs): open(path, 'w').close()
    def __enter__(self): return self
    def __exit__(self, *exc): return False
    def createDimension(self, name, size): pass
    def createVariable(self, *args, **kwargs): return _Var()
'''


@pytest.mark.parametrize('tool', ['adcircpy', 'thalassa', 'xarray'])
def test_metadata_and_eager_materialization(tmp_path, tool):
    (tmp_path / (tool + '-1.2.3.dist-info')).mkdir()
    (tmp_path / (tool + '-1.2.3.dist-info') / 'METADATA').write_text(f'Name: {tool}\nVersion: 1.2.3\n')
    fake = '''
class Dataset:
    sizes = {'node': 4, 'face': 2, 'nele': 2}
    def load(self):
        self.loaded = True
        return self
    def __getitem__(self, key):
        assert self.loaded, 'lazy access before load'
        raise KeyError(key)
class Warm:
    def load(self): return self
    def close(self): pass
def open_dataset(path, *args, **kwargs): return Warm() if str(path).endswith('warm.nc') else Dataset()
'''
    if tool != 'adcircpy': (tmp_path / 'netCDF4.py').write_text(NETCDF4_STUB)
    if tool == 'thalassa':
        (tmp_path / 'xarray.py').write_text('')
        (tmp_path / tool).mkdir(); (tmp_path / tool / '__init__.py').write_text(''); (tmp_path / tool / 'api.py').write_text(fake)
    elif tool == 'adcircpy':
        (tmp_path / (tool + '.py')).write_text('class AdcircMesh:\n    @staticmethod\n    def open(*args, **kwargs): return object()\n')
    else: (tmp_path / (tool + '.py')).write_text(fake)
    env = {**os.environ, 'PYTHONPATH': str(tmp_path) + os.pathsep + os.environ['PYTHONPATH']}
    p = subprocess.run([sys.executable, str(ADAPTERS / (tool + '_read.py')), '--mesh', 'not-a-source', '--fair'], env=env, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    data = json.loads(p.stdout)
    assert data['version'] == '1.2.3'
    assert data['materialization'] is True
    assert data['output']['verification_state'] == 'incomplete'
    assert 'mapping' in data['output']['reason'].lower()


def test_chilmesh_parse_preparation_does_not_read_input():
    from benchmarks.competitors.w2 import measures
    events = []
    class Mesh:
        @classmethod
        def read_from_fort14(cls, path, **kw): events.append('parse'); return cls()
        n_verts, n_elems = 4, 2
    prepared = measures.prepare_fair_operation('parse', 'synthetic', mesh_class=Mesh)
    assert events == []
    prepared['operation']()
    assert events == ['parse']


@pytest.mark.parametrize('operation,expected', [('locate-tree-build', ['parse', 'tree']), ('locate-query', ['parse', 'tree', 'query', 'query'])])
def test_locate_separate_region_work(operation, expected):
    import numpy as np
    from benchmarks.competitors.w2 import measures
    events = []
    class Mesh:
        n_verts, n_elems = 3, 1
        points = np.array([[0., 0., 1.], [1., 0., 2.], [0., 1., 3.]])
        connectivity_list = np.array([[0, 1, 2]])
        @classmethod
        def read_from_fort14(cls, path, **kw): events.append('parse'); return cls()
        def _build_spatial_indices(self): events.append('tree')
        def find_element(self, point): events.append('query'); return 0
    prepared = measures.prepare_fair_operation(operation, 'synthetic', mesh_class=Mesh, n_queries=2)
    assert events == (['parse', 'tree'] if operation == 'locate-query' else ['parse'])
    prepared['operation']()
    assert events == expected


@pytest.mark.parametrize('tool,dim', [('xarray', 'nele'), ('thalassa', 'face')])
def test_adapter_control_pipe_ack_precedes_hashes(tmp_path, fair_fixture, tool, dim):
    import time
    from benchmarks.competitors.reference import build_reference
    from benchmarks.competitors.verification import verify_output
    from benchmarks.competitors.protocol import ProtocolState
    reference = build_reference(fair_fixture['mesh'], fair_fixture['output'])
    distribution = tmp_path / (tool + '-1.2.3.dist-info'); distribution.mkdir()
    (distribution / 'METADATA').write_text(f'Name: {tool}\nVersion: 1.2.3\n')
    fake = '''
class Values:
    def __init__(self, value): self.value = value
    def tolist(self): return self.value
class Var:
    def __init__(self, value): self.values = Values(value)
class Dataset:
    sizes = {'node': 4, 'nele': 2, 'face': 2}
    def load(self):
        import deferred_marker
        self.loaded = True
        return self
    def __getitem__(self, key):
        assert self.loaded
        return Var({'node_id':[10,20,40,80], 'element_id':[7,9],
                    'x':[0.,1.,1.,0.], 'y':[-0.,0.,1.,1.], 'depth':[1.,2.,3.,4.],
                    'lon':[0.,1.,1.,0.], 'lat':[-0.,0.,1.,1.],
                    'arity':[3,4], 'connectivity':[[10,20,40,0],[10,20,40,80]],
                    'face_nodes':[[0,1,2,-1],[0,1,2,3]]}[key])
class Warm:
    def load(self): return self
    def close(self): pass
def open_dataset(path, *args, **kwargs): return Warm() if str(path).endswith('warm.nc') else Dataset()
'''
    (tmp_path / 'deferred_marker.py').write_text('')
    (tmp_path / 'netCDF4.py').write_text(NETCDF4_STUB)
    if tool == 'thalassa':
        (tmp_path / 'xarray.py').write_text('')
        (tmp_path / tool).mkdir(); (tmp_path / tool / '__init__.py').write_text(''); (tmp_path / tool / 'api.py').write_text(fake)
    else: (tmp_path / (tool + '.py')).write_text(fake)
    incoming_read, incoming_write = os.pipe(); outgoing_read, outgoing_write = os.pipe()
    result_path = tmp_path / 'result.json'
    with result_path.open('w') as result, os.fdopen(incoming_write, 'w', buffering=1) as send, os.fdopen(outgoing_read) as receive:
        env = {**os.environ, 'PYTHONPATH': str(tmp_path) + os.pathsep + os.environ['PYTHONPATH'],
               'CHILMESH_CONTROL_IN_FD': str(incoming_read), 'CHILMESH_CONTROL_OUT_FD': str(outgoing_write),
               'CHILMESH_RESULT_FD': str(result.fileno()), 'CHILMESH_REPEAT_ID': 'test'}
        p = subprocess.Popen([sys.executable, str(ADAPTERS / (tool + '_read.py')), '--mesh', 'unread-source', '--fair'],
                             env=env, pass_fds=(incoming_read, outgoing_write, result.fileno()), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        os.close(incoming_read); os.close(outgoing_write)
        try:
            state = ProtocolState('test')
            for event in ('READY', 'START', 'END'):
                msg = json.loads(receive.readline()); state.accept(msg)
                assert msg['event'] == event
                if event in ('READY', 'END'):
                    assert result_path.read_text() == ''
                    reply = {**msg, 'event': 'GO' if event == 'READY' else 'ACK', 'monotonic_s': time.perf_counter()}
                    state.accept(reply); send.write(json.dumps(reply) + '\n')
            stdout, stderr = p.communicate(timeout=10)
            assert p.returncode == 0, stderr
            assert state.complete and stdout == ''
        finally:
            if p.poll() is None: p.kill(); p.wait()
    data = json.loads(result_path.read_text())
    assert data['version'] == '1.2.3' and data['new_modules'] == ['deferred_marker']
    assert verify_output(data['output'], reference)['verification_state'] == 'verified'
    assert data['output']['node_sha256'] == reference['canonical_node_sha256']


@pytest.mark.parametrize('operation,output', [('locate-query', {'element_indices': [], 'n_queries': 0}), ('locate-tree-build', {'tree_built': False}), ('stations', {'n_sampled': 0, 'max_abs_diff': None})])
def test_operation_completion_rejects_empty_or_missing_work(operation, output):
    from benchmarks.competitors.w2.measures import verify_fair_operation
    import types
    prepared = {'mesh': types.SimpleNamespace(n_elems=2), 'query_sha256': 'a' * 64}
    assert verify_fair_operation(operation, prepared, output)['verification_state'] == 'failed'


def test_render_and_write_require_persistent_completed_artifacts(tmp_path):
    from benchmarks.competitors.w2.measures import verify_fair_operation
    path = tmp_path / 'output'
    assert verify_fair_operation('render', {}, {'path': str(path)})['verification_state'] == 'incomplete'
    path.write_bytes(b'not a PNG')
    assert verify_fair_operation('render', {}, {'path': str(path)})['verification_state'] == 'failed'
    assert verify_fair_operation('write', {}, {'path': str(path)})['sha256']


def test_render_and_stations_missing_prerequisite_before_start(tmp_path):
    import types
    from benchmarks.competitors.w2.measures import prepare_fair_operation
    class Mesh:
        @classmethod
        def read_from_fort14(cls, *args, **kw): return types.SimpleNamespace()
    for operation in ('render', 'stations'):
        with pytest.raises(ValueError, match='maxele'):
            prepare_fair_operation(operation, 'synthetic', mesh_class=Mesh, aux_dir=tmp_path)


def test_chilmesh_protocol_materialized_output_after_ack(fair_fixture, monkeypatch):
    import io
    import types
    import numpy as np
    from benchmarks.competitors import protocol
    from benchmarks.competitors.w2 import measures
    incoming = io.StringIO(''.join(json.dumps(dict(protocol_version=1, repeat_id='r', event=event,
        monotonic_s=stamp, clock='perf_counter', clock_metadata=protocol.clock_metadata())) + '\n'
        for event, stamp in [('GO', 1.), ('ACK', 4.)]))
    controls, results = io.StringIO(), io.StringIO()
    child = protocol.ChildProtocol('r', incoming, controls, results, clock=iter([0., 2., 3.]).__next__)
    monkeypatch.setattr(protocol, 'protocol_from_environment', lambda: child)
    monkeypatch.setattr(measures.importlib.metadata, 'version', lambda name: '1.2.3' if name == 'chilmesh' else pytest.fail(name))
    data = fair_fixture['adapter']
    mesh = types.SimpleNamespace(n_verts=4, n_elems=2, node_ids=data['node_ids'], element_ids=data['element_ids'],
        points=np.array(data['coordinates']), connectivity_list=[[0, 1, 2], [0, 1, 2, 3]])
    def load(*args, **kw):
        assert child.state.events[-1]['event'] == 'START'
        return mesh
    monkeypatch.setattr(measures, '_load', load)
    extract = measures.chilmesh_output
    def checked(output):
        assert child.state.complete
        return extract(output)
    monkeypatch.setattr(measures, 'chilmesh_output', checked)
    assert measures.fair_main(types.SimpleNamespace(op='parse', mesh='never-read-source', n_queries=1, aux_dir=None, stations=None, output=None)) == 0
    actual = json.loads(results.getvalue())
    assert actual['version'] == '1.2.3' and actual['materialization'] is True
    assert actual['actual_backend']['operation_implementation']
    from benchmarks.competitors.reference import build_reference
    from benchmarks.competitors.verification import verify_output
    reference = build_reference(fair_fixture['mesh'], fair_fixture['output'])
    assert verify_output(actual['output'], reference)['verification_state'] == 'verified'
