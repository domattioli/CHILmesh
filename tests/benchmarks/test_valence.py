"""Valence wrapper: streaming download, not-served handling, rights filter, no fallback."""
from __future__ import annotations

import hashlib
import io
import json
import types
from pathlib import Path

import pytest

from benchmarks.competitors import valence as V


class _FakeResp(io.BytesIO):
    def __init__(self, data):
        super().__init__(data)
        self.max_read = 0

    def read(self, n=-1):
        self.max_read = max(self.max_read, n if n and n > 0 else 10**12)
        return super().read(n)

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class _NotFound(Exception):
    pass


class _Unavailable(Exception):
    pass


class _Auth(Exception):
    pass


def _fake_registry(digest="", *, hash_exc=None):
    def fetch_hash(name):
        if hash_exc:
            raise hash_exc
        return digest

    return types.SimpleNamespace(
        RegistryNotFoundError=_NotFound,
        RegistryUnavailableError=_Unavailable,
        RegistryAuthError=_Auth,
        RegistryError=Exception,
        _config=lambda: ("https://valence.example", "tok"),
        _check_name=lambda n: n,
        fetch_hash=fetch_hash,
        list_meshes=lambda: ["Shinnecock"],
    )


def test_streaming_download_is_chunked_and_verified(tmp_path):
    payload = b"x" * (3 * 1024 + 7)
    reg = _fake_registry(hashlib.sha256(payload).hexdigest())
    resp = _FakeResp(payload)
    calls = []

    def opener(req, timeout=None):
        calls.append(req.full_url)
        return resp

    res = V.fetch_mesh_file("Shinnecock", tmp_path, registry=reg, opener=opener, chunk_bytes=1024)
    assert res["status"] == "ok"
    assert res["sha256"] == hashlib.sha256(payload).hexdigest()
    assert (tmp_path / "Shinnecock.14").read_bytes() == payload
    assert resp.max_read <= 1024  # never reads the whole payload at once
    assert calls == ["https://valence.example/meshes/Shinnecock"]


def test_hash_mismatch_is_failed_and_file_removed(tmp_path):
    reg = _fake_registry("0" * 64)
    res = V.fetch_mesh_file("Shinnecock", tmp_path, registry=reg,
                            opener=lambda req, timeout=None: _FakeResp(b"abc"))
    assert res["status"] == "failed"
    assert "hash" in res["reason"]
    assert not (tmp_path / "Shinnecock.14").exists()


@pytest.mark.parametrize("exc", [_NotFound("nope"), _Unavailable("down")])
def test_not_served_or_unavailable_is_not_run_without_fallback(tmp_path, exc):
    reg = _fake_registry(hash_exc=exc)
    opened = []

    def opener(req, timeout=None):  # must never be reached
        opened.append(req)
        raise AssertionError("fallback download attempted")

    res = V.fetch_mesh_file("STOFS", tmp_path, registry=reg, opener=opener)
    assert res["status"] == "not-run"
    assert res["reason"] == "mesh not yet served by Valence"
    assert opened == []
    assert list(tmp_path.iterdir()) == []


def test_missing_token_is_not_run_with_distinct_reason(tmp_path):
    reg = _fake_registry(hash_exc=_Auth("no token"))
    res = V.fetch_mesh_file("Shinnecock", tmp_path, registry=reg, opener=None)
    assert res["status"] == "not-run"
    assert "VALENCE_TOKEN" in res["reason"]


def test_module_has_no_noaa_or_local_clone_fallback():
    src = open(V.__file__, encoding="utf-8").read().lower()
    assert "noaa" not in src
    assert "valence_data_dir" not in src


def test_rights_manual_file_excludes_unpermitted_meshes(tmp_path):
    rights = tmp_path / "rights.json"
    rights.write_text(json.dumps({
        "Shinnecock": {"publish_results": True, "note": "ok"},
        "WNAT_Hagen": {"publish_results": False, "note": "no"},
    }))
    allowed, excluded = V.apply_rights(["Shinnecock", "WNAT_Hagen", "ENPAC2003"], V.load_rights(rights))
    assert allowed == ["Shinnecock"]
    assert {e["mesh"] for e in excluded} == {"WNAT_Hagen", "ENPAC2003"}
    assert all(e["reason"] for e in excluded)


def test_rights_file_missing_excludes_everything(tmp_path):
    allowed, excluded = V.apply_rights(["Shinnecock"], V.load_rights(tmp_path / "absent.json"))
    assert allowed == []
    assert excluded[0]["mesh"] == "Shinnecock"


def test_fair_reference_bound_to_exact_asset(fair_fixture):
    p = fair_fixture['mesh']; raw = p.read_bytes(); digest = hashlib.sha256(raw).hexdigest()
    fetched = {'status': 'ok', 'path': str(p), 'sha256': digest, 'bytes': len(raw)}
    out = V.reference_from_fetch('synthetic', 'v1', 'asset-1', fetched, fair_fixture['output'])
    assert out['exact_version'] == 'v1' and out['valence_asset'] == 'asset-1'
    assert out['served_sha256'] == out['raw_file_sha256'] == digest
    assert V.reference_from_fetch('x', 'v1', 'asset', {'status': 'not-run', 'reason': 'not served'}, fair_fixture['output'])['status'] == 'not-run'


def test_conversion_missing_proof_and_costs(fair_fixture):
    from benchmarks.competitors.campaign import prepare_conversion, validate_conversion
    assert prepare_conversion(None, 'xarray', fair_fixture['output'])['status'] == 'not-run'
    assert validate_conversion({'conversion_s': 2., 'read_s': 3.})['verification_state'] == 'incomplete'


@pytest.mark.parametrize('reader,dimension', [('xarray', 'nele'), ('thalassa', 'nele')])
def test_independent_conversion_roundtrip_and_lineage(fair_fixture, tmp_path, reader, dimension):
    import numpy as np
    from benchmarks.competitors.reference import build_reference
    from benchmarks.competitors.scripts.convert_fair_netcdf import convert
    from benchmarks.competitors.campaign import prepare_conversion, validate_conversion
    reference = build_reference(fair_fixture['mesh'], fair_fixture['output'])
    stores = {}
    class Variable:
        def __init__(self, dtype, shape): self.data = np.zeros(shape, dtype=dtype); self.dtype = self.data.dtype
        def __setitem__(self, key, value): self.data[key] = value
        def __getitem__(self, key): return self.data[key]
    class Dataset:
        def __init__(self, path, mode, **kw):
            self.path = path
            if mode == 'w':
                self.dimensions, self.variables = {}, {}; stores[path] = self
                Path(path).write_bytes(b'synthetic NetCDF container')
            else:
                self.dimensions = stores[path].dimensions; self.variables = stores[path].variables
        def createDimension(self, name, size): self.dimensions[name] = range(size)
        def createVariable(self, name, dtype, dims, fill_value):
            assert fill_value is False
            self.variables[name] = Variable(dtype, [len(self.dimensions[d]) for d in dims])
        def __getitem__(self, key): return self.variables[key]
        def set_auto_mask(self, value): assert value is False
        def __enter__(self): return self
        def __exit__(self, *args): pass
    fake = types.SimpleNamespace(Dataset=Dataset, __version__='1.2.3', __netcdf4libversion__='4.9.0', __hdf5libversion__='1.14.0')
    def runner(argv, **kw):
        assert argv[2:4] == ['--reference', argv[3]]
        lineage = convert(json.loads(Path(argv[3]).read_text()), reader, argv[-1], netcdf=fake)
        return types.SimpleNamespace(returncode=0, stdout=json.dumps(lineage), stderr='')
    conversion = prepare_conversion(reference, reader, tmp_path / 'convert', runner=runner,
                                    converter_evidence={'lock_sha256': 'a' * 64, 'compatible': True}, clock=iter([10., 12.]).__next__)
    assert conversion['status'] == 'ok', conversion
    assert conversion['contract']['dimensions'] == {'node': 4, 'nele': 2, 'max_face_nodes' if reader == 'xarray' else 'nvertex': 4}
    assert conversion['source_content_verification']['node_sha256'] == reference['canonical_node_sha256']
    assert conversion['conversion_s'] == 2. and conversion['read_s'] is None and conversion['end_to_end_s'] is None
    with_read = validate_conversion({**conversion, 'read_s': 3.}, reference)
    assert with_read['end_to_end_s'] == 5.
    assert validate_conversion({**conversion, 'source_sha256': '0' * 64}, reference)['verification_state'] == 'incomplete'
    corrupt = {**conversion, 'source_content_verification': {**conversion['source_content_verification'], 'element_sha256': '0' * 64}}
    assert validate_conversion(corrupt, reference)['verification_state'] == 'incomplete'
    assert validate_conversion({**conversion, 'native_versions': {}}, reference)['verification_state'] == 'incomplete'
    path = tmp_path / 'bad.nc'
    Path(reference['nodes_path']).write_text('tampered')
    with pytest.raises(ValueError, match='hash'): convert(reference, reader, path, netcdf=fake)


@pytest.mark.parametrize('change', [None, {}, {'publish_results': True}, {'publish_results': False},
    {'exact_version':'wrong'}, {'mesh':'wrong'}, {'citation':'https://example.org/issues/1'},
    {'ruling_date':'yesterday'}, {'scope':''}, {'evidence_sha256':''}, {'ambiguous':True}])
def test_fair_rights_default_deny_and_exact_authority(change):
    from _fair_saved import ruling
    entry=None if change is None else dict(ruling(),**change)
    if change == {}: entry={}
    if change == {'publish_results':True}: entry=change
    result=V.fair_rights_decision('STOFS','v1',entry)
    assert not result['allowed'] and result['reasons']


def test_fair_rights_admits_only_exact_dated_ruling():
    from _fair_saved import ruling
    decision=V.fair_rights_decision('STOFS','v1',ruling())
    assert decision['allowed'] and decision['citation'].endswith('/123') and decision['ruling_date']=='2026-11-10'
    for name in ['Shinnecock','WNAT_Hagen','ENPAC2003','STOFS']:
        assert not V.fair_rights_decision(name,'v1',None)['allowed']


def _fake_netcdf(corrupt=None):
    import numpy as np
    stores = {}
    class Variable:
        def __init__(self, name, dtype, shape): self.name = name; self.data = np.zeros(shape, dtype=dtype); self.dtype = self.data.dtype
        def __setitem__(self, key, value):
            if corrupt and self.name == corrupt[0] and key[0] == corrupt[1]: value = np.asarray(value) + 1
            self.data[key] = value
        def __getitem__(self, key): return self.data[key]
    class Dataset:
        def __init__(self, path, mode, **kw):
            if mode == 'w':
                self.dimensions, self.variables = {}, {}; stores[path] = self
                Path(path).write_bytes(b'synthetic NetCDF container')
            else: self.dimensions, self.variables = stores[path].dimensions, stores[path].variables
        def createDimension(self, name, size): self.dimensions[name] = range(size)
        def createVariable(self, name, dtype, dims, fill_value):
            assert fill_value is False
            self.variables[name] = Variable(name, dtype, [len(self.dimensions[d]) for d in dims])
        def __getitem__(self, key): return self.variables[key]
        def set_auto_mask(self, value): pass
        def __enter__(self): return self
        def __exit__(self, *args): pass
    return types.SimpleNamespace(Dataset=Dataset, __version__='1.2.3', __netcdf4libversion__='4.9.0', __hdf5libversion__='1.14.0'), stores


def test_thalassa_layout_converts_and_validates(fair_fixture, tmp_path):
    from benchmarks.competitors.reference import build_reference
    from benchmarks.competitors.scripts.convert_fair_netcdf import convert
    from benchmarks.competitors.campaign import validate_conversion
    reference = build_reference(fair_fixture['mesh'], fair_fixture['output'])
    fake, stores = _fake_netcdf()
    lineage = convert(reference, 'thalassa', tmp_path / 'ok.nc', netcdf=fake)
    ds = stores[str(tmp_path / 'ok.nc')]
    assert set(ds.dimensions) == {'node', 'nele', 'nvertex'} and 'connectivity' not in ds.variables
    assert ds['element'].data.min() >= 0 and ds['element'].data.max() <= 4
    assert lineage['contract']['indexing'] == 'one-based-row-index' and lineage['contract']['fill_value'] == 0
    assert lineage['source_content_verification']['node_sha256'] == reference['canonical_node_sha256']
    full = {**lineage, 'conversion_s': 1., 'read_s': 1., 'environment_lock_sha256': 'a' * 64, 'argv': ['c'], 'options': {'r': 1}}
    assert validate_conversion(full, reference)['verification_state'] == 'verified'
    wrong = {**full, 'contract': {**lineage['contract'], 'indexing': 'original-file-ids'}}
    assert validate_conversion(wrong, reference)['verification_state'] == 'incomplete'


def test_thalassa_layout_wrong_element_index_fails(fair_fixture, tmp_path):
    from benchmarks.competitors.reference import build_reference
    from benchmarks.competitors.scripts.convert_fair_netcdf import convert
    reference = build_reference(fair_fixture['mesh'], fair_fixture['output'])
    fake, _ = _fake_netcdf(corrupt=('element', 0))
    with pytest.raises(ValueError):
        convert(reference, 'thalassa', tmp_path / 'bad.nc', netcdf=fake)


def test_thalassa_output_extractor_with_padded_quad():
    from benchmarks.competitors.protocol import thalassa_output
    class V:
        def __init__(self, v): self.values = types.SimpleNamespace(tolist=lambda: v)
    class DS:
        sizes = {'node': 4, 'face': 2, 'max_no_vertices': 4}
        def __init__(self, **over): self.data = {'node_id': [10, 20, 40, 80], 'element_id': [7, 9], 'lon': [0., 1., 1., 0.],
                'lat': [0., 0., 1., 1.], 'depth': [1., 2., 3., 4.], 'arity': [3, 4], 'face_nodes': [[0, 1, 2, -1], [0, 1, 2, 3]], **over}
        def __getitem__(self, k): return V(self.data[k])
    out = thalassa_output(DS())
    assert out['verification_state'] == 'materialized' and out['connectivity'] == [[0, 1, 2], [0, 1, 2, 3]]
    assert out['connectivity_indexing'] == 'row-offset' and out['index_to_node_id'] == [10, 20, 40, 80]
    assert out['n_elems'] == 2 and 'face_nodes = element - 1' in out['mapping_provenance']['kind']
    assert thalassa_output(DS(arity=[3., 4.]))['verification_state'] == 'incomplete'
    assert thalassa_output(DS(face_nodes=[[0, 1, 2, 0], [0, 1, 2, 3]]))['verification_state'] == 'incomplete'


def test_dest_name_not_double_suffixed(tmp_path):
    payload = b"mesh"
    reg = _fake_registry(hashlib.sha256(payload).hexdigest())
    res = V.fetch_mesh_file("Shinnecock_coarse.14", tmp_path, registry=reg, opener=lambda req, timeout=None: _FakeResp(payload))
    assert res["status"] == "ok"
    assert (tmp_path / "Shinnecock_coarse.14").read_bytes() == payload
    assert not (tmp_path / "Shinnecock_coarse.14.14").exists()
