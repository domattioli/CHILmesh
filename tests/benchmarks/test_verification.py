import copy
import pytest


def evidence(fair_fixture, output=None):
    from benchmarks.competitors.reference import build_reference
    from benchmarks.competitors.verification import verify_output
    r = build_reference(fair_fixture['mesh'], fair_fixture['output'])
    data = copy.deepcopy(fair_fixture['adapter'] if output is None else output)
    data.setdefault('mapping_provenance', {'kind': 'original-file-ids', 'evidence': 'synthetic API'})
    return verify_output(data, r)


def test_actual_materialized_output(fair_fixture):
    assert evidence(fair_fixture)['verification_state'] == 'verified'


@pytest.mark.parametrize('field,value', [('n_nodes', 0), ('n_elems', None), ('node_ids', [10, 10, 40, 80]), ('element_ids', [8, 9]), ('coordinates', [[0., 0., 1.], [1., 0., 2.], [1., 1., 3.], [0., 1., 4.]]), ('connectivity', [[40, 20, 10], [10, 20, 40, 80]])])
def test_corruption_and_missing(fair_fixture, field, value):
    out = copy.deepcopy(fair_fixture['adapter']); out[field] = value
    assert evidence(fair_fixture, out)['verification_state'] != 'verified'


@pytest.mark.parametrize('column', [0, 1, 2])
def test_coordinate_or_depth_corruption(fair_fixture, column):
    out = copy.deepcopy(fair_fixture['adapter']); out['coordinates'][1][column] += .1
    assert evidence(fair_fixture, out)['verification_state'] == 'failed'


def test_explicit_offset_map_and_no_inferred_map(fair_fixture):
    out = copy.deepcopy(fair_fixture['adapter']); out['connectivity'] = [[0, 1, 2], [0, 1, 2, 3]]
    out['connectivity_indexing'] = 'row-offset'; out['index_to_node_id'] = [10, 20, 40, 80]
    assert evidence(fair_fixture, out)['verification_state'] == 'verified'


def test_unavailable_map_incomplete(fair_fixture):
    out = copy.deepcopy(fair_fixture['adapter']); out['connectivity_indexing'] = 'row-offset'
    assert evidence(fair_fixture, out)['verification_state'] == 'incomplete'


def test_float32_loss(fair_fixture):
    import numpy as np
    p = fair_fixture['mesh']; p.write_text(p.read_text().replace('20 1.0', '20 1.00000001', 1))
    out = copy.deepcopy(fair_fixture['adapter']); out['coordinates'][1][0] = float(np.float32(1.00000001))
    assert evidence(fair_fixture, out)['verification_state'] == 'failed'


def test_verifier_never_reads_source(fair_fixture, monkeypatch):
    from benchmarks.competitors.reference import build_reference
    from benchmarks.competitors.verification import verify_output
    r = build_reference(fair_fixture['mesh'], fair_fixture['output']); fair_fixture['mesh'].unlink()
    out = {**fair_fixture['adapter'], 'mapping_provenance': {'kind': 'original-file-ids'}}
    assert verify_output(out, r)['verification_state'] == 'verified'


def test_written_artifact_checked_against_independent_reference(fair_fixture, tmp_path):
    from benchmarks.competitors.reference import build_reference
    from benchmarks.competitors.verification import verify_written_output
    reference = build_reference(fair_fixture['mesh'], fair_fixture['output'])
    written = tmp_path / 'written.14'; written.write_bytes(fair_fixture['mesh'].read_bytes())
    assert verify_written_output(written, reference)['verification_state'] == 'verified'
    written.write_text(written.read_text().replace('20 1.0', '20 2.0'))
    assert verify_written_output(written, reference)['verification_state'] == 'failed'
    written.unlink()
    assert verify_written_output(written, reference)['verification_state'] == 'incomplete'
