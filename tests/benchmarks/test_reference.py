"""Independent canonical goldens; no timed reader supplies expected values."""
import ast
import hashlib
from pathlib import Path

import pytest


def build(fair_fixture, **kw):
    from benchmarks.competitors.reference import build_reference
    return build_reference(fair_fixture['mesh'], fair_fixture['output'], chunk_size=1, **kw)


def test_reference_golden_chunks_and_binding(fair_fixture):
    r = build(fair_fixture, name='synthetic', version='v1', served_sha256=hashlib.sha256(fair_fixture['mesh'].read_bytes()).hexdigest())
    nodes = b'10 0.0 -0.0 1.0\n20 1.0 0.0 2.0\n40 1.0 1.0 3.0\n80 0.0 1.0 4.0'
    elems = b'7 3 10 20 40\n9 4 10 20 40 80'
    assert Path(r['nodes_path']).read_bytes() == nodes
    assert Path(r['elements_path']).read_bytes() == elems
    assert r['canonical_node_sha256'] == hashlib.sha256(nodes).hexdigest()
    assert r['canonical_element_sha256'] == hashlib.sha256(elems).hexdigest()
    assert r['raw_file_sha256'] == hashlib.sha256(fair_fixture['mesh'].read_bytes()).hexdigest()
    assert (r['header_node_count'], r['header_element_count']) == (4, 2)
    assert r['verification_state'] == 'verified'
    assert r['sort_chunks'] >= 6


@pytest.mark.parametrize('before,after', [('0.0', '0e0'), ('10 0.0', '10 0.000000000000000000'), ('20 1.0', '20 1e0')])
def test_decimal_spellings(fair_fixture, before, after):
    path = fair_fixture['mesh']; path.write_text(path.read_text().replace(before, after, 1))
    assert build(fair_fixture)['header_node_count'] == 4


@pytest.mark.parametrize('before,after', [('20 1', '10 1'), ('7 3 10 20 40', '7 3 10 20 99'), ('7 3 10 20 40', '7 5 10 20 40'), ('20 1', '20 nan'), ('20 1', '20 inf'), ('2 4', '2 5')])
def test_reject_bad_content(fair_fixture, before, after):
    from benchmarks.competitors.reference import ReferenceError
    p = fair_fixture['mesh']; p.write_text(p.read_text().replace(before, after, 1))
    with pytest.raises(ReferenceError): build(fair_fixture)


def test_precision_and_digest_mismatch(fair_fixture):
    from benchmarks.competitors.reference import ReferenceError
    with pytest.raises(ReferenceError, match='digest'): build(fair_fixture, served_sha256='0' * 64)
    p = fair_fixture['mesh']; p.write_text(p.read_text().replace('20 1.0', '20 1.00000000000000001', 1))
    r = build(fair_fixture)
    assert r['precision_token_count_gt17'] == 1
    assert r['verification_state'] == 'incomplete'


def test_reorder_and_winding(fair_fixture, tmp_path):
    r = build(fair_fixture)
    p = fair_fixture['mesh']; rows = p.read_text().splitlines(); rows[2:6] = reversed(rows[2:6]); p.write_text('\n'.join(rows) + '\n')
    from benchmarks.competitors.reference import build_reference
    other = build_reference(p, tmp_path / 'other', chunk_size=2)
    assert r['canonical_node_sha256'] == other['canonical_node_sha256']
    p.write_text(p.read_text().replace('7 3 10 20 40', '7 3 40 20 10'))
    changed = build_reference(p, tmp_path / 'changed')
    assert r['canonical_element_sha256'] != changed['canonical_element_sha256']


def test_reference_import_independence():
    from benchmarks.competitors import reference
    tree = ast.parse(Path(reference.__file__).read_text())
    imports = [n.module if isinstance(n, ast.ImportFrom) else a.name for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom)) for a in (n.names if isinstance(n, ast.Import) else [None])]
    assert not any(x and any(t in x for t in ('chilmesh', 'adcircpy', 'thalassa', 'numpy', 'xarray')) for x in imports)
