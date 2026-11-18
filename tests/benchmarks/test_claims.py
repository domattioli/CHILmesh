import copy
import pytest
from _fair_saved import claim_cells, proposal


def test_claim_exact_arithmetic_and_bounded_evidence():
    from benchmarks.competitors.claims import evaluate_claim
    result=evaluate_claim(proposal('a','b'),claim_cells())
    assert result['accepted'] and result['value']==20 and result['failed_predicates']==[]
    for text in ['STOFS','v1','rss_growth_bytes','20','a','b','original','1','n=3','spread','frozen campaign']:
        assert text in result['sentence']
    assert result['source_record_hashes']==['a'*64]


@pytest.mark.parametrize('field,value,predicate', [
 ('equivalence','CHILmesh-only','equivalence'),('verification','failed','verification'),
 ('complete',False,'complete'),('stable',False,'stable'),('conditions_valid',False,'conditions'),
 ('rights_cleared',False,'rights'),('input_format','NetCDF','input'),('input_sha256','other','input'),
 ('contract_id','other','contract'),('operation','write','operation'),('exact_version','v2','mesh'),
 ('manifest_sha256','other','resources'),('source_record_hashes',[],'evidence')])
def test_claim_rejects_every_cell_gate(field,value,predicate):
    from benchmarks.competitors.claims import evaluate_claim
    cells=claim_cells();cells['a'][field]=value
    r=evaluate_claim(proposal('a','b'),cells)
    assert not r['accepted'] and predicate in r['failed_predicates'] and r['sentence'] is None


@pytest.mark.parametrize('change,predicate', [({'metric':'rss_whole_process_peak_bytes','unit':'seconds'},'metric'),
 ({'denominator_metric':'rss_whole_process_peak_bytes'},'metric'),({'direction':'higher'},'direction'),
 ({'formula':'B/A'},'formula'),({'template':'40–60% less memory everywhere'},'template'),
 ({'resource_condition_reference':''},'resources'),({'denominator_cell_id':'missing'},'cells'),
 ({'sentence':'CHILmesh is always faster'},'template')])
def test_claim_rejects_unsupported_proposals(change,predicate):
    from benchmarks.competitors.claims import evaluate_claim
    r=evaluate_claim(dict(proposal('a','b'),**change),claim_cells())
    assert not r['accepted'] and predicate in r['failed_predicates']


def test_zero_denominator_and_all_metric_instability_reject():
    from benchmarks.competitors.claims import evaluate_claim
    cells=claim_cells();cells['b']['metrics']['rss_growth_bytes']['median']=0
    assert 'denominator' in evaluate_claim(proposal('a','b'),cells)['failed_predicates']
    cells=claim_cells();cells['b']['metrics']['operation_elapsed_s']['unstable']=True
    assert 'stable' in evaluate_claim(proposal('a','b'),cells)['failed_predicates']


@pytest.mark.parametrize('field,value,predicate',[('tool_version','unknown','versions'),('tool_version','>=1','versions'),
    ('source_record_hashes',['not-a-hash'],'evidence'),('series_id','other','resources')])
def test_claim_exact_versions_hashes_and_series(field,value,predicate):
    from benchmarks.competitors.claims import evaluate_claim
    cells=claim_cells();cells['a'][field]=value
    result=evaluate_claim(proposal('a','b'),cells)
    assert not result['accepted'] and predicate in result['failed_predicates']
