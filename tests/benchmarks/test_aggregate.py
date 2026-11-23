import json
import pytest
from test_records import successful_fair_data


def attempts(n=5):
    rows=[]
    for i in range(n):
        r=successful_fair_data();r.update(repeat_id=str(i),repeat_index=i);rows.append(r)
    return rows


@pytest.mark.parametrize('values,spread,state,unstable', [([95,100,105],10,'finite',False),([94,100,105],11,'finite',True),([0,0,0],0,'finite',False),([0,0,1],None,'infinite',True)])
def test_metric_arithmetic_boundaries(values,spread,state,unstable):
    from benchmarks.competitors.aggregate import summarize_metric
    s=summarize_metric(list(zip(map(str,range(len(values))),values)),unit='bytes')
    assert s['n']==len(values) and s['min']==min(values) and s['max']==max(values)
    assert s['absolute_spread']==max(values)-min(values)
    assert s['relative_spread_pct']==spread and s['spread_state']==state and s['unstable']==unstable
    json.dumps(s,allow_nan=False)


@pytest.mark.parametrize('mesh,n',[('Shinnecock',5),('STOFS',3)])
def test_complete_counts_and_all_metric_stability(mesh,n):
    from benchmarks.competitors.aggregate import aggregate_series
    rows=attempts(n)
    s=aggregate_series(rows,mesh_class='stofs' if mesh=='STOFS' else 'small')
    assert s['complete'] and s['stable'] and len(s['metrics'])==3
    rows[-1]['measurement']['rss_growth_bytes']=999
    assert not aggregate_series(rows,mesh_class='stofs' if mesh=='STOFS' else 'small')['stable']
    assert not aggregate_series(rows[:-1],mesh_class='stofs' if mesh=='STOFS' else 'small')['complete']


def test_failures_remain_excluded_and_never_zero():
    from benchmarks.competitors.aggregate import aggregate_series
    rows=attempts();rows[2]['outcome']='failed';rows[2]['reason']='timeout';rows[2]['measurement']['operation_elapsed_s']=None
    s=aggregate_series(rows,mesh_class='small')
    assert not s['complete'] and s['excluded_attempts']==[{'repeat_id':'2','reason':'timeout'}]
    assert s['metrics']['operation_elapsed_s']['n']==4
    assert all(x['value'] is not None for x in s['metrics']['operation_elapsed_s']['raw_values'])


@pytest.mark.parametrize('key',['series_id','environment_id','reference_id','contract_id','backend','manifest_sha256','host_allocation'])
def test_no_identity_pooling(key):
    from benchmarks.competitors.aggregate import aggregate_series
    rows=attempts();rows[-1][key]='different'
    with pytest.raises(ValueError): aggregate_series(rows,mesh_class='small')


def test_replacement_lineage_requires_rationale():
    from benchmarks.competitors.aggregate import aggregate_series
    with pytest.raises(ValueError): aggregate_series(attempts(),mesh_class='small',previous_series_id='old')
    assert aggregate_series(attempts(),mesh_class='small',previous_series_id='old',replacement_rationale='complete rerun')['previous_series_id']=='old'


@pytest.mark.parametrize('key',['input_sha256','input_format','conversion','actual_backend'])
def test_no_input_or_actual_backend_pooling(key):
    from benchmarks.competitors.aggregate import aggregate_series
    rows=attempts()
    if key=='actual_backend':
        for r in rows: r['adapter_result']={'actual_backend':{'implementation':'python'}}
        rows[-1]['adapter_result']['actual_backend']={'implementation':'cpp'}
    else:
        for r in rows: r[key]='same'
        rows[-1][key]='different'
    with pytest.raises(ValueError): aggregate_series(rows,mesh_class='small')
