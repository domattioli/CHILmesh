import copy
import pytest
from _fair_saved import saved,persist,ruling,proposal
from benchmarks.competitors import campaign as C


def test_full_synthetic_matrix_saved_hashes_and_redaction(tmp_path):
    from benchmarks.competitors.report.render import fair_report
    from benchmarks.competitors.w2.tools import FAIR_OPERATIONS
    m,rows,refs=saved(tmp_path,meshes=('Shinnecock','WNAT_Hagen','ENPAC2003','STOFS'),operations=FAIR_OPERATIONS)
    persist(tmp_path,m,rows,refs)
    result=fair_report(tmp_path,published=True,rights={'STOFS':ruling()})
    assert len(result['cells'])==6*len(FAIR_OPERATIONS)
    assert len(result['excluded'])==3
    cells=list(result['cells'].values());a=next(c for c in cells if c['tool']=='chilmesh' and c['operation']=='parse');b=next(c for c in cells if c['tool']=='adcircpy' and c['operation']=='parse')
    p=proposal(a['cell_id'],b['cell_id']);p['resource_condition_reference']=m['manifest_sha256']
    result=fair_report(tmp_path,published=True,rights={'STOFS':ruling()},proposals=[p])
    assert result['claims'][0]['accepted'] and result['claims'][0]['value']==0
    assert result['claims'][0]['source_record_hashes']
    p['template']='40–60% less memory'
    assert not fair_report(tmp_path,published=True,rights={'STOFS':ruling()},proposals=[p])['claims'][0]['accepted']


def test_controller_abort_preserves_active_and_every_unstarted(tmp_path):
    from test_campaign import ready_manifest,runtime
    m=ready_manifest()
    def abort(entry): raise KeyboardInterrupt()
    rows=C.execute_schedule(m,tmp_path,execute=abort,host_manifest=m['host_manifest'],runtime_manifest=runtime(),
        lock_check=lambda e:True,cgroup_check=lambda p:{'memory_events':{'oom_kill':0}},lock_path=tmp_path/'lock')
    assert len(rows)==3 and rows[0]['outcome']=='failed'
    assert all(r['outcome']=='not-run' for r in rows[1:])
    from benchmarks.competitors.records import read_record_hashes
    assert read_record_hashes(tmp_path/'evidence'/'record-hashes.json')==rows


@pytest.mark.parametrize('gate',['cache','load','runtime','memory','time'])
def test_invalid_conditions_stay_descriptive_never_claims(tmp_path,gate):
    from benchmarks.competitors.report.render import fair_report
    m,rows,refs=saved(tmp_path,tools=('chilmesh','adcircpy'))
    rows[0].update(outcome='invalid',reason=gate+' invalid')
    rows[0]['measurement_validity'][gate]={'valid':False,'reason':'synthetic invalidation'}
    persist(tmp_path,m,rows,refs)
    r=fair_report(tmp_path,published=True,rights={'STOFS':ruling()})
    cell=r['cells'][rows[0]['cell_id']]
    assert not cell['complete'] and not cell['conditions_valid']
    assert cell['metrics']['rss_growth_bytes']['n']==2


@pytest.mark.parametrize('kind',['empty','same-count'])
def test_corrupt_output_saved_series_rejected_by_real_gates(tmp_path,fair_fixture,kind):
    from benchmarks.competitors.reference import build_reference
    from benchmarks.competitors.verification import verify_output
    from benchmarks.competitors.report.render import fair_report
    ref=build_reference(fair_fixture['mesh'],fair_fixture['output'])
    output=copy.deepcopy(fair_fixture['adapter']);output['mapping_provenance']={'kind':'original IDs'}
    if kind=='empty': output.update(n_nodes=0,n_elems=0,node_ids=[],element_ids=[],coordinates=[],connectivity=[])
    else: output['coordinates'][1][0]+=1
    proof=verify_output(output,ref)
    assert proof['verification_state']=='failed'
    m,rows,refs=saved(tmp_path,tools=('chilmesh','adcircpy'))
    rows[0].update(outcome='invalid',verification='failed',reason='actual output corrupt',output_verification=proof)
    persist(tmp_path,m,rows,refs)
    cells=fair_report(tmp_path,published=True,rights={'STOFS':ruling()})['cells']
    p=proposal(rows[0]['cell_id'],rows[1]['cell_id']);p['resource_condition_reference']=m['manifest_sha256']
    from benchmarks.competitors.claims import evaluate_claim
    result=evaluate_claim(p,cells)
    assert not result['accepted'] and {'verification','complete'} <= set(result['failed_predicates'])


@pytest.mark.parametrize('values,complete,stable', [([95,100,105],True,True),([94,100,105],True,False),([100,100],False,True)])
def test_saved_repeat_spread_boundaries(tmp_path,values,complete,stable):
    from benchmarks.competitors.report.render import fair_report
    from test_records import fair_data
    m,rows,refs=saved(tmp_path,tools=('chilmesh',))
    for i,row in enumerate(rows):
        if i>=len(values):
            blank=fair_data()
            for key in ('measurement','protocol','termination','cache','load','output_verification','measurement_validity'): row[key]=blank[key]
            row.update(outcome='not-run',reason='missing tool',verification='incomplete')
        else:
            # Change growth while preserving sample-backed peak and baseline.
            row['measurement'].update(rss_growth_bytes=values[i],rss_operation_peak_bytes=100+values[i],rss_whole_process_peak_bytes=100+values[i])
            row['measurement']['samples'][1]['rss_bytes']=100+values[i]
    persist(tmp_path,m,rows,refs)
    cell=next(iter(fair_report(tmp_path)['cells'].values()))
    assert cell['complete']==complete and cell['stable']==stable


@pytest.mark.parametrize('artifact',['references','samples','stdout','stderr'])
def test_saved_reference_and_sidecar_hashes_checked_even_when_excluded(tmp_path,artifact):
    from benchmarks.competitors.report.render import fair_report
    from benchmarks.competitors.records import RecordError
    m,rows,refs=saved(tmp_path,tools=('chilmesh',));persist(tmp_path,m,rows,refs)
    if artifact=='references': path=tmp_path/'evidence/references.json'
    elif artifact=='samples': path=tmp_path/rows[0]['measurement']['samples_path']
    else: path=tmp_path/rows[0]['termination'][artifact+'_path']
    path.write_bytes(path.read_bytes()+b' tamper')
    with pytest.raises(RecordError): fair_report(tmp_path,published=True)


def test_full_matrix_stofs_termination_inventory_and_separate_conversion(tmp_path):
    from benchmarks.competitors.report.render import fair_report
    from benchmarks.competitors.campaign import validate_conversion
    from benchmarks.competitors.w2.tools import FAIR_OPERATIONS
    from test_records import fair_data
    m,rows,refs=saved(tmp_path,meshes=('Shinnecock','WNAT_Hagen','ENPAC2003','STOFS'),operations=FAIR_OPERATIONS)
    for row in rows:
        ref=refs[row['reference_id']]
        if row['operation']=='parse-netcdf' and row['tool'] in ('thalassa','xarray'):
            lineage=dict(source_reference_id=ref['reference_id'],source_sha256=ref['raw_file_sha256'],converted_sha256='2'*64,
                converter_revision='fair-netcdf-v1',converter_sha256='3'*64,netcdf4_version='1.0',native_versions={'netcdf':'1','hdf5':'1'},
                environment_lock_sha256='4'*64,argv=['synthetic-converter'],options={'reader':row['tool']},reader=row['tool'],
                source_to_row_maps={'nodes':[10,20,40,80],'elements':[7,9]},source_to_converted_id_map={'independently_reopened':True},
                source_content_verification=dict(row['output_verification'],verification_state='verified'),conversion_s=2,read_s=row['measurement']['operation_elapsed_s'],
                contract=dict(dimensions={'node':4,'nele':2,'max_face_nodes' if row['tool']=='xarray' else 'nvertex':4},
                    variables={'node_id':'i8','x':'f8','y':'f8','depth':'f8','element_id':'i8','arity':'i8',
                               'connectivity' if row['tool']=='xarray' else 'element':'i8'},
                    fill_value=0,indexing='original-file-ids' if row['tool']=='xarray' else 'one-based-row-index'))
            row['conversion']=validate_conversion(lineage,ref)
            assert row['conversion']['verification_state']=='verified'
            row['input_sha256']='2'*64
        if ref['name']!='STOFS': continue
        triggers={('chilmesh','parse'):'timeout',('adcircpy','parse'):'kernel-oom',
                  ('thalassa','parse-netcdf'):'rss-limit',('xarray','parse-netcdf'):'signal-unknown'}
        trigger=triggers.get((row['tool'],row['operation']))
        if trigger:
            row.update(outcome='failed',reason=trigger)
            row['termination'].update(trigger=trigger,returncode=-9,signal=9,limits={'wall_s':1800,'rss_bytes':48*1024**3})
            row['termination']['missing'].pop('signal',None)
            if trigger=='kernel-oom':
                row['termination']['memory_events_delta']['oom_kill']=1;row['termination']['memory_events_after']['oom_kill']=1
        if row['tool']=='chilmesh' and row['operation']=='write':
            blank=fair_data()
            for k in ('measurement','protocol','termination','cache','load','output_verification','measurement_validity'): row[k]=blank[k]
            row.update(outcome='not-run',reason='supported tool unavailable',verification='incomplete')
    persist(tmp_path,m,rows,refs)
    result=fair_report(tmp_path,published=True,rights={'STOFS':ruling()})
    for term in ['timeout','rss-limit','kernel-oom','signal-unknown','not-run','not-comparable','conversion_s','read_s','end_to_end_s','2.01']:
        assert term in result['markdown']
    for tool in m['tools']: assert tool in result['markdown']
    assert len(result['cells'])==60 and len(result['excluded'])==3


def test_active_abort_uses_sampler_saved_partial_evidence(tmp_path,fair_fixture):
    import json
    from benchmarks.competitors.w2.run import measure_fair_attempt
    from benchmarks.competitors.reference import build_reference
    from test_records import successful_fair_data
    from test_campaign import runtime
    reference=build_reference(fair_fixture['mesh'],fair_fixture['output'])
    identity=successful_fair_data();identity['tool_version']='1'
    partial=successful_fair_data()
    partial.update(status='failed',reason='controller abort',returncode=-9,adapter_result={})
    partial['termination'].update(returncode=-9,signal=9,trigger='signal-unknown')
    partial['termination']['missing'].pop('signal',None)
    (tmp_path/'controller.json').write_text(json.dumps(partial))
    def interrupted(*a,**k): raise KeyboardInterrupt()
    row=measure_fair_attempt(identity,reference,['synthetic'],mesh_class='small',expected_runtime=runtime(),artifact_dir=tmp_path,
        launcher=interrupted,cache_reset=lambda:partial['cache'],load_evidence=lambda:partial['load'])
    assert row['outcome']=='failed' and row['reason']=='controller abort'
    assert row['measurement']['last_rss_bytes']==140
    assert row['termination']['trigger']=='signal-unknown' and row['termination']['returncode']==-9


def test_abort_finalization_preserves_prepared_cell_identity_after_success(tmp_path,fair_fixture,monkeypatch):
    from benchmarks.competitors.w2 import run,tools
    from benchmarks.competitors.reference import build_reference
    from benchmarks.competitors.records import check_fair_coverage
    from benchmarks.competitors.aggregate import aggregate_series
    from test_campaign import ready_manifest,runtime,locked_environment
    from test_records import successful_fair_data
    base=ready_manifest();env=locked_environment();env.update(environment_id='chilmesh',tool_distribution='chilmesh',
        distributions={'chilmesh':'1','dependency':'2'},direct_dependencies=['chilmesh'])
    m=C.make_manifest('abort',['STOFS'],['parse'],['chilmesh'],
        **{k:base[k] for k in C.MANIFEST_EVIDENCE if k!='environment_manifests'},environment_manifests=[env])
    ref=build_reference(fair_fixture['mesh'],fair_fixture['output'],name='STOFS',version='v1')
    identities={};commands={}
    for e in m['schedule']:
        row=successful_fair_data();row.update({k:e[k] for k in ('cell_id','series_id','repeat_id','repeat_index','execution_order','tool','operation','backend')})
        row.update(campaign_id=m['campaign_id'],manifest_sha256=m['manifest_sha256'],reference_id=ref['reference_id'],
            contract_id=tools.operation_contract('chilmesh','parse')['contract_id'],tool_version='1',environment_id='chilmesh',
            host_allocation={k:base['host_manifest'][k] for k in ('affinity','cpu_quota','cpuset','memory_limits')},argv=['/venv/bin/python','parse'])
        identities[e['repeat_id']]=row;commands[e['repeat_id']]=row['argv']
    calls=[]
    def launch(cmd,**kw):
        calls.append(cmd)
        if len(calls)>1: raise KeyboardInterrupt()
        row=successful_fair_data();row.update(status='ok',reason='',returncode=0,adapter_result=dict(version='1',runtime=runtime(),
            new_modules=[],materialization=True,output={**fair_fixture['adapter'],'mapping_provenance':{'kind':'IDs'}}))
        return row
    def load():
        load=successful_fair_data()['load'];load['allocation']=next(iter(identities.values()))['host_allocation'];return load
    rows=run.run_fair_campaign(m,identities,{'STOFS':ref},commands,tmp_path/'run',host_probe=lambda:m['host_manifest'],runtime_probe=runtime,
        launcher=launch,cache_reset=lambda:successful_fair_data()['cache'],load_evidence=load,
        lock_check=lambda e:True,cgroup_check=lambda p:{'memory_events':{'oom_kill':0}},lock_path=tmp_path/'lock')
    assert [row['outcome'] for row in rows]==['ok','failed','not-run']
    check_fair_coverage(rows,[])
    summary=aggregate_series(rows,mesh_class='stofs')
    assert not summary['complete'] and summary['metrics']['rss_growth_bytes']['n']==1


def test_campaign_conversion_transport_binds_each_repeat_source(tmp_path,monkeypatch):
    from test_campaign import ready_manifest,locked_environment,runtime
    from benchmarks.competitors.w2 import run
    base=ready_manifest();env=locked_environment();env.update(environment_id='thalassa')
    m=C.make_manifest('conversion-transport',['STOFS'],['parse-netcdf'],['thalassa'],
        **{k:base[k] for k in C.MANIFEST_EVIDENCE if k!='environment_manifests'},environment_manifests=[env])
    identities={};commands={};conversions={};seen=[]
    for e in m['schedule']:
        row=C.unstarted_record(m,e,'prepared synthetic input')
        row.update(environment_id='thalassa',tool_version='1',host_allocation={k:base['host_manifest'][k] for k in ('affinity','cpu_quota','cpuset','memory_limits')})
        identities[e['repeat_id']]=row;commands[e['repeat_id']]=['/venv/bin/python','synthetic']
        conversions[e['repeat_id']]={'source_reference_id':'source-'+str(e['repeat_index'])}
    def measured(identity,reference,command,**kw):
        seen.append(kw['conversion']);return identity
    monkeypatch.setattr(run,'measure_fair_attempt',measured)
    run.run_fair_campaign(m,identities,{'STOFS':{'reference_id':'synthetic'}},commands,tmp_path/'run',conversions=conversions,
        host_probe=lambda:m['host_manifest'],runtime_probe=runtime,lock_check=lambda e:True,
        cgroup_check=lambda p:{'memory_events':{'oom_kill':0}},lock_path=tmp_path/'lock')
    assert seen==[conversions[e['repeat_id']] for e in m['schedule']]
