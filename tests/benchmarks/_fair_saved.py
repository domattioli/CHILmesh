"""Saved synthetic campaigns; no benchmark process or privileged operations."""
import copy
import hashlib
import json
from benchmarks.competitors import campaign as C, records as R
from benchmarks.competitors.w2.tools import FAIR_TOOLS, FAIR_OPERATIONS, operation_contract
from test_records import successful_fair_data


def ruling(name='STOFS', version='v1'):
    return dict(mesh=name, exact_version=version, publish_results=True,
                citation='https://github.com/Valence/Valence/issues/123', ruling_date='2026-11-10',
                scope='publication of benchmark results', evidence_sha256='a'*64)


def saved(tmp_path, meshes=('STOFS',), operations=('parse',), tools=FAIR_TOOLS):
    operations=[(op,operation_contract('chilmesh',op)['backend_selected']) for op in operations]
    manifest=C.make_manifest('saved',meshes,operations,tools)
    from _grids import FAIR_FORT14
    from benchmarks.competitors.reference import build_reference
    rows=[]; references={}; mesh_refs={}
    for mesh in meshes:
        raw=tmp_path/(mesh+'.14');raw.write_text(mesh+'\n'+FAIR_FORT14.split('\n',1)[1])
        ref=build_reference(raw,tmp_path/('ref-'+mesh),name=mesh,version='v1')
        mesh_refs[mesh]=ref;references[ref['reference_id']]=ref
    for entry in manifest['schedule']:
        contract=operation_contract(entry['tool'],entry['operation'])
        ref=mesh_refs[entry['mesh']];reference_id=ref['reference_id']
        row=successful_fair_data()
        row.update({k:entry[k] for k in ('cell_id','series_id','repeat_id','repeat_index','execution_order','tool','operation','backend')})
        row.update(campaign_id=manifest['campaign_id'],manifest_sha256=manifest['manifest_sha256'],reference_id=reference_id,
                   contract_id=contract['contract_id'],equivalence=contract['equivalence'],input_format=contract['input_format'],input_sha256=ref['raw_file_sha256'])
        expected=dict(n_nodes=ref['header_node_count'],n_elems=ref['header_element_count'],node_sha256=ref['canonical_node_sha256'],element_sha256=ref['canonical_element_sha256'])
        row['output_verification'].update(expected,expected=expected.copy())
        if contract['equivalence']=='not-comparable':
            from test_records import fair_data
            blank=fair_data()
            for k in R.EVIDENCE_FIELDS: row[k]=blank[k]
            row.update(outcome='not-comparable',reason=contract['reason'],verification='not-applicable',measurement_validity=blank['measurement_validity'])
        rows.append(R.make_fair_record(**row))
    return manifest,rows,references


def persist(path, manifest, rows, references):
    evidence=path/'evidence';evidence.mkdir(parents=True)
    (evidence/'campaign.json').write_text(json.dumps(manifest))
    raw=(json.dumps(references,sort_keys=True,allow_nan=False)+'\n').encode()
    (evidence/'references.json').write_bytes(raw)
    (evidence/'reference-hashes.json').write_text(json.dumps({'references.json':hashlib.sha256(raw).hexdigest()}))
    paths=[]
    for row in rows:
        for section,key,payload in [('measurement','samples',json.dumps(row['measurement']['samples']).encode()),
            ('termination','stdout',str(row['termination']['stdout_tail']).encode()),
            ('termination','stderr',str(row['termination']['stderr_tail']).encode())]:
            if row[section][key+'_path'] is not None:
                sidecar=path/'raw'/row['repeat_id']/(key+'.json' if key=='samples' else key+'.txt')
                sidecar.parent.mkdir(parents=True,exist_ok=True);sidecar.write_bytes(payload)
                row[section][key+'_path']=str(sidecar.relative_to(path))
                row[section][key+'_sha256']=hashlib.sha256(payload).hexdigest()
        p=path/'attempts'/(row['repeat_id']+'.json');R.write_fair_record(row,p);paths.append(p)
    R.write_record_hashes(paths,evidence/'record-hashes.json')


def proposal(a,b):
    return dict(claim_id='claim',template='less',numerator_cell_id=a,denominator_cell_id=b,
                metric='rss_growth_bytes',unit='bytes',direction='lower',mesh='STOFS',exact_version='v1',
                operation='parse',formula='100*(B-A)/B',resource_condition_reference='frozen campaign')


def claim_cells():
    from benchmarks.competitors.aggregate import aggregate_series
    from test_aggregate import attempts
    cells={}
    for name,value in [('a',40),('b',50)]:
        rows=attempts(3)
        for r in rows: r['measurement']['rss_growth_bytes']=value
        cells[name]=dict(aggregate_series(rows,mesh_class='stofs'),cell_id=name,tool=name,tool_version='1',
            mesh='STOFS',exact_version='v1',operation='parse',contract_id='contract',backend='default',
            input_format='fort.14',input_sha256='1'*64,reference_id='ref',equivalence='like-for-like',
            verification='verified',conditions_valid=True,rights_cleared=True,resource_condition_reference='frozen campaign',
            series_id='original',source_record_hashes=['a'*64],manifest_sha256='b'*64,host_id='host',host_allocation={'cpu':1},
            environment_id=name,thread_environment=C.THREAD_ENV,options={})
    return cells
