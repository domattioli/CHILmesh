"""Renderer: tables per workflow and track, failed and not-run rows keep their reasons."""
from __future__ import annotations

import pytest

from benchmarks.competitors import records as R
from benchmarks.competitors.report import render


def _rec(tool, status="ok", reason="", metrics=None, workflow="w2", track=None, cache="cold", mesh="Synth"):
    return R.make_record(
        tool=tool, tool_version="1", license="MIT", run_mode="subprocess",
        mesh=dict(name=mesh, version="v1", sha256="0" * 64, provenance="p", crs="cartesian",
                  units="m", n_nodes=3, n_elems=1),
        workflow=workflow, track=track, metric_values=metrics or {}, status=status, reason=reason,
        cache=cache, command=["x"], wall_s=1.0, peak_rss_bytes=1, repeat_index=0,
        machine={}, chilmesh_git_sha="abc")


def test_table_has_row_per_tool_and_column_per_metric():
    recs = [_rec("a", metrics={"parse_s": R.metric(1.25), "adj_s": R.metric(2.5)}),
            _rec("b", metrics={"parse_s": R.metric(3.0),
                               "adj_s": R.metric(None, na=True, reason="no adjacency API")})]
    md = render.render_markdown(recs, tools=["a", "b"], metrics=["parse_s", "adj_s"])
    assert "| a |" in md and "| b |" in md
    assert "parse_s" in md and "adj_s" in md
    assert "1.25" in md
    assert "not applicable" in md and "no adjacency API" in md


def test_failed_and_not_run_rows_show_reason():
    recs = [_rec("a", metrics={"parse_s": R.metric(1.0)}),
            _rec("b", status="failed", reason="OOM-killed (SIGKILL)"),
            _rec("c", status="not-run", reason="mesh not yet served by Valence")]
    md = render.render_markdown(recs, tools=["a", "b", "c"], metrics=["parse_s"])
    assert "failed" in md and "OOM-killed (SIGKILL)" in md
    assert "not run" in md and "mesh not yet served by Valence" in md


def test_missing_cell_fails_loudly():
    recs = [_rec("a", metrics={"parse_s": R.metric(1.0)})]
    with pytest.raises(R.MissingCellError):
        render.render_markdown(recs, tools=["a"], metrics=["parse_s", "adj_s"])


def test_separate_sections_per_workflow_and_track():
    recs = [_rec("a", workflow="w1", track=1, metrics={"q": R.metric(0.9)}),
            _rec("a", workflow="w1", track=2, metrics={"q": R.metric(0.8)})]
    md = render.render_sections(recs, tools=["a"], metrics=["q"])
    assert "w1 / track 1" in md and "w1 / track 2" in md


def test_non_cold_runs_are_flagged_and_excluded_from_published_table():
    recs = [_rec("a", cache="warm", metrics={"parse_s": R.metric(1.0)})]
    md = render.render_markdown(recs, tools=["a"], metrics=["parse_s"], published=True)
    assert "warm" in md and "excluded" in md.lower()


def test_invalid_status_is_shown():
    recs = [_rec("a", status="invalid", reason="boundary moved", metrics={"parse_s": R.metric(1.0)})]
    md = render.render_markdown(recs, tools=["a"], metrics=["parse_s"])
    assert "invalid" in md and "boundary moved" in md


def test_w2_metrics_include_station_max_abs_diff_and_render_covers_it():
    from benchmarks.competitors.w2 import run, tools

    assert "station_max_abs_diff" in tools.METRICS
    na = run._na_secondary("adcircpy")
    assert na["station_max_abs_diff"]["na"] and na["station_max_abs_diff"]["reason"]
    full = {m: R.metric(1.0) for m in tools.METRICS if m != "station_max_abs_diff"}
    recs = [_rec("chilmesh", metrics=full),
            _rec("chilmesh", metrics={"station_max_abs_diff": R.metric(0.25)}),
            _rec("adcircpy", metrics={"parse_s": R.metric(1.0), **na})]
    md = render.render_markdown(recs, tools=["chilmesh", "adcircpy"], metrics=list(tools.METRICS))
    assert "station_max_abs_diff" in md and "0.25" in md


def test_w2_station_not_run_reason_shown_in_table():
    from benchmarks.competitors.w2 import tools

    full = {m: R.metric(1.0) for m in tools.METRICS if m != "station_max_abs_diff"}
    recs = [_rec("chilmesh", metrics=full),
            _rec("chilmesh", status="not-run", reason="NOAA inputs not fetched (maxele)")]
    md = render.render_markdown(recs, tools=["chilmesh"], metrics=list(tools.METRICS))
    assert "NOAA inputs not fetched (maxele)" in md


def test_exclusions_listed():
    md = render.render_markdown([_rec("a", metrics={"m": R.metric(1.0)})], tools=["a"], metrics=["m"],
                                exclusions=[{"mesh": "WNAT_Hagen", "reason": "no rights entry"}])
    assert "WNAT_Hagen" in md and "no rights entry" in md


def test_fair_report_complete_grouped_hashed_and_rights_first(tmp_path):
    from _fair_saved import saved, persist, ruling
    m,rows,refs=saved(tmp_path,meshes=('STOFS','WNAT_Hagen'),operations=('parse','parse-netcdf','write'))
    # Even a diagnostic must never leak an excluded mesh's numeric evidence.
    for row in rows:
        if refs[row['reference_id']]['name']=='WNAT_Hagen':
            row['termination']['stdout_tail']='SECRET_NUMBER_87654321'
            row['termination'].get('missing',{}).pop('stdout_tail',None)
    persist(tmp_path,m,rows,refs)
    result=render.fair_report(tmp_path,published=True,rights={'STOFS':ruling()})
    text=result['markdown']
    for needle in ['Native like-for-like','conversion / NetCDF','CHILmesh-only','not-comparable','STOFS','123','sha256','raw','median','growth','10 ms']:
        assert needle in text
    assert 'SECRET_NUMBER_87654321' not in str(result)
    assert len(result['cells'])==18
    assert result['excluded']==[{'mesh':'WNAT_Hagen','reason':'publication rights unavailable or invalid'}]
    assert all(c['mesh']=='STOFS' for c in result['cells'].values())


def test_fair_report_missing_attempt_or_tamper_fails_before_rights(tmp_path):
    from _fair_saved import saved, persist
    m,rows,refs=saved(tmp_path)
    persist(tmp_path,m,rows[:-1],refs)
    with pytest.raises(R.MissingCellError): render.fair_report(tmp_path,published=True)
    p=next((tmp_path/'attempts').glob('*.json'));p.write_text(p.read_text()+' ')
    with pytest.raises(R.RecordError): render.fair_report(tmp_path)


def test_fair_report_failures_and_incomplete_descriptive_series(tmp_path):
    from _fair_saved import saved, persist, ruling
    m,rows,refs=saved(tmp_path)
    triggers=['timeout','rss-limit','kernel-oom','signal-unknown']
    for i,row in enumerate(rows):
        if row['outcome']=='ok':
            row.update(outcome='failed',reason=triggers[i%4])
            row['termination'].update(trigger=triggers[i%4],returncode=-9,signal=9)
            row['termination']['missing'].pop('signal',None)
            if triggers[i%4]=='kernel-oom':
                row['termination']['memory_events_delta']['oom_kill']=1;row['termination']['memory_events_after']['oom_kill']=1
    persist(tmp_path,m,rows,refs)
    r=render.fair_report(tmp_path,published=True,rights={'STOFS':ruling()})
    for tool in m['tools']: assert tool in r['markdown']
    assert 'incomplete' in r['markdown'] and 'timeout' in r['markdown'] and 'limits' in r['markdown']


def test_public_claim_diagnostics_contain_no_excluded_numbers_hashes_or_free_text(tmp_path):
    from _fair_saved import saved,persist,ruling,proposal
    import json
    m,rows,refs=saved(tmp_path,meshes=('STOFS','WNAT_Hagen'),tools=('chilmesh','adcircpy'))
    persist(tmp_path,m,rows,refs)
    private=[r for r in rows if refs[r['reference_id']]['name']=='WNAT_Hagen']
    p=proposal(private[0]['cell_id'],private[1]['cell_id'])
    p.update(claim_id='PRIVATE_VALUE_1234567',sentence='PRIVATE_VALUE_1234567',mesh='WNAT_Hagen')
    r=render.fair_report(tmp_path,published=True,rights={'STOFS':ruling()},proposals=[p])
    text=json.dumps(r)
    assert 'PRIVATE_VALUE_1234567' not in text
    index=json.loads((tmp_path/'evidence/record-hashes.json').read_text())
    for entry in index['records']:
        if entry['repeat_id'] in {row['repeat_id'] for row in private}: assert entry['sha256'] not in text
    assert not r['publication_passed']
