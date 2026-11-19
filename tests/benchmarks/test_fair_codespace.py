"""Stub executable only: no installations, meshes, network or cache reset."""
import os
import subprocess
from pathlib import Path
import pytest

SCRIPT=Path(__file__).resolve().parents[2]/'benchmarks/competitors/scripts/run_fair_codespace.sh'


@pytest.mark.parametrize('mode,optin,expected', [(None,False,0),('--preflight',False,0),('--run',False,2),('--run',True,0)])
def test_codespace_explicit_modes_and_base_python(tmp_path,mode,optin,expected):
    python=tmp_path/'python';log=tmp_path/'argv'
    python.write_text('#!/bin/sh\nprintf "%s\\n" "$@" >> "$CALL_LOG"\nexit 0\n');python.chmod(0o755)
    env={**os.environ,'CHILMESH_FAIR_PYTHON':str(python),'CALL_LOG':str(log)}
    env.pop('CHILMESH_RUN_BENCH',None)
    if optin: env['CHILMESH_RUN_BENCH']='1'
    argv=['bash',str(SCRIPT),'--manifest','frozen.json','--out',str(tmp_path/'out')]
    if mode: argv.append(mode)
    p=subprocess.run(argv,env=env,text=True,capture_output=True)
    assert p.returncode==expected,p.stderr
    if expected==0:
        args=log.read_text();assert 'fair-preflight' in args and 'frozen.json' in args
        assert ('fair-run' in args)==(mode=='--run')
    else: assert not log.exists()


def test_codespace_stops_on_failed_capability_preflight(tmp_path):
    python=tmp_path/'python';python.write_text('#!/bin/sh\nexit 1\n');python.chmod(0o755)
    p=subprocess.run(['bash',str(SCRIPT),'--run','--manifest','missing','--out',str(tmp_path/'out')],
        env={**os.environ,'CHILMESH_RUN_BENCH':'1','CHILMESH_FAIR_PYTHON':str(python)},capture_output=True)
    assert p.returncode==1


def test_codespace_serial_output_lock_and_evidence_archive(tmp_path):
    import tarfile
    python=tmp_path/'python';log=tmp_path/'argv';out=tmp_path/'out'
    python.write_text('#!/bin/sh\nprintf "%s\\n" "$@" >> "$CALL_LOG"\nif [ "$3" = fair-run ]; then\n mkdir -p "$7/evidence" "$7/attempts" "$7/raw"\n echo frozen > "$7/evidence/campaign.json"\n echo attempt > "$7/attempts/a.json"\n echo samples > "$7/raw/samples.json"\nfi\n')
    python.chmod(0o755)
    env={**os.environ,'CHILMESH_RUN_BENCH':'1','CHILMESH_FAIR_PYTHON':str(python),'CALL_LOG':str(log)}
    argv=['bash',str(SCRIPT),'--run','--manifest','immutable.json','--out',str(out)]
    (out/'.controller-lock').mkdir(parents=True)
    assert subprocess.run(argv,env=env,capture_output=True).returncode==1
    assert not log.exists()
    (out/'.controller-lock').rmdir()
    p=subprocess.run(argv,env=env,capture_output=True,text=True)
    assert p.returncode==0,p.stderr
    with tarfile.open(out/'evidence-archive.tar.gz') as archive:
        assert {'evidence/campaign.json','attempts/a.json','raw/samples.json'} <= set(archive.getnames())
    assert log.read_text().count('immutable.json')==2
    assert not (out/'.controller-lock').exists()
