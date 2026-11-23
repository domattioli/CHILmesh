"""Bounded child lifecycle; no benchmark or external tool needed."""
import os
import sys
import time
import pytest
from benchmarks.competitors import measure as M


def test_limits_strict_exceeds():
    assert M.limit_trigger(10, 100, wall_s=10, rss_bytes=100) is None
    assert M.limit_trigger(10.001, 100, wall_s=10, rss_bytes=100) == 'timeout'
    assert M.limit_trigger(10, 101, wall_s=10, rss_bytes=100) == 'rss-limit'


@pytest.mark.parametrize('code', ['pass', 'raise SystemExit(3)'])
def test_missing_protocol_early_exit_is_bounded(code, tmp_path):
    r = M.run_measured([sys.executable,'-c',code], timeout_s=2, repeat_id='r', artifact_dir=tmp_path)
    assert r['status'] == 'failed' and 'protocol' in r['reason']
    assert r['termination']['reap_diagnostics']['reaped']
    assert r['protocol']['process_group_id'] == r['protocol']['pid']


def test_protocol_result_not_stdout_tail_and_verification_outside_region(tmp_path):
    code = '''from benchmarks.competitors.protocol import protocol_from_environment
import time
p=protocol_from_environment(); p.ready()
with p.region(): time.sleep(.06)
time.sleep(.03)
p.result(dict(payload='x'*20000, materialization=True))
'''
    r=M.run_measured([sys.executable,'-c',code],timeout_s=3,repeat_id='r',artifact_dir=tmp_path,tail_bytes=10)
    assert r['status']=='ok', r
    assert len(r['adapter_result']['payload'])==20000
    assert r['measurement']['operation_elapsed_s'] < r['wall_s']
    assert r['termination']['reap_diagnostics']['reaped']


@pytest.mark.skipif(sys.platform != 'linux', reason='Linux procfs lifecycle')
@pytest.mark.parametrize('trigger',['timeout','rss-limit','abort'])
def test_group_kill_reaps_child_and_kills_descendant(tmp_path, trigger):
    pidfile=tmp_path/'pid'
    code="import subprocess,time,os; p=subprocess.Popen(['sleep','10']); open(%r,'w').write(str(p.pid)); time.sleep(10)" % str(pidfile)
    def sample(pid, **kw):
        s=M.read_rss_sample(pid, **kw)
        if trigger=='abort' and pidfile.exists(): raise KeyboardInterrupt
        if trigger=='rss-limit': s['rss_bytes']=101 if pidfile.exists() and pidfile.read_text() else 0
        return s
    if trigger=='abort':
        with pytest.raises(KeyboardInterrupt):
            M.run_measured([sys.executable,'-c',code],timeout_s=.3,sample_reader=sample,artifact_dir=tmp_path)
    else:
        r=M.run_measured([sys.executable,'-c',code],timeout_s=.3,rss_limit_bytes=100 if trigger=='rss-limit' else None,sample_reader=sample,artifact_dir=tmp_path)
        assert r['termination']['trigger']==trigger
        assert r['termination']['reap_diagnostics']['reaped'] and r['measurement']['samples']
    child=int(pidfile.read_text())
    status='/proc/%s/status'%child
    # SIGKILL delivery to the reparented descendant is asynchronous; wait for it.
    def gone():
        try:
            with open(status) as fh: return '\nState:\tZ' in fh.read()
        except FileNotFoundError: return True
    deadline=time.monotonic()+5
    while not gone() and time.monotonic()<deadline: time.sleep(.02)
    assert gone()


def test_launch_error_has_explicit_partial_evidence(tmp_path):
    r=M.run_measured(['/nonexistent/us2-binary'],timeout_s=1,repeat_id='r',artifact_dir=tmp_path)
    assert r['termination']['trigger']=='launch-error'
    assert r['measurement']['rss_growth_bytes'] is None
    assert r['measurement']['missing']['rss_growth_bytes']


def test_malformed_protocol_kills_and_reaps(tmp_path):
    code="import os,time; os.write(int(os.environ['CHILMESH_CONTROL_OUT_FD']),b'not-json\\n'); time.sleep(5)"
    r=M.run_measured([sys.executable,'-c',code],timeout_s=2,repeat_id='r',artifact_dir=tmp_path)
    assert r['status']=='failed' and 'protocol' in r['reason']
    assert r['termination']['reap_diagnostics']['reaped'] and r['termination']['kill_diagnostics']


def test_controller_abort_closes_inherited_descriptors_and_reaps(tmp_path, monkeypatch):
    captured={};original=M.subprocess.Popen
    def popen(*args,**kwargs):
        proc=original(*args,**kwargs);captured.update(proc=proc,fds=list(kwargs['pass_fds']))
        assert kwargs['start_new_session'] is True
        return proc
    monkeypatch.setattr(M.subprocess,'Popen',popen)
    def abort(pid): raise KeyboardInterrupt
    with pytest.raises(KeyboardInterrupt):
        M.run_measured([sys.executable,'-c','import time;time.sleep(5)'],timeout_s=2,repeat_id='r',sample_reader=abort,artifact_dir=tmp_path)
    assert captured['proc'].returncode is not None
    for fd in captured['fds']:
        with pytest.raises(OSError): os.fstat(fd)


@pytest.mark.parametrize('delta,attributable,expected',[(0,True,'signal-unknown'),(1,False,'signal-unknown'),(1,True,'kernel-oom')])
def test_cgroup_counters_integrate_termination(tmp_path,delta,attributable,expected):
    def counters(pid,previous=None):
        return dict(cgroup='/fake',memory_events={'oom_kill':delta if previous else 0},memory_peak_bytes=123,memory_peak_scope='shared cgroup lifetime')
    code='import os,signal;os.kill(os.getpid(),signal.SIGKILL)'
    r=M.run_measured([sys.executable,'-c',code],timeout_s=2,artifact_dir=tmp_path,cgroup_reader=counters,oom_attributable=attributable)
    assert r['termination']['trigger']==expected
    assert r['termination']['memory_events_delta']=={'oom_kill':delta}
    assert r['termination']['memory_peak_scope']=='shared cgroup lifetime'
