"""One subprocess lifecycle and Linux 10 ms RSS observer for every tool.

Legacy wall/peak keys retain their v1 meaning. Fair sampled peaks stop at END;
wait4 lifetime RSS remains a separately labelled diagnostic, never a fallback.
"""
from __future__ import annotations
import hashlib
import json
import os
import selectors
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from .cache import utc_now
from .protocol import ProtocolState, ProtocolError, clock_metadata

_POLL_S = .010


def _tail(path, n):
    with open(path, 'rb') as fh:
        fh.seek(0, os.SEEK_END)
        fh.seek(max(0, fh.tell()-n))
        return fh.read().decode('utf-8', errors='replace')[-n:]


def _rss_bytes(value):
    return int(value) if sys.platform == 'darwin' else int(value)*1024


def read_rss_sample(pid, *, proc_root='/proc', clock=time.perf_counter):
    root = Path(proc_root)/str(pid)
    sample = dict(monotonic_s=clock(), pid=pid, rss_bytes=None, descendants=[], missing_reason='')
    try:
        status = root.joinpath('status').read_text()
        line = next(line for line in status.splitlines() if line.startswith('VmRSS:'))
        fields = line.split()
        if fields[2] != 'kB': raise ValueError('unexpected VmRSS unit')
        sample['rss_bytes'] = int(fields[1])*1024
        descendants = set()
        for task in root.joinpath('task').iterdir():
            descendants.update(int(p) for p in task.joinpath('children').read_text().split())
        sample['descendants'] = sorted(descendants)
    except (OSError, ValueError, StopIteration) as exc:
        sample['missing_reason'] = f'procfs unavailable: {exc}'
    return sample


def summarize_samples(samples, *, launch_s, start_s, end_s):
    samples = list(samples)
    reasons = []
    markers = start_s is not None and end_s is not None and launch_s <= start_s <= end_s
    before = [s for s in samples if launch_s <= s['monotonic_s'] < start_s] if markers else []
    region = [s for s in samples if start_s <= s['monotonic_s'] <= end_s] if markers else []
    whole = [s for s in samples if launch_s <= s['monotonic_s'] <= end_s] if markers else []
    gaps = [b['monotonic_s']-a['monotonic_s'] for a,b in zip(samples,samples[1:])]
    if not markers: reasons.append('missing or unordered region markers')
    if not before: reasons.append('missing baseline sample')
    if not region: reasons.append('missing in-region sample')
    if any(s['rss_bytes'] is None or s['descendants'] or s.get('missing_reason') for s in samples): reasons.append('missing RSS or descendants in child-only observation')
    if any(g > .020 + 1e-12 or g < 0 for g in gaps): reasons.append('sample gap exceeds 20ms or out of order')
    if whole and (whole[0]['monotonic_s']-launch_s > .020 or end_s-whole[-1]['monotonic_s'] > .020): reasons.append('uncovered boundary')
    baseline = before[-1]['rss_bytes'] if before else None
    def peak(rows):
        values = [s['rss_bytes'] for s in rows if s['rss_bytes'] is not None]
        return max(values) if values else None
    operation = peak(region)
    result = dict(operation_elapsed_s=end_s-start_s if markers else None, process_wall_s=None,
                  rss_baseline_bytes=baseline, rss_operation_peak_bytes=operation,
                  rss_whole_process_peak_bytes=peak(whole),
                  rss_growth_bytes=max(0,operation-baseline) if baseline is not None and operation is not None else None,
                  last_rss_bytes=samples[-1]['rss_bytes'] if samples else None, wait4_lifetime_maxrss_bytes=None,
                  policy_version='procfs-v1', interval_s=.010, scope='child-only', coverage=not reasons,
                  coverage_reasons=reasons, gaps=gaps, sample_count=len(samples), samples=samples,
                  samples_path=None, samples_sha256=None)
    result['missing'] = {k: 'unavailable observation; no lifetime RSS substitution' for k,v in result.items() if v is None}
    return result


def save_samples(samples, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = (json.dumps(samples, sort_keys=True, allow_nan=False)+'\n').encode()
    with path.open('xb') as stream: stream.write(raw)
    return path, hashlib.sha256(raw).hexdigest()


def limit_trigger(elapsed, rss, *, wall_s, rss_bytes):
    if elapsed > wall_s: return 'timeout'
    if rss_bytes is not None and rss is not None and rss > rss_bytes: return 'rss-limit'
    return None


def termination_trigger(returncode, *, before=None, after=None, attributable=False, limit=None):
    if limit: return limit
    if returncode == -signal.SIGKILL and attributable and before is not None and after is not None and after.get('oom_kill',0) > before.get('oom_kill',0): return 'kernel-oom'
    if returncode is not None and returncode < 0: return 'signal-unknown'
    if returncode: return 'nonzero-exit'
    return 'none'


def run_measured(cmd, *, timeout_s, env=None, cwd=None, tail_bytes=4000,
                 rss_limit_bytes=None, repeat_id=None, artifact_dir=None,
                 sample_reader=read_rss_sample, cgroup_reader=None, oom_attributable=False,
                 load_observer=None):
    """Launch a fresh group; preserve partial evidence and reap on every exit.

    repeat_id enables the fair control protocol. result.json is an untruncated
    inherited file channel, independent of stdout. No cache operation occurs here.
    """
    from .machine import child_cgroup_snapshot
    cgroup_reader = cgroup_reader or child_cgroup_snapshot
    directory = Path(artifact_dir) if artifact_dir is not None else Path(tempfile.mkdtemp(prefix='chilmesh-measure-'))
    directory.mkdir(parents=True, exist_ok=True)
    out_path, err_path, result_path = (directory/n for n in ('stdout.txt','stderr.txt','result.json'))
    result = dict(command=list(cmd), argv=list(cmd), status='failed', reason='', returncode=None,
                  wall_s=None, peak_rss_bytes=None, stdout_tail='', stderr_tail='')
    protocol = dict(pid=None, process_group_id=None, ready_s=None, start_s=None, end_s=None, ack_s=None,
                    clock=None, receipt_timestamps={}, new_modules=None, materialization=None)
    termination = dict(launch_utc=None, exit_utc=None, launch_s=None, exit_s=None, returncode=None, signal=None,
                       limits={'wall_s':timeout_s,'rss_bytes':rss_limit_bytes}, trigger='launch-error', cgroup=None,
                       memory_events_before=None, memory_events_after=None, memory_events_delta=None,
                       memory_peak_bytes=None, memory_peak_scope=None, kill_diagnostics=[], reap_diagnostics={},
                       stdout_path=str(out_path), stderr_path=str(err_path))
    samples=[]; proc=None; fds=[]; state=ProtocolState(repeat_id) if repeat_id else None
    sel=selectors.DefaultSelector(); buffer=b''; limit=None; protocol_error=''; usage=None; before={}; after={}
    launch=time.perf_counter(); termination.update(launch_s=launch,launch_utc=utc_now())
    try:
        with out_path.open('xb') as out, err_path.open('xb') as err, result_path.open('x+b') as structured:
            child_env=dict(os.environ if env is None else env)
            inherited=[]
            if state:
                child_in, parent_out=os.pipe(); parent_in, child_out=os.pipe()
                fds.extend([child_in,parent_out,parent_in,child_out])
                os.set_blocking(parent_in,False);sel.register(parent_in,selectors.EVENT_READ)
                child_env.update(CHILMESH_CONTROL_IN_FD=str(child_in),CHILMESH_CONTROL_OUT_FD=str(child_out),
                                 CHILMESH_RESULT_FD=str(structured.fileno()),CHILMESH_REPEAT_ID=repeat_id)
                inherited=[child_in,child_out,structured.fileno()]
            try:
                proc=subprocess.Popen(list(cmd),stdout=out,stderr=err,env=child_env,cwd=cwd,
                                      start_new_session=True,pass_fds=inherited)
            except OSError as exc:
                result['reason']=f'could not start command: {exc}'
                termination['launch_s']=None;termination['launch_utc']=None
            if proc is not None:
                protocol.update(pid=proc.pid,process_group_id=proc.pid)
                if state:
                    for fd in (child_in,child_out): os.close(fd);fds.remove(fd)
                before=cgroup_reader(proc.pid)
                termination.update(cgroup=before.get('cgroup'),memory_events_before=before.get('memory_events'))
                next_sample=launch
                def sample():
                    s=sample_reader(proc.pid);samples.append(s);return s
                def send(event):
                    msg=dict(protocol_version=1,repeat_id=repeat_id,event=event,monotonic_s=time.perf_counter(),
                             clock='perf_counter',clock_metadata=state.events[0]['clock_metadata'])
                    state.accept(msg,received_s=msg['monotonic_s'])
                    os.write(parent_out,(json.dumps(msg)+'\n').encode())
                    protocol['receipt_timestamps'][event]=msg['monotonic_s']
                    if event=='ACK': protocol['ack_s']=msg['monotonic_s']
                while True:
                    pid,status,usage_now=os.wait4(proc.pid,os.WNOHANG)
                    if pid:
                        usage=usage_now;proc.returncode=os.waitstatus_to_exitcode(status)
                        termination['reap_diagnostics']={'reaped':True};break
                    now=time.perf_counter()
                    if now>=next_sample:
                        s=sample();next_sample=now+_POLL_S
                        if s['missing_reason']:
                            # Zombie between exit and reap has no VmRSS; drop that post-exit sample.
                            pid,status,usage_now=os.wait4(proc.pid,os.WNOHANG)
                            if pid:
                                samples.pop();usage=usage_now;proc.returncode=os.waitstatus_to_exitcode(status)
                                termination['reap_diagnostics']={'reaped':True};break
                        limit=limit_trigger(now-launch,s['rss_bytes'],wall_s=timeout_s,rss_bytes=rss_limit_bytes)
                    else: limit=limit_trigger(now-launch,None,wall_s=timeout_s,rss_bytes=rss_limit_bytes)
                    if load_observer is not None: load_observer.poll()
                    if limit: break
                    if state:
                        for key,_ in sel.select(0):
                            chunk=os.read(key.fd,65536)
                            if not chunk: sel.unregister(key.fd);continue
                            buffer+=chunk
                            while b'\n' in buffer:
                                line,buffer=buffer.split(b'\n',1)
                                try:
                                    msg=json.loads(line);state.accept(msg,received_s=time.perf_counter())
                                    event=msg['event'];protocol['receipt_timestamps'][event]=state.receipts[-1]
                                    protocol[event.lower()+'_s']=msg['monotonic_s'];protocol['clock']=msg['clock_metadata']
                                    if event=='READY': sample();send('GO')
                                    if event=='END': send('ACK')
                                except (ValueError,TypeError,KeyError) as exc:
                                    protocol_error=f'protocol failure: {exc}';break
                        if protocol_error: break
                    time.sleep(min(.002,max(0,next_sample-time.perf_counter())))
                after=cgroup_reader(proc.pid, previous=before)
                if limit or protocol_error:
                    _kill_group(proc,termination)
                    _,status,usage=os.wait4(proc.pid,0);proc.returncode=os.waitstatus_to_exitcode(status)
                    termination['reap_diagnostics']={'reaped':True}
                result.update(returncode=proc.returncode,wall_s=time.perf_counter()-launch,
                              peak_rss_bytes=_rss_bytes(usage.ru_maxrss) if usage else None)
                termination.update(exit_s=launch+result['wall_s'],exit_utc=utc_now(),returncode=proc.returncode,
                                   signal=-proc.returncode if proc.returncode<0 else None,
                                   memory_events_after=after.get('memory_events'),memory_peak_bytes=after.get('memory_peak_bytes'),
                                   memory_peak_scope=after.get('memory_peak_scope'))
                b,a=before.get('memory_events'),after.get('memory_events')
                if b is not None and a is not None and set(b)==set(a) and all(a[k]>=b[k] for k in b):
                    termination['memory_events_delta']={k:a[k]-b[k] for k in b}
                termination['trigger']=termination_trigger(proc.returncode,before=b,after=a,attributable=oom_attributable,limit=limit)
                if state:
                    try:
                        state.finish()
                        structured.seek(0);result['adapter_result']=json.load(structured)
                        protocol.update(new_modules=result['adapter_result'].get('new_modules'),materialization=result['adapter_result'].get('materialization'))
                    except (ValueError,TypeError) as exc: protocol_error=protocol_error or f'protocol/result incomplete: {exc}'
                trigger=termination['trigger']
                if limit: result['reason']=f'{limit} after {result["wall_s"]:g} s'
                elif protocol_error: result['reason']=protocol_error
                elif proc.returncode==-9: result['reason']='killed by SIGKILL; '+('kernel OOM evidenced' if trigger=='kernel-oom' else 'cause unknown')
                elif proc.returncode<0: result['reason']=f'killed by signal {-proc.returncode}; cause unknown'
                elif proc.returncode: result['reason']=f'exit code {proc.returncode}'
                else: result['status']='ok'
    finally:
        if proc is not None and proc.returncode is None:
            _kill_group(proc,termination)
            try:
                _,status,usage=os.wait4(proc.pid,0);proc.returncode=os.waitstatus_to_exitcode(status)
                termination['reap_diagnostics']={'reaped':True}
            except ChildProcessError: pass
        # Also remove descendants left behind by a normally exited direct child.
        if proc is not None:
            try: os.killpg(proc.pid,signal.SIGKILL)
            except ProcessLookupError: pass
        sel.close()
        for fd in fds:
            os.close(fd)
        # A process in exit teardown has no VmRSS before it can be reaped. Such samples after the
        # END marker observe no live address space; drop them instead of failing coverage.
        while (samples and proc is not None and proc.returncode is not None and protocol['end_s'] is not None
               and samples[-1]['rss_bytes'] is None and samples[-1]['missing_reason'] and samples[-1]['monotonic_s']>protocol['end_s']):
            samples.pop()
        m=summarize_samples(samples,launch_s=launch,start_s=protocol['start_s'],end_s=protocol['end_s'])
        m.update(process_wall_s=result['wall_s'],wait4_lifetime_maxrss_bytes=result['peak_rss_bytes'])
        path,digest=save_samples(samples,directory/'samples.json')
        m.update(samples_path=str(path),samples_sha256=digest)
        m['missing']={k:'unavailable observation' for k,v in m.items() if v is None}
        for name,path in (('stdout',out_path),('stderr',err_path)):
            if path.exists():
                result[name+'_tail']=_tail(path,tail_bytes)
                termination[name+'_tail']=result[name+'_tail']
                with path.open('rb') as stream: termination[name+'_sha256']=hashlib.file_digest(stream,'sha256').hexdigest()
        for section in (protocol,termination): section['missing']={k:'unavailable evidence' for k,v in section.items() if v is None}
        result.update(measurement=m,protocol=protocol,termination=termination)
        if sys.exc_info()[0] is not None:
            result.update(status='failed',reason='controller abort',returncode=proc.returncode if proc else None,
                          wall_s=time.perf_counter()-launch)
            termination.update(trigger='signal-unknown' if proc and proc.returncode and proc.returncode<0 else 'nonzero-exit',
                               exit_s=launch+result['wall_s'],exit_utc=utc_now(),returncode=result['returncode'],
                               signal=-proc.returncode if proc and proc.returncode and proc.returncode<0 else None)
            m['process_wall_s']=result['wall_s'];m['missing'].pop('process_wall_s',None)
            termination['missing']={k:'unavailable evidence' for k,v in termination.items() if v is None}
        with (directory/'controller.json').open('x') as stream:
            json.dump(result,stream,sort_keys=True,allow_nan=False)
    return result


def _kill_group(proc, termination):
    try:
        os.killpg(proc.pid,signal.SIGKILL)
        termination['kill_diagnostics'].append({'process_group_id':proc.pid,'signal':signal.SIGKILL})
    except ProcessLookupError: pass
