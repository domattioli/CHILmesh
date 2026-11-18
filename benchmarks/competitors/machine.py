"""Machine and environment description recorded with every run (FR-009)."""
from __future__ import annotations

import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]


def _ram_bytes() -> Optional[int]:
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        return None


def _git_sha() -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, capture_output=True,
                             text=True, timeout=10, check=False)
        return out.stdout.strip() or "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def chilmesh_git_sha() -> str:
    return _git_sha()


def describe_machine(tool_versions: Optional[dict] = None) -> dict:
    return {
        "cpu": platform.processor() or platform.machine(),
        "cpu_count": os.cpu_count(),
        "ram_bytes": _ram_bytes(),
        "os": f"{platform.system()} {platform.release()}",
        "python": sys.version.split()[0],
        "git_sha": _git_sha(),
        "tool_versions": dict(tool_versions or {}),
    }


# Frozen campaign evidence. No privileged operation occurs in this module.
import hashlib
import time
from contextlib import contextmanager

DETECTOR_POLICY = {'cpu_percentage_points': 5, 'io_bytes_per_s': 1048576}
IDLE_SERVICE_ALLOWLIST = ()


def resolve_cgroup(cgroup_text, mountinfo):
    """Resolve cgroup v2 using the mount root (including delegated subtrees)."""
    path = next((s.split(':',2)[2] for s in cgroup_text.splitlines() if s.startswith('0::')), None)
    if path is None: raise ValueError('unavailable cgroup v2 membership')
    for line in mountinfo.splitlines():
        left, sep, right = line.partition(' - ')
        fields = left.split()
        if sep and right.split()[0] == 'cgroup2':
            root, mount = fields[3:5]
            if path == root or path.startswith(root.rstrip('/')+'/'):
                suffix = path[len(root):].lstrip('/')
                return path, str(Path(mount)/suffix)
    raise ValueError('unresolved cgroup mount')


def cgroup_snapshot(path):
    path = Path(path)
    result = dict(cgroup=str(path), memory_events=None, memory_peak_bytes=None, memory_peak_scope=None, missing={})
    try:
        result['memory_events'] = {k:int(v) for k,v in (line.split() for line in (path/'memory.events').read_text().splitlines())}
        if 'oom_kill' not in result['memory_events']: raise ValueError('oom_kill unavailable')
    except (OSError,ValueError) as exc:
        result['memory_events']=None;result['missing']['memory_events']=str(exc)
    try:
        result['memory_peak_bytes']=int((path/'memory.peak').read_text())
        result['memory_peak_scope']='cgroup lifetime; potentially shared, not child-only'
    except (OSError,ValueError) as exc:
        result['missing']['memory_peak_bytes']=str(exc)
        result['missing']['memory_peak_scope']='memory.peak unavailable'
    return result


def child_cgroup_snapshot(pid, *, previous=None, proc_root='/proc'):
    if previous and previous.get('cgroup'): return cgroup_snapshot(previous['cgroup'])
    try:
        root=Path(proc_root)/str(pid)
        _, mount=resolve_cgroup((root/'cgroup').read_text(),(root/'mountinfo').read_text())
        return cgroup_snapshot(mount)
    except (OSError,ValueError) as exc:
        return dict(cgroup=None,memory_events=None,memory_peak_bytes=None,memory_peak_scope=None,missing={'cgroup':str(exc)})


def host_identity(*, proc_root='/proc'):
    """Retain raw allocation facts; unavailable fields remain explicit nulls."""
    root=Path(proc_root);raw={};missing={}
    paths={'boot':root/'sys/kernel/random/boot_id','container':root/'self/cgroup','mountinfo':root/'self/mountinfo',
           'cpuinfo':root/'cpuinfo','meminfo':root/'meminfo'}
    for key,path in paths.items():
        try: raw[key]=path.read_text()
        except OSError as exc: raw[key]=None;missing[key]=str(exc)
    # Retain stable host facts only; free memory and CPU MHz change each run.
    if raw['cpuinfo'] is not None:
        raw['cpuinfo']='\n'.join(line for line in raw['cpuinfo'].splitlines() if line.startswith(('model name', 'Hardware', 'CPU architecture', 'CPU implementer', 'CPU part')))
    if raw['meminfo'] is not None:
        raw['meminfo']='\n'.join(line for line in raw['meminfo'].splitlines() if line.startswith('MemTotal:'))
    host_id=hashlib.sha256((raw['boot']+'\n'+raw['container']).encode()).hexdigest() if raw['boot'] and raw['container'] else None
    try: cgroup_path,mount=resolve_cgroup(raw['container'] or '', raw['mountinfo'] or '')
    except ValueError: cgroup_path=mount=None
    # Keep only the cgroup mount lines: other mounts (overlay snapshots, bind mounts) change while a container
    # runs and are not part of the CPU/memory allocation.
    if raw['mountinfo'] is not None:
        raw['mountinfo']='\n'.join(line for line in raw['mountinfo'].splitlines() if ' - cgroup' in line)
    for key in ('cpu.max','cpuset.cpus.effective','memory.max','memory.high'):
        try: raw[key]=(Path(mount)/key).read_text().strip() if mount else None
        except OSError: raw[key]=None
        if raw[key] is None: missing[key]='allocation counter unavailable'
    try: affinity=sorted(os.sched_getaffinity(0))
    except AttributeError: affinity=None
    quota=raw['cpu.max'];available=None
    if affinity is not None and quota:
        q,period=quota.split();available=min(len(affinity),float(q)/float(period)) if q!='max' else len(affinity)
    model=next((line.split(':',1)[1].strip() for line in (raw['cpuinfo'] or '').splitlines() if line.startswith('model name')),None)
    if model is None:
        # arm64 /proc/cpuinfo has no "model name"; identify the core by implementer and part.
        ids={k:next((line.split(':',1)[1].strip() for line in (raw['cpuinfo'] or '').splitlines() if line.startswith(k)),None)
             for k in ('CPU implementer','CPU part')}
        if all(ids.values()): model=f"{platform.machine()} implementer {ids['CPU implementer']} part {ids['CPU part']}"
    ram=next((int(line.split()[1])*1024 for line in (raw['meminfo'] or '').splitlines() if line.startswith('MemTotal:')),None)
    return dict(host_id=host_id,cpu_model=model,logical_count=os.cpu_count(),affinity=affinity,available_cores=available,
                physical_ram_bytes=ram,kernel=platform.release(),os=platform.system(),cgroup_path=cgroup_path,cgroup_mount=mount,
                cpu_quota={'raw':quota} if quota else None,cpuset=raw['cpuset.cpus.effective'],
                memory_limits={'max':raw['memory.max'],'high':raw['memory.high']} if raw['memory.max'] else None,
                raw_evidence=raw,missing=missing)


@contextmanager
def host_lock(path):
    """An advisory host lock; all campaign controllers must share this path."""
    import fcntl
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('a+') as stream:
        try: fcntl.flock(stream,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as exc: raise RuntimeError('overlapping benchmark campaign') from exc
        try: yield
        finally: fcntl.flock(stream,fcntl.LOCK_UN)


def load_snapshot(*, proc_root='/proc', host_probe=host_identity, clock=time.perf_counter):
    root=Path(proc_root);host=host_probe();missing=[]
    result=dict(monotonic_s=clock(),host_id=host['host_id'],allocation={k:host[k] for k in ('affinity','cpu_quota','cpuset','memory_limits')},
                processes={},complete=True,missing=missing)
    for key,file in (('proc_stat','stat'),('loadavg','loadavg'),('diskstats','diskstats')):
        try: result[key]=(root/file).read_text()
        except OSError as exc: result[key]=None;missing.append(str(exc))
    hz=os.sysconf('SC_CLK_TCK')
    try:
        values=list(map(int,result['proc_stat'].splitlines()[0].split()[1:9]))
        result['system_cpu_s']=sum(values[i] for i in (0,1,2,5,6,7))/hz
        # Keep raw diskstats. Per-process accounting below is the attribution
        # source; aggregate disk evidence is checked for unexplained activity.
        rows=[f for f in (line.split() for line in result['diskstats'].splitlines()) if len(f)>=14]
        names={f[2] for f in rows}
        # Whole devices only: a partition row repeats the sectors of its parent device.
        def partition(name):
            base=name.rstrip('0123456789')
            return base!=name and (base in names or (base.endswith('p') and base[:-1] in names))
        result['system_io_bytes']=sum((int(f[5])+int(f[9]))*512 for f in rows if not partition(f[2]))
    except (TypeError,ValueError,IndexError): missing.append('unreadable system accounting')
    try: entries=list(root.iterdir())
    except OSError: entries=[];missing.append('process inventory unavailable')
    for entry in entries:
        if not entry.name.isdigit(): continue
        try:
            stat=(entry/'stat').read_text();tail=stat[stat.rfind(')')+2:].split()
            io={k:int(v) for k,v in (line.split(':') for line in (entry/'io').read_text().splitlines())}
            result['processes'][entry.name]=dict(name=stat[stat.find('(')+1:stat.rfind(')')],starttime=int(tail[19]),
                                                cpu_s=(int(tail[11])+int(tail[12]))/hz,io_bytes=io['read_bytes']+io['write_bytes'])
        except FileNotFoundError: missing.append(f'process {entry.name} disappeared during inventory')
        except (OSError,ValueError,IndexError,KeyError) as exc: missing.append(f'process {entry.name}: {exc}')
    result['complete']=not missing and host['host_id'] is not None and not host.get('missing')
    return result


# Nominal interval is 0.100 s; a real sleep overshoots by a few ms. Operator ruling 2026-11-12.
LOAD_CADENCE_TOLERANCE_S=.150


def evaluate_load(observations, *, launch_s, host_id, allocation, allowlist=IDLE_SERVICE_ALLOWLIST, allowed_pids=()):
    observations=list(observations);unknown=[];contaminated=[]
    before=[o for o in observations if o['monotonic_s']<=launch_s]
    if not before or launch_s-before[0]['monotonic_s'] < 1-1e-9: unknown.append('missing one-second baseline')
    if not observations: unknown.append('missing observations')
    for o in observations:
        if not o.get('complete') or not all(o.get(k) is not None for k in ('proc_stat','loadavg','diskstats','processes','system_cpu_s','system_io_bytes')): unknown.append('unavailable load accounting')
        if o.get('host_id')!=host_id or o.get('allocation')!=allocation: contaminated.append('host/allocation drift')
    allowed={str(i) for i in allowed_pids};window=[]
    for a,b in zip(observations,observations[1:]):
        dt=b['monotonic_s']-a['monotonic_s']
        if dt<=0 or dt>LOAD_CADENCE_TOLERANCE_S+1e-8: unknown.append('load cadence gap');continue
        total_cpu=total_io=0
        for pid,p in b['processes'].items():
            old=a['processes'].get(pid)
            if not old or old['starttime']!=p['starttime']:
                if str(pid) not in allowed and p['name'] not in allowlist: contaminated.append('new foreign process')
                # An expected new process started inside this interval: all of its usage so far is attributed.
                elif not old: total_cpu+=p['cpu_s'];total_io+=p['io_bytes']
                continue
            cpu=p['cpu_s']-old['cpu_s'];io=p['io_bytes']-old['io_bytes']
            if cpu<0 or io<0: unknown.append('accounting counter reset');continue
            total_cpu+=cpu;total_io+=io
            if str(pid) in {str(i) for i in allowed_pids} or p['name'] in allowlist: continue
            if cpu/dt > .05+1e-9 or io/dt > 1048576+1e-6: contaminated.append(f'foreign PID {pid} CPU/I/O')
        # The measured child and the controller are expected to exit; only a foreign exit hides accounting.
        if any(pid not in allowed and a['processes'][pid]['name'] not in allowlist for pid in set(a['processes'])-set(b['processes'])): unknown.append('disappeared process accounting unavailable')
        # An expected process (the measured child) exited inside this interval: its final CPU and I/O are in the
        # system counters but no longer readable per process, so this one interval cannot be attributed either way.
        if set(a['processes'])-set(b['processes']): continue
        if all(k in a and k in b for k in ('system_cpu_s','system_io_bytes')):
            # CPU accounting ticks are 10 ms, so one 0.1 s interval resolves no better than 10 points.
            # Judge unattributed activity over a trailing window of at least one second.
            window.append((dt,b['system_cpu_s']-a['system_cpu_s']-total_cpu,b['system_io_bytes']-a['system_io_bytes']-total_io))
            while sum(w[0] for w in window[1:])>=1: window.pop(0)
            span=sum(w[0] for w in window)
            if span>=1-1e-9 and (sum(w[1] for w in window)/span > .05+1e-9 or sum(w[2] for w in window)/span > 1048576+1e-6):
                unknown.append('unattributable system activity')
    return dict(pre_run_s=launch_s-before[0]['monotonic_s'] if before else 0,interval_s=.100,observations=observations,
                allowlist=list(allowlist),detector_policy=dict(DETECTOR_POLICY),host_id=host_id,allocation=allocation,
                verdict='contaminated' if contaminated else ('unknown' if unknown else 'clean'),reasons=sorted(set(contaminated+unknown)))


class LoadObserver:
    def __init__(self, *, reader=load_snapshot, clock=time.perf_counter, sleep=time.sleep):
        self.reader,self.clock,self.sleep=reader,clock,sleep
        self.observations=[];self.next_s=None

    def poll(self):
        now=self.clock()
        if self.next_s is None or now>=self.next_s-1e-9:
            self.observations.append(self.reader());self.next_s=now+.100

    def baseline(self):
        start=self.clock();self.poll()
        while self.clock()-start < 1:
            self.sleep(min(.100,1-(self.clock()-start)));self.poll()
        return self.observations
