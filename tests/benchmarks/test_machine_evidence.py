import copy
import pytest


def observation(t, cpu=0, io=0, **changes):
    d=dict(monotonic_s=t,host_id='host',allocation={'memory':'8G'},proc_stat='cpu 1',loadavg='0',diskstats='0',
           processes={'99':dict(starttime=1,name='foreign',cpu_s=cpu,io_bytes=io)}, complete=True,
           system_cpu_s=cpu,system_io_bytes=io)
    d.update(changes); return d


def test_frozen_load_thresholds_baseline_and_unknown():
    from benchmarks.competitors.machine import evaluate_load
    clean=[observation(i/10) for i in range(12)]
    assert evaluate_load(clean,launch_s=1.,host_id='host',allocation={'memory':'8G'})['verdict']=='clean'
    for changes in ({'host_id':'other'}, {'allocation':{}}, {'complete':False}):
        data=copy.deepcopy(clean); data[-1].update(changes)
        assert evaluate_load(data,launch_s=1.,host_id='host',allocation={'memory':'8G'})['verdict']!='clean'
    assert evaluate_load(clean[2:],launch_s=1.,host_id='host',allocation={'memory':'8G'})['verdict']=='unknown'
    for cpu,io in [(.00501,0),(0,104858)]:
        data=copy.deepcopy(clean); data[-1]=observation(1.1,cpu,io)
        assert evaluate_load(data,launch_s=1.,host_id='host',allocation={'memory':'8G'})['verdict']=='contaminated'
    data=copy.deepcopy(clean); data[-1]=observation(1.1,.005,104857.6)
    assert evaluate_load(data,launch_s=1.,host_id='host',allocation={'memory':'8G'})['verdict']=='clean'
    assert evaluate_load(data,launch_s=1.,host_id='host',allocation={'memory':'8G'},allowlist=['foreign'])['verdict']=='clean'


def test_cgroup_mount_resolution_counters_scope_and_missing(tmp_path):
    from benchmarks.competitors.machine import resolve_cgroup, cgroup_snapshot
    cgroup='0::/parent/job\n'
    mounts='29 23 0:26 /parent /sys/fs/cgroup rw - cgroup2 cgroup rw\n'
    assert resolve_cgroup(cgroup,mounts)==('/parent/job','/sys/fs/cgroup/job')
    p=tmp_path/'cg';p.mkdir();(p/'memory.events').write_text('oom 1\noom_kill 2\n');(p/'memory.peak').write_text('123')
    r=cgroup_snapshot(p)
    assert r['memory_events']=={'oom':1,'oom_kill':2} and r['memory_peak_bytes']==123
    assert 'cgroup' in r['memory_peak_scope']
    (p/'memory.events').unlink()
    assert cgroup_snapshot(p)['memory_events'] is None
    assert cgroup_snapshot(p)['missing']


def test_overlap_lock_is_exclusive_and_released(tmp_path):
    from benchmarks.competitors.machine import host_lock
    with host_lock(tmp_path/'host.lock'):
        with pytest.raises(RuntimeError):
            with host_lock(tmp_path/'host.lock'): pass
    with host_lock(tmp_path/'host.lock'): pass


def test_load_observer_one_second_baseline_100ms():
    from benchmarks.competitors.machine import LoadObserver
    now=[0.]
    def sleep(dt): now[0]+=dt
    observer=LoadObserver(reader=lambda:observation(now[0]),clock=lambda:now[0],sleep=sleep)
    observer.baseline()
    assert len(observer.observations)>=11 and now[0]>=1.
    assert max(b['monotonic_s']-a['monotonic_s'] for a,b in zip(observer.observations,observer.observations[1:]))<=.100000001


def test_host_identity_stable_despite_free_memory_or_cpu_frequency(tmp_path,monkeypatch):
    from benchmarks.competitors.machine import host_identity
    root=tmp_path/'proc';cg=tmp_path/'cg';cg.mkdir()
    for name,value in [('cpu.max','200000 100000'),('cpuset.cpus.effective','0-3'),('memory.max','8589934592'),('memory.high','max')]: (cg/name).write_text(value)
    (root/'self').mkdir(parents=True);(root/'sys/kernel/random').mkdir(parents=True)
    (root/'sys/kernel/random/boot_id').write_text('fake-boot')
    (root/'self/cgroup').write_text('0::/job\n')
    (root/'self/mountinfo').write_text(f'29 23 0:26 /job {cg} rw - cgroup2 cgroup rw\n')
    (root/'cpuinfo').write_text('model name: synthetic\ncpu MHz: 1000\n')
    (root/'meminfo').write_text('MemTotal: 100000 kB\nMemFree: 100 kB\n')
    monkeypatch.setattr('benchmarks.competitors.machine.os.sched_getaffinity',lambda _:set(range(4)),raising=False)
    first=host_identity(proc_root=root)
    assert first['available_cores']==2 and first['physical_ram_bytes']==102400000 and not first['missing']
    (root/'meminfo').write_text('MemTotal: 100000 kB\nMemFree: 200 kB\n')
    (root/'cpuinfo').write_text('model name: synthetic\ncpu MHz: 3000\n')
    assert host_identity(proc_root=root)==first
    (cg/'memory.max').write_text('4294967296')
    assert host_identity(proc_root=root)['memory_limits']!=first['memory_limits']
