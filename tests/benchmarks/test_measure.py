"""Subprocess runner: wall time, peak RSS, timeout and OOM as failed with cause."""
from __future__ import annotations

import sys

import pytest

from benchmarks.competitors import machine
from benchmarks.competitors.measure import run_measured

pytestmark = pytest.mark.skipif(sys.platform == "win32", reason="harness targets POSIX")


def test_success_records_wall_time_and_rss():
    code = "import time; b = bytearray(40*1024*1024); b[::4096] = b'x'*len(b[::4096]); time.sleep(0.2); print('hi')"
    res = run_measured([sys.executable, "-c", code], timeout_s=30)
    assert res["status"] == "ok"
    assert res["returncode"] == 0
    assert res["wall_s"] >= 0.2
    assert res["peak_rss_bytes"] >= 40 * 1024 * 1024
    assert "hi" in res["stdout_tail"]


def test_timeout_is_failed_with_cause_and_peak_reached():
    res = run_measured([sys.executable, "-c", "import time; time.sleep(60)"], timeout_s=0.5)
    assert res["status"] == "failed"
    assert "timeout" in res["reason"]
    assert res["peak_rss_bytes"] > 0


def test_sigkill_is_failed_cause_unknown():
    code = "import os, signal; os.kill(os.getpid(), signal.SIGKILL)"
    res = run_measured([sys.executable, "-c", code], timeout_s=30)
    assert res["status"] == "failed"
    assert "SIGKILL" in res["reason"] and "unknown" in res["reason"]
    assert res["termination"]["trigger"] == "signal-unknown"


def test_nonzero_exit_is_failed_with_stderr_tail():
    res = run_measured([sys.executable, "-c", "import sys; sys.stderr.write('boom'); sys.exit(3)"], timeout_s=30)
    assert res["status"] == "failed"
    assert res["returncode"] == 3
    assert "boom" in res["stderr_tail"]
    assert "exit code 3" in res["reason"]


def test_missing_executable_is_failed_not_exception():
    res = run_measured(["/nonexistent/binary-xyz"], timeout_s=5)
    assert res["status"] == "failed"
    assert res["reason"]


def test_tail_is_bounded():
    res = run_measured([sys.executable, "-c", "print('x'*100000)"], timeout_s=30, tail_bytes=500)
    assert len(res["stdout_tail"]) <= 500


def test_machine_description_fields():
    desc = machine.describe_machine(tool_versions={"numpy": "x"})
    for key in ("cpu", "ram_bytes", "os", "python", "git_sha", "tool_versions"):
        assert key in desc
    assert desc["tool_versions"] == {"numpy": "x"}


def _samples(values, start=1.):
    return [dict(monotonic_s=start+i*.01, pid=42, rss_bytes=v, descendants=[]) for i,v in enumerate(values)]


def test_region_windows_exclude_imports_and_verification(tmp_path):
    from benchmarks.competitors.measure import summarize_samples, save_samples
    samples = _samples([900, 100, 120, 160, 110, 9999])
    m = summarize_samples(samples, launch_s=1., start_s=1.015, end_s=1.045)
    assert m['rss_baseline_bytes'] == 100
    assert m['rss_operation_peak_bytes'] == 160
    assert m['rss_whole_process_peak_bytes'] == 900
    assert m['rss_growth_bytes'] == 60 and m['coverage']
    assert m['operation_elapsed_s'] == pytest.approx(.03)
    path, digest = save_samples(samples, tmp_path/'samples.json')
    import hashlib
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    with pytest.raises(FileExistsError): save_samples(samples, path)


@pytest.mark.parametrize('case', ['empty', 'baseline', 'region', 'gap', 'missing', 'descendant', 'markers'])
def test_sampler_missing_never_uses_wait4(case):
    from benchmarks.competitors.measure import summarize_samples
    samples = _samples([100, 90, 80, 70])
    start, end = 1.005, 1.025
    if case == 'empty': samples = []
    if case == 'baseline': start = 1.
    if case == 'region': start, end = 1.001, 1.002
    if case == 'gap': samples[2]['monotonic_s'] = 1.04; samples[3]['monotonic_s'] = 1.05
    if case == 'missing': samples[1].update(rss_bytes=None, missing_reason='gone')
    if case == 'descendant': samples[1]['descendants'] = [43]
    if case == 'markers': end = None
    m = summarize_samples(samples, launch_s=1., start_s=start, end_s=end)
    assert not m['coverage'] and m['coverage_reasons']
    assert m['wait4_lifetime_maxrss_bytes'] is None


def test_growth_is_clipped_and_procfs_uses_bytes(tmp_path):
    from benchmarks.competitors.measure import summarize_samples, read_rss_sample
    m = summarize_samples(_samples([100, 90, 80]), launch_s=1., start_s=1.005, end_s=1.025)
    assert m['rss_growth_bytes'] == 0
    p = tmp_path/'42'; (p/'task'/'42').mkdir(parents=True)
    (p/'status').write_text('VmRSS: 123 kB\n')
    (p/'task'/'42'/'children').write_text('43 44')
    s = read_rss_sample(42, proc_root=tmp_path, clock=lambda: 1.)
    assert s['rss_bytes'] == 123*1024 and s['descendants'] == [43,44]


def test_sigkill_requires_attributable_counter_delta():
    from benchmarks.competitors.measure import termination_trigger
    assert termination_trigger(-9, before={'oom_kill':1}, after={'oom_kill':1}, attributable=True) == 'signal-unknown'
    assert termination_trigger(-9, before={'oom_kill':1}, after={'oom_kill':2}, attributable=False) == 'signal-unknown'
    assert termination_trigger(-9, before={'oom_kill':1}, after={'oom_kill':2}, attributable=True) == 'kernel-oom'
    assert termination_trigger(-9, before={'oom_kill':1}, after={'oom_kill':2}, attributable=True, limit='rss-limit') == 'rss-limit'
