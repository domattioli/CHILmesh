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


def test_sigkill_is_failed_oom_suspected():
    code = "import os, signal; os.kill(os.getpid(), signal.SIGKILL)"
    res = run_measured([sys.executable, "-c", code], timeout_s=30)
    assert res["status"] == "failed"
    assert "SIGKILL" in res["reason"] and "OOM" in res["reason"]


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
