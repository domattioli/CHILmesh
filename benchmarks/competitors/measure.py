"""Run one command in a fresh subprocess; record wall time, peak RSS and failure cause.

Peak RSS comes from ``os.wait4`` for that exact child (POSIX), so it is not
polluted by earlier children. Linux reports ``ru_maxrss`` in KiB, macOS in bytes.
"""
from __future__ import annotations

import os
import signal
import subprocess
import sys
import tempfile
import time
from typing import Optional, Sequence

_POLL_S = 0.02


def _tail(path: str, n: int) -> str:
    with open(path, "rb") as fh:
        fh.seek(0, os.SEEK_END)
        size = fh.tell()
        fh.seek(max(0, size - n))
        return fh.read().decode("utf-8", errors="replace")[-n:]


def _rss_bytes(ru_maxrss: int) -> int:
    return int(ru_maxrss) if sys.platform == "darwin" else int(ru_maxrss) * 1024


def run_measured(cmd: Sequence[str], *, timeout_s: float, env: Optional[dict] = None,
                 cwd: Optional[str] = None, tail_bytes: int = 4000) -> dict:
    """Run ``cmd``; never raise. Returns a dict ready to merge into a run record."""
    out_f = tempfile.NamedTemporaryFile(delete=False)
    err_f = tempfile.NamedTemporaryFile(delete=False)
    result = {"command": list(cmd), "status": "failed", "reason": "", "returncode": None,
              "wall_s": None, "peak_rss_bytes": None, "stdout_tail": "", "stderr_tail": ""}
    try:
        t0 = time.perf_counter()
        try:
            proc = subprocess.Popen(list(cmd), stdout=out_f, stderr=err_f, env=env, cwd=cwd)
        except OSError as exc:
            result["reason"] = f"could not start command: {exc}"
            return result
        timed_out = False
        while True:
            pid, status, rusage = os.wait4(proc.pid, os.WNOHANG)
            if pid:
                break
            if time.perf_counter() - t0 > timeout_s:
                timed_out = True
                proc.kill()
                pid, status, rusage = os.wait4(proc.pid, 0)
                break
            time.sleep(_POLL_S)
        result["wall_s"] = time.perf_counter() - t0
        result["peak_rss_bytes"] = _rss_bytes(rusage.ru_maxrss)
        if os.WIFSIGNALED(status):
            sig = os.WTERMSIG(status)
            result["returncode"] = -sig
        else:
            result["returncode"] = os.WEXITSTATUS(status)
        if timed_out:
            result["reason"] = f"timeout after {timeout_s:g} s (peak RSS {result['peak_rss_bytes']} bytes)"
        elif result["returncode"] == -signal.SIGKILL:
            result["reason"] = (f"killed by SIGKILL, OOM suspected (peak RSS {result['peak_rss_bytes']} bytes)")
        elif result["returncode"] < 0:
            result["reason"] = f"killed by signal {-result['returncode']}"
        elif result["returncode"] != 0:
            result["reason"] = f"exit code {result['returncode']}"
        else:
            result["status"] = "ok"
    finally:
        out_f.close()
        err_f.close()
        result["stdout_tail"] = _tail(out_f.name, tail_bytes)
        result["stderr_tail"] = _tail(err_f.name, tail_bytes)
        for name in (out_f.name, err_f.name):
            try:
                os.unlink(name)
            except OSError:
                pass
    return result
