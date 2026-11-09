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
