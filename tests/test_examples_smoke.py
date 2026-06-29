"""Smoke-test example scripts to prevent silent rot.

Executes each example script in an isolated subprocess to prevent in-process
mutations of the cached test mesh from breaking downstream tests.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples"


@pytest.mark.parametrize("script", [
    "01_quickstart.py",
    "02_fort14_roundtrip.py",
    "03_smoothing.py",
    "04_spatial_queries.py",
])
def test_example_script_runs(script: str, tmp_path: Path) -> None:
    """Execute example script in isolated subprocess without error."""
    path = EXAMPLES_DIR / script
    if not path.exists():
        pytest.skip(f"{script} not present")

    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"

    result = subprocess.run(
        [sys.executable, str(path)],
        cwd=str(tmp_path),
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
    )

    assert result.returncode == 0, (
        f"{script} failed (exit {result.returncode}):\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )
