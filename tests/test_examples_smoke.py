"""Smoke-test example scripts to prevent silent rot.

Executes each example script end-to-end via runpy.
"""
from __future__ import annotations

import runpy
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples"


@pytest.mark.parametrize("script", [
    "01_quickstart.py",
    "02_fort14_roundtrip.py",
    "03_smoothing.py",
    "04_spatial_queries.py",
])
def test_example_script_runs(script: str) -> None:
    """Execute example script without error."""
    path = EXAMPLES_DIR / script
    if not path.exists():
        pytest.skip(f"{script} not present")
    runpy.run_path(str(path), run_name="__main__")
