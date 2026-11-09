"""W2 tool table: which tools run, how they launch, and which secondary metrics they cannot produce."""
from __future__ import annotations

import os
from pathlib import Path

ADAPTER_DIR = Path(__file__).resolve().parents[1] / "adapters"

# Metrics every tool must answer (value, failure reason or justified na).
PRIMARY = "parse_s"
SECONDARY = ("adj_cpp_s", "adj_python_s", "tree_build_s", "query_per_s", "render_s",
             "station_max_abs_diff")
METRICS = (PRIMARY,) + SECONDARY

# Human-readable metric descriptions. Metric names stay fixed.
METRIC_DESCRIPTIONS = {
    "adj_cpp_s": "mesh construction with adjacencies and layers, C++ backend (CHILMESH_BACKEND=cpp)",
    "adj_python_s": "mesh construction with adjacencies and layers, Python backend (CHILMESH_BACKEND=python)",
}

LICENSE_UNVERIFIED = "unverified (research.md T002)"

# Competitor tools run as subprocess adapters in their own environments.
COMPETITORS = {
    "adcircpy": {"env": "CHILMESH_BENCH_PY_ADCIRCPY", "script": ADAPTER_DIR / "adcircpy_read.py",
                 "input": "fort.14", "license": LICENSE_UNVERIFIED,
                 "na": "no equivalent public API for this operation (parse only)"},
    "thalassa": {"env": "CHILMESH_BENCH_PY_THALASSA", "script": ADAPTER_DIR / "thalassa_read.py",
                 "input": "fort.14", "license": LICENSE_UNVERIFIED,
                 "na": "no equivalent public API for this operation (parse only)"},
    "xarray": {"env": "CHILMESH_BENCH_PY_XARRAY", "script": ADAPTER_DIR / "xarray_read.py",
               "input": "netcdf", "license": LICENSE_UNVERIFIED,
               "na": "xarray is a generic array library with no mesh-operation API"},
}

XARRAY_FORT14_REASON = "xarray reads NetCDF; the mesh is served as fort.14, so it is not run"


def competitor_python(tool: str):
    """Interpreter for a competitor environment, or None when its env var is unset."""
    return os.environ.get(COMPETITORS[tool]["env"])


def all_tools() -> list:
    return ["chilmesh"] + sorted(COMPETITORS)
