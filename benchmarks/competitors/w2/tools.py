"""W2 tool table: which tools run, how they launch, and which secondary metrics they cannot produce."""
from __future__ import annotations

import os
from pathlib import Path

ADAPTER_DIR = Path(__file__).resolve().parents[1] / "adapters"

# Metrics every tool must answer (value, failure reason or justified na).
PRIMARY = "parse_s"
SECONDARY = ("adj_cpp_s", "adj_python_s", "tree_build_s", "query_per_s", "render_s",
             "station_max_abs_diff", "write_fort14_s", "roundtrip_ok", "full_init_s",
             "import_s", "install_mb", "peak_rss_per_melem", "interp_max_abs_diff")
METRICS = (PRIMARY,) + SECONDARY

# Human-readable metric descriptions. Metric names stay fixed.
METRIC_DESCRIPTIONS = {
    "adj_cpp_s": "mesh construction with adjacencies and layers, C++ backend (CHILMESH_BACKEND=cpp)",
    "adj_python_s": "mesh construction with adjacencies and layers, Python backend (CHILMESH_BACKEND=python)",
    "write_fort14_s": "time to write a fort.14 file",
    "roundtrip_ok": "whether nodes, elements and boundaries survive a fort.14 round trip",
    "full_init_s": "full CHILmesh construction with layers and adjacencies",
    "import_s": "fresh subprocess wall time for import chilmesh",
    "install_mb": "installed distribution file size in MiB",
    "peak_rss_per_melem": "peak RSS in bytes per million elements for each subprocess step",
    "interp_max_abs_diff": "station interpolation maximum absolute difference",
}

LICENSE_UNVERIFIED = "unverified (research.md T002)"

# Competitor tools run as subprocess adapters in their own environments.
COMPETITORS = {
    "adcircpy": {"env": "CHILMESH_BENCH_PY_ADCIRCPY", "script": ADAPTER_DIR / "adcircpy_read.py",
                 "input": "fort.14", "license": LICENSE_UNVERIFIED,
                 "na": "no equivalent public API for this operation (parse only)"},
    "thalassa": {"env": "CHILMESH_BENCH_PY_THALASSA", "script": ADAPTER_DIR / "thalassa_read.py",
                 "input": "netcdf", "license": LICENSE_UNVERIFIED,
                 "na": "no equivalent public API for this operation (parse only)"},
    "xarray": {"env": "CHILMESH_BENCH_PY_XARRAY", "script": ADAPTER_DIR / "xarray_read.py",
               "input": "netcdf", "license": LICENSE_UNVERIFIED,
               "na": "xarray is a generic array library with no mesh-operation API"},
}


def competitor_python(tool: str):
    """Interpreter for a competitor environment, or None when its env var is unset."""
    return os.environ.get(COMPETITORS[tool]["env"])


def all_tools() -> list:
    return ["chilmesh"] + sorted(COMPETITORS)

# Fair mode inventories all tools without changing the legacy W2 scheduling API.
FAIR_TOOLS = ('chilmesh', 'adcircpy', 'thalassa', 'xarray', 'oceanmesh', 'OceanMesh2D')
FAIR_OPERATIONS = ('parse', 'parse-netcdf', 'adjacency-cpp', 'adjacency-python',
                   'locate-tree-build', 'locate-query', 'render', 'stations', 'write', 'full-init')
NETCDF_ONLY_REASON = '{tool} reads NetCDF; native fort.14 work is not comparable'


def operation_contract(tool, operation, *, options=None):
    """Bind operation, backend, inputs and auxiliary/query/field options by hash."""
    import hashlib
    import json
    if tool not in FAIR_TOOLS or operation not in FAIR_OPERATIONS:
        raise ValueError('unknown fair tool/operation')
    backend = operation.removeprefix('adjacency-') if operation.startswith('adjacency-') else 'default'
    fmt = 'NetCDF' if operation == 'parse-netcdf' else 'fort.14'
    equivalence, reason = 'not-comparable', 'no equivalent public operation/output contract'
    if tool in ('oceanmesh', 'OceanMesh2D'):
        reason = 'W1-only tool; no W2 operation contract' + ('; Octave runtime is not Python equivalent' if tool == 'OceanMesh2D' else '')
    elif operation == 'parse' and tool in ('thalassa', 'xarray'):
        reason = NETCDF_ONLY_REASON.format(tool=tool)
    elif operation == 'parse' and tool in ('chilmesh', 'adcircpy'):
        equivalence, reason = 'like-for-like', 'requires verified actual native output and original ID mapping'
    elif operation == 'parse-netcdf' and tool in ('thalassa', 'xarray'):
        equivalence, reason = 'like-for-like', 'conversion / NetCDF only; verified fair-netcdf-v1 lineage required; no native ratio'
    elif tool == 'chilmesh' and operation != 'parse-netcdf':
        equivalence, reason = 'CHILmesh-only', 'no proven competitor work/output equivalence'
    contract = dict(tool=tool, operation=operation, backend_selected=backend, actual_backend=None,
                    input_format=fmt, options=dict(options or {}), equivalence=equivalence, reason=reason,
                    measured_work=operation, required_output='canonical mesh content' if operation.startswith('parse') else operation,
                    materialization='complete before END', prerequisites=[] if operation.startswith('parse') else ['verified parsed input'],
                    output_verification='actual output after END ACK')
    contract['tools_covered'] = ['chilmesh', 'adcircpy'] if operation == 'parse' and tool in ('chilmesh', 'adcircpy') else [tool]
    contract['canonicalization_revision'] = 'binary64-repr-v1'
    identity = {key: value for key, value in contract.items() if key not in ('tool', 'actual_backend', 'reason')}
    contract['contract_id'] = hashlib.sha256(json.dumps(identity, sort_keys=True, allow_nan=False).encode()).hexdigest()
    return contract


def operation_inventory(meshes):
    return [dict(mesh=mesh, **operation_contract(tool, op))
            for mesh in meshes for tool in FAIR_TOOLS for op in FAIR_OPERATIONS]
