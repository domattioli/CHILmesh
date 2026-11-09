# Competitor benchmark harness

Measures CHILmesh against adcircpy, thalassa, xarray, oceanmesh and OceanMesh2D.
The harness is not part of the sdist (`prune benchmarks`) and needs Python 3.11 or newer.

Two workflows run:

- **W2**: mesh read, mesh construction with adjacencies and layers (C++ vs Python backend), point location and render.
  The metrics `adj_cpp_s` and `adj_python_s` time construction with adjacencies and layers. Building the adjacencies is pure
  Python on every backend; the C++ backend only accelerates layerization.
- **W1**: pre-run mesh QA and repair. Track 1 freezes topology and pins the boundary. Track 2 allows a full clean.

## Setup

1. Install CHILmesh in the main environment: `pip install -e .`
2. Create one environment per competitor from `envs/requirements-<tool>.txt`. Competitors are never added to `pyproject.toml`.
   GPL tools (adcircpy, oceanmesh) run only as standalone subprocess scripts from `adapters/`. CHILmesh never imports them.
3. Point the harness at each environment's interpreter:
   `CHILMESH_BENCH_PY_ADCIRCPY`, `CHILMESH_BENCH_PY_THALASSA`, `CHILMESH_BENCH_PY_XARRAY`, `CHILMESH_BENCH_PY_OCEANMESH`.
   A tool with no variable set is recorded as `not-run`, with the reason.
   OceanMesh2D runs through GNU Octave and needs `octave` on `PATH`.
4. Set `VALENCE_TOKEN`. Meshes come only from Valence through `chilmesh.registry`.
   A mesh that Valence does not serve yet is recorded as `not-run`: "mesh not yet served by Valence". There is no local or NOAA fallback.
5. Create `rights.json` by hand. A mesh with no entry that allows publishing is excluded from the published tables.

## Run

Real runs are opt-in. Without `CHILMESH_RUN_BENCH=1` the commands exit with code 2.

```
CHILMESH_RUN_BENCH=1 python -m benchmarks.competitors w2 --mesh WNAT_Hagen --out results/w2 --cache cold
CHILMESH_RUN_BENCH=1 python -m benchmarks.competitors w1 --mesh Shinnecock --out results/w1
python -m benchmarks.competitors report --records results/w2 --workflow w2 --published --rights rights.json
```

## NOAA auxiliary inputs (`--aux-dir`)

Station sampling and the render metric read these files from the `--aux-dir` directory:

| File | Content |
|---|---|
| `coops_stations.json` | CO-OPS station list (`{"stations": [{"lat", "lng"}, ...]}`) |
| `points.cwl.nc` | Station time series. Variables `x`, `y`, `zeta`. The reference per station is the max over time of `zeta`. |
| `maxele.63.nc` | STOFS max-elevation field, one value per mesh node, variable `zeta_max` |

`station_max_abs_diff` is the max absolute difference between the `maxele.63.nc` field sampled at the stations
(inverse-distance over the containing element) and the peak `zeta` from `points.cwl.nc`. The render metric draws
the same field. Mesh nodal z (bathymetry) is never used for either.

- `coops_stations.json` or `points.cwl.nc` absent: station sampling is `not-run`, "NOAA inputs not fetched".
- `maxele.63.nc` absent: station sampling is `not-run`, "NOAA inputs not fetched (maxele)". The render cell is `not applicable` with the same reason.
- `maxele.63.nc` present but without `zeta_max`, or with a node count that differs from the mesh: the record is `failed`, with the reason.

## Cold cache

A published timing must be a cold-cache timing. Before each run, flush the file cache:

- Linux: `sync; echo 3 | sudo tee /proc/sys/vm/drop_caches`
- macOS: `sudo purge`

Pass `--cache cold` only after flushing. Warm and unknown runs are flagged and left out of the published table.

## Reproduction by a second person

A result counts as reproduced when a second person, on their own machine, follows this runbook and gets records
with the same status for each tool and mesh. Each record stores the machine description, tool versions, the mesh
SHA-256 and the CHILmesh git SHA, so differences can be traced.

## Records

Each run writes one JSON record (`schema_version` 1) with status `ok`, `failed`, `not-run` or `invalid`.
A cell that a tool cannot produce is `not applicable` with a stated reason. The report step fails if any
tool x metric cell has no value, no failure reason and no justified `not applicable`.
