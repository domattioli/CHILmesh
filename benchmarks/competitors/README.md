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

## Fair comparison (schema v2)

The legacy commands above preserve historical W1/W2 behavior. Their cache labels,
first-value tables and schema-v1 records cannot establish fair-comparison eligibility.
Use the [fair-comparison quickstart](../../specs/007-fair-comparison/quickstart.md)
for the offline checks and operator handoff. No real campaign or publication is
established by synthetic tests.

```bash
python -m benchmarks.competitors fair-preflight --manifest campaign.json --out results/fair
# Later, with operator authorization and all prerequisites satisfied:
CHILMESH_RUN_BENCH=1 python -m benchmarks.competitors fair-run --manifest campaign.json --out results/fair
python -m benchmarks.competitors fair-report --campaign results/fair
python -m benchmarks.competitors fair-report --campaign results/fair --published --rights rights.json --claims claims.json
```

Preflight is read-only: no downloads, measured children or cache reset. It checks
the frozen manifest, physical locks, common Python 3.11 realpath/digest, host and
allocation, readable cgroup OOM counters, reset permission and load observability.
Real exact direct/transitive/native locks and successful install/import compatibility
must be supplied; placeholder locks block execution. Every environment and the
parent use the same declared base Python. Four thread variables (OMP, OPENBLAS,
MKL and NUMEXPR thread counts) must equal 1.

Preparation is explicit. Beside `campaign.json`, provide `evidence/prepared.json`
with `identities` and `commands` keyed by frozen repeat ID, `references` keyed by
manifest mesh alias, and optional `conversions` keyed by repeat ID. A tool-keyed
conversion remains supported for a single source; repeat keys bind distinct sources
in multi-mesh campaigns. These are the `run_fair_campaign` API inputs, with full
RepeatRecord identity and MeshReference/ConversionReference objects, not invented
pins or source-count defaults. Acquire inputs only through `valence.fetch_mesh_file`
and `reference_from_fetch`; use `campaign.prepare_conversion` for independent
reader-specific conversion. Prepare references, inputs, queries and conversions
before scheduling. The CLI does not download a replacement source or invent
preparation when this file is absent: it writes the full scheduled matrix as
reasoned not-run/not-comparable outcomes and exits 1.

For the Codespace wrapper, declare the same base interpreter and default to preflight:

```bash
CHILMESH_FAIR_PYTHON=/absolute/base/python3.11 bash benchmarks/competitors/scripts/run_fair_codespace.sh --manifest campaign.json --out results/fair
CHILMESH_RUN_BENCH=1 CHILMESH_FAIR_PYTHON=/absolute/base/python3.11 bash benchmarks/competitors/scripts/run_fair_codespace.sh --run --manifest campaign.json --out results/fair
```

The wrapper never installs packages or publishes. It serializes output controllers;
Python also holds a host-wide campaign lock across output directories. Run mode
retains `evidence/`, immutable `attempts/`, `raw/` sample/diagnostic files and a new
`evidence-archive.tar.gz`. Existing attempts and archives are never overwritten.
A replacement requires a new manifest/series/output, the superseded series ID and
a rationale; no success selection across series. Five repeats are required for
small meshes, three for STOFS. A failed or invalid attempt makes its series incomplete;
all three metrics must have range/median <=10% to support claims.

Each fresh child follows one protocol: a one-second load baseline; `sync` and write
`3` to drop_caches with ordered status/timestamp and lower `Cached:` evidence; launch;
imports; READY/GO; START/materialized operation/END; ACK; output verification.
A `--cache cold` label cannot satisfy this protocol. Load is observed every 100 ms;
unrelated CPU >5 percentage points of one core or I/O >1 MiB/s invalidates evidence.
Unknown load/cache/runtime evidence also invalidates comparison.

Child-only VmRSS samples use a 10 ms cadence; gaps >20 ms, descendants or missing
baseline/region samples invalidate memory. Whole-process sampled peak covers launch
through END, including imports. Growth is max(0, region peak minus pre-START baseline),
an upper bound affected by allocator retention and mapped pages, not total mesh
allocation. Shorter spikes can be missed. Verification is outside both windows;
wait4 lifetime peak is a separate diagnostic. Watchdogs cover startup through exit:
600 s / 8 GiB for small meshes, 1800 s / 48 GiB for STOFS. RSS-limit, timeout,
cgroup-evidenced kernel OOM and unknown signal stay distinct failures. Interruptions
retain active partial evidence and every unstarted cell.

Fair reports recompute saved hashes, full declared schedule coverage and exact
series summaries offline. Reference metadata and sample/stdout/stderr sidecars are
hash checked before rights filtering. Native fort.14, conversion / NetCDF,
CHILmesh-only and not-comparable sections remain distinct. Thalassa/xarray have no
native fort.14 comparison; oceanmesh/OceanMesh2D remain W1-only with explicit reasons.
NetCDF conversion cost, measured read cost and their end-to-end sum stay separate.

All public results default to denied rights, including WNAT_Hagen and the other
three meshes. An entry must transcribe an actual permitting Valence ruling; the
following illustrates its format and grants no permission:

```json
{
  "STOFS": {
    "mesh": "STOFS",
    "exact_version": "EXACT-SERVED-VERSION",
    "publish_results": true,
    "citation": "https://github.com/Valence/Valence/issues/123",
    "ruling_date": "2026-11-10",
    "scope": "publication of benchmark results",
    "evidence_sha256": "SHA256-OF-RETAINED-AUTHORITATIVE-RULING"
  }
}
```

Citations identify a GitHub Valence repository/owner issue or PR (`pull/N`); review
its authority and retain the ruling evidence hash. Missing, denied, ambiguous,
undated or version-mismatched rulings exclude all measurements, including failure
numbers and claim diagnostics. Public exclusions contain only mesh names/reasons.
Admitted meshes show the citation/date/scope/hash. Legacy `apply_rights` remains
unchanged for historical workflows.

Claims are a JSON list of proposals with `claim_id`, `template: "less"`, numerator
and denominator cell IDs, `mesh`, `exact_version`, `operation`, `metric`, `unit`,
`direction: "lower"`, `formula: "100*(B-A)/B"` and
`resource_condition_reference` equal to the manifest SHA-256. Supported metrics are
`operation_elapsed_s` (seconds), `rss_whole_process_peak_bytes` and `rss_growth_bytes`
(bytes). The denominator is B and must be positive. Both cells must be native,
like-for-like, verified, complete, stable across every metric, condition-valid and
rights-cleared. Accepted bounded sentences cite versions, medians, spreads, repeats,
series, resource scope and raw hashes. Free prose and historical memory claims
cannot be certified. Exit 0 means valid command artifacts; recorded tool failures
alone are not controller errors. Exit 2 means malformed arguments/missing opt-in;
exit 1 means preflight, coverage, integrity, publication or claim-gate failure.
