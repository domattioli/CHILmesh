<!-- provenance: author=domattioli model=claude-sonnet-5-5 effort=medium date=2026-11-08 skill=speckit-plan repo=CHILmesh session=2d051659 -->
# Implementation Plan: Competitor Head-to-Head Benchmarks (W2 then W1)

**Branch**: current (no switch) | **Date**: 2026-11-08 | **Spec**: [spec.md](spec.md)

## Summary

Portable, CLI-driven harness under `benchmarks/competitors/` (excluded from sdist and wheel). Phase 1 builds the core and W2. Phase 2 runs W1 (Track 1 headline, Track 2 secondary). Phase 3 makes doc fixes 1, 4, 5. Phase 4 reports and closes. Meshes come only through the Valence client. GPL tools run only as subprocesses. Every cell is one JSON run record, status `ok | failed | not-run | invalid`. Real runs wait on Valence serving the meshes.

## Research Findings (evidence is file:line; unknowns marked)

- R1. Valence client exists, stdlib only: `src/chilmesh/registry.py` (`list_meshes`, `fetch_hash`, `fetch_mesh_bytes`, `download_mesh`, `fetch_mesh`, L152-230). Needs `VALENCE_TOKEN`; base URL via `VALENCE_API_URL`. SHA-256 verified per download (L174-183). Errors: `RegistryNotFoundError`, `RegistryUnavailableError`, `RegistryAuthError`.
- R2. No provenance or rights endpoint (only `/meshes`, `/meshes/{n}`, `/meshes/{n}/hash`). `fetch_mesh_bytes` holds the whole payload in RAM; 1.73 GB STOFS needs a streaming download in the harness (not a change to `registry.py`). Whether Valence exposes source/version/rights metadata: unknown.
- R3. `docs/BENCHMARK.md:259,272-276`: STOFS known out-of-environment (needs hosted runner, Valence #77). Whether Valence now serves STOFS: unknown.
- R4. `benchmarks/` does not exist; `scripts/benchmark*.py` are CHILmesh-vs-itself single files. `scripts/benchmark_wnat_hagen.py:4-9,48-60` reads a local Valence clone (`VALENCE_DATA_DIR`); the new harness must NOT copy that (conflicts with "Valence API only"). `MANIFEST.in:28` has `prune scripts`; add `prune benchmarks`. No existing parse or render benchmark (`docs/BENCHMARK.md:47,79-86`).
- R5. #276 / FR-011: `read_fort14_raw` uses `fh.readlines()` (`src/chilmesh/fort14_io.py:298`) and dict-of-tuples (L105-159). At 12.8M nodes / 24.9M elements this likely exceeds 32 GB; measured peak unknown. ENPAC parse 2.11 s at 273k nodes (`BENCHMARK.md:114`); scaled ~100 s is an estimate only. Verdict comes from W2. Reader change is a follow-up PR, not part of this feature.
- R6. `find_element`: cKDTree over centroids, k=20 candidates, Python point-in-element loop (`src/chilmesh/CHILmesh.py:1019-1043`); the tree is built in the constructor when `build_spatial_indices` is true (L420-422, L971-974), not lazily; `find_element` only queries it. Report tree build separately from per-query time; use inside-mesh points or label outside-mesh cost.
- R7. `courant_number` / `cfl_gate` at `src/chilmesh/quality.py:360-470`. "Largest stable dt" is derived: `dt_max = courant_max / max C(dt=1)` (C is linear in dt); helper lives in the harness. CRS spherical for ENPAC/WNAT; Shinnecock units unknown.
- R8. Smoothers: `smooth` (`CHILmesh.py:1623`, `fixed_boundary=True`), `smooth_mesh` (L1672), `angle_based_smoother` (L1773), `direct_smoother` (L2276). Constitution II: coordinates move, topology constant, so Track 1 is natural for CHILmesh.
- R9. Doc defects: `BENCHMARK.md:222` (WNAT "~4M DOF", true ~105k) and `:245` (STOFS "~4M DOF", true ~25.6M); `README.md:8` range 8.6x-14.7x mixes meshes (14.7x WNAT, 8.6x ENPAC; v1.4.1 ENPAC ~14.8x per `BENCHMARK.md:36-38,51`); Python 3.10 at `README.md:27` and `:46` (badge L18, req L62, `pyproject.toml:12` already 3.11+). `BENCHMARK.md:5` header says v1.2.2 (stale, out of scope).
- R10. Constitution gates that bite: I, II, III, IV, V, VII, VIII; Governance "Medium decision" may need a `design` issue. License: PolyForm Noncommercial (`pyproject.toml:8`); GPL stays separate by process boundary only.
- R11. Unknown from repo: Octave + OceanMesh2D compatibility; adcircpy/Thalassa licenses and versions; oceanmesh frozen-topology mode; Valence rights flags; source of `points.cwl.nc` and CO-OPS station list.

## Technical Context

Python 3.11+ harness; stdlib + numpy + psutil/resource; chilmesh + chilmesh-cpp. Competitors in separate envs (never in CHILmesh `pyproject.toml`). Meshes cached and JSON records written outside the repo; only reduced tables committed. Tests: pytest, competitor-free (schema, status, displacement check, scaling fit, GPL scan); real runs opt-in behind `CHILMESH_RUN_BENCH=1`. Target: remote Linux ~32 GB. Repro within 10 percent (SC-003): median of 3 repeats (N=1 for STOFS parse, flagged). Cold cache attempted and recorded (`cold | warm | unknown`).

## Constitution Check

Pass with conditions: I (unit tests on synthetic meshes, real runs opt-in), III (tests first), V (record CRS per mesh), VII (no public API change), VIII (docs match data). Governance: classified Medium (adds a top-level `benchmarks/` tree, no public API change). Operator exception 2026-11-08: no `design` issue and no 3-day feedback window for this feature. Constitution II/III: metric and displacement tests cover tri, quad and mixed meshes and compare boundary coordinates with `np.array_equal`. FR-007/SC-005: subprocess only, scan test, `prune benchmarks`. No violations.

## Project Structure

```text
specs/006-competitor-benchmarks/{spec,plan,research,data-model,quickstart,tasks}.md, contracts/, checklists/
benchmarks/competitors/
  README.md  cli.py  records.py  machine.py  valence.py  measure.py
  metrics/{quality,courant,displacement}.py
  w2/{tools,scaling,run}.py   w1/{tools,track1,track2,run}.py
  adapters/{adcircpy_read,oceanmesh_clean,thalassa_read,xarray_read}.py  oceanmesh2d.m
  report/render.py   envs/requirements-*.txt
tests/benchmarks/
docs/BENCHMARK.md  README.md  MANIFEST.in   # Phase 3 edits
```

Decision: `benchmarks/competitors/`, not `scripts/` (needs package, adapters, schema). GPL-importing adapters run only via subprocess in the competitor's env; chilmesh never imports them. Whether adapter scripts count as "GPL code in the repo" is a licensing question to confirm.

## Design

Run record fields: schema_version, tool, tool_version, license, run_mode, mesh (name, version, valence sha256, provenance reference, crs and units, n_nodes, n_elems), workflow, track, metric_values (null + `na: true` for not applicable), peak_rss_bytes, wall_s, repeat_index, cache, machine, command, chilmesh_git_sha, status, reason, stdout/stderr tail. Renderer emits a cell for every tool x metric and fails loudly on a missing one (SC-001). Display labels "not run" and "not applicable" map to status `not-run` and `na: true`; a justified `na` cell satisfies SC-001.

Cache policy: the runbook gives cold-cache steps; the run records `cache` as `cold | warm | unknown`. Only `cold` runs enter the published comparison; others are listed flagged. Each cell repeats 3 times (STOFS parse N=1, flagged). Timing boundaries for point location: tree construction (constructor) is timed apart from the 1,000 `find_element` queries, never both in one number.

Measurement: fresh subprocess per tool and metric; peak RSS via `resource.getrusage` or `/usr/bin/time -v`; timeout and OOM become `failed` records with peak reached.

Valence: wrapper over `chilmesh.registry` plus streaming SHA-256 download. Not found or unavailable gives a `not-run` record "mesh not yet served by Valence". No local-clone or NOAA fallback. Rights: manual `rights.json` in results dir; meshes without a permitting entry are excluded and listed (FR-008).

### Phase 0: resolve unknowns, write research.md
Confirm Valence holdings and metadata; source of `points.cwl.nc` and CO-OPS stations; competitor licenses and versions; Octave probe; default smoother. Exit gate: each competitor row is runnable, not-run-with-reason, or needs-decision.

### Phase 1: W2 harness (first, FR-001)
Tests first (schema, status, machine capture). Implement records, machine, measure, valence, renderer. Adapters: CHILmesh `read_from_fort14`, adcircpy, Thalassa, xarray/netCDF4 (input format per tool recorded). Adjacency C++ vs Python (confirm active backend via `backend_info()`). 1,000-query point location with build reported separately, outside-mesh count. Station sampling vs `points.cwl.nc` (max abs diff). Render time of max-elevation field (timeout and `failed` acceptable). Scaling test at >=3 sizes (Shinnecock 5,780; WNAT 98,365; ENPAC 531,680; STOFS 24.9M; optional synthetic): fit time vs log n and vs n, O(log n) supported only if log fit wins; note cache effects. Publish W2 section with #276 verdict. Gate: SC-001.

### Phase 2: W1 (after W2 published)
Track 1 (headline): frozen topology; CHILmesh smoother, OceanMesh2D via Octave, oceanmesh subprocess; Shinnecock, WNAT_Hagen, ENPAC2003; metrics per spec US2.1 computed by one shared `metrics/quality.py` on each tool's output; boundary displacement exactly 0 else `invalid` (also assert node and element counts unchanged). Track 2 (secondary): full clean, plus vertices deleted; note CHILmesh deletes none. Octave failure gives `not-run` with stderr tail. Gate: SC-002.

### Phase 3: doc fixes (independent)
Fix 1: `BENCHMARK.md:222` reword (not only number swap; direct solve is feasible at ~105k DOF, check L263-270), `:245` to ~25.6M DOF. Fix 4: `README.md:8` cite mesh and version per figure. Fix 5: remove `README.md:27` callout, rewrite `:46`, grep other "3.10", add a grep test. Gate: SC-006.

### Phase 4: report and closure
Render tables, add W2/W1 sections to `BENCHMARK.md`, GPL scan of repo, wheel and sdist (SC-005), repeat one cell within 10 percent (SC-003), verify the FR-014 W3 exclusion (no W3 code or doc), update #276 and #301, then move specs to DomI as the final step.

## Risks
Valence not serving meshes (runs block; unit tests proceed on synthetic meshes). Harness RAM from whole-file download (streaming fix). STOFS adjacency or render may OOM or take hours (timeouts, `failed`, N=1). Metric mismatch across tools (one shared implementation). Format rounding may break displacement == 0. License boundary is process-only.

## Decisions (2026-11-08, operator)

- D1. NOAA data (`points.cwl.nc`, CO-OPS station list): "Valence only" covers meshes. The harness fetches these directly from NOAA, recording URL, SHA-256 and fetch date. Separately, file an issue to decide whether they belong in domain data (Valence, or a `.chil` extension table in CHILmesh). Ownership per memory: CHILmesh owns the container and extension tables; Valence owns fingerprints. File where it fits after reading both repos.
- D2. W1 reports all three CHILmesh smoothers (`smooth` with `fixed_boundary=True`, `angle_based_smoother`, `direct_smoother`) on each Track 1 and Track 2 mesh.
- D3. README claim (fix 4) is about instantiation and mesh operations, not smoothing. Rewrite `README.md:8` to say exactly what was measured, with mesh and version per figure; do not reuse the blended 8.6x-14.7x range. Confirm the source of each figure in `docs/BENCHMARK.md` before wording.

## Open Questions (remaining)
1. Does Valence serve all four meshes, and what provenance and rights metadata?
2. oceanmesh and OceanMesh2D function mapping for Tracks 1 and 2.
3. Fixed dt, depth source, Shinnecock CRS for Courant.
4. How is displacement == 0 checked across competitor file formats?
5. Does SC-006 cover CHANGELOG history?
6. Do GPL-importing adapter scripts count as GPL code in the repo; is a design issue needed?
7. What does "raw xarray/netCDF4" read, fort.14 or NetCDF?
