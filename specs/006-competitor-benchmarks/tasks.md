<!-- provenance: author=domattioli model=claude-sonnet-5-5 effort=medium date=2026-11-08 skill=speckit-tasks repo=CHILmesh session=2d051659 -->
# Tasks: Competitor Head-to-Head Benchmarks (W2 then W1)

**Input**: spec.md, plan.md (this directory). Tests are included (plan: tests first).
**Format**: `[ID] [P?] [Story] Description` with exact paths. `[P]` = parallel (different files). T001-T006 run serially (all append to `research.md`).
**Order**: US1 (W2) before US2 (W1) per FR-001. US3 (docs) is independent and can run any time.
**Real runs** (STOFS, Shinnecock, WNAT, ENPAC) are opt-in (`CHILMESH_RUN_BENCH=1`) and wait on Valence. Tasks marked RUN need the remote machine and are not part of the code-complete gate.

## Phase 0: Research (blocks design details)

- [ ] T001 Ask Valence which of Shinnecock, WNAT_Hagen, ENPAC2003, STOFS v3.1 it serves and what provenance/rights metadata exists; write `specs/006-competitor-benchmarks/research.md`
- [ ] T002 Record license and version for adcircpy, Thalassa, oceanmesh, OceanMesh2D in `research.md`
- [ ] T003 Probe GNU Octave + OceanMesh2D on a small mesh; record result in `research.md`
- [ ] T004 Decide fixed dt, depth source, Shinnecock CRS for Courant; record in `research.md`
- [ ] T005 Map oceanmesh and OceanMesh2D calls to Track 1 (frozen topology) and Track 2 (full clean) in `research.md`
- [ ] T006 Read Valence and CHILmesh `.chil` ownership notes; file one issue on whether NOAA `points.cwl.nc` and CO-OPS station list belong in domain data (Valence or a `.chil` extension table). Draft via `write-like-scientist`; show draft to operator before posting (decision D1)

## Phase 1: Setup

- [x] T007 Create `benchmarks/competitors/` tree per plan (README.md, cli.py, records.py, machine.py, valence.py, measure.py, metrics/, w2/, w1/, adapters/, report/, envs/) and `tests/benchmarks/`
- [x] T008 [P] Add `prune benchmarks` to `MANIFEST.in`
- [x] T009 [P] Add per-tool `benchmarks/competitors/envs/requirements-*.txt` (competitor deps never in `pyproject.toml`)

## Phase 2: Foundational (blocks US1 and US2)

- [x] T010 [P] Tests `tests/benchmarks/test_records.py`: run-record schema (incl. mesh version, provenance reference, CRS and units), status enum, `na: true` metrics, loud failure on a missing tool x metric cell, justified `na` accepted (SC-001)
- [x] T011 [P] Tests `tests/benchmarks/test_measure.py`: subprocess runner records wall time, peak RSS, timeout and OOM as `failed` with cause
- [x] T012 [P] Tests `tests/benchmarks/test_valence.py`: not-found / unavailable gives `not-run` "mesh not yet served by Valence"; no local-clone or NOAA mesh fallback
- [x] T013 Implement `records.py` (mesh version, provenance, CRS fields; cache `cold|warm|unknown`), `machine.py` (CPU, RAM, OS, python, git sha, tool versions), `measure.py` to pass T010-T011
- [x] T014 Implement `valence.py`: wrapper over `chilmesh.registry` plus streaming SHA-256 download (no whole-file RAM copy); rights via manual `rights.json`, excluded meshes listed (FR-008)
- [x] T015 [P] Implement `report/render.py`: tables per workflow and track, failed and not-run rows with reasons
- [x] T016 [P] Test `tests/benchmarks/test_gpl_scan.py`: no `import oceanmesh|adcircpy` outside `adapters/`, no GPL text in repo (SC-005)
- [x] T017 [P] Implement `noaa.py`: direct NOAA fetch for `points.cwl.nc` and station list, recording URL, SHA-256, fetch date (D1)

**Checkpoint**: T010-T016 green on a synthetic mesh.

## Phase 3: User Story 1 - W2 STOFS post-processing (P1)

- [x] T018 [P] [US1] Tests `tests/benchmarks/test_w2_scaling.py`: log-n versus linear-n fit picks correctly on synthetic timings
- [x] T019 [P] [US1] Tests for point-location harness: inside-mesh points, outside-mesh count reported, tree construction (constructor) timed apart from 1,000 `find_element` queries, no double counting
- [x] T020 [US1] Implement `w2/tools.py` adapters: CHILmesh `read_from_fort14`, plus subprocess adapters `adapters/adcircpy_read.py`, `adapters/thalassa_read.py`, `adapters/xarray_read.py`; record input format per tool
- [x] T021 [US1] Implement adjacency build timing C++ and Python; confirm active backend via `backend_info()`
  - Note: time_adjacency now includes layerization; adjacency alone is pure Python on every backend
- [x] T022 [US1] Implement 1,000-query point-location benchmark at CO-OPS stations (FR-002) (wired; real run pending T027)
- [x] T023 [US1] Implement station sampling versus `points.cwl.nc`, max absolute difference (wired; real run pending T027). Note: samples `maxele.63.nc` `zeta_max` against peak `zeta`, not nodal bathymetry; also T022 render metric uses that field
- [x] T024 [US1] Implement `plot_field` render timing with timeout and `failed` record
- [x] T025 [US1] Implement `w2/scaling.py`: query time at 3 or more sizes (Shinnecock, WNAT_Hagen, ENPAC2003, STOFS, optional synthetic); fit and report both (FR-013, SC-004) (wired; real run pending T027)
- [x] T026 [US1] Implement `w2/run.py` and CLI `python -m benchmarks.competitors w2`
- [ ] T027 [US1] RUN on remote machine: W2 on STOFS v3.1 (cold cache verified and recorded, N=1 parse flagged, N=3 others); non-cold runs are flagged and excluded from the published comparison; write records
- [ ] T028 [US1] Publish W2 section in `docs/BENCHMARK.md`, including #276 verdict on `fort14_io.py:298` (FR-011); if a reader change is needed, file a follow-up issue, do not change the reader here

**Checkpoint**: SC-001, SC-004. W2 results published before any W1 result (FR-001).

## Phase 4: User Story 2 - W1 mesh QA, repair, timestep gate (P2)

- [x] T029 [P] [US2] Tests `tests/benchmarks/test_displacement.py` on tri, quad and mixed synthetic meshes: boundary compared with `np.array_equal`, exactly 0 passes, any move gives `invalid`; node and element counts unchanged in Track 1 (FR-006)
- [x] T030 [P] [US2] Tests `tests/benchmarks/test_metrics.py`: shared quality metrics (min, mean minus 3 sigma versus 0.75, percent angles under 10 degrees, percent irregular valence) and `dt_max = courant_max / max C(dt=1)` on tri, quad and mixed hand-checked meshes
- [x] T031 [US2] Implement `metrics/quality.py`, `metrics/courant.py`, `metrics/displacement.py`
- [x] T032 [US2] Implement `w1/tools.py`: CHILmesh with all three smoothers (D2); `adapters/oceanmesh_clean.py` (subprocess); `adapters/oceanmesh2d.m` via Octave, `not-run` plus stderr tail on failure
- [x] T033 [US2] Implement `w1/track1.py` (headline, frozen topology) and `w1/track2.py` (full clean, vertices deleted, note CHILmesh deletes none)
- [x] T034 [US2] Implement `w1/run.py` and CLI `python -m benchmarks.competitors w1`
- [ ] T035 [US2] RUN: Track 1 then Track 2 on Shinnecock, WNAT_Hagen, ENPAC2003 once Valence serves them
- [ ] T036 [US2] Publish W1 section in `docs/BENCHMARK.md`; Track 1 as headline (SC-002)

## Phase 5: User Story 3 - Doc fixes (P3, independent)

- [x] T037 [P] [US3] `docs/BENCHMARK.md:222`: WNAT is 52,774 nodes, about 105k DOF; reword sentence (direct solve is feasible at this size; check L263-270). `:245`: STOFS about 25.6M DOF
- [x] T038 [P] [US3] `README.md:8`: rewrite speedup claim to what it measures (instantiation and mesh operations), one mesh and version per figure, sourced from `docs/BENCHMARK.md:36-38,51` (D3)
- [x] T039 [P] [US3] `README.md:27` remove 3.10 callout, `:46` rewrite; grep README, docs, CONTRIBUTING, workflows for other "3.10"; add a grep test (SC-006); decide CHANGELOG scope
- [x] T040 [US3] Verify each edited figure against `docs/BENCHMARK.md` rows

## Phase 6: Polish

- [x] T041 [P] `benchmarks/competitors/README.md` runbook for the remote machine (env setup, `VALENCE_TOKEN`, cold-cache steps and verification, commands)
- [ ] T042 GPL scan of repo, built wheel and sdist (SC-005)
- [ ] T043 Repeat published cells from their recorded commands, median of 3, confirm within 10 percent (SC-003); document a second-person repro procedure in the runbook and ask a second person to repeat one cell; state the STOFS N=1 exception
- [ ] T044 Verify FR-014: grep shows no W3 (partition, tidal hindcast) code or doc added; update #276 and #301 with verdicts and links
- [ ] T045 (plan Phase 4 final step) Move specs to DomI: `bash <DomI>/skills/speckit-pipeline/scripts/move_specs_to_domi.sh --apply` (last task)

## Governance

Operator exception 2026-11-08: no `design` issue or 3-day window for this feature (plan Constitution Check).

## Dependencies

T001-T005 inform T031-T033 details. T007 before T010. T013-T014 before US1/US2. US1 before US2 publication. T037-T040 independent. T027, T035 need Valence and the remote machine.
