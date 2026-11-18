<!-- provenance: author=domattioli model=unknown effort=medium date=2026-11-10 skill=speckit-plan repo=CHILmesh session=dispatch-sol-78a08f; manual stamp, exact provider model ID unavailable; transport rung=sol -->
# Implementation Plan: Fair Head-to-Head Benchmark

**Branch**: `feat/301-competitor-benchmarks` (existing; no change) | **Date**: 2026-11-10 | **Spec**: [spec.md](spec.md)
**Input**: `specs/007-fair-comparison/spec.md`, including Clarifications and Resolved decisions.
**Status**: Design complete; no implementation, measurement, rights approval, or publication performed.

## Summary

Extend the existing W2 subprocess harness with independently verified outputs, one region-aware RSS sampler, immutable schema-v2 evidence, and separate equivalence, measurement, stability, rights and claims gates. Keep every planned cell and failed attempt, including all six tools on STOFS. CHILmesh winning is not an acceptance condition. Replace no historical records and carry forward no historical memory claim.

## Technical Context

**Language/Version**: Python 3.11, one real base executable for parent and all Python venvs; installation compatibility remains a later preflight gate, not an assumed success.
**Primary Dependencies**: Standard-library subprocess, selectors, hashlib, importlib.metadata, statistics and Linux procfs/cgroup evidence; existing NumPy/CHILmesh workers, isolated adcircpy, Thalassa, xarray environments, and the independent stdlib+netCDF4 conversion script (research.md D7). No GPL imports into CHILmesh or the controller.
**Storage**: Append-only campaign manifests, mesh references, JSON repeat records, hashed sample/diagnostic sidecars, aggregate and claim decisions. References and sidecars are outside the record-directory `*.json` glob.
**Testing**: pytest, synthetic tiny fort.14 fixtures, injected clocks/procfs/cgroup readers/launchers; bounded Linux subprocess lifecycle tests. No real meshes or network for default tests.
**Target Platform**: One Linux 16-core Codespace; record actual 67424980992-byte-class RAM, kernel, CPU affinity/quota and memory limits rather than round labels. Host matching requires a frozen campaign identity and unchanged enforced allocation.
**Project Type**: Internal benchmark library plus CLI and a Codespace shell entry point; no public CHILmesh API changes.
**Performance Goals**: Auditable cost, not a target speedup; 10 ms shared RSS samples; at least five valid cold repeats per small-mesh cell and three per STOFS cell.
**Constraints**: Small meshes 600 s / 8 GiB RSS; STOFS 1800 s / 48 GiB RSS; serial children; OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS, NUMEXPR_NUM_THREADS all `1`. Sampled limits can miss sub-interval spikes.
**Scale/Scope**: Four meshes; six tools; parse, adjacency by selected backend, locate tree build, locate query, render, stations, write and full-init cells. Import/install diagnostics remain context, not new comparison operations. No W1/W3, backend fixes or #309.

## Constitution Check

| Principle | Before research | After design |
|---|---|---|
| I library first / isolated tests | Pass: reusable benchmark helpers | Pass: helpers under benchmarks/competitors, thin CLI |
| II immutability; V coordinates; IX mixed elements | Pass: no mesh mutation/API change | Pass: preserve IDs, winding, tri/quad arity and opaque coordinates |
| III tests first / regressions | Pass: test tasks precede implementation | Pass: synthetic adversarial checks plus legacy W1 regression isolation |
| IV correctness before performance | Pass: independent reference mandatory | Pass: invalid verification never yields ratios |
| VI formats; VII API stability; VIII documentation | Pass: benchmark-only changes | Pass: conversion separated; v1 evidence preserved and ineligible for v2 claims |
| X MATLAB parity | Not affected | Not affected: no geometry algorithms changed |

No constitutional exception is needed. Test-first implementation does not imply approval of any new public API. The constitution's historical CLAUDE.md reference is not an instruction to recreate that file.

## Project Structure

Documentation: `specs/007-fair-comparison/{plan,research,data-model,quickstart,tasks}.md` and `.specify/workflow-state.json`. CLI and record contracts are embedded in data-model.md because this dispatch permits no contracts/ files. Skip consumer agent-context updates; no DomI writes or final spec-move task.

Existing integration points:

- `benchmarks/competitors/measure.py`: replace the 20 ms wait-only loop with common 10 ms VmRSS sampling, marker handling, deadlines and process-group cleanup; retain wait4 only as a separately labelled cross-check.
- `records.py`: versioned validators/readers/writers, missing flags, evidence references and full matrix coverage instead of first matching metric.
- `machine.py`: host/runtime/build/load evidence; `valence.py`: retain Valence-only source rule and strengthen exact-version publication rights.
- `w2/run.py`: separate operation records rather than one bundled tool/repeat record; propagate adapter counts and checksum evidence, never `_mesh_meta` zeros; immutable series scheduling and explicit outcome rows.
- `w2/measures.py`, `adapters/{adcircpy_read,thalassa_read,xarray_read}.py`: common region protocol, complete materialization, distribution versions, output verification after region end.
- `w2/tools.py`: six-tool operation matrix; correct Thalassa to NetCDF-only; inventory oceanmesh/OceanMesh2D without creating W1 adapters.
- `cli.py`, `report/render.py`: campaign controls, aggregate/claims commands, rights-filtered reports with no excluded numeric data.

Proposed reusable modules: `reference.py` (independent fort.14 parser), `verification.py` (materialized-output normalization), `protocol.py` (standard-library-only child markers, importable with repo root on sys.path), `cache.py`, `campaign.py` (manifest, pins, schedule), `aggregate.py`, `claims.py`. Proposed independent converter: `benchmarks/competitors/scripts/convert_fair_netcdf.py`, stdlib+netCDF4 over reference-parser output with the reader-specific dimensions/variables and independent content verification in research.md D7. Proposed entry point: `benchmarks/competitors/scripts/run_fair_codespace.sh`; environment locks under `benchmarks/competitors/envs/`. Tests extend `tests/benchmarks/`.

## Phase 0: Research decisions

[research.md](research.md) resolves holes 7–12, fixes the interval and observation policies, and distinguishes decisions from untested host capabilities. No design clarification remains open. Exact dependency locks, tool API proof, cgroup access, cold-reset permission and rights citations are inputs to preflight; absent inputs block eligibility and produce explicit outcomes.

## Phase 1: Execution design

1. Prepare a frozen campaign manifest: host/allocation, one Python 3.11 digest and build, all dependency/native-runtime pins, sampler revision, 10 ms interval, limits, all cell contracts, repeat counts and deterministic order. Acquire only through Valence and independently build references before repeat execution. No source hash, reference parse or input conversion runs between a repeat's cache reset and its operation.
2. Classify each cell before launch. Native fort.14 parse can be like-for-like only for CHILmesh/adcircpy with proven equivalent outputs. Thalassa/xarray use conversion / NetCDF cells under research.md D7: an independent harness script, reader-specific dimension/variable contracts, verified source/converted content and hashes, and separate conversion/read/end-to-end costs. They never enter native fort.14 ratios. oceanmesh/OceanMesh2D are not-comparable W1-only inventory cells; Octave additionally fails Python parity. CHILmesh-only operations have explicit output obligations and backend identity; unavailable C++ remains explicit, with no backend repair.
3. Hold an exclusive campaign lock. Record one second of pre-run load observations, then `sync`, `/proc/meminfo` before, drop_caches write of `3`, and `/proc/meminfo` after. Cold requires successful commands and lower `Cached:`; otherwise cache unknown and invalid, not a trusted user label. Launch a fresh process session only after this flush. No application cache or prior operation result survives to a new child.
4. Competitor launches set PYTHONPATH to REPO_ROOT so adapters can import the standard-library-only `benchmarks.competitors.protocol`; GPL scan tests restrict adapter harness imports to that module. Child performs metadata collection, imports and allowed prerequisites without touching the parse input. It sends READY on a dedicated inherited control pipe. Parent obtains an immediate baseline sample, acknowledges GO, and continues fixed-period sampling. Child emits START and END with perf_counter timestamps bracketing the actual operation and full materialization. Parent timestamps receipt and uses same-host monotonic clock identity to classify samples by child timestamps. Child waits for END acknowledgement before checksum verification, so verification cannot contaminate the sampled windows. Missing/malformed markers, no pre-start baseline or no in-region samples invalidate memory. For non-parse operations, input parse and output-contract prerequisites precede START and must be verified; input parsing for parse never moves there. Each locate suboperation gets its own process/region rather than reusing an already built tree.
5. Parent records VmRSS bytes and descendants every 10 ms and load every 100 ms. Baseline is last sample before START; region peak is max within START..END; whole-process sampled peak is max launch..END; growth is max(0, region peak minus baseline). A gap greater than 20 ms across an observed window is flagged as incomplete coverage and invalidates memory; actual timestamps and scheduler jitter remain visible. Spikes shorter than 10 ms remain a stated limitation. Any descendant invalidates child-only memory equivalence, but all descendants are still killed on abort. `wait4.ru_maxrss` includes verification and is diagnostic only, never substituted for either sampled peak.
6. Shared safety watchdog measures wall from launch, uses 600/1800 s bounds, and samples total child RSS against 8/48 GiB. Imports and verification remain outside reported operation time but inside this process watchdog. A limit breach kills the fresh process group, reaps the child and retains partial evidence; RSS-limit termination is not asserted to be kernel OOM. cgroup memory.events before/after plus memory.peak where available support kernel OOM attribution; unexplained SIGKILL stays unknown. Unreadable cgroup evidence is explicit, prevents OOM attribution, and blocks campaign preflight for eligible measurements.
7. Verify every parse from the actual adapter materialized output after END. Independent reference uses raw-file header and binary64 canonical nodes/elements, not CHILmesh. Preserve arbitrary file IDs using explicit adapter mapping; storage offsets alone cannot establish IDs. Missing or ambiguous mapping is verification incomplete. `sys.modules` changes across region invalidate isolated operation comparison; preload only imports, then start a fully disclosed replacement series. Never warm the measured input by parsing it first.
8. Write each attempt atomically to a new, non-overwriting path. On controller interruption, finalize active attempts and synthesize reasoned not-run records for planned unstarted cells. Aggregate only named complete series; retain all original attempts and replacement rationales. All three comparison metrics must have relative spread <=10% for cell stability. Failed values remain null, not zero.
9. Public rendering first applies exact-version rights filtering, including mandatory WNAT_Hagen exclusion while no ruling exists. Then show permitted outcomes, raw repeat references, medians/spreads and caveats separately for native, conversion, CHILmesh-only and not-comparable sections. Claims come only from structured proposals passed through claims.py; generic user prose is not automatically certified.

## Validation and delivery

Implement tests before each production component; use mocks for `/proc`, cgroup, cache writes, metadata and tool APIs. Include real bounded process-group cleanup on Linux, without real meshes. End-to-end synthetic campaigns must reject corruption, deferred imports, runtime mismatch, contamination, cache failure, short/missing samples, mixed series, excluded rights and invalid claims. Recalculate aggregate/claim arithmetic from saved records offline. The later Codespace script defaults to preflight and requires explicit opt-in to execute; no campaign is run in this planning dispatch. [quickstart.md](quickstart.md) supplies future commands and expected outcomes.

## Complexity Tracking

No violations. Dedicated control pipe and hashed sample sidecars add bounded complexity needed to exclude post-region verification and retain raw evidence without truncating JSON to stdout tails. Keep v1 workflow compatibility in a separate validator/legacy path; never silently upgrade historical records to eligible v2 evidence.
