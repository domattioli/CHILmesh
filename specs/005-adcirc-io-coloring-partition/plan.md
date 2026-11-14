<!-- provenance: author=domattioli model=claude-opus-5-5 effort=medium date=2026-10-12 skill=speckit-pipeline repo=CHILmesh session=4e33ac09-e474-4ebc-a2dd-027a9e0a0442 -->
# Implementation Plan: ADCIRC I/O, element coloring, and mesh partitioning (#289 B1-B3)

**Branch**: `development` (no feature branch, single-branch policy) | **Date**: 2026-10-12 | **Spec**: [spec.md](spec.md)
**Input**: `specs/005-adcirc-io-coloring-partition/spec.md` (clarified 2026-10-12)

## Summary

Add to CHILmesh, in three ordered stages: (B1) a public `.2dm` writer with nodestrings, node-sharing element coloring, and a partition quality report, plus a partitioner only if SegmentAnyMesh's seeded watershed fails the balance gate; (B2) readers for ADCIRC ASCII and NetCDF output and a field-on-mesh plot; (B3) fort.22 best track and OWI readers and writers, and fort.19/20/24 readers and writers. A local, non-CI script checks one pinned real file per format against the readers. Each stage closes only on its own gate.

## Technical Context

**Language/Version**: Python >=3.10 [code pyproject.toml:8]
**Primary Dependencies**: numpy>=1.23, scipy>=1.10, matplotlib>=3.6 (existing required) [code pyproject.toml:37-39]; new optional extras `netcdf = ["netCDF4"]` and, only if the CHILmesh partitioner path is taken, `partition = ["pymetis"]`
**Storage**: Files (ADCIRC text formats, SMS .2dm, NetCDF)
**Testing**: pytest; synthetic fixtures written by test helpers; real-file gate is a standalone script, not collected by pytest
**Target Platform**: Linux, macOS, Windows (pure Python)
**Project Type**: Library
**Performance Goals**: Single-step read of a 1M-node global output file within 3.0x that step's array `nbytes`, peak `tracemalloc` (SC-005)
**Constraints**: No network in the test suite; no third-party ADCIRC file in the repo; importing chilmesh never needs netCDF4 or pymetis
**Scale/Scope**: About 8 new public function groups, 7 new modules, 2 scripts (gate, exploration driver), 1 manifest

## Constitution Check

| Principle | Status | Note |
|---|---|---|
| I Library-first | PASS | Each format is its own module with its own public functions, re-exported from `chilmesh/__init__.py` |
| II Mesh immutability | PASS | Coloring, partitioning and readers do not mutate the mesh |
| III Test-first | PASS | Tasks put tests before implementation for every public function |
| IV Correctness over speed | PASS | Coloring is greedy and checked; balance threshold explicit |
| V Coordinate agnosticism | PASS | Readers return coordinates as read; OWI grids keep their own lat/lon fields as data |
| VI Format pluralism | AMEND (no MUST conflict) | VI's 2DM sentence is descriptive state, not a normative rule: it records that write omits nodestrings (#228). This feature adds them, which advances the principle. The sentence is updated when US1 ships (task in B1 gate) so docs match code (VIII) [code .specify/memory/constitution.md VI] |
| VII API stability | PASS | `_write_2dm` stays as alias of the new public method (FR-003); new names are additive |
| VIII Docs = contract | PASS | Each public function gets a docstring with a tested example |
| IX Mixed elements | PASS | .2dm, coloring, partition and plotting tests cover tri-only, quad-only, mixed |

Post-design re-check: unchanged (one AMEND, justified).

## Project Structure

### Documentation (this feature)

```text
specs/005-adcirc-io-coloring-partition/
├── spec.md, plan.md, research.md, data-model.md, quickstart.md
├── contracts/api.md
├── checklists/requirements.md, checklists/io-quality.md
└── tasks.md
```

### Source Code (repository root)

```text
src/chilmesh/
├── CHILmesh.py              # public write_to_2dm (+ nodestrings); reader reads NS cards
├── coloring.py              # NEW: color_elements, check_coloring
├── partition.py             # NEW: partition_report; partition_mesh only if exploration fails
├── adcirc_output_io.py      # NEW: ASCII fort.61/62/63/64, maxele/maxvel readers
├── adcirc_netcdf_io.py      # NEW: NetCDF readers, lazy import of netCDF4
├── fort22_io.py             # NEW: best track (NWS 8/20) read/write
├── owi_io.py                # NEW: fort.221/222 read/write
├── boundary_forcing_io.py   # NEW: fort.19/20/24 read/write
├── chilplotting.py          # NEW function plot_field (builds on plot_filled :298)
└── __init__.py              # re-exports
scripts/
├── check_adcirc_testsuite.py        # NEW real-file gate (not in CI)
├── adcirc_testsuite_manifest.json   # NEW: url, size, sha256, format, stage per file
└── explore_watershed_partition.py   # NEW: B1 exploration driver (reads segmentanymesh from a local clone path arg)
tests/
├── fixtures/synthetic_adcirc.py     # NEW: writers of synthetic files from the format description
├── fixtures/atcf/                   # NEW: one NOAA ATCF best-track file + SOURCE.md
└── test_2dm_nodestrings.py, test_coloring.py, test_partition_report.py,
    test_adcirc_output_ascii.py, test_adcirc_output_netcdf.py, test_plot_field.py,
    test_fort22.py, test_owi.py, test_boundary_forcing.py
```

**Structure Decision**: Single library; one module per format family, matching the existing `fort13_io.py` / `fort14_io.py` / `fort15_io.py` pattern [code src/chilmesh/].

## Stage plan

| Stage | Order inside stage | Exit gate |
|---|---|---|
| B1 | US1 .2dm → US2 coloring → US3 report → exploration → (partitioner or decision record) | pytest green; exploration table recorded in research.md; constitution VI amendment drafted |
| B2 | pin B2 real files (manifest entries; the script itself is built in the Foundational phase) → US4 ASCII → US5 NetCDF → US6 plot | pytest green; gate script passes all B2 formats |
| B3 | start fort.20 search, pin other B3 real files → US7 fort.22 + OWI → US8 fort.19/20/24 (in parallel with the search) → fort.20 gate | real fort.20 pinned; pytest green; gate script passes all B3 formats (and all earlier ones) |

## Complexity Tracking

| Item | Why needed | Simpler alternative rejected because |
|---|---|---|
| Constitution VI amendment | Feature adds nodestrings to 2DM write | Leaving text as is would make docs contradict code (VIII) |
| Conditional partitioner | Operator decision ties location to exploration result | Building both up front wastes work in one branch |
