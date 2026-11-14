<!-- provenance: author=domattioli model=claude-opus-5-5 effort=medium date=2026-10-12 skill=speckit-pipeline repo=CHILmesh session=4e33ac09-e474-4ebc-a2dd-027a9e0a0442 -->
# Tasks: ADCIRC I/O, element coloring, and mesh partitioning (#289 B1-B3)

**Input**: spec.md, plan.md, research.md, data-model.md, contracts/api.md
**Tests**: Required (constitution III, TDD for public API). Tests are written first and must fail before implementation.
**Order**: Stage B1 → B2 → B3. A stage starts only after the previous stage's exit gate task is checked.

## Phase 1: Setup

- [x] T001 Add `netcdf = ["netCDF4"]` optional extra in pyproject.toml, and add `netCDF4` to the `dev` extra so CI runs the NetCDF tests (the `netcdf` extra itself stays netCDF4 only)
- [x] T002 Add `.adcirc_testsuite/` to .gitignore
- [x] T003 [P] Create synthetic fixture writer module skeleton in tests/fixtures/synthetic_adcirc.py

## Phase 2: Foundational

- [x] T004 Create scripts/adcirc_testsuite_manifest.json schema (name, format, stage, url, lfs, sha256, size) with empty entries list
- [x] T005 Implement download, LFS-pointer rejection, sha256 check, format→reader dispatch table, `--stage` and `--cache` options in scripts/check_adcirc_testsuite.py (never executes downloaded files; readers not yet registered are reported as FAIL)

## Stage B1

### US1 .2dm writer with nodestrings (P1)

- [x] T006 [US1] Write failing tests for public write, NS round-trip (tri-only, quad-only, mixed, no boundary) and `_write_2dm` alias in tests/test_2dm_nodestrings.py
- [x] T007 [US1] Add public `write_to_2dm(filename, *, nodestrings=True)` emitting NS cards, keep `_write_2dm` alias, in src/chilmesh/CHILmesh.py
- [x] T008 [US1] Parse NS cards into boundary segments in `read_from_2dm` in src/chilmesh/CHILmesh.py
- [x] T009 [US1] Re-export and docstring example; update CHANGELOG.md

### US2 Element coloring (P1)

- [x] T010 [P] [US2] Write failing tests (no node-sharing pair same color, deterministic, bound max_degree+1, single element, tri/quad/mixed) in tests/test_coloring.py
- [x] T011 [US2] Implement `color_elements` and `check_coloring` using Vert2Elem in src/chilmesh/coloring.py
- [x] T012 [US2] Re-export in src/chilmesh/__init__.py

### US3 Partition report and partitioner decision (P2)

- [x] T013 [P] [US3] Write failing tests for `partition_report` (hand-checked small mesh: sizes, imbalance, cut edges, pieces; disconnected part; k bounds) in tests/test_partition_report.py
- [x] T014 [US3] Implement `partition_report` in src/chilmesh/partition.py and re-export
- [x] T015 [US3] Write exploration driver taking a local segmentanymesh clone path argument in scripts/explore_watershed_partition.py (3 connected meshes, k=4 and 8, farthest-point seeds; reads segmentanymesh only)
- [x] T016 [US3] Run exploration; record table (mesh, k, imbalance, pieces) and PASS/FAIL against <=1.25 and one-piece rule in specs/005-adcirc-io-coloring-partition/research.md R3
- [x] T017 [US3] If T016 FAIL: write failing tests for `partition_mesh` spectral (including the SC-004 threshold, <=1.25 and one piece, at k=4 and 8 on the three FR-006 meshes), k-range refusal, and the metis-missing message in tests/test_partition_mesh.py, then implement in src/chilmesh/partition.py and add `partition = ["pymetis"]` extra in pyproject.toml; if spectral also misses the threshold, STOP and report to operator. If T016 PASS: record decision "partitioner in SegmentAnyMesh" in research.md and skip
- [x] T017a Find a real SMS `.2dm` from any public source (adcirc-testsuite first, then other public repositories); pin it by SHA-256 at a fixed-commit URL in scripts/adcirc_testsuite_manifest.json (stage B1) with source and license noted, and register `read_from_2dm` in the gate script; if none is found, STOP and report to operator (B1 cannot close)

- [x] T018 **B1 EXIT GATE**: full pytest green; SC-003 checked on all bundled meshes; T016 result recorded; real `.2dm` pinned and `python scripts/check_adcirc_testsuite.py --stage B1` exits 0; constitution VI 2DM sentence updated in .specify/memory/constitution.md (nodestrings now written)

## Stage B2

- [x] T019 Pin one real file per B2 format (fort.61, 62, 63, 64, maxele.63, maxvel.63, fort.63.nc, fort.64.nc, maxele.63.nc, maxvel.63.nc, fort.61.nc, fort.62.nc) from adcirc/adcirc-testsuite in scripts/adcirc_testsuite_manifest.json (text via raw.githubusercontent.com, LFS via media.githubusercontent.com)

### US4 ASCII output readers (P1)

- [x] T020 [P] [US4] Add synthetic writers for full and sparse scalar/vector global and station files, plus truncated variant, in tests/fixtures/synthetic_adcirc.py
- [x] T021 [US4] Write failing tests (values, sparse default, -99999→NaN, node-count mismatch message, single step, truncation warning with dropped count, maxele/maxvel two-set file returns time-of-maximum as second record) in tests/test_adcirc_output_ascii.py
- [x] T022 [US4] Implement `read_adcirc_output` with offset index and single-step seek in src/chilmesh/adcirc_output_io.py; re-export; register in gate script
- [x] T023 [US4] Measure SC-005 single-step memory on a synthetic 1M-node file (peak `tracemalloc` during the call vs 3.0 x returned array `nbytes`); record numbers and PASS/FAIL in research.md

### US5 NetCDF readers (P2)

- [x] T024 [US5] Write failing tests (same meaning as ASCII, single step lazy, ImportError names extra, `import chilmesh` works without netCDF4) in tests/test_adcirc_output_netcdf.py, skipping data tests when netCDF4 absent
- [x] T025 [US5] Implement `read_adcirc_output_nc` with function-level import in src/chilmesh/adcirc_netcdf_io.py; re-export; register in gate script

### US6 Plot field on mesh (P2)

- [x] T026 [P] [US6] Write failing headless tests (colorbar present, NaN blank, length mismatch message, tri/quad/mixed) in tests/test_plot_field.py
- [x] T027 [US6] Implement `plot_field` on top of plot_filled in src/chilmesh/chilplotting.py

- [x] T028 **B2 EXIT GATE**: full pytest green; `python scripts/check_adcirc_testsuite.py --stage B2` exits 0 with every B2 format pinned

## Stage B3

### Real fort.20 (blocking for B3 close)

- [x] T029 Locate a real fort.20 (redistribution-free source or a stable URL the script fetches at run time); record source, licence status and URL in research.md R9
- [x] T030 Pin the fort.20 and one real file per other B3 format (fort.19, fort.22 NWS 8, fort.22 NWS 20 if present, fort.221, fort.222, fort.24) in scripts/adcirc_testsuite_manifest.json
- [x] T031 **fort.20 GATE** (checked immediately before T039; US7 and US8 do not wait on it): manifest has a fort.20 entry with sha256; if no source found, STOP and report to operator (B3 cannot close)

### US7 fort.22 best track and OWI (P1)

- [x] T032 [P] [US7] Add one NOAA ATCF best-track file and tests/fixtures/atcf/SOURCE.md (URL, retrieval date, US government work)
- [x] T033 [US7] Write failing tests (typed fields, untyped columns kept verbatim, blank fields as missing, round-trip for NWS 8 and NWS 20, unsupported NWS error) in tests/test_fort22.py
- [x] T034 [US7] Implement best-track read/write for NWS 8 and 20 in src/chilmesh/fort22_io.py; re-export; register in gate script
- [x] T035 [P] [US7] Write failing tests with synthetic OWI 221/222 (grid header, U/V blocks, round-trip) in tests/test_owi.py
- [x] T036 [US7] Implement `read_owi` / `write_owi` in src/chilmesh/owi_io.py; re-export; register in gate script

### US8 fort.19/20/24 (P2)

- [x] T037 [US8] Write failing round-trip and mesh-count tests for fort.19, fort.20, fort.24 with synthetic files in tests/test_boundary_forcing.py; confirm the fort.20 IBTYPE set in FR-017 against the ADCIRC fort.14 documentation and record it in research.md R8
- [x] T038 [US8] Implement read/write for fort.19, fort.20, fort.24 in src/chilmesh/boundary_forcing_io.py; re-export; register in gate script

- [x] T039 **B3 EXIT GATE**: full pytest green offline; `python scripts/check_adcirc_testsuite.py --stage all` exits 0; no adcirc-testsuite file in git (SC-007)

## Final Phase: Polish

- [x] T040 [P] Update README format table and docs for new readers/writers and extras; every new public function has a docstring example linked to a test (constitution VIII)
- [x] T041 Update CHANGELOG.md with all new public names
- [x] T042 Move specs to DomI: `bash ~/Projects/DomI/skills/speckit-pipeline/scripts/move_specs_to_domi.sh --apply` (after code is committed; runs last)

## Dependencies

Setup → Foundational → B1 (US1, US2, US3 in any order; T016 before T017; T017a any time in B1) → T018 → B2 (T019, then US4 → US5; US6 after US4) → T028 → B3 (T029 and T030 start the stage; US7 and US8 run while the fort.20 search continues; T031 before T039) → T039 → Polish → T042.

## Parallel examples

- B1: T006, T010, T013 together (different test files).
- B2: T020 and T026 together.
- B3: T032 and T035 together; T029 (fort.20 search) alongside US7/US8.

## Strategy

MVP = B1 US1 + US2. Each stage ships as its own increment after its gate.
