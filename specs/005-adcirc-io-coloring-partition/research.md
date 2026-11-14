<!-- provenance: author=domattioli model=claude-opus-5-5 effort=medium date=2026-10-12 skill=speckit-pipeline repo=CHILmesh session=4e33ac09-e474-4ebc-a2dd-027a9e0a0442 -->
# Research: #289 B1-B3

## R1 .2dm nodestrings
- Decision: write one `NS` card sequence per boundary segment; the last node ID of a nodestring is negative (SMS convention); reader parses NS cards into boundary segments.
- Rationale: SMS 2DM format convention; current writer omits them [code src/chilmesh/CHILmesh.py:2759].
- Alternatives: separate sidecar file (rejected, not readable by SMS).

## R2 Coloring
- Decision: greedy coloring over the element node-sharing graph (built from Vert2Elem [code src/chilmesh/CHILmesh.py:685]), elements visited in a fixed order (largest degree first, ties by element ID) so results are deterministic. Bound: max node-sharing degree + 1.
- Rationale: operator chose vertex-sharing rule (Clarifications Q4); greedy meets the bound by construction.
- Alternatives: DSatur (fewer colors, slower; may be added later as option).

## R3 Partitioner decision path
- Decision: run `select_labels` (equidistant multi-source watershed) [code domattioli/segmentanymesh prototype/click_select.py:186] with k seeds chosen by farthest-point sampling on element centroids, on 3 connected meshes (bundled `donut` ~276, `structured` ~660, `block_o` ~5,200 elements [code src/chilmesh/examples.py:43-59]; `annulus` substitutes for any found disconnected), k=4 and 8. Pass: every run imbalance <= 1.25 and every part one piece.
- If pass: record decision; partitioner work moves to SegmentAnyMesh (separate task, out of this feature).
- If fail: CHILmesh `partition_mesh(method="spectral")` using scipy recursive spectral bisection of the element dual graph; `method="metis"` via optional pymetis.
- Result table (scripts/explore_watershed_partition.py, T016; imbalance over the k labeled parts; max pieces over the k parts; unassigned = contested/unreached elements left by `select_labels`, which are in no part):

```
mesh          k  imbalance  max_pieces  unassigned
donut         4      1.038           1          10
donut         8      1.184           1          26
structured    4      1.333           1           0
structured    8      1.679           1          36
block_o       4      1.181           1          37
block_o       8      1.260           1          71
```
- Result: **FAIL** (pass needs every imbalance <= 1.25 and one piece per part). structured k=4 (1.333), structured k=8 (1.679) and block_o k=8 (1.260) exceed 1.25. Every labeled part is one piece, but `select_labels` also leaves 10-71 contested elements unassigned in 5 of 6 runs. 
- T017 spectral result (`partition_mesh(method="spectral")`, measured):

```
mesh          k  imbalance  max_pieces
donut         4      1.000           1
donut         8      1.014           1
structured    4      1.000           1
structured    8      1.006           1
block_o       4      1.000           1
block_o       8      1.000           1
```
- Decision: **partitioner in CHILmesh (spectral)**. All 6 runs meet imbalance <= 1.25 with one piece per part and every element assigned. METIS stays optional (`partition` extra).

## R4 ASCII output formats
- Decision: one parser for the shared layout: header line, then `NDSETS NP DT NSPOOL IRTYPE` line, then per set a `TIME IT` line (or `TIME IT NNONDEF DEFVAL` for sparse) and node records. IRTYPE 1 = scalar, 2 = vector. Station files have the same layout with NSTA in place of NP. maxele/maxvel: same layout, usually 1 or 2 sets (second set is time of max).
- Single-step read: scan record offsets, then seek; memory one step.
- Truncation: drop incomplete last set, warn with count dropped (Clarifications Q5).
- Fill value -99999 → NaN.

## R5 NetCDF
- Decision: `netcdf = ["netCDF4"]` only (Clarifications Q3); import inside functions; variables `zeta`, `u-vel`/`v-vel`, `zeta_max`, `time`; slice by index without loading the whole variable.

## R6 fort.22 best track
- Decision: ATCF b-deck comma-separated records; NWS 8 and NWS 20 share the columns, NWS 20 adds extra columns used by the generalized asymmetric vortex. Typed fields per FR-015; all other columns kept as untyped text in record order (`extra`) and written back unchanged, so the version-dependent column set is a non-issue for round-trip. Writer re-emits fixed widths so read-write-read is equal.
- Fixture: one NOAA ATCF file (US government work), source URL and date in `tests/fixtures/atcf/SOURCE.md`.

## R7 OWI fort.221/222
- Decision: header line with start and end time; per snapshot a line with iLat, iLong, DX, DY, SWLat, SWLon, date; values in 8 per line fixed width; 222 has U then V blocks.

## R8 fort.19/20/24
- fort.19: first line ETIMINC, then per time step one value per elevation-specified boundary node.
- fort.20: first line FTIMINC, then per time step one value per flux-specified boundary node (non-periodic normal-flow IBTYPE set per FR-017; confirmed in T037).
- T037 IBTYPE confirmation: adcirc/adcirc docs (main at e8b62a70db39dbc8875785ee999a2214be75302b), `docs/technical_reference/input_files/fort20.rst`, states fort.20 is read for "specified non-zero normal flow" boundaries `IBTYPE=2, 12, or 22` with NFFR=0. `docs/user_guide/model_configuration/boundary_conditions/flux_specified_boundaries.rst` lists 2, 12, 22 as external nonzero-flow types (102, 112, 122 are their baroclinic variants and are not named for fort.20); IBTYPE 32 and 52 do not appear in either document. Decision: matching set is {2, 12, 22} (`FORT20_IBTYPES`), narrower than the FR-017 text, which named 32 and "52 if listed"; neither is listed. FR-017 wording should be updated to match.
- fort.19 (`fort19.rst`): ETIMINC, then NETA values per step, one per line. fort.20: FTIMINC, then NFLBN values per step, one per line (real files carry trailing `!` comments; the reader takes the first token).
- fort.24 (`fort24.rst`): per constituent four header lines (alpha line, frequency, the literal `1`, constituent name; ADCIRC skips them) then `JN amp phase` per node. The data-model fields `nodal_factor` / `eq_arg` do not exist in the file: the third line is stored as `nodal_factor` (value 1) and there is no equilibrium argument. Real file `global-tide-2d/fort.24` has one constituent (M2) and 27330 nodes.
- fort.24: per constituent: alpha line, frequency, nodal factor and equilibrium argument line, name, then one `node amp phase` line per node.

## R9 Real-file gate
- Decision: manifest JSON lists one or more real files per format (url, sha256, size, format, stage). Script downloads to gitignored `.adcirc_testsuite/`, uses raw.githubusercontent.com for text and media.githubusercontent.com for LFS objects, rejects a file whose head is an LFS pointer, verifies sha256, parses with the matching reader, exits non-zero on any failure. Never executes files.
- T017a outcome (B1): no `.2dm` found in adcirc-testsuite; pinned `estuary.2dm` from pwcazenave/PyFVCOM at commit f98e11234d8f92bc45c497b18477978d482ab9d8 (SMS-written, MIT, 79700 bytes, sha256 82ada67c...a0996): 1570 E3T, 845 ND, 3 NS cards (one nodestring of 25 nodes), plus SMS cards BEG, GM, BCP, MAT, BED, BEF, END, NUM, SI, DY, TD, TU that `read_from_2dm` ignores. `read_from_2dm` registered in the gate script; `--stage B1` exits 0.
- T029 outcome (B3 fort.20): adcirc-testsuite has no fort.20. Pinned instead `input/meshes/LA_v12g-WithUpperAtch/20170527.fort.20` from StormSurgeLive/asgs at commit c6376f94e3d02e824883fca88937f736bffb6cbb, licence GPL-3.0 (repository licence), URL https://raw.githubusercontent.com/StormSurgeLive/asgs/c6376f94e3d02e824883fca88937f736bffb6cbb/input/meshes/LA_v12g-WithUpperAtch/20170527.fort.20, 2494 bytes, sha256 01f8a8e3c51af63708dd10d77c18159a84f9b20c989672771d10e215bf46ba5a, not Git LFS. Content: `31536000  ! number of seconds in a year` then 84 value lines with `!` comments. Fetched at run time by the gate script, never vendored. Verified by download in B3.
- T030 B3 pins (adcirc-testsuite at aeba15e unless noted): fort.19 `adcirc_ideal_channel-woffset-2d-parallel/fort.19` (490 B, raw); fort.22 NWS 8 `adcirc_irene_holland_legacy/fort.22` (fort.15 NWS=8, 17716 B, raw); fort.22 NWS 20 `adcirc_nws30_wlcorrection/NWS_20_fort.22` (8034 B, raw); fort.221/222 `adcirc_katrina-2d/` (5.7 MB and 11.3 MB, Git LFS, media URL); fort.24 `adcirc_global-tide-2d/fort.24` (972808 B, raw). Hashes in the manifest. Vendored exception (FR-020): `tests/fixtures/atcf/bal092011.dat`, NOAA ATCF b-deck for Irene 2011, 17617 B, see its SOURCE.md.
- Coverage: the testsuite has fort.24, fort.221/222, fort.61/62/64, maxvel and NetCDF variants, but no fort.20 [web adcirc/adcirc-testsuite tree listing, coordinator]. fort.20 source search is a B3 task (candidates: ADCIRC example problems page, published case-study data, NOAA/USGS public archives) [inference].

## R10 B2 real files and SC-005 measurement
- T019 pins (adcirc/adcirc-testsuite at commit aeba15eb5e665f2eaf21defe7ecf593a276ccc46, all Git LFS, media.githubusercontent.com): ASCII fort.61/62/63/64, maxele.63, maxvel.63 from `adcirc/adcirc_quarterannular-2d/control/`; NetCDF fort.61/62/63/64.nc, maxele.63.nc, maxvel.63.nc from `adcirc/adcirc_quarterannular-2d-netcdf/control/`. All 12 formats have a real file. Quarter annular grid: 63 nodes, 3 stations, 824 time sets; maxele/maxvel hold 2 sets. sha256 and size are in the manifest. ASCII and NetCDF readings of the same quantity agree to 6e-11 with identical NaN masks [measured, coordinator run].
- Real-file layout facts found: line 2 carries `NDSETS NP DT NSPOOL IRTYPE FileFmtVersion: n` (trailing text ignored); station files use the global layout with NSTA in place of NP; the second maxele/maxvel set repeats the end time in its header.
- NetCDF files store no iteration count, so `OutputSeries.iterations` is -1 for them.
- T023 SC-005 (synthetic 1M-node file, `read_adcirc_output(path, step=last)`, peak `tracemalloc` during the call vs returned `values.nbytes`; macOS, Python 3.14, numpy 2.5.3):

  | case | returned nbytes | peak bytes | ratio | limit | result |
  |---|---|---|---|---|---|
  | scalar, 1 step file | 8,000,000 | 16,133,065 | 2.017 | 3.0 | PASS |
  | scalar, 3 step file | 8,000,000 | 16,132,859 | 2.017 | 3.0 | PASS |
  | vector, 1 step file | 16,000,000 | 32,132,489 | 2.008 | 3.0 | PASS |

  Overall SC-005: PASS. The 2x comes from the per-set work array plus the stacked return copy; the text is parsed in 20,000-line slabs, so it adds no full-file buffer. Enforced in `tests/test_adcirc_output_ascii.py::test_single_step_memory_1m_nodes`.

## Open risks
- No real fort.20 may exist with clear redistribution or fetchable URL → B3 blocked (by operator design).
- NWS 20 column set differs between ADCIRC versions [inference].
- 1M-node single-step memory target depends on text seek efficiency. Measure: peak `tracemalloc` during the call vs 3.0 x returned array `nbytes` (SC-005).
