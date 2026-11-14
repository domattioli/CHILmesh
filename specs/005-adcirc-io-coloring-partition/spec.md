<!-- provenance: author=domattioli model=claude-opus-5-5 effort=medium date=2026-10-12 skill=speckit-pipeline repo=CHILmesh session=4e33ac09-e474-4ebc-a2dd-027a9e0a0442 -->
# Feature Specification: ADCIRC I/O, element coloring, and mesh partitioning (CHILmesh #289 B1-B3)

**Feature Branch**: `development` (CHILmesh single-branch policy; no feature branch)
**Created**: 2026-10-12
**Status**: Draft
**Input**: CHILmesh issue #289 "Catalogue: ADCIRC file I/O still missing from CHILmesh" and its binding comment "Decisions (operator grill, 2026-10-12)", scope B1-B3. [issue]

## Overview and stage order

This specification covers three ordered, iterative stages. Each stage must be complete and pass its own gates before the next stage starts. A later stage may reuse the earlier stages' results, but no earlier stage depends on a later one.

| Stage | Contents | Gate to leave the stage |
|---|---|---|
| B1 | Public `.2dm` writer with nodestrings; element coloring; mesh partitioner prototype with a partition quality report | B1 user stories pass; partitioner decision recorded (SegmentAnyMesh or CHILmesh) |
| B2 | ASCII output readers (fort.61, fort.62, fort.63, fort.64, maxele.63, maxvel.63); NetCDF output readers as an optional extra; plotting of output fields on the mesh | B2 user stories pass; real-file gate passes for every B2 format |
| B3 | Locate a real fort.20; fort.22 meteorological forcing (best track NWS 8/20; OWI fort.221/222 grids); fort.19, fort.20, fort.24 boundary forcing, read and write | Real fort.20 located and pinned; B3 user stories pass; real-file gate passes for every B3 format |

Out of scope (B4, deferred by the operator): structured editing of fort.15 and hot start files fort.67/fort.68. Also out of scope: running ADCIRC or `adcprep`, or wrapping Fortran. [issue]

## Clarifications

### Session 2026-10-12

- Q: Real-file gate: only pinned files gate, or a real file for every format? → A: A real file for every format. The script pins one real file per format. adcirc-testsuite has real fort.24, OWI fort.221/222, fort.61/62/64, maxvel and NetCDF variants (some under Git LFS); only fort.20 has none, so locating a real fort.20 (any redistribution-free or fetch-at-runtime source) is an explicit B3 task with its own exit gate. [issue] [web adcirc/adcirc-testsuite tree listing]
- Q: Balance threshold for the partitioner decision? → A: Imbalance ratio (max part size / mean part size) at most 1.25 at k=4 and k=8 on 3 meshes, AND every part is one connected piece.
- Q: NetCDF optional extra? → A: `netcdf = ["netCDF4"]` only; xarray not used.
- Q: Coloring adjacency rule? → A: Vertex-sharing only: elements of one color share no node.
- Q: Truncated ADCIRC output? → A: Return the complete records and warn, stating how many records were dropped.

### Operator follow-up 2026-10-12

- Q: `.2dm` real-file gate when adcirc-testsuite has no `.2dm`? → A: A real SMS `.2dm` MUST be found from any public source (not only adcirc-testsuite), pinned by SHA-256 at a fixed commit URL, and parsed by `scripts/check_adcirc_testsuite.py` before Stage B1 closes, the same rule as fort.20. The file is fetched at run time and not vendored. If none is found, B1 cannot close and the operator decides.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Write a mesh to SMS `.2dm` with boundaries (Stage B1, Priority: P1)

A modeler who built or edited a mesh in CHILmesh saves it as `.2dm` so that SMS and other tools can open it with its boundaries intact.

**Why this priority**: A private writer already exists and is tested, so this is the smallest step that closes a listed gap. [code src/chilmesh/CHILmesh.py:2755] [code tests/test_2dm_roundtrip.py]

**Independent Test**: Write a mesh with boundary segments to `.2dm`, read it back, and compare nodes, elements, and nodestrings.

**Acceptance Scenarios**:

1. **Given** a mixed triangle and quadrilateral mesh, **When** the user calls the public `.2dm` write method, **Then** reading the file back gives the same node coordinates (within the written precision), the same elements, and the same element types.
2. **Given** a mesh with boundary segments, **When** it is written to `.2dm`, **Then** each boundary segment appears as one nodestring with the same node order, and reading it back restores those nodestrings.
3. **Given** existing code that calls the private writer, **When** the public method is added, **Then** that code keeps working (the private name stays as an alias).

---

### User Story 2 - Color elements so no two neighbors share a color (Stage B1, Priority: P1)

A developer who runs per-element work in parallel (for example smoothing or assembly) asks CHILmesh for an element coloring, so that elements of one color can be processed at the same time without conflicts.

**Why this priority**: Coloring stays in CHILmesh in every partitioner outcome, and it reuses the existing vertex-to-element adjacency. [issue] [code src/chilmesh/CHILmesh.py:92,685]

**Independent Test**: Color several bundled meshes and check that no two elements sharing a node share a color.

**Acceptance Scenarios**:

1. **Given** any valid mesh, **When** coloring is requested, **Then** every element gets exactly one color and no two elements of the same color share a node.
2. **Given** the same mesh and the same options, **When** coloring is run twice, **Then** the result is identical.
3. **Given** a mesh, **When** coloring is requested, **Then** the number of colors is reported and is at most the largest number of elements that touch any one element through shared nodes, plus one.

---

### User Story 3 - Split a mesh into balanced subdomains and see how good the split is (Stage B1, Priority: P2)

A modeler asks for a mesh to be split into k subdomains and receives a label per element plus a quality report (size balance, cut edges, connectivity of each part).

**Why this priority**: This is a prototype. The operator decision requires an exploration first: test whether SegmentAnyMesh's seeded watershed over element adjacency, run with k seeds, gives balanced subdomains. If it does, the partitioner lives in SegmentAnyMesh (which depends on CHILmesh). If it does not, CHILmesh gets a SciPy-based partitioner with an optional pymetis extra. The quality report stays in CHILmesh in both cases. [issue] [code domattioli/segmentanymesh prototype/click_select.py:186 `select_labels`, "Equidistant multi-source watershed labeling"]

**Independent Test**: Run the quality report on any element labeling (including a hand-made one) and check its numbers by hand on a small mesh. Run the exploration on at least three meshes of different size and shape.

**Acceptance Scenarios**:

1. **Given** a mesh and an element labeling, **When** the quality report runs, **Then** it reports per-part element counts, the imbalance ratio (largest part over mean part size), the number of cut edges, and the number of connected pieces per part.
2. **Given** the exploration result, **When** the watershed approach meets the balance threshold on the exploration meshes, **Then** the decision "partitioner in SegmentAnyMesh" is recorded and CHILmesh ships only the report; **otherwise** CHILmesh ships its own partitioner.
3. **Given** CHILmesh ships its own partitioner, **When** pymetis is not installed, **Then** the SciPy-based method works, and asking for the pymetis method gives a clear message that names the extra to install.

---

### User Story 4 - Read ADCIRC ASCII output aligned to mesh nodes (Stage B2, Priority: P1)

A modeler who ran ADCIRC loads fort.63, fort.64, maxele.63, maxvel.63 (global, node-aligned) and fort.61, fort.62 (station) results into arrays, ready for analysis and plotting with the matching mesh.

**Why this priority**: Reading results is the most common post-run task. [issue]

**Independent Test**: Read synthetic files written from the format description and compare to the known values; then read the real test-suite files through the real-file gate.

**Acceptance Scenarios**:

1. **Given** a global output file and its mesh, **When** the user reads it, **Then** the result holds the time values and one value per node per time step (two components for velocity), in mesh node order.
2. **Given** a sparse-format record (only non-default nodes listed), **When** it is read, **Then** missing nodes take the stated default (dry) value.
3. **Given** a station file, **When** it is read, **Then** the result holds the time values and one value per station per time step.
4. **Given** an output file whose node count does not match the mesh, **When** the user reads it with that mesh, **Then** the reader raises an error that names both counts.

---

### User Story 5 - Read ADCIRC NetCDF output when the optional extra is installed (Stage B2, Priority: P2)

A modeler with fort.63.nc, fort.64.nc or maxele.63.nc reads them the same way as the ASCII files, and can read one time step without loading the whole file.

**Why this priority**: Most modern runs write NetCDF, but the dependency must stay optional. [issue]

**Independent Test**: Read a small synthetic NetCDF file; check that a single time slice can be read alone; check the message when the extra is missing.

**Acceptance Scenarios**:

1. **Given** the NetCDF extra is installed, **When** the user reads a NetCDF output file, **Then** the result has the same shape and meaning as the ASCII reader's result for the same data.
2. **Given** a large file, **When** the user asks for one time step, **Then** only that step is loaded.
3. **Given** the NetCDF extra is not installed, **When** the user calls the NetCDF reader, **Then** a clear message names the extra to install, and importing CHILmesh itself still works.

---

### User Story 6 - Plot an output field on the mesh (Stage B2, Priority: P2)

A modeler draws a maximum-elevation map or a single time step of elevation on the mesh with one call, with a color bar and with dry nodes shown as missing.

**Why this priority**: It turns the readers into something the user can see; the existing filled-plot code is the base. [code src/chilmesh/chilplotting.py:298]

**Independent Test**: Plot a synthetic node field on a bundled mesh without a display, and check that the figure contains the field and a color bar.

**Acceptance Scenarios**:

1. **Given** a mesh and a node-aligned field from a B2 reader, **When** the plot entry point is called, **Then** a figure is produced with a color bar and the elements colored by the field.
2. **Given** a field with dry (missing) values, **When** it is plotted, **Then** those areas are left blank, not drawn as a number.
3. **Given** a field whose length does not equal the node count, **When** it is plotted, **Then** an error names both lengths.

---

### User Story 7 - Read and write fort.22 meteorological forcing (Stage B3, Priority: P1)

A modeler reads a best-track fort.22 (NWS 8 and NWS 20) or OWI wind and pressure grids (fort.221 and fort.222), checks or edits the values, and writes them back.

**Why this priority**: Storm forcing is the most used forcing input in coastal studies. [issue]

**Independent Test**: Read the NOAA ATCF best-track fixture and synthetic OWI files; write them; read them back; compare.

**Acceptance Scenarios**:

1. **Given** a best-track fort.22, **When** it is read, **Then** each record gives time, position, central pressure, maximum wind and wind radii as typed fields.
2. **Given** a parsed best-track file, **When** it is written without changes, **Then** reading the output again gives equal records.
3. **Given** OWI fort.221 and fort.222 files, **When** they are read, **Then** each time snapshot gives its grid definition and the pressure or wind component arrays.
4. **Given** a fort.22 with an NWS value this feature does not support, **When** it is read, **Then** an error names the unsupported format.

---

### User Story 8 - Read and write boundary forcing fort.19, fort.20, fort.24 (Stage B3, Priority: P2)

A modeler reads and writes non-periodic elevation (fort.19), non-periodic flux (fort.20) and self-attraction and load tide (fort.24) files that match their mesh's boundaries.

**Why this priority**: These inputs complete the set of boundary forcing a user touches around a run. [issue]

**Independent Test**: Round-trip synthetic files for each type; read the real fort.19 test-suite file through the real-file gate.

**Acceptance Scenarios**:

1. **Given** a fort.19 file, **When** it is read, **Then** the result holds the time increment and one value per open-boundary node per time step.
2. **Given** a fort.20 file, **When** it is read, **Then** the result holds one value per flux-boundary node per time step.
3. **Given** a fort.24 file, **When** it is read, **Then** the result holds, per constituent, its name, frequency and nodal factors, and amplitude and phase per node.
4. **Given** any of these files read and written without change, **When** read again, **Then** the values are equal.

### Edge Cases

- Mesh with no boundary segments: `.2dm` writing produces no nodestrings and does not fail.
- Mesh with a single element: coloring gives one color; partitioning into k > 1 parts is refused with a clear message.
- k larger than the number of elements: refused with a clear message.
- Disconnected mesh: the partition report counts pieces per part correctly; coloring still works. A mesh with more disconnected pieces than k cannot meet the one-piece rule; the exploration meshes are connected.
- ADCIRC output file cut short (fewer records than the header says): the reader returns only the complete records and warns with the number of records dropped; a partly written last record is dropped, never returned.
- Output written with a fill value for dry nodes (-99999): treated as missing.
- fort.22 best-track lines with blank or short fields: parsed as missing, not as zero.
- Very large output files: global readers must not need more than one time step in memory when the user asks for a single step.
- No real fort.20 found: B3 cannot close; the search result is recorded and the operator decides (FR-023a).

## Requirements *(mandatory)*

### Functional Requirements

**Stage B1**

- **FR-001**: CHILmesh MUST provide a public method to write a mesh to `.2dm`, including nodes, triangles, quadrilaterals, and one nodestring per boundary segment.
- **FR-002**: The `.2dm` reader MUST read nodestrings back so that write then read keeps boundary segments.
- **FR-003**: The existing private `.2dm` writer name MUST keep working.
- **FR-004**: CHILmesh MUST provide element coloring in which no two elements that share a node have the same color, and whose result is the same on every run with the same input.
- **FR-005**: CHILmesh MUST provide a partition quality report for any element labeling: per-part size, imbalance ratio, cut-edge count, and connected pieces per part.
- **FR-006**: An exploration MUST measure SegmentAnyMesh's seeded watershed with k seeds on at least three connected meshes (the bundled `donut` (~276 elements), `structured` (~660) and `block_o` (~5,200) examples; the driver checks each is one connected piece and substitutes `annulus` for any that is not), at k=4 and k=8, and record the imbalance ratio and pieces per part for each mesh and k, before any partitioner is written. The approach passes only if every run has imbalance ratio at most 1.25 and every part is one connected piece.
- **FR-007**: If the exploration does not meet the balance threshold, CHILmesh MUST provide a partitioner that works with its required dependencies alone, and MAY use pymetis when that optional extra is installed.
- **FR-008**: If the exploration meets the threshold, the partitioner MUST NOT be added to CHILmesh, and SegmentAnyMesh is not changed by this feature (work there is a separate task).

**Stage B2**

- **FR-009**: CHILmesh MUST read ASCII fort.61, fort.62, fort.63, fort.64, maxele.63 and maxvel.63, in full and sparse record forms, and return time values plus node-aligned (or station-aligned) arrays. For maxele.63 and maxvel.63 every data set present is returned; when a second set exists it is the time of the maximum per node and is returned as the second record.
- **FR-010**: Global readers MUST check the node count against a given mesh and raise an error that names both counts when they differ.
- **FR-011**: Global readers MUST allow reading a single time step without reading all steps into memory.
- **FR-011a**: When an ASCII output file ends early, readers MUST return only complete records and emit a warning that states how many records were dropped.
- **FR-012**: CHILmesh MUST read NetCDF fort.63.nc, fort.64.nc, maxele.63.nc, maxvel.63.nc, fort.61.nc and fort.62.nc when the optional extra `netcdf` (netCDF4 only) is installed, with single-step reading, and with results of the same meaning as the ASCII readers. No other NetCDF file is in scope.
- **FR-013**: Importing CHILmesh MUST NOT need the NetCDF extra; calling a NetCDF reader without it MUST give a message naming the extra.
- **FR-014**: The plotting module MUST provide an entry point that draws a node-aligned field on the mesh with a color bar and blank dry areas.

**Stage B3**

- **FR-015**: CHILmesh MUST read and write best-track fort.22 for NWS 8 and NWS 20. Typed fields are the ATCF b-deck columns for time, latitude, longitude, maximum wind, central pressure, isotach and its four radii, and radius of maximum wind; every other column (including the NWS 20 extra columns) is kept as untyped text in record order and written back unchanged, so the writer needs no per-version column list.
- **FR-016**: CHILmesh MUST read and write OWI fort.221 (pressure) and fort.222 (wind) grid files.
- **FR-017**: CHILmesh MUST read and write fort.19, fort.20 and fort.24. When a mesh is given, the per-step value count MUST equal the matching boundary node count and a mismatch raises an error naming both counts. Matching nodes: fort.19 uses the fort.14 elevation-specified (open) boundary nodes, in file order; fort.20 uses the fort.14 flux-specified boundary nodes with IBTYPE values 2, 12 and 22, the set the ADCIRC fort.20 documentation lists [web: adcirc/adcirc e8b62a7 fort20.rst, flux_specified_boundaries.rst]; fort.24 uses all mesh nodes.
- **FR-018**: Every B3 writer MUST round-trip: read, write without change, read again gives equal values.

**Fixtures and real-file gate (all stages)**

- **FR-019**: Repository tests MUST use synthetic files written from the format descriptions. No third-party ADCIRC file is copied into the repository or into Valence. [issue]
- **FR-020**: The one exception is a single small best-track fort.22 fixture (one storm, under 100 KB) taken directly from NOAA's ATCF archive, with its source URL and retrieval date recorded next to the file. [issue]
- **FR-021**: `scripts/check_adcirc_testsuite.py` MUST pin at least one real file per supported format (fort.61, 62, 63, 64, maxele.63, maxvel.63, each NetCDF variant read, fort.19, fort.20, fort.22 best track, fort.221, fort.222, fort.24, .2dm). It MUST download the needed `adcirc/adcirc-testsuite` files (and, for `.2dm`, a real SMS file from any public source) into a gitignored folder, check each against a pinned SHA-256, and parse each with the matching reader. Small text files come from raw.githubusercontent.com; Git LFS files come from media.githubusercontent.com. [issue] [inference: LFS pointers are not the file]
- **FR-022**: The real-file script MUST NOT run in CI. It is a success gate for the pipeline: any pinned real file that fails to download-check or parse blocks completion of the stage that owns its file type. [issue]
- **FR-023**: Downloaded files are data only; the script MUST never execute them.
- **FR-023a**: Before B3 closes, a real fort.20 file MUST be located (any source that is redistribution-free or that the script can fetch at run time), pinned in the manifest, and parsed by the gate. B3 cannot close without it. If no such source is found, work on B3 closure stops and the question goes to the operator; the fort.22, fort.221/222, fort.19 and fort.24 work may still be implemented and tested meanwhile.
- **FR-024**: No redistribution request goes to the ADCIRC maintainers within this feature. [issue]

### Key Entities

- **Mesh**: existing CHILmesh object: nodes, elements, boundary segments, adjacencies.
- **Nodestring**: ordered list of node IDs for one boundary segment in `.2dm`.
- **Element coloring**: one color index per element, with the color count.
- **Partition labeling**: one part index per element, for k parts.
- **Partition quality report**: per-part sizes, imbalance ratio, cut edges, pieces per part.
- **Output time series**: time values plus per-node or per-station values per step (one or two components); dry values marked missing.
- **Storm track**: ordered best-track records (time, position, pressure, wind, radii).
- **Gridded met snapshot**: grid definition plus pressure or wind arrays at one time.
- **Boundary forcing series**: per-boundary-node values per time step (fort.19, fort.20) or per-constituent amplitude and phase per node (fort.24).
- **Real-file manifest**: list of test-suite files with URL, size and pinned SHA-256.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: All 8 user stories pass their acceptance scenarios with synthetic fixtures.
- **SC-002**: Every format in the FR-021 list has at least one pinned real file, and 100% of pinned real files for a stage's formats pass the hash check and parse without error before that stage closes. `.2dm` MUST have a pinned real SMS file from any public source (adcirc-testsuite or elsewhere) that passes the hash check and parses before B1 closes; if none is found, B1 cannot close and the operator decides.
- **SC-003**: Coloring on every bundled CHILmesh example mesh has zero adjacent elements with the same color and uses no more colors than the largest node-sharing neighbor count plus one.
- **SC-004**: On the 3 exploration meshes (FR-006), at k = 4 and k = 8, the imbalance ratio is at most 1.25 and every part is one connected piece. Evidence depends on the FR-006 outcome: if the exploration passes, the recorded exploration table (research.md R3) is the evidence and no CHILmesh partitioner exists to test; if it fails, the CHILmesh partitioner's required-dependency method is tested against the same threshold on the same meshes. If that method also misses the threshold, work stops and the operator decides.
- **SC-005**: Reading one time step from a global output file of a 1-million-node mesh allocates at most 3.0 times that step's array size (returned array `nbytes`), measured as peak `tracemalloc` allocation during the single-step call on a synthetic file. Pass or fail is the comparison of those two numbers.
- **SC-006**: `.2dm`, fort.22, fort.19, fort.20 and fort.24 round-trips give equal values in 100% of fixture cases.
- **SC-007**: No file from adcirc-testsuite is present in the repository at the end of the feature, and the full test suite passes without network access.

## Assumptions

- Synthetic fixtures can be written from the ADCIRC documentation descriptions of each format.
- The `.2dm` reader currently ignores nodestrings; adding them to the reader is part of FR-002. [inference]
- SegmentAnyMesh is private and read-only for this feature; the exploration runs its code against CHILmesh meshes outside both repositories' source trees.
- Stage order is strict: B2 starts after B1's gate, B3 after B2's gate.
