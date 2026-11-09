<!-- provenance: author=domattioli model=claude-sonnet-5-5 effort=medium date=2026-11-08 skill=speckit-specify repo=CHILmesh session=2d051659 -->
# Feature Specification: Competitor Head-to-Head Benchmarks (W2 then W1)

**Feature Branch**: n/a (stays on current branch by operator instruction)
**Created**: 2026-11-08
**Status**: Draft
**Input**: GitHub issue domattioli/CHILmesh#301, Workflows 1 and 2. Workflow 3 deferred.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Operational-scale post-processing on STOFS-2D-Global (Priority: P1)

A coastal modeler wants to know whether CHILmesh can load, index, query and plot the NOAA STOFS-2D-Global ADCIRC mesh (12,785,004 nodes, 24,875,336 elements) and how it compares with adcircpy, Thalassa and raw xarray/netCDF4 on the same machine.

**Why this priority**: Needs no paid license or solver. Tests the largest public claims (memory, O(log n) point location). Closes #276. Most likely place for failure (1.73 GB fort.14 read), so run first.

**Independent Test**: Run the W2 benchmark on the STOFS grid and outputs; a results table with every metric for every tool is produced.

**Acceptance Scenarios**:

1. **Given** the STOFS grid and a cold file cache, **When** each tool reads the grid by its documented default path, **Then** parse time and peak memory are recorded per tool, or a recorded failure with the cause.
2. **Given** a loaded mesh, **When** CHILmesh builds adjacency with the C++ and Python paths, **Then** both build times are recorded.
3. **Given** the CO-OPS station list, **When** 1,000 point-location queries run, **Then** total and per-query time are recorded for each tool that supports it.
4. **Given** the STOFS `points.cwl.nc` file, **When** values are sampled at the stations, **Then** the maximum absolute difference from the file is recorded.
5. **Given** a maximum-elevation field, **When** CHILmesh renders it, **Then** render time is recorded.

---

### User Story 2 - Pre-run mesh QA, repair and timestep gate (Priority: P2)

An ADCIRC user wants to know whether CHILmesh's quality, smoothing and Courant/timestep tools match or beat OceanMesh2D (MATLAB or Octave) and oceanmesh (Python) before a run.

**Why this priority**: Every ADCIRC user does this step, but it needs a MATLAB/Octave seat and carries a license boundary, so it follows W2.

**Independent Test**: Run Track 1 on Shinnecock Inlet; the quality, Courant and cost metrics table is produced for each tool.

**Acceptance Scenarios**:

1. **Given** Shinnecock Inlet (3,070 nodes, 5,780 elements), WNAT_Hagen and ENPAC2003, **When** each tool smooths with topology frozen (Track 1), **Then** per-mesh element quality (minimum, mean minus 3 sigma, against the 0.75 bar), percent of angles below 10 degrees, percent irregular valence, maximum Courant number at fixed dt, largest stable dt, wall time and peak memory are recorded.
2. **Given** any Track 1 run, **When** it finishes, **Then** boundary-node displacement is 0 for every tool, or the run is marked invalid.
3. **Given** the same meshes, **When** each tool performs a full clean (Track 2), **Then** the same metrics plus count of vertices deleted are recorded and reported as secondary, with the note that CHILmesh does not delete vertices.

---

### User Story 3 - Correct disputed documentation claims (Priority: P3)

A reader of the README and benchmark docs wants figures that match measured data.

**Why this priority**: Small, independent of any benchmark run.

**Independent Test**: Read the three changed passages; each matches the cited evidence.

**Acceptance Scenarios**:

1. **Given** docs/BENCHMARK.md, **When** corrected, **Then** WNAT is stated as 52,774 nodes (about 105k DOF) and STOFS as 12.8M nodes (about 25.6M DOF).
2. **Given** README.md, **When** corrected, **Then** the speedup range cites one mesh per figure with its version, and the v1.4.1 ENPAC figure (about 14.8x) is not outside the stated range.
3. **Given** README.md, **When** corrected, **Then** the Python requirement reads 3.11+ everywhere.

---

### Edge Cases

- STOFS read exhausts memory or exceeds the time limit: record the failure, peak memory reached and the cause; do not silently drop the tool.
- A competitor cannot run (no MATLAB seat, Octave incompatibility, SMS unavailable): mark the cell "not run" with the reason; the report still publishes.
- A competitor has no equivalent for a metric: mark "not applicable".
- Mesh rights do not allow publishing results for a mesh: that mesh is excluded and the exclusion is listed.
- STOFS grid version differs (v1.1 has 12,784,991 nodes): the version used is stated in results.
- Warm cache contaminates a run: runs repeat with documented cache clearing; a run without it is flagged.

## Clarifications

### Session 2026-11-08

- Q: Which machine runs the STOFS benchmark? -> A: A machine other than this Mac or a cloud host; the harness MUST be portable and script-driven, and MUST NOT assume local paths or this machine.
- Q: How do meshes reach the harness? -> A: Valence API only. The harness MUST fetch meshes (with provenance) from Valence and block, with a clear message, when a mesh is not yet served. No committed mesh copies, no direct NOAA mesh download. Auxiliary NOAA data (`points.cwl.nc`, CO-OPS station list) is exempt: the harness fetches it directly and records URL, SHA-256 and fetch date (plan D1).
- Q: How is OceanMesh2D run? -> A: GNU Octave attempt; on failure the cell is "not run" with the reason. oceanmesh (Python, separate process) always runs. No MATLAB seat assumed.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The benchmark MUST run W2 and publish its results before W1 results are published.
- **FR-002**: W2 MUST measure parse time, peak memory, adjacency build time (C++ and Python), 1,000-query point-location time, render time, and maximum absolute difference from `points.cwl.nc` at the stations.
- **FR-003**: W2 MUST compare against adcircpy, Thalassa and raw xarray/netCDF4 on the same machine, with a cold cache, each using its documented default path.
- **FR-004**: W1 MUST report Track 1 (frozen topology) as the headline result and Track 2 (full clean) as secondary.
- **FR-005**: W1 MUST record every metric listed in User Story 2, scenario 1, for CHILmesh, OceanMesh2D and oceanmesh; SMS is optional.
- **FR-006**: Every Track 1 run MUST verify boundary-node displacement is 0 and mark violations invalid.
- **FR-007**: GPL-licensed competitors (oceanmesh, adcircpy) MUST run only as separate processes; no GPL code is imported into or distributed with CHILmesh.
- **FR-008**: Benchmark meshes MUST come from Valence only, with provenance; the spec MUST NOT publish results for a mesh whose rights do not allow it.
- **FR-009**: Each result MUST record machine description, tool versions, mesh version, and the command used, so a reader can repeat it.
- **FR-010**: Failures and not-run cells MUST appear in the results with a reason.
- **FR-011**: The W2 outcome MUST state whether the fort.14 reader needs a code change (closes #276).
- **FR-012**: The three doc corrections in User Story 3 MUST be made.
- **FR-013**: The benchmark MUST test the "O(log n) point location" claim by recording query time across mesh sizes from about 5k to 25M elements.
- **FR-014**: Workflow 3 (partition and tidal hindcast skill) is out of scope.

### Key Entities

- **Benchmark mesh**: a named mesh with node/element counts, source, version, and rights status (Shinnecock, WNAT_Hagen, ENPAC2003, STOFS-2D-Global).
- **Tool**: a CHILmesh or competitor entry with version, license, and run mode (in-process or separate process).
- **Run record**: tool, mesh, track, metrics, machine, command, status (ok, failed, not run, invalid).
- **Results report**: tables per workflow and track, with failures and caveats.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: The W2 report has a metric value, a recorded failure or not-run reason, or a justified "not applicable" for every tool and every metric in FR-002. Display labels ("not run", "not applicable") map to machine statuses `not-run` and `na` in the run record.
- **SC-002**: The W1 Track 1 report covers all three meshes, with zero boundary displacement verified for every valid run.
- **SC-003**: A second person can repeat any published result from the recorded machine, versions and command and land within 10 percent on timing.
- **SC-004**: Point-location query time is reported at 3 or more mesh sizes between 5k and 25M elements, so the O(log n) claim is supported or refuted.
- **SC-005**: No GPL source file appears in the CHILmesh repository or wheel, checked by a repository scan.
- **SC-006**: After the doc fixes, no README or BENCHMARK passage contradicts the measured node count, DOF or speedup figures, and no "Python 3.10" mention remains.

## Assumptions

- W2 runs on a separate, larger machine or cloud host (about 32 GB RAM assumed); the harness is portable.
- STOFS grid v3.1 (same size and counts as v2.1) is used; the version is recorded.
- Valence ingestion of Shinnecock, WNAT_Hagen, ENPAC2003 and STOFS is an open dependency in Valence. Benchmarks fetch meshes only through Valence, so real runs wait on it.
- No MATLAB seat is assumed. OceanMesh2D gets a GNU Octave attempt; if the attempt fails, the cell is "not run" with the stderr tail as reason (see Clarifications). oceanmesh (Python) always runs as a separate process.
- Whether Valence rights allow publishing benchmark results is not yet checked.
- Not checked: OceanMesh2D GMD paper details, SMS thresholds, adcircpy license and version.
