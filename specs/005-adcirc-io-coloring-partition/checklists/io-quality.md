<!-- provenance: author=domattioli model=claude-opus-5-5 effort=medium date=2026-10-12 skill=speckit-pipeline repo=CHILmesh session=4e33ac09-e474-4ebc-a2dd-027a9e0a0442 -->
# Requirements Quality Checklist: ADCIRC I/O, coloring, partitioning

**Purpose**: Test the written requirements (not the code) for completeness, clarity and consistency
**Created**: 2026-10-12
**Feature**: [spec.md](../spec.md)

## Completeness
- [x] CHK001 Are the exact columns read for NWS 20 best track specified, beyond "extra columns"? [Gap, FR-015; resolved: typed set named, other columns verbatim]
- [x] CHK002 Is the meaning of "matching boundary nodes" for fort.19/20 defined in terms of fort.14 IBTYPE values? [Gap, US8; resolved in FR-017, IBTYPE set confirmed in T037]
- [x] CHK003 Is a real-file requirement stated for every supported format? [FR-021, FR-023a]
- [x] CHK004 Is the second dataset of maxele/maxvel (time of maximum) required to be returned? [Gap, FR-009; resolved: yes, record index 1]
- [x] CHK005 Are NetCDF station files (fort.61.nc, fort.62.nc) explicitly in or out of scope? [Ambiguity, FR-012; resolved: in scope, listed]

## Clarity
- [x] CHK006 Is "balanced" quantified? [SC-004, ≤1.25 and one piece]
- [x] CHK007 Is the coloring adjacency rule single and explicit? [FR-004, node-sharing]
- [x] CHK008 Is "about 3 times" in SC-005 precise enough to pass or fail? [Clarity, SC-005; resolved: 3.0x nbytes, peak tracemalloc]
- [x] CHK009 Are the 3 exploration meshes named? [Clarity, FR-006; resolved: donut, structured, block_o]

## Consistency
- [x] CHK010 Do the truncation edge case and FR-011a agree? [Consistency]
- [x] CHK011 Does constitution VI (2DM write has no nodestrings) conflict with FR-001 until amended? [Conflict, plan Constitution Check; resolved: VI sentence is descriptive state, no MUST violated; updated at T018]
- [x] CHK012 Does the stage table gate match FR-023a (fort.20 before B3 close)? [Consistency]

## Acceptance criteria and scenarios
- [x] CHK013 Does every user story have an independent test? [Coverage]
- [x] CHK014 Is the behavior when no real fort.20 can be found stated as a stop for operator decision? [Exception flow, FR-023a; resolved: spec FR-023a and edge case now state it]

## Dependencies and assumptions
- [x] CHK015 Is SegmentAnyMesh read-only and out of the change set? [FR-008]
- [x] CHK016 Is the NOAA ATCF fixture's size limit ("small") quantified? [Clarity, FR-020; resolved: one storm, under 100 KB]
