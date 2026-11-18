<!-- provenance: author=domattioli model=unknown effort=medium date=2026-11-10 skill=speckit-specify repo=CHILmesh session=dispatch-sol-6f5ce8; manual metadata, provenance script unavailable in authorized context -->
# Feature Specification: Fair Head-to-Head Benchmark

**Feature Branch**: Not created or changed (operator instruction)  
**Created**: 2026-11-10  
**Status**: Clarified; plan and tasks authored 2026-11-10  
**Issue**: Follow-up to competitor benchmark work (#301)  
**Input**: Operator mandate for a fair, auditable 1-for-1 comparison on one 16-core, 62 GB Codespace; four meshes; verified outputs; repeatable measurements; rights and claims gates.

## Clarifications

### Session 2026-11-10

Operator decisions (2026-11-10): budgets per class (see Resolved decisions 1); checksum binary64 repr (2); frozen interpreter Python 3.11, one base binary for all venvs, provisional until an install test passes for adcircpy, thalassa, xarray and CHILmesh (FR-007); cold-cache flush happens before the child launches, with `/proc/meminfo` evidence, and FR-009 is read as "after the previous repeat and before the child launches"; every child runs single-threaded with OMP/OPENBLAS/MKL/NUMEXPR thread variables set to 1, recorded per repeat (FR-006).

Resolved from repository evidence. Each bullet names the file read. The operator answered all five open questions the same day; see the decisions above.

- Q: Which host is "the 16-core, 62 GB Codespace", and what do the records already show for it? → A: The host that produced `.scratch/results/w2stofs/` records: `os` "Linux 6.8.0-1064-azure", `cpu_count` 16, `ram_bytes` 67424980992 (62.8 GiB; advertised as 64 GB class), harness Python 3.14.2 (`.scratch/results/w2stofs/w2_STOFS2DGlobal_v2.1.14_chilmesh_0_0.json`). FR-006 host identity MUST record these raw values plus kernel, cgroup memory and CPU limits, and a stable host identifier; "62 GB" is the GiB figure, not a rounded GB figure.
- Q: Were the prior adcircpy counts really zero, meaning an empty parse? → A: No. The adcircpy adapter printed correct counts (Shinnecock 3070/5780, WNAT_Hagen 52774/98365, ENPAC2003 272913/531680; `stdout_tail` in `.scratch/results/w2cold/w2_*_adcircpy_0_1.json`), matching CHILmesh. The zeros came from the harness: `benchmarks/competitors/w2/run.py:46-48` defaults `n_nodes=0, n_elems=0` and `run.py:155-163` never copies the adapter's counts into the record. FR-015 therefore requires that counts and checksums in a repeat record come from the tool's materialized output as emitted by the adapter, never from a harness default; a record whose counts are not traceable to adapter output is invalid. The Assumptions bullet about zero counts is corrected accordingly.
- Q: What does the current harness measure for memory, and does it meet FR-010/FR-011? → A: It does not. `benchmarks/competitors/measure.py:47-57` takes one `ru_maxrss` from `os.wait4` per child: a whole-process peak with no baseline, no region samples and no sampling interval. It is kept only as a kernel-accounted whole-process cross-check. The FR-010 shared sampler is an external poller in the parent that reads the child's resident set (Linux `/proc/<pid>/status` VmRSS) at a fixed interval and records every sample with a timestamp; the same code path already launches CHILmesh workers and competitor adapters (`run.py:33-40`, `run.py:155-156`), so one sampler serves every tool.
- Q: Can operation RSS growth be measured across the parse call for a child-process tool such as adcircpy? → A: Yes, because every adapter is harness-owned code (`benchmarks/competitors/adapters/adcircpy_read.py`, `thalassa_read.py`, `xarray_read.py`, `w2/measures.py`). The adapter emits a region-start marker immediately before the operation call and a region-end marker immediately after complete materialization, each with its `time.perf_counter()` stamp, on a channel the parent reads (stderr line or pipe). The parent maps markers to sample indices: baseline = last sample before region-start; operation-region peak = max sample between markers; whole-process peak = max of all samples from launch to region-end. Verification (counts, checksums) runs in the child after region-end, so its memory is outside both windows (FR-011). Limitation to disclose: a sampled peak misses spikes shorter than the interval, and RSS growth includes interpreter allocator retention and file-backed pages a tool maps, so growth is an upper bound on mesh memory. An opaque tool with no harness-owned adapter (OceanMesh2D through Octave) cannot emit markers and is not-comparable for memory.
- Q: Is import exclusion already in place, and how are deferred imports detected? → A: Imports already precede the timer: `adcircpy_read.py:19-21` imports adcircpy before `t0`; `w2/measures.py:259-262` resolves `from chilmesh import CHILmesh` before `t0`; `thalassa_read.py:18-21` and `xarray_read.py:18-20` likewise. Detection of deferred imports (Edge Cases): the adapter snapshots `sys.modules` keys at region-start and region-end and reports the difference as `new_modules`; a non-empty list marks the repeat invalid for isolated operation comparison, as the edge case already states.
- Q: How is elapsed operation time measured uniformly? → A: In the child, `time.perf_counter()` around the operation call only, which every adapter already does (`adcircpy_read.py:21-23`, `measures.py:261-263`, `thalassa_read.py:21-23`, `xarray_read.py:20-23`). The parent's whole-process wall time (`measure.py:56`) is recorded separately and labelled as including startup and imports. FR-011 elapsed operation time means the in-child figure.
- Q: Is the RSS scope child-only or child plus descendants? → A: Child-only. Every adapter is a single Python process. The sampler also records the child's descendant pids (Linux `/proc/<pid>/task/*/children`) at each sample; any descendant invalidates the memory result for that repeat (FR-010), since `ru_maxrss` from `wait4` and VmRSS of one pid do not cover it.
- Q: Does any record prove a cold cache today? → A: No. Records carry only the label `cache: cold|warm|unknown` (`benchmarks/competitors/records.py:16`), set from the `--cache` flag (`cli.py:25`); the README asks the operator to flush by hand (`benchmarks/competitors/README.md:54-61`). FR-009 evidence on the campaign host (Linux) is: the harness itself runs `sync` then writes `3` to `/proc/sys/vm/drop_caches`, and records the command, exit code, timestamp, and `/proc/meminfo` `Cached:` before and after (after must be lower). A repeat without that evidence is `cache: unknown` and invalid. Cache-reset ordering is resolved: flush before child launch.
- Q: Which interpreter does each tool use today? → A: Different ones. Local records show the harness on Python 3.14.7 (`machine.python`) while the adcircpy and thalassa venvs run Python 3.12 (`.scratch/cenvs/adcircpy/lib/python3.12/...` in `stderr_tail` of `.scratch/results/w2/w2_WNAT_Hagen.14_adcircpy_0_1.json`; thalassa likewise in `w2cold`). On the Codespace, adcircpy and thalassa were never configured (`reason: no competitor environment configured` in `w2stofs`). The adcircpy import on local 3.12 failed once with `ModuleNotFoundError: No module named 'distutils'` from an unpinned geopandas (`.scratch/results/w2/w2_WNAT_Hagen.14_adcircpy_0_1.json`). `machine.py:34-43` records only `sys.version` of the parent, never the child's build or binary digest. FR-007 evidence is therefore recorded in each child: `sys.version`, `sys.implementation`, `sysconfig.get_platform()`, and the SHA-256 of `os.path.realpath(sys.executable)`; every venv is created from the same base binary with `<base python> -m venv`, so realpath resolves to one file. Interpreter version is resolved: Python 3.11.
- Q: Does adcircpy expose a version for the pin record? → A: Not through `__version__`: every ok adcircpy record shows `"version": "unknown"` (`adcircpy_read.py:26` falls back). FR-007 pins MUST come from `importlib.metadata.version()` and a frozen requirements lock per environment, not from module attributes. `envs/requirements-adcircpy.txt` has no pins today.
- Q: Does Thalassa read fort.14? → A: No. On a fort.14 path `thalassa.api.open_dataset` raises `ValueError: did not find a match in any of xarray's currently installed IO backends` (`.scratch/results/w2cold/w2_ShinnecockInlet_coarse.14_thalassa_0_2.json`). Thalassa and xarray are both NetCDF-only; FR-004 applies to both, and the prior "thalassa failed: exit code 1" cells are reclassified not-comparable on fort.14 rather than failed.
- Q: Which tools have a W2 read adapter at all? → A: Only chilmesh, adcircpy, thalassa and xarray (`benchmarks/competitors/w2/tools.py:32-42`). oceanmesh and OceanMesh2D are W1-only tools with no W2 adapter; under FR-001 and FR-016 every W2 cell for them is `not-comparable` with reason "no mesh-read operation in scope; W1-only tool", never a blank.
- Q: What is the W2 operation inventory in the harness? → A: `parse`, `adjacency` (cpp and python backend variants), `locate` (tree build and query), `render`, `stations`, `write`, `full-init` (`w2/measures.py:249`). `write`, `full-init` and both adjacency variants are CHILmesh-only cells (FR-005); competitors answer only `parse` (`tools.py:35,38`). FR-002's inventory is this list.
- Q: What are the Valence mesh names and the STOFS version? → A: `ShinnecockInlet_coarse.14`, `WNAT_Hagen.14`, `EasternPacific_ENPAC2003.14`, `STOFS2DGlobal_v2.1.14` (record file names under `.scratch/results/w2stofs/`). STOFS served is v2.1 with 12,785,004 nodes and 24,875,336 elements, not the v3.1 assumed in 006; FR-014 binds that version.
- Q: What STOFS cost is already measured on the campaign host? → A: CHILmesh parse 96.17 s in-child, 96.43 s whole-process wall, `ru_maxrss` peak 9,950,588,928 bytes (9.27 GiB), cold label, repeat count 1 (`w2stofs/w2_STOFS2DGlobal_v2.1.14_chilmesh_0_0.json`). This is the evidence base for the budget options in Q1; it is not a result of this campaign.
- Q: Is a local `rights.json` entry enough for FR-018? → A: No. `benchmarks/competitors/valence.py:91-102` admits any mesh whose entry has `publish_results: true`. FR-018 requires the entry to also carry the ruling citation (Valence issue or PR number and date) and the exact mesh version it covers; the report prints the citation next to each admitted mesh, and an entry without a citation excludes the mesh.
- Q: How does the current harness decide OOM? → A: It infers it from SIGKILL (`measure.py:65-66`, "OOM suspected"). FR-016 evidence for OOM is the host cgroup `memory.events` `oom_kill` counter delta across the run (plus `memory.peak` where present); SIGKILL without that delta is "killed by signal, cause unknown".
- Q: How are the numeric values stored by CHILmesh's own reader, for the checksum policy? → A: As IEEE-754 binary64 via Python `float()` of each token (`src/chilmesh/fort14_io.py:123-125`); no tool retains the source decimal text. This bounds the options in Q2.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Trust each comparison (Priority: P1)

A coastal modeler wants to compare CHILmesh and competitors doing the same work on the same mesh, with evidence that both tools loaded the actual mesh correctly.

**Why this priority**: Timing an empty or incomplete result, or comparing different operations, cannot establish an advantage.

**Independent Test**: Review a candidate comparison and deliberately substitute an empty parse, altered connectivity, and a converted input. Each is either rejected or placed in its correctly labelled category without producing a like-for-like advantage claim.

**Acceptance Scenarios**:

1. **Given** identical source bytes and equivalent native mesh-read operations, **When** counts and node/element checksums match the independent reference, **Then** the cell is verified; eligibility for claims still depends on measurement, stability and rights gates.
2. **Given** a nonempty source mesh, **When** adcircpy reports zero nodes or elements, or a tool returns correct counts but altered coordinates or connectivity, **Then** the cell is invalid and cannot support an advantage claim.
3. **Given** Thalassa or xarray reading NetCDF derived from fort.14, **When** results are reported, **Then** they appear only in a labelled conversion section, never in the fort.14 like-for-like table.
4. **Given** a tool without an equivalent adjacency, query or render operation, **When** the operation matrix is reviewed, **Then** it says not-comparable with a reason; a capability exclusive to CHILmesh is labelled CHILmesh-only.

---

### User Story 2 - Separate operation cost from runtime overhead (Priority: P1)

A benchmark reviewer wants repeated timing and memory measurements that isolate the operation while retaining whole-process memory context.

**Why this priority**: Different import overhead must not masquerade as a mesh-memory advantage.

**Independent Test**: Audit records from tools with different import footprints and a repeat series containing an outlier. Confirm import exclusion, both RSS measures, identical measurement methods, and the unstable flag without needing to trust the scoreboard prose.

**Acceptance Scenarios**:

1. **Given** different tool environments, **When** their interpreter version or build identity differs, **Then** affected measurements are ineligible for comparison until corrected and rerun.
2. **Given** completed imports and operation prerequisites, **When** the operation runs, **Then** elapsed time excludes startup, imports and verification; RSS growth uses the immediately preceding RSS baseline, while whole-process peak is separately labelled as including startup and imports.
3. **Given** Shinnecock Inlet, WNAT_Hagen or ENPAC2003, **When** a successful cell is summarized, **Then** at least five cold-cache repeats and their median and spread are visible; STOFS requires at least three.
4. **Given** a metric with relative spread of exactly 10%, **When** eligibility is checked, **Then** it passes the spread threshold; above 10% is unstable and cannot support a claim for that metric.
5. **Given** a failed cache reset, concurrent benchmark, unrelated workload or missing memory samples, **When** a run is reviewed, **Then** its measurements cannot enter a valid comparison summary.

---

### User Story 3 - See complete results without unsupported claims (Priority: P2)

A reader wants the largest mesh's failures to remain visible and every published advantage sentence to cite eligible evidence from a mesh cleared for publication.

**Why this priority**: Hiding failures or rights exclusions would distort the benchmark even with accurate successful runs.

**Independent Test**: Review a report containing a STOFS timeout, a missing tool, an unstable memory series and WNAT_Hagen. Confirm explicit outcomes for every tool, no WNAT published row, and rejection of unsupported advantage sentences.

**Acceptance Scenarios**:

1. **Given** a STOFS operation that exceeds its fixed time or memory limit, **When** its child terminates, **Then** a failed record names timeout or out-of-memory, with available measurements and diagnostic evidence, rather than dropping the tool or treating failure as a zero.
2. **Given** WNAT_Hagen with no Valence rights ruling, **When** a published table is produced, **Then** no WNAT_Hagen measurement appears; the exclusion and reason remain visible.
3. **Given** an advantage sentence citing an invalid, unstable, non-equivalent or rights-ineligible cell, **When** publication is checked, **Then** publication of that sentence is blocked with the failed gate identified.
4. **Given** all required valid repeats and permitted rights for both tools, **When** an advantage sentence is accepted, **Then** it identifies the mesh, operation, metric, versions and supporting table cells, with arithmetic reproducible from the recorded medians.

### Edge Cases

- Counts match but node coordinates, depth, identifiers or element connectivity differ: checksum verification fails; the result is invalid.
- A checksum reference is unavailable, cannot be independently established, or normalization is unresolved: verification is incomplete and publication of a comparison is blocked.
- A tool defers parsing until first access: the measured read must force complete materialization before the operation ends; lazy handle creation is not an equivalent read.
- Imports trigger deferred imports during the operation: repeat after completing those imports outside the region; otherwise mark the run invalid for isolated operation comparison.
- An operation's peak occurs during imports: whole-process peak retains it; operation growth still uses only samples in the operation region.
- Zero median with nonzero range: relative spread is infinite and unstable; all-zero memory growth has zero spread but cannot be a denominator for a percentage advantage.
- One repeat fails among otherwise successful repeats: retain it; the cell is incomplete for claims until a new, fully disclosed complete series meets the repetition rule. No selective removal of slow repeats.
- An OOM kill gives no reliable cause: record failed with unknown termination cause and available diagnostics; do not infer OOM from exit code alone.
- Octave cannot share the Python interpreter: record its distinct pinned runtime and not-comparable classification; do not silently waive runtime parity.
- Valence does not serve a mesh or a tool cannot start: record not-run with reason, never an omitted row; this does not fulfill a successful comparison.

## Requirements *(mandatory)*

### Functional Requirements

#### Scope and equivalence

- **FR-001**: The campaign MUST cover Shinnecock Inlet (existing harness alias Shinnecock), WNAT_Hagen, ENPAC2003 and STOFS-2D-Global in internal records. It MUST inventory CHILmesh, adcircpy, oceanmesh, OceanMesh2D, Thalassa and xarray, as named by the existing harness. Inclusion in the inventory does not establish equivalence or publication rights.
- **FR-002**: Before measurement, each tool × mesh × operation cell MUST declare input format, input identity, prerequisites, work performed, required output, completion condition and classification: like-for-like, CHILmesh-only or not-comparable. The operation inventory covers parse, adjacency (per backend), locate tree-build, locate query, render, stations, write and full-init. A capability mapping MUST justify equivalence of both work and output; matching metric names alone is insufficient. Unknown capability remains not-comparable pending evidence.
- **FR-003**: A like-for-like comparison MUST use the same input bytes and format, numeric precision policy, operation boundaries, prerequisite state, query/field data and output contract. A lighter operation MUST NOT be compared with a combined operation; adjacency/layer construction cannot be presented as parse alone. Backend identity and actual work MUST be recorded for CHILmesh variants; variants performing equivalent work are separate cells.
- **FR-004**: Thalassa and xarray MUST appear only through a section labelled conversion / NetCDF, outside every fort.14 like-for-like table. It MUST name the source fort.14 SHA-256, converted-file SHA-256, converter/version/options, conversion verification, and separate conversion cost and NetCDF read cost. It MUST show end-to-end conversion-plus-read cost separately if used, and MUST NOT describe a NetCDF read-only figure as a fort.14 read advantage. No direct fort.14 ratio against these tools is eligible for claims.
- **FR-005**: Non-equivalent operations MUST show not-comparable and the reason; CHILmesh-exclusive operations MUST show CHILmesh-only. Neither class may enter a head-to-head speed or memory ratio. Unsupported cells MUST remain in the inventory rather than disappearing from reports.

#### Shared execution conditions

- **FR-006**: All campaign measurements MUST run on ONE identified 16-core, 62 GB Codespace, tools one at a time, without other workload. Records MUST capture host identity, CPU, available cores, RAM and container limits, OS/kernel, thread settings, and execution order. The same resource limits and thread policy MUST apply to each like-for-like cell. An overlapping tool, other user job, detected unrelated CPU/I/O load, or changed resource allocation invalidates the affected run; load observations before and throughout each run MUST be retained. Ordinary idle system services are not another workload.
- **FR-007**: Every Python tool MUST use the same interpreter version AND build, evidenced by executable digest, full version/build identity and architecture. Separate dependency environments are allowed only with that identical interpreter build. Every tool, direct and transitive dependency set, native runtime, converter and CHILmesh revision MUST be pinned and recorded with exact versions and environment identity before running. A missing or mismatching pin blocks comparison eligibility. Non-Python tools MUST record their distinct runtime/version and remain not-comparable to Python tools under this campaign's strict interpreter rule.
- **FR-008**: Tools MUST remain isolated in fresh child processes for each repeat; GPL competitors MUST not be imported into CHILmesh. Startup, all necessary imports, input acquisition/conversion, and independently prepared operation prerequisites MUST finish before the timed and operation-memory region. Input parsing required by the operation MUST NOT be hidden in preparation. Region end MUST represent complete required output, including eager materialization. Output verification occurs after the region and MUST not inflate measured cost.
- **FR-009**: Every valid repeat MUST have documented, successful cold file-cache reset after all hashing, reference parsing and conversion, and before the child process launches; imports run inside the child after the flush and stay outside the timed region. Merely using a fresh process or labelling a run cold is insufficient. Each repeat MUST start without retained tool/application caches or prior operation state; failure to establish cold conditions makes it invalid. Mesh fetch and cache-reset duration MUST be excluded from operation timing.

#### Uniform time and memory measurement

- **FR-010**: Every tool MUST use the SAME external child-process RSS measurement method and SAME sampler version, sampling interval, process scope, units and boundary synchronization. The record MUST state whether scope is child-only or child plus descendants; descendant use MUST be captured consistently or invalidate equivalence. Tools MUST NOT substitute self-reported memory, a different peak source, or interpreter-specific allocation tracing for RSS. Sampling coverage and limitations MUST be disclosed, including that a sampled peak may miss shorter spikes.
- **FR-011**: Each repeat MUST report elapsed operation time and three RSS quantities: baseline immediately before operation start; whole-process peak from child launch through operation completion (including startup/imports); and operation-region peak. It MUST also report operation RSS growth = max(0, operation-region peak − baseline). Only region samples contribute to operation-region peak. Verification memory MUST be excluded from both reported peaks by ending their observation windows before verification. Missing baseline, samples or region boundaries invalidate the memory result. Whole-process peak MUST be labelled as including runtime/import overhead; operation growth MUST NOT be described as total mesh memory or an allocation count.
- **FR-012**: Each small-mesh cell (Shinnecock Inlet, WNAT_Hagen, ENPAC2003) MUST have at least five valid cold-cache repeats; each STOFS cell MUST have at least three. The campaign MUST preregister repeat counts, execution ordering and the same per-operation time and memory limits across tools before measurement. All attempts, failures and reruns MUST remain visible. A supported cell with failed or invalid attempts is not claims-eligible merely because successes were selectively retained; a replacement complete series requires an explicit rationale and retains the original series.
- **FR-013**: For elapsed time, whole-process peak RSS and operation RSS growth separately, the cell MUST report raw repeat values, sample count, median, minimum, maximum, absolute spread (maximum − minimum), and relative spread (100 × absolute spread / median). Relative spread >10% MUST flag that metric unstable; equality to 10% passes. Zero median and zero range gives 0%; zero median and positive range gives infinite spread. No failed, invalid or not-run value may be replaced with zero or included in a numerical summary. An incomplete series may show a descriptive summary of its successful repeats, labelled incomplete and never eligible. A cell is stable only when every reported comparison metric is stable; an unstable cell cannot support any advantage sentence.

#### Parse verification and outcome completeness

- **FR-014**: Each mesh reference MUST bind Valence provenance, mesh version, raw-file SHA-256, fort.14 header node/element counts, and independent canonical node and element SHA-256 checksums. Header counts MUST be read from the bytes identified by that file digest. Reference checksums MUST be established independently of any timed tool, never adopted from CHILmesh simply because it is the subject tool.
- **FR-015**: For EVERY parse repeat, actual materialized output MUST match both reference header counts AND both canonical checksums. Nodes MUST cover identifiers, coordinates and depth; elements MUST cover identifiers, type/arity and connectivity. Canonicalization MUST reconcile documented storage indexing without losing identity, ignore only irrelevant record ordering, and detect missing, duplicated or altered nodes/elements. A declared, shared numeric normalization policy MUST apply to every tool and the reference. Missing evidence, zero counts for a nonempty reference, count mismatch, or checksum mismatch makes the cell invalid and excludes it from ratios, claims and successful-repeat counts. Verification of derived conversion inputs MUST also prove the same source mesh content. Non-parse operations MUST use verified mesh inputs and meet the shared output contract in FR-002.
- **FR-016**: Every planned tool × mesh × operation cell MUST retain an explicit outcome: ok, invalid, failed, not-run, or not-comparable (with CHILmesh-only where applicable), with reason and evidence. Execution outcome, equivalence class, verification, stability and publication eligibility MUST remain distinct properties. In particular, STOFS MUST record an outcome for EVERY inventoried tool, even if none succeeds. Actual OOM or timeout MUST be failed, with reason, configured limits, elapsed time, available last/peak RSS, termination diagnostics and missing-measurement flags; it MUST NOT become not-run, a blank, or a numeric speedup. Unsupported operations and unavailable tools need explicit reasons and cannot masquerade as successful STOFS coverage.
- **FR-017**: Records MUST preserve campaign/cell/repeat identifiers, commands/options, mesh and auxiliary-data hashes, tool/environment pins, operation contracts, host and load evidence, cache-reset evidence, region boundaries, raw measurements, verification evidence, failures, aggregation inputs and exclusions. A reviewer MUST be able to trace every displayed number and eligibility decision to immutable raw records and reproduce its arithmetic without running the benchmark.

#### Rights and claims gates

- **FR-018**: Meshes MUST be acquired from Valence with provenance, retaining the prior no-fallback rule. Only a documented Valence rights ruling explicitly permitting publication of benchmark results for the exact mesh/version may admit a mesh into published tables. Missing, ambiguous, denied or version-mismatched rights MUST exclude it. A locally written allow flag alone is insufficient. WNAT_Hagen currently has NO ruling and MUST stay out of published tables, charts, numerical summaries and advantage sentences; it may remain in internal records. Reports MUST list exclusions and reasons without leaking excluded measurements. Rights for the other three meshes MUST be checked, never presumed.
- **FR-019**: Every proposed published sentence asserting a speed or memory advantage MUST pass a claims gate: both supporting cells are like-for-like, verified, complete, stable, measured under FR-006–FR-013 and rights-eligible. Each sentence MUST cite its cells, mesh/version, operation, metric, tool versions, repeat counts, medians and spread. Ratios MUST use the same named metric and unit on each side, with explicit denominator and direction; zero denominators prohibit percentage/ratio claims. Memory prose MUST distinguish whole-process peak from operation growth. No claim may generalize a result beyond its measured meshes, operations, versions and resource conditions.
- **FR-020**: The prior “40 to 60 percent less memory than adcircpy” statement MUST be treated as unsupported background, not carried into this campaign's conclusions or used as its expected outcome. Fresh eligible cells are required for any new claim. Equal performance, worse performance and universal failure are acceptable results; there is no required CHILmesh win. Reporting MUST preserve eligible failures and caveats even when they weaken the narrative.

### Key Entities

- **Campaign**: One host and pinned runtime manifest, predeclared operation matrix, limits, order, repeat counts and measurement policy.
- **Mesh reference**: Named/versioned Valence asset, byte digest, header counts, canonical content checksums and publication-rights ruling.
- **Operation contract**: Same work, inputs, prerequisites and output obligations used to justify or reject equivalence.
- **Repeat record**: One tool/mesh/operation attempt, raw metrics, verification evidence, environment/cache/load evidence and outcome.
- **Comparison cell**: Complete repeat series and metric summaries, with separate equivalence, validity, stability and rights states.
- **Claim**: Proposed sentence linked to supporting cells, calculation and gate decision.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Internal coverage accounts for 100% of the four-mesh, six-tool operation inventory with a measured result or explicit outcome and reason; STOFS has zero silently omitted tools, including after timeout or memory exhaustion.
- **SC-002**: 100% of successful comparison cells meet the five-repeat small-mesh or three-repeat STOFS minimum and expose medians, ranges and relative spreads. Every >10% spread is visibly unstable and excluded from claims.
- **SC-003**: 100% of cells admitted as verified have matching source counts and both content checksums for every parse repeat. Empty outputs, same-count content corruption and missing verification are rejected in acceptance review.
- **SC-004**: 100% of eligible comparisons show matching host/runtime conditions, import-excluded operation timing, both clearly distinguished memory views and one shared measurement policy; no mismatched environment or contaminated run contributes an eligible cell.
- **SC-005**: Publication review finds zero fort.14 like-for-like rows for NetCDF-only tools, zero ratios from non-equivalent operations, and zero published measurements for meshes lacking a permitting Valence ruling. WNAT_Hagen has zero published measurement rows while its ruling is absent.
- **SC-006**: A reviewer can independently recalculate 100% of published aggregates and advantage figures from cited raw records; every deliberately introduced invalid, unstable, non-equivalent or rights-ineligible advantage sentence is blocked.

## Assumptions

- This feature strengthens the W2 evidence contract from [006](../006-competitor-benchmarks/spec.md) and the [harness README](../../benchmarks/competitors/README.md). Conflicting older host, repeat, measurement, rights or equivalence rules are superseded for this campaign only; older artifacts are not edited here.
- The operator's account of the prior scoreboard is the problem statement, not newly verified measurement evidence: whole-process RSS mixed runtime/import overhead, a 144–249 MB gap persisted as mesh size grew about 90-fold, environments differed, and adcircpy record counts were zero. Clarification (Session 2026-11-10) traced the zeros to the harness record builder, not to an empty parse; the adapter output carried correct counts. Those observations motivate gates; they establish no new advantage.
- “62 GB” and “16-core” identify the required Codespace class. Actual advertised and enforced byte/core limits must be recorded; this spec does not guess a GB/GiB conversion or exact CPU model.
- All six harness tools remain visible in coverage, including non-equivalent tools; strict Python build parity applies to Python comparisons. OceanMesh2D's Octave runtime remains explicitly non-comparable under that rule rather than being silently exempted.
- Relative spread means range divided by median, chosen conservatively so a single slow or high-memory repeat remains visible. This is a spread flag, not a confidence interval or proof of statistical significance.
- Tool/runtime version numbers, numeric normalization and resource-limit values are not invented here. Exact pins and policies must be settled and frozen before any campaign can produce eligible results.
- Valence availability, independent reference generation, approved publication rights and cold-cache control on the target Codespace are dependencies. An unavailable dependency produces explicit blocked/not-run records, not a substitute source or a claimed result.

## Resolved decisions (clarify 2026-11-10)

1. **Resource budgets** (resolved, operator 2026-11-10): small meshes 600 s wall and 8 GiB RSS per operation; STOFS 1800 s wall and 48 GiB RSS per operation. Same limits for every tool. The parent kills the child when a sampled RSS or the wall clock passes the limit; the record names the limit that fired. Disclose that spikes shorter than the sampling interval can pass unseen.
2. **Checksum normalization** (resolved, operator 2026-11-10): binary64 shortest round-trip form. Node line `id x y depth`, element line `id arity n1 n2 ...`, 1-based file ids, values via `repr(float(token))`, ascending id, `\n` separated, SHA-256 over UTF-8, nodes and elements separately. The reference comes from an independent minimal parser, not CHILmesh. A token with more than 17 significant digits makes the cell "verification incomplete".

## Out of Scope

Workflow 1 (including smoothing, repair and quality comparisons), Workflow 3, the C++ backend fix, and the fort.14 IBTYPEE bug (#309) are excluded. This feature has specification, plan, research, data-model, quickstart, tasks and review artifacts; implementation, benchmark execution and publication remain outside this documentation remediation.
