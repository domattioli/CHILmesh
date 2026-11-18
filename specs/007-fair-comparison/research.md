<!-- provenance: author=domattioli model=unknown effort=medium date=2026-11-10 skill=speckit-plan repo=CHILmesh session=dispatch-sol-78a08f; manual stamp, exact provider model ID unavailable; transport rung=sol -->
# Research: Fair comparison

Status: repository-only design research. No network, package installs, benchmark runs or host privilege probes. Decisions below are implementation choices unless explicitly resolved by the operator in spec.md. No untested capability is treated as available.

## Repository evidence

- `measure.py:16,31-76`: 20 ms wait4 loop, direct-child timeout kill, SIGKILL labelled suspected OOM; no region samples.
- `records.py:14-23,63-99,101-126`: schema v1; only four statuses; first matching tool/metric wins, losing series identity.
- `w2/run.py:33-58,87-150,155-174,210-239`: common launch path but bundled operation evidence, default zero counts, competitor output counts not propagated, arbitrary cache label and repeated worker scheduling.
- `w2/tools.py:32-47`: only four W2 tools, Thalassa incorrectly marked fort.14, no six-tool inventory.
- `adapters/adcircpy_read.py:19-28`, `thalassa_read.py:18-27`, `xarray_read.py:18-28`: imports before time; module-attribute versions; only xarray explicitly forces load.
- `w2/measures.py:249-289`: CHILmesh worker operation inventory; initial mesh read before operation dispatch; locate combines build and query, requiring split regions/processes for unambiguous costs.
- `valence.py:33-80,91-102`: streaming Valence-only acquisition exists; rights gate currently accepts a bare publish_results flag.
- `report/render.py:35-71`, `cli.py:18-73`: publication currently filters cache label and optionally rights, not verification or complete series.
- Existing tests: `tests/benchmarks/test_measure.py`, `test_records.py`, `test_w2_wiring.py`, `test_bench_cli.py`, `test_report_render.py`, `test_valence.py`, `test_gpl_scan.py` provide regression integration points.

## D1 — Independent reference and canonical bytes (FR-014/015)

**Decision**: Minimal standard-library reference parser, importing no timed tool, independently reads raw fort.14 bytes, binds file digest and header counts, validates unique IDs/arity/connectivity, sorts by original numeric IDs, and hashes nodes and elements separately. Canonical line fields are single-space separated; records are joined by LF with no trailing LF (design choice fixing the unspecified terminal separator); UTF-8 bytes. Node numeric fields use repr(float(token)); connectivity order/winding is retained. Preserve signed zero; reject nonfinite values and malformed tokens. Count significant decimal digits from the mantissa, excluding leading zeros and exponent; conservatively retain trailing digits, and mark >17-digit input verification incomplete. Reference schema records this algorithm revision. Use bounded chunks and external merge sorting outside measurement for STOFS rather than retaining tens of millions of Python records.

**Rationale**: Header/count-only verification misses corruption; a timed tool cannot certify its own reference. Explicit IDs and immutable references prevent row-order and index drift. Adapter normalization may share the canonical byte encoder but cannot establish the source reference using CHILmesh.

**Alternatives considered**: Raw text hashes reject harmless decimal notation changes; decimal rounding loses binary64 differences; using CHILmesh as reference violates independence. Unknown source identity mappings never default to +1 unless the mapping is proven contiguous.

## D2 — One sampler and region protocol (FR-008/010/011; hole 11)

**Decision**: One Linux child-only VmRSS sampler in measure.py, fixed 0.010 s interval, same version for every tool. Dedicated control pipe READY/GO/START/END/ACK, child and parent perf_counter stamps with clock metadata. START immediately precedes operation; END follows complete output. END acknowledgement gates subsequent verification. Before GO take a baseline; periodic samples classify by timestamps, not delayed pipe-receipt time. Record every sample including descendant IDs and gaps; >0.020 s gaps or absent in-window samples fail memory coverage. Whole-process sampled peak stops at END, unlike the separately retained wait4 lifetime maximum. Preserve child output in a structured result artifact; stdout/stderr tails are diagnostics only.

**Rationale**: Approximately 20 samples for the previously observed 0.2 s parse, without claiming sub-interval precision. The handshake prevents verification racing the observer; gap checking exposes sampling stalls. A child finishing between ticks has invalid sampled memory rather than an invented peak.

**Alternatives considered**: wait4 alone lacks baseline/region boundaries and includes verification. tracemalloc is not RSS. A different sampler per adapter defeats equivalence. 20 ms gives fewer samples; 1 ms increases controller overhead.

## D3 — Lifecycle, limits and OOM (holes 8/9)

**Decision**: Popen with start_new_session=True, process-group SIGKILL on timeout/RSS breach, reap and close descriptors on every path. Watchdog launch-to-exit 600 s/8 GiB for small meshes, 1800 s/48 GiB for STOFS; elapsed operation time remains START..END. Use strict exceeds checks for sampled RSS/deadline. Resolve cgroup from `/proc/<pid>/cgroup` and mountinfo, snapshot memory.events before/after and memory.peak if available. Require readable memory.events in campaign preflight; unavailable evidence yields explicit blocked/not-run cells. Label kernel OOM only with oom_kill delta and attributable child termination; shared-cgroup noise or concurrent load makes attribution uncertain and the measurement invalid. A sampler kill is rss-limit, SIGKILL without evidence is unknown, never OOM by default.

**Rationale**: A child-only kill leaks descendants. Kernel counters provide cause evidence while sampler limits implement the chosen unprivileged enforcement. memory.peak may be historical in a shared cgroup and must be labelled scope/lifetime rather than as child peak.

**Alternatives considered**: cgroup memory.max enforcement needs unverified write permissions; RLIMIT_AS is virtual address space, not RSS. Assuming any SIGKILL means OOM is false. Unbounded import/verification hangs are prevented by the shared wall watchdog.

## D4 — Cache and host contamination (FR-006/009; hole 10)

**Decision**: Cache reset before child launch, per resolved operator decision; record ordered sync/write argv or operation, exit statuses, UTC/monotonic timestamps, Cached KiB before/after, converted bytes and launch timestamp. Require Cached_after < Cached_before. Permission failure is an invalid cold attempt, never a labelled success. Precompute hashes/references/conversion first, and do not open the measured file between flush and START. New child means no retained application state. Acquire exclusive host campaign lock; record process inventory, CPU and I/O deltas, loadavg, proc/stat, proc/diskstats and own process accounting for one second before launch and every 100 ms through exit. Freeze idle-service allowlist from preflight, never auto-whitelist later jobs. Any second benchmark/user job, host identity/allocation change, or unexplained activity invalidates the repeat. Preregister a conservative detector: unrelated CPU >5 percentage points of one core or unrelated I/O >1 MiB/s in a 100 ms window is contamination; unattributable evidence is unknown and ineligible. Record thresholds and observations; thresholds are operational detectors, not proof of absolute host idleness.

**Rationale**: loadavg alone includes the measured child and cannot distinguish an unrelated job; corroborate per-process and total counters. Shared reset sequence prevents tool-specific page-cache treatment. FR-009 uses the same before-launch ordering.

**Alternatives considered**: --cache cold is only a label; manual flush lacks per-repeat evidence. Flush after imports evicts large native-runtime pages and contradicts the selected protocol. Lower Cached is a required campaign criterion even where flush succeeded but background activity obscures it.

## D5 — Build parity and environment pins (FR-007; hole 7)

**Decision**: Parent and every Python child derive from one Python 3.11 real base binary. Record realpath SHA-256, full sys.version, implementation identity, sysconfig platform, architecture and executable; equality applies to digest/build/architecture, not merely '3.11'. Use importlib.metadata.version for every adapter's distribution and freeze all direct/transitive versions, native libraries, converter and CHILmesh revision plus lock digest. Install/import-test all four environments before freezing; a failure blocks the campaign, never silently changes interpreter. Force/record four thread variables to 1 at launch and in child. Proposed exact lock versions are generated only after later compatibility testing, not invented now.

**Rationale**: adcircpy lacks reliable __version__; per-tool venv paths alone prove neither interpreter nor dependency parity. No network in this dispatch, so compatibility is not claimed.

**Alternatives considered**: Parent-only machine.python misses child drift; accepting same minor version permits different builds; forcing current 3.14 contradicts the resolved 3.11 choice.

## D6 — Schema, aggregation, rights and claims (hole 12)

**Decision**: Schema v2 records one operation/backend/attempt; separate status, equivalence, verification, measurement validity, stability and publication rights. Explicit legacy v1 reader remains for existing workflows but cannot grant fair-campaign eligibility. Append attempts, hash raw evidence, prevent overwrite. Group by frozen campaign/cell/series and exact input/contract/env identities; no mixing successes across series. For each metric expose raw values, count, median/min/max/range/relative spread, with >10% unstable. Cell unstable if any of three metrics unstable. Encode infinite spread as null plus spread_state=infinite to keep strict JSON. Replacement series requires rationale and links to original failed/invalid series; no cherry-picking.

**Rationale**: v1 cannot represent the required evidence and its first-match rendering silently discards repeat variation. All-zero growth has zero spread but prohibits division by zero. Rights approval is separate from correctness.

**Alternatives considered**: Optional v1 fields cannot establish a mandatory evidence contract. A bare rights allow flag lacks a ruling citation and mesh version. Omitting a failed STOFS tool distorts coverage. Free-form advantage prose without structured evidence cannot be automatically certified.

## D7 — Independent NetCDF conversion (FR-004/015)

**Decision**: A small harness-owned standard-library + netCDF4 script, `benchmarks/competitors/scripts/convert_fair_netcdf.py`, writes NetCDF from the independent reference parser's validated output, never from CHILmesh code or a timed reader. Converter contract revision is `fair-netcdf-v1`; record the script SHA-256/source revision, exact netCDF4/native-library versions and environment lock hash. Exact dependency versions are frozen after compatibility testing, not invented here. Proposed invocation: `<base-python> benchmarks/competitors/scripts/convert_fair_netcdf.py --reference <reference-output> --reader <xarray|thalassa> --output <new-file.nc>`. The reference output supplies the validated D1 node/element records and source digest, not merely their checksums. Record full argv/options, source and converted hashes, ID maps and separate conversion/read/end-to-end costs before any repeat flush.

**Dimension/variable contract**: Produce reader-specific files: xarray requires `node=N`, `nele=E` (`adapters/xarray_read.py:25-26`); Thalassa requires `node=N`, `element=E` (`adapters/thalassa_read.py:26`). Let `elem_dim` denote that reader's element dimension, with `max_face_nodes` equal to the largest source arity. Both files contain signed int64 `node_id(node)`, binary64 `x(node)`, `y(node)`, `depth(node)`, signed int64 `element_id(elem_dim)`, `arity(elem_dim)` and `connectivity(elem_dim,max_face_nodes)`. Connectivity contains original 1-based file node IDs in original winding order, padded with an explicit zero fill value beyond arity; IDs can be noncontiguous and are never inferred from row offsets. Rows follow ascending original IDs. Record dimension/variable names, dtypes, fill/index conventions and source-to-row maps in ConversionReference. The current adapters only count dimensions; their future verification paths must expose all these variables after eager materialization. Prove the pinned Thalassa API preserves this contract in preflight; otherwise mark the conversion incomplete, without assuming an untested API works.

**Verification and independence**: Reopen the written file outside measurement using netCDF4, reconstruct canonical node/element records from the declared variables/maps, and match counts and both D1 checksums against the independent reference. Each reader's actual materialized output is also verified after END. A CHILmesh-written file is not accepted for this campaign: conversion is independent harness code written against reference-parser output. Dimension-only matches or unverified mappings are incomplete and cannot support claims. Conversion remains separate from native fort.14 comparisons.

**Rationale**: A named independent writer avoids having the subject tool manufacture a competitor's input. Reader-specific dimensions satisfy the two existing count lookups without assuming that their element dimension names are interchangeable; explicit content variables enable verification beyond counts.

## Design completion and runtime prerequisites

No unanswered design question remains. Later runtime prerequisites: exact compatible locks, Valence availability/version evidence, access to cache reset and cgroup counters, idle host/allocation identity, proven adapter APIs/ID maps and authoritative rights rulings. Missing prerequisites yield explicit outcomes, never assumed measurements or publication permission. WNAT_Hagen remains internally inventoried and publicly excluded while its ruling is absent. No external research or measurement is needed to author these tasks.
