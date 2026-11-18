<!-- provenance: author=domattioli model=unknown effort=medium date=2026-11-10 skill=speckit-plan repo=CHILmesh session=dispatch-sol-78a08f; manual stamp, exact provider model ID unavailable; transport rung=sol -->
# Data model and internal contracts

All contracts below are proposed v2 interfaces, not currently implemented commands. Persist strict UTF-8 JSON: no NaN/Infinity, bytes as integer counts, durations in seconds, UTC timestamps plus monotonic timestamps for ordering. Unknown quantities are null with an explicit missing flag/reason, never fabricated zero. Reference/sampler/canonicalization revisions are mandatory identity inputs.

## Campaign

Fields: schema_version=2, campaign_id, manifest_sha256, created_at, chilmesh_git_sha, host_id, host_manifest, runtime_manifest, environment_manifests, sampler_policy, load_policy, cache_policy, canonicalization_revision, ordered_cells, schedule, limits_by_mesh_class, repeats_by_class, previous_campaign_id (optional).

Host manifest: stable campaign host identifier (hash of boot identity plus container identity, without credentials), CPU model, logical count, affinity/available cores, physical RAM bytes, kernel/OS, cgroup path/mount and CPU quota, cpuset and memory limits. Preserve raw values and evidence; do not equate advertised RAM to cgroup allowance. Allocation/host changes invalidate continuation under the same manifest.

Runtime manifest: full sys.version, sys.implementation fields, sysconfig platform, architecture, executable path/realpath, executable SHA-256. Every Python child and parent must match the frozen Python 3.11 build. Non-Python runtime is recorded distinctly and never passes Python equivalence.

Environment manifest: tool distribution name and exact importlib.metadata version, environment identifier, lock path/hash, complete direct/transitive distributions, native runtime/library versions, converter identity/options, executable and source revision where appropriate. Unknown/missing pins block eligibility. Actual install/import compatibility is a later preflight result, not a design assumption.

Policies: 0.010 s RSS cadence, child-only scope, byte units, fixed sampler version; >0.020 s sample-gap invalidation; 0.100 s load observations and 1 s pre-run window; four thread environment variables fixed to 1; cache reset before child launch. Limits: small 600 s watchdog/8589934592 RSS bytes, STOFS 1800 s/51539607552 bytes; counts >=5 small, >=3 STOFS. Watchdog begins at process launch, operation elapsed excludes other work.

Freeze before first measurement; changes create a new manifest/campaign rather than mutating evidence. Serial deterministic order: declared mesh order, operation/backend order, repeat index, fixed tool order. Include unsupported cells in matrix even though no process launches. Replacement series names superseded series and a rationale, and keeps all attempts visible.

## MeshReference and ConversionReference

MeshReference: reference_id, name, aliases, exact version, Valence asset/provenance/served digest, raw_file_sha256, raw_size_bytes, header_node_count, header_element_count, canonical_node_sha256, canonical_element_sha256, canonicalization_revision, parser_revision, precision_token_count_gt17, malformed/nonfinite/duplicate/missing evidence, verification_state, evidence hashes, rights_ruling_id.

Header and canonical counts derive from the same raw bytes identified by the digest. Canonical records use original 1-based file identifiers, not row positions; node `id x y depth`, element `id arity n1 n2 ...`, single spaces, ascending numeric IDs, LF between records, no terminal LF; repr(float(token)) for each binary64 numeric value. Keep winding and signed zero. Reject duplicate IDs, invalid references and unsupported arity; mixed triangles/quads must pass. >17 significant digits makes verification incomplete. Bounded external sorting may create temporary chunks outside timed regions; parser imports no competitor or CHILmesh parser.

ConversionReference (research.md D7, `fair-netcdf-v1`): source_reference_id, source_sha256, converted_sha256, converter script revision/SHA-256, netCDF4/native-library versions and environment lock hash, full invocation/options, converted format, reader-specific dimension/variable/dtype/fill/index contract, source-to-row maps, source-to-converted ID map evidence, source-content verification, conversion_s and read_s as separate costs, end_to_end_s with explicit definition. D7 conversion uses only independent reference-parser output, never CHILmesh. xarray files use node/nele and Thalassa files node/element; both preserve node_id/x/y/depth and element_id/arity/connectivity as specified in D7. Reopen with netCDF4 and match counts and both canonical hashes against the independent reference; each reader must also verify its actual materialized output. Converted content must match reference, not just dimensions. Missing fields or unverifiable mapping make conversion invalid/incomplete. NetCDF rows are conversion / NetCDF only; no native fort.14 ratio.

## OperationContract and ComparisonCell

OperationContract: contract_id/hash, operation, backend_selected, actual_backend/work evidence, input_format/digests, auxiliary/query/field hashes and options, prerequisites, measured_work, required_output, materialization condition, output_verification, equivalence_class, reason/evidence, tools covered.

Operations: parse; adjacency-cpp and adjacency-python; locate-tree-build and locate-query as distinct cells; render; stations; write; full-init. Native parse equivalence is permitted only after evidence for CHILmesh/adcircpy. Other CHILmesh operations are CHILmesh-only absent proven competitor work/output equivalence. Backend labels never imply C++ acceleration; record actual implementation and mark unavailable backend explicitly. Thalassa/xarray native cells and both W1-only tools' W2 cells are not-comparable with reason; conversion contracts are distinct. Import/install context is outside the comparison operation matrix.

Cell key: campaign_id + mesh_reference_id + tool + environment_id + contract_id + backend + series_id. Cell stores ordered attempt IDs, required repeat count, classification, complete flag, metric summaries, overall stable flag, rights/claim eligibility and all rejection reasons. A missing attempt is a coverage failure, not an inferred not-run success.

## RepeatRecord v2

Required fields:

- Identity: schema_version=2, campaign_id/manifest hash, cell_id, series_id, repeat_id, repeat_index, execution_order, reference_id, contract_id, tool/environment/version, CHILmesh revision, exact argv/options and thread environment, host identity/allocation.
- Orthogonal decisions: outcome (`ok|invalid|failed|not-run|not-comparable`), reason, equivalence (`like-for-like|CHILmesh-only|not-comparable`), verification (`verified|failed|incomplete|not-applicable`), measurement_validity with per-time/memory/cache/load/runtime flags and reasons, publication eligibility. Stability is computed on the series, never copied from a single repeat.
- Child protocol: pid, process_group_id, READY/START/END/ACK timestamps, child clock metadata, receipt timestamps, deferred new_modules, materialization/output contract result.
- Measurement: operation_elapsed_s, process_wall_s, rss_baseline_bytes, rss_operation_peak_bytes, rss_whole_process_peak_bytes (launch..END), rss_growth_bytes, last_rss_bytes, wait4_lifetime_maxrss_bytes separately labelled; policy version, interval, scope, coverage/gaps/count, missing flags, raw sample sidecar path/hash.
- Samples: each has monotonic timestamp, child pid, VmRSS bytes or missing reason, descendant PID list. Baseline is final pre-START sample. Region peak uses only START..END samples; whole-process sampled peak only launch..END. Growth=max(0,region_peak-baseline). Missing baseline/region sample/boundary or any descendants invalidates memory.
- Cache: ordered sync/drop_caches operations and exit codes, Cached before/after (raw KiB and normalized bytes), UTC and monotonic timestamps, reset/launch order and success. A user's --cache value never supplies evidence. Unknown cache is invalid for campaign comparison.
- Load: pre-run and throughout-run proc/stat/loadavg/diskstats plus process CPU/I/O/accounting snapshots or hashed sidecar, allowlist, detector policy and verdict. Overlap, foreign workload, allocation change or unknown evidence invalidates the affected run.
- Verification: actual output n_nodes/n_elems, node/element hashes, reference expected values, adapter mapping provenance, normalization revision and precise mismatch/missing reasons. Counts must come from materialized output; no default zeros and no re-reading source as a substitute for tool output.
- Termination: launch/exit times, returncode/signal, configured limits, trigger (`timeout|rss-limit|kernel-oom|signal-unknown|nonzero-exit|launch-error|none`), cgroup identity and memory.events counters/delta, memory.peak optional and scoped, kill/reap diagnostics, stdout/stderr artifact hashes and bounded tails.

Successful process exit is not enough for outcome ok: evidence failures change outcome to invalid, retaining process return code. Resource breach or actual tool failure yields failed, with partial values and missing flags; it is never not-run. Unavailable tool/input/preflight yields not-run with reason; unsupported operation yields not-comparable. Aborted active children become failed; planned unstarted cells become not-run. A single failure or invalid attempt makes that series incomplete for claims. All records are append-only atomic creations; conflicting IDs/overwrite attempts fail. Hash files after finalization and preserve a manifest of raw record hashes for reproducible offline review.

### Legacy v1 compatibility

Dispatch read/validate by schema version, rejecting unsupported versions explicitly. Preserve legacy v1 writer behavior for W1 and old CLI paths; fair campaigns explicitly use v2 creation. Legacy readers may display v1 as historical, never promote absent evidence to v2 defaults. Do not rewrite previous records. Sample/reference/manifest/aggregate/claim files reside in distinct subdirectories so existing read_records globs do not confuse artifact types.

## Aggregate and eligibility

Each named metric has unit, raw repeat ID/value pairs, included/excluded attempt IDs and reasons, n, median, min, max, absolute_spread=max-min, relative_spread_pct=100*(max-min)/median. All-zero series => 0%; zero median and positive spread => null relative_spread_pct plus spread_state=infinite, unstable=true. Exactly 10% passes; >10% fails. Compute summaries for operation_elapsed_s, rss_whole_process_peak_bytes and rss_growth_bytes separately. For incomplete series, observed successful-value summaries may be displayed as descriptive/incomplete with all failures, never eligible. Do not pool environments/contracts/inputs/backends/series. Cell stable requires all three comparison metrics stable; whole-process peak and operation growth remain different named metrics.

## RightsRuling and Claim

RightsRuling: mesh name/exact version, permitting publication-results decision, authoritative Valence issue/PR citation, ruling date, evidence reference/hash and scope. Missing, ambiguous, denied or mismatched rulings exclude. A boolean allow flag alone fails. WNAT_Hagen has no ruling and must remain excluded until one exists; other meshes are never presumed cleared.

Claim proposal: claim_id, sentence template, supporting cell/series IDs, mesh/version, operation/contract, named metric/unit, exact tool versions, repeat counts, medians/spreads, numerator/denominator cell IDs, direction, formula, resource-condition reference. Gate result: accepted/rejected, all failed predicates, source record hashes and computed value. For 'A uses X% less than B', X=100*(B_median-A_median)/B_median with B explicitly denominator; require B>0. Ratio B/A requires A>0. Same metric/unit only, no whole-RSS-to-growth division. Both cells must be like-for-like, verified, complete, stable across every comparison metric, valid cache/load/runtime/measurement and rights-cleared. Conversion/native crossings and CHILmesh-only cells reject. Reject unsupported generalization; render only bounded generated claim templates, not arbitrary certified prose. Every accepted sentence cites cells and scope. FR-020 historical memory sentence is not a template or expected result.

## Internal CLI and child contracts

Proposed additions preserve existing w1/w2/report behavior via explicit fair-campaign mode:

- `python -m benchmarks.competitors fair-preflight --manifest PATH --out DIR`: offline-capable validation, runtime/host/cache/cgroup capability results; no measured child, no downloads by default. Exact-version rights are needed for publication, not internal coverage.
- `python -m benchmarks.competitors fair-run --manifest PATH --out DIR`: requires CHILMESH_RUN_BENCH=1 and passing execution preflight; Valence-only acquisition, frozen scheduling, append-only records. Unsupported/missing prerequisites still produce full-matrix outcomes; failed preflight returns nonzero after writing its evidence/outcomes.
- `python -m benchmarks.competitors fair-report --campaign DIR [--published --rights PATH] [--claims PATH]`: offline recomputation, verifies raw hashes and coverage, rights-filtered publication. Internal report clearly labelled internal. Published mode defaults deny without rights. Rejected proposed claims return nonzero and produce rejection reasons; unrelated eligible tables may be retained for review.

Exit codes: 0 command completed with valid artifacts (a recorded tool failure alone is not a controller error); 2 malformed args/missing opt-in; 1 preflight/coverage/integrity/claim/publication-gate failure. No output path overwrites existing raw attempts.

Dedicated inherited control pipe uses newline JSON protocol_version, repeat_id, event, monotonic_s, clock identity. READY waits for GO; START/END bracket operation; END waits for ACK; only then verification/result JSON is written to an untruncated designated result channel/artifact. Validate sequence, IDs, clock and timestamps; missing or duplicate boundaries invalidate evidence. Stdout/stderr are not trusted as a complete result transport.

## Report contract

Check internal full-matrix coverage before rights filtering. Public output then excludes rights-ineligible numeric data before aggregation/rendering, including failure timing, peaks, charts, summary totals and claim diagnostics. Exclusion list contains names/reasons only. Group separately by mesh/version, operation/backend, format/equivalence and series, never merely workflow/track. Preserve permitted STOFS failures and not-run/not-comparable tools with reasons, limits and available diagnostics. Show native like-for-like, conversion / NetCDF, CHILmesh-only and not-comparable sections; visible caveats for sampling, import-inclusive whole-process RSS, growth as upper bound rather than total mesh allocation, failures, incomplete series and instability. Each displayed number links to hashed raw records and deterministic arithmetic inputs. Rights citations accompany admitted meshes.
