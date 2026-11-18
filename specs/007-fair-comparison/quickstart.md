<!-- provenance: author=domattioli model=unknown effort=medium date=2026-11-10 skill=speckit-plan repo=CHILmesh session=dispatch-sol-78a08f; manual stamp, exact provider model ID unavailable; transport rung=sol -->
# Quickstart: implementing and reviewing fair comparison

The fair CLI and Codespace wrapper are implemented. Offline synthetic checks below validate their contracts; they do not establish real-host acceptance, compatible installations, publication rights or measured advantages. This implementation session ran no real benchmarks, installs, network requests or cache resets.

## 1. Develop using synthetic evidence

Use a disposable local test environment and synthetic tiny meshes; no Valence assets or real competitor installations are needed for default tests. Keep temporary HOME/cache/output under the test temporary directory. Tests must first fail for their intended missing behavior, then pass after implementation.

```bash
python -m pytest tests/benchmarks/test_reference.py tests/benchmarks/test_verification.py tests/benchmarks/test_records.py
python -m pytest tests/benchmarks/test_protocol.py tests/benchmarks/test_measure.py tests/benchmarks/test_measure_lifecycle.py tests/benchmarks/test_campaign.py tests/benchmarks/test_cache.py tests/benchmarks/test_machine_evidence.py
python -m pytest tests/benchmarks/test_fair_adapters.py tests/benchmarks/test_w2_wiring.py tests/benchmarks/test_aggregate.py
python -m pytest tests/benchmarks/test_claims.py tests/benchmarks/test_valence.py tests/benchmarks/test_report_render.py tests/benchmarks/test_bench_cli.py tests/benchmarks/test_fair_codespace.py tests/benchmarks/test_fair_campaign.py
python -m pytest tests/benchmarks
```

Expected: fake procfs/cache/cgroup data makes all tests offline and unprivileged; Linux lifecycle tests use only bounded synthetic child processes, including a descendant and timeout. No real download or cache-reset command is executed by tests. Existing W1 records and report behavior remain compatible; schema-v1 evidence cannot become fair-campaign eligible.

## 2. Prepare the later Codespace campaign

This step belongs to the operator-authorized execution session. Use one 16-core, approximately 62.8 GiB host; record actual kernel, RAM bytes, CPU/cgroup limits and stable host identity. The historical host values in spec.md identify the intended class, not proof that any current machine matches.

Create every venv from the same selected Python 3.11 base binary, including the harness. Export `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1` in the parent environment too. Confirm install/import compatibility for CHILmesh, adcircpy, Thalassa and xarray before freezing dependency locks. Resolve tool versions via importlib.metadata; record all dependencies/native runtime versions, lock hashes, full build/platform and real executable digest per child. Keep GPL imports confined to their subprocess adapters. Do not silently change Python if compatibility fails.

The script is `benchmarks/competitors/scripts/run_fair_codespace.sh`. Its default mode validates/produces preflight evidence without measuring or downloading. It accepts an explicit manifest path; runtime/environment lock paths and tool executables are declared there. Credentials are supplied through the established Valence runtime configuration and never copied into manifests, command evidence or logs.

```bash
CHILMESH_FAIR_PYTHON=/absolute/base/python3.11 bash benchmarks/competitors/scripts/run_fair_codespace.sh --preflight --manifest .scratch/fair/campaign.json --out .scratch/fair/run-001
```

The manifest is authored/frozen using data-model.md: four exact mesh identities/versions from Valence, all six tools and operations, equivalence/output contracts, reference parser revision, deterministic serial schedule, 5 small/3 STOFS repeats, 600 s/8 GiB and 1800 s/48 GiB budgets, 10 ms child-only RSS policy, 100 ms load observation policy, thread variables all 1, and cache-reset evidence requirements. Include independent source references and verified conversion lineage prepared before repeat flush using the harness stdlib+netCDF4 converter and reader-specific contract in research.md D7; no CHILmesh-written conversion is accepted. Freeze exact pins only after compatibility succeeds. A preflight capability check may inspect reset permissions; an actual cache flush happens only in explicit run mode.

Supply `evidence/prepared.json` beside the input manifest: `identities` and `commands` keyed by frozen repeat ID, `references` keyed by manifest mesh alias, and optional `conversions` keyed by repeat ID (legacy tool-keyed fallback for one source). These are the existing `run_fair_campaign` prepared inputs; see the [harness README](../../benchmarks/competitors/README.md#fair-comparison-schema-v2) for the transport and rights formats. The CLI never invents preparation or fetches a fallback. Acquisition/reference/conversion must be prepared through the Valence-only harness APIs before repeats.

Missing interpreter pins, cgroup memory.events, host identity, load observability or reset capability blocks eligible execution and creates a reasoned preflight report plus full planned-cell not-run outcomes. Rights are a separate publication prerequisite. Do not presume any of the three other meshes is cleared; WNAT_Hagen currently has no ruling.

## 3. Execute only when authorized later

```bash
CHILMESH_RUN_BENCH=1 CHILMESH_FAIR_PYTHON=/absolute/base/python3.11 bash benchmarks/competitors/scripts/run_fair_codespace.sh --run --manifest .scratch/fair/campaign.json --out .scratch/fair/run-001
```

Script verifies its single-host lock and invokes fair-run serially. Each attempt: observe load; sync; record Cached before; write 3 to drop_caches; record lower Cached after; launch new session with four thread vars set to 1; import/prepare; baseline/markers; materialize; END acknowledgement; verify output; persist evidence. No file hashing or reference parser runs between flush and parse region. A reset failure yields unknown cache and invalid evidence; no 'cold' flag can override it. The host's existing permissions must support the protocol; no fallback to warm runs.

Outputs under the campaign directory: `evidence/campaign.json`, `evidence/preflight.json`, `evidence/references.json`, `evidence/reference-hashes.json`, `evidence/record-hashes.json`, immutable `attempts/` and `raw/` sample/diagnostic files. The wrapper also retains `evidence-archive.tar.gz` after run mode, including failure finalization. The offline report recomputes aggregates and claim decisions and prints Markdown; it does not overwrite or cache verdicts in raw attempts. Reference metadata, sample/diagnostic files and records have hashes. A timeout/RSS breach kills and reaps the process group, saves available elapsed/last/peak values and missing flags, and retains the failed attempt. A cgroup oom_kill delta is required to name kernel OOM. Every STOFS tool has an explicit outcome even if all executions fail. Controller abort finalizes the active repeat and explains unstarted cells. Do not restart into the same IDs or select only successful attempts; a replacement series requires a new ID and rationale with links to the original series.

## 4. Review offline and test publication gates

```bash
python -m benchmarks.competitors fair-report --campaign .scratch/fair/run-001
python -m benchmarks.competitors fair-report --campaign .scratch/fair/run-001 --published --rights .scratch/fair/rights.json --claims .scratch/fair/claims.json
```

The first report is internal. The second must deny rights by default, require a permitting exact-version Valence issue/PR ruling with date, list WNAT_Hagen's exclusion without its measurements, and retain permitted failures. It checks full internal coverage first, then excludes rights-ineligible data from published tables, summaries, charts and claim diagnostics. Every admitted mesh displays its ruling citation. Native fort.14 and conversion / NetCDF are separate; CHILmesh-only and unsupported cells have reasons, no ratios.

Acceptance review uses saved synthetic evidence to deliberately introduce:

1. Empty output, duplicate IDs, changed coordinates/depth/winding/connectivity, mixed tri/quad errors and >17-digit tokens: reject verification or mark incomplete.
2. Interpreter digest/build drift, missing pins, deferred imports, lazy materialization, descendants, absent baseline/markers/samples, sampling gaps, failed reset and foreign CPU/I/O: reject affected comparisons.
3. Four small/two STOFS repeats, one failed repeat, mixed-series successes, zero denominators and relative spreads 10% versus >10%: exact 10 passes; incomplete/unstable cells and division by zero reject claims.
4. STOFS timeout, sampler RSS limit, true cgroup OOM evidence, unexplained SIGKILL and missing tool: preserve distinct outcomes, limits, diagnostics and missing flags; no zero substitution or omissions.
5. WNAT or another mesh with missing/denied/version-mismatched/bare-boolean rights; unsupported advantage text; native-versus-NetCDF ratio: reject publication/claim and list the failed predicates without leaking excluded numbers.

Independently recompute raw n/median/min/max/range/relative spread for operation time, sampled whole-process RSS and growth. Match all displayed arithmetic and accepted claim percentages to cited record hashes. Whole-process RSS includes startup/imports through END; operation growth is an upper bound affected by allocator retention/file mappings, not total mesh memory. Sub-10-ms spikes remain unobserved. The old 40–60% statement supplies no expected outcome. Equal, worse or universally failing results satisfy honest-reporting acceptance.

## Completion boundary

Complete tests, documentation and script validation; hand off evidence for review. This feature plan does not include moving specs to DomI, implementing a backend, W1/W3, #309, committing, pushing or publishing. The operator handles spec relocation separately.


## Verified offline handoff — 2026-11-11

US3 and polish implementation is present; task checkboxes were not changed.
The dispatch had no Skill tool, so caveman reporting was applied by hand.
Preparation transport was unspecified: the implementation documents the
`evidence/prepared.json` seam above. Per-repeat conversion keys prevent a single
tool-keyed conversion from silently referring to the wrong mesh in a multi-mesh
campaign. All code and test writes stayed under the authorized benchmark paths.

Actual acceptance commands (repository root; no network, installations, real cache
reset or real meshes):

```text
.scratch/venv-cpp/bin/python -m pytest tests/benchmarks -q
673 passed, 6 skipped in 24.62s; exit 0
bash -n benchmarks/competitors/scripts/run_fair_codespace.sh
(no output); exit 0
```

The benchmark suite includes GPL isolation, legacy W1/schema-v1 regressions and
deterministic offline arithmetic. Three Linux-only lifecycle cases and three
unavailable-native-backend cases are skipped on this host. Linux lifecycle and
compatible native-backend execution remain operator prerequisites.
The original 602 passing tests remain in the passing suite; 71 new tests cover
rights/claims/report/CLI/script and saved-campaign acceptance. Synthetic temporary
campaigns retain `attempts/`, reference metadata/hashes and sample/diagnostic files
under pytest temporary directories for independent review.

Test-first and intermediate evidence (commands below all use
`.scratch/venv-cpp/bin/python -m pytest`; paths abbreviated to `tests/benchmarks/`):

| Arguments / stage | Actual output excerpt | Exit |
|---|---|---|
| `test_claims.py test_fair_campaign.py test_fair_codespace.py test_valence.py test_report_render.py test_bench_cli.py -q`, before implementation | 24 failed; missing claims/renderer; abort propagated KeyboardInterrupt | 2 |
| `test_valence.py test_report_render.py test_bench_cli.py test_fair_codespace.py -q --tb=line`, before implementation | 22 failed, 28 passed; absent rights/CLI/script APIs | 1 |
| Same six-file focused set, `-q --tb=short`, first implementation | 1 failed, 79 passed; synthetic adjacency backend needed explicit manifest tuple | 1 |
| Full `tests/benchmarks -q --tb=short`, after fixture correction | 654 passed, 6 skipped | 0 |
| `test_fair_campaign.py -q --tb=line -k saved_reference_and_sidecar`, first/second red runs | 4 failed each; first exposed missing synthetic sidecars, second confirmed absent integrity checks | 1, 1 |
| `test_fair_campaign.py test_report_render.py test_bench_cli.py -q --tb=short`, after integrity checks | 37 passed | 0 |
| `test_fair_campaign.py -q --tb=short`, expanded saved matrix | 18 passed | 0 |
| `test_claims.py -q --tb=line -k exact_versions`, before stricter claim gates | 4 failed; unknown versions, malformed hashes, cross-series claims accepted | 1 |
| `test_claims.py test_fair_codespace.py test_report_render.py -q --tb=short`, after fixes | 46 passed | 0 |
| `test_fair_campaign.py -q --tb=short -k prepared_cell_identity`, red then green | cell identity collision after abort; then 1 passed | 1, 0 |
| Six-file focused set, `-q --tb=short`, after abort fix | 98 passed | 0 |
| `test_fair_campaign.py -q --tb=short -k conversion_transport`, red | per-repeat conversions arrived as None | 1 |
| `test_fair_campaign.py -q --tb=short`, after transport fix | 20 passed | 0 |
| Six-file focused set, `-q --tb=short`, final | 99 passed | 0 |
| Full `tests/benchmarks -q`, final | 673 passed, 6 skipped | 0 |

An initial attempted test-log redirect to `/private/tmp/us3-red.log` was denied
before pytest launched (exit 1); this was an incidental log write, not a required
operation. Tests then ran directly without that redirect. Shell syntax validation
also passed independently. A no-allow snapshot verification returned exit 1 and
listed only authorized changed benchmark files. The final scope check uses the
supplied snapshot verifier with the three allowed paths:

```bash
python3 .scratch/agents-inc-runs/20261111T125745Z-sol-731977/pre_dispatch_snapshot.py verify .scratch/agents-inc-runs/20261111T125745Z-sol-731977/snapshot.json --allow benchmarks/competitors/ --allow tests/benchmarks/ --allow specs/007-fair-comparison/quickstart.md
```

Scope result: `UNCHANGED outside allowed paths`; exit 0. Protected pre-existing
changes were preserved.

### Success-criterion review

| Criterion | Retained synthetic evidence and boundary |
|---|---|
| SC-001 | `test_fair_campaign.py` saves four meshes × six tools × ten operation contracts. STOFS timeout, RSS limit, evidenced OOM, unknown signal, not-run and unsupported outcomes remain. Abort tests retain partial active evidence and preserve cell identity for unstarted repeats. |
| SC-002 | Saved-series tests recompute all three metrics; exactly 10% passes, 11% is unstable, two STOFS observations remain incomplete. Existing aggregate tests cover five/three counts and zero/infinite spread rules. |
| SC-003 | Independent tiny reference parser plus materialized-output verification rejects empty/same-count corruption before saved-series claim gates. Missing conversion proof cannot establish verified conversion eligibility. |
| SC-004 | Cache, load, runtime, memory and timing invalidations stay descriptive/ineligible. Existing shared-controller tests cover native/competitor launches; new tests retain abort sampler evidence. No real-host compatibility or idleness is inferred. |
| SC-005 | Public return objects and Markdown omit excluded measurements, failure tails, raw hashes and free-form claim diagnostics. Only names/reasons remain; admitted exact-version rights citations accompany cells. Native/conversion/CHILmesh-only/unsupported sections remain distinct. |
| SC-006 | Offline hash checks reject record/reference/sidecar tampering before rights filtering. Accepted bounded claims reproduce exact arithmetic from cited medians/repeats; failed predicates reject unstable/incomplete/non-equivalent/invalid/rights-ineligible supports and historical prose. |

### Remaining operator prerequisites

No real-host acceptance, campaign, publication, package installation or rights
ruling was performed. Before authorized execution, establish Python 3.11 common
binary compatibility and full real locks/native/converter pins; prove the pinned
adapter APIs and ID maps; obtain exact Valence assets/versions and independent
source/conversion proofs; freeze prepared identities/commands; verify Linux host
allocation, cgroup OOM attribution, cache-reset permissions/lower-Cached evidence,
load observability and serial host locking. Obtain authoritative exact-version
publication rulings for every admitted mesh. WNAT_Hagen stays excluded while no
ruling exists; no other mesh is presumed cleared. Run the Linux lifecycle tests
and the complete acceptance commands again on the intended Codespace before any
operator-authorized campaign. Existing placeholder locks deliberately block
execution. No backend/W1/W3/#309 changes, spec relocation, commits or pushes are
part of this handoff.
