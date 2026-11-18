<!-- provenance: author=domattioli model=claude-fable-5-1 effort=medium date=2026-11-10 skill=speckit-analyze repo=CHILmesh session=dispatch-fable-9655eb; manual stamp, doc_provenance.py not reachable from this sandbox -->
# Specification Analysis Report: 007-fair-comparison

Date: 2026-11-10. Analyze cycle 1. Read-only pass over spec.md, plan.md, research.md, data-model.md, quickstart.md, tasks.md, checklists/requirements.md, `.specify/memory/constitution.md` and the current code under `benchmarks/competitors/`. No artifact other than this report was edited. No benchmark, test, install or network action was run.

Binding operator decisions checked against every artifact: budgets 600 s / 8 GiB (small) and 1800 s / 48 GiB (STOFS); checksum `repr(float(token))` binary64; one Python 3.11 base binary; cache flush before child launch; OMP/OPENBLAS/MKL/NUMEXPR thread variables = 1. All five appear consistently in spec.md:14,173-174, plan.md:14,21, research.md D1/D4/D5, data-model.md:12,16,24 and tasks.md T004/T022/T028/T032. No artifact contradicts a decision. The findings below are about stale wording, ordering, code mismatches and underspecified items.

Constitution: no MUST principle is violated. Principle III (test-first) is the only one with a finding (O1, MEDIUM, ordering inside tasks.md rather than a missing gate).

## Findings

| ID | Category | Severity | Location(s) | Summary | Recommendation |
|----|----------|----------|-------------|---------|----------------|
| S1 | Stale wording | HIGH | spec.md:6, spec.md:171, spec.md:178; checklists/requirements.md:18-19,29,46 | Spec still says an obsolete draft status, heads the resolved decisions an obsolete open-decisions heading, and its Out of Scope paragraph says "Runtime limits and checksum policy remain for clarify" and "does not ... advance clarify, plan, tasks, analyze or implement phases". The checklist leaves three items unchecked for the same two now-resolved decisions. plan.md:35 already has to explain this away. | Set Status to "Clarified; plan and tasks authored 2026-11-10". Rename section to "Resolved decisions (clarify 2026-11-10)". Delete the last two sentences of Out of Scope. Check the three checklist boxes with the decision citation. |
| S2 | Contradiction | HIGH | spec.md:119 (FR-009) vs spec.md:14, research.md:44-46 (D4), plan.md:61 | FR-009 text still says the cold reset happens "after preparation and before the operation". Operator decision and D4 put the flush before child launch, so imports and in-child preparation run after the flush. research.md:46 openly says D4 "supersedes the older literal FR-009 ordering". A requirement line that its own plan supersedes is an audit hazard. | Rewrite FR-009 to: "documented, successful cold file-cache reset after all hashing, reference parsing and conversion, and before the child process launches; imports run inside the child after the flush and stay outside the timed region." |
| S3 | Stale wording | MEDIUM | spec.md:25 (the obsolete cache-order question), spec.md:26 (the obsolete interpreter-version question) | Two clarification bullets still defer to questions that spec.md:14 records as answered. | Replace both tails with the decision: "resolved: flush before child launch" and "resolved: Python 3.11". |
| S4 | Inconsistency | MEDIUM | spec.md:109 (FR-002) vs spec.md:30, plan.md:22, data-model.md:32 | FR-002's operation inventory lists read, adjacency/layers, point location, station sampling and rendering. The clarification bullet, plan and data model also inventory `write` and `full-init`, and split locate into tree-build and query. FR-002 is the normative list. | Extend FR-002's list to: parse, adjacency (per backend), locate tree-build, locate query, render, stations, write, full-init. |
| S5 | Underspecification | MEDIUM | spec.md:111 (FR-004), data-model.md:26, tasks.md:47 (T018) | FR-004 and ConversionReference require "converter/version/options" but no artifact names a converter or says who writes the NetCDF. T018 puts "conversion preparation" in campaign.py with no converter identity. The xarray adapter expects dimension names `node` and `nele` (xarray_read.py:25-26) and Thalassa expects `node`/`element` (thalassa_read.py:26); the converted file's variable and dimension contract is nowhere stated. If CHILmesh writes the NetCDF, the subject tool produces a competitor's input, which D1's independence rationale should at least address. | Add a research decision D7: name the converter (tool, version, invocation), the NetCDF variable/dimension contract each reader needs, and state that converted content is verified against the independent reference (data-model.md:26) so a CHILmesh-written file is acceptable or not. Reference D7 from T018. |
| O1 | Ordering | MEDIUM | tasks.md:6, tasks.md:72 (T025), tasks.md:94 (T039), constitution III | tasks.md:6 promises "For each phase, write tests and demonstrate their intended failures before implementing." T025 is a test task placed after the US2 implementation tasks T026–T033, and T039 is "Add and pass" an acceptance test after all US3 implementation. Both are tests written after the code they test. | Either move T025 into the US2 "Tests first" block (before T026) and T039 into the US3 "Tests first" block, or amend tasks.md:6 to say that integration/acceptance tests T025 and T039 are the exception and explain why. |
| C1 | Code mismatch | MEDIUM | tasks.md:68 (T030) vs tests/benchmarks | T030 says "Update historical 'OOM suspected' test expectation in `tests/benchmarks/test_measure.py`". A search of tests/benchmarks for "OOM suspected" finds no match. The string exists only in measure.py:66. The only reason-string assertion is `timeout after 0.1 s` in tests/benchmarks/test_w1_fixes.py:26 (an uncommitted file). | Change T030 to "add a test that SIGKILL without an oom_kill delta is `signal-unknown`, and confirm no existing test asserts 'OOM suspected'". |
| C2 | Underspecification | MEDIUM | plan.md:62, tasks.md:25 (T006), tasks.md:44 (T015), w2/run.py:35 vs w2/run.py:156 | Adapters must import the new control protocol, but competitor children launch with `env=dict(os.environ)` and no PYTHONPATH (run.py:155-156), while the CHILmesh worker gets REPO_ROOT on PYTHONPATH (run.py:35). `benchmarks/competitors/__init__.py` is a one-line docstring today, so a package import is light, but no task says the competitor launch adds the repo to `sys.path`, nor that `protocol.py` must stay free of CHILmesh and NumPy imports so a GPL venv can load it. test_gpl_scan.py:46-54 forbids the harness importing adapters; the reverse direction (adapters importing harness) has no rule. | In T006 state: protocol.py is standard-library only and importable with the repo root on `sys.path`. In T015 or T033 state: competitor launches set PYTHONPATH to REPO_ROOT. Add a test_gpl_scan case that adapters import only `benchmarks.competitors.protocol`. |
| C3 | Stale wording | LOW | w2/tools.py:44, w2/run.py:147-149, tasks.md:46 (T017) | `XARRAY_FORT14_REASON` ("xarray reads NetCDF...") is the reason emitted for any tool whose input is not fort.14. After T017 marks Thalassa NetCDF-only, Thalassa cells would carry a reason that names xarray. | T017: rename to `NETCDF_ONLY_REASON` and format the tool name into it. |
| R1 | Stale wording | LOW | research.md:9, research.md:14, research.md:15 | Line citations overrun the files: records.py cited as 102-136 (file ends at 127), valence.py as 91-115 (ends 103), report/render.py as 35-88 (ends 72). | Correct to records.py:101-126, valence.py:91-102, render.py:35-71. |
| Q1 | Inconsistency | LOW | quickstart.md:11-15 vs tasks.md:57 (T021) | The quickstart's test command list omits `tests/benchmarks/test_measure_lifecycle.py`, which T021 creates. The final `python -m pytest tests/benchmarks` line still covers it. | Add the file to the second command line in quickstart.md. |
| Q2 | Ambiguity | LOW | spec.md:97 vs spec.md:125 (FR-012) vs data-model.md:50 | Edge case says one failed repeat makes the cell "incomplete for claims until a new, fully disclosed complete series meets the repetition rule". FR-012 and data-model.md:50 agree, but FR-013 (spec.md:126) says "No failed ... value may be ... included in a numerical summary" while data-model.md:58 allows "descriptive/incomplete" summaries of the successful values of an incomplete series. Both can be true (descriptive, never eligible) but the spec never says descriptive summaries are permitted. | Add one sentence to FR-013: "An incomplete series may show a descriptive summary of its successful repeats, labelled incomplete and never eligible." |

No duplicate requirements found. No unresolved placeholder tokens (TODO, TKTK, ???) found in the six artifacts.

## Coverage summary

Mapping from tasks.md:121 was checked requirement by requirement against task text. Every FR and every buildable SC has at least one task.

| Requirement | Has task? | Task IDs | Notes |
|---|---|---|---|
| FR-001 | yes | T011, T017, T019, T043 | six-tool inventory; tools.py:32-42 has three today |
| FR-002 | yes | T011, T017 | see S4 for the inventory list |
| FR-003 | yes | T016, T017, T019 | backend identity in T016 |
| FR-004 | yes | T010, T012, T015, T018, T042 | see S5 for converter identity |
| FR-005 | yes | T011, T017, T042 | |
| FR-006 | yes | T023, T029, T033 | host lock, load, allocation |
| FR-007 | yes | T004, T007, T032 | 3.11 digest parity |
| FR-008 | yes | T003, T006, T015, T016 | |
| FR-009 | yes | T022, T028, T033 | see S2 for FR text |
| FR-010 | yes | T020, T026 | |
| FR-011 | yes | T020, T026, T027 | |
| FR-012 | yes | T004, T024, T031, T033 | |
| FR-013 | yes | T024, T031, T025 | |
| FR-014 | yes | T008, T013 | |
| FR-015 | yes | T009, T014, T019 | |
| FR-016 | yes | T002, T005, T021, T027, T030, T043 | |
| FR-017 | yes | T002, T005, T026, T036, T039 | |
| FR-018 | yes | T034, T040, T042 | |
| FR-019 | yes | T035, T041, T042 | |
| FR-020 | yes | T035, T041, T046 | |
| SC-001 | yes | T011, T043, T039 | |
| SC-002 | yes | T024, T031, T039 | |
| SC-003 | yes | T009, T014, T039 | |
| SC-004 | yes | T004, T023, T025 | |
| SC-005 | yes | T034, T036, T042 | |
| SC-006 | yes | T025, T039, T048 | |

Unmapped tasks: none. T001 (fixtures), T046–T048 (docs, regression run, handoff) are infrastructure and map to constitution III and VIII rather than an FR.

Code-vs-task check: every file a task names exists today or is declared new in plan.md:51. Verified existing: measure.py, records.py, machine.py, valence.py, cli.py, report/render.py, w2/run.py, w2/measures.py, w2/tools.py, adapters/{adcircpy_read,thalassa_read,xarray_read}.py, envs/, tests/benchmarks/{conftest.py,_grids.py,test_records.py,test_measure.py,test_w2_wiring.py,test_valence.py,test_report_render.py,test_bench_cli.py}. Declared new: protocol.py, campaign.py, cache.py, reference.py, verification.py, aggregate.py, claims.py, scripts/run_fair_codespace.sh, envs/README-fair.md, envs/*-fair.lock, and eight new test files. Task claims about current code that were verified true: 20 ms wait4 loop (measure.py:16,47-55); SIGKILL labelled "OOM suspected" (measure.py:65-66); `_mesh_meta` zero defaults (run.py:46); competitor counts never copied into the record (run.py:157-163); one locate call reused for tree-build and query (run.py:113-121, measures.py:272-277); Thalassa marked fort.14 (tools.py:37); `publish_results: true` is the whole rights gate (valence.py:96); parent-only Python version (machine.py:40); `read_records` globs `*.json` (records.py:98); `cell_for` returns the first matching record (records.py:106-113); xarray forces load, Thalassa does not (xarray_read.py:22, thalassa_read.py:22). The one false claim is C1.

## Metrics

- Total requirements: 20 FR + 6 SC = 26
- Total tasks: 48
- Coverage: 26/26 (100%)
- Ambiguity count: 1 (Q2)
- Duplication count: 0
- Critical issues: 0
- High: 2 (S1, S2); Medium: 6 (S3, S4, S5, O1, C1, C2); Low: 4 (C3, R1, Q1, Q2)

## Next actions

No CRITICAL finding. Implementation may proceed after the two HIGH edits, which are wording-only and take minutes:

1. Edit spec.md for S1, S2, S3, S4 and Q2 (status line, section title, Out of Scope tail, FR-009 text, FR-002 list, FR-013 sentence) and check the three checklist boxes.
2. Add research.md D7 (converter identity and NetCDF contract) for S5 and cite it from T018.
3. Edit tasks.md for O1 (move or justify T025/T039), C1 (T030 wording), C2 (protocol import path in T006/T015/T033), C3 (T017 reason constant).
4. Fix the three line citations in research.md (R1) and the quickstart test list (Q1).
5. Rerun this analysis once; expected result zero findings above LOW.

Extension hooks: `.specify/extensions.yml` was not read in this sandbox; no hook output is claimed.

## Remediation

Applied 2026-11-10; findings and source line locations above describe the original review snapshot. Task IDs above now refer to the renumbered backlog. No implementation or benchmark evidence is claimed.

- S1: Set clarified status and resolved-decision heading, removed obsolete Out of Scope sentences and plan workaround, and checked the three decision-related checklist items with citations.
- S2: Aligned FR-009 and D4 on hashing/reference/conversion before flush, flush before launch, and imports in the child outside timing.
- S3: Replaced deferred cache-order/interpreter questions with the resolved before-launch/Python 3.11 decisions.
- S4: Made FR-002 inventory explicit for parse, per-backend adjacency, separate locate build/query, render, stations, write and full-init.
- S5: Added D7 naming the independent harness stdlib+netCDF4 converter, version/invocation, reader-specific dimensions and content variables, and independent-reference verification. Linked T018, plan, data model and quickstart; CHILmesh-written conversion is excluded.
- O1: Moved the US2 integration test (now T025) and US3 acceptance test (now T039) into tests-first blocks, required intended failures before implementation and passing integration checkpoints afterward, and renumbered all task references.
- C1: Replaced the nonexistent historical assertion update with a SIGKILL/no-oom_kill-delta signal-unknown test and an existing-assertion check (now T030).
- C2: Specified stdlib-only protocol imports, competitor PYTHONPATH=REPO_ROOT, and an adapter-to-harness GPL scan case in T006/T010/T015 and the plan.
- C3: Required tool-named NETCDF_ONLY_REASON in T017 for tools.py and run.py.
- R1: Corrected records.py, valence.py and report/render.py line ranges in research evidence.
- Q1: Added test_measure_lifecycle.py to the second quickstart pytest command.
- Q2: Explicitly permitted descriptive successful-repeat summaries of incomplete series in FR-013, labelled incomplete and never eligible.

All five operator decisions remain: small 600 s/8 GiB, STOFS 1800 s/48 GiB; binary64 repr checksum; one Python 3.11 base binary; flush before child launch; single-threaded thread environment. Verification is document/scope validation only; campaign prerequisites remain later work.
