caveman: ultra
GRILL: none
[verified] No Skill tool available in dispatch tool inventory; caveman rules applied by hand, no real Skill invocation.
[verified] Evidence paths below: S = specs/006-competitor-benchmarks/spec.md; P = specs/006-competitor-benchmarks/plan.md; T = specs/006-competitor-benchmarks/tasks.md; Q = specs/006-competitor-benchmarks/checklists/requirements.md; C = .specify/memory/constitution.md. Ranges are inclusive.

## Findings

| ID | Severity | Location | Summary | Recommended artifact edit |
|---|---|---|---|---|
| F01 | MEDIUM | S:75,89; P:54,76; T:17,33 | [inferred] FR-008 restricts meshes, so D1 is compatible with its literal scope; clarification's unqualified “no direct NOAA download path” and plan's “No ... NOAA fallback” leave auxiliary-data scope ambiguous. | Amend spec clarification and plan Valence paragraph: mesh downloads only via Valence; NOAA station/output data explicitly exempt, with D1 provenance. |
| F02 | MEDIUM | S:76,120; P:63; T:13,58 | [verified] Clarification requires an Octave attempt and no assumed MATLAB seat; assumption still allows an untested Octave run and mentions MATLAB availability. | Replace assumption with clarified Octave attempt/failure rule; retain separate oceanmesh execution requirement. |
| F03 | HIGH | S:110; P:26,69; T:75 | [inferred] SC-003 requires a second person to reproduce any published result within 10%; T043/plan verify one cell, without an independent-person requirement. | Expand validation tasks to the specified repeatability coverage and independent reproducer, or explicitly obtain a narrower SC-003. Define the STOFS N=1 exception consistently. |
| F04 | HIGH | C:118,138; P:30,86; T:11-77 | [inferred] LESSON-CANDIDATE: plan leaves Medium governance unresolved while declaring no violations; tasks contain no classification/design-issue/feedback-window gate. Constitution alignment cannot be certified. | Add a plan classification decision and task gate before implementation; if Medium, record design issue, tradeoffs and 3-day window. No constitution edit proposed. |
| F05 | MEDIUM | S:68,84; P:26; T:48,73 | [inferred] Spec requires cold-cache comparisons and flagged contaminated runs; plan allows warm/unknown records without specifying publication eligibility or repeat enforcement. | Define cold-cache verification, repeat/flag policy, and comparable-result eligibility in plan and T027/T041. |
| F06 | MEDIUM | S:89-90; P:30,50,54; T:11,29-30 | [inferred] Plan promises per-mesh CRS and spec requires mesh version/provenance; enumerated run-record mesh fields contain only name/hash/counts. No task explicitly persists mesh version or CRS in each record. | Add mesh version, provenance reference and CRS/unit fields to plan schema and T010/T013/T014; preserve rights exclusions. |
| F07 | MEDIUM | T:5,11-15 | [verified] T002-T005 marked parallel “different files,” but all write research.md. | Remove parallel markers or specify separate owned fragments followed by a serial merge task. |
| F08 | LOW | S:95; P:68-69; T:77 | [verified] FR-014 has no explicit task mapping; final T045 moves specs to DomI but has no spec item or plan-phase counterpart. | Add W3 exclusion verification to closure tasks and add the final relocation step to plan Phase 4. Preserve T045 as last task. |
| F09 | MEDIUM | P:17; src/chilmesh/CHILmesh.py:420-422,971-974,1033-1043 | [verified] R6's “lazy build” is inaccurate: constructor calls the index builder when enabled; find_element queries an existing tree. | Correct R6 and specify timing boundaries for construction versus queries in plan/T019, avoiding double counting. |
| F10 | MEDIUM | C:93-96,198; P:30; T:55-57 | [inferred] LESSON-CANDIDATE: synthetic metric/displacement tests lack explicit tri/quad/mixed coverage and bit-exact boundary comparison required by constitution. | Extend T029/T030 and constitution-check paragraph with all three element combinations and np.array_equal boundary checks; document applicability of parity principle C:100-107. |
| F11 | LOW | S:65,101,108; P:8,50; T:26 | [verified] Spec says “not run,” plan uses “not-run”; spec permits not-applicable metrics but SC-001 only accepts values or failure reasons, while plan uses null plus na=true. | Define display labels versus machine status and explicitly allow justified not-applicable cells in SC-001/schema acceptance. |

## FR/SC coverage map

[inferred] “Covered” means a named implementation/validation task exists, not completion. All tasks remain unchecked (T:11-77). Partial denotes a listed finding; boundary-only denotes no explicit task.

| Spec item and evidence | Task evidence | Coverage |
|---|---|---|
| FR-001 S:82 | T028/T036; T:49-51,62,81 | Covered: W2 publication before W1 |
| FR-002 S:83 | T011/T020-T024/T027; T:27,41-48 | Covered |
| FR-003 S:84 | T020/T027/T041; T:41,48,73 | Partial F05 |
| FR-004 S:85 | T033/T036; T:59,62 | Covered |
| FR-005 S:86; S:39 | T004/T030-T035/T011/T013; T:14,27,29,56-61 | Covered via shared metrics and measurement |
| FR-006 S:87 | T029/T031/T033; T:55,57,59 | Covered; constitution refinement F10 |
| FR-007 S:88 | T009/T016/T020/T032/T042; T:22,32,41,58,74 | Covered |
| FR-008 S:89 | T001/T012/T014; T:11,28,30 | Partial F01/F06 |
| FR-009 S:90 | T010/T013/T014/T041/T043; T:26,29-30,73,75 | Partial F06 |
| FR-010 S:91 | T010-T015/T024/T032; T:26-31,45,58 | Covered |
| FR-011 S:92 | T028/T044; T:49,76 | Covered |
| FR-012 S:93 | T037-T040; T:66-69 | Covered |
| FR-013 S:94 | T018/T025; T:39,46 | Covered |
| FR-014 S:95 | No explicit task; workflow task boundaries T:37,53,64,71 | Boundary-only F08 |
| SC-001 S:108 | T010/T015/T027/T028; T:26,31,48-51 | Partial F11 |
| SC-002 S:109 | T029/T035/T036; T:55,61-62 | Covered |
| SC-003 S:110 | T043; T:75 | Partial F03 |
| SC-004 S:111 | T018/T025; T:39,46 | Covered |
| SC-005 S:112 | T008/T016/T042; T:21,32,74 | Covered; licensing interpretation remains open P:46,86 |
| SC-006 S:113 | T037-T040; T:66-69 | Covered |

## Reverse trace and requested checks

| Tasks (source lines) | Spec/plan trace |
|---|---|
| T001-T005 (T:11-15) | [verified] Research P:56-57; FR-005/008 S:86,89 |
| T006 (T:16) | [verified] D1 P:76 |
| T007-T009 (T:20-22) | [verified] Structure/environment P:26,34-46 |
| T010-T016 (T:26-32) | [verified] Core schema, measurement, Valence, renderer, GPL P:30,50-60 |
| T017 (T:33) | [verified] D1 P:76 |
| T018-T028 (T:39-49) | [verified] W2 P:59-60; FR-001/002/003/011/013 S:82-84,92,94 |
| T029-T036 (T:55-62) | [verified] W1 P:62-63; D2 P:77; FR-004/005/006 S:85-87 |
| T037-T040 (T:66-69) | [verified] Doc fixes P:65-66; D3 P:78; FR-012 S:93 |
| T041 (T:73) | [verified] Portability and recorded execution S:74,90; plan README P:37 |
| T042-T044 (T:74-76) | [verified] Closure P:68-69 |
| T045 (T:77) | [verified] No explicit spec/plan trace; dispatch explicitly requires final DomI move; F08 |

[verified] D1 maps to T006/T017 (P:76; T:16,33), D2 to T032 (P:77; T:58), D3 to T038 (P:78; T:67). W3 remains excluded (S:95; T:37-77). W2-before-W1 publication is explicit (T:6,51,81); DomI relocation is last (T:77).
[inferred] Constitution: library/test isolation and stable API intent align (C:18-21,75-78; P:26,30); frozen topology and exact boundary intent align (C:27-31,198; T:55,59). Governance, CRS persistence and mixed-element verification need F04/F06/F10. ADCIRC-only benchmark scope does not itself remove library format support (S:13,31; C:64-69).

## Three research spot-checks

| Claim | Repository evidence | Result |
|---|---|---|
| R4 prune scripts (P:15) | MANIFEST.in:28 | [verified] Matches exactly |
| R5 readlines (P:16) | src/chilmesh/fort14_io.py:296-298 | [verified] Matches exactly; memory estimate not validated |
| R6 index/query behavior (P:17) | src/chilmesh/CHILmesh.py:420-422,971-974,1033-1043 | [verified] Centroid cKDTree, capped 20 candidates and loop match; lazy-build claim fails (F09) |

[inferred] Counts: 11 findings = 0 CRITICAL, 2 HIGH, 7 MEDIUM, 2 LOW (F01-F11). 20 requirements = 14 FR + 6 SC (S:82-95,108-113); 19/20 have task associations, FR-014 boundary-only. 45 tasks (T:11-77), 44 trace to spec/plan, T045 dispatch-only. Decisions 3/3 mapped; source spot-checks 2 full matches, 1 partial (tables above).

## OUT OF SCOPE / INCOMPLETE:

[verified] Starting assumptions: use dispatch's five-step analyze procedure; inspect local artifacts only; no network or benchmark execution. Unchanged. No separate installed analyze procedure claimed (dispatch TASK/CONSTRAINTS).
[verified] Only file created: specs/006-competitor-benchmarks/analysis.md (report-write command). No spec/plan/task/source edits intended; scope gate below verifies changed-file allowance.
[inferred] Artifact review complete; implementation readiness needs F03/F04 resolved and remaining findings triaged. Rights, competitor licensing/runtime availability, governance classification and measured performance remain unverified (S:119-122; P:80-87); no external research authorized. LESSON-CANDIDATE F04/F10 relayed in final response; no canon edits.
Gate evidence pending execution; final response records actual commands, output and exit codes.
WORKERS SPAWNED: 0
