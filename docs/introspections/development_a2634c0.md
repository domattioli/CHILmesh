<!-- Session handoff + corpus entry. Caveman style. -->

# Session Handoff — CHILmesh · development_a2634c0 · 2026-06-09

**Task:** #203 (skeletonize non-termination) + #202 verification
**Phase:** bugfix / hardening
**Progress:** #203 closed (fixed+guarded+tested); #202 commented (problem 1 fixed, problem 2 no-repro)
**Branch:** development
**Duration:** ~30 min
**Tool failures:** 0
**Outcome:** complete

## Pre-flight

- branch_policy_conflict: caught_and_resolved   <!-- harness put session on claude/amazing-cray-ii32ht; switched to development per branching.md -->
- mcp_scope_gap: no
- label_scheme_mismatch: no

## What worked (top 3, with evidence)

1. Verify-before-fix paid off → both #203 and #202 turned out already-fixed on `development` HEAD (#155 `3492f03`, #202 `e56ac9e`). Reproducing live (Block_O default read 0.11s, n_layers=9; 10-elem mixed 0.0013s; 400-mesh fuzz 0 hangs) prevented re-fixing solved bugs.
2. git-log triage of recent commits surfaced the already-landed fixes fast (grep `skelet|backend|#202`), reframing the task from "fix hang" → "harden + close stale reports".
3. Root-caused the discrepancy: reporter envs (QuADMesh/MADMESHing) ran STALE sibling `../CHILmesh` checkouts predating the fixes → false bug reports. Real deliverable = defensive guard + mixed-element regression test (donut test only covered pure-tri).

## What didn't (top 3, with evidence)

1. Could not reproduce the exact #203 hang — `Mixed_Test.14` is a QuADMesh fixture, absent from CHILmesh repo. Spent effort fuzzing + analyzing the loop before concluding it already terminates.
2. Fresh container shipped no numpy/scipy — venv build tax paid again (recurring; #148 / ensure-test-venv pattern).
3. Two issues (#202, #203) were near-duplicates of work already on `development` — symptom of cross-repo stale-checkout reporting.

## Recurring frictions (from local corpus)

- **Cross-repo stale sibling checkout → false bug reports.** QuADMesh/MADMESHing install chilmesh from `../CHILmesh`; if that checkout lags `development`, downstream sessions file bugs against fixes that already shipped. Observed twice this session (#202, #203).
- Harness injects `claude/*` session branch every run (CLAUDE.md branch-sprawl log).
- venv rebuild tax on fresh container (≥5 sessions).

## Pain → skill table

| Pain | Severity | DomI issue | Saved-min/session |
|---|---|---|---|
| stale sibling checkout → duplicate/false downstream bug reports | med | (candidate: sibling-pin-check) | 15 |
| venv rebuild tax on fresh container | low | ensure-test-venv (exists) | 2 |

## Next steps

- #202 problem 2 (optional `pip install -e` cpp build-hook) left for operator — slow path no longer reproduces, so low urgency.
- #203 direction #2 (dispatch read path to `chilmesh_cpp.full_init` when `CPP_AVAILABLE`) → perf enhancement, lives with C++ availability work in #163.
- Consider a downstream guard: sibling consumers (QuADMesh #76, MADMESHing) should pin/verify chilmesh checkout == `development` HEAD before filing perf/hang bugs.

## Open questions

- Should pure-Python remain the source-install default now that it's fast (Block_O 0.11s), with cpp as opt-in build? (#202 problem 2)
