<!-- Session handoff + corpus entry. Caveman style. -->

# Session Handoff — CHILmesh · development_b764338 · 2026-06-09

**Task:** #199 admesh-gmsh-io draft rehome/cleanup + #200 DomI sync
**Phase:** implementation
**Progress:** complete — both issues closed; major pre-existing test-redness surfaced (operator-only)
**Branch:** development
**Duration:** ~25 min
**Tool failures:** 0
**Outcome:** complete (assigned work) + research (test-suite finding)

## Pre-flight

- branch_policy_conflict: caught_and_resolved   <!-- harness put session on claude/cool-ptolemy-30lytp; switched to development -->
- mcp_scope_gap: no
- label_scheme_mismatch: no

## What worked (top 3, with evidence)

1. #199 disambiguated fast by checking BOTH siblings: ADMESH already has native `src/admesh/gmsh.py`; CHILmesh already has native `src/chilmesh/gmsh_io.py` (23/23 tests green) → draft provably superseded → safe delete, no architecture decision needed.
2. Drift resolved cheaply: CHILmesh consumes DomI via plugins (no vendored skills) → sync = pin refresh only (`074e4a0 → 69fdeb7`), not a full artifact pull.
3. Test-failure triage via signature-dedup (`sed` normalize + `uniq -c`) collapsed 237 failures into ~8 root causes in one run → identified the daily-maintenance divergence pattern, avoided rabbit-hole.

## What didn't (top 3, with evidence)

1. `development` test suite is badly red (~237 failed / 596 passed / 87 errors, block_o excluded) — NOT this session's doing (commits touch only `.domi-pin` + `.planning/`). Burned ~6 tool calls proving it pre-existing before I could trust the gmsh delete.
2. Fresh container shipped no chilmesh/pytest — venv build tax paid again (#148 pattern).
3. Ready-queue thin again: only #199 cleanly executable; #155 env-blocked, #167 GPU-blocked, #187 awaits operator, brainstorming issues need decisions. Same shape as prior sessions.

## Recurring frictions (from local corpus)

- Harness injects `claude/*` session branch every run — observed in many prior sessions (CLAUDE.md branch-sprawl log).
- `daily-maintenance ↔ development` divergence — observed in 2026-06-08 corpus + CLAUDE.md note; now manifests as a red test suite on development.

## Pain → skill table

| Pain | Severity | DomI issue | Saved-min/session |
|---|---|---|---|
| venv rebuild tax on fresh container | low | ensure-test-venv (exists) | 2 |
| sync-from-domi plugin not loaded at container start → manual pin refresh | low | #114 (known) | 3 |

## Pain corpus (machine-readable)

```yaml
session_id: development_b764338
repo: CHILmesh
branch: development
date: 2026-06-09
duration_min: 25
issue_worked: "#199, #200"
phase: implementation
outcome: complete

tool_failure_count: 0
workarounds:
  - "sync-from-domi plugin absent → refreshed .domi-pin by hand from local DomI origin/main"
  - "fresh container had no chilmesh/pytest → built .venv, pip install -e . pytest"

pre_flight:
  branch_policy_conflict: true
  mcp_scope_gap: false
  label_scheme_mismatch: false
  notes: "harness branch claude/cool-ptolemy-30lytp → switched to development per CLAUDE.md precedence"

pain_points:
  - pain: "development test suite ~237 failed (API drift: elem_quality tuple API, copy, paths_on_outer_vertices, advancing_front_*, write_to_fort14, get_layer, bbox vs bounding_box, fort14 int('1.0') parse)"
    frequency: recurring-across-sessions
    severity: high
    evidence: "pytest tests/ -k 'not block_o' → 237 failed/596 passed/87 errors; signatures dedup'd to ~8 root causes; matches daily-maintenance 45-behind note in CLAUDE.md"
    existing_skill_should_have_caught_it: "no"
    missing_skill_would_have_prevented_it: "no — operator-only branch reconciliation, explicitly out-of-scope for autonomous sessions"
    domi_issue: null
    saved_time_estimate_min: 0
    tokens_wasted: 8000
  - pain: "venv rebuild tax on fresh container"
    frequency: recurring-across-sessions
    severity: low
    evidence: "no chilmesh/pytest on import; pip install -e . needed before gate"
    existing_skill_should_have_caught_it: "ensure-test-venv"
    missing_skill_would_have_prevented_it: "no"
    domi_issue: "#148"
    saved_time_estimate_min: 2
    tokens_wasted: 1000

actions_taken:
  votes_cast: []
  new_requests_filed: []
  closed_issues_flagged_for_reopen: []
  introspect_design_proposal_on_9: false

introspection_meta:
  what_worked: "cross-sibling check to disambiguate a rehome chore into a safe delete"
  what_was_hard: "distinguishing pre-existing repo-wide redness from my own change before trusting a delete"
```

## Next session — pick up here

1. [ ] OPERATOR: reconcile `daily-maintenance → development` — development is missing the API the tests/chilplotting/bridge depend on (`elem_quality` tuple form, `copy`, `paths_on_outer_vertices`, `advancing_front_*`, `write_to_fort14`, `get_layer`, `bounding_box` key). High regression risk — operator-only per CLAUDE.md.
2. [ ] Once reconciled, re-run full `pytest tests/` (incl. block_o) to confirm green before any `development → main` merge.
3. [ ] #155 / #167 / #168 stay env-blocked (STOFS file needs hosted runner; GPU device absent).

**Files to read first:**
- `.claude/CLAUDE.md` — branch policy + daily-maintenance divergence note (2026-06-08).
- `src/chilmesh/CHILmesh.py` — `element_quality` (1174) vs callers' `elem_quality`; the API-drift epicenter.
- PR #195 body — ⚠ operator-finding block enumerates the 8 root causes.

**Context to remember:**
- This session touched ONLY `.domi-pin` + `.planning/` — it did not cause the test redness.
- gmsh I/O fully owned by CHILmesh via `src/chilmesh/gmsh_io.py` (23/23); the admesh draft is gone.

## Routing decisions taken this session

- Votes on existing skill-proposal issues: 0 (probation active, #203)
- New requests filed: 0
- Closed issues flagged for reopen: 0
- Comments on DomI #9: 0
- PR description updated: yes (#195)

---
_Written via `introspect@DomI` v1.3 from CHILmesh (inline fallback — plugin not loaded at container start, #114). Caveman style. Pairs with `handoff@DomI` v1.0._
