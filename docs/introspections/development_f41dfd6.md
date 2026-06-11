<!-- Session handoff + corpus entry. Caveman style. -->

# Session Handoff — CHILmesh · development_f41dfd6 · 2026-06-09

**Task:** routine #199/#200 → operator-directed CI-badge greening
**Phase:** implementation
**Progress:** complete — CI PR lane 992 passed / 0 failed (deterministic ×4), coverage 82%
**Branch:** development
**Duration:** ~90 min (extended, two passes)
**Tool failures:** 0
**Outcome:** complete

## Pre-flight

- branch_policy_conflict: caught_and_resolved   <!-- harness claude/cool-ptolemy-30lytp → development -->
- mcp_scope_gap: no
- label_scheme_mismatch: no

## What worked (top 3, with evidence)

1. Signature-dedup of failures (`sed` normalize + `uniq -c`) collapsed 237 → ~8 root causes → found they were a dropped-API regression, not 237 bugs.
2. `git log --all -S "def <method>"` located commit `fa50a89` holding the full API → faithful restore, zero reinvention (8 methods).
3. Root-caused the xdist flake to `write_to_fort14` mutating `self.grid_name` (poisoned shared cache via direct `examples.annulus()` callers) → fixed the writer, not the symptom; reverted a wrong conftest-copy hypothesis.

## What didn't (top 3, with evidence)

1. First flake hypothesis wrong: added conftest copy-on-handout → still flaked (direct callers bypass fixtures). Cost ~2 validation rounds before finding the real `write_to_fort14` mutation.
2. caveman plugin un-loadable mid-session (DomI #114) — operator flagged twice; only fixable thing was marketplace owner-case for next container.
3. fresh container venv tax again (#148): no chilmesh/pytest/xdist/cov until installed.

## Recurring frictions (from local corpus)

- DomI contract plugins not loaded at container start — DomI #114, every session.
- Shared mutable fixture cache + in-place mutators → nondeterministic xdist failures (latent; `fresh_mesh` exists but direct `examples.X()` callers bypass all guards).

## Pain → skill table

| Pain | Severity | DomI issue | Saved-min/session |
|---|---|---|---|
| plugins not loaded at container start | med | #114 | 5 |
| venv rebuild tax | low | #148 | 2 |

## Pain corpus (machine-readable)

```yaml
session_id: development_f41dfd6
repo: CHILmesh
branch: development
date: 2026-06-09
duration_min: 90
issue_worked: "#199, #200, CI-badge (no issue)"
phase: implementation
outcome: complete

tool_failure_count: 0
workarounds:
  - "sync-from-domi plugin absent → refreshed .domi-pin by hand"
  - "caveman plugin absent → emulated ultra; fixed marketplace owner-case for next container"
  - "fresh container → built .venv, pip install -e . pytest pytest-xdist pytest-cov"

pre_flight:
  branch_policy_conflict: true
  mcp_scope_gap: false
  label_scheme_mismatch: false
  notes: "harness claude/cool-ptolemy-30lytp → development per CLAUDE.md precedence"

pain_points:
  - pain: "CI badge red ~237 failures = dropped CHILmesh API (impl regressed behind tests + chilplotting/bridge/cli consumers)"
    frequency: once
    severity: high
    evidence: "pytest -n auto -m 'not slow' 237 failed → restored 8 methods from fa50a89 → 0 failed; coverage 82%"
    existing_skill_should_have_caught_it: "no"
    missing_skill_would_have_prevented_it: "a pre-merge CI-status gate on branch consolidation merges would have caught the dropped API at #195 merge time"
    domi_issue: null
    saved_time_estimate_min: 0
    tokens_wasted: 0
  - pain: "xdist-ordering flake: write_to_fort14 mutated self.grid_name, poisoning shared fixture cache"
    frequency: recurring-this-session
    severity: medium
    evidence: "test_find_element_center_of_mesh[annulus] assert -1>=0; nondeterministic; fixed by save/restore in writer"
    existing_skill_should_have_caught_it: "no"
    missing_skill_would_have_prevented_it: "no"
    domi_issue: null
    saved_time_estimate_min: 0
    tokens_wasted: 4000
  - pain: "caveman plugin not loadable mid-session"
    frequency: recurring-across-sessions
    severity: medium
    evidence: "Skill /caveman:caveman → Unknown skill x2; operator annoyed; #114"
    existing_skill_should_have_caught_it: "no"
    missing_skill_would_have_prevented_it: "no — container-start declarative load only"
    domi_issue: "#114"
    saved_time_estimate_min: 5
    tokens_wasted: 1500

actions_taken:
  votes_cast: []
  new_requests_filed: []
  closed_issues_flagged_for_reopen: []
  introspect_design_proposal_on_9: false

introspection_meta:
  what_worked: "git -S history archaeology to restore a dropped API faithfully instead of reinventing"
  what_was_hard: "a nondeterministic xdist flake whose poisoner was a direct cached-mesh caller bypassing fixture guards"
```

## Next session — pick up here

1. [ ] Confirm `python-package.yml` push-lane run on `f41dfd6` green on GitHub (local proxy only locally; full suite incl block_o/slow runs in CI).
2. [ ] Operator: rolling PR #195 `development → main` merge (operator-only) + the separate `daily-maintenance` (#182, 45-commit) reconciliation still pending.
3. [ ] Low-pri: #187 lexicon naming (awaits operator), #155 STOFS/GPU env-blocked.

**Files to read first:**
- `src/chilmesh/CHILmesh.py` — restored API (`elem_quality`, `write_to_fort14`, advancing-front, `copy`, `get_layer`, `paths_on_outer_vertices`).
- `tests/test_fort14_boundary_types.py` — boundary round-trip contract.
- PR #195 body §6 — green-state summary.

**Context to remember:**
- Implementation dispatched to Haiku subagents (cavecrew-builder, caveman-prefixed) per CLAUDE.md; two ≤6-line test-fix edits (writer no-mutate, conftest revert) done inline — deviation noted.
- Restore source-of-truth = commit `fa50a89`.

## Routing decisions taken this session

- Votes on existing skill-proposal issues: 0 (probation #203)
- New requests filed: 0
- Closed issues flagged for reopen: 0
- Comments on DomI #9: 0
- PR description updated: yes (#195)

---
_Written via `introspect@DomI` v1.3 inline (plugin not loaded — #114). Pairs with `handoff@DomI` doc `.claude/handoffs/2026-06-09-chilmesh-ci-green.md` (operator explicitly requested both this session). Caveman style._
