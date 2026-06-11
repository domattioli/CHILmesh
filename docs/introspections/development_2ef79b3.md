<!-- Session handoff + corpus entry. Caveman style. -->

# Session Handoff — CHILmesh · development_2ef79b3 · 2026-06-08

**Task:** #187 layer-lexicon scope creep
**Phase:** planning
**Progress:** decision artifact shipped; rename deferred to operator naming pick
**Branch:** development
**Duration:** ~20 min
**Tool failures:** 0
**Outcome:** partial (unblocked — awaiting operator decision)

## Pre-flight

- branch_policy_conflict: caught_and_resolved   <!-- harness put session on claude/modest-archimedes-3ypj3g; switched to development -->
- mcp_scope_gap: no
- label_scheme_mismatch: no

## What worked (top 3, with evidence)

1. Hour→repo routing + issue-queue triage cut quickly to executable work (22Z → CHILmesh; queue scan filtered #155 recently-worked, #163 blocked, #167 GPU-env-blocked, #168 rec done).
2. Lexicon inventory via grep counts surfaced root cause fast (`_skeletonize` = MATLAB `meshLayers` port, mislabeled as medial-axis) — counts: layer 290 / skeletoniz 48 / peel 7 / front 1.
3. Recognized rename = public-API break → shipped decision doc instead of unilateral rename (Principle VII compliance).

## What didn't (top 3, with evidence)

1. Ready-queue thin: of 7 open issues, most genuinely-actionable ones env-blocked (GPU #167, STOFS file #155) or awaiting operator (#142/#187/#153 brainstorming). Routine had to settle for a planning deliverable.
2. Fresh container shipped no test deps — venv build tax paid again (#148 pattern); `pip install -e .` ~OK but cost a step.
3. `#168` rec#1 already done in prior session — wasted one verification grep to confirm before skipping.

## Recurring frictions (from local corpus)

- Harness injects `claude/*` session branch every run — observed in many prior sessions (CLAUDE.md branch-sprawl log).
- Big-mesh (WNAT/STOFS) benchmark perpetually env-blocked — observed across #155 thread (≥5 sessions).

## Pain → skill table

| Pain | Severity | DomI issue | Saved-min/session |
|---|---|---|---|
| venv rebuild tax on fresh container | low | ensure-test-venv (exists) | 2 |
| thin ready-queue → planning-only sessions | low | n/a (operator triage) | — |

## Pain corpus (machine-readable)

```yaml
session_id: development_2ef79b3
repo: CHILmesh
branch: development
date: 2026-06-08
duration_min: 20
issue_worked: "#187"
phase: planning
outcome: partial

tool_failure_count: 0
workarounds:
  - "switched off harness claude/* branch to development at bootstrap"

pre_flight:
  branch_policy_conflict: true
  mcp_scope_gap: false
  label_scheme_mismatch: false
  notes: "harness branch claude/modest-archimedes-3ypj3g overridden to development per CLAUDE.md"

pain_points:
  - pain: "ready-queue thin; most actionable issues env-blocked or operator-gated"
    frequency: recurring-across-sessions
    severity: low
    evidence: "7 open issues; #155 recently-worked, #163 blocked, #167 GPU-no-device, #168 done, #142/#153/#187 brainstorming"
    existing_skill_should_have_caught_it: no
    missing_skill_would_have_prevented_it: no
    domi_issue: null
    saved_time_estimate_min: 0
    tokens_wasted: unknown

actions_taken:
  votes_cast: []                # probation #203 — voting suspended
  new_requests_filed: []        # probation #203
  closed_issues_flagged_for_reopen: []
  introspect_design_proposal_on_9: false

introspection_meta:
  what_worked: "inventory-by-grep → root-cause naming conflation in 2 calls"
  what_was_hard: "deciding to NOT rename (API stability) vs operator's 'revise top to bottom' ask"
```

## Next session — pick up here

- If operator checked the 3 boxes in `docs/LEXICON_PROPOSAL.md` → execute the mechanical rename sweep (deprecation-shim public symbols, rename privates) as a single `type: refactor` PR.
- Else next CHILmesh slot: #142 validator-promotion vote, or #155 STOFS only if a hosted runner / iterative solver lands.
