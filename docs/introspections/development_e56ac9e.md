<!-- Session handoff + corpus entry. Caveman style. v1.3 — corpus canonical, not a GitHub comment. -->

# Session Handoff — CHILmesh · development_e56ac9e · 2026-06-09

**Task:** #202 — backend_info() falsely reports cpp/rust available when no compiled extension built
**Phase:** implementation
**Progress:** 100% of problem 1 (misreporting); problem 2 (no-cpp-build perf cliff) left open
**Branch:** development
**Duration:** ~35 min
**Tool failures:** 0
**Outcome:** complete (bounded scope)

## Pre-flight

- branch_policy_conflict: caught_and_resolved   <!-- harness injected `claude/amazing-cray-pnjv8z`; switched to `development` per CLAUDE.md + branching.md -->
- mcp_scope_gap: no
- label_scheme_mismatch: no

## What worked (top 3, with evidence)

1. Reuse-over-reinvent — fix gated backend_info() on existing `CPP_AVAILABLE`/`RUST_AVAILABLE` flags already correct in backends/*.py (per #163). No new detection logic. (`e56ac9e` diff: -10/+7 lines)
2. Re-verify env claim (#223) — reproduced bug in THIS container (`{'selected':'cpp'}` with no `.so`) before fixing; did not trust issue body alone. Confirmed flags `False` + `chilmesh_cpp.__file__ is None` (namespace stub).
3. Fix surfaced latent false-pass — `test_backend_equivalence.py` had `RUST_AVAILABLE = ... or True`; the rust test passed only because backend_info() was lying. Caught when fix turned it red; corrected the skipif source. Net: a real regression test now exists.

## What didn't (top 3, with evidence)

1. chilmesh not pip-installed in fresh container → ~90s venv build before any validation could run (`pip install -e ".[dev]"`). Per-session tax (#148).
2. caveman plugin not loaded at container start (DomI #114) → emulated ultra from SKILL.md. Recurring across every cloud routine.
3. Docstring example for backend_info() still shows `>>> info['selected']` → `'cpp'` — now misleading in a no-ext env. Left unchanged (not doctested; out of bounded scope). Minor debt.

## Recurring frictions (from local corpus)

- caveman plugin absent at container start — observed in many prior sessions (DomI #114)
- fresh-container venv tax before pytest gate — observed across CHILmesh/QuADMESH routines (#148)

## Pain → skill table

| Pain | Severity | DomI issue | Saved-min/session |
|---|---|---|---|
| Capability flag lies on namespace-stub import | low | n/a (repo bug, fixed) | — |
| venv build tax before gate | low | ensure-test-venv exists | ~2 |

## Pain corpus (machine-readable)

```yaml
session_id: development_e56ac9e
repo: CHILmesh
branch: development
date: 2026-06-09
duration_min: 35
issue_worked: "#202"
phase: implementation
outcome: complete

tool_failure_count: 0
workarounds:
  - "manual venv build (.venv) + editable install before pytest gate could run"

pre_flight:
  branch_policy_conflict: true
  mcp_scope_gap: false
  label_scheme_mismatch: false
  notes: "harness branch claude/amazing-cray-pnjv8z -> switched to development per CLAUDE.md/branching.md"

pain_points:
  - pain: "capability-detection helper trusted bare `import name` success, which passes for an importable-but-empty namespace-package stub"
    frequency: once
    severity: low
    evidence: "backend_info() reported selected:cpp with no .so; backends/*.py already had hasattr-gated flags"
    existing_skill_should_have_caught_it: false
    missing_skill_would_have_prevented_it: false
    domi_issue: null
    saved_time_estimate_min: 0
    tokens_wasted: 0
  - pain: "fresh container has no chilmesh/numpy/pytest; ~90s venv build before validation"
    frequency: recurring-across-sessions
    severity: low
    evidence: "pip install -e .[dev] required before any pytest run; #148"
    existing_skill_should_have_caught_it: true
    missing_skill_would_have_prevented_it: false
    domi_issue: null
    saved_time_estimate_min: 2
    tokens_wasted: 0

actions_taken:
  votes_cast: []          # probation active (#191/#203) — no voting
  new_requests_filed: []  # probation active — pain routed to corpus tokens_wasted only
  closed_issues_flagged_for_reopen: []
  introspect_design_proposal_on_9: false

introspection_meta:
  what_worked: "reproduce-then-fix; reuse existing flags; bounded scope kept diff to 3 files"
  what_was_hard: "fix correctly turned a previously-green-but-wrong test red — had to trace it to an `or True` skipif bug, not a regression"
```

## Next session — pick up here

- #202 problem 2 still OPEN: source/editable install ships no compiled cpp ext → silent slow pure-Python (Block_O init >200s). Options: (a) PEP 517 build hook builds ext on `pip install -e`, (b) document explicit build step, (c) de-pathologize pure-Python skeletonization on Block_O. Needs operator direction on which.
- #155 (priority: now) recently worked (<24h) — skipped this session per recency filter; ready when window clears.
- Optional tidy: update backend_info() docstring example (`>>> info['selected']` → no longer always 'cpp').
