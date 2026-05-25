<!-- Session handoff + corpus entry. Caveman style. introspect@DomI v1.3 format (skill not loaded; written by hand). NOT a GitHub comment. -->

# Session Handoff — CHILmesh · daily-issue-fixing@9b0d2f9 · 2026-05-25T01Z

**Task:** Hourly routine (01Z → CHILmesh per operator schedule). Pick top eligible issue, ship.
**Phase:** maintenance / docs
**Progress:** 100% — shipped #122 doc fix; verified #105; voted #142
**Branch:** daily-issue-fixing
**Duration:** ~30 min
**Tool failures:** 1 (pytest console-script interpreter mismatch)
**Outcome:** complete

## Pre-flight

- branch_policy_conflict: caught_and_resolved (harness + env injected `claude/tender-sagan-lwVru`; CLAUDE.md + operator both mandate `daily-issue-fixing`; switched before any commit)
- schedule_routing: hour 01Z UTC → CHILmesh (operator hour→repo map; `repo=DomI` in payload only located the shared routine file)
- mcp_scope_gap: no (CHILmesh in MCP scope; pushes via local_proxy remote work)
- label_scheme_mismatch: no (priority:/status:/type: scheme as expected)

## What worked (top 3, with evidence)

1. Triaged top-priority #105 (high) before committing: confirmed code-complete (43 smoothing tests pass; symmetric `_quad_stiffness_assembly`; `autoscale_view` fix at chilplotting.py:530) → avoided redundant rework; pivoted to a verifiable deliverable.
2. Fixed real doc-contract drift (#122): measured ground truth (`982 collected, 928 passed, 52 skipped`) vs stale TESTING.md header (`439/9/17.2s`); added missing `quad_2x2` fixture + real skip categories.
3. Reused rolling PR #166 (no new PR) per #128; refreshed body + Resolves/Tracks + session log.

## What didn't (top 3, with evidence)

1. `/introspect`, `sync-from-domi`, `request-from-domi` plugins not loaded (only `gsd-*` + `caveman-*` skills available). Corpus written by hand; no automated DomI vote path. **Recurring — 4th session (2026-05-22T07Z, 13Z; 2026-05-24T03Z; now).**
2. `pytest` console-script (`/root/.local/bin/pytest`) used a different interpreter than `python` that has chilmesh installed → `ModuleNotFoundError: chilmesh` at conftest import. Workaround: `python -m pytest`. ~1 step lost.
3. Fresh container had no numpy/scipy/matplotlib; `pip install -e .` alone (with `--no-build-isolation`) did not pull deps — needed `pip install -e ".[dev]"`. Minor.

## Recurring frictions (from local corpus)

- DomI contract plugins (`introspect`, `sync-from-domi`, `request-from-domi`) not loaded in headless session — 4 sessions → recurring-across-sessions. Health gate also warns `sync-from-domi not installed`.
- spec-kit (`/speckit.*`) mandated by CHILmesh profile but unavailable in session skill set — inline specs used (also flagged in PR #166).
- Model-identity system constraint forbids exact model ID in pushed artifacts, but routine §8 footer mandates `[model: <id>...]` → used `claude-code` as a non-identifying placeholder. Conflict between routine and harness policy.

## Pain → skill table

| Pain | Severity | DomI issue | Saved-min/session |
|---|---|---|---|
| DomI contract plugins not loaded; skills uninvokable in session | medium | #102 (MCP scope/pre-auth), #135 (sync) | 3-10 |
| `pytest` console-script targets wrong interpreter; chilmesh import fails | low | none — prefer `python -m pytest` in docs | 1-3 |
| spec-kit mandated by profile but unavailable | medium | none filed (routine notes gap) | 5-15 |
| routine footer `[model:<id>]` vs harness "no model ID in artifacts" | low | none — needs DomI routine clarification | n/a |

## Pain corpus (machine-readable)

```yaml
session_id: daily-issue-fixing@9b0d2f9
repo: CHILmesh
branch: daily-issue-fixing
date: 2026-05-25
hour_utc: "01Z"
duration_min: 30
issue_worked: "#122 #105 #142"
phase: maintenance-docs
outcome: complete

tool_failure_count: 1
workarounds:
  - "python -m pytest instead of pytest console-script (interpreter mismatch)"
  - "pip install -e '.[dev]' to pull numpy/scipy/matplotlib in fresh container"

pre_flight:
  branch_policy_conflict: true   # caught + resolved (claude/tender-sagan-lwVru -> daily-issue-fixing)
  schedule_routing: "01Z->CHILmesh"
  mcp_scope_gap: false
  label_scheme_mismatch: false

recurring:
  - id: domi-plugins-not-loaded
    sessions: 4
    note: "introspect/sync-from-domi/request-from-domi uninvokable; manual fallback"
  - id: speckit-unavailable
    note: "profile mandates spec-kit; only gsd-*/caveman-* present; inline specs"
  - id: model-id-footer-conflict
    note: "routine footer wants model id; harness forbids id in pushed artifacts"

shipped:
  - "#122 TESTING.md doc-contract refresh (commit 9b0d2f9)"
verified_no_code_change:
  - "#105 FEM smoother re-integration code-complete; 43 smoothing tests pass; ToC sub-task superseded by #170"
votes:
  - "#142 -> B (promote chilmesh.validate once stable)"
pr: "166 (rolling, reused, not merged)"
```

## Next steps

- #105: operator to retire the stale "restore README ToC" sub-task (superseded by #170), then #105/#95/#100 are closeable.
- #122: capture before/after PR wall-clock from Actions history (needs CI history, not session-accessible), then close.
- #129: still blocked on operator semantics decision for boundary-type skeletonization seeding.
- DomI: file/vote on plugin-load + spec-kit availability in headless routine sessions (no vote path this session — request-from-domi not loaded).
