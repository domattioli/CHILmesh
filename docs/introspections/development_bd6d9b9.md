<!-- Session handoff + corpus entry. Caveman style. -->

# Session Handoff — CHILmesh · development_bd6d9b9 · 2026-06-10

**Task:** #204 — ctor opt-out flags (build_spatial_indices, validate)
**Phase:** implementation
**Progress:** 100% — primary ask shipped + closed
**Branch:** development
**Duration:** ~25 min
**Tool failures:** 1 (raw routine-file fetch 404 → fell back to local DomI checkout)
**Outcome:** complete

## Pre-flight

- branch_policy_conflict: caught_and_resolved (harness injected `claude/amazing-cray-33dpmk`; CLAUDE.md precedence → checked out `development`)
- mcp_scope_gap: no
- label_scheme_mismatch: no

## What worked (top 3, with evidence)

1. Hour→repo routing + profile lookup clean (hour 09 → CHILmesh; profile knob row + §6 confirmed branch=development, gate=`pytest tests/`).
2. Issue triage avoided dead ends — #155 (priority:now) correctly identified as env-blocked (global 1.73GB STOFS mesh + GPU solver), skipped; picked #204 as high-confidence additive win.
3. Haiku subagent one-shot implemented the edit faithfully (20 ins/7 del, additive, defaults preserved); gate green 931 passed.

## What didn't (top 3, with evidence)

1. Routine raw fetch `raw.githubusercontent.com/.../22cc443/.../claude_routine_instructions.md` → 404 (commit-pinned path); recovered via local `/home/user/DomI/.claude/` copy (NOT a stop condition per routine).
2. caveman plugin not loaded at container start (DomI #114) → `Unknown skill` on `/caveman:caveman`; emulated ultra from SKILL.md inline.
3. Fresh container shipped no test venv (numpy/scipy/pytest/editable chilmesh absent) — had to build `.venv` before the gate could run (#148 tax).

## Recurring frictions (from local corpus)

- caveman plugin absent at container start — observed in multiple prior sessions (#114).
- test-venv build tax on pytest-gated repos — observed recurrently (#148).

## Pain corpus (machine-readable)

```yaml
session_id: development_bd6d9b9
repo: CHILmesh
branch: development
date: 2026-06-10
duration_min: 25
issue_worked: "#204"
phase: implementation
outcome: complete
tool_failure_count: 1
workarounds:
  - "routine raw fetch 404 at pinned commit -> read local DomI/.claude copy"
  - "caveman /skill unknown -> emulate ultra from SKILL.md"
  - "no test venv -> python -m venv .venv && pip install -e .[dev]"
pains:
  - desc: "caveman + DomI contract plugins not loaded at container start"
    severity: low
    domi_issue: "#114"
    tokens_wasted: "~1 extra skill call + 1 line"
  - desc: "pytest-gated repo ships no venv; per-session build tax"
    severity: med
    domi_issue: "#148"
    tokens_wasted: "~1 build step (~30s wall)"
notes: "probation active (#191/#203) -> no skill-request votes/filing; pains routed to corpus only."
next_steps:
  - "#204 second hotspot: cache/decrement degree array in layer_paths.py:173 (degree_remaining 2.85M calls) — separate layer_paths change."
  - "QuADMesh side: wire identify_edges_in_layer to pass build_spatial_indices=False, validate=False (downstream)."
```
