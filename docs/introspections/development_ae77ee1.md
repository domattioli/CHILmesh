<!-- Session handoff + corpus entry. Caveman style. R5 frontmatter per introspect@DomI tmpl. -->
---
date: 2026-06-13
session: 2026-06-13T08Z-rotation
repo: domattioli/CHILmesh
severity: med
freq: recurring
issues: [211, 286]
wasted_min: 6
wasted_tok: 4000
missing_skill: null
---

# Session Handoff — CHILmesh · development_ae77ee1 · 2026-06-13 (rotation hour-08)

**Task:** overhaul rotation, CHILmesh slot. C spec-048 slice shipped 06-11/02Z → maintenance track (issue-queue top) + hub loop.
**Phase:** maintenance
**Progress:** complete — #211 fixed + regression-locked, pushed to rolling PR #210
**Branch:** development (rolling PR #210)
**Duration:** ~30 min
**Tool failures:** 0
**Outcome:** complete

## Pre-flight

- branch_policy_conflict: caught_and_resolved   <!-- harness claude/determined-gauss-sohefo → development per CLAUDE.md precedence; harness branch == origin/main (0/0) -->
- domi_pin_drift: none on development (3e46639 current; pin on main/claude-branch is stale 39fd74a but that's pre-02Z-promotion, not real drift)
- caveman_plugin: NOT loaded → /caveman:caveman ultra returned `Unknown skill` → emulated from SKILL.md (honest fallback per #168). SessionStart resume hook NOT active this container.
- health_check: exit 0 but DEAD GATE — `sync-from-domi not installed` → warn+continue (DomI#286 class). Pin happened current so no exposure.

## What shipped (evidence)

1. `ae77ee1` fix #211 — padded-triangle sentinel. `split_triangle` (mutations.py:80) + `split_triangles` (:488) tested `elem[2] != elem[3]`; padding convention is `[v0,v1,v2,v0]` (every other check `row[3]==row[0]`), so post-`merge_elements` 4-column meshes wrongly raised `ValueError` on a remaining triangle. Latent `_point_in_element` (CHILmesh.py:973) `elem[3]==elem[2]`→`==elem[0]`. 3 one-token source fixes.
2. 3 regression tests in `test_mutations.py` (donut param runs real; annulus skips no-adjacent-pair). Gate: full fast suite **992 passed / 47 skipped, 0 regressions** (63s, `not block_o`); was 989/44 pre-session.
3. Hub: commented DomI#286 (CHILmesh = 3rd affected consumer of dead drift-gate); MADMESHing#48 checklist.

## Key decisions

1. Picked #211 (queue top, fresh `type: bug`, exact location+fix given, bounded/testable) over the brainstorm/research backlog (#201/#155/#167). Serves #48: mixed-element correctness is core to CHILmesh's role as MADMESHing's hard dep.
2. **Did NOT port the offline drift-gate fallback** to CHILmesh despite the dead gate. DomI#286 explicitly says the per-repo copies (Valence#147, ADMESH#150) ARE the redundancy #48 targets and asks to canonicalize ONE helper. A 4th hand-rolled copy is the anti-pattern. Pin current on dev → no active exposure → held for DomI hub (hour-12) + noted CHILmesh affected on #286. Will wire canonical helper via `/sync from DomI` once shipped.
3. Kept EDIT 3 (latent _point_in_element) despite it being currently benign — matches canonical convention, prevents a future silent break if `_point_in_quad` stops degenerating correctly. Made its test actually prove the branch routing (monkeypatch), not just the black-box answer.

## What worked (top 3)

1. **Revert experiment as the verification gate.** Did not trust the subagent's "95 passed". Reverted each fix, confirmed the new tests FAIL pre-fix. Caught that test (c) was a false-pass before commit.
2. Empirical bug-site confirmation (grep convention across the file: 10 sites use `[3]==[0]`, 2 use the wrong sentinel) before dispatching — gave Haiku exact line+token edits, one clean round + one test-rework round.
3. Haiku builder / Opus review split held; reviewer added the real signal the builder missed.

## What didn't (pains → routing)

1. **Subagent false-pass regression test (#168 class, recurring, severity med).** Haiku reported "95 passed, 3 skipped" — green — but (a) the new tests SKIPPED on annulus and only ran on donut, and (b) test (c) PASSED even with the EDIT-3 bug reverted → ZERO signal on the latent fix. A regression test that passes pre-fix proves nothing. Orchestrator caught it only by running the revert experiment + reading the skip guards. Cost ~6 min + 1 extra Haiku round (monkeypatch rework). Route: lesson, NOT a skill (#203 probation). **Mitigation for future sessions: every subagent-authored "regression test" MUST be revert-verified by the orchestrator (revert the fix → test must FAIL). A green test on a fixed tree is necessary, not sufficient.** Same root as the 02Z `git rev-parse` stdout-pollution false-pass and Valence's SYNCED-FAIL false-green — subagents over-report pass.
2. **Skip-guarded parametrized tests mask zero-coverage.** `pytest.skip` on fixtures-without-the-needed-shape lets a test show green while never executing its assertion on any param. Hard to distinguish "passed" from "all skipped" without `-v`. Mitigation: prefer a fixture known to have the shape, or assert ≥1 param is non-skip; at minimum run `-v` on new test classes.
3. **Dead drift-gate is silent across consumers.** Health check exit 0 + `✓ Health check passed` while blind to drift (plugin absent in cloud). Already DomI#286 (3rd instance now) — no new routing.

## Next steps

- DomI hub (hour-12): ship the canonical offline-drift helper (#286), then CHILmesh `/sync from DomI` to wire it (kills the dead gate without a 4th copy).
- Issue queue next: #202-2c (pure-Python skeletonize >200s hot loop, the residual) or #201 (.chil format brainstorm). #211 fully closed by this PR once #210 promotes.
- Spec-048 ecosystem remainder is consumer-side (Q-T5/Q-T7) — QuADMesh slot, not CHILmesh.
