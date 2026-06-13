<!-- Session handoff + corpus entry. Caveman style. R5 frontmatter per introspect@DomI tmpl. -->
---
date: 2026-06-13
session: 2026-06-13T14Z-rotation
repo: domattioli/CHILmesh
severity: med
freq: recurring
issues: [202, 48]
wasted_min: 0
wasted_tok: 0
missing_skill: null
---

# Session Handoff — CHILmesh · development_f06a0ad · 2026-06-13 (rotation hour-14)

**Task:** overhaul rotation, CHILmesh slot. C spec-048 slice shipped prior slots → maintenance track (issue-queue top) + hub loop.
**Phase:** maintenance
**Progress:** complete — #202 problem-2c debunked + fixed, pushed to rolling PR #210
**Branch:** development (rolling PR #210)
**Duration:** ~35 min
**Tool failures:** 1 minor (pytest absent → pip install; chilmesh not pre-installed → editable install)
**Outcome:** complete

## Pre-flight

- branch_policy_conflict: caught_and_resolved — harness `claude/determined-gauss-ogyo6q` → `development` per CLAUDE.md precedence. Harness branch == origin/main baseline.
- domi_pin_drift: none on `development` (3e46639 = DomI main HEAD; pin on the harness branch reads stale 39fd74a but that's pre-02Z-promotion, not real drift).
- caveman_plugin: NOT loaded → `/caveman:caveman ultra` returned `Unknown skill` → emulated from SKILL.md (honest fallback per #168). SessionStart resume hook NOT active.
- health_check: exit 0 but DEAD GATE — `sync-from-domi not installed` warn+continue (DomI#286 class). Pin already current → no exposure.

## What shipped (evidence)

1. `f06a0ad` fix #202 problem-2c. The 02Z slow-path `UserWarning` hardcoded "Block_O ~5k elems can exceed 200s". Re-measured: Block_O full pure-Python init = **0.24s** (n_layers=9), `_skeletonize`=0.018s. Synthetic scaling LINEAR (4k→0.05s · 16k→0.23s · 60k→1.0s) — no O(n²) anywhere. Claim ~1000× wrong.
2. Edits (1 src file + 1 test comment + 1 new test): corrected warning text (now "linear, ~1s/60k elems") + raised over-eager 2k→**50k** elem threshold; vectorized `_get_centroids` (measured 0.135s hot spot → bit-identical fancy-index mean, incl. padded `[v0,v1,v2,v0]`); `tests/test_skeletonize_perf.py` Block_O 30s regression tripwire; refreshed stale conftest comment ("~30s / O(n²)" → reality).
3. Gate: full suite **1069 passed / 56 skipped, 0 regressions** (100s). Verified centroid bit-identity (max|diff|=0.0) + warning fire@60k / no-fire@5k independently before commit.
4. Hub/#48: #202 evidence comment (recommend MADMESHing drop `MADMESHING_RUN_BLOCK_O=1` gate-skip); #48 checklist.

## Key decisions

1. Picked #202-2c (08Z "next steps" nominee, queue top, on-mission for #48: it's *why* MADMESHing env-gates Block_O) over the brainstorm/research backlog (#201/#155/#167).
2. **Re-measured before trusting the claim.** The ">200s" was inherited from the issue body (2026-06-09) into the 02Z warning text AND the conftest comment AND the 08Z handoff "next steps" — three layers deep, never re-measured. First action was a cProfile, which killed the premise in one shot. Pivoted slice from "fix a 200s hot loop" to "correct a 1000×-wrong perf claim" — the real bug.
3. **Left #202 OPEN** for its real half (problem 2: editable/source install ships no compiled cpp extension; `CPP_AVAILABLE` wiring, tracked w/ #163). Only the perf-cliff sub-problem is resolved. Did not over-close.
4. **Did NOT port DomI#286 offline drift-gate** (4th hand-rolled copy = the anti-pattern #48 targets). Pin current → no exposure → held for DomI hub.

## What worked (top 3)

1. **Profile-first killed a stale premise.** ~5 min of cProfile + a 3-point scaling test turned "deep O(n²) hot-loop fix (risky/deep)" into "fix a wrong string + a threshold (low-risk)". Cheapest possible course-correction.
2. **Orchestrator-side independent verification.** Did not trust Haiku's "applied cleanly" — ran a reference-loop centroid diff (proved bit-identity) + a warning fire/no-fire harness at the new threshold. Caught nothing wrong this time, but the check is the gate, not the trust.
3. Haiku builder / Fable-5 review split: exact old/new strings given → one clean round, zero rework.

## What didn't (pains → routing)

1. **Stale perf claims propagate across artifacts without re-measurement (recurring, severity med).** "Block_O >200s" originated in #202's body, then got copied into a *user-facing `UserWarning`* (02Z), a conftest comment, and a handoff "next steps" — a 1000×-wrong number shipped to users because each consumer trusted the prior layer instead of running the 5-min profile. Same family as the 08Z subagent-false-pass pain: **assertions about runtime/behavior get propagated, not verified.** Route: lesson, NOT a skill (#203 probation). Mitigation for future sessions: **any perf/runtime claim baked into shipped text (warnings, docs, comments) MUST cite a re-measurement with date+commit; treat an un-dated magnitude claim as unverified.** A cheap profile beats inheriting a number.
2. **Dead drift-gate still silent** (DomI#286, now Nth instance). Health exit 0 + `✓` while blind to drift (plugin absent in cloud). No new routing — held for #286 canonical helper.
3. Container ships no chilmesh/pytest — every CHILmesh routine slot pays the editable-install + `pip install pytest` tax. Minor, but recurring across slots; a cached venv / setup step would save ~30s/slot. (Same class as QuADMesh `scripts/dev_setup.sh`.)

## Next steps

- **MADMESHing M slot:** drop `MADMESHING_RUN_BLOCK_O=1` gate-skip + un-skip the two Block_O tests — basis cleared (init sub-second). MADMESHing-repo edit, not CHILmesh.
- DomI hub: ship canonical offline-drift helper (#286) → CHILmesh `/sync from DomI` to wire it (kills the dead gate without a 4th copy).
- #202 stays open for cpp source-build wiring (problem 2, w/ #163). Queue next after that: #201 (.chil format brainstorm) or #198 (hero gif).
- #202-2c fully closed once #210 promotes to main.
