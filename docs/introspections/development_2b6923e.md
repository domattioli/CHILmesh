<!-- Session handoff + corpus entry. Caveman style. R5 frontmatter per introspect@DomI tmpl. -->
---
date: 2026-06-15
session: 2026-06-15T10Z-rotation
repo: domattioli/CHILmesh
severity: low
freq: recurring
issues: [211, 48]
wasted_min: 4
wasted_tok: 1500
missing_skill: null
---

# Session Handoff — CHILmesh · development_2b6923e · 2026-06-15 (rotation hour-10)

**Task:** rotation maintenance track — C spec-048 slice already shipped 06-11→; issue-queue work
**Phase:** maintenance
**Progress:** complete — pin synced, quad insert_vertex crash fixed + regression-tested
**Branch:** development (rolling PR #210)
**Duration:** ~30 min
**Tool failures:** 0
**Outcome:** complete

## Pre-flight
- branch_policy_conflict: caught_and_resolved   <!-- harness claude/ecstatic-maxwell-4cr9ax → development per CLAUDE.md precedence; origin/development fetched fresh (didn't exist locally at start) -->
- domi_pin_drift: caught_and_synced   <!-- 3e46639 → 69b073d via sibling-clone update_pin.sh (REPO_ROOT override; CHILmesh has no local skill copy) -->
- caveman_plugin: NOT loaded → Unknown skill → emulated from CLAUDE.md (DomI#268; expected cold-path race)

## What shipped (evidence)

1. `2b6923e^` chore: sync DomI pin 3e46639 → 69b073d (manifest 9d57d1f verified vs DomI main).
2. `2b6923e` fix #211 tail — quad/mixed `insert_vertex` crashed (`np.vstack` ValueError, 3-col new elem vs 4-col mesh) and dropped a quad's 4th edge from cavity detection (wrong sentinel `elem[2]!=elem[3]` + dead ternary `elem[:3] if … else elem[:3]`). Routed ring selection through canonical `_is_triangle()`, added `_ring_to_edges()` (closed 3-/4-edge ring), padded new tris `[v1,v2,new,v1]` on 4-col meshes. Regression test `test_insert_vertex_quad_mesh`.
3. Gate: `pytest tests/ -k "not block_o"` → 998 passed / 47 skipped, 0 regressions; `test_mutations.py` 96 passed.

## Key decisions

1. Picked the #211 tail over priority:now #155 (lifecycle benchmark): #155 is OOM/timeout-bound at scale (env-bound, in-progress); the insert_vertex crash is bounded, reproducible, verifiable, and serves #48 mixed-element robustness directly.
2. Did NOT close #211 — it rides rolling PR #210, not yet on main; closing is on-merge. Same for #202 (problem-1 fixed/locked, 2c debunked, problem-2 build-hook operator-gated).
3. Used existing `_is_triangle()` helper rather than re-deriving a sentinel → cannot drift from the package padding convention again (this is exactly how #211 left these three sites broken).

## What worked (top 3)

1. Empirical repro first — built a 2×2 quad mesh, hit the `np.vstack` ValueError live before touching code → reframed from "dead-ternary cleanup" to "the quad path fully crashes", and proved the regression test fails pre-fix.
2. Sibling-clone pin sync via `REPO_ROOT=/home/user/CHILmesh bash /home/user/DomI/.../update_pin.sh` — no network, no gh, 1 call.
3. Haiku builder for the 1-file fix + test; orchestrator reviewed the full `git diff` and independently re-ran the gate (998 pass) before commit.

## What didn't (pains → routing)

1. **"Fixed" issue had an untested tail that stayed broken.** #211 (`ae77ee1`) patched `mutations.py:80/:488` + `_point_in_element` but missed the *same* wrong-sentinel + dead-ternary in `insert_vertex`/`_get_edge_set` — because the original scan keyed on the literal `"is not a triangle"` strings, and these sites phrase the check differently (`elem[2] != elem[3]` inside a ternary). A grep for the *string* missed them; a grep for the *pattern* `elem\[2\].*elem\[3\]` would have caught all sites at once. Cost: the bug survived a "complete" fix ~2 days. Severity low, freq recurring (sentinel-convention bugs cluster). Route: lesson, NOT a skill (#203 probation). Mitigation: when fixing a convention/sentinel bug, grep the *expression shape* repo-wide (`elem[2]`/`elem[3]` index pairs), not the error string, to find every co-located site.
2. **Quad path was never tested.** `insert_vertex` had 5 tests, all on 3-col triangle fixtures → a hard crash on any quad mesh shipped undetected. Mirrors the #211 root note ("existing tests only cover 3-column pure-triangle fixtures"). No new routing — covered by the standing mixed-element-coverage gap; this session adds one quad case.

## Next steps

- On PR #210 merge: #211 closes (now fully patched incl. insert_vertex), #202 closes for problem-1/2c (problem-2 build-hook stays open, operator-gated).
- Mixed-element test coverage is thin for mutations beyond split/merge — a quad-fixture parametrize sweep over `insert_vertex`/`delete`/`flip` would surface more latent 3-col assumptions. Future C maintenance slot.
- Spec-048 ecosystem remainder unchanged: M-T4 (MADMESHing `quality.py` benchmark → D2) verdict already KEEP-fast-path (06-13 05Z); no open C ecosystem item.
