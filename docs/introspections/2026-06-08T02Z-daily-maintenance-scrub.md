---
session: 2026-06-08T02Z
repo: CHILmesh (scheduled hour-02) — operator-directed cross-repo task
model: claude-opus-4-8 (+ sonnet/haiku subagents)
branch: development
---

# Session: scrub `daily-maintenance` → `development` across the fleet (#196)

## What happened
- Scheduled CHILmesh routine. Picked branch-policy doc-drift fix: CHILmesh `.claude/CLAUDE.md`
  still named the deprecated `daily-maintenance` as sole working branch (Doc v1.4) 3 wks after
  the #196 migration → observed `development` 8 ahead / **45 behind** `daily-maintenance`
  (real fixes #173/#189/#190/#192 stranded on the deprecated branch).
- Operator then directed (mid-session): "stop referring to daily-maintenance in the
  documentation; apply to the other 7 repos." Scope expanded to fleet-wide active-doc scrub.

## What shipped (6 repos, all on `development`, pushed to rolling PRs)
- **CHILmesh** (#194): `.claude/CLAUDE.md` branch policy + Lessons-Learned; `.specify/speckit-constitution.md`.
- **DomI** (#220): governing docs (`branching.md`, `specs/domi-constitution.md`, TEST-AUDIT,
  code-quality-plan, MANIFEST active skill desc, TOOL_CONFIGURATIONS, IMPLEMENTATION_PLAN),
  10 skill SKILL.md, and **scripts + smoke tests + fixtures in lockstep** (pr-base-validator,
  git-preflight, deploy-state-check, session-resume, git-push-fallback, merge-gracefully,
  introspect, onstart). All smokes green (session-resume 18/18 after fixture fix).
- **ADMESH** (#139): `.claude/CLAUDE.md`, constitution, CONTRIBUTING ×2, docs links, routine copy.
- **QuADMesh** (#80): local routine copy (14 stale directives).
- **Valence** (#131): local routine copy + docs/CONTRIBUTING.
- **MADMESHing** (#53): local routine copy + PROJECT-VISION.
- Health-Wealth-Assistant + XNAT-Interact: **0 refs**, no change.

## Key decisions
- **History preserved, active directives scrubbed.** Did NOT rewrite `docs/introspections/*`,
  `.claude/handoffs/*`, `.planning/*`, `docs/sessions/*`, MANIFEST version-history changelog,
  shipped-spec design records, benchmark.md metric logs, or historical incident writeups in
  SKILL.md. Kept deprecation explainers (branching.md migration section) that intentionally
  name the old branch.
- **Branch reconciliation (daily-maintenance → development merge) flagged operator-only** — 45-commit
  divergence, high regression risk; out of scope for an autonomous session.
- Stale local `.claude/claude_routine_instructions.md` copies fixed by surgical sed (branch refs
  only); they remain stale in OTHER dims (token budgets etc.) vs DomI upstream → flagged for a
  proper `sync from DomI` of that file.

## Pains (tokens_wasted; probation #203 active → corpus only, no skill-request filing)
- pain: subagent completion reports were INACCURATE (claimed files unchanged that were edited,
  and vice versa) → had to review every git diff manually before committing. tokens_wasted: ~4000.
  fix-idea: subagents should emit `git diff --name-only` as ground truth, not prose recall.
- pain: running `scripts/instructions_on_start.sh` in DomI silently REGENERATED
  `docs/introspections/PAIN_MATRIX.md` (auto-gen side-effect, date + wiped a manual section —
  itself the known pain #215) → polluted the diff, nearly committed unrelated churn. Reverted.
  tokens_wasted: ~1500. fix-idea: health check should not mutate tracked files, or `--check` mode.
- pain: a sonnet subagent rewrote a `daily-maintenance` ref inside a MANIFEST historical
  CHANGELOG entry (v2.18 incident description) → falsified history; caught on review + reverted.
  tokens_wasted: ~1000. fix-idea: explicit "version-history/changelog blocks are immutable" rule.
- pain: `grep -r` skipped the per-repo local routine files (tracked regular files) on first sweep
  while the term clearly existed → near-miss on the highest-impact files (14 directives each).
  tokens_wasted: ~800. fix-idea: when a count and a list disagree, re-grep explicit paths.

## Flagged follow-ups (operator / future session)
1. **Reconcile `daily-maintenance` → `development`** in each repo (operator merge; 45-commit gap in CHILmesh).
2. **Full re-sync** of consumer `.claude/claude_routine_instructions.md` from DomI upstream (still stale beyond branch name).
3. **CI workflow triggers** still list `daily-maintenance` on `on: push: branches`: MADMESHing `.github/workflows/benchmark.yml`, Valence `fetch-wnat-hero.yml` + `feature_request.yml`, ADMESH `mkdocs.yml`/`docs.yml`. Left as config (not docs); trivial trigger-list cleanup, deferred.
4. Historical `specs/**` design records across ADMESH/Valence/QuADMesh mention daily-maintenance as past branch — intentionally NOT rewritten (history).

## Next steps
- CHILmesh substantive backlog untouched this session (operator task took priority). Top
  open: #155 (world-mesh benchmark, infra-blocked), #184 (Gmsh label round-trip, deferred sub-scope).
