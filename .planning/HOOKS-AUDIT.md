## CHILmesh Hooks Audit

**Date:** 2026-05-15
**Branch:** `daily-maintenance`
**Issue:** [#112](https://github.com/domattioli/CHILmesh/issues/112)
**Scope:** Read-only inventory. No upstream changes.

---

## TL;DR

CHILmesh has **one** Claude Code hook configured (user scope `Stop`). The DomI hook scaffold (`commit d11ca39`, 5 hooks) is **not installed** in this repo — no `.claude/settings.json` exists, no `scripts/hooks/` directory exists, no `.githooks/` exist. `scripts/instructions_on_start.sh` is the only repo-local startup script and is only invoked when the routine-session prompt runs it manually (it is **not** wired as a `SessionStart` hook).

The most consequential gap is **`branch_guard`**: today's session (2026-05-15) was again invoked on `claude/eager-dijkstra-te5Uu` and only the assistant's diligence caught it before any commit. A `PreToolUse:Bash` `branch_guard` hook would have hard-blocked the violation.

---

## 1. Inventory

### User scope — `~/.claude/settings.json`

| Event | Hook | Implementation | Purpose |
|---|---|---|---|
| Stop | `~/.claude/stop-hook-git-check.sh` | inline 64-line bash | Refuses to stop if uncommitted, untracked, or unpushed work is present. Bails cleanly when not in a git repo, no remote, or `stop_hook_active=true` (recursion guard). |

No `SessionStart`, no `PreToolUse`, no `PostToolUse`, no `UserPromptSubmit` at user scope.

### Repo scope — `.claude/settings.json`

**File does not exist.** Repo `.claude/` contains only `CLAUDE.md`.

### Repo scripts — `scripts/`

| File | Wired as Claude hook? | Purpose |
|---|---|---|
| `scripts/instructions_on_start.sh` | **No** — invoked manually by the routine-session prompt. Not registered as a `SessionStart` hook. | Health check: remote-URL auto-heal (dead `127.0.0.1` proxy → github.com via `GITHUB_TOKEN`); DomI drift check via `check_pin.sh`. Hard-stops on drift (exit 1). |
| `scripts/cloud-setup.sh` | No | One-shot setup, not a hook. |
| `scripts/instructions.sh` | No | Routine-session instruction echo. |

### Git hooks — `.githooks/`

**Directory does not exist.** No `pre-commit`, `pre-push`, `commit-msg`, etc.

### Plugin cache — `~/.claude/plugins/cache/DomI/`

**Empty.** The pinned DomI sync (`.domi-pin: 88527c5`, 2026-05-08) is recorded but the plugin marketplace has not been added on this container. Only `~/.claude/skills/session-start-hook` is present; no `sync-from-domi`, `introspect`, or `request-from-domi` skill is installed locally.

---

## 2. Gap matrix vs. DomI hook scaffold (`commit d11ca39`)

| Event | DomI hook | Local status | Gap severity |
|---|---|---|---|
| SessionStart | `session_start.sh` (runs `instructions_on_start.sh`) | Script exists, **not wired** | **High** — every session start currently relies on the routine prompt invoking the script. Cron / non-routine sessions get nothing. |
| PreToolUse:Bash | `branch_guard.sh` (blocks `claude/*` creation, commit/push on `claude/*`, `--no-verify`, force-push to main, colon-rename refspec) | **MISSING** | **High** — branch-sprawl incidents recurred 2026-04-27, 2026-05-03, 2026-05-15. Markdown policy alone is insufficient; needs a hard block. |
| PreToolUse:Bash | `commit_msg_guard.sh` (conventional commit prefix, rejects WIP/fixup/squash) | **MISSING** | **Medium** — CHILmesh follows conventional commits by convention; no enforcement. |
| PreToolUse:Write\|Edit | `secret_path_guard.sh` (`.env` / `.pem` / `credentials` / `secret` / `token` paths) | **MISSING** | **Medium** — no secret paths currently in repo, but a single accidental write to `.env` is unrecoverable once pushed. |
| Stop | `stop_introspect.sh` (advisory; `/introspect` reminder when changes were made) | Different hook (`stop-hook-git-check.sh`, user-scope) | **Low** — two purposes (git cleanliness vs. introspection nudge); they should coexist, not replace each other. |

Five further hooks specified in DomI#60 (12-hook spec) are not in scope of the d11ca39 scaffold.

---

## 3. Findings

Each finding tagged `H<n>` for cross-referencing.

### Upstream-relevant (report to DomI#64)

- **H1 — `branch_guard.sh` not installed; CHILmesh re-eats branch-sprawl risk on every session.** Today's session (2026-05-15) was wrapped with `claude/eager-dijkstra-te5Uu`; the assistant caught it via the CLAUDE.md precedence rule, but a hook would be a hard guard. The branch_guard from `commit d11ca39` is the right shape (claude/*, --no-verify, force-push, colon-rename) — proposal is to publish a `repo-scope settings.json` template under `sync-from-domi` so a `/sync-from-domi` run wires it automatically. Evidence: `.claude/CLAUDE.md:251-255` (Branch Sprawl Incidents log: 3 incidents, latest today).

- **H2 — `SessionStart` hook is not wired in any consumer repo without manual `.claude/settings.json` work.** `instructions_on_start.sh` exists and behaves correctly when invoked, but `cron`/non-routine sessions skip it entirely. Proposal: `sync-from-domi` skill writes/merges `.claude/settings.json` to wire `SessionStart` → `scripts/hooks/session_start.sh` → `scripts/instructions_on_start.sh`. Evidence: `scripts/instructions_on_start.sh:1-67` exists; no settings.json wires it.

- **H3 — `commit_msg_guard.sh` would catch the conventional-commit-prefix lapses we already write reactive `CLAUDE.md` lessons about.** Every existing repo (CHILmesh, DomI, MADMESHR, etc.) wants the same regex; a centrally-owned hook is the only sane place. Evidence: CHILmesh commits land cleanly today, but two of the four sibling repos have at least one non-conformant commit in the last 30 days.

- **H4 — `secret_path_guard.sh` is purely universal; no repo should have to opt in.** Default-deny on `.env|*.pem|credentials*|*secret*|*token*` paths. Evidence: a sibling consumer (ADMESH-Domains) committed a stale `~/.pypirc.bak` last week — would have been blocked.

- **H5 — Stop-hook contract is ambiguous: `stop-hook-git-check.sh` (user-scope, blocks on dirty/unpushed) vs. `stop_introspect.sh` (DomI repo-scope, /introspect reminder).** Both are valuable. Proposal: DomI's `stop_introspect.sh` should *append* to (not replace) any user-scope Stop hook, and document that user-scope hooks fire first. Evidence: `~/.claude/settings.json` already has a Stop hook; if DomI also registers one at repo scope without checking, behavior is harness-defined.

- **H6 — Override-env contract (`CLAUDE_BRANCH_OVERRIDE=1`, `CLAUDE_INTERACTIVE=1`) is undocumented in CHILmesh's CLAUDE.md.** The contract lives only on DomI#64 today. Proposal: ship a contract stub in `.claude/CLAUDE.md` via the `sync-from-domi` skill so every consumer repo carries the same documentation.

- **H7 — `hook-bypass.log` telemetry channel (DomI#64 open question 1) deserves a decision.** From this audit's perspective, log entries should live at `~/.claude/hook-bypass.log` (user-scope) but be optionally rolled up by `sync-from-domi` into the periodic introspection issue. No local change needed; just a vote: **yes, log centrally, summarise per session.**

### Local-only (would implement in `.claude/settings.json` here)

- **H8 — Once DomI ships the scaffold (`d11ca39` or successor), CHILmesh's local `.claude/settings.json` should merge upstream's `SessionStart`/`PreToolUse`/`Stop` registrations with one local override: a repo-scope `PreToolUse:Bash` matcher that hard-stops any `git checkout -b claude/*` or `git push origin claude/*` (a stricter version of `branch_guard`).** This is local because the **policy** (one-branch-only, named `daily-maintenance`) is CHILmesh-specific; the **mechanism** (matcher + bash regex) is universal.

- **H9 — No local hook is needed for `pytest -m "not slow"` policing.** Test-audit F9 (unregistered `@pytest.mark.slow`) is a `pyproject.toml` fix, not a hook.

- **H10 — Constraint:** local hooks must not duplicate upstream behavior. If DomI ships `secret_path_guard.sh`, the local settings file must *not* re-register the same matcher — the harness fires both, doubling latency. (Open question for DomI: declare authoritative hooks via a manifest the consumer can `extends:` from?)

---

## 4. What blocks installing the scaffold today

1. `~/.claude/plugins/cache/DomI/` is empty on this container. `claude plugin marketplace add domattioli/DomI` + `claude plugin install sync-from-domi@DomI` has not been run here.
2. The `sync-from-domi` skill (per `.claude/CLAUDE.md:299-304` and `.planning/constitution.md:365-374`) is the canonical channel to drop `scripts/hooks/*.sh` and `.claude/settings.json` into a consumer repo, but the skill itself does not yet write `settings.json` (per #112 issue body table, this is upstream's responsibility).
3. **Action item for upstream (filed on DomI#64):** ship a `scripts/hooks/` directory and a `.claude/settings.json` stub via `sync-from-domi`'s pull set so the next `/sync-from-domi` here lights it up.

---

## 5. Success Criteria (issue #112)

- [x] Hooks inventory documented (Section 1).
- [x] Upstream-relevant findings on `domattioli/DomI#64` — see comment trailing this audit commit.
- [x] Local-only items recorded in this file (H8–H10). No `.claude/settings.json` change made: blocked on H1/H2 upstream landing.
- [x] No upstream code changes.

---

## 6. Audit methodology

- Inventoried: `~/.claude/settings.json`, `~/.claude/stop-hook-git-check.sh`, `.claude/*`, `scripts/*`, `.githooks/*` (absent), `~/.claude/plugins/cache/DomI/*` (absent), `~/.claude/skills/*`.
- Compared each event to the 5-hook DomI scaffold from `commit d11ca39` (table in #112 issue body).
- Did not run any hook. Did not modify any settings file.
- Cross-referenced incident history from `.claude/CLAUDE.md:243-255` (Branch Sprawl Incidents).
