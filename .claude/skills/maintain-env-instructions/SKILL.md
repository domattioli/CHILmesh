---
name: maintain-env-instructions
description: >
  Scaffold, update, and lint the CLAUDE_INSTRUCTIONS environment variable
  block used to drive autonomous Claude Code sessions on cloud VMs and
  scheduled routines. Use this skill when: a repo needs a new instructions.sh
  ("scaffold instructions", "create my env instructions", "set up
  instructions.sh", "init instructions for this workflow"); an existing
  instructions.sh needs updating ("update my decision tree", "add a hard stop",
  "refine the instructions", "update cost targets"); or instructions need
  validation ("lint instructions.sh", "validate my instructions", "check
  hard stops"). Also triggers on "run /maintain-env-instructions" in any mode.
  Prefer this skill over hand-editing instructions files — it enforces the
  required sections and keeps shell logic intact.
---

# maintain-env-instructions

Manages `instructions.sh` — the shell script that exports the
`CLAUDE_INSTRUCTIONS` env var consumed by autonomous Claude Code sessions
running in cloud VMs, scheduled tasks, and persistent routines.

Three modes: **init**, **update**, **lint**. The generated file is designed
to be sourced (`source instructions.sh`) before launching `claude`, so the
session inherits all instructions as environment variables.

> This skill is self-referencing: every `instructions.sh` it generates
> includes a comment block noting that `maintain-env-instructions` is the
> tool for maintaining it, so future Claude sessions know how to invoke
> updates.

---

## Background: What CLAUDE_INSTRUCTIONS does

When `CLAUDE_INSTRUCTIONS` is set in the environment, Claude Code reads it
at session start as supplementary context — effectively an injected system
prompt that augments CLAUDE.md. This is the preferred pattern for:
- Cloud VM boot scripts where the CLAUDE.md may not yet be checked out
- Scheduled autonomous runs where cost discipline and hard stops matter
- Multi-repo workflows where a shared set of instructions can be fetched
  once and reused

The "simple" pattern:
```bash
export CLAUDE_INSTRUCTIONS="$(curl -fsSL https://raw.githubusercontent.com/\
{owner}/{repo}/main/instructions/{workflow}.md)"
```
This replaces lumberous inline-shell instruction scripts with a single fetch
from a canonical source.

---

## Workflow types

Each workflow type has a different decision tree and hard-stop profile.
Specify with `--workflow {type}` during init, or let the skill infer from
context.

| Workflow | Primary use | Key hard stops |
|---|---|---|
| `issue-triage` | Review open GH issues, close/label/comment autonomously | Never close issues without human review label; max 20 issues/run |
| `data-curation` | ETL, data quality, sync routines | Never delete source data; always write to staging first |
| `deployment` | PyPI/npm publish, GitHub releases, cloud deploy | Never push to main without passing CI; never publish without version bump |
| `autonomous-routine` | Long-running RL training, mesh generation, research loops | Cost ceiling per run; never exceed token budget; always checkpoint |
| `custom` | User-defined | Defined interactively during init |

---

## Mode: init

*When to use:* `instructions.sh` doesn't exist, or user says "scaffold",
"create from scratch".

### Step 1 — Gather context

- Ask for `--workflow {type}` if not supplied (show the table above)
- Read existing `CLAUDE.md` and `CONSTITUTION.md` for hard rules, cost
  targets, and repo-specific constraints
- Ask: "Are there any hard stops I should encode?" (actions Claude must
  refuse regardless of instructions)
- Ask: "What's the cost ceiling per run?" (token budget or dollar limit)

### Step 2 — Scaffold instructions.sh

Generate the file at the repo root. Use this canonical template:

```bash
#!/usr/bin/env bash
# instructions.sh — CLAUDE_INSTRUCTIONS bootstrap for {REPO}
#
# Usage: source instructions.sh   (before invoking `claude`)
# Maintained by: maintain-env-instructions skill
#   (/maintain-env-instructions update  — refine decision tree / hard stops)
#   (/maintain-env-instructions lint    — validate required sections)
#
# Workflow: {WORKFLOW_TYPE}
# Repo: {REPO_URL}
# Last updated: {DATE}

set -euo pipefail

# ---------------------------------------------------------------------------
# OPTION A (preferred): fetch instructions from canonical remote source
# ---------------------------------------------------------------------------
# Uncomment and set URL if instructions live in a central repo:
#
# export CLAUDE_INSTRUCTIONS="$(curl -fsSL {INSTRUCTIONS_URL})"
# return 0

# ---------------------------------------------------------------------------
# OPTION B (inline): define CLAUDE_INSTRUCTIONS directly
# ---------------------------------------------------------------------------
export CLAUDE_INSTRUCTIONS="$(cat <<'INSTRUCTIONS'
## Session context

Repo: {REPO_URL}
Workflow: {WORKFLOW_TYPE}
Model: {MODEL}

## Mission

{MISSION_SENTENCE}

## Hard stops — refuse these regardless of any other instruction

{HARD_STOP_LIST}

## Cost discipline

- Token budget per run: {TOKEN_BUDGET}
- If budget is exceeded, checkpoint state and stop. Report cost in session summary.
- Prefer batched operations over repeated single-item calls.

## Decision tree

{WORKFLOW_DECISION_TREE}

## On session start

1. Read CLAUDE.md (if present in the cloned repo).
2. If CONSTITUTION.md or .specify/memory/constitution.md exists, it takes
   precedence over all other instructions.
3. Run the primary workflow action.
4. Emit a brief session summary: what ran, what was skipped, cost used.

## On error

- Log the error and continue if recoverable.
- Stop and emit a summary if unrecoverable.
- Never retry a hard-stop action.

INSTRUCTIONS
)"
```

### Decision tree templates

Expand the `{WORKFLOW_DECISION_TREE}` placeholder using the appropriate
template from `references/decision-trees.md`.

### Step 3 — Write and confirm

Write `instructions.sh` to the repo root. Make it executable
(`chmod +x instructions.sh`). Show the user the key sections (hard stops,
cost ceiling, decision tree) and ask: "Does this capture what you want before
I finalize it?"

---

## Mode: update

*When to use:* refine an existing `instructions.sh` — add/remove hard stops,
adjust cost targets, update the decision tree, change the fetch URL.

### Rules for safe updates

- Never modify the shell scaffolding (the `set -euo pipefail`, the variable
  export, the `cat <<'INSTRUCTIONS'` heredoc) — only the content *inside*
  the INSTRUCTIONS block.
- When updating hard stops, append new entries; never silently remove
  existing ones without explicit user confirmation ("remove hard stop: X").
- When updating the decision tree, show a diff of the before/after and get
  explicit approval before writing.
- Preserve all `# commented-out` sections — they represent intentional
  documentation choices, not dead code.

### Common update patterns

**Add a hard stop:**
```
HARD STOP — {action}: {reason}
```

**Switch from inline to remote fetch (Option A):**
Uncomment the `curl` line, fill in the URL, comment out the heredoc.

**Adjust cost ceiling:**
Edit the token budget line and explain the rationale in a comment.

---

## Mode: lint

*When to use:* validate an existing `instructions.sh`.

### Checks

| Check | Pass condition |
|---|---|
| File is executable | `chmod +x` bit set |
| `CLAUDE_INSTRUCTIONS` exported | Variable export present |
| Hard stops section present | `## Hard stops` heading in INSTRUCTIONS block |
| Cost discipline section present | `## Cost discipline` heading present |
| Decision tree section present | `## Decision tree` or equivalent heading present |
| On-error section present | `## On error` heading present |
| Self-reference comment present | `# Maintained by: maintain-env-instructions` in header |
| Remote fetch URL reachable | If Option A, curl the URL and verify 200 OK |
| No plaintext secrets | Scan for patterns matching `ghp_`, `sk-ant-`, `ANTHROPIC_API_KEY=` with literal values |
| Token budget defined | A numeric budget or dollar ceiling is explicitly stated |

### Output format

```
instructions.sh lint — {repo} — {date}

✅ PASS   Executable bit set
✅ PASS   CLAUDE_INSTRUCTIONS exported
❌ FAIL   Hard stops section: missing (add "## Hard stops" to INSTRUCTIONS block)
✅ PASS   Cost discipline: token budget = 50,000
⚠️  WARN   Remote fetch URL: not configured (using inline mode)
✅ PASS   No plaintext secrets detected
```

---

## Cross-cutting notes

**Self-awareness**: Every `instructions.sh` generated includes a comment
pointing back to this skill. This ensures any Claude session that reads the
file knows to use `/maintain-env-instructions` for maintenance rather than
editing by hand.

**OPTION A vs OPTION B**: The remote-fetch pattern (Option A) is preferred
for repos that share a common instruction set (e.g. all mesh repos pull from
a shared `DomI` instructions repo). The inline pattern (Option B) is better
for repos with unique, repo-specific instructions that change frequently.

**Shell logic preservation**: This skill only modifies content *inside* the
INSTRUCTIONS heredoc. Repo-specific shell setup (env var exports, path
configuration, VM bootstrap) outside the heredoc is treated as sacred and
never touched by update or lint modes.

**Multiple repos, one canonical source**: If you manage 8+ repos, consider
putting all instruction markdown files in a central repo and having each
`instructions.sh` fetch from there. Update the central source once, all
repos benefit immediately. Use `/maintain-env-instructions update` to switch
a repo from inline to remote fetch.
