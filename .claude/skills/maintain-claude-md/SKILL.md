---
name: maintain-claude-md
description: >
  Scaffold, update, and lint CLAUDE.md files across Claude Code repos.
  Use this skill whenever: (1) a repo needs a new CLAUDE.md from scratch
  ("init my CLAUDE.md", "scaffold CLAUDE.md", "set up Claude config for this
  repo"); (2) a session ends and lessons should be recorded ("add a lesson",
  "update CLAUDE.md with what we learned", "log a lesson to CLAUDE.md",
  "append lesson", "record this in CLAUDE.md"); (3) an existing CLAUDE.md
  needs quality-checking ("lint CLAUDE.md", "validate my CLAUDE.md", "check
  for stale references in CLAUDE.md"). Also trigger when user asks to "run
  /maintain-claude-md" in any mode.
---

# maintain-claude-md

Manages the lifecycle of CLAUDE.md — the operational reference file Claude
Code reads at every session start. Three modes: **init**, **update**,
**lint**. Determine which mode the user wants from context, or ask if
ambiguous.

> This skill is self-referencing: every CLAUDE.md it scaffolds includes a
> note that `maintain-claude-md` is the tool used to evolve it, so future
> sessions know how to invoke maintenance.

---

## Mode: init

*When to use:* repo has no CLAUDE.md, or user says "scaffold", "create from
scratch", "reinitialize".

### Step 1 — Introspect the repo

Read the following files (skip if absent, note what's missing):
- `README.md` or `pyproject.toml` / `package.json` — project purpose, stack
- `CONSTITUTION.md` or `.specify/memory/constitution.md` — governing rules
- `PROJECT_PLAN.md` — roadmap, current phase
- Existing `.claude/settings.json` — model, permissions
- Any `specs/` directory — active feature branches

Run `git log --oneline -10` to get a sense of recent activity.

### Step 2 — Ask the user (only what you can't infer)

Keep the interview to ≤5 questions. You likely already know the stack,
language, and purpose from Step 1. Focus on gaps:

1. **Mission** — one sentence: what does this repo do?
2. **Workflow type** — how does Claude primarily help here?
   (issue-triage / data-curation / deployment / autonomous-routine / other)
3. **Key commands** — install, test, build, release (if not in README)
4. **External dependencies** — APIs, databases, cloud environments, required
   env vars that Claude needs to know about
5. **Hard rules** — anything Claude must never do in this repo

### Step 3 — Scaffold the CLAUDE.md

Use the canonical structure below. Populate each section from your
introspection and the interview. Leave explicit `<!-- TODO: fill in -->` 
markers for anything the user hasn't provided — don't invent.

```markdown
# CLAUDE.md

Operational reference for Claude Code sessions on {PROJECT}.

**Read these at every session start (in order):**
{CONSTITUTION_PATH} → {PROJECT_PLAN_PATH} → CLAUDE.md
If CLAUDE.md contradicts the constitution, the constitution wins.

---

## Doc router — for task X, read doc Y

| I need to... | Primary doc | Section |
|---|---|---|
| Know the mission / hard rules | {CONSTITUTION_PATH} | — |
| Run install / test / release | CLAUDE.md (this file) | § Commands |
| Understand architecture | CLAUDE.md | § Architecture |
| Look up a lesson learned | CLAUDE.md | § Lessons learned |
| Find the active spec | specs/{ACTIVE_SPEC}/ | spec.md / plan.md |

---

## Project overview

{MISSION_SENTENCE}
Stack: {STACK}. Workflow: {WORKFLOW_TYPE}.

<!-- maintained-by: maintain-claude-md skill (use /maintain-claude-md update to add lessons, /maintain-claude-md lint to validate) -->

---

## Commands

\`\`\`bash
# Install
{INSTALL_CMD}

# Test
{TEST_CMD}

# Release / deploy
{RELEASE_CMD}
\`\`\`

---

## Architecture

{ASCII_DIAGRAM_OR_FILE_MAP}

Key files:
| File/dir | Role |
|---|---|
| {FILE} | {ROLE} |

---

## Key patterns

{CODE_IDIOMS_OR_CONVENTIONS}

---

## Lessons learned

<!-- Append entries with /maintain-claude-md update. Do not edit manually. -->
<!-- Format: ### YYYY-MM-DD — {topic} -->

(No entries yet.)

---

## Session history

<!-- One-liner per session, newest first. Updated by /maintain-claude-md update. -->

(No entries yet.)
```

### Step 4 — Write and confirm

Write the file to `CLAUDE.md` at the repo root. Show the user a diff-style
summary of what was generated. Ask: "Does this look right? Anything to add
before I write it?"

---

## Mode: update

*When to use:* user wants to append a lesson, session summary, or other
new content without clobbering existing CLAUDE.md.

### What to append and where

**Lessons learned** (`## Lessons learned` section):
```
### {YYYY-MM-DD} — {topic}
{2–5 bullet points. Be specific: what failed, what worked, why it matters.}
```

**Session history** (`## Session history` section):
```
- {YYYY-MM-DD}: {one-liner summary of what shipped/changed}
```

**Speckit active feature** (if a new spec branch was just started):
Update the `<!-- SPECKIT START/END -->` block near the top:
```markdown
<!-- SPECKIT START -->
Active spec-kit feature: `{SPEC_ID}-{name}` (branch `{SPEC_ID}-{name}`).
For context, read: `specs/{SPEC_ID}-{name}/plan.md`
<!-- SPECKIT END -->
```

### Rules for safe appends

- Never rewrite existing sections, only append to them.
- Preserve all existing content verbatim — including typos, formatting
  choices, and comment blocks.
- When appending to `## Lessons learned`, insert the new entry *after* the
  section heading and the `<!-- Append entries... -->` comment, before any
  existing entries (newest-first ordering).
- After writing, show the user the lines you added (a brief diff).

---

## Mode: lint

*When to use:* user wants to validate an existing CLAUDE.md.

### Checks to run

Run all checks, report as a table: **check | status | detail**.

| Check | How to verify |
|---|---|
| Required sections present | `## Doc router`, `## Project overview`, `## Commands`, `## Architecture`, `## Lessons learned` |
| Doc router links resolve | For each path in the table, check the file exists |
| Constitution path correct | Verify `.specify/memory/constitution.md` or `CONSTITUTION.md` exists |
| Speckit block accurate | If `<!-- SPECKIT START/END -->` present, verify the referenced spec folder exists |
| No broken skill references | Scan for `.claude/skills/*/SKILL.md` paths — verify each exists |
| Commands work | Run each command snippet in a dry-run / syntax check (don't execute if destructive) |
| Lessons are dated | Every `###` heading in `## Lessons learned` should start `### YYYY-MM-DD` |
| No stale session history | Flag entries > 90 days old as "review: may be stale" |
| Self-reference comment present | `<!-- maintained-by: maintain-claude-md skill -->` exists in `## Project overview` |

### Output format

```
CLAUDE.md lint report — {repo} — {date}

✅ PASS   Required sections: all 5 present
⚠️  WARN   Doc router: specs/002-feature/plan.md not found (spec may be closed)
❌ FAIL   Constitution path: CONSTITUTION.md missing (.specify/memory/constitution.md also missing)
✅ PASS   Skills: all referenced paths exist
✅ PASS   Lessons: 4 entries, all dated
⚠️  WARN   Session history: 2 entries older than 90 days
```

Exit with a suggested fix for each ❌ and ⚠️.

---

## Cross-cutting notes

**Self-awareness**: Every CLAUDE.md this skill generates includes a
`<!-- maintained-by -->` comment so future sessions know to invoke
`/maintain-claude-md` for upkeep. This closes the loop: the skill's output
teaches Claude how to re-invoke the skill.

**Extensibility**: The structure above is the canonical base. Individual
repos may add sections (e.g. `## Domain insights`, `## Release skills`).
The lint mode only enforces the base set; extra sections are never flagged.

**Conflict resolution**: If a repo's CONSTITUTION.md or `.specify/memory/
constitution.md` contradicts what this skill would scaffold, the constitution
wins. Surface the conflict to the user before writing anything.
